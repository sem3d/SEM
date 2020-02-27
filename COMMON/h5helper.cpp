
#include "h5helper.h"
#include "sem_types.h"
#include <cassert>
#include <vector>
#include <cstdlib>

using std::vector;

#undef  _MIN
#define _MIN(x,y) ((x) < (y) ? (x) : (y))

#undef  _MAX
#define _MAX(x,y) ((x) > (y) ? (x) : (y))


#define H5H_MSG_SIZE 64
static herr_t h5h_err_handler(hid_t stack, void* client_data)
{
    exit(1);
}

static void* old_handler;
static void* old_data;
void h5h_set_errhandler()
{
    H5Eget_auto(H5E_DEFAULT, (H5E_auto_t*)&old_handler, &old_data);
    H5Eset_auto(H5E_DEFAULT, h5h_err_handler, NULL);
}

void h5h_restore_errhandler()
{
    H5Eset_auto(H5E_DEFAULT, (H5E_auto_t)old_handler, old_data);
}

void h5h_create_attr(hid_t parent, const char* name, int value)
{
    hid_t attr_id, space_id;
    space_id = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate(parent, name, H5T_STD_I32LE, space_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    H5Sclose(space_id);
}

void h5h_create_attr(hid_t parent, const char* name, double value)
{
    hid_t attr_id, space_id;
    space_id = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr_id);
    H5Sclose(space_id);
}

#define _CHUNK_SZ  (256*1024)

hid_t h5h_dset_prop(hsize_t d0, hsize_t d1)
{
    // heuristic to find a good chunk size of about 256k
    hsize_t chunk[2];
    if (d0*d1<_CHUNK_SZ) {
	chunk[0] = d0;
	chunk[1] = d1;
    } else {
	chunk[0] = _MIN(d0, _CHUNK_SZ);
	chunk[1] = _MAX(1, _MIN(d1, _CHUNK_SZ/chunk[0]));
    }
    hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
    if (d0*d1>128) {
	H5Pset_deflate(prop_id, 5);
	H5Pset_chunk(prop_id, 2, chunk);
    }
    return prop_id;
}

hid_t h5h_dset_prop(hsize_t d0)
{
    // heuristic to find a good chunk size of about 256k
    hsize_t chunk[1];
    chunk[0] = _MIN(d0, _CHUNK_SZ);
    hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
    if (d0>128) {
	H5Pset_deflate(prop_id, 5);
	H5Pset_chunk(prop_id, 1, chunk);
    }
    return prop_id;
}


/** WRITE DSET **/

// Use specialized template functions to get HDF5 type from T (to avoid C++11 is_same).
template<typename T> hid_t h5h_get_type();
template<> hid_t h5h_get_type<double>() {return H5T_NATIVE_DOUBLE;}
template<> hid_t h5h_get_type<int>() {return H5T_NATIVE_INT;}
template<> hid_t h5h_get_type<int64_t>() {return H5T_NATIVE_INT64;}

// Template : 1D array
template<typename T> void h5h_write_dset(hid_t parent, const char* name, int d0, const T* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[1];
    dims[0] = d0;
    if (d0==0) return;
    space_id = H5Screate_simple(1, dims, dims);
    prop_id = h5h_dset_prop(d0);
    dset_id = H5Dcreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, h5h_get_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}
// Instanciate functions from template.
template void h5h_write_dset<double>(hid_t parent, const char* name, int d0, const double* arr);
template void h5h_write_dset<int>(hid_t parent, const char* name, int d0, const int* arr);
template void h5h_write_dset<index_t>(hid_t parent, const char* name, int d0, const index_t* arr);

template<typename T> void h5h_write_dset(hid_t parent, const char* name, const std::vector<T>& arr)
{
    h5h_write_dset(parent, name, arr.size(), &arr[0]);
}
// Instanciate functions from template.
template void h5h_write_dset<double>(hid_t parent, const char* name, const std::vector<double>& arr);
template void h5h_write_dset<int>(hid_t parent, const char* name, const std::vector<int>& arr);
template void h5h_write_dset<index_t>(hid_t parent, const char* name, const std::vector<index_t>& arr);

// Templete : 1D empty array

template<typename T> void h5h_write_dset_empty(hid_t parent, const char* name, int d0, const T* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[1];
    dims[0] = d0;
    //if (d0==0) return;
    space_id = H5Screate_simple(1, dims, NULL);
    prop_id = h5h_dset_prop(d0);
    dset_id = H5Dcreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, h5h_get_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}
// Instanciate functions from template.
template void h5h_write_dset_empty<double>(hid_t parent, const char* name, int d0, const double* arr);
template void h5h_write_dset_empty<int>(hid_t parent, const char* name, int d0, const int* arr);
template void h5h_write_dset_empty<index_t>(hid_t parent, const char* name, int d0, const index_t* arr);

template<typename T> void h5h_write_dset_empty(hid_t parent, const char* name, const std::vector<T>& arr)
{
    h5h_write_dset_empty(parent, name, arr.size(), &arr[0]);
}
// Instanciate functions from template.
template void h5h_write_dset_empty<double>(hid_t parent, const char* name, const std::vector<double>& arr);
template void h5h_write_dset_empty<int>(hid_t parent, const char* name, const std::vector<int>& arr);
template void h5h_write_dset_empty<index_t>(hid_t parent, const char* name, const std::vector<index_t>& arr);
// Template : 2D array
template<typename T> void h5h_write_dset_2d(hid_t parent, const char* name, int d0, int d1, const T* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[2];
    dims[0] = d0;
    dims[1] = d1;
    if (d0*d1==0) return;
    space_id = H5Screate_simple(2, dims, dims);
    prop_id = h5h_dset_prop(d0, d1);
    dset_id = H5Dcreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, h5h_get_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}
// Instanciate functions from template.
template void h5h_write_dset_2d<double>(hid_t parent, const char* name, int d0, int d1, const double* arr);
template void h5h_write_dset_2d<int>(hid_t parent, const char* name, int d0, int d1, const int* arr);
template void h5h_write_dset_2d<index_t>(hid_t parent, const char* name, int d0, int d1, const index_t* arr);

template<typename T> void h5h_write_dset_2d(hid_t parent, const char* name, int d1, const vector<T>& v)
{
    h5h_write_dset_2d(parent, name, v.size()/d1, d1,  &v[0]);
}
// Instanciate functions from template.
template void h5h_write_dset_2d<double>(hid_t parent, const char* name, int d1, const std::vector<double>& arr);
template void h5h_write_dset_2d<int>(hid_t parent, const char* name, int d1, const std::vector<int>& arr);
template void h5h_write_dset_2d<index_t>(hid_t parent, const char* name, int d1, const std::vector<index_t>& arr);

// Template write dataset 2d
template<typename T> void h5h_write_dset_2d_empty(hid_t parent, const char* name, int d0, int d1, const T* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[2];
    dims[0] = d0;
    dims[1] = d1;
    space_id = H5Screate_simple(2, dims, NULL);
    prop_id = h5h_dset_prop(d0, d1);
    dset_id = H5Dcreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, h5h_get_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}
// Instanciate functions from template.
template void h5h_write_dset_2d_empty<double>(hid_t parent, const char* name, int d0, int d1, const double* arr);
template void h5h_write_dset_2d_empty<int>(hid_t parent, const char* name, int d0, int d1, const int* arr);
template void h5h_write_dset_2d_empty<index_t>(hid_t parent, const char* name, int d0, int d1, const index_t* arr);

template<typename T> void h5h_write_dset_2d_empty(hid_t parent, const char* name, int d1, const vector<T>& v)
{
    h5h_write_dset_2d_empty(parent, name, v.size()/d1, d1,  &v[0]);
}
// Instanciate functions from template.
template void h5h_write_dset_2d_empty<double>(hid_t parent, const char* name, int d1, const std::vector<double>& arr);
template void h5h_write_dset_2d_empty<int>(hid_t parent, const char* name, int d1, const std::vector<int>& arr);
template void h5h_write_dset_2d_empty<index_t>(hid_t parent, const char* name, int d1, const std::vector<index_t>& arr);

/** ATTRIBUTES **/

int h5h_read_attr_int(hid_t dset_id, const char* attrname)
{
    //hsize_t dims[1] = {1,};
    int res;
    hid_t attid = H5Aopen(dset_id, attrname, H5P_DEFAULT);
    //hid_t spcid = H5Aget_space(attid);
    //TODO check dims of attr
    H5Aread(attid, H5T_NATIVE_INT, &res);
    H5Aclose(attid);
    return res;
}

void h5h_write_attr_int(hid_t dset_id, const char* attrname, int val)
{
    //hsize_t dims[1] = {1,};
    hid_t spcid = H5Screate(H5S_SCALAR);
    hid_t attid = H5Acreate2(dset_id, attrname, H5T_STD_I64LE, spcid, H5P_DEFAULT, H5P_DEFAULT);

    //TODO check dims of attr
    H5Awrite(attid, H5T_NATIVE_INT, &val);
    H5Aclose(attid);
    H5Sclose(spcid);
}

int h5h_get_dset1d_size(hid_t dset_id)
{
    hsize_t dims[1], maxdims[1];
    hid_t space_id = H5Dget_space(dset_id);
    int ndims = H5Sget_simple_extent_ndims(space_id);
    assert(ndims==1);
    H5Sget_simple_extent_dims(space_id, dims, maxdims);
    H5Sclose(space_id);
    return dims[0];
}

void h5h_get_dset2d_size(hid_t dset_id, hsize_t& d1, hsize_t& d2)
{
    hsize_t dims[2], maxdims[2];
    hid_t space_id = H5Dget_space(dset_id);
    int ndims = H5Sget_simple_extent_ndims(space_id);
    assert(ndims==1 || ndims==2);
    H5Sget_simple_extent_dims(space_id, dims, maxdims);
    H5Sclose(space_id);
    d1 = (ndims == 1) ?       1 : dims[0];
    d2 = (ndims == 1) ? dims[0] : dims[1];
}

void h5h_read_dset_Nx2(hid_t g, const char* dname, vector<double>& v, vector<double>& w)
{
    hid_t dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    hsize_t dim1 = -1, dim2 = -1;
    h5h_get_dset2d_size(dset_id, dim1, dim2);

    v.resize(dim1); w.resize(dim1);
    hid_t memspace_id = H5Screate_simple(1, &dim1, NULL);
    hsize_t startmem[1] = {0}; hsize_t countmem[1] = {dim1};
    H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, startmem, NULL, countmem, NULL);

    hid_t filespace_id = H5Dget_space(dset_id);
    hsize_t countfile[2] = {dim1, 1}; hsize_t startfile[2] = {0, 0};
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, startfile, NULL, countfile, NULL); // Get X, mask Y
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, &v[0]);
    startfile[1] = 1;
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, startfile, NULL, countfile, NULL); // Get Y, mask X
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, &w[0]);

    H5Sclose(memspace_id);
    H5Sclose(filespace_id);
    H5Dclose(dset_id);
}

void h5h_read_dset_Nx3(hid_t g, const char* dname, vector<double>& u, vector<double>& v, vector<double>& w)
{
    hid_t dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    hsize_t dim1 = -1, dim2 = -1;
    h5h_get_dset2d_size(dset_id, dim1, dim2);

    u.resize(dim1);
    v.resize(dim1);
    w.resize(dim1);
    hid_t memspace_id = H5Screate_simple(1, &dim1, NULL);
    hsize_t startmem[1] = {0}; hsize_t countmem[1] = {dim1};
    H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, startmem, NULL, countmem, NULL);

    hid_t filespace_id = H5Dget_space(dset_id);
    hsize_t countfile[2] = {dim1, 1}; hsize_t startfile[2] = {0, 0};

    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, startfile, NULL, countfile, NULL); // Get X, mask Y
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, &u[0]);

    startfile[1] = 1;
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, startfile, NULL, countfile, NULL); // Get Y, mask X
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, &v[0]);

    startfile[1] = 2;
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, startfile, NULL, countfile, NULL); // Get Y, mask X
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, &w[0]);

    H5Sclose(memspace_id);
    H5Sclose(filespace_id);
    H5Dclose(dset_id);
}
void h5h_read_dset(hid_t g, const char* dname, vector<int>& v)
{
    hid_t dset_id;
    int dim;
    dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    dim = h5h_get_dset1d_size(dset_id);
    v.resize(dim);
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    H5Dclose(dset_id);
}

void h5h_read_dset_2d(hid_t g, const char* dname, int& d0, int& d1, vector<double>& data)
{
    hid_t dset_id;
    hsize_t dim0, dim1;

    dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    h5h_get_dset2d_size(dset_id, dim0, dim1);
    d0 = dim0;
    d1 = dim1;
    data.resize(dim0*dim1);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    H5Dclose(dset_id);
}

void h5h_read_dset_2d(hid_t g, const char* dname, int& d0, int& d1, vector<int>& data)
{
    hid_t dset_id;
    hsize_t dim0, dim1;

    dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    h5h_get_dset2d_size(dset_id, dim0, dim1);
    d0 = dim0;
    d1 = dim1;
    data.resize(dim0*dim1);
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    H5Dclose(dset_id);
}

void h5h_read_dset_2d(hid_t g, const char* dname, int& d0, int& d1, vector<int64_t>& data)
{
    hid_t dset_id;
    hsize_t dim0, dim1;

    dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    h5h_get_dset2d_size(dset_id, dim0, dim1);
    d0 = dim0;
    d1 = dim1;
    data.resize(dim0*dim1);
    H5Dread(dset_id, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    H5Dclose(dset_id);
}


#undef _MIN
#undef _MAX


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
