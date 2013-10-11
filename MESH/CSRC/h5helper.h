//

#ifndef _H5HELPER_H_
#define _H5HELPER_H_

// hdf5.h may  include mpi.h which can include mpicxx.h and we don't want that since
// we would need to link against libmpicxx ...
#define OMPI_SKIP_MPICXX
#include "hdf5.h"

#undef  _MIN
#define _MIN(x,y) ((x) < (y) ? (x) : (y))

#undef  _MAX
#define _MAX(x,y) ((x) > (y) ? (x) : (y))


static inline void h5h_create_attr(hid_t parent, const char* name, int value)
{
    hid_t attr_id, space_id;
    space_id = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate(parent, name, H5T_STD_I32LE, space_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    H5Sclose(space_id);
}

static inline void h5h_create_attr(hid_t parent, const char* name, double value)
{
    hid_t attr_id, space_id;
    space_id = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr_id);
    H5Sclose(space_id);
}

#define _CHUNK_SZ  (256*1024)

static inline hid_t h5h_dset_prop(int d0, int d1)
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

static inline hid_t h5h_dset_prop(int d0)
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

static inline void h5h_write_dset(hid_t parent, const char* name, int d0, int d1, double* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[2];
    dims[0] = d0;
    dims[1] = d1;
    space_id = H5Screate_simple(2, dims, dims);
    prop_id = h5h_dset_prop(d0, d1);
    dset_id = H5Dcreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}

static inline void h5h_write_dset(hid_t parent, const char* name, int d0, int d1, int* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[2];
    dims[0] = d0;
    dims[1] = d1;
    space_id = H5Screate_simple(2, dims, dims);
    prop_id = h5h_dset_prop(d0, d1);
    dset_id = H5Dcreate(parent, name, H5T_STD_I32LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}

static inline void h5h_write_dset(hid_t parent, const char* name, int d0, double* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[1];
    dims[0] = d0;
    space_id = H5Screate_simple(1, dims, dims);
    prop_id = h5h_dset_prop(d0);
    dset_id = H5Dcreate(parent, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}

static inline void h5h_write_dset(hid_t parent, const char* name, std::vector<double>& arr)
{
    h5h_write_dset(parent, name, arr.size(), &arr[0]);
}

static inline void h5h_write_dset(hid_t parent, const char* name, int d0, int* arr)
{
    hid_t dset_id, space_id, prop_id;
    hsize_t dims[1];
    dims[0] = d0;
    space_id = H5Screate_simple(1, dims, dims);
    prop_id = h5h_dset_prop(d0);
    dset_id = H5Dcreate(parent, name, H5T_STD_I32LE, space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
    H5Dclose(dset_id);
    H5Pclose(prop_id);
    H5Sclose(space_id);
}

#undef _MIN
#undef _MAX

#endif
// Local Variables:
// mode: c++
// coding: utf-8
// c-file-style: "stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 tw=80 smartindent */
