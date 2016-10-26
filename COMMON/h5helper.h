/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
//

#ifndef _H5HELPER_H_
#define _H5HELPER_H_

// hdf5.h may  include mpi.h which can include mpicxx.h and we don't want that since
// we would need to link against libmpicxx ...
#define OMPI_SKIP_MPICXX
#include <hdf5.h>
#include <vector>

void h5h_set_errhandler();
void h5h_restore_errhandler();

void  h5h_create_attr(hid_t parent, const char* name, int value);
void  h5h_create_attr(hid_t parent, const char* name, double value);

hid_t h5h_dset_prop(hsize_t d0, hsize_t d1);
hid_t h5h_dset_prop(hsize_t d0);

void  h5h_write_dset(hid_t parent, const char* name, int d0, const double* arr);
void  h5h_write_dset(hid_t parent, const char* name, const std::vector<double>& arr);
void  h5h_write_dset_2d(hid_t parent, const char* name, int d0, int d1, const double* arr);
void  h5h_write_dset_2d(hid_t parent, const char* dname, int d1, const std::vector<double>& v);

void  h5h_write_dset(hid_t parent, const char* name, int d0, const int* arr);
void  h5h_write_dset(hid_t parent, const char* dname, const std::vector<int>& v);

void  h5h_write_dset_2d(hid_t parent, const char* name, int d0, int d1, const int* arr);
void  h5h_write_dset_2d(hid_t parent, const char* dname, int d1, const std::vector<int>& v);



int   h5h_read_attr_int(hid_t dset_id, const char* attrname);
void  h5h_write_attr_int(hid_t dset_id, const char* attrname, int val);
int   h5h_get_dset1d_size(hid_t dset_id);
void  h5h_get_dset2d_size(hid_t dset_id, hsize_t& d1, hsize_t& d2);

void  h5h_read_dset_Nx2(hid_t g, const char* dname, std::vector<double>& v, std::vector<double>& w);
void  h5h_read_dset_Nx3(hid_t g, const char* dname,
                       std::vector<double>& u, std::vector<double>& v, std::vector<double>& w);
void  h5h_read_dset(hid_t g, const char* dname, std::vector<int>& v);
void  h5h_read_dset_2d(hid_t g, const char* dname, int& d0, int& d1, std::vector<double>& data);
void  h5h_read_dset_2d(hid_t g, const char* dname, int& d0, int& d1, std::vector<int>& data);



#undef _MIN
#undef _MAX

#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
