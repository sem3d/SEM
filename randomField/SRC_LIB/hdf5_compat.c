/* hdf5_compat.c
 * Weak fallback symbols for parallel HDF5 Fortran routines
 * when linking with serial HDF5 on developer machines.
 */

void h5pset_dxpl_mpio_f_(long long *prp_id, int *comm, int *hdferr) __attribute__((weak));
void h5pset_dxpl_mpio_f_(long long *prp_id, int *comm, int *hdferr) {
    (void)prp_id;
    (void)comm;
    if (hdferr) *hdferr = 0;
}

void h5pset_fapl_mpio_f_(long long *prp_id, int *comm, int *info, int *hdferr) __attribute__((weak));
void h5pset_fapl_mpio_f_(long long *prp_id, int *comm, int *info, int *hdferr) {
    (void)prp_id;
    (void)comm;
    (void)info;
    if (hdferr) *hdferr = 0;
}
