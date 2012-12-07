!>
!!\file sem_hdf5.f90
!!\brief Contient les définitions utiles pour la lecture/ecriture de fichiers hdf5
!!\author L. A.
!!\version 1.0
!!\date 2011.11
!!
!<

module sem_hdf5
    use HDF5

    logical :: h5_initialized = .false.

    integer(HID_T) :: h5_gzip_prop_1d
    integer(HID_T) :: h5_gzip_prop_2d

    integer, parameter :: chunk0 = 1024
    integer, parameter :: chunk1 = 1024

    integer(HSIZE_T), DIMENSION(1:2) :: chunk = (/chunk0,chunk1/)

contains

    !> Structure des fichiers HDF5 utilises :
    !!
    !! Protection/reprise: (Numero, champ)
    !! Elements/
    !! x Offets : Un tableau (Nelem x Nfield) d'offset contenant la position dans
    !!            le tableau de numero nfield des donnees de l'element Nelem
    !! 1 : Veloc
    !! 2 : Veloc1
    !! 2 : Veloc2
    !! 2 : Veloc3
    !! 3 : Displ
    !! 4 : EpsilonVol
    !! 5 : R_vol
    !! 6 : R_xx
    !! 6 : R_yy
    !! 6 : R_xy
    !! 6 : R_xz
    !! 6 : R_yz
    !! 7 : EpsilonDev_xx
    !! 7 : EpsilonDev_yy
    !! 7 : EpsilonDev_xy
    !! 7 : EpsilonDev_xz
    !! 7 : EpsilonDev_yz
    !! 8 : Stress
    subroutine init_hdf5()

        integer :: hdferr

        if (h5_initialized) then
            return
        end if

        CALL h5open_f(hdferr)
        CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5_gzip_prop_2d, hdferr)
        CALL h5pset_deflate_f(h5_gzip_prop_2d, 5, hdferr)
    end subroutine init_hdf5

    subroutine create_dset(parent, name, dtype, nmax, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer, INTENT(IN) :: nmax
        integer(HSIZE_T), dimension(1) :: dims, chunk
        integer(HID_T) :: space_id, prop_id
        integer :: hdferr

        dims(1) = nmax
        chunk(1) = min(nmax, 256*1024)
        call h5screate_simple_f(1, dims, space_id, hdferr, dims)
        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        if (nmax.gt.128) then
            call h5pset_deflate_f(prop_id, 5, hdferr)
            call h5pset_chunk_f(prop_id, 1, chunk, hdferr)
            call h5pset_shuffle_f(prop_id, hdferr)
        end if
        call h5dcreate_f(parent, name, dtype, space_id, dset_id, hdferr, prop_id)
        !write(*,*) "h5dcreate: ", name, hdferr, dims, chunk
        call h5pclose_f(prop_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine create_dset

    subroutine create_dset_2d(parent, name, dtype, d1, d2, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer, INTENT(IN) :: d1, d2
        integer(HSIZE_T), dimension(2) :: dims, chunk
        integer(HID_T) :: space_id, prop_id
        integer :: hdferr

        dims(1) = d1
        dims(2) = d2
        chunk(1) = min(d1, 256*1024)
        chunk(2) = max(1, min(d2, 256*1024/chunk(1)))
        call h5screate_simple_f(2, dims, space_id, hdferr, dims)
        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        if ((d1*d2).gt.128) then
            call h5pset_deflate_f(prop_id, 5, hdferr)
            call h5pset_chunk_f(prop_id, 2, chunk, hdferr)
        end if
        call h5dcreate_f(parent, name, dtype, space_id, dset_id, hdferr, prop_id)
        !write(*,*) "h5dcreate: ", name, hdferr, dims, chunk
        call h5pclose_f(prop_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine create_dset_2d

    subroutine read_dset_1d_real(parent, name, data)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        real*8, dimension(:), allocatable, intent(out) :: data
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(1) :: dims, maxdims
        integer :: hdferr, ndims

        call h5dopen_f(parent, name, dset_id, hdferr)
        !write(*,*) "Reading attr: ", name, dset_id, hdferr
        call h5dget_space_f(dset_id, space_id, hdferr)
        !write(*,*) "Get type: ", space_id, hdferr
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        !write(*,*) "Dims: ", dims, maxdims, hdferr
        allocate(data(1:dims(1)))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine read_dset_1d_real

    subroutine write_attr_int(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        integer*4, intent(in) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = 1
        !write(*,*) "save_attr_int: ", attr, value
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(dset, attr, H5T_STD_I32LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_STD_I32LE, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_attr_int

    subroutine write_attr_real(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        real*8, intent(in) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = 1
        !write(*,*) "save_attr_real: ", attr, value
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(dset, attr, H5T_NATIVE_DOUBLE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_attr_real

    subroutine read_attr_int(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        integer*4, intent(out) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        dims(1) = 1
        call h5aopen_f(dset, attr, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        ! TODO: check h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        call h5aread_f(attr_id, H5T_STD_I32LE, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine read_attr_int

    subroutine read_attr_real(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        real*8, intent(out) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        dims(1) = 1
        call h5aopen_f(dset, attr, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        ! TODO: check h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine read_attr_real

end module sem_hdf5
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
