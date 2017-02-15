!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
    use constants, only : fpp
    implicit none

    logical :: h5_initialized = .false.

    integer(HID_T) :: h5_gzip_prop_1d
    integer(HID_T) :: h5_gzip_prop_2d

    integer, parameter :: chunk0 = 1024
    integer, parameter :: chunk1 = 1024

    integer(HSIZE_T), DIMENSION(1:2) :: chunk = (/chunk0,chunk1/)

#ifdef SINGLEPRECISION
#define H5T_REAL  H5T_NATIVE_REAL
#else
#define H5T_REAL  H5T_NATIVE_DOUBLE
#endif
    interface create_dset
       module procedure create_dset_i4, create_dset_i8
    end interface create_dset

    interface create_dset_2d
       module procedure create_dset_2d_i4, create_dset_2d_i8
    end interface create_dset_2d

    interface write_dataset
        module procedure write_dataset_d1, write_dataset_d2, write_dataset_i1, write_dataset_i2, &
            write_dataset_d3
    end interface write_dataset

    interface read_dataset
        module procedure read_dset_1d_int, read_dset_2d_int, read_dset_1d_real, read_dset_2d_real, &
            read_dset_3d_real
    end interface read_dataset

    interface append_dataset_2d
       module procedure append_dataset_2d_r, append_dataset_2d_i
    end interface append_dataset_2d
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

    subroutine create_dset_i8(parent, name, dtype, nmax, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer(HSIZE_T), INTENT(IN) :: nmax
        integer(HSIZE_T), dimension(1) :: dims, chunk, maxdims
        integer(HID_T) :: space_id, prop_id
        integer :: hdferr

        dims(1) = nmax
        maxdims(1) = nmax
        if (nmax==H5S_UNLIMITED_F) then
            chunk(1) = 256
            dims(1) = 0
        else
            chunk(1) = min(nmax, 256*1024_HSIZE_T)
        end if
        call h5screate_simple_f(1, dims, space_id, hdferr, dims)
        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        if (nmax.gt.128 .or. nmax==H5S_UNLIMITED_F) then
            call h5pset_deflate_f(prop_id, 5, hdferr)
            call h5pset_chunk_f(prop_id, 1, chunk, hdferr)
            call h5pset_shuffle_f(prop_id, hdferr)
        end if
        call h5dcreate_f(parent, name, dtype, space_id, dset_id, hdferr, prop_id)
        !write(*,*) "h5dcreate: ", name, hdferr, dims, chunk
        call h5pclose_f(prop_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine create_dset_i8

    subroutine create_dset_i4(parent, name, dtype, nmax, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer, INTENT(IN) :: nmax

        call create_dset_i8(parent, name, dtype, int(nmax,HSIZE_T), dset_id)
    end subroutine create_dset_i4

    subroutine create_dset_2d_i8(parent, name, dtype, d1, d2, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer(HSIZE_T), INTENT(IN) :: d1, d2
        integer(HSIZE_T), dimension(2) :: dims, chunk, maxdims
        integer(HID_T) :: space_id, prop_id
        integer :: hdferr

        dims(1) = d1
        dims(2) = d2
        maxdims(1) = d1
        maxdims(2) = d2
        chunk(1) = min(d1, 256*1024_HSIZE_T)
        if (d2==H5S_UNLIMITED_F) then
            chunk(2) = 64
            dims(2) = 0
        else
            chunk(2) = max(1_HSIZE_T, min(d2, int(256*1024/chunk(1),HSIZE_T)))
        endif
        call h5screate_simple_f(2, dims, space_id, hdferr, maxdims)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5screate KO"
        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pcreate KO"
        if ((d1*d2).gt.128 .or. d2==H5S_UNLIMITED_F) then
            call h5pset_chunk_f(prop_id, 2, chunk, hdferr)
            if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_chunk KO"
            if (dtype/=H5T_IEEE_F32LE .and. dtype/=H5T_IEEE_F64LE) then
                ! Les donnees en float donnent un taux de compression bas pour un
                ! cout de calcul eleve
                call h5pset_deflate_f(prop_id, 5, hdferr)
                if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_deflate KO"
                call h5pset_shuffle_f(prop_id, hdferr)
                if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_shuffle KO"
            endif
        end if
        call h5dcreate_f(parent, name, dtype, space_id, dset_id, hdferr, prop_id)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5dcreate KO"
        call h5pclose_f(prop_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5sclose KO"
    end subroutine create_dset_2d_i8

    subroutine create_dset_3d_i8(parent, name, dtype, dims, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer(HSIZE_T), INTENT(INOUT), dimension(3) :: dims
        integer(HSIZE_T), dimension(3) :: chunk, maxdims
        integer(HID_T) :: space_id, prop_id
        integer :: hdferr

        maxdims = dims
        chunk(1) = min(dims(1), 256_HSIZE_T)
        chunk(2) = min(dims(2), 1024_HSIZE_T)
        if (dims(3)==H5S_UNLIMITED_F) then
            chunk(3) = 64
            dims(3) = 0
        else
            chunk(3) = max(1_HSIZE_T, min(dims(3), int(256*1024/(chunk(1)*chunk(2)),HSIZE_T)))
        endif
        call h5screate_simple_f(3, dims, space_id, hdferr, maxdims)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5screate KO"
        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pcreate KO"
        if ((dims(1)*dims(2)*dims(3)).gt.128 .or. dims(3)==H5S_UNLIMITED_F) then
            call h5pset_chunk_f(prop_id, 3, chunk, hdferr)
            if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_chunk KO"
            if (dtype/=H5T_IEEE_F32LE .and. dtype/=H5T_IEEE_F64LE) then
                ! Les donnees en float donnent un taux de compression bas pour un
                ! cout de calcul eleve
                call h5pset_deflate_f(prop_id, 5, hdferr)
                if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_deflate KO"
                call h5pset_shuffle_f(prop_id, hdferr)
                if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pset_shuffle KO"
            endif
        end if
        call h5dcreate_f(parent, name, dtype, space_id, dset_id, hdferr, prop_id)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5dcreate KO"
        call h5pclose_f(prop_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5pclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "create_dset_2d_i8 : h5sclose KO"
    end subroutine create_dset_3d_i8

    subroutine create_dset_2d_i4(parent, name, dtype, d1, d2, dset_id)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent, dtype
        integer(HID_T), INTENT(OUT) :: dset_id
        integer, INTENT(IN) :: d1, d2

        call create_dset_2d_i8(parent, name, dtype, int(d1,HSIZE_T), int(d2,HSIZE_T), dset_id)
    end subroutine create_dset_2d_i4

    subroutine read_dset_1d_real(parent, name, data, ibase)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        real(fpp), dimension(:), allocatable, intent(out) :: data
        integer, intent(in), optional :: ibase
        !
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(1) :: dims, maxdims
        integer :: hdferr, i0
        i0 = 1
        if (present(ibase)) then
          i0 = ibase
        end if

        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_real : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_real : h5dgetspace KO"
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr .ne. 1) stop "read_dset_1d_real : h5sgetdim KO"
        allocate(data(i0:dims(1)+i0-1))
        call h5dread_f(dset_id, H5T_REAL, data, dims, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_real : h5dread KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_real : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_real : h5sclose KO"
    end subroutine read_dset_1d_real

    subroutine read_dset_1d_int(parent, name, data, ibase)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        integer, dimension(:), allocatable, intent(out) :: data
        integer, intent(in), optional :: ibase
        !
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(1) :: dims, maxdims
        integer :: hdferr, i0
        i0 = 1
        if (present(ibase)) then
          i0 = ibase
        end if

        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_int : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_int : h5dgetspace KO"
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr .ne. 1) stop "read_dset_1d_int : h5sgetdim KO"
        allocate(data(i0:dims(1)+i0-1))
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_int : h5dread KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_int : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_1d_int : h5sclose KO"
    end subroutine read_dset_1d_int

    subroutine read_dims(parent, name, dims)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        integer(HSIZE_T), dimension(:), allocatable, intent(out) :: dims
        !
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), allocatable, dimension(:) :: maxdims
        integer :: hdferr, ndims
        !
        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dims_3d : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dims_3d : h5dgetspace KO"
        call h5sget_simple_extent_ndims_f(space_id, ndims, hdferr)
        if (ndims<0) stop "read_dims_3d: error"
        allocate(dims(ndims), maxdims(ndims))
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine read_dims

    subroutine read_dset_2d_real(parent, name, data, ibase)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        real(fpp), dimension(:,:), allocatable, intent(out) :: data
        integer, intent(in), optional :: ibase
        !
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(2) :: dims, maxdims
        integer :: hdferr, i0, i1
        i0 = 1
        i1 = 1
        if (present(ibase)) then
          i0 = ibase
          i1 = ibase
        end if

        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_real : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_real : h5dgetspace KO"
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr .ne. 2) stop "read_dset_2d_real : h5sgetdim KO "
        allocate(data(i0:dims(1)+i0-1, i1:dims(2)+i1-1))
        call h5dread_f(dset_id, H5T_REAL, data, dims, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_real : h5dread KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_real : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_real : h5sclose KO"
    end subroutine read_dset_2d_real

    subroutine read_dset_3d_real(parent, name, data, ibase)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        real(fpp), dimension(:,:,:), allocatable, intent(out) :: data
        integer, intent(in), optional :: ibase
        !
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(3) :: dims, maxdims
        integer :: hdferr, i0, i1,i2
        i0 = 1
        i1 = 1
        i2 = 1
        if (present(ibase)) then
          i0 = ibase
          i1 = ibase
          i2 = ibase
        end if

        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5dgetspace KO"
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr .ne. 3) stop "read_dset_3d_real : h5sgetdim KO "
        allocate(data(i0:dims(1)+i0-1, i1:dims(2)+i1-1, i2:dims(3)+i2-1))
        call h5dread_f(dset_id, H5T_REAL, data, dims, hdferr)
        if (hdferr .ne. 0) then
            write(*,*) "read_dset_3d_real : h5dread KO", dims
            stop 1
        endif
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_3d_real : h5sclose KO"
    end subroutine read_dset_3d_real

    subroutine read_subset_3d_real(parent, name, imin, imax, data)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        integer, intent(in), dimension(0:2) :: imin, imax
        real(fpp), dimension(:,:,:), allocatable, intent(out) :: data
        !
        integer(HID_T) :: dset_id, space_id, memspace_id
        integer(HSIZE_T), dimension(3) :: start, count
        integer :: hdferr

        start = imin
        count = imax-imin+[1,1,1]
        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5dgetspace KO"
        call H5Sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, start, count, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5sgetdim KO "
        allocate(data(imin(0):imax(0), imin(1):imax(1), imin(2):imax(2)))
        call H5Screate_simple_f(3, count, memspace_id, hdferr)
        call h5dread_f(dset_id, H5T_REAL, data, count, hdferr, memspace_id, space_id)
        if (hdferr .ne. 0) then
            write(*,*) "read_subset_3d_real : h5dread KO", imin, imax
            stop 1
        end if
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5sclose KO"
        call h5sclose_f(memspace_id, hdferr)
        if (hdferr .ne. 0) stop "read_subset_3d_real : h5sclose KO"
    end subroutine read_subset_3d_real

    subroutine read_dset_2d_int(parent, name, data, ibase)
        use HDF5
        character(len=*), INTENT(IN) :: name
        integer(HID_T), INTENT(IN) :: parent
        integer, dimension(:,:), allocatable, intent(out) :: data
        integer, intent(in), optional :: ibase
        !
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(2) :: dims, maxdims
        integer :: hdferr, i0, i1
        i0 = 1
        i1 = 1
        if (present(ibase)) then
          i0 = ibase
          i1 = ibase
        end if

        call h5dopen_f(parent, name, dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_int : h5dopen KO"
        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_int : h5dgetspace KO"
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr .ne. 2) stop "read_dset_2d_int : h5sgetdim KO"
        allocate(data(i0:dims(1)+i0-1, i1:dims(2)+i1-1))
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_int : h5dread KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_int : h5dclose KO"
        call h5sclose_f(space_id, hdferr)
        if (hdferr .ne. 0) stop "read_dset_2d_int : h5sclose KO"
    end subroutine read_dset_2d_int

!<<<<<<< HEAD
!    subroutine write_attr_int(dset, attr, value)
!        use HDF5
!        integer(HID_T), intent(in) :: dset
!        character(len=*), intent(in) :: attr
!        integer, intent(in) :: value
!        integer :: hdferr
!        integer(HID_T) :: attr_id, space_id
!        integer(HSIZE_T), dimension(1) :: dims
!
!        dims(1) = 1
!        !write(*,*) "save_attr_int: ", attr, value
!        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
!        call h5acreate_f(dset, attr, H5T_STD_I32LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
!        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)
!        call h5aclose_f(attr_id, hdferr)
!        call h5sclose_f(space_id, hdferr)
!    end subroutine write_attr_int
     subroutine write_attr_int(dset, attr, value)
         use HDF5
         integer(HID_T), intent(in) :: dset
         character(len=*), intent(in) :: attr
         integer, intent(in) :: value
         integer :: hdferr
         integer(HID_T) :: attr_id, space_id
         integer(HSIZE_T), dimension(1) :: dims
         integer :: value_tmp
         value_tmp = value
         
         dims(1) = 1
         !write(*,*) "save_attr_int: ", attr, value
         call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
         call h5acreate_f(dset, attr, H5T_STD_I32LE, space_id, attr_id,hdferr, H5P_DEFAULT_F)
         call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value_tmp, dims,hdferr)
         call h5aclose_f(attr_id, hdferr)
         call h5sclose_f(space_id, hdferr)
     end subroutine write_attr_int
!=======
!    subroutine write_attr_int(dset, attr, value)
!        use HDF5
!        integer(HID_T), intent(in) :: dset
!        character(len=*), intent(in) :: attr
!        integer, intent(in) :: value
!        integer :: hdferr
!        integer(HID_T) :: attr_id, space_id
!        integer(HSIZE_T), dimension(1) :: dims
!        integer :: value_tmp
!
!        value_tmp = value
!        dims(1) = 1
!        !write(*,*) "save_attr_int: ", attr, value
!        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
!        call h5acreate_f(dset, attr, H5T_STD_I32LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
!        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value_tmp, dims, hdferr)
!        call h5aclose_f(attr_id, hdferr)
!        call h5sclose_f(space_id, hdferr)
!    end subroutine write_attr_int
!>>>>>>> 9a6b45d5e8874b8bcbb9288394b6f44f17e9d77e

    subroutine write_attr_real(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        real(fpp), intent(in) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = 1
        !write(*,*) "save_attr_real: ", attr, value
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(dset, attr, H5T_NATIVE_DOUBLE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_REAL, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_attr_real

    subroutine write_attr_bool(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        logical, intent(in) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        integer :: tmp_val

        dims(1) = 1
        !write(*,*) "save_attr_int: ", attr, value
        if (value) then
            tmp_val = 1
        else
            tmp_val = 0
        end if
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(dset, attr, H5T_STD_I8LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, tmp_val, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_attr_bool

    subroutine read_attr_int(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        integer, intent(out) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        dims(1) = 1
        call h5aopen_f(dset, attr, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        ! TODO: check h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine read_attr_int

    subroutine read_attr_bool(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        logical, intent(out) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        integer :: ival
        dims(1) = 1
        call h5aopen_f(dset, attr, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        ! TODO: check h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, ival, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
        if (ival==0) then
            value = .false.
        else
            value = .true.
        end if
    end subroutine read_attr_bool

    subroutine read_attr_real(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        real(fpp), intent(out) :: value
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        dims(1) = 1
        call h5aopen_f(dset, attr, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        ! TODO: check h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        call h5aread_f(attr_id, H5T_REAL, value, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine read_attr_real

    subroutine write_dataset_d1(parent, name, arr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: parent
        character(len=*), intent(in) :: name
        real(fpp), dimension(:), intent(in) :: arr
        !
        integer(HSIZE_T), dimension(1) ::  dims
        integer(HID_T) :: dset_id
        integer :: hdferr
        dims(1) = size(arr,1)
        call create_dset(parent, name, H5T_IEEE_F64LE, dims(1), dset_id)
        call h5dwrite_f(dset_id, H5T_REAL, arr, dims, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d1 : h5dwrite KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d1 : h5dclose KO"
    end subroutine write_dataset_d1

    subroutine write_dataset_d2(parent, name, arr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: parent
        character(len=*), intent(in) :: name
        real(fpp), dimension(:,:), intent(in) :: arr
        !
        integer(HSIZE_T), dimension(2) ::  dims
        integer(HID_T) :: dset_id
        integer :: hdferr
        dims(1) = size(arr,1)
        dims(2) = size(arr,2)
        call create_dset_2d_i8(parent, name, H5T_IEEE_F64LE, dims(1), dims(2), dset_id)
        call h5dwrite_f(dset_id, H5T_REAL, arr, dims, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d2 : h5dwrite KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d2 : h5dclose KO"
    end subroutine write_dataset_d2

    subroutine write_dataset_d3(parent, name, arr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: parent
        character(len=*), intent(in) :: name
        real(fpp), dimension(:,:,:), intent(in) :: arr
        !
        integer(HSIZE_T), dimension(3) ::  dims
        integer(HID_T) :: dset_id
        integer :: hdferr
        dims(1) = size(arr,1)
        dims(2) = size(arr,2)
        dims(3) = size(arr,3)
        call create_dset_3d_i8(parent, name, H5T_IEEE_F64LE, dims, dset_id)
        call h5dwrite_f(dset_id, H5T_REAL, arr, dims, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d3 : h5dwrite KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_d3 : h5dclose KO"
    end subroutine write_dataset_d3

    subroutine write_dataset_i1(parent, name, arr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: parent
        character(len=*), intent(in) :: name
        integer, dimension(:), intent(in) :: arr
        !
        integer(HSIZE_T), dimension(1) ::  dims
        integer(HID_T) :: dset_id
        integer :: hdferr
        dims(1) = size(arr,1)
        call create_dset(parent, name, H5T_STD_I32LE, dims(1), dset_id)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, arr, dims, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_i1 : h5dwrite KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_i1 : h5dclose KO"
    end subroutine write_dataset_i1

    subroutine write_dataset_i2(parent, name, arr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: parent
        character(len=*), intent(in) :: name
        integer, dimension(:,:), intent(in) :: arr
        !
        integer(HSIZE_T), dimension(2) ::  dims
        integer(HID_T) :: dset_id
        integer :: hdferr
        dims(1) = size(arr,1)
        dims(2) = size(arr,2)
        call create_dset_2d_i8(parent, name, H5T_STD_I32LE, dims(1), dims(2), dset_id)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, arr, dims, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_i2 : h5dwrite KO"
        call h5dclose_f(dset_id, hdferr)
        if (hdferr .ne. 0) stop "write_dataset_i2 : h5dclose KO"
    end subroutine write_dataset_i2

    subroutine append_dataset_2d_r(dset_id, arr, hdferr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: dset_id
        real(fpp), dimension(:,:), intent(in) :: arr
        integer, intent(out) :: hdferr
        !
        integer(HSIZE_T), dimension(2) ::  dims, maxdims, offset, dsize
         integer(HID_T) :: memspace, filespace
        !
        dims(1) = size(arr, 1)
        dims(2) = size(arr, 2)
        call H5Screate_simple_f(2, dims, memspace, hdferr)
        call H5Dget_space_f(dset_id, filespace, hdferr)
        call H5Sget_simple_extent_dims_f(filespace, dims, maxdims, hdferr)
        !write(*,*) "RFile dims:", dims, maxdims
        call H5Sclose_f(filespace, hdferr)

        dsize(1) = size(arr,1)
        dsize(2) = dims(2) + size(arr, 2)
        call H5Dextend_f(dset_id, dsize, hdferr)

        call H5Dget_space_f(dset_id, filespace, hdferr)
        offset(1) = 0
        offset(2) = dims(2)
        dims(1) = size(arr, 1)
        dims(2) = size(arr, 2)
        !write(*,*) "RW: offset:", offset
        !write(*,*) "RW: dims:", dims
        call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, hdferr)

        call H5Dwrite_f(dset_id, H5T_REAL, arr, dims, hdferr, memspace, filespace)
        call H5Sclose_f(filespace, hdferr)
        call H5Sclose_f(memspace, hdferr)
    end subroutine append_dataset_2d_r

    subroutine append_dataset_2d_i(dset_id, arr, hdferr)
        use HDF5
        implicit none
        integer(HID_T), intent(in) :: dset_id
        integer, dimension(:,:), intent(in) :: arr
        integer, intent(out) :: hdferr
        !
        integer(HSIZE_T), dimension(2) ::  dims, maxdims, offset, dsize
        integer(HID_T) :: memspace, filespace
        !
        dims(1) = size(arr, 1)
        dims(2) = size(arr, 2)
        call H5Screate_simple_f(2, dims, memspace, hdferr)
        call H5Dget_space_f(dset_id, filespace, hdferr)
        call H5Sget_simple_extent_dims_f(filespace, dims, maxdims, hdferr)
        call H5Sclose_f(filespace, hdferr)

        dsize(1) = size(arr,1)
        dsize(2) = dims(2) + size(arr, 2)
        call H5Dextend_f(dset_id, dsize, hdferr)

        call H5Dget_space_f(dset_id, filespace, hdferr)
        offset(1) = 0
        offset(2) = dims(2)
        dims(1) = size(arr, 1)
        dims(2) = size(arr, 2)
        call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, hdferr)

        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, arr, dims, hdferr, memspace, filespace)
        call H5Sclose_f(filespace, hdferr)
        call H5Sclose_f(memspace, hdferr)
    end subroutine append_dataset_2d_i

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_attr_real_vec(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        real(fpp), dimension(:), intent(out) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims, max_dims

        dims(1) = size(data)
        call h5aopen_f(file_id, attr_name, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdferr) !Get the size of the attribute
        call h5aread_f(attr_id, H5T_REAL, data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine read_attr_real_vec



end module sem_hdf5

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
