module hdf5_RF

    use displayCarvalhol
    use math_RF
    use hdf5
    use mpi
    use write_Log_File

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_h5attr_bool(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        logical, intent(in) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        integer :: tmp_data

        dims(1) = 1
        if (data) then
            tmp_data = 1
        else
            tmp_data = 0
        end if
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(file_id, attr_name, H5T_STD_I8LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, tmp_data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_h5attr_bool

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_h5attr_int(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        integer, intent(in) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        integer :: tmp_data

        tmp_data = data

        dims(1) = 1
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(file_id, attr_name, H5T_STD_I32LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, tmp_data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_h5attr_int

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine write_h5attr_int_long(file_id, attr_name, data)
!        use HDF5
!        integer(HID_T), intent(in) :: file_id
!        character(len=*), intent(in) :: attr_name
!        integer(kind=8), intent(in) :: data
!        integer :: hdferr
!        integer(HID_T) :: attr_id, space_id
!        integer(HSIZE_T), dimension(1) :: dims
!
!        !NOT WORKING, h5awrite_f doesn't accept integer(kind=8)
!
!        dims(1) = 1
!        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
!        call h5acreate_f(file_id, attr_name, H5T_STD_I64LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
!        call h5awrite_f(attr_id, H5T_STD_I64LE, int(data,4), dims, hdferr)
!        call h5aclose_f(attr_id, hdferr)
!        call h5sclose_f(space_id, hdferr)
!    end subroutine write_h5attr_int_long

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_h5attr_real(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        double precision, intent(in) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = 1
        call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
        call h5acreate_f(file_id, attr_name, H5T_NATIVE_DOUBLE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_h5attr_real

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_h5attr_int_vec(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        integer, dimension(:), intent(in) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        integer, dimension(size(data)) :: tmp_data

        tmp_data = data
        dims(1) = size(data)
        call h5screate_f(H5S_SIMPLE_F, space_id, hdferr)
        call h5sset_extent_simple_f(space_id, 1, dims, dims, hdferr) !Set the size of the attribute
        call h5acreate_f(file_id, attr_name, H5T_STD_I32LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, tmp_data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine write_h5attr_int_vec

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine write_h5attr_int_long_vec(file_id, attr_name, data)
!        use HDF5
!        integer(HID_T), intent(in) :: file_id
!        character(len=*), intent(in) :: attr_name
!        integer(kind=8), dimension(:), intent(in) :: data
!        integer :: hdferr
!        integer(HID_T) :: attr_id, space_id
!        integer(HSIZE_T), dimension(1) :: dims
!
!        !NOT WORKING, h5awrite_f doesn't accept integer(kind=8)
!
!        dims(1) = size(data)
!        call h5screate_f(H5S_SIMPLE_F, space_id, hdferr)
!        call h5sset_extent_simple_f(space_id, 1, dims, dims, hdferr) !Set the size of the attribute
!        call h5acreate_f(file_id, attr_name, H5T_STD_I64LE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
!        call h5awrite_f(attr_id, H5T_STD_I64LE, int(data), dims, hdferr)
!        call h5aclose_f(attr_id, hdferr)
!        call h5sclose_f(space_id, hdferr)
!
!    end subroutine write_h5attr_int_long_vec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_h5attr_real_vec(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        double precision, dimension(:), intent(in) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = size(data)
        call h5screate_f(H5S_SIMPLE_F, space_id, hdferr)
        call h5sset_extent_simple_f(space_id, 1, dims, dims, hdferr) !Set the size of the attribute
        call h5acreate_f(file_id, attr_name, H5T_NATIVE_DOUBLE, space_id, attr_id, hdferr, H5P_DEFAULT_F)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine write_h5attr_real_vec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_h5attr_int(dset, attr, value)
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

    end subroutine read_h5attr_int

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_h5attr_bool(dset, attr, value)
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
    end subroutine read_h5attr_bool

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_h5attr_real(dset, attr, value)
        use HDF5
        integer(HID_T), intent(in) :: dset
        character(len=*), intent(in) :: attr
        double precision, intent(out) :: value
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
    end subroutine read_h5attr_real

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_h5attr_real_vec(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        double precision, dimension(:), intent(out) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims, max_dims

        dims(1) = size(data)
        call h5aopen_f(file_id, attr_name, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdferr) !Get the size of the attribute
        call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine read_h5attr_real_vec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_h5attr_int_vec(file_id, attr_name, data)
        use HDF5
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: attr_name
        integer, dimension(:), intent(out) :: data
        integer :: hdferr
        integer(HID_T) :: attr_id, space_id
        integer(HSIZE_T), dimension(1) :: dims, max_dims

        dims(1) = size(data)
        call h5aopen_f(file_id, attr_name, attr_id, hdferr)
        call h5aget_space_f(attr_id, space_id, hdferr)
        call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdferr) !Get the size of the attribute
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine read_h5attr_int_vec

end module hdf5_RF
