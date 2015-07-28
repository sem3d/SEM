
module msnapshots
    use sdomain
    use hdf5
    use sem_hdf5
    use semdatafiles
    use mpi
    use mfields
    use sem_c_config, only : sem_mkdir
    use constants, only : M_1_3
    implicit none
contains

    subroutine grp_write_real_2d(Tdomain, parent_id, name, dim1, dim2, data, ntot_nodes)
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim1, dim2
        real, dimension(0:dim1-1,0:dim2-1), intent(in) :: data
        character(len=*), INTENT(IN) :: name
        integer, intent(out) :: ntot_nodes
        !
        integer(HID_T) :: dset_id
        real, dimension(:,:), allocatable :: all_data
        integer, dimension(:), allocatable :: displs, counts
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr

        if (Tdomain%output_rank==0) then
            allocate(displs(0:Tdomain%nb_output_procs-1))
            allocate(counts(0:Tdomain%nb_output_procs-1))
        end if
        call MPI_Gather(dim2, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, Tdomain%comm_output, ierr)

        ntot_nodes = 0
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = dim1*ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
                counts(n)=counts(n)*dim1 ! for the next gatherv
            end do
            allocate(all_data(0:dim1-1,0:ntot_nodes-1))
        end if

        call MPI_Gatherv(data, dim1*dim2, MPI_DOUBLE_PRECISION, all_data, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = dim1
            dims(2) = ntot_nodes
            call create_dset_2d(parent_id, name, H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
            deallocate(all_data)
            deallocate(displs)
            deallocate(counts)
        end if
    end subroutine grp_write_real_2d


    subroutine grp_write_fields(Tdomain, parent_id, dim2, displ, veloc, accel, press, &
        eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, eps_dev_xy, eps_dev_xz, eps_dev_yz, &
        sig_dev_xx, sig_dev_yy, sig_dev_zz, sig_dev_xy, sig_dev_xz, sig_dev_yz, P_energy, S_energy, ntot_nodes)

        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim2
        real, dimension(0:2,0:dim2-1), intent(in) :: displ, veloc, accel
        real, dimension(0:dim2-1), intent(in) :: press
        real, dimension(0:dim2-1), intent(in) :: eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, &
            eps_dev_xy, eps_dev_xz, eps_dev_yz
        real, dimension(0:dim2-1), intent(in) :: sig_dev_xx, sig_dev_yy, sig_dev_zz, sig_dev_xy, &
            sig_dev_xz, sig_dev_yz, P_energy, S_energy
        integer, intent(out) :: ntot_nodes
        !
        integer(HID_T) :: dset_id
        real, dimension(:,:), allocatable :: all_data_2d
        real, dimension(:), allocatable :: all_data_1d
        integer, dimension(:), allocatable :: displs, counts
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr


        if (Tdomain%output_rank==0) then
            allocate(displs(0:Tdomain%nb_output_procs-1))
            allocate(counts(0:Tdomain%nb_output_procs-1))
        end if
        call MPI_Gather(dim2, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, Tdomain%comm_output, ierr)

        ! 1D FIELDS
        ntot_nodes = 0
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
            end do
            allocate(all_data_1d(0:ntot_nodes-1))
        end if
        ! PRESSION
        call MPI_Gatherv(press, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "press", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! VOL_EPS
        call MPI_Gatherv(eps_vol, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "eps_vol", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! EPS_DEV_XX
        call MPI_Gatherv(eps_dev_xx, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "eps_dev_xx", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! EPS_DEV_YY
        call MPI_Gatherv(eps_dev_yy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "eps_dev_yy", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! EPS_DEV_ZZ
        call MPI_Gatherv(eps_dev_zz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "eps_dev_zz", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! EPS_DEV_XY
        call MPI_Gatherv(eps_dev_xy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "eps_dev_xy", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! EPS_DEV_XZ
        call MPI_Gatherv(eps_dev_xz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "eps_dev_xz", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! EPS_DEV_YZ
        call MPI_Gatherv(eps_dev_yz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "eps_dev_yz", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if

        ! SIG_DEV_XX
        call MPI_Gatherv(sig_dev_xx, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "sig_dev_xx", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! SIG_DEV_YY
        call MPI_Gatherv(sig_dev_yy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "sig_dev_yy", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! SIG_DEV_ZZ
        call MPI_Gatherv(sig_dev_zz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "sig_dev_zz", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! SIG_DEV_XY
        call MPI_Gatherv(sig_dev_xy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "sig_dev_xy", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! SIG_DEV_XZ
        call MPI_Gatherv(sig_dev_xz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "sig_dev_xz", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! SIG_DEV_YZ
        call MPI_Gatherv(sig_dev_yz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "sig_dev_yz", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! P_ENERGY
        call MPI_Gatherv(P_energy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "P_energy", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! S_ENERGY
        call MPI_Gatherv(S_energy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "S_energy", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! 3D FIELDS
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = displs(n)*3
                counts(n) = counts(n)*3
            end do
            allocate(all_data_2d(0:2,0:ntot_nodes-1))
        end if
        ! VELOCITY
        call MPI_Gatherv(veloc, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = 3
            dims(2) = ntot_nodes
            call create_dset_2d(parent_id, "veloc", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! DISPLACEMENT
        call MPI_Gatherv(displ, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = 3
            dims(2) = ntot_nodes
            call create_dset_2d(parent_id, "displ", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        ! ACCELERATION
        call MPI_Gatherv(accel, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = 3
            dims(2) = ntot_nodes
            call create_dset_2d(parent_id, "accel", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        if (Tdomain%output_rank==0) then
            deallocate(all_data_1d)
            deallocate(all_data_2d)
            deallocate(displs)
            deallocate(counts)
        end if

    end subroutine grp_write_fields


    subroutine grp_write_int_2d(Tdomain, parent_id, name, dim1, dim2, data, ntot_nodes)
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim1, dim2
        integer, dimension(:,:), intent(in) :: data
        character(len=*), INTENT(IN) :: name
        integer, intent(out) :: ntot_nodes
        !
        integer(HID_T) :: dset_id
        integer, dimension(:,:), allocatable :: all_data
        integer, dimension(:), allocatable :: displs, counts
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr

        if (Tdomain%output_rank==0) then
            allocate(displs(0:Tdomain%nb_output_procs-1))
            allocate(counts(0:Tdomain%nb_output_procs-1))
        end if
        call MPI_Gather(dim2, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, Tdomain%comm_output, ierr)

        ntot_nodes = 0
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = dim1*ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
                counts(n)=counts(n)*dim1 ! for the next gatherv
            end do
            allocate(all_data(dim1,0:ntot_nodes-1))
        end if

        call MPI_Gatherv(data, dim1*dim2, MPI_INTEGER, all_data, counts, displs, &
            MPI_INTEGER, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = dim1
            dims(2) = ntot_nodes
            call create_dset_2d(parent_id, name, H5T_STD_I32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, all_data, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
            deallocate(all_data)
            deallocate(displs)
            deallocate(counts)
        end if
    end subroutine grp_write_int_2d

    subroutine grp_write_int_1d(Tdomain, parent_id, name, dim1, data, ntot_nodes)
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim1
        integer, dimension(:), intent(in) :: data
        character(len=*), INTENT(IN) :: name
        integer, intent(OUT) :: ntot_nodes
        !
        integer(HID_T) :: dset_id
        integer, dimension(:), allocatable :: all_data
        integer, dimension(:), allocatable :: displs, counts
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr

        if (Tdomain%output_rank==0) then
            allocate(displs(0:Tdomain%nb_output_procs-1))
            allocate(counts(0:Tdomain%nb_output_procs-1))
        end if
        call MPI_Gather(dim1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, Tdomain%comm_output, ierr)

        ntot_nodes = 0
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
            end do
            allocate(all_data(0:ntot_nodes-1))
        end if

        call MPI_Gatherv(data, dim1, MPI_INTEGER, all_data, counts, displs, &
            MPI_INTEGER, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, name, H5T_STD_I32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, all_data, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
            deallocate(all_data)
            deallocate(displs)
            deallocate(counts)
        end if
    end subroutine grp_write_int_1d

    subroutine grp_write_real_1d(Tdomain, parent_id, name, dim1, data, ntot_nodes)
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim1
        real, dimension(:), intent(in) :: data
        character(len=*), INTENT(IN) :: name
        integer, intent(OUT) :: ntot_nodes
        !
        integer(HID_T) :: dset_id
        real, dimension(:), allocatable :: all_data
        integer, dimension(:), allocatable :: displs, counts
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr

        if (Tdomain%output_rank==0) then
            allocate(displs(0:Tdomain%nb_output_procs-1))
            allocate(counts(0:Tdomain%nb_output_procs-1))
        end if
        call MPI_Gather(dim1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, Tdomain%comm_output, ierr)

        ntot_nodes = 0
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
            end do
            allocate(all_data(0:ntot_nodes-1))
        end if

        call MPI_Gatherv(data, dim1, MPI_DOUBLE_PRECISION, all_data, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, name, H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
            deallocate(all_data)
            deallocate(displs)
            deallocate(counts)
        end if
    end subroutine grp_write_real_1d


    subroutine compute_saved_elements(Tdomain, irenum, nnodes)
        type (domain), intent (INOUT):: Tdomain
        integer, allocatable, dimension(:), intent(out) :: irenum ! maps Iglobnum to file node number
        integer, intent(out) :: nnodes
        integer :: n, i, j, k, ngllx, nglly, ngllz, ig, gn, ne
        !
        integer :: group, count
        integer :: ierr

        group = Tdomain%rank/Tdomain%ngroup

        allocate(irenum(0:Tdomain%n_glob_points-1))

        irenum = -1
        ig = 0
        ne = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            ne = ne + 1
            do k = 0,ngllz - 1
                do j = 0,nglly - 1
                    do i = 0,ngllx - 1
                        gn = Tdomain%specel(n)%Iglobnum(i,j,k)
                        if (irenum(gn) == -1) then
                            irenum(gn) = ig
                            ig = ig + 1
                        end if
                    end do
                end do
            end do
        end do
        nnodes = ig

        if (.not. allocated(Tdomain%output_nodes)) then
            allocate(Tdomain%output_nodes(0:Tdomain%nb_output_procs-1))
            allocate(Tdomain%output_nodes_offset(0:Tdomain%nb_output_procs-1))
        end if

        call MPI_Allgather(nnodes, 1, MPI_INTEGER, Tdomain%output_nodes, 1, MPI_INTEGER, Tdomain%comm_output, ierr)

        count = 0
        do i=0,Tdomain%nb_output_procs-1
            Tdomain%output_nodes_offset(i) = count
            count = count + Tdomain%output_nodes(i)
        end do

    end subroutine compute_saved_elements

    subroutine create_dir_sorties(Tdomain, isort)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer,intent(in) :: isort
        character(Len=MAX_FILE_SIZE) :: temp
        integer :: code, ierr, rg

        rg = Tdomain%rank
        if (rg==0) then
            call semname_snap_result_dir(isort, temp)
            ierr = sem_mkdir(trim(adjustl(temp)))
        end if
        call mpi_barrier(Tdomain%communicateur, code)
    end subroutine create_dir_sorties

    !>
    !! Ecrit la geometrie pour les sorties dans un fichier HDF5
    !! Le format est destine a etre relu par paraview / xdmf
    !<
    subroutine write_snapshot_geom(Tdomain)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer :: rg
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        integer :: hdferr, code, ierr
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes
        integer, dimension(0:Tdomain%nb_procs-1) :: nodes_per_proc
        integer :: group

        rg = Tdomain%rank
        group = rg/Tdomain%ngroup

        call init_hdf5()
        if (rg==0) then
            ierr = sem_mkdir(trim(adjustl(path_results)))
        end if
        call mpi_barrier(Tdomain%communicateur, code)
        call semname_snap_geom_file(group, fnamef)

        if (Tdomain%output_rank==0) then
            call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        else
            fid = -1
        endif

        call mpi_barrier(Tdomain%communicateur, code)

        call compute_saved_elements(Tdomain, irenum, nnodes)

        call write_global_nodes(Tdomain, fid, irenum, nnodes)

        call write_elem_connectivity(Tdomain, fid, irenum)

        call write_constant_fields(Tdomain, fid, irenum, nnodes)

        if (Tdomain%output_rank==0) then
            call h5fclose_f(fid, hdferr)
        endif

        call mpi_gather(nnodes, 1, MPI_INTEGER, nodes_per_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (rg==0) call write_master_xdmf(Tdomain, nodes_per_proc)
    end subroutine write_snapshot_geom


    subroutine write_global_nodes(Tdomain, fid, irenum, nnodes)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, intent(in) :: nnodes
        integer, dimension(:), intent(in), allocatable :: irenum
        !
        real, dimension(:,:), allocatable :: nodes
        integer :: n
        integer :: nnodes_tot
        !
        allocate(nodes(0:2,0:nnodes-1))
        do n = 0, Tdomain%n_glob_points-1
            if (irenum(n)>=0) then
                nodes(:,irenum(n)) = Tdomain%GlobCoord(:,n)
            end if
        end do

        call grp_write_real_2d(Tdomain, fid, "Nodes", 3, nnodes, nodes, nnodes_tot)
        deallocate(nodes)
    end subroutine write_global_nodes

    subroutine write_elem_connectivity(Tdomain, fid, irenum)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        !
        integer :: ngllx, nglly, ngllz
        integer(HSIZE_T), dimension(2) :: dims
        integer, dimension(:,:), allocatable :: data
        integer, dimension(:), allocatable :: mat, iglobnum, proc
        integer, dimension(3,0:Tdomain%n_elem-1) :: ngll
        integer :: count, ig, nglobnum
        integer :: i, j, k, n, nb_elem
        integer :: nb_elem_tot, nb_gll_tot, nglob_tot
        integer :: ioffset

        ! First we count the number of hexaedrons
        count = 0
        nglobnum = 0
        k = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll(1,k) = Tdomain%specel(n)%ngllx
            ngll(2,k) = Tdomain%specel(n)%nglly
            ngll(3,k) = Tdomain%specel(n)%ngllz
            ! Max number of global points (can count elem vertices twice)
            nglobnum = nglobnum + ngll(1,k)*ngll(2,k)*ngll(3,k)
            ! Number of subelements
            count = count+(ngll(1,k)-1)*(ngll(2,k)-1)*(ngll(3,k)-1)
            k = k + 1
        enddo
        nb_elem = k

        call grp_write_int_2d(Tdomain, fid, "NGLL", 3, nb_elem, ngll, nb_gll_tot)

        allocate( data(1:8,0:count-1))
        allocate( mat(0:count-1))
        allocate( proc(0:count-1))
        allocate( iglobnum(0:nglobnum-1))
        dims(1) = 8
        dims(2) = count
        count = 0
        ig = 0
        ioffset = Tdomain%output_nodes_offset(Tdomain%output_rank)
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-2
                do j = 0,nglly-2
                    do i = 0,ngllx-2
                        data(1,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+0,j+0,k+0))
                        data(2,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+1,j+0,k+0))
                        data(3,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+1,j+1,k+0))
                        data(4,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+0,j+1,k+0))
                        data(5,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+0,j+0,k+1))
                        data(6,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+1,j+0,k+1))
                        data(7,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+1,j+1,k+1))
                        data(8,count) = ioffset+irenum(Tdomain%specel(n)%Iglobnum(i+0,j+1,k+1))
                        mat(count) = Tdomain%specel(n)%mat_index
                        proc(count) = Tdomain%rank
                        count=count+1
                    end do
                end do
            end do
            do k = 0,ngllz - 1
                do j = 0,nglly - 1
                    do i = 0,ngllx - 1
                        iglobnum(ig) = ioffset + irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        ig = ig + 1
                    end do
                end do
            end do
        end do

        call grp_write_int_2d(Tdomain, fid, "Elements", 8, count, data, nb_elem_tot)
        Tdomain%n_hexa = nb_elem_tot
        call grp_write_int_1d(Tdomain, fid, "Iglobnum", nglobnum, iglobnum, nglob_tot)
        call grp_write_int_1d(Tdomain, fid, "Material", count, mat, nb_elem_tot)
        call grp_write_int_1d(Tdomain, fid, "Proc", count, proc, nb_elem_tot)

        deallocate(mat)
        deallocate(proc)
        deallocate(iglobnum)
        deallocate(data)
    end subroutine write_elem_connectivity

    subroutine save_field_h5(Tdomain, isort)

        use sdomain
        use forces_aniso
        use assembly
        implicit none

        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: isort
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        real, dimension(:,:),allocatable :: displ, veloc, accel
        real, dimension(:), allocatable :: press
        real, dimension(:,:,:,:),allocatable :: field_displ, field_veloc, field_accel
        real, dimension(:,:,:),allocatable :: field_press
        integer, dimension(:), allocatable :: valence
        integer :: hdferr
        integer :: ngllx, nglly, ngllz, idx, mat_idx
        integer :: i, j, k, n
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes, group, nnodes_tot

        ! start modif
        !        real, dimension (:,:), allocatable  :: htprimex, hprimey, hprimez
        !        real, dimension (:,:,:),allocatable :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
        real, dimension(:), allocatable     :: eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, &
            eps_dev_xy, eps_dev_xz, eps_dev_yz
        real, dimension(:), allocatable     :: sig_dev_xx, sig_dev_yy, sig_dev_zz, sig_dev_xy, &
            sig_dev_xz, sig_dev_yz
        real, dimension(:), allocatable     :: P_energy, S_energy
        ! end modif

        call create_dir_sorties(Tdomain, isort)

        call compute_saved_elements(Tdomain, irenum, nnodes)

        allocate(displ(0:2,0:nnodes-1))
        allocate(veloc(0:2,0:nnodes-1))
        allocate(accel(0:2,0:nnodes-1))
        allocate(press(0:nnodes-1))
        allocate(valence(0:nnodes-1))

        ngllx = 0
        nglly = 0
        ngllz = 0
        valence(:) = 0
        veloc(:,:) = 0
        accel(:,:) = 0
        !        eps_vol(:) = 0
        !        eps_dev_xx(:) = 0
        !        eps_dev_yy(:) = 0
        !        eps_dev_zz(:) = 0
        !        eps_dev_xy(:) = 0
        !        eps_dev_xz(:) = 0
        !        eps_dev_yz(:) = 0

        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            if (ngllx /= Tdomain%specel(n)%ngllx .or. &
                nglly /= Tdomain%specel(n)%nglly .or. &
                ngllz /= Tdomain%specel(n)%ngllz) then
                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz
                if (allocated(field_displ)) deallocate(field_displ)
                if (allocated(field_veloc)) deallocate(field_veloc)
                if (allocated(field_accel)) deallocate(field_accel)
                if (allocated(field_press)) deallocate(field_press)
                allocate(field_displ(0:ngllx-1,0:nglly-1,0:ngllz-1,3))
                allocate(field_veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,3))
                allocate(field_accel(0:ngllx-1,0:nglly-1,0:ngllz-1,3))
                allocate(field_press(0:ngllx-1,0:nglly-1,0:ngllz-1))

                ! start modif
            !                if (allocated(DXX)) deallocate(DXX)
            !                if (allocated(DXY)) deallocate(DXY)
            !                if (allocated(DXZ)) deallocate(DXZ)
            !                if (allocated(DYX)) deallocate(DYX)
            !                if (allocated(DYY)) deallocate(DYY)
            !                if (allocated(DYZ)) deallocate(DYZ)
            !                if (allocated(DZX)) deallocate(DZX)
            !                if (allocated(DZY)) deallocate(DZY)
            !                if (allocated(DZZ)) deallocate(DZZ)
            !                if (allocated(hTprimex)) deallocate(hTprimex)
            !                if (allocated(hprimey)) deallocate(hprimey)
            !                if (allocated(hprimez)) deallocate(hprimez)
            !                allocate(DXX(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DXY(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DXZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DYX(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DYY(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DYZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DZX(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DZY(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(DZZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
            !                allocate(hTprimex(0:ngllx-1,0:ngllx-1))
            !                allocate(hprimey(0:nglly-1,0:nglly-1))
            !                allocate(hprimez(0:ngllz-1,0:ngllz-1))
                ! end modif
            endif
            call gather_elem_displ(Tdomain, n, field_displ)
            call gather_elem_veloc(Tdomain, n, field_veloc)
            call gather_elem_accel(Tdomain, n, field_accel)
            call gather_elem_press(Tdomain, n, field_press)

            ! start modif
            mat_idx=Tdomain%specel(n)%mat_index
            !            hTprimex=Tdomain%sSubDomain(mat_idx)%hTprimex
            !            hprimey=Tdomain%sSubDomain(mat_idx)%hprimey
            !            hprimez=Tdomain%sSubDomain(mat_idx)%hprimez

            !            if(Tdomain%specel(n)%solid .and. (.not. Tdomain%specel(n)%PML))then   ! SOLID PART OF THE DOMAIN
            !                call physical_part_deriv(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,Tdomain%specel(n)%InvGrad,field_displ(:,:,:,0),DXX,DYX,DZX)
            !                call physical_part_deriv(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,Tdomain%specel(n)%InvGrad,field_displ(:,:,:,1),DXY,DYY,DZY)
            !                call physical_part_deriv(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,Tdomain%specel(n)%InvGrad,field_displ(:,:,:,2),DXZ,DYZ,DZZ)
            !            endif

            ! end modif

            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        valence(idx) = valence(idx)+1
                        displ(:,idx) = field_displ(i,j,k,:)
                        veloc(:,idx) = veloc(:,idx)+field_veloc(i,j,k,:)
                        accel(:,idx) = accel(:,idx)+field_accel(i,j,k,:)
                        press(idx) = field_press(i,j,k)
                    !                        ! start modif
                    !                        if (Tdomain%specel(n)%solid .and. (.not. Tdomain%specel(n)%PML)) then
                    !                            eps_vol(idx) = DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k)
                    !                            eps_dev_xx(idx) = DXX(i,j,k) - eps_vol(idx) * 0.333333333333333333333333333333d0
                    !                            eps_dev_yy(idx) = DYY(i,j,k) - eps_vol(idx) * 0.333333333333333333333333333333d0
                    !                            eps_dev_zz(idx) = DZZ(i,j,k) - eps_vol(idx) * 0.333333333333333333333333333333d0
                    !                            eps_dev_xy(idx) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                    !                            eps_dev_xz(idx) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                    !                            eps_dev_yz(idx) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                    !                        else
                    !                            eps_vol(idx) = 0
                    !                            eps_dev_xx(idx) = 0
                    !                            eps_dev_yy(idx) = 0
                    !                            eps_dev_zz(idx) = 0
                    !                            eps_dev_xy(idx) = 0
                    !                            eps_dev_xz(idx) = 0
                    !                            eps_dev_yz(idx) = 0
                    !                        endif
                    !                         ! end modif
                    end do
                end do
            end do
        end do

        ! start modif
        allocate(eps_vol(0:nnodes-1))
        allocate(eps_dev_xx(0:nnodes-1))
        allocate(eps_dev_yy(0:nnodes-1))
        allocate(eps_dev_zz(0:nnodes-1))
        allocate(eps_dev_xy(0:nnodes-1))
        allocate(eps_dev_xz(0:nnodes-1))
        allocate(eps_dev_yz(0:nnodes-1))

        allocate(sig_dev_xx(0:nnodes-1))
        allocate(sig_dev_yy(0:nnodes-1))
        allocate(sig_dev_zz(0:nnodes-1))
        allocate(sig_dev_xy(0:nnodes-1))
        allocate(sig_dev_xz(0:nnodes-1))
        allocate(sig_dev_yz(0:nnodes-1))

        allocate(P_energy(0:nnodes-1))
        allocate(S_energy(0:nnodes-1))

        call compute_stress_strain(Tdomain, irenum, eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, &
            eps_dev_xy, eps_dev_xz, eps_dev_yz, sig_dev_xx, sig_dev_yy, sig_dev_zz, &
            sig_dev_xy, sig_dev_xz, sig_dev_yz, P_energy, S_energy)
        ! end modif

        ! normalization
        do i = 0,nnodes-1
            if (valence(i)/=0) then
                veloc(0:2,i) = veloc(0:2,i)/valence(i)
                accel(0:2,i) = accel(0:2,i)/valence(i)
            else
                write(*,*) "Elem",i," non traite"
            end if
        end do

        if (Tdomain%output_rank==0) then
            group = Tdomain%rank/Tdomain%ngroup
            call semname_snap_result_file(group, isort, fnamef)
            call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        else
            fid = -1
        endif
        call grp_write_fields(Tdomain, fid, nnodes, displ, veloc, accel, press, &
            eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, eps_dev_xy, eps_dev_xz, eps_dev_yz, &
             sig_dev_xx, sig_dev_yy, sig_dev_zz, sig_dev_xy, sig_dev_xz, sig_dev_yz, P_energy, &
             S_energy, nnodes_tot)

        if (Tdomain%output_rank==0) then
            call h5fclose_f(fid, hdferr)
            call write_xdmf(Tdomain, group, isort, nnodes_tot)
        endif
        deallocate(displ,veloc,accel,press,valence, &
            eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, eps_dev_xy, eps_dev_xz, eps_dev_yz, &
            sig_dev_xx, sig_dev_yy, sig_dev_zz, sig_dev_xy, sig_dev_xz, sig_dev_yz, P_energy, S_energy)

        if (allocated(field_displ)) deallocate(field_displ)
        if (allocated(field_veloc)) deallocate(field_veloc)
        if (allocated(field_accel)) deallocate(field_accel)
        if (allocated(field_press)) deallocate(field_press)
        !        if (allocated(DXX)) deallocate(DXX)
        !        if (allocated(DXY)) deallocate(DXY)
        !        if (allocated(DXZ)) deallocate(DXZ)
        !        if (allocated(DYX)) deallocate(DYX)
        !        if (allocated(DYY)) deallocate(DYY)
        !        if (allocated(DYZ)) deallocate(DYZ)
        !        if (allocated(DZX)) deallocate(DZX)
        !        if (allocated(DZY)) deallocate(DZY)
        !        if (allocated(DZZ)) deallocate(DZZ)
        !        if (allocated(hTprimex)) deallocate(hTprimex)
        !        if (allocated(hprimey)) deallocate(hprimey)
        !        if (allocated(hprimez)) deallocate(hprimez)

        call mpi_barrier(Tdomain%communicateur, hdferr)

    end subroutine save_field_h5

    subroutine write_master_xdmf(Tdomain, nodes_per_proc)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: n_procs, nelem, n_groups
        character (len=MAX_FILE_SIZE) :: fnamef
        integer, dimension(0:Tdomain%nb_procs-1), intent(in) :: nodes_per_proc
        integer :: group
        n_procs = Tdomain%nb_procs
        n_groups = (n_procs+Tdomain%ngroup-1)/Tdomain%ngroup
        nelem = Tdomain%n_elem
        call semname_xdmf_master(fnamef)

        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write(61,"(a)") '<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
        write(61,"(a)") '<Domain>'
        write(61,"(a)") '<Grid CollectionType="Spatial" GridType="Collection">'
        !!! XXX: recuperer le nom par semname_*
        do group=0,n_groups-1
            write(61,"(a,I4.4,a)") '<xi:include href="mesh.',group,'.xmf"/>'
        end do

        !do rg=0,n_procs-1
        !    if (nodes_per_proc(rg).gt.0) then
        !        write(61,"(a,I4.4,a)") '<xi:include href="mesh.',rg,'.xmf"/>'
        !    endif
        !end do
        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
        close(61)

    end subroutine write_master_xdmf

    subroutine write_xdmf(Tdomain, group, isort, nnodes)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer, intent(in) :: group, isort, nnodes
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: i, nn, ne
        real :: time
        call semname_xdmf(group, fnamef)

        nn = nnodes
        ne = Tdomain%n_hexa
        if (ne==0) then
            return
        end if
        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a,I4.4,a)") '<Grid CollectionType="Temporal" GridType="Collection" Name="space.',group,'">'
        write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mat" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',group,'.h5:/Material</DataItem>'
        write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Proc" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',group,'.h5:/Proc</DataItem>'

        write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',group,'.h5:/Mass</DataItem>'
        write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Jac" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',group,'.h5:/Jac</DataItem>'
        time = 0
        do i=1,isort
            write(61,"(a,I4.4,a,I4.4,a)") '<Grid Name="mesh.',i,'.',group,'">'
            write(61,"(a,F20.10,a)") '<Time Value="', time,'"/>'
            write(61,"(a,I9,a)") '<Topology Type="Hexahedron" NumberOfElements="',ne,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Int" Dimensions="',ne,' 8">'
            write(61,"(a,I4.4,a)") 'geometry',group,'.h5:/Elements'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Topology>'
            write(61,"(a)") '<Geometry Type="XYZ">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a)") 'geometry',group,'.h5:/Nodes'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Geometry>'
            ! DISPL
            write(61,"(a,I9,a)") '<Attribute Name="Displ" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/displ'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! VELOC
            write(61,"(a,I9,a)") '<Attribute Name="Veloc" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/veloc'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! ACCEL
            write(61,"(a,I9,a)") '<Attribute Name="Accel" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/accel'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! PRESSURE
            write(61,"(a,I9,a)") '<Attribute Name="Pressure" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/press'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! VOLUMETRIC STRAIN
            write(61,"(a,I9,a)") '<Attribute Name="eps_vol" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_vol'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! start modifs
            ! EPS_DEV_XX
            write(61,"(a,I9,a)") '<Attribute Name="eps_dev_xx" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_xx'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! EPS_DEV_XX
            write(61,"(a,I9,a)") '<Attribute Name="eps_dev_yy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_yy'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! EPS_DEV_ZZ
            write(61,"(a,I9,a)") '<Attribute Name="eps_dev_zz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_zz'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! EPS_DEV_XY
            write(61,"(a,I9,a)") '<Attribute Name="eps_dev_xy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_xy'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! EPS_DEV_XZ
            write(61,"(a,I9,a)") '<Attribute Name="eps_dev_xz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_xz'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! EPS_DEV_YZ
            write(61,"(a,I9,a)") '<Attribute Name="eps_dev_yz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_yz'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! SIG_DEV_XX
            write(61,"(a,I9,a)") '<Attribute Name="sig_dev_xx" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_xx'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! SIG_DEV_XX
            write(61,"(a,I9,a)") '<Attribute Name="sig_dev_yy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_yy'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! SIG_DEV_ZZ
            write(61,"(a,I9,a)") '<Attribute Name="sig_dev_zz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_zz'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! SIG_DEV_XY
            write(61,"(a,I9,a)") '<Attribute Name="sig_dev_xy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_xy'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! SIG_DEV_XZ
            write(61,"(a,I9,a)") '<Attribute Name="sig_dev_xz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_xz'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! SIG_DEV_YZ
            write(61,"(a,I9,a)") '<Attribute Name="sig_dev_yz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_yz'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! P_ENERGY
            write(61,"(a,I9,a)") '<Attribute Name="P_energy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/P_energy'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! S_ENERGY
            write(61,"(a,I9,a)") '<Attribute Name="S_energy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/S_energy'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! end modifs

            ! DOMAIN
            write(61,"(a)") '<Attribute Name="Domain" Center="Grid" AttributeType="Scalar" Dimensions="1">'
            write(61,"(a,I4,a)") '<DataItem Format="XML" Datatype="Int"  Dimensions="1">',group,'</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! MAT
            write(61,"(a,I9,a)") '<Attribute Name="Mat" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',group, &
                '"]/DataItem[@Name="Mat"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! PROC
            write(61,"(a,I9,a)") '<Attribute Name="Proc" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',group, &
                '"]/DataItem[@Name="Proc"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! MASS
            write(61,"(a,I9,a)") '<Attribute Name="Mass" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',group, &
                '"]/DataItem[@Name="Mass"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I9,a)") '<Attribute Name="Jac" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',group, &
                '"]/DataItem[@Name="Jac"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '</Grid>'
            ! XXX inexact pour l'instant
            time = time+Tdomain%TimeD%time_snapshots
        end do
        write(61,"(a)") '</Grid>'
        close(61)
    end subroutine write_xdmf

    subroutine write_constant_fields(Tdomain, fid, irenum, nnodes)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        integer, intent(in) :: nnodes
        !
        real, dimension(:),allocatable :: mass, jac
        integer :: ngllx, nglly, ngllz, idx
        integer :: i, j, k, n, nnodes_tot


        allocate(mass(0:nnodes-1))
        allocate(jac(0:nnodes-1))
        ! mass
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        mass(idx) = Tdomain%specel(n)%MassMat(i,j,k)
                    end do
                end do
            end do
        end do
        ! jac
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        jac(idx) = Tdomain%specel(n)%Jacob(i,j,k)
                    end do
                end do
            end do
        end do


        do n = 0,Tdomain%n_face-1
            ngllx = Tdomain%sface(n)%ngll1
            nglly = Tdomain%sface(n)%ngll2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    idx = irenum(Tdomain%sface(n)%Iglobnum_Face(i,j))
                    if (idx>=0) mass(idx) = Tdomain%sface(n)%MassMat(i,j)
                end do
            end do
        end do

        do n = 0,Tdomain%n_edge-1
            ngllx = Tdomain%sedge(n)%ngll
            do i = 1,ngllx-2
                idx = irenum(Tdomain%sedge(n)%Iglobnum_Edge(i))
                if (idx>=0) mass(idx) = Tdomain%sedge(n)%MassMat(i)
            end do
        end do

        do n = 0,Tdomain%n_vertex-1
            idx = irenum(Tdomain%svertex(n)%Iglobnum_Vertex)
            if (idx>=0) mass(idx) = Tdomain%svertex(n)%MassMat
        end do

        call grp_write_real_1d(Tdomain, fid, "Mass", nnodes, mass, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Jac", nnodes, jac, nnodes_tot)
        deallocate(mass,jac)

    end subroutine write_constant_fields

    subroutine compute_stress_strain(Tdomain, irenum, eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, &
        eps_dev_xy, eps_dev_xz, eps_dev_yz, sig_dev_xx, sig_dev_yy, sig_dev_zz, &
        sig_dev_xy, sig_dev_xz, sig_dev_yz, P_energy, S_energy)

        type(domain), intent(in) :: Tdomain

        integer :: ngllx, nglly, ngllz, idx, mat_idx, nnodes, i, j, k, n, n_solid

        real, dimension(:,:,:,:), allocatable :: field_displ
        real, dimension(:,:,:), allocatable :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
        real, dimension(:,:), allocatable :: hTprimex, hprimey, hprimez
        real, dimension(:), allocatable, intent(INOUT) :: eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, &
            eps_dev_xy, eps_dev_xz, eps_dev_yz
        real, dimension(:), allocatable, intent(INOUT) :: sig_dev_xx, sig_dev_yy, sig_dev_zz, &
            sig_dev_xy, sig_dev_xz, sig_dev_yz
        real, dimension(:), allocatable, intent(INOUT) :: P_energy, S_energy
        integer, dimension(:), allocatable, intent(IN) :: irenum
        logical :: aniso, solid
        real :: xmu, xlambda, xkappa, x2mu, xlambda2mu, onemSbeta, onemPbeta
        real :: eps_trace

        ngllx = 0
        nglly = 0
        ngllz = 0

        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            if (ngllx /= Tdomain%specel(n)%ngllx .or. &
                nglly /= Tdomain%specel(n)%nglly .or. &
                ngllz /= Tdomain%specel(n)%ngllz) then
                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz

                if (allocated(field_displ)) deallocate(field_displ)

                allocate(field_displ(0:ngllx-1,0:nglly-1,0:ngllz-1,3))

                if (allocated(DXX)) deallocate(DXX)
                if (allocated(DXY)) deallocate(DXY)
                if (allocated(DXZ)) deallocate(DXZ)
                if (allocated(DYX)) deallocate(DYX)
                if (allocated(DYY)) deallocate(DYY)
                if (allocated(DYZ)) deallocate(DYZ)
                if (allocated(DZX)) deallocate(DZX)
                if (allocated(DZY)) deallocate(DZY)
                if (allocated(DZZ)) deallocate(DZZ)
                if (allocated(hTprimex)) deallocate(hTprimex)
                if (allocated(hprimey)) deallocate(hprimey)
                if (allocated(hprimez)) deallocate(hprimez)
                allocate(DXX(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DXY(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DXZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DYX(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DYY(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DYZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DZX(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DZY(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(DZZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(hTprimex(0:ngllx-1,0:ngllx-1))
                allocate(hprimey(0:nglly-1,0:nglly-1))
                allocate(hprimez(0:ngllz-1,0:ngllz-1))

            endif

            call gather_elem_displ(Tdomain, n, field_displ)

            solid=Tdomain%specel(n)%solid
            n_solid=Tdomain%n_sls
            mat_idx=Tdomain%specel(n)%mat_index
            hTprimex=Tdomain%sSubDomain(mat_idx)%hTprimex
            hprimey=Tdomain%sSubDomain(mat_idx)%hprimey
            hprimez=Tdomain%sSubDomain(mat_idx)%hprimez
            aniso=Tdomain%aniso

            ! displacement gradient
            if((solid) .and. (.not. Tdomain%specel(n)%PML))then
                call physical_part_deriv(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,Tdomain%specel(n)%InvGrad,field_displ(:,:,:,0),DXX,DYX,DZX)
                call physical_part_deriv(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,Tdomain%specel(n)%InvGrad,field_displ(:,:,:,1),DXY,DYY,DZY)
                call physical_part_deriv(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,Tdomain%specel(n)%InvGrad,field_displ(:,:,:,2),DXZ,DYZ,DZZ)
            endif

            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1

                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))

                        if ((solid) .and. (.not. Tdomain%specel(n)%PML)) then
                            ! strain tensor
                            eps_trace       = DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k)
                            eps_vol(idx)    = eps_trace
                            eps_dev_xx(idx) = DXX(i,j,k) - eps_trace / 3
                            eps_dev_yy(idx) = DYY(i,j,k) - eps_trace / 3
                            eps_dev_zz(idx) = DZZ(i,j,k) - eps_trace / 3
                            eps_dev_xy(idx) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                            eps_dev_xz(idx) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                            eps_dev_yz(idx) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                            ! stress tensor
                            ! TODO stress tensor for anisotropic and attenuation
                            if (aniso) then
                            else
                                xmu     = Tdomain%specel(n)%Mu(i,j,k)
                                xlambda = Tdomain%specel(n)%Lambda(i,j,k)
                                xkappa  = Tdomain%specel(n)%Kappa(i,j,k)

                                if (n_solid>0) then
                                    onemSbeta=Tdomain%specel(n)%sl%onemSbeta(i,j,k)
                                    onemPbeta=Tdomain%specel(n)%sl%onemPbeta(i,j,k)
                                    !  mu_relaxed -> mu_unrelaxed
                                    xmu    = xmu * onemSbeta
                                    !  kappa_relaxed -> kappa_unrelaxed
                                    xkappa = xkappa * onemPbeta
                                endif
                                x2mu       = 2. * xmu
                                xlambda2mu = xlambda + x2mu

                                sig_dev_xx(idx) = xlambda2mu * DXX(i,j,k) + xlambda * (DYY(i,j,k) + DZZ(i,j,k))
                                sig_dev_yy(idx) = xlambda2mu * DYY(i,j,k) + xlambda * (DXX(i,j,k) + DZZ(i,j,k))
                                sig_dev_zz(idx) = xlambda2mu * DZZ(i,j,k) + xlambda * (DXX(i,j,k) + DYY(i,j,k))
                                sig_dev_xy(idx) = xmu * (DXY(i,j,k) + DYX(i,j,k))
                                sig_dev_xz(idx) = xmu * (DXZ(i,j,k) + DZX(i,j,k))
                                sig_dev_yz(idx) = xmu * (DYZ(i,j,k) + DZY(i,j,k))

                                P_energy(idx) = .5 * xlambda2mu * eps_trace**2
                                S_energy(idx) = xmu * ((DXX(i,j,k) - eps_trace/ 3)**2 + &
                                                       (DYY(i,j,k) - eps_trace/ 3)**2 + &
                                                       (DZZ(i,j,k) - eps_trace/ 3)**2 + &
                                                        DXY(i,j,k)**2 + DYX(i,j,k)**2 + &
                                                        DXZ(i,j,k)**2 + DZX(i,j,k)**2 + &
                                                        DYZ(i,j,k)**2 + DZY(i,j,k)**2)


                            endif
                        else
                            eps_vol(idx) = 0
                            eps_dev_xx(idx) = 0
                            eps_dev_yy(idx) = 0
                            eps_dev_zz(idx) = 0
                            eps_dev_xy(idx) = 0
                            eps_dev_xz(idx) = 0
                            eps_dev_yz(idx) = 0
                            sig_dev_xx(idx) = 0
                            sig_dev_yy(idx) = 0
                            sig_dev_zz(idx) = 0
                            sig_dev_xy(idx) = 0
                            sig_dev_xz(idx) = 0
                            sig_dev_yz(idx) = 0
                            P_energy(idx)   = 0
                            S_energy(idx)   = 0

                        endif
                    end do
                end do
            end do
        end do

        if (allocated(field_displ)) deallocate(field_displ)
        if (allocated(hTprimex)) deallocate(hTprimex)
        if (allocated(hprimey)) deallocate(hprimey)
        if (allocated(hprimez)) deallocate(hprimez)
        if (allocated(DXX)) deallocate(DXX)
        if (allocated(DXY)) deallocate(DXY)
        if (allocated(DXZ)) deallocate(DXZ)
        if (allocated(DYX)) deallocate(DYX)
        if (allocated(DYY)) deallocate(DYY)
        if (allocated(DYZ)) deallocate(DYZ)
        if (allocated(DZX)) deallocate(DZX)
        if (allocated(DZY)) deallocate(DZY)
        if (allocated(DZZ)) deallocate(DZZ)

    end subroutine compute_stress_strain

end module msnapshots
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
