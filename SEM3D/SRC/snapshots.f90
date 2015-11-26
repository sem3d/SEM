!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module msnapshots
    use sdomain
    use hdf5
    use sem_hdf5
    use semdatafiles
    use mpi
    use mfields
    use deriv3d
    use sem_c_config, only : sem_mkdir
    use constants
    implicit none

    type :: output_var_t
        real, dimension(:,:), allocatable :: displ, veloc, accel
        real, dimension(:), allocatable :: press
        real, dimension(:), allocatable :: eps_vol, eps_dev_xx, eps_dev_yy, eps_dev_zz, &
            eps_dev_xy, eps_dev_xz, eps_dev_yz
        real, dimension(:), allocatable :: sig_dev_xx, sig_dev_yy, sig_dev_zz, sig_dev_xy, &
            sig_dev_xz, sig_dev_yz, P_energy, S_energy
        real, dimension(:), allocatable :: eps_dev_pl_xx, eps_dev_pl_yy, eps_dev_pl_zz, &
            eps_dev_pl_xy, eps_dev_pl_xz, eps_dev_pl_yz
    end type output_var_t

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

    subroutine grp_write_fields(Tdomain, parent_id, dim2, out_variables, outputs, ntot_nodes)

        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim2
        type(output_var_t), intent(in) :: outputs
        integer, dimension(0:8), intent(in) :: out_variables
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

        ! P_ENERGY
        if (out_variables(OUT_ENERGYP) == 1) then
            call MPI_Gatherv(outputs%P_energy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "P_energy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! S_ENERGY
        if (out_variables(OUT_ENERGYS) == 1) then
            call MPI_Gatherv(outputs%S_energy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "S_energy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! VOL_EPS
        if (out_variables(OUT_EPS_VOL) == 1) then
            call MPI_Gatherv(outputs%eps_vol, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "eps_vol", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! PRESSION
        if (out_variables(OUT_PRESSION) == 1) then
            call MPI_Gatherv(outputs%press, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "press", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if

        ! EPS_DEV
        if (out_variables(OUT_EPS_DEV) == 1) then
            ! EPS_DEV_XX
            call MPI_Gatherv(outputs%eps_dev_xx, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "eps_dev_xx", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_YY
            call MPI_Gatherv(outputs%eps_dev_yy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "eps_dev_yy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_ZZ
            call MPI_Gatherv(outputs%eps_dev_zz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "eps_dev_zz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_XY
            call MPI_Gatherv(outputs%eps_dev_xy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "eps_dev_xy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_XZ
            call MPI_Gatherv(outputs%eps_dev_xz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "eps_dev_xz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_YZ
            call MPI_Gatherv(outputs%eps_dev_yz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "eps_dev_yz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if

            if (Tdomain%nl_flag==1) then
                ! EPS_DEV_PL_XX
                call MPI_Gatherv(outputs%eps_dev_pl_xx, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_nodes
                    call create_dset(parent_id, "eps_dev_pl_xx", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_YY
                call MPI_Gatherv(outputs%eps_dev_pl_yy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_nodes
                    call create_dset(parent_id, "eps_dev_pl_yy", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_ZZ
                call MPI_Gatherv(outputs%eps_dev_pl_zz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_nodes
                    call create_dset(parent_id, "eps_dev_pl_zz", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_XY
                call MPI_Gatherv(outputs%eps_dev_pl_xy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_nodes
                    call create_dset(parent_id, "eps_dev_pl_xy", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_XZ
                call MPI_Gatherv(outputs%eps_dev_pl_xz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_nodes
                    call create_dset(parent_id, "eps_dev_pl_xz", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_YZ
                call MPI_Gatherv(outputs%eps_dev_pl_yz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_nodes
                    call create_dset(parent_id, "eps_dev_pl_yz", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
            end if
        end if

        ! SIG_DEV
        if (out_variables(OUT_STRESS_DEV) == 1) then
            ! SIG_DEV_XX
            call MPI_Gatherv(outputs%sig_dev_xx, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "sig_dev_xx", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_YY
            call MPI_Gatherv(outputs%sig_dev_yy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "sig_dev_yy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_ZZ
            call MPI_Gatherv(outputs%sig_dev_zz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "sig_dev_zz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_XY
            call MPI_Gatherv(outputs%sig_dev_xy, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "sig_dev_xy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_XZ
            call MPI_Gatherv(outputs%sig_dev_xz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "sig_dev_xz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_YZ
            call MPI_Gatherv(outputs%sig_dev_yz, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "sig_dev_yz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
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
        if (out_variables(OUT_VITESSE) == 1) then
            call MPI_Gatherv(outputs%veloc, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = 3
                dims(2) = ntot_nodes
                call create_dset_2d(parent_id, "veloc", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! DISPLACEMENT
        if (out_variables(OUT_DEPLA) == 1) then
            call MPI_Gatherv(outputs%displ, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = 3
                dims(2) = ntot_nodes
                call create_dset_2d(parent_id, "displ", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if

        ! ACCELERATION
        if (out_variables(OUT_ACCEL) == 1) then
            call MPI_Gatherv(outputs%accel, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = 3
                dims(2) = ntot_nodes
                call create_dset_2d(parent_id, "accel", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
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

    subroutine compute_saved_elements(Tdomain, irenum, nnodes, domains)
        type (domain), intent (INOUT):: Tdomain
        integer, allocatable, dimension(:), intent(out) :: irenum ! maps Iglobnum to file node number
        integer, intent(out) :: nnodes
        integer :: n, i, j, k, ngllx, nglly, ngllz, ig, gn, ne
        !
        integer :: count
        integer :: ierr, domain_type, imat
        integer, dimension(:), allocatable, intent(out) :: domains

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

        allocate(domains(0:ig-1))
        domain_type = 0
        domains = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            imat = Tdomain%specel(n)%mat_index
            domain_type = get_domain(Tdomain%sSubDomain(imat))

            do k = 0,ngllz - 1
                do j = 0,nglly - 1
                    do i = 0,ngllx - 1
                        gn = Tdomain%specel(n)%Iglobnum(i,j,k)
                        if (domain_type>domains(irenum(gn))) then
                            domains(irenum(gn)) = domain_type
                        endif
                    end do
                end do
            end do
        end do

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
        integer, dimension(:), allocatable :: domains

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

        call compute_saved_elements(Tdomain, irenum, nnodes, domains)

        call write_global_nodes(Tdomain, fid, irenum, nnodes)

        call write_elem_connectivity(Tdomain, fid, irenum)

        call write_constant_fields(Tdomain, fid, irenum, nnodes, domains)

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

    subroutine allocate_fields(nnodes, out_flags, fields, nl_flag)

        type (output_var_t), intent(inout) :: fields
        integer, dimension(0:8), intent(in) :: out_flags
        integer, intent(in) :: nnodes, nl_flag
        integer             :: flag_gradU
        flag_gradU =  out_flags(OUT_ENERGYP)+out_flags(OUT_ENERGYS)+&
            out_flags(OUT_EPS_VOL)+out_flags(OUT_EPS_DEV)+&
            out_flags(OUT_STRESS_DEV)

        ! ALLOCATE FIELDS
        if (out_flags(OUT_ENERGYP   ) == 1) allocate(fields%P_energy(0:nnodes-1))
        if (out_flags(OUT_ENERGYS   ) == 1) allocate(fields%S_energy(0:nnodes-1))
        if (out_flags(OUT_EPS_VOL   ) == 1) allocate(fields%eps_vol(0:nnodes-1))
        if (out_flags(OUT_PRESSION  ) == 1) allocate(fields%press(0:nnodes-1))
        if (out_flags(OUT_DEPLA     ) == 1) allocate(fields%displ(0:2,0:nnodes-1))
        if (out_flags(OUT_VITESSE   ) == 1) allocate(fields%veloc(0:2,0:nnodes-1))
        if (out_flags(OUT_ACCEL     ) == 1) allocate(fields%accel(0:2,0:nnodes-1))
        if (flag_gradU/=0) then
            allocate(fields%eps_dev_xx(0:nnodes-1))
            allocate(fields%eps_dev_yy(0:nnodes-1))
            allocate(fields%eps_dev_zz(0:nnodes-1))
            allocate(fields%eps_dev_xy(0:nnodes-1))
            allocate(fields%eps_dev_xz(0:nnodes-1))
            allocate(fields%eps_dev_yz(0:nnodes-1))
            if (nl_flag==1) then
                allocate(fields%eps_dev_pl_xx(0:nnodes-1))
                allocate(fields%eps_dev_pl_yy(0:nnodes-1))
                allocate(fields%eps_dev_pl_zz(0:nnodes-1))
                allocate(fields%eps_dev_pl_xy(0:nnodes-1))
                allocate(fields%eps_dev_pl_xz(0:nnodes-1))
                allocate(fields%eps_dev_pl_yz(0:nnodes-1))
            end if
            if (out_flags(OUT_STRESS_DEV) == 1) then
                allocate(fields%sig_dev_xx(0:nnodes-1))
                allocate(fields%sig_dev_yy(0:nnodes-1))
                allocate(fields%sig_dev_zz(0:nnodes-1))
                allocate(fields%sig_dev_xy(0:nnodes-1))
                allocate(fields%sig_dev_xz(0:nnodes-1))
                allocate(fields%sig_dev_yz(0:nnodes-1))
            end if
        end if
        ! INITIALIZE FIELDS
        if (out_flags(OUT_ENERGYP ) == 1) fields%P_energy = 0.
        if (out_flags(OUT_ENERGYS ) == 1) fields%S_energy = 0.
        if (out_flags(OUT_EPS_VOL ) == 1) fields%eps_vol = 0.
        if (out_flags(OUT_PRESSION) == 1) fields%press = 0.
        if (out_flags(OUT_DEPLA   ) == 1) fields%displ = 0.
        if (out_flags(OUT_VITESSE ) == 1) fields%veloc = 0.
        if (out_flags(OUT_ACCEL   ) == 1) fields%accel = 0.

        if (flag_gradU/=0) then
            fields%eps_dev_xx = 0.
            fields%eps_dev_yy = 0.
            fields%eps_dev_zz = 0.
            fields%eps_dev_xy = 0.
            fields%eps_dev_xz = 0.
            fields%eps_dev_yz = 0.
            if (nl_flag==1) then
                fields%eps_dev_pl_xx = 0.
                fields%eps_dev_pl_yy = 0.
                fields%eps_dev_pl_zz = 0.
                fields%eps_dev_pl_xy = 0.
                fields%eps_dev_pl_xz = 0.
                fields%eps_dev_pl_yz = 0.
            endif
            if (out_flags(OUT_STRESS_DEV) == 1) then
                fields%sig_dev_xx = 0.
                fields%sig_dev_yy = 0.
                fields%sig_dev_zz = 0.
                fields%sig_dev_xy = 0.
                fields%sig_dev_xz = 0.
                fields%sig_dev_yz = 0.
            endif
        end if

    end subroutine allocate_fields

    subroutine allocate_second_fields(domain_type, out_flags, nl_flag, ngllx, nglly, ngllz, &
        field_displ, field_veloc, field_accel, field_press, field_stress, field_eps_pl, field_phi, field_vphi, &
        DXX, DYY, DZZ, DXY, DYX, DXZ, DZX, DYZ, DZY)

        integer, intent(in) :: domain_type, nl_flag, ngllx, nglly, ngllz
        integer, dimension(0:8), intent(in) :: out_flags
        real, dimension(:,:,:,:), allocatable, intent(inout) :: field_displ, field_veloc, field_accel
        real, dimension(:,:,:,:), allocatable, intent(inout) :: field_stress, field_eps_pl
        real, dimension(:,:,:),   allocatable, intent(inout) :: field_press
        real, dimension(:,:,:),   allocatable, intent(inout) :: field_phi, field_vphi
        real, dimension(:,:,:),   allocatable, intent(inout) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
        integer                                              :: flag_gradU

        flag_gradU =  out_flags(OUT_ENERGYP)+out_flags(OUT_ENERGYS)+&
            out_flags(OUT_EPS_VOL)+out_flags(OUT_EPS_DEV)+&
            out_flags(OUT_STRESS_DEV)
        if (allocated(field_veloc))  deallocate(field_veloc)
        if (allocated(field_accel))  deallocate(field_accel)
        if (allocated(field_press))  deallocate(field_press)
        if (allocated(field_displ))  deallocate(field_displ)
        if (allocated(field_stress)) deallocate(field_stress)
        if (allocated(field_eps_pl)) deallocate(field_eps_pl)
        if (allocated(field_phi))    deallocate(field_phi)
        if (allocated(field_vphi))   deallocate(field_vphi)
        if (allocated(DXX)) then
            deallocate(DXX)
            deallocate(DXY)
            deallocate(DXZ)
            deallocate(DYX)
            deallocate(DYY)
            deallocate(DYZ)
            deallocate(DZX)
            deallocate(DZY)
            deallocate(DZZ)
        end if

        ! ALLOCATION
        select case(domain_type)
            case (DM_SOLID,DM_SOLID_PML) !  SOLID PART OF THE DOMAIN
                if (out_flags(OUT_DEPLA)==1 .or. flag_gradU/=0) allocate(field_displ(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                if (out_flags(OUT_VITESSE)==1) allocate(field_veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                if (out_flags(OUT_ACCEL)==1) allocate(field_accel(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                if (out_flags(OUT_PRESSION)==1) allocate(field_press(0:ngllx-1,0:nglly-1,0:ngllz-1))

                if (flag_gradU /=0) then
                    allocate(DXX(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DXY(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DXZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DYX(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DYY(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DYZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DZX(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DZY(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(DZZ(0:ngllx-1,0:nglly-1,0:ngllz-1))
                endif
                if (nl_flag == 1) then
                    if (flag_gradU /= 0 .or. out_flags(OUT_PRESSION) == 1) then
                        allocate(field_stress(0:ngllx-1,0:nglly-1,0:ngllz-1,0:5))
                    end if
                    if (flag_gradU /=0) then
                        allocate(field_eps_pl(0:ngllx-1,0:nglly-1,0:ngllz-1,0:5))
                    end if
                end if

            case (DM_FLUID,DM_FLUID_PML) ! FLUID PART OF THE DOMAIN
                allocate(field_phi(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(field_vphi(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(field_veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                allocate(field_accel(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                allocate(field_press(0:ngllx-1,0:nglly-1,0:ngllz-1))
        end select

        ! INITIALIZATION
        if (domain_type == DM_SOLID_PML) then
            if (out_flags(OUT_DEPLA   ) == 1) field_displ = 0
            if (out_flags(OUT_ACCEL   ) == 1) field_accel = 0
            if (out_flags(OUT_PRESSION) == 1) field_press = 0

            if (flag_gradU /=0) then
                DXX = 0
                DYY = 0
                DZZ = 0
                DXY = 0
                DYX = 0
                DXZ = 0
                DZX = 0
                DYZ = 0
                DZY = 0
            else if (flag_gradU == 0 .and. nl_flag == 0) then
                field_stress = 0
                field_eps_pl = 0
            end if
        end if

    end subroutine allocate_second_fields

    subroutine deallocate_fields(out_flags, fields, nl_flag)
        integer, dimension(0:8), intent(in) :: out_flags
        type (output_var_t), intent(inout) :: fields
        integer, intent(in) :: nl_flag
        
        if (out_flags(OUT_ENERGYP   ) == 1) deallocate(fields%P_energy)
        if (out_flags(OUT_ENERGYS   ) == 1) deallocate(fields%S_energy)
        if (out_flags(OUT_EPS_VOL   ) == 1) deallocate(fields%eps_vol)
        if (out_flags(OUT_PRESSION  ) == 1) deallocate(fields%press)
        if (out_flags(OUT_DEPLA     ) == 1) deallocate(fields%displ)
        if (out_flags(OUT_VITESSE   ) == 1) deallocate(fields%veloc)
        if (out_flags(OUT_ACCEL     ) == 1) deallocate(fields%accel)
        if (allocated(fields%eps_dev_xx)) then
            deallocate(fields%eps_dev_xx)
            deallocate(fields%eps_dev_yy)
            deallocate(fields%eps_dev_zz)
            deallocate(fields%eps_dev_xy)
            deallocate(fields%eps_dev_xz)
            deallocate(fields%eps_dev_yz)
        end if
        if (allocated(fields%eps_dev_pl_xx)) then
            deallocate(fields%eps_dev_pl_xx)
            deallocate(fields%eps_dev_pl_yy)
            deallocate(fields%eps_dev_pl_zz)
            deallocate(fields%eps_dev_pl_xy)
            deallocate(fields%eps_dev_pl_xz)
            deallocate(fields%eps_dev_pl_yz)
        end if
        
        if (allocated(fields%sig_dev_xx)) then
            deallocate(fields%sig_dev_xx)
            deallocate(fields%sig_dev_yy)
            deallocate(fields%sig_dev_zz)
            deallocate(fields%sig_dev_xy)
            deallocate(fields%sig_dev_xz)
            deallocate(fields%sig_dev_yz)
        end if
    end subroutine deallocate_fields

    subroutine save_field_h5(Tdomain, isort)
        use sdomain
        use forces_aniso

        implicit none

        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: isort
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        real, dimension(:,:,:,:), allocatable :: field_displ, field_veloc, field_accel
        real, dimension(:,:,:,:), allocatable :: field_eps_pl, field_stress
        real, dimension(:,:,:),   allocatable :: field_press
        real, dimension(:,:,:),   allocatable :: field_phi, field_vphi

        integer :: domain_type
        integer, dimension(:), allocatable :: valence
        integer :: hdferr
        integer :: ngllx, nglly, ngllz, idx
        integer :: i, j, k, n
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes, group, nnodes_tot
        integer, dimension(:), allocatable :: domains
        type(Element), pointer :: el
        type(subdomain), pointer :: sub_dom_mat
        type(output_var_t) :: out_fields
        integer :: n_solid

        real, dimension(:,:,:), allocatable   :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
        real    :: xmu, xlambda, xkappa, x2mu, xlambda2mu, onemSbeta, onemPbeta

        integer, dimension(0:8) :: out_variables
        integer                 :: flag_gradU


        out_variables(0:8) = Tdomain%out_variables(0:8)
        flag_gradU = out_variables(OUT_ENERGYP)+out_variables(OUT_ENERGYS)+&
            out_variables(OUT_EPS_VOL)+out_variables(OUT_EPS_DEV)+&
            out_variables(OUT_STRESS_DEV)
        n_solid = Tdomain%n_sls
        call create_dir_sorties(Tdomain, isort)

        call compute_saved_elements(Tdomain, irenum, nnodes, domains)

        call allocate_fields(nnodes, out_variables, out_fields, Tdomain%nl_flag)
        
        allocate(valence(0:nnodes-1))

        valence(:) = 0

        ngllx = 0
        nglly = 0
        ngllz = 0

        do n = 0,Tdomain%n_elem-1
            el => Tdomain%specel(n)
            sub_dom_mat => Tdomain%sSubdomain(el%mat_index)
            domain_type = get_domain(sub_dom_mat)
            if (.not. el%OUTPUT) cycle
            if (ngllx /= el%ngllx .or. nglly /= el%nglly .or. ngllz /= el%ngllz) then
                ngllx = el%ngllx
                nglly = el%nglly
                ngllz = el%ngllz
                call allocate_second_fields(domain_type, out_variables, Tdomain%nl_flag, ngllx, nglly, ngllz, &
                    field_displ, field_veloc, field_accel, field_press, field_stress, field_eps_pl, &
                    field_phi, field_vphi, DXX, DYY, DZZ, DXY, DYX, DXZ, DZX, DYZ, DZY)
            endif
            !GATHER FIELDS
            select case(domain_type)
                case (DM_SOLID) ! SOLID PART OF THE DOMAIN
                    if (out_variables(OUT_DEPLA  ) == 1) call gather_field(el, field_displ, Tdomain%champs0%Depla)
                    if (out_variables(OUT_VITESSE) == 1) call gather_field(el, field_veloc, Tdomain%champs0%Veloc)
                    if (out_variables(OUT_ACCEL  ) == 1) call gather_field(el, field_accel, Tdomain%champs0%Forces)
                    
                    if (Tdomain%nl_flag == 1) then ! NON LINEAR
                        if (out_variables(OUT_PRESSION) == 1 .or. flag_gradU /=0) then
                            call gather_field(el, field_stress, Tdomain%champs0%Stress)
                            do k = 0,ngllz-1
                                do j = 0,nglly-1
                                    do i = 0,ngllx-1
                                        if (out_variables(OUT_PRESSION) == 1) then
                                            field_press(i,j,k) = sum(field_stress(i,j,k,1:3))
                                        end if
                                        if (flag_gradU /= 0) then
                                            call gather_field(el, field_eps_pl, Tdomain%champs0%Epsilon_pl)
                                        end if
                                    end do
                                end do
                            end do
                        end if
                    else    ! ELASTIC
                        write(*,*) 'allocate_displ', allocated(field_displ)
                        if (out_variables(OUT_PRESSION) == 1) then
                            call pressure_solid(ngllx,nglly,ngllz,sub_dom_mat%htprimex, sub_dom_mat%hprimey,sub_dom_mat%hprimez, &
                                el%InvGrad,field_displ, el%Lambda, el%Mu, field_press)
                        end if
                        if (flag_gradU/=0) then
                            call grad_displ_solid(ngllx,nglly,ngllz,sub_dom_mat%htprimex,sub_dom_mat%hprimey,sub_dom_mat%hprimez, &
                                el%InvGrad,field_displ,dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz)
                        endif
                    end if
                write(*,*) 'allocate', allocated(dxx)
                case (DM_SOLID_PML) ! SOLID PML PART OF THE DOMAIN
                    call gather_field_pml(el, field_veloc, Tdomain%champs0%VelocPML)
                case (DM_FLUID) ! FLUID PART OF THE DOMAIN
                    call gather_field_fluid(el, field_phi, Tdomain%champs0%Phi)
                    call fluid_velocity(ngllx,nglly,ngllz,sub_dom_mat%htprimex, &
                        sub_dom_mat%hprimey,sub_dom_mat%hprimez, el%InvGrad,el%density,field_phi,field_veloc)
                    call gather_field_fluid(el, field_vphi, Tdomain%champs0%VelPhi)
                    call fluid_velocity(ngllx,nglly,ngllz,sub_dom_mat%htprimex, &
                        sub_dom_mat%hprimey,sub_dom_mat%hprimez, el%InvGrad,el%density,field_vphi,field_accel)
                    field_press = -field_vphi
                case (DM_FLUID_PML) !FLUID PML PART OF THE DOMAIN
                    call gather_field_fpml(el, field_phi, Tdomain%champs0%fpml_Phi)
                    call fluid_velocity(ngllx,nglly,ngllz,sub_dom_mat%htprimex, &
                        sub_dom_mat%hprimey,sub_dom_mat%hprimez, el%InvGrad,el%density,field_phi,field_veloc)
                    call gather_field_fpml(el, field_vphi, Tdomain%champs0%fpml_VelPhi)
                    call fluid_velocity(ngllx,nglly,ngllz,sub_dom_mat%htprimex, &
                        sub_dom_mat%hprimey,sub_dom_mat%hprimez, el%InvGrad,el%density,field_vphi,field_accel)
                    field_press = -field_vphi
            end select
            ! COMPUTING SNAPSHOTS
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        idx = irenum(el%Iglobnum(i,j,k))
                        if (domains(idx)==domain_type) then
                            valence(idx) = valence(idx)+1
                            select case(domain_type)
                                case (DM_SOLID,DM_SOLID_PML) ! SOLID PART OF THE DOMAIN
                                    if (out_variables(OUT_DEPLA)   ==1) out_fields%displ(:,idx) = field_displ(i,j,k,:)
                                    if (out_variables(OUT_VITESSE) ==1) out_fields%veloc(:,idx) = out_fields%veloc(:,idx) + field_veloc(i,j,k,:)
                                    if (out_variables(OUT_ACCEL)   ==1) out_fields%accel(:,idx) = out_fields%accel(:,idx) + field_accel(i,j,k,:)
                                    if (out_variables(OUT_PRESSION)==1) out_fields%press(idx)   = field_press(i,j,k)
                                    if (flag_gradU /= 0) then
                                        if (Tdomain%nl_flag==1) then ! NON LINEAR
                                            call stress_strain_nl(out_variables, out_fields, &
                                                DXX(i,j,k), DYY(i,j,k), DZZ(i,j,k), &
                                                DXY(i,j,k), DYX(i,j,k), DXZ(i,j,k), &
                                                DZX(i,j,k), DYZ(i,j,k), DZY(i,j,k), &
                                                field_stress(i,j,k,:),field_eps_pl(i,j,k,:), idx)
                                        else ! ELASTIC
                                            if (Tdomain%aniso) then
                                            else
                                                xmu     = el%Mu(i,j,k)
                                                xlambda = el%Lambda(i,j,k)
                                                xkappa  = el%Kappa(i,j,k)
                                                if (n_solid>0) then
                                                    onemSbeta=el%sl%onemSbeta(i,j,k)
                                                    onemPbeta=el%sl%onemPbeta(i,j,k)
                                                    !  mu_relaxed -> mu_unrelaxed
                                                    xmu    = xmu * onemSbeta
                                                    !  kappa_relaxed -> kappa_unrelaxed
                                                    xkappa = xkappa * onemPbeta
                                                endif
                                                x2mu       = 2. * xmu
                                                xlambda2mu = xlambda + x2mu

                                                call stress_strain_el(out_variables, out_fields, DXX(i,j,k), DYY(i,j,k), DZZ(i,j,k), &
                                                    DXY(i,j,k), DYX(i,j,k), DXZ(i,j,k), DZX(i,j,k), DYZ(i,j,k), DZY(i,j,k), &
                                                    xmu, x2mu, xlambda2mu, idx)
                                            end if
                                        end if
                                    end if

                                case (DM_FLUID,DM_FLUID_PML)
                                    out_fields%veloc(:,idx) = out_fields%veloc(:,idx) + field_veloc(i,j,k,:)
                                    out_fields%accel(:,idx) = out_fields%accel(:,idx) + field_accel(i,j,k,:)
                                    out_fields%press(idx)   = field_press(i,j,k)
                            end select

                        end if
                    enddo
                enddo
            enddo
        enddo

        ! normalization
        do i = 0,nnodes-1
            if (valence(i)/=0) then
                if (out_variables(OUT_VITESSE)==1) then
                    out_fields%veloc(0:2,i) = out_fields%veloc(0:2,i)/valence(i)
                end if
                if (out_variables(OUT_ACCEL)==1) then
                    out_fields%accel(0:2,i) = out_fields%accel(0:2,i)/valence(i)
                end if
            end if
        enddo

        if (Tdomain%output_rank==0) then
            group = Tdomain%rank/Tdomain%ngroup
            call semname_snap_result_file(group, isort, fnamef)
            call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        else
            fid = -1
        endif

        call grp_write_fields(Tdomain, fid, nnodes, out_variables, out_fields, nnodes_tot)

        if (Tdomain%output_rank==0) then
            call h5fclose_f(fid, hdferr)
            call write_xdmf(Tdomain, group, isort, nnodes_tot, out_variables)
        endif

        deallocate(valence)

        call deallocate_fields(out_variables, out_fields, Tdomain%nl_flag)

        if (allocated(field_displ))  deallocate(field_displ)
        if (allocated(field_veloc))  deallocate(field_veloc)
        if (allocated(field_accel))  deallocate(field_accel)
        if (allocated(field_press))  deallocate(field_press)
        if (allocated(field_stress)) deallocate(field_stress)
        if (allocated(field_eps_pl)) deallocate(field_eps_pl)
        if (allocated(field_phi))    deallocate(field_phi)
        if (allocated(field_vphi))   deallocate(field_vphi)
        if (allocated(DXX)) then
            deallocate(DXX)
            deallocate(DXY)
            deallocate(DXZ)
            deallocate(DYX)
            deallocate(DYY)
            deallocate(DYZ)
            deallocate(DZX)
            deallocate(DZY)
            deallocate(DZZ)
        end if
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

    subroutine write_xdmf(Tdomain, group, isort, nnodes, out_variables)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer, intent(in) :: group, isort, nnodes
        integer, dimension(0:8), intent(in) :: out_variables
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

        write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" Datatype="Float" Precision="8"  Dimensions="',nn, &
            '">geometry',group,'.h5:/Mass</DataItem>'
        write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Jac" Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn, &
            '">geometry',group,'.h5:/Jac</DataItem>'
        write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Dom" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',group,'.h5:/Dom</DataItem>'
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
            if (out_variables(OUT_DEPLA) == 1) then
                write(61,"(a,I9,a)") '<Attribute Name="Displ" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/displ'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! VELOC
            if (out_variables(OUT_VITESSE) == 1) then
                write(61,"(a,I9,a)") '<Attribute Name="Veloc" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/veloc'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! ACCEL
            if (out_variables(OUT_ACCEL) == 1) then
                write(61,"(a,I9,a)") '<Attribute Name="Accel" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/accel'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! PRESSURE
            if (out_variables(OUT_PRESSION) == 1) then
                write(61,"(a,I9,a)") '<Attribute Name="Pressure" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/press'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! VOLUMETRIC STRAIN
            if (out_variables(OUT_EPS_VOL) == 1) then
                write(61,"(a,I9,a)") '<Attribute Name="eps_vol" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_vol'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            if (out_variables(OUT_EPS_DEV) == 1) then
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
                if (Tdomain%nl_flag==1) then
                    ! EPS_DEV_PL_XX
                    write(61,"(a,I9,a)") '<Attribute Name="eps_dev_pl_xx" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_xx'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_XX
                    write(61,"(a,I9,a)") '<Attribute Name="eps_dev_pl_yy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_yy'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_ZZ
                    write(61,"(a,I9,a)") '<Attribute Name="eps_dev_pl_zz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_zz'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_XY
                    write(61,"(a,I9,a)") '<Attribute Name="eps_dev_pl_xy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_xy'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PLXZ
                    write(61,"(a,I9,a)") '<Attribute Name="eps_dev_pl_xz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_xz'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_YZ
                    write(61,"(a,I9,a)") '<Attribute Name="eps_dev_pl_yz" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_yz'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                end if
            end if

            if (out_variables(OUT_STRESS_DEV) == 1) then
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
            end if

            if (out_variables(OUT_ENERGYP) == 1) then
                ! P_ENERGY
                write(61,"(a,I9,a)") '<Attribute Name="P_energy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/P_energy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if

            if (out_variables(OUT_ENERGYS) == 1) then
                ! S_ENERGY
                write(61,"(a,I9,a)") '<Attribute Name="S_energy" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/S_energy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if

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
            write(61,"(a,I9,a)") '<Attribute Name="Dom" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',group, &
                '"]/DataItem[@Name="Dom"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '</Grid>'
            ! XXX inexact pour l'instant
            time = time+Tdomain%TimeD%time_snapshots
        end do
        write(61,"(a)") '</Grid>'
        close(61)
    end subroutine write_xdmf

    subroutine write_constant_fields(Tdomain, fid, irenum, nnodes, domains)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        integer, intent(in) :: nnodes
        integer, dimension(:), allocatable, intent(in) :: domains
        !
        real, dimension(:),allocatable :: mass, jac
        integer :: ngllx, nglly, ngllz, idx
        integer :: i, j, k, n, nnodes_tot
        integer :: domain_type, imat


        allocate(mass(0:nnodes-1))
        allocate(jac(0:nnodes-1))
        mass = 0d0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            imat = Tdomain%specel(n)%mat_index
            domain_type = get_domain(Tdomain%sSubDomain(imat))
            select case(domain_type)
                case (DM_SOLID)
                    do k = 0,ngllz-1
                        do j = 0,nglly-1
                            do i = 0,ngllx-1
                                idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                                if (domains(idx)==domain_type) then
                                    mass(idx) = Tdomain%MassMatSol(Tdomain%specel(n)%Idom(i,j,k))
                                endif
                            end do
                        end do
                    end do
                case (DM_SOLID_PML)
                    do k = 0,ngllz-1
                        do j = 0,nglly-1
                            do i = 0,ngllx-1
                                idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                                if (domains(idx)==domain_type) then
                                    mass(idx) = Tdomain%MassMatSolPML(Tdomain%specel(n)%Idom(i,j,k))
                                endif
                            end do
                        end do
                    end do
                case (DM_FLUID)
                    do k = 0,ngllz-1
                        do j = 0,nglly-1
                            do i = 0,ngllx-1
                                idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                                if (domains(idx)==domain_type) then
                                    mass(idx) = Tdomain%MassMatFlu(Tdomain%specel(n)%Idom(i,j,k))
                                endif
                            end do
                        end do
                    end do
                case (DM_FLUID_PML)
                    do k = 0,ngllz-1
                        do j = 0,nglly-1
                            do i = 0,ngllx-1
                                idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                                if (domains(idx)==domain_type) then
                                    mass(idx) = Tdomain%MassMatFluPML(Tdomain%specel(n)%Idom(i,j,k))
                                endif
                            end do
                        end do
                    end do
            end select
        end do
        if (Tdomain%ngll_pmls>0) deallocate(Tdomain%MassMatSolPml)
        if (Tdomain%ngll_pmlf>0) deallocate(Tdomain%MassMatFluPml)
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

        call grp_write_real_1d(Tdomain, fid, "Mass", nnodes, mass, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Jac", nnodes, jac, nnodes_tot)
        call grp_write_int_1d(Tdomain, fid, "Dom", nnodes, domains, nnodes_tot)
        deallocate(mass,jac)

    end subroutine write_constant_fields

    subroutine stress_strain_el(out_flags, fields, DXX, DYY, DZZ, DXY, DYX, DXZ, DZX, DYZ, DZY, &
        xmu, x2mu, xlambda2mu, idx)

        integer, dimension(0:8), intent(in) :: out_flags
        type (output_var_t), intent(inout) :: fields
        real, intent(in) :: DXX, DYY, DZZ, DXY, DYX, DXZ, DZX, DYZ, DZY, xmu, x2mu, xlambda2mu
        real    :: eps_trace
        integer :: idx

        eps_trace = DXX + DYY + DZZ

        if (out_flags(OUT_ENERGYP) == 1) then ! P_energy
            fields%P_energy(idx) = .5 * xlambda2mu * eps_trace**2
        end if

        if (out_flags(OUT_ENERGYS) == 1) then ! S_energy
            fields%S_energy(idx) =  xmu/2 * (&
                DXY**2 + DYX**2 &
                + DXZ**2 + DZX**2 &
                + DYZ**2 + DZY**2 &
                - 2 * DXY * DYX - 2 * DXZ * DZX - 2 * DYZ * DZY)
        end if

        if (out_flags(OUT_EPS_VOL) == 1) then ! volumetric strain
            fields%eps_vol(idx)    = eps_trace
        end if

        if (out_flags(OUT_EPS_DEV) == 1) then ! deviatoric strain
            fields%eps_dev_xx(idx) = DXX - eps_trace * M_1_3
            fields%eps_dev_yy(idx) = DYY - eps_trace * M_1_3
            fields%eps_dev_zz(idx) = DZZ - eps_trace * M_1_3
            fields%eps_dev_xy(idx) = 0.5 * (DXY + DYX)
            fields%eps_dev_xz(idx) = 0.5 * (DZX + DXZ)
            fields%eps_dev_yz(idx) = 0.5 * (DZY + DYZ)
        end if

        if (out_flags(OUT_STRESS_DEV) == 1) then ! deviatoric stress state
            fields%sig_dev_xx(idx) = x2mu * (DXX - eps_trace * M_1_3)
            fields%sig_dev_yy(idx) = x2mu * (DYY - eps_trace * M_1_3)
            fields%sig_dev_zz(idx) = x2mu * (DZZ - eps_trace * M_1_3)
            fields%sig_dev_xy(idx) = xmu  * (DXY + DYX)
            fields%sig_dev_xz(idx) = xmu  * (DXZ + DZX)
            fields%sig_dev_yz(idx) = xmu  * (DYZ + DZY)
        end if
    end subroutine stress_strain_el

    subroutine stress_strain_nl(out_flags, fields, dUxx, dUyy, dUzz, dUxy, dUyx, dUxz, dUzx, dUyz, dUzy, &
        stress, strain_pl, idx)

        integer, dimension(0:8), intent(in) :: out_flags
        type (output_var_t), intent(inout) :: fields
        real, intent(in) :: dUxx, dUyy, dUzz, dUxy, dUyx, dUxz, dUzx, dUyz, dUzy
        real, dimension(0:5), intent(in) :: stress, strain_pl
        real    :: eps_trace, eps_trace_pl, pressure
        integer :: idx

        eps_trace    = dUxx + dUyy + dUzz
        eps_trace_pl = 0
        pressure     = sum(stress(0:2)) * M_1_3

        ! A VOIR
        if (out_flags(OUT_ENERGYP) == 1) then ! P_energy
            fields%P_energy(idx) = .5 * pressure * eps_trace
        end if
        ! A VOIR
        if (out_flags(OUT_ENERGYS) == 1) then ! S_energy
            fields%S_energy(idx) = 0
        end if

        if (out_flags(OUT_EPS_VOL) == 1) then ! volumetric strain
            fields%eps_vol(idx) = eps_trace
        end if

        if (out_flags(OUT_EPS_DEV) == 1) then ! deviatoric strain
            fields%eps_dev_xx(idx) = dUxx - eps_trace * M_1_3
            fields%eps_dev_yy(idx) = dUyy - eps_trace * M_1_3
            fields%eps_dev_zz(idx) = dUzz - eps_trace * M_1_3
            fields%eps_dev_xy(idx) = 0.5 * (dUxy + dUyx)
            fields%eps_dev_xz(idx) = 0.5 * (dUzx + dUxz)
            fields%eps_dev_yz(idx) = 0.5 * (dUyz + dUzy)
            fields%eps_dev_pl_xx(idx) = strain_pl(0) - eps_trace_pl * M_1_3
            fields%eps_dev_pl_yy(idx) = strain_pl(1) - eps_trace_pl * M_1_3
            fields%eps_dev_pl_zz(idx) = strain_pl(2) - eps_trace_pl * M_1_3
            fields%eps_dev_pl_xy(idx) = strain_pl(3)
            fields%eps_dev_pl_xz(idx) = strain_pl(4)
            fields%eps_dev_pl_yz(idx) = strain_pl(5)
        end if

        if (out_flags(OUT_STRESS_DEV) == 1) then ! deviatoric stress state
            fields%sig_dev_xx(idx) = stress(0)-pressure
            fields%sig_dev_yy(idx) = stress(1)-pressure
            fields%sig_dev_zz(idx) = stress(2)-pressure
            fields%sig_dev_xy(idx) = stress(3)
            fields%sig_dev_xz(idx) = stress(4)
            fields%sig_dev_yz(idx) = stress(5)
        end if

    end subroutine stress_strain_nl

end module msnapshots

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
