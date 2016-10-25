!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP

module msnapshots
    use sdomain
    use hdf5
    use sem_hdf5
    use semdatafiles
    use mpi
    use deriv3d
    use sem_c_config, only : sem_mkdir
    use constants
    implicit none
#include "index.h"

!integrate

    type :: output_var_t
        ! Nodes fields
        real(fpp), dimension(:,:), allocatable :: displ, veloc, accel
        real(fpp), dimension(:)  , allocatable :: press_n
        ! Cell fields
        real(fpp), dimension(:)  , allocatable   :: press_c, eps_vol
        real(fpp), dimension(:,:), allocatable   :: eps_dev, sig_dev
        real(fpp), dimension(:), allocatable     :: P_energy, S_energy
        real(fpp), dimension(:,:), allocatable   :: eps_dev_pl
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

    subroutine grp_write_fields(Tdomain, parent_id, dim2_n, dim2_c,&
        out_variables, outputs, ntot_nodes, ntot_cells)

        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim2_c,dim2_n ! Nb (N)odes and (C)ells
        type(output_var_t), intent(in) :: outputs
        integer, dimension(0:), intent(in) :: out_variables
        integer, intent(out) :: ntot_nodes,ntot_cells
        !
        integer(HID_T) :: dset_id
        real, dimension(:,:), allocatable :: all_data_2d_n
        real, dimension(:), allocatable :: all_data_1d_n, all_data_1d_c
        integer, dimension(:), allocatable :: displs_n, displs_c
        integer, dimension(:), allocatable :: counts_n, counts_c
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr

        if (Tdomain%output_rank==0) then
            allocate(displs_n(0:Tdomain%nb_output_procs-1))
            allocate(displs_c(0:Tdomain%nb_output_procs-1))
            allocate(counts_n(0:Tdomain%nb_output_procs-1))
            allocate(counts_c(0:Tdomain%nb_output_procs-1))
        end if
        call MPI_Gather(dim2_n, 1, MPI_INTEGER, counts_n, 1, MPI_INTEGER, 0, Tdomain%comm_output, ierr)
        call MPI_Gather(dim2_c, 1, MPI_INTEGER, counts_c, 1, MPI_INTEGER, 0, Tdomain%comm_output, ierr)
        ! 1D FIELDS-NOEUD
        ntot_nodes = 0
        ntot_cells = 0
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs_n(n)   = ntot_nodes
                displs_c(n)   = ntot_cells
                ntot_nodes    = ntot_nodes + counts_n(n)
                ntot_cells    = ntot_cells + counts_c(n)
            end do
            allocate(all_data_1d_n(0:ntot_nodes-1))
            allocate(all_data_1d_c(0:ntot_cells-1))
        end if

        ! P_ENERGY
        if (out_variables(OUT_ENERGYP) == 1) then
            call MPI_Gatherv(outputs%P_energy, dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, MPI_DOUBLE_PRECISION, 0,&
                Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "P_energy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! S_ENERGY
        if (out_variables(OUT_ENERGYS) == 1) then
            call MPI_Gatherv(outputs%S_energy, dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c,MPI_DOUBLE_PRECISION, 0, &
                Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "S_energy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! VOL_EPS
        if (out_variables(OUT_EPS_VOL) == 1) then
            call MPI_Gatherv(outputs%eps_vol, dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c,MPI_DOUBLE_PRECISION, 0,&
                Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "eps_vol", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! PRESSION
        if (out_variables(OUT_PRESSION) == 1) then
            call MPI_Gatherv(outputs%press_c, dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, MPI_DOUBLE_PRECISION, 0,&
                Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "press_elem", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            call MPI_Gatherv(outputs%press_n, dim2_n, MPI_DOUBLE_PRECISION, &
                all_data_1d_n, counts_n, displs_n, MPI_DOUBLE_PRECISION, 0,&
                Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_nodes
                call create_dset(parent_id, "press_gll", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_n, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if

        ! EPS_DEV
        if (out_variables(OUT_EPS_DEV) == 1) then
            ! EPS_DEV_XX
            call MPI_Gatherv(outputs%eps_dev(0,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "eps_dev_xx", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_YY
            call MPI_Gatherv(outputs%eps_dev(1,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "eps_dev_yy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_ZZ
            call MPI_Gatherv(outputs%eps_dev(2,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "eps_dev_zz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_XY
            call MPI_Gatherv(outputs%eps_dev(3,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "eps_dev_xy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_XZ
            call MPI_Gatherv(outputs%eps_dev(4,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "eps_dev_xz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! EPS_DEV_YZ
            call MPI_Gatherv(outputs%eps_dev(5,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "eps_dev_yz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if

            if (out_variables(OUT_EPS_DEV_PL) == 1) then
                ! EPS_DEV_PL_XX
                call MPI_Gatherv(outputs%eps_dev_pl(0,:), dim2_c, MPI_DOUBLE_PRECISION, &
                    all_data_1d_c, counts_c, displs_c, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_cells
                    call create_dset(parent_id, "eps_dev_pl_xx", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_YY
                call MPI_Gatherv(outputs%eps_dev_pl(1,:), dim2_c, MPI_DOUBLE_PRECISION, &
                    all_data_1d_c, counts_c, displs_c, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_cells
                    call create_dset(parent_id, "eps_dev_pl_yy", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_ZZ
                call MPI_Gatherv(outputs%eps_dev_pl(2,:), dim2_c, MPI_DOUBLE_PRECISION, &
                    all_data_1d_c, counts_c, displs_c, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_cells
                    call create_dset(parent_id, "eps_dev_pl_zz", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_XY
                call MPI_Gatherv(outputs%eps_dev_pl(3,:), dim2_c, MPI_DOUBLE_PRECISION, &
                    all_data_1d_c, counts_c, displs_c, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_cells
                    call create_dset(parent_id, "eps_dev_pl_xy", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_XZ
                call MPI_Gatherv(outputs%eps_dev_pl(4,:), dim2_c, MPI_DOUBLE_PRECISION, &
                    all_data_1d_c, counts_c, displs_c, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_cells
                    call create_dset(parent_id, "eps_dev_pl_xz", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
                ! EPS_DEV_PL_YZ
                call MPI_Gatherv(outputs%eps_dev_pl(5,:), dim2_c, MPI_DOUBLE_PRECISION, &
                    all_data_1d_c, counts_c, displs_c, &
                    MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
                if (Tdomain%output_rank==0) then
                    dims(1) = ntot_cells
                    call create_dset(parent_id, "eps_dev_pl_yz", H5T_IEEE_F32LE, dims(1), dset_id)
                    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                end if
            end if
        end if

        if (out_variables(OUT_STRESS_DEV) == 1) then
            ! SIG_DEV_XX
            call MPI_Gatherv(outputs%sig_dev(0,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "sig_dev_xx", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_YY
            call MPI_Gatherv(outputs%sig_dev(1,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "sig_dev_yy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_ZZ
            call MPI_Gatherv(outputs%sig_dev(2,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "sig_dev_zz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_XY
            call MPI_Gatherv(outputs%sig_dev(3,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "sig_dev_xy", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_XZ
            call MPI_Gatherv(outputs%sig_dev(4,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "sig_dev_xz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
            ! SIG_DEV_YZ
            call MPI_Gatherv(outputs%sig_dev(5,:), dim2_c, MPI_DOUBLE_PRECISION, &
                all_data_1d_c, counts_c, displs_c, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = ntot_cells
                call create_dset(parent_id, "sig_dev_yz", H5T_IEEE_F32LE, dims(1), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d_c, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if

        ! 3D FIELDS
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs_n(n) = displs_n(n)*3
                counts_n(n) = counts_n(n)*3
            end do
            allocate(all_data_2d_n(0:2,0:ntot_nodes-1))
        end if

        ! VELOCITY
        if (out_variables(OUT_VITESSE) == 1) then
            call MPI_Gatherv(outputs%veloc, 3*dim2_n, MPI_DOUBLE_PRECISION, &
                all_data_2d_n, counts_n, displs_n, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = 3
                dims(2) = ntot_nodes
                call create_dset_2d(parent_id, "veloc", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d_n, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        ! DISPLACEMENT
        if (out_variables(OUT_DEPLA) == 1) then
            call MPI_Gatherv(outputs%displ, 3*dim2_n, MPI_DOUBLE_PRECISION, &
                all_data_2d_n, counts_n, displs_n, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = 3
                dims(2) = ntot_nodes
                call create_dset_2d(parent_id, "displ", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d_n, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if

        ! ACCELERATION
        if (out_variables(OUT_ACCEL) == 1) then
            call MPI_Gatherv(outputs%accel, 3*dim2_n, MPI_DOUBLE_PRECISION, &
                all_data_2d_n, counts_n, displs_n, &
                MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
            if (Tdomain%output_rank==0) then
                dims(1) = 3
                dims(2) = ntot_nodes
                call create_dset_2d(parent_id, "accel", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d_n, dims, hdferr)
                call h5dclose_f(dset_id, hdferr)
            end if
        end if
        if (Tdomain%output_rank==0) then
            deallocate(all_data_1d_n)
            deallocate(all_data_2d_n)
            deallocate(all_data_1d_c)
            deallocate(displs_n)
            deallocate(counts_n)
            deallocate(displs_c)
            deallocate(counts_c)
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

    subroutine compute_saved_elements(Tdomain, irenum, nnodes, ncells, domains)
        type (domain), intent (INOUT):: Tdomain
        integer, allocatable, dimension(:), intent(out) :: irenum ! maps Iglobnum to file node number
        integer, intent(out) :: nnodes, ncells
        integer :: n, i, j, k, ngll, ig, gn, ne
        !
        integer :: count
        integer :: ierr, domain_type, imat
        integer, dimension(:), allocatable, intent(out) :: domains
        ! on sauvegarde nsubelements a la place de le calculer
        ncells = Tdomain%n_hexa_local
        allocate(irenum(0:Tdomain%n_glob_points-1))

        irenum = -1
        ig = 0
        ne = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = 0
            select case (Tdomain%specel(n)%domain)
                 case (DM_SOLID)
                     ngll = Tdomain%sdom%ngll
                 case (DM_FLUID)
                     ngll = Tdomain%fdom%ngll
                 case (DM_SOLID_PML)
                     ngll = Tdomain%spmldom%ngll
                 case (DM_FLUID_PML)
                     ngll = Tdomain%fpmldom%ngll
            end select
            ne = ne + 1
            do k = 0,ngll - 1
                do j = 0,ngll - 1
                    do i = 0,ngll - 1
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
            ngll = 0
            select case (Tdomain%specel(n)%domain)
                 case (DM_SOLID)
                     ngll = Tdomain%sdom%ngll
                 case (DM_FLUID)
                     ngll = Tdomain%fdom%ngll
                 case (DM_SOLID_PML)
                     ngll = Tdomain%spmldom%ngll
                 case (DM_FLUID_PML)
                     ngll = Tdomain%fpmldom%ngll
            end select
            imat = Tdomain%specel(n)%mat_index
            domain_type = Tdomain%specel(n)%domain

            do k = 0,ngll - 1
                do j = 0,ngll - 1
                    do i = 0,ngll - 1
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
        integer :: nsubelements
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

        call compute_saved_elements(Tdomain, irenum, nnodes, nsubelements, domains)

        call write_global_nodes(Tdomain, fid, irenum, nnodes)

        call write_elem_connectivity(Tdomain, fid, irenum)

        call write_constant_fields(Tdomain, fid, irenum, nnodes, domains)

        if (Tdomain%output_rank==0) then
            call h5fclose_f(fid, hdferr)
        endif

        call mpi_gather(nnodes, 1, MPI_INTEGER, nodes_per_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (rg==0) call write_master_xdmf(Tdomain)
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

    !------------------------------------------------------
    subroutine write_elem_energy(Tdomain, fid)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid

        integer :: ngllx, nglly, ngllz

        real, dimension(:), allocatable :: En_S_int, En_P_int
        integer :: count
        integer :: i, j, k, n

        allocate( En_P_int(0:Tdomain%n_hexa_local-1))
        allocate( En_S_int(0:Tdomain%n_hexa_local-1))

        count = 0
        do n = 0,Tdomain%n_elem-1

            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = 0
            nglly = 0
            ngllz = 0

            select case (Tdomain%specel(n)%domain)
                 case (DM_SOLID)
                     ngllx = Tdomain%sdom%ngll
                 case (DM_FLUID)
                     ngllx = Tdomain%fdom%ngll
                 case (DM_SOLID_PML)
                     ngllx = Tdomain%spmldom%ngll
                 case (DM_FLUID_PML)
                     ngllx = Tdomain%fpmldom%ngll
                 case default
                     stop "unknown domain"
            end select

            nglly = ngllx
            ngllz = ngllx

            do k = 0,ngllz-2
                do j = 0,nglly-2
                    do i = 0,ngllx-2
                        !!Output of the integral for every element
                        En_S_int(count) = Tdomain%specel(n)%En_S_int !Output of the integral
                        En_P_int(count) = Tdomain%specel(n)%En_P_int
                        !!Output of the average for every sub-element
                        !En_S_int(count) = Tdomain%specel(n)%En_S_avg(i,j,k)
                        !En_P_int(count) = Tdomain%specel(n)%En_P_avg(i,j,k)
                        count=count+1
                    end do
                end do
            end do

        end do

        call grp_write_real_1d(Tdomain, fid, "En_S_int", count, En_S_int, Tdomain%n_hexa)
        call grp_write_real_1d(Tdomain, fid, "En_P_int", count, En_P_int, Tdomain%n_hexa)
    end subroutine write_elem_energy

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
            select case (Tdomain%specel(n)%domain)
                 case (DM_SOLID)
                     ngll(1,k) = Tdomain%sdom%ngll
                     ngll(2,k) = Tdomain%sdom%ngll
                     ngll(3,k) = Tdomain%sdom%ngll
                 case (DM_FLUID)
                     ngll(1,k) = Tdomain%fdom%ngll
                     ngll(2,k) = Tdomain%fdom%ngll
                     ngll(3,k) = Tdomain%fdom%ngll
                 case (DM_SOLID_PML)
                     ngll(1,k) = Tdomain%spmldom%ngll
                     ngll(2,k) = Tdomain%spmldom%ngll
                     ngll(3,k) = Tdomain%spmldom%ngll
                 case (DM_FLUID_PML)
                     ngll(1,k) = Tdomain%fpmldom%ngll
                     ngll(2,k) = Tdomain%fpmldom%ngll
                     ngll(3,k) = Tdomain%fpmldom%ngll
                 case default
                     stop "unknown domain"
            end select
            ! Max number of global points (can count elem vertices twice)
            nglobnum = nglobnum + ngll(1,k)*ngll(2,k)*ngll(3,k)
            ! Number of subelements
            count = count+(ngll(1,k)-1)*(ngll(2,k)-1)*(ngll(3,k)-1)
            k = k + 1
        enddo
        nb_elem = k
        Tdomain%n_hexa_local = count

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
            ngllx = 0; nglly = 0; ngllz = 0;
            select case (Tdomain%specel(n)%domain)
                 case (DM_SOLID)
                     ngllx = Tdomain%sdom%ngll
                     nglly = Tdomain%sdom%ngll
                     ngllz = Tdomain%sdom%ngll
                 case (DM_FLUID)
                     ngllx = Tdomain%fdom%ngll
                     nglly = Tdomain%fdom%ngll
                     ngllz = Tdomain%fdom%ngll
                 case (DM_SOLID_PML)
                     ngllx = Tdomain%spmldom%ngll
                     nglly = Tdomain%spmldom%ngll
                     ngllz = Tdomain%spmldom%ngll
                 case (DM_FLUID_PML)
                     ngllx = Tdomain%fpmldom%ngll
                     nglly = Tdomain%fpmldom%ngll
                     ngllz = Tdomain%fpmldom%ngll
                 case default
                     stop "unknown domain"
            end select
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

    subroutine allocate_fields(nnodes, nsubelements, out_flags, fields, nl_flag)
        type (output_var_t), intent(inout) :: fields
        integer, dimension(0:N_OUT_VARS-1), intent(in) :: out_flags
        integer, intent(in) :: nnodes,nsubelements
        logical, intent(in) :: nl_flag
        logical :: flag_gradU

        if (nl_flag) then
            flag_gradU = (out_flags(OUT_ENERGYP)     + &
                          out_flags(OUT_ENERGYS)     + &
                          out_flags(OUT_EPS_VOL)) /= 0
        else
            flag_gradU = (out_flags(OUT_PRESSION)    + &
                          out_flags(OUT_ENERGYP)     + &
                          out_flags(OUT_ENERGYS)     + &
                          out_flags(OUT_EPS_VOL)     + &
                          out_flags(OUT_EPS_DEV)     + &
                          out_flags(OUT_STRESS_DEV)) /= 0
        endif
        ! sortie par noeud
        if (out_flags(OUT_DEPLA     ) == 1) allocate(fields%displ(0:2,0:nnodes-1))
        if (out_flags(OUT_VITESSE   ) == 1) allocate(fields%veloc(0:2,0:nnodes-1))
        if (out_flags(OUT_ACCEL     ) == 1) allocate(fields%accel(0:2,0:nnodes-1))
        if (out_flags(OUT_PRESSION  ) == 1) allocate(fields%press_n(0:nnodes-1))
        ! sortie par element
        if (out_flags(OUT_ENERGYP   ) == 1) allocate(fields%P_energy(0:nsubelements-1))
        if (out_flags(OUT_ENERGYS   ) == 1) allocate(fields%S_energy(0:nsubelements-1))
        if (out_flags(OUT_EPS_VOL   ) == 1) allocate(fields%eps_vol(0:nsubelements-1))
        if (out_flags(OUT_PRESSION  ) == 1) allocate(fields%press_c(0:nsubelements-1))
        if (out_flags(OUT_EPS_DEV   ) == 1) allocate(fields%eps_dev(0:5,0:nsubelements-1))
        if (out_flags(OUT_STRESS_DEV) == 1) allocate(fields%sig_dev(0:5,0:nsubelements-1))
        if (out_flags(OUT_EPS_DEV_PL) == 1) allocate(fields%eps_dev_pl(0:5,0:nsubelements-1))
        ! initialize
        if (out_flags(OUT_DEPLA     ) == 1) fields%displ      = 0.
        if (out_flags(OUT_VITESSE   ) == 1) fields%veloc      = 0.
        if (out_flags(OUT_ACCEL     ) == 1) fields%accel      = 0.
        if (out_flags(OUT_ENERGYP   ) == 1) fields%P_energy   = 0.
        if (out_flags(OUT_ENERGYS   ) == 1) fields%S_energy   = 0.
        if (out_flags(OUT_EPS_VOL   ) == 1) fields%eps_vol    = 0.
        if (out_flags(OUT_PRESSION  ) == 1) fields%press_c    = 0.
        if (out_flags(OUT_PRESSION  ) == 1) fields%press_n    = 0.
        if (out_flags(OUT_EPS_DEV   ) == 1) fields%eps_dev    = 0.
        if (out_flags(OUT_STRESS_DEV) == 1) fields%sig_dev    = 0.
        if (out_flags(OUT_EPS_DEV_PL) == 1) fields%eps_dev_pl = 0.
        !
        return
        !
    end subroutine allocate_fields

    subroutine deallocate_fields(out_flags, fields)
        integer, dimension(0:10), intent(in) :: out_flags
        type (output_var_t), intent(inout) :: fields

        if (out_flags(OUT_DEPLA     ) == 1) deallocate(fields%displ)
        if (out_flags(OUT_VITESSE   ) == 1) deallocate(fields%veloc)
        if (out_flags(OUT_ACCEL     ) == 1) deallocate(fields%accel)
        if (out_flags(OUT_ENERGYP   ) == 1) deallocate(fields%P_energy)
        if (out_flags(OUT_ENERGYS   ) == 1) deallocate(fields%S_energy)
        if (out_flags(OUT_EPS_VOL   ) == 1) deallocate(fields%eps_vol)
        if (out_flags(OUT_PRESSION  ) == 1) deallocate(fields%press_c)
        if (out_flags(OUT_PRESSION  ) == 1) deallocate(fields%press_n)
        if (out_flags(OUT_EPS_DEV   ) == 1) deallocate(fields%eps_dev)
        if (out_flags(OUT_STRESS_DEV) == 1) deallocate(fields%sig_dev)
        if (out_flags(OUT_EPS_DEV_PL) == 1) deallocate(fields%eps_dev_pl)
        !
        return
        !
    end subroutine deallocate_fields

    subroutine integrate_on_element(ngll,jac,GLLw,input_field,output_integral)
        implicit none
        ! intent IN
        integer,intent(in) :: ngll
        real(fpp), dimension(0:ngll-1), intent(in) :: GLLw
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1),intent(in):: jac
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1),intent(in):: input_field
        ! intent INOUT
        !integer, intent(inout) :: count_subel
        !real(fpp), dimension(0:nsubelements-1), intent(inout) :: output_field
        real(fpp), intent(out) :: output_integral
        !
        integer :: i,j,k
        real(fpp) :: Whei, mult

        output_integral = 0.0d0
        ! calcul integrale
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    Whei = GLLw(i)*GLLw(j)*GLLw(k)
                    mult = Whei*jac(i,j,k)
                    output_integral = output_integral + mult*input_field(i,j,k)
                enddo
            enddo
        enddo

        return
    end subroutine integrate_on_element

    subroutine apply_integrated_value_on_output(ngll, input_integral, output_field, count_subel)

        implicit none
        !INPUT
        integer,intent(in) :: ngll
        real(fpp), intent(in) :: input_integral
        !OUTPUT
        integer, intent(inout) :: count_subel
        real(fpp), dimension(0:), intent(inout) :: output_field
        !LOCAL
        integer :: i,j,k


        do k = 0,ngll-2
            do j = 0,ngll-2
                do i = 0,ngll-2
                    ! Output of the integral for every element
                    output_field(count_subel) = input_integral !Integrated Value
                    count_subel=count_subel+1
                end do
            end do
        end do
    end subroutine apply_integrated_value_on_output

    subroutine evaluate_cell_centers(ngll, gllc, start, input_field, output_field)
        ! intent IN
        integer, intent(in)                             :: ngll
        real(fpp), intent(in), dimension(0:,0:,0:)      :: input_field
        real(fpp), intent(in), dimension(0:)            :: gllc
        integer, intent(in) :: start
        ! intent INOUT
        real(fpp), dimension(0:), intent(inout) :: output_field
        !
        integer         :: i, j, k, idx
        real(fpp)       :: xi, eta, zeta
        !

        idx = 0
        do k=0,ngll-2
            zeta = .5d0 * (gllc(k)+gllc(k+1))
            do j=0,ngll-2
                eta = .5d0 * (gllc(j)+gllc(j+1))
                do i = 0,ngll-2
                    xi = .5d0 * (gllc(i)+gllc(i+1))
                    output_field(start+idx) = evaluate_field(ngll,gllc,xi,eta,zeta,input_field)
                    idx = idx + 1
                end do
            end do
        end do
        !
        return
        !
    end subroutine evaluate_cell_centers

    function evaluate_field(ngll,gllc,xi,eta,zeta,field) result(r)
        ! intent IN
        integer,   intent(in)                           :: ngll
        real(fpp)                                       :: xi, eta, zeta
        real(fpp), intent(in), dimension(0:)            :: gllc
        real(fpp), intent(in), dimension(0:,0:,0:)      :: field
        !intent OUT
        real(fpp) :: r, weight
        real(fpp), dimension(0:ngll-1) :: outx, outy, outz
        integer :: i,j,k
        !
        r = 0.0d0
        do i=0,ngll-1
            call pol_lagrange(ngll,gllc,i,xi  ,outx(i)) ! P_i(xi)
            call pol_lagrange(ngll,gllc,i,eta ,outy(i)) ! P_i(eta)
            call pol_lagrange(ngll,gllc,i,zeta,outz(i)) ! P_i(zeta)
        end do

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    weight = outx(i)*outy(j)*outz(k)
                    r = r + (field(i,j,k) * weight)
                end do
            end do
        end do
        !
        return
        !
    end function evaluate_field


    subroutine save_field_h5(Tdomain, isort)
        use sdomain
        use dom_solid
        use dom_fluid
        use dom_solidpml
        use dom_fluidpml

        implicit none

        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: isort
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        integer :: domain_type
        integer, dimension(:), allocatable :: valence
        integer :: hdferr
        integer :: ngll
        integer :: i, j, k, n, m, ind
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes, ncells, group, nnodes_tot,nelements_tot
        integer, dimension(:), allocatable :: domains
        type(Element), pointer :: el
        type(subdomain), pointer :: sub_dom_mat
        type(output_var_t) :: out_fields
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev,eps_dev_pl
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        integer :: bnum, ee
        real, dimension(:), allocatable :: GLLc ! GLLw
        real(fpp), dimension(:,:,:), allocatable :: jac

        integer, dimension(0:size(Tdomain%out_variables)-1) :: out_variables
        logical :: nl_flag

        integer :: cell_start

        nl_flag=Tdomain%nl_flag

        out_variables(:) = Tdomain%out_variables(:)

        call create_dir_sorties(Tdomain, isort)
        call compute_saved_elements(Tdomain, irenum, nnodes, ncells, domains)
        call allocate_fields(nnodes, ncells, Tdomain%out_variables, out_fields, nl_flag)
        allocate(valence(0:nnodes-1))

        valence(:) = 0
        ngll = 0
        cell_start = 0
        do n = 0,Tdomain%n_elem-1
            el => Tdomain%specel(n)
            sub_dom_mat => Tdomain%sSubdomain(el%mat_index)
            if (.not. el%OUTPUT ) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            domain_type = Tdomain%specel(n)%domain
            select case(domain_type)
                case (DM_SOLID)
                    if (Tdomain%sdom%PlaneW%Exist) call compute_planeW_Exafield(el%lnum,Tdomain%TimeD%rtime,Tdomain)
                    call get_solid_dom_var(Tdomain%sdom, el%lnum, out_variables,    &
                        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol,&
                        eps_dev, sig_dev, nl_flag, eps_dev_pl)
                case (DM_FLUID)
                    call get_fluid_dom_var(Tdomain%fdom, el%lnum, out_variables,        &
                        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
                case (DM_SOLID_PML)
                    call get_solidpml_dom_var(Tdomain%spmldom, el%lnum, out_variables,           &
                        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
                case (DM_FLUID_PML)
                  call get_fluidpml_dom_var(Tdomain%fpmldom, el%lnum, out_variables,           &
                  fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
                case default
                  stop "unknown domain"
            end select

            bnum = el%lnum/VCHUNK
            ee = mod(el%lnum,VCHUNK)
            allocate(jac(0:ngll-1,0:ngll-1,0:ngll-1))
            jac (:,:,:) = 0.0d0
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        ind = irenum(el%Iglobnum(i,j,k))
                        valence(ind) = valence(ind)+1
                        select case (domain_type)
                            case (DM_SOLID)
                                jac(i,j,k) = Tdomain%sdom%Jacob_   (i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
                                jac(i,j,k) = Tdomain%spmldom%Jacob_(i,j,k,bnum,ee)
                            case (DM_FLUID)
                                jac(i,j,k) = Tdomain%fdom%Jacob_   (i,j,k,bnum,ee)
                            case (DM_FLUID_PML)
                                jac(i,j,k) = Tdomain%fpmldom%Jacob_(i,j,k,bnum,ee)
                            case default
                                stop "unknown domain"
                        end select

                        ! sortie par noeud
                        if (out_variables(OUT_DEPLA)   == 1)  out_fields%displ(0:2,ind) = fieldU(i,j,k,0:2)
                        if (out_variables(OUT_VITESSE) == 1)  out_fields%veloc(0:2,ind) = fieldV(i,j,k,0:2)
                        if (out_variables(OUT_ACCEL)   == 1)  out_fields%accel(0:2,ind) = fieldA(i,j,k,0:2)
                        if (out_variables(OUT_PRESSION) == 1) out_fields%press_n(ind) = fieldP(i,j,k)
                    enddo
                enddo
            enddo
            ! sortie de la valeur au centre par sub-element
            !call domain_gllw(Tdomain, domain_type, GLLw)
            call domain_gllc(Tdomain, Tdomain%specel(n)%domain, GLLc)

            if (out_variables(OUT_PRESSION) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, fieldP, out_fields%press_c)
            endif
            if (out_variables(OUT_ENERGYP) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, P_energy, out_fields%P_energy)
            endif
            if (out_variables(OUT_ENERGYS) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, S_energy, out_fields%S_energy)
            endif
            if (out_variables(OUT_EPS_VOL) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, eps_vol, out_fields%eps_vol)
            endif
            if (out_variables(OUT_EPS_DEV) == 1) then
                do m = 0,OUT_VAR_DIMS_3D(OUT_EPS_DEV)-1
                    call evaluate_cell_centers(ngll, GLLc, cell_start, eps_dev(:,:,:,m), out_fields%eps_dev(m,:))
                end do
            endif
            if (out_variables(OUT_EPS_DEV_PL) == 1) then
                do m = 0,OUT_VAR_DIMS_3D(OUT_EPS_DEV_PL)-1
                    call evaluate_cell_centers(ngll, GLLc, cell_start, eps_dev_pl(:,:,:,m), out_fields%eps_dev_pl(m,:))
                end do
            endif
            if (out_variables(OUT_STRESS_DEV) == 1) then
                do m = 0,OUT_VAR_DIMS_3D(OUT_STRESS_DEV)-1
                    call evaluate_cell_centers(ngll, GLLc, cell_start, sig_dev(:,:,:,m), out_fields%sig_dev(m,:))
                end do
            endif
            cell_start = cell_start + (ngll-1)**3
            if(allocated(jac)) deallocate(jac)
            if(allocated(GLLc)) deallocate(GLLc)
            !if(allocated(GLLw)) deallocate(GLLw)
        enddo
        if(allocated(fieldU))       deallocate(fieldU)
        if(allocated(fieldV))       deallocate(fieldV)
        if(allocated(fieldA))       deallocate(fieldA)
        if(allocated(fieldP))       deallocate(fieldP)
        if(allocated(P_energy))     deallocate(P_energy)
        if(allocated(S_energy))     deallocate(S_energy)
        if(allocated(eps_vol))      deallocate(eps_vol)
        if(allocated(eps_dev))      deallocate(eps_dev)
        if(allocated(eps_dev_pl))   deallocate(eps_dev_pl)
        if(allocated(sig_dev))      deallocate(sig_dev)

        if (Tdomain%output_rank==0) then
            group = Tdomain%rank/Tdomain%ngroup
            call semname_snap_result_file(group, isort, fnamef)
            call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        else
            fid = -1
        endif
        call grp_write_fields(Tdomain, fid, nnodes, ncells, &
            out_variables, out_fields, nnodes_tot, nelements_tot)
        if (Tdomain%output_rank==0) then
            call h5fclose_f(fid, hdferr)
            call write_xdmf(Tdomain, group, isort, nnodes_tot, out_variables)
        endif

        deallocate(valence)
        call deallocate_fields(out_variables, out_fields)
        call mpi_barrier(Tdomain%communicateur, hdferr)

    end subroutine save_field_h5

    subroutine write_master_xdmf(Tdomain)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: n_procs, nelem, n_groups
        character (len=MAX_FILE_SIZE) :: fnamef
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
            write(61,"(a,I4.4,a)") '<xi:include href="mesh.',group,'.xmf" xpointer="xpointer(//Xdmf/Domain/Grid)"/>'
        end do

        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
        close(61)

    end subroutine write_master_xdmf

    subroutine write_xdmf(Tdomain, group, isort,&
        nnodes, out_variables)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer, intent(in) :: group, isort, nnodes
        integer, dimension(0:), intent(in) :: out_variables
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
        write(61,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write(61,"(a)") '<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
        write(61,"(a)") '<Domain>'
        write(61,"(a,I4.4,a)") '<Grid CollectionType="Temporal" GridType="Collection" Name="space.',group,'">'

        time = 0

        do i=1,isort
            write(61,"(a,I4.4,a,I4.4,a)") '<Grid Name="mesh.',i,'.',group,'">'
            write(61,"(a,F20.10,a)") '<Time Value="', time,'"/>'
            write(61,"(a,I9,a)") '<Topology TopologyType="Hexahedron" NumberOfElements="',ne,'">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Int" Dimensions="',ne,' 8">'
            write(61,"(a,I4.4,a)") 'geometry',group,'.h5:/Elements'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Topology>'
            write(61,"(a)") '<Geometry GeometryType="XYZ">'
            write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a)") 'geometry',group,'.h5:/Nodes'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Geometry>'
            ! DISPL
            if (out_variables(OUT_DEPLA) == 1) then
                write(61,"(a)") '<Attribute Name="Displ" Center="Node" AttributeType="Vector">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn,' 3">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/displ'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! VELOC
            if (out_variables(OUT_VITESSE) == 1) then
                write(61,"(a)") '<Attribute Name="Veloc" Center="Node" AttributeType="Vector">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn,' 3">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/veloc'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! ACCEL
            if (out_variables(OUT_ACCEL) == 1) then
                write(61,"(a)") '<Attribute Name="Accel" Center="Node" AttributeType="Vector">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn,' 3">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/accel'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! PRESSURE
            if (out_variables(OUT_PRESSION) == 1) then
                write(61,"(a)") '<Attribute Name="Press_gll" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/press_gll'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            if (out_variables(OUT_PRESSION) == 1) then
                write(61,"(a)") '<Attribute Name="Press_elem" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/press_elem'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! VOLUMETRIC STRAIN
            if (out_variables(OUT_EPS_VOL) == 1) then
                write(61,"(a)") '<Attribute Name="eps_vol" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_vol'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! DEVIATORIC STRAIN
            if (out_variables(OUT_EPS_DEV) == 1) then
                ! EPS_DEV_XX
                write(61,"(a)") '<Attribute Name="eps_dev_xx" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_xx'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! EPS_DEV_XX
                write(61,"(a)") '<Attribute Name="eps_dev_yy" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_yy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! EPS_DEV_ZZ
                write(61,"(a)") '<Attribute Name="eps_dev_zz" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_zz'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! EPS_DEV_XY
                write(61,"(a)") '<Attribute Name="eps_dev_xy" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_xy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! EPS_DEV_XZ
                write(61,"(a)") '<Attribute Name="eps_dev_xz" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_xz'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! EPS_DEV_YZ
                write(61,"(a)") '<Attribute Name="eps_dev_yz" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_yz'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                if (Tdomain%nl_flag) then
                    ! EPS_DEV_PL_XX
                    write(61,"(a)") '<Attribute Name="eps_dev_pl_xx" Center="Cell" AttributeType="Scalar">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="4" Dimensions="',ne,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_xx'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_XX
                    write(61,"(a)") '<Attribute Name="eps_dev_pl_yy" Center="Cell" AttributeType="Scalar">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="4" Dimensions="',ne,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_yy'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_ZZ
                    write(61,"(a)") '<Attribute Name="eps_dev_pl_zz" Center="Cell" AttributeType="Scalar">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="4" Dimensions="',ne,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_zz'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_XY
                    write(61,"(a)") '<Attribute Name="eps_dev_pl_xy" Center="Cell" AttributeType="Scalar">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="4" Dimensions="',ne,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_xy'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_XZ
                    write(61,"(a)") '<Attribute Name="eps_dev_pl_xz" Center="Cell" AttributeType="Scalar">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="4" Dimensions="',ne,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_xz'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                    ! EPS_DEV_PL_YZ
                    write(61,"(a)") '<Attribute Name="eps_dev_pl_yz" Center="Cell" AttributeType="Scalar">'
                    write(61,"(a,I9,a)") '<DataItem Format="HDF" Datatype="Float" Precision="4" Dimensions="',ne,'">'
                    write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/eps_dev_pl_yz'
                    write(61,"(a)") '</DataItem>'
                    write(61,"(a)") '</Attribute>'
                end if
            end if
            ! DEVIATORIC STRESS
            if (out_variables(OUT_STRESS_DEV) == 1) then
                ! SIG_DEV_XX
                write(61,"(a)") '<Attribute Name="sig_dev_xx" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_xx'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! SIG_DEV_XX
                write(61,"(a)") '<Attribute Name="sig_dev_yy" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_yy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! SIG_DEV_ZZ
                write(61,"(a)") '<Attribute Name="sig_dev_zz" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_zz'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! SIG_DEV_XY
                write(61,"(a)") '<Attribute Name="sig_dev_xy" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_xy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! SIG_DEV_XZ
                write(61,"(a)") '<Attribute Name="sig_dev_xz" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_xz'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
                ! SIG_DEV_YZ
                write(61,"(a)") '<Attribute Name="sig_dev_yz" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/sig_dev_yz'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! P_ENERGY
            if (out_variables(OUT_ENERGYP) == 1) then
                write(61,"(a)") '<Attribute Name="P_energy" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/P_energy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! S_ENERGY
            if (out_variables(OUT_ENERGYS) == 1) then
                write(61,"(a)") '<Attribute Name="S_energy" Center="Cell" AttributeType="Scalar">'
                write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="',ne,'">'
                write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',group,'.h5:/S_energy'
                write(61,"(a)") '</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! DOMAIN
            write(61,"(a)") '<Attribute Name="Domain" Center="Grid" AttributeType="Scalar">'
            write(61,"(a,I4,a)") '<DataItem Format="XML" NumberType="Int"  Dimensions="1">',group,'</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! MAT
            write(61,"(a)") '<Attribute Name="Mat" Center="Cell" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mat" Format="HDF" NumberType="Int"  Dimensions="',ne, &
                '">geometry',group,'.h5:/Material</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! PROC
            write(61,"(a)") '<Attribute Name="Proc" Center="Cell" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Proc" Format="HDF" NumberType="Int"  Dimensions="',ne, &
                '">geometry',group,'.h5:/Proc</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! MASS
            write(61,"(a)") '<Attribute Name="Mass" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" NumberType="Float" Precision="4"  Dimensions="',nn, &
                '">geometry',group,'.h5:/Mass</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! ALPHA/DUMPSX
            write(61,"(a)") '<Attribute Name="Alpha" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Alpha" Format="HDF" NumberType="Float" Precision="4"  Dimensions="',nn, &
                '">geometry',group,'.h5:/Alpha</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! JACOBIAN
            write(61,"(a)") '<Attribute Name="Jac" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Jac" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                '">geometry',group,'.h5:/Jac</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! DENSITY
            write(61,"(a)") '<Attribute Name="Dens" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Dens" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                '">geometry',group,'.h5:/Dens</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! LAMBDA
            write(61,"(a)") '<Attribute Name="Lamb" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Lamb" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                '">geometry',group,'.h5:/Lamb</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! MU
            write(61,"(a)") '<Attribute Name="Mu" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mu" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                '">geometry',group,'.h5:/Mu</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! KAPPA
            write(61,"(a)") '<Attribute Name="Kappa" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Kappa" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                '">geometry',group,'.h5:/Kappa</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! DOMAIN
            write(61,"(a)") '<Attribute Name="Dom" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Dom" Format="HDF" NumberType="Int"  Dimensions="',nn, &
                '">geometry',group,'.h5:/Dom</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '</Grid>'
            ! XXX inexact pour l'instant
            time = time+Tdomain%TimeD%time_snapshots
        end do
        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
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
        real, dimension(:),allocatable :: mass, jac, dumpsx
        real, dimension(:),allocatable :: dens, lamb, mu, kappa
        integer :: ngll, idx
        integer :: i, j, k, n, lnum, nnodes_tot, bnum, ee
        integer :: domain_type, imat
        real(fpp) :: dx, dy, dz, dt

        allocate(mass(0:nnodes-1))
        allocate(dumpsx(0:nnodes-1))
        allocate(jac(0:nnodes-1))
        allocate(dens(0:nnodes-1))
        allocate(lamb(0:nnodes-1))
        allocate(mu(0:nnodes-1))
        allocate(kappa(0:nnodes-1))

        mass = 0d0
        dumpsx = 0d0
        dens = 0d0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            imat = Tdomain%specel(n)%mat_index
            lnum = Tdomain%specel(n)%lnum
            bnum = lnum/VCHUNK
            ee = mod(lnum,VCHUNK)

            domain_type = Tdomain%specel(n)%domain
            select case(domain_type)
            case (DM_SOLID)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (domains(idx)==domain_type) then
                                mass(idx) = Tdomain%sdom%MassMat(Tdomain%sdom%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case (DM_SOLID_PML)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (domains(idx)==domain_type) then
                                mass(idx) = Tdomain%spmldom%MassMat(Tdomain%spmldom%Idom_(i,j,k,bnum,ee))
                                dt = 2d0*Tdomain%TimeD%dtmin
#ifdef CPML
                                dumpsx(idx) = 0.
#else
                                dx = ((1d0/Tdomain%spmldom%PMLDumpSx_(i,j,k,1,bnum,ee))-1.)/dt
                                dy = ((1d0/Tdomain%spmldom%PMLDumpSy_(i,j,k,1,bnum,ee))-1.)/dt
                                dz = ((1d0/Tdomain%spmldom%PMLDumpSz_(i,j,k,1,bnum,ee))-1.)/dt
                                dumpsx(idx) = dx+dy+dz
#endif
                            endif
                        end do
                    end do
                end do
            case (DM_FLUID)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (domains(idx)==domain_type) then
                                mass(idx) = Tdomain%fdom%MassMat(Tdomain%fdom%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case (DM_FLUID_PML)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (domains(idx)==domain_type) then
                                mass(idx) = Tdomain%fpmldom%MassMat(Tdomain%fpmldom%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case default
                stop "unknown domain"
            end select
        enddo
        ! We need that for computing accel
!        if (Tdomain%spmldom%ngll>0) deallocate(Tdomain%spmldom%MassMat)
!        if (Tdomain%fpmldom%ngll>0) deallocate(Tdomain%fpmldom%MassMat)
        ! jac
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                jac(idx) = Tdomain%sdom%Jacob_   (i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
                                jac(idx) = Tdomain%spmldom%Jacob_(i,j,k,bnum,ee)
                            case (DM_FLUID)
                                jac(idx) = Tdomain%fdom%Jacob_   (i,j,k,bnum,ee)
                            case (DM_FLUID_PML)
                                jac(idx) = Tdomain%fpmldom%Jacob_(i,j,k,bnum,ee)
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
        end do

        ! dens
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                dens(idx) = Tdomain%sdom%Density_        (i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
#ifdef CPML
                                dens(idx) = Tdomain%spmldom%sSubDomain(Tdomain%specel(n)%mat_index)%DDensity
#else
                                dens(idx) = Tdomain%spmldom%Density_     (i,j,k,bnum,ee)
#endif
                            case (DM_FLUID)
                                dens(idx) = 1.0D0/Tdomain%fdom%IDensity_ (i,j,k,bnum,ee)
                            case (DM_FLUID_PML)
#ifdef CPML
                                dens(idx) = 0. ! Tdomain%fpmldom%sSubDomain(Tdomain%specel(n)%mat_index)%DDensity ! TODO
#else
                                dens(idx) = Tdomain%fpmldom%Density_     (i,j,k,bnum,ee)
#endif
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
        end do

        ! lamb
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                lamb(idx) = Tdomain%sdom%Lambda_        (i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
#ifdef CPML
                                lamb(idx) = Tdomain%spmldom%sSubDomain(Tdomain%specel(n)%mat_index)%DLambda
#else
                                lamb(idx) = Tdomain%spmldom%Lambda_     (i,j,k,bnum,ee)
#endif
                            case (DM_FLUID)
                                lamb(idx) = Tdomain%fdom%Lambda_        (i,j,k,bnum,ee)
                            case (DM_FLUID_PML)
                                lamb(idx) = Tdomain%fpmldom%Lambda_     (i,j,k,bnum,ee)
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
        end do

        ! mu
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                mu(idx) = Tdomain%sdom%Mu_(i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
#ifdef CPML
                                mu(idx) = Tdomain%spmldom%sSubDomain(Tdomain%specel(n)%mat_index)%DMu
#else
                                mu(idx) = Tdomain%spmldom%Mu_(i,j,k,bnum,ee)
#endif
                            case (DM_FLUID)
                                mu(idx) = -1d0
                            case (DM_FLUID_PML)
                                mu(idx) = -1d0
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
        end do

        ! kappa
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                kappa(idx) = Tdomain%sdom%Kappa_(i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
                                kappa(idx) = -1d0
                            case (DM_FLUID)
                                kappa(idx) = -1d0
                            case (DM_FLUID_PML)
                                kappa(idx) = -1d0
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
        end do

        call grp_write_real_1d(Tdomain, fid, "Mass", nnodes, mass, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Alpha", nnodes, dumpsx, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Jac", nnodes, jac, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Dens", nnodes, dens, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Lamb", nnodes, lamb, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Mu", nnodes, mu, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Kappa", nnodes, kappa, nnodes_tot)
        call grp_write_int_1d(Tdomain, fid, "Dom", nnodes, domains, nnodes_tot)
        deallocate(mass,jac)
        deallocate(dens, lamb, mu, kappa)

    end subroutine write_constant_fields

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
