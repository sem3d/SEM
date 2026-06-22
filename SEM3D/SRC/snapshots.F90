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
    use msnapdata, only : output_var_t
    implicit none
#include "index.h"
contains
#define DEBUG_CPML
#ifdef SINGLEPRECISION
#define MPI_REAL_FPP MPI_FLOAT
#else
#define MPI_REAL_FPP MPI_DOUBLE_PRECISION
#endif

    subroutine grp_write_real_2d(outputs, parent_id, name, dim1, dim2, data, ntot_nodes)
        type (output_var_t), intent (INOUT):: outputs
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim1, dim2
        real(fpp), dimension(0:dim1-1,0:dim2-1), intent(in) :: data
        character(len=*), INTENT(IN) :: name
        integer, intent(out) :: ntot_nodes
        !
        integer(HID_T) :: dset_id
        real(fpp), dimension(:,:), allocatable :: all_data
        integer, dimension(:), allocatable :: displs, counts
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr

        if (outputs%rank==0) then
            allocate(displs(0:outputs%nprocs-1))
            allocate(counts(0:outputs%nprocs-1))
        end if
        call MPI_Gather(dim2, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, outputs%comm, ierr)

        ntot_nodes = 0
        if (outputs%rank==0) then
            do n=0, outputs%nprocs-1
                displs(n) = dim1*ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
                counts(n)=counts(n)*dim1 ! for the next gatherv
            end do
            allocate(all_data(0:dim1-1,0:ntot_nodes-1))
        end if

        call MPI_Gatherv(data, dim1*dim2, MPI_REAL_FPP, all_data, counts, displs, &
            MPI_REAL_FPP, 0, outputs%comm, ierr)
        if (outputs%rank==0) then
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

    subroutine write_1d_var_c(outputs, parent_id, fname, field)
        type(output_var_t), intent(inout) :: outputs
        integer(HID_T), intent(in) :: parent_id
        character(len=*), intent(in) :: fname
        real(fpp), intent(in), dimension(:) :: field
        !
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: dset_id
        integer :: hdferr, ierr

        call MPI_Gatherv(field, outputs%ncells, MPI_REAL_FPP, &
            outputs%all_data_1d_c, outputs%counts_c, outputs%displs_c, MPI_REAL_FPP, 0,&
            outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = outputs%ntot_cells
            call create_dset(parent_id, trim(fname), H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, outputs%all_data_1d_c, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
    end subroutine write_1d_var_c

    subroutine write_2d_var_vecc(outputs, parent_id, fname, field)
        type(output_var_t), intent(inout) :: outputs
        integer(HID_T), intent(in) :: parent_id
        character(len=*), intent(in) :: fname
        real(fpp), intent(in), dimension(:,:) :: field
        !
        integer(HSIZE_T), dimension(2) :: dims
        integer(HID_T) :: dset_id
        integer :: hdferr, ierr

        call MPI_Gatherv(field, outputs%ncells, MPI_REAL_FPP, &
            outputs%all_data_2d_c, outputs%counts_c, outputs%displs_c, MPI_REAL_FPP, 0,&
            outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = 3
            dims(2) = outputs%ntot_cells
            call create_dset_2d(parent_id, trim(fname), H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, outputs%all_data_2d_c, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
    end subroutine write_2d_var_vecc

    subroutine write_1d_var_n(outputs, parent_id, fname, field)
        type(output_var_t), intent(inout) :: outputs
        integer(HID_T), intent(in) :: parent_id
        character(len=*), intent(in) :: fname
        real(fpp), intent(in), dimension(:) :: field
        !
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: dset_id
        integer :: hdferr, ierr

        call MPI_Gatherv(field, outputs%nnodes, MPI_REAL_FPP, &
            outputs%all_data_1d_n, outputs%counts_n, outputs%displs_n, MPI_REAL_FPP, 0,&
            outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = outputs%ntot_nodes
            call create_dset(parent_id, trim(fname), H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, outputs%all_data_1d_n, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
    end subroutine write_1d_var_n

    subroutine write_2d_var_vecn(outputs, parent_id, fname, field)
        type(output_var_t), intent(inout) :: outputs
        integer(HID_T), intent(in) :: parent_id
        character(len=*), intent(in) :: fname
        real(fpp), intent(in), dimension(:,:) :: field
        !
        integer(HSIZE_T), dimension(2) :: dims
        integer(HID_T) :: dset_id
        integer :: hdferr, ierr

        call MPI_Gatherv(field, 3*outputs%nnodes, MPI_REAL_FPP, &
            outputs%all_data_2d_n, outputs%counts2d_n, outputs%displs2d_n, MPI_REAL_FPP, 0,&
            outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = 3
            dims(2) = outputs%ntot_nodes
            call create_dset_2d(parent_id, trim(fname), H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, outputs%all_data_2d_n, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
    end subroutine write_2d_var_vecn

    subroutine grp_write_fields(outputs, parent_id, out_variables)
        type(output_var_t), intent(inout) :: outputs
        integer(HID_T), intent(in) :: parent_id
        integer, dimension(0:), intent(in) :: out_variables
        !

        ! P_ENERGY
        if (out_variables(OUT_ENERGYP) == 1) then
            call write_1d_var_c(outputs, parent_id, "P_energy", outputs%P_energy)
        end if
        ! S_ENERGY
        if (out_variables(OUT_ENERGYK) == 1) then
            call write_1d_var_c(outputs, parent_id, "K_energy", outputs%K_energy)
        end if
        ! VOL_EPS
        if (out_variables(OUT_EPS_VOL) == 1) then
            call write_1d_var_c(outputs, parent_id, "eps_vol", outputs%eps_vol)
        end if
        ! PRESSION
        if (out_variables(OUT_PRESSION) == 1) then
            call write_1d_var_c(outputs, parent_id, "press_elem", outputs%press_c)
            call write_1d_var_n(outputs, parent_id, "press_gll", outputs%press_n)
        end if

        ! DUDX
        if (out_variables(OUT_DUDX) == 1) then
            call write_1d_var_n(outputs, parent_id, "dUxdx", outputs%dUdX(0,:))
            call write_1d_var_n(outputs, parent_id, "dUxdy", outputs%dUdX(1,:))
            call write_1d_var_n(outputs, parent_id, "dUxdz", outputs%dUdX(2,:))
            call write_1d_var_n(outputs, parent_id, "dUydx", outputs%dUdX(3,:))
            call write_1d_var_n(outputs, parent_id, "dUydy", outputs%dUdX(4,:))
            call write_1d_var_n(outputs, parent_id, "dUydz", outputs%dUdX(5,:))
            call write_1d_var_n(outputs, parent_id, "dUzdx", outputs%dUdX(6,:))
            call write_1d_var_n(outputs, parent_id, "dUzdy", outputs%dUdX(7,:))
            call write_1d_var_n(outputs, parent_id, "dUzdz", outputs%dUdX(8,:))
        end if

        ! EPS_DEV
        if (out_variables(OUT_EPS_DEV) == 1) then
            call write_1d_var_c(outputs, parent_id, "eps_dev_xx", outputs%eps_dev(0,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_yy", outputs%eps_dev(1,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_zz", outputs%eps_dev(2,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_xy", outputs%eps_dev(3,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_xz", outputs%eps_dev(4,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_yz", outputs%eps_dev(5,:))
        end if
        if (out_variables(OUT_EPS_DEV_PL) == 1) then
            call write_1d_var_c(outputs, parent_id, "eps_dev_pl_xx", outputs%eps_dev_pl(0,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_pl_yy", outputs%eps_dev_pl(1,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_pl_zz", outputs%eps_dev_pl(2,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_pl_xy", outputs%eps_dev_pl(3,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_pl_xz", outputs%eps_dev_pl(4,:))
            call write_1d_var_c(outputs, parent_id, "eps_dev_pl_yz", outputs%eps_dev_pl(5,:))
        end if

        if (out_variables(OUT_STRESS_DEV) == 1) then
            call write_1d_var_c(outputs, parent_id, "sig_dev_xx", outputs%sig_dev(0,:))
            call write_1d_var_c(outputs, parent_id, "sig_dev_yy", outputs%sig_dev(1,:))
            call write_1d_var_c(outputs, parent_id, "sig_dev_zz", outputs%sig_dev(2,:))
            call write_1d_var_c(outputs, parent_id, "sig_dev_xy", outputs%sig_dev(3,:))
            call write_1d_var_c(outputs, parent_id, "sig_dev_xz", outputs%sig_dev(4,:))
            call write_1d_var_c(outputs, parent_id, "sig_dev_yz", outputs%sig_dev(5,:))
        end if

        if (out_variables(OUT_ENERGYD) == 1) then
            call write_1d_var_c(outputs, parent_id, "L_energy", outputs%D_energy(0,:))
            call write_1d_var_c(outputs, parent_id, "S_energy", outputs%D_energy(1,:))
            call write_1d_var_c(outputs, parent_id, "R_energy", outputs%D_energy(2,:))
        end if

        ! VELOCITY
        if (out_variables(OUT_VITESSE) == 1) then
            call write_2d_var_vecn(outputs, parent_id, "veloc", outputs%veloc)
        end if
        ! DISPLACEMENT
        if (out_variables(OUT_DEPLA) == 1) then
            call write_2d_var_vecn(outputs, parent_id, "displ", outputs%displ)
        end if

        ! ACCELERATION
        if (out_variables(OUT_ACCEL) == 1) then
            call write_2d_var_vecn(outputs, parent_id, "accel", outputs%accel)
        end if

#ifdef CPML
#ifdef DEBUG_CPML
!        write(*,*) "R1X:", maxval(outputs%R1_x), shape(outputs%R1_x), outputs%nnodes
        call write_2d_var_vecn(outputs, parent_id, "R1_0", outputs%R1_0)
        call write_2d_var_vecn(outputs, parent_id, "R1_1", outputs%R1_1)
        call write_2d_var_vecn(outputs, parent_id, "R1_2", outputs%R1_2)
        call write_2d_var_vecn(outputs, parent_id, "R2_0_dX", outputs%R2_0_dX)
        call write_2d_var_vecn(outputs, parent_id, "R2_0_dY", outputs%R2_0_dY)
        call write_2d_var_vecn(outputs, parent_id, "R2_0_dZ", outputs%R2_0_dZ)
        call write_2d_var_vecn(outputs, parent_id, "R2_1_dX", outputs%R2_1_dX)
        call write_2d_var_vecn(outputs, parent_id, "R2_1_dY", outputs%R2_1_dY)
        call write_2d_var_vecn(outputs, parent_id, "R2_1_dZ", outputs%R2_1_dZ)
        call write_2d_var_vecn(outputs, parent_id, "R2_2_dX", outputs%R2_2_dX)
        call write_2d_var_vecn(outputs, parent_id, "R2_2_dY", outputs%R2_2_dY)
        call write_2d_var_vecn(outputs, parent_id, "R2_2_dZ", outputs%R2_2_dZ)
! DEBUG_CPML
#endif

        !call write_2d_var_vecn(outputs, parent_id, "FDump", outputs%FDump)
        !call write_2d_var_vecn(outputs, parent_id, "FMasU", outputs%FMasU)
        !call write_2d_var_vecn(outputs, parent_id, "Fint" , outputs%Fint )
#endif
    end subroutine grp_write_fields

    subroutine grp_write_int_2d(outputs, parent_id, name, dim1, dim2, data, ntot_nodes)
        type (output_var_t), intent (INOUT):: outputs
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

        if (outputs%rank==0) then
            allocate(displs(0:outputs%nprocs-1))
            allocate(counts(0:outputs%nprocs-1))
        end if
        call MPI_Gather(dim2, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, outputs%comm, ierr)

        ntot_nodes = 0
        if (outputs%rank==0) then
            do n=0, outputs%nprocs-1
                displs(n) = dim1*ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
                counts(n)=counts(n)*dim1 ! for the next gatherv
            end do
            allocate(all_data(dim1,0:ntot_nodes-1))
        end if

        call MPI_Gatherv(data, dim1*dim2, MPI_INTEGER, all_data, counts, displs, &
            MPI_INTEGER, 0, outputs%comm, ierr)
        if (outputs%rank==0) then
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

    subroutine grp_write_int_1d(outputs, parent_id, name, dim1, data, ntot_nodes)
        type (output_var_t), intent (INOUT):: outputs
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

        if (outputs%rank==0) then
            allocate(displs(0:outputs%nprocs-1))
            allocate(counts(0:outputs%nprocs-1))
        end if
        call MPI_Gather(dim1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, outputs%comm, ierr)

        ntot_nodes = 0
        if (outputs%rank==0) then
            do n=0, outputs%nprocs-1
                displs(n) = ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
            end do
            allocate(all_data(0:ntot_nodes-1))
        end if


        call MPI_Gatherv(data, dim1, MPI_INTEGER, all_data, counts, displs, &
            MPI_INTEGER, 0, outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, name, H5T_STD_I32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, all_data, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
            deallocate(all_data)
            deallocate(displs)
            deallocate(counts)
        end if
    end subroutine grp_write_int_1d

    subroutine grp_write_real_1d(outputs, parent_id, name, dim1, data, ntot_nodes)
        type (output_var_t), intent (INOUT):: outputs
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim1
        real(fpp), dimension(:), intent(in) :: data
        character(len=*), INTENT(IN) :: name
        integer, intent(OUT) :: ntot_nodes
        !
        integer(HID_T) :: dset_id
        real(fpp), dimension(:), allocatable :: all_data
        integer, dimension(:), allocatable :: displs, counts
        integer :: n
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, ierr

        if (outputs%rank==0) then
            allocate(displs(0:outputs%nprocs-1))
            allocate(counts(0:outputs%nprocs-1))
        end if
        call MPI_Gather(dim1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, outputs%comm, ierr)

        ntot_nodes = 0
        if (outputs%rank==0) then
            do n=0, outputs%nprocs-1
                displs(n) = ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
            end do
            allocate(all_data(0:ntot_nodes-1))
        end if

        call MPI_Gatherv(data, dim1, MPI_REAL_FPP, all_data, counts, displs, &
            MPI_REAL_FPP, 0, outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, name, H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
            deallocate(all_data)
            deallocate(displs)
            deallocate(counts)
        end if
    end subroutine grp_write_real_1d

    subroutine compute_saved_elements(Tdomain, outputs)
        type (domain), intent (INOUT):: Tdomain
        type (output_var_t), intent(inout) :: outputs
        integer :: n, i, j, k, ngll, ig, gn, ne
        !
        integer :: domain_type, imat
        ! on sauvegarde ncells a la place de le calculer
        allocate(outputs%irenum(0:Tdomain%n_glob_points-1))

        outputs%irenum = -1
        outputs%ncells = 0
        ig = 0
        ne = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            ne = ne + 1
            do k = 0,ngll - 1
                do j = 0,ngll - 1
                    do i = 0,ngll - 1
                        gn = Tdomain%specel(n)%Iglobnum(i,j,k)
                        if (outputs%irenum(gn) == -1) then
                            outputs%irenum(gn) = ig
                            ig = ig + 1
                        end if
                    end do
                end do
            end do
            outputs%ncells = outputs%ncells + (ngll-1)**3
        end do
        outputs%nelems = ne
        outputs%nnodes = ig

        allocate(outputs%domains(0:ig-1))
        domain_type = 0
        outputs%domains = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            imat = Tdomain%specel(n)%mat_index
            domain_type = Tdomain%specel(n)%domain

            do k = 0,ngll - 1
                do j = 0,ngll - 1
                    do i = 0,ngll - 1
                        gn = Tdomain%specel(n)%Iglobnum(i,j,k)
                        if (domain_type>outputs%domains(outputs%irenum(gn))) then
                            outputs%domains(outputs%irenum(gn)) = domain_type
                        endif
                    end do
                end do
            end do
        end do

        call compute_offsets(outputs)

    end subroutine compute_saved_elements

    subroutine compute_offsets(outputs)
        type(output_var_t), intent(inout) :: outputs
        !
        integer :: nnodes, ncells
        integer :: ntot_nodes, ntot_cells
        integer :: ierr, n
        nnodes = outputs%nnodes
        ncells = outputs%ncells
        ! Needed on all proc to compute displacement offset for node numbering
        allocate(outputs%displs_n(0:outputs%nprocs-1))
        allocate(outputs%counts_n(0:outputs%nprocs-1))
        if (outputs%rank==0) then
            allocate(outputs%displs_c(0:outputs%nprocs-1))
            allocate(outputs%counts_c(0:outputs%nprocs-1))
            allocate(outputs%displs2d_n(0:outputs%nprocs-1))
            allocate(outputs%counts2d_n(0:outputs%nprocs-1))
        end if
        call MPI_AllGather(nnodes, 1, MPI_INTEGER, outputs%counts_n, 1, &
            MPI_INTEGER, outputs%comm, ierr)
        call MPI_Gather(ncells, 1, MPI_INTEGER, outputs%counts_c, 1, &
            MPI_INTEGER, 0, outputs%comm, ierr)
        ! 1D FIELDS-NOEUD
        ntot_nodes = 0
        ntot_cells = 0
        do n=0, outputs%nprocs-1
            outputs%displs_n(n)   = ntot_nodes
            ntot_nodes    = ntot_nodes + outputs%counts_n(n)
        end do
        if (outputs%rank==0) then
            do n=0, outputs%nprocs-1
                outputs%displs_c(n)   = ntot_cells
                ntot_cells    = ntot_cells + outputs%counts_c(n)
            end do
            outputs%displs2d_n = 3*outputs%displs_n
            outputs%counts2d_n = 3*outputs%counts_n
            allocate(outputs%all_data_1d_n(0:ntot_nodes-1))
            allocate(outputs%all_data_1d_c(0:ntot_cells-1))
            allocate(outputs%all_data_2d_n(0:2,0:ntot_nodes-1))
            allocate(outputs%all_data_2d_c(0:2,0:ntot_cells-1))
            outputs%ntot_nodes = ntot_nodes
            outputs%ntot_cells = ntot_cells
        else
            outputs%ntot_nodes = -1
            outputs%ntot_cells = -1
        end if
    end subroutine compute_offsets


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
    subroutine write_snapshot_geom(Tdomain, outputs)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        type (output_var_t), intent(INOUT) :: outputs
        integer :: rg
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        integer :: hdferr, code, ierr
        integer :: group, subgroup
        integer, dimension(0:Tdomain%nb_procs-1) :: nodes_per_proc
        !- Create subdomains communicators
        rg = Tdomain%rank
        group = rg/Tdomain%ngroup
        subgroup = mod(rg, Tdomain%ngroup)
        call MPI_Comm_split(Tdomain%communicateur, group, subgroup, outputs%comm, ierr)
        call MPI_Comm_size(outputs%comm, outputs%nprocs,  code)
        call MPI_Comm_rank(outputs%comm, outputs%rank, code)
        outputs%group = group

        call init_hdf5()
        if (rg==0) then
            ierr = sem_mkdir(trim(adjustl(path_results)))
        end if
        call MPI_Barrier(Tdomain%communicateur, code)
        call semname_snap_geom_file(group, fnamef)

        if (outputs%rank==0) then
            call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        else
            fid = -1
        endif

        call mpi_barrier(Tdomain%communicateur, code)
        
        call compute_saved_elements(Tdomain, outputs)
        
        call write_global_nodes(Tdomain, fid, outputs)
        
        call write_elem_connectivity(Tdomain, fid, outputs)
        
        call allocate_fields(outputs, Tdomain%out_var_snap, Tdomain%nl_flag)
        
        call write_constant_fields(Tdomain, fid, outputs)
        
        if (outputs%rank==0) then
            call h5fclose_f(fid, hdferr)
        endif

        call mpi_gather(outputs%nnodes, 1, MPI_INTEGER, nodes_per_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (rg==0) call write_master_xdmf(Tdomain, nodes_per_proc)
    end subroutine write_snapshot_geom

    subroutine write_global_nodes(Tdomain, fid, outputs)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        type(output_var_t), intent(inout) :: outputs
        !
        real(fpp), dimension(:,:), allocatable :: nodes
        integer :: n
        !
        allocate(nodes(0:2,0:outputs%nnodes-1))
        do n = 0, Tdomain%n_glob_points-1
            if (outputs%irenum(n)>=0) then
                nodes(:,outputs%irenum(n)) = Tdomain%GlobCoord(:,n)
            end if
        end do

        call write_2d_var_vecn(outputs, fid, "Nodes", nodes)
!        call grp_write_real_2d(outputs, fid, "Nodes", 3, outputs%nnodes, nodes, nnodes_tot)
        deallocate(nodes)
    end subroutine write_global_nodes

    subroutine write_elem_connectivity(Tdomain, fid, outputs)
        implicit none
        type (domain), intent (INOUT) :: Tdomain
        type (output_var_t), intent(INOUT) :: outputs
        integer(HID_T), intent(in) :: fid
        !
        integer :: ngll
        integer(HSIZE_T), dimension(2) :: dims
        integer, dimension(:,:), allocatable :: data
        integer, dimension(:), allocatable :: mat, proc
        integer :: count, ig, nglobnum
        integer :: i, j, k, n, nb_elem
        integer :: nb_elem_tot
        integer :: ioffset

        ! First we count the number of hexaedrons
        count = 0
        nglobnum = 0
        k = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            ! Max number of global points (can count elem vertices twice)
            nglobnum = nglobnum + ngll**3
            ! Number of subelements
            count = count+(ngll-1)**3
            k = k + 1
        enddo
        nb_elem = k

        allocate( data(1:8,0:count-1))
        allocate( mat(0:count-1))
        allocate( proc(0:count-1))
        dims(1) = 8
        dims(2) = count
        count = 0
        ig = 0
        ioffset = outputs%displs_n(outputs%rank)
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            do k = 0,ngll-2
                do j = 0,ngll-2
                    do i = 0,ngll-2
                        data(1,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+0,j+0,k+0))
                        data(2,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+1,j+0,k+0))
                        data(3,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+1,j+1,k+0))
                        data(4,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+0,j+1,k+0))
                        data(5,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+0,j+0,k+1))
                        data(6,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+1,j+0,k+1))
                        data(7,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+1,j+1,k+1))
                        data(8,count) = ioffset+outputs%irenum(Tdomain%specel(n)%Iglobnum(i+0,j+1,k+1))
                        mat(count) = Tdomain%specel(n)%mat_index
                        proc(count) = Tdomain%rank
                        count=count+1
                    end do
                end do
            end do
        end do

        call grp_write_int_2d(outputs, fid, "Elements", 8, count, data, nb_elem_tot)
        Tdomain%n_hexa = nb_elem_tot
        call grp_write_int_1d(outputs, fid, "Material", count, mat, nb_elem_tot)
        call grp_write_int_1d(outputs, fid, "Proc", count, proc, nb_elem_tot)

        deallocate(mat)
        deallocate(proc)
        deallocate(data)
    end subroutine write_elem_connectivity

    subroutine allocate_fields(outputs, out_flags, nl_flag)
        type (output_var_t), intent(inout) :: outputs
        integer, dimension(0:OUT_LAST), intent(in) :: out_flags
        logical, intent(in) :: nl_flag
        !
        logical :: flag_gradU
        integer :: nnodes,ncells
        nnodes = outputs%nnodes
        ncells = outputs%ncells

        if (nl_flag) then
            flag_gradU = (out_flags(OUT_ENERGYP)     + &
                          out_flags(OUT_ENERGYK)     + &
                          out_flags(OUT_DUDX)        + &
                          out_flags(OUT_EPS_VOL)) /= 0
        else
            flag_gradU = (out_flags(OUT_PRESSION)    + &
                          out_flags(OUT_ENERGYP)     + &
                          out_flags(OUT_ENERGYK)     + &
                          out_flags(OUT_DUDX)        + &
                          out_flags(OUT_EPS_VOL)     + &
                          out_flags(OUT_EPS_DEV)     + &
                          out_flags(OUT_STRESS_DEV)  + &
                          out_flags(OUT_ENERGYD)   )/= 0
        endif
        ! sortie par noeud
        if (out_flags(OUT_DEPLA     ) == 1) allocate(outputs%displ(0:2,0:nnodes-1))
        if (out_flags(OUT_VITESSE   ) == 1) allocate(outputs%veloc(0:2,0:nnodes-1))
        if (out_flags(OUT_ACCEL     ) == 1) allocate(outputs%accel(0:2,0:nnodes-1))
        if (out_flags(OUT_PRESSION  ) == 1) allocate(outputs%press_n(0:nnodes-1))
        if (out_flags(OUT_DUDX      ) == 1) allocate(outputs%dUdX(0:8,0:nnodes-1))
        ! sortie par element
        if (out_flags(OUT_ENERGYP   ) == 1) allocate(outputs%P_energy(0:ncells-1))
        if (out_flags(OUT_ENERGYK   ) == 1) allocate(outputs%K_energy(0:ncells-1))
        if (out_flags(OUT_ENERGYD   ) == 1) allocate(outputs%D_energy(0:2,0:ncells-1))
        if (out_flags(OUT_EPS_VOL   ) == 1) allocate(outputs%eps_vol(0:ncells-1))
        if (out_flags(OUT_PRESSION  ) == 1) allocate(outputs%press_c(0:ncells-1))
        if (out_flags(OUT_EPS_DEV   ) == 1) allocate(outputs%eps_dev(0:5,0:ncells-1))
        if (out_flags(OUT_STRESS_DEV) == 1) allocate(outputs%sig_dev(0:5,0:ncells-1))
        if (out_flags(OUT_EPS_DEV_PL) == 1) allocate(outputs%eps_dev_pl(0:5,0:ncells-1))
        ! initialize
        if (out_flags(OUT_DEPLA     ) == 1) outputs%displ      = 0.
        if (out_flags(OUT_VITESSE   ) == 1) outputs%veloc      = 0.
        if (out_flags(OUT_ACCEL     ) == 1) outputs%accel      = 0.
        if (out_flags(OUT_ENERGYP   ) == 1) outputs%P_energy   = 0.
        if (out_flags(OUT_ENERGYK   ) == 1) outputs%K_energy   = 0.
        if (out_flags(OUT_DUDX      ) == 1) outputs%dUdX       = 0.
        if (out_flags(OUT_EPS_VOL   ) == 1) outputs%eps_vol    = 0.
        if (out_flags(OUT_PRESSION  ) == 1) outputs%press_c    = 0.
        if (out_flags(OUT_PRESSION  ) == 1) outputs%press_n    = 0.
        if (out_flags(OUT_EPS_DEV   ) == 1) outputs%eps_dev    = 0.
        if (out_flags(OUT_STRESS_DEV) == 1) outputs%sig_dev    = 0.
        if (out_flags(OUT_EPS_DEV_PL) == 1) outputs%eps_dev_pl = 0.
        if (out_flags(OUT_ENERGYD   ) == 1) outputs%D_energy   = 0.
        !
#ifdef CPML
        allocate(outputs%R1_0   (0:2,0:(nnodes-1)))
        allocate(outputs%R1_1   (0:2,0:(nnodes-1)))
        allocate(outputs%R1_2   (0:2,0:(nnodes-1)))
        allocate(outputs%R2_0_dX(0:2,0:(nnodes-1)))
        allocate(outputs%R2_0_dY(0:2,0:(nnodes-1)))
        allocate(outputs%R2_0_dZ(0:2,0:(nnodes-1)))
        allocate(outputs%R2_1_dX(0:2,0:(nnodes-1)))
        allocate(outputs%R2_1_dY(0:2,0:(nnodes-1)))
        allocate(outputs%R2_1_dZ(0:2,0:(nnodes-1)))
        allocate(outputs%R2_2_dX(0:2,0:(nnodes-1)))
        allocate(outputs%R2_2_dY(0:2,0:(nnodes-1)))
        allocate(outputs%R2_2_dZ(0:2,0:(nnodes-1)))
        outputs%R1_0 = 0.0
        outputs%R1_1 = 0.0
        outputs%R1_2 = 0.0
        outputs%R2_0_dX = 0.0
        outputs%R2_0_dY = 0.0
        outputs%R2_0_dZ = 0.0
        outputs%R2_1_dX = 0.0
        outputs%R2_1_dY = 0.0
        outputs%R2_1_dZ = 0.0
        outputs%R2_2_dX = 0.0
        outputs%R2_2_dY = 0.0
        outputs%R2_2_dZ = 0.0
        allocate(outputs%FDump(0:2,0:(nnodes-1)))
        allocate(outputs%FMasU(0:2,0:(nnodes-1)))
        allocate(outputs%Fint (0:2,0:(nnodes-1)))
        outputs%FDump = 0.
        outputs%FMasU = 0.
        outputs%Fint  = 0.
#endif
        return
        !
    end subroutine allocate_fields

    subroutine deallocate_fields(out_flags, fields)
        integer, dimension(0:OUT_LAST), intent(in) :: out_flags
        type (output_var_t), intent(inout) :: fields

        if (out_flags(OUT_DEPLA     ) == 1) deallocate(fields%displ)
        if (out_flags(OUT_VITESSE   ) == 1) deallocate(fields%veloc)
        if (out_flags(OUT_ACCEL     ) == 1) deallocate(fields%accel)
        if (out_flags(OUT_ENERGYP   ) == 1) deallocate(fields%P_energy)
        if (out_flags(OUT_ENERGYK   ) == 1) deallocate(fields%K_energy)
        if (out_flags(OUT_DUDX      ) == 1) deallocate(fields%dUdX)
        if (out_flags(OUT_EPS_VOL   ) == 1) deallocate(fields%eps_vol)
        if (out_flags(OUT_PRESSION  ) == 1) deallocate(fields%press_c)
        if (out_flags(OUT_PRESSION  ) == 1) deallocate(fields%press_n)
        if (out_flags(OUT_EPS_DEV   ) == 1) deallocate(fields%eps_dev)
        if (out_flags(OUT_STRESS_DEV) == 1) deallocate(fields%sig_dev)
        if (out_flags(OUT_EPS_DEV_PL) == 1) deallocate(fields%eps_dev_pl)
        if (out_flags(OUT_ENERGYD   ) == 1) deallocate(fields%D_energy)
        !
#ifdef CPML
        deallocate(fields%R1_0)
        deallocate(fields%R1_1)
        deallocate(fields%R1_2)
        deallocate(fields%R2_0_dX)
        deallocate(fields%R2_0_dY)
        deallocate(fields%R2_0_dZ)
        deallocate(fields%R2_1_dX)
        deallocate(fields%R2_1_dY)
        deallocate(fields%R2_1_dZ)
        deallocate(fields%R2_2_dX)
        deallocate(fields%R2_2_dY)
        deallocate(fields%R2_2_dZ)
        deallocate(fields%FDump)
        deallocate(fields%FMasU)
        deallocate(fields%Fint )
#endif
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


    subroutine save_field_h5(Tdomain, isort, outputs)
        use sdomain
        use dom_solid
        use dom_solid_dg
        use dom_fluid
        use dom_fluid_aniso
        use dom_solidpml
        use dom_fluidpml

        implicit none
        type (domain), intent (INOUT):: Tdomain
        type(output_var_t), intent(INOUT) :: outputs
        integer, intent(in) :: isort
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        integer :: domain_type
        integer, dimension(:), allocatable :: valence
        integer :: hdferr
        integer :: ngll, oldngll
        integer :: i, j, k, n, nn, m, ind
        integer :: nnodes
        type(Element), pointer :: el
        type(subdomain), pointer :: sub_dom_mat
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, K_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: D_energy
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev, eps_dev_pl, dUdX
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        integer :: bnum, ee
        real(fpp), dimension(:), allocatable :: GLLc ! GLLw
        real(fpp), dimension(:,:,:), allocatable :: jac

        integer, dimension(0:OUT_LAST) :: out_variables
        logical :: nl_flag

        integer :: cell_start

        nl_flag=Tdomain%nl_flag
        nnodes = outputs%nnodes
        out_variables(:) = Tdomain%out_var_snap(:)
        call create_dir_sorties(Tdomain, isort)
        allocate(valence(0:nnodes-1))

        valence(:) = 0
        oldngll = 0

        cell_start = 0
#ifdef CPML
        outputs%R1_0 = 0.
        outputs%R1_1 = 0.
        outputs%R1_2 = 0.
        outputs%R2_0_dX = 0.
        outputs%R2_1_dX = 0.
        outputs%R2_2_dX = 0.
        outputs%R2_0_dY = 0.
        outputs%R2_1_dY = 0.
        outputs%R2_2_dY = 0.
        outputs%R2_0_dZ = 0.
        outputs%R2_1_dZ = 0.
        outputs%R2_2_dZ = 0.
#endif
        do n = 0,Tdomain%n_elem-1
            el => Tdomain%specel(n)
            sub_dom_mat => Tdomain%sSubdomain(el%mat_index)
            if (.not. el%OUTPUT ) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            if (ngll/=oldngll .or. oldngll==0) then
                ! Allocate everything here, it simpler and not that costly...
                if (oldngll/=0) then
                    deallocate(fieldP,fieldU,fieldV,fieldA)
                    deallocate(eps_vol,eps_dev,sig_dev,dUdX)
                    deallocate(P_energy,K_energy,D_energy)
                endif
                allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                allocate(dUdX(0:ngll-1,0:ngll-1,0:ngll-1,0:8))
                allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(K_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(D_energy(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                allocate(eps_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                ! tot energy 5
                allocate(eps_dev_pl(0:ngll-1,0:ngll-1,0:ngll-1,0:6))
                allocate(sig_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                oldngll = ngll
            endif
            domain_type = Tdomain%specel(n)%domain
            select case(domain_type)
                case (DM_SOLID_CG)
!                    if (Tdomain%sdom%PlaneW%Exist) then
!                        call compute_planeW_Exafield(el%lnum,Tdomain%TimeD%rtime,Tdomain,0)
!                    end if
                    call get_solid_dom_var(Tdomain%sdom, el%lnum, out_variables,    &
                        fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol,&
                        eps_dev, sig_dev, dUdX, nl_flag, eps_dev_pl)
                case (DM_SOLID_DG)
                    call get_solid_dg_dom_var(Tdomain%sdomdg, el%lnum, out_variables, &
                        fieldV, fieldA, eps_vol, eps_dev)
                case (DM_FLUID_CG)
                    call get_fluid_dom_var(Tdomain%fdom, el%lnum, out_variables,        &
                        fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol, eps_dev, &
                        sig_dev, dUdX)
                case (DM_FLUID_CG_ANISO)
                    call get_fluid_aniso_dom_var(Tdomain%fanisodom, el%lnum, out_variables, &
                        fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol, eps_dev, &
                        sig_dev, dUdX)
                case (DM_SOLID_CG_PML)
                    call get_solidpml_dom_var(Tdomain%spmldom, el%lnum, out_variables,           &
                        fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol, eps_dev, sig_dev)
#ifdef CPML
                    call get_solidpml_rfields(Tdomain, Tdomain%spmldom, n, el%lnum, outputs)
#endif
                case (DM_FLUID_CG_PML)
                    call get_fluidpml_dom_var(Tdomain%fpmldom, el%lnum, out_variables,           &
                    fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol, eps_dev, sig_dev)
#ifdef CPML
                    call get_fluidpml_rfields(Tdomain, Tdomain%fpmldom, n, el%lnum, outputs)
#endif
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
                        ind = outputs%irenum(el%Iglobnum(i,j,k))
                        valence(ind) = valence(ind)+1
                        select case (domain_type)
                            case (DM_SOLID_CG)
                                jac(i,j,k) = Tdomain%sdom%Jacob_   (i,j,k,bnum,ee)
                            case (DM_SOLID_DG)
                                jac(i,j,k) = Tdomain%sdomdg%Jacob_ (i,j,k,bnum,ee)
                            case (DM_SOLID_CG_PML)
                                jac(i,j,k) = Tdomain%spmldom%Jacob_(i,j,k,bnum,ee)
                            case (DM_FLUID_CG)
                                jac(i,j,k) = Tdomain%fdom%Jacob_      (i,j,k,bnum,ee)
                            case (DM_FLUID_CG_ANISO)
                                jac(i,j,k) = Tdomain%fanisodom%Jacob_ (i,j,k,bnum,ee)
                            case (DM_FLUID_CG_PML)
                                jac(i,j,k) = Tdomain%fpmldom%Jacob_(i,j,k,bnum,ee)
                            case default
                                stop "unknown domain"
                        end select

                        ! sortie par noeud
                        if (out_variables(OUT_DEPLA)    == 1) outputs%displ(0:2,ind) = fieldU(i,j,k,0:2)
                        if (out_variables(OUT_VITESSE)  == 1) outputs%veloc(0:2,ind) = fieldV(i,j,k,0:2)
                        if (out_variables(OUT_ACCEL)    == 1) outputs%accel(0:2,ind) = fieldA(i,j,k,0:2)
                        if (out_variables(OUT_PRESSION) == 1) outputs%press_n(ind)   = fieldP(i,j,k)
                        if (out_variables(OUT_DUDX) == 1) outputs%dUdX(0:8,ind) = dUdX(i,j,k,0:8)
!                       ajout dUdX
                        !if (domain_type==DM_SOLID_DG) then
!                            if (out_variables(OUT_EPS_VOL)    == 1) outputs%eps_vol(ind) = eps_vol(i,j,k)
!                            if (out_variables(OUT_EPS_DEV)    == 1) outputs%eps_dev(0:5,ind) = eps_dev(i,j,k,0:5)
!                        end if
                    enddo
                enddo
            enddo
            ! sortie de la valeur au centre par sub-element
            !call domain_gllw(Tdomain, domain_type, GLLw)
            call domain_gllc(Tdomain, Tdomain%specel(n)%domain, GLLc)

            if (out_variables(OUT_PRESSION) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, fieldP, outputs%press_c)
            endif
            if (out_variables(OUT_ENERGYP) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, P_energy, outputs%P_energy)
            endif
            if (out_variables(OUT_ENERGYK) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, K_energy, outputs%K_energy)
            endif
            if (out_variables(OUT_EPS_VOL) == 1) then !  .and. domain_type/=DM_SOLID_DG
                call evaluate_cell_centers(ngll, GLLc, cell_start, eps_vol, outputs%eps_vol)
            endif
            if (out_variables(OUT_EPS_DEV) == 1) then !  .and. domain_type/=DM_SOLID_DG
                do m = 0,OUT_VAR_DIMS_3D(OUT_EPS_DEV)-1
                    call evaluate_cell_centers(ngll, GLLc, cell_start, eps_dev(:,:,:,m), outputs%eps_dev(m,:))
                end do
            endif
            if (out_variables(OUT_EPS_DEV_PL) == 1) then
                do m = 0,OUT_VAR_DIMS_3D(OUT_EPS_DEV_PL)-1
                    call evaluate_cell_centers(ngll, GLLc, cell_start, eps_dev_pl(:,:,:,m), outputs%eps_dev_pl(m,:))
                end do
            endif
            if (out_variables(OUT_STRESS_DEV) == 1) then
                do m = 0,OUT_VAR_DIMS_3D(OUT_STRESS_DEV)-1
                    call evaluate_cell_centers(ngll, GLLc, cell_start, sig_dev(:,:,:,m), outputs%sig_dev(m,:))
                end do
            endif
            if (out_variables(OUT_ENERGYD) == 1) then
                do m = 0,OUT_VAR_DIMS_3D(OUT_ENERGYD)-1
                    call evaluate_cell_centers(ngll, GLLc, cell_start, D_energy(:,:,:,m), outputs%D_energy(m,:))
                end do
            endif
            
            cell_start = cell_start + (ngll-1)**3
            if(allocated(jac)) deallocate(jac)
            if(allocated(GLLc)) deallocate(GLLc)
            !if(allocated(GLLw)) deallocate(GLLw)
        enddo
        if (oldngll/=0) then
            deallocate(fieldP,fieldU,fieldV,fieldA)
            deallocate(eps_vol,eps_dev,sig_dev,dUdX)
            deallocate(P_energy,K_energy,D_energy)
        endif

        if (outputs%rank==0) then
            call semname_snap_result_file(outputs%group, isort, fnamef)
            call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        else
            fid = -1
        endif
        call grp_write_fields(outputs, fid, out_variables)
        if (outputs%rank==0) then
            call h5fclose_f(fid, hdferr)
            call write_xdmf(Tdomain, isort, outputs, out_variables)
        endif

        deallocate(valence)

        call mpi_barrier(Tdomain%communicateur, hdferr)

    end subroutine save_field_h5

    subroutine write_master_xdmf(Tdomain, nodes_per_proc)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: n_procs, nelem, n_groups
        character (len=MAX_FILE_SIZE) :: fnamef
        integer, intent(in), dimension(0:Tdomain%nb_procs-1) :: nodes_per_proc
        integer :: group, k, group_nodes
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

        do group=0,n_groups-1
            ! check if the group does any output
            group_nodes = 0
            do k=group*Tdomain%ngroup, min(n_procs-1, (group+1)*Tdomain%ngroup-1)
                group_nodes = group_nodes + nodes_per_proc(k)
            end do
            if (group_nodes /= 0) then
                write(61,"(a,I4.4,a)") '<xi:include href="mesh.',group,'.xmf" xpointer="xpointer(//Xdmf/Domain/Grid)"/>'
            end if

        end do

        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
        close(61)

    end subroutine write_master_xdmf

    subroutine write_xdmf_attr_scalar_nodes(aname, nn, i, group, adata)
        character(len=*) :: aname, adata
        integer :: nn, i, group
        !
        write(61,"(a,a,a,a,a,a)") '<Attribute Name="', trim(aname), '" Center="Node" AttributeType="Scalar">'
        write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="', nn, '">'
        write(61,"(a,I4.4,a,I4.4,a,a)") 'Rsem', i, '/sem_field.', group, '.h5:/', adata
        write(61,"(a)") '</DataItem>'
        write(61,"(a)") '</Attribute>'
    end subroutine write_xdmf_attr_scalar_nodes

    subroutine write_xdmf_attr_vector_nodes(aname, nn, i, group, adata)
        character(len=*) :: aname, adata
        integer :: nn, i, group
        !
        write(61,"(a,a,a,a,a,a)") '<Attribute Name="', trim(aname), '" Center="Node" AttributeType="Vector">'
        write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="', nn,' 3">'
        write(61,"(a,I4.4,a,I4.4,a,a)") 'Rsem', i, '/sem_field.', group, '.h5:/', adata
        write(61,"(a)") '</DataItem>'
        write(61,"(a)") '</Attribute>'
    end subroutine write_xdmf_attr_vector_nodes

    subroutine write_xdmf_attr_scalar_cells(aname, ne, i, group, adata)
        character(len=*) :: aname, adata
        integer :: ne, i, group
        !
        write(61,"(a,a,a,a,a,a)") '<Attribute Name="', trim(aname), '" Center="Cell" AttributeType="Scalar">'
        write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="', ne, '">'
        write(61,"(a,I4.4,a,I4.4,a,a)") 'Rsem', i, '/sem_field.', group, '.h5:/', adata
        write(61,"(a)") '</DataItem>'
        write(61,"(a)") '</Attribute>'
    end subroutine write_xdmf_attr_scalar_cells

    subroutine write_xdmf_attr_vector_cells(aname, ne, i, group, adata)
        character(len=*) :: aname, adata
        integer :: ne, i, group
        !
        write(61,"(a,a,a,a,a,a)") '<Attribute Name="', trim(aname), '" Center="Cell" AttributeType="Vector">'
        write(61,"(a,I9,a)") '<DataItem Format="HDF" NumberType="Float" Precision="4" Dimensions="', ne, '">'
        write(61,"(a,I4.4,a,I4.4,a,a)") 'Rsem', i, '/sem_field.', group, '.h5:/', adata
        write(61,"(a)") '</DataItem>'
        write(61,"(a)") '</Attribute>'
    end subroutine write_xdmf_attr_vector_cells

    subroutine write_xdmf(Tdomain, isort, outputs, out_variables)
        implicit none
        type (domain), intent (IN):: Tdomain
        type(output_var_t), intent(in) :: outputs
        integer, intent(in) :: isort
        integer, dimension(0:), intent(in) :: out_variables
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: i, j, nn, ne, group
        real(fpp) :: time
        character(len=11), dimension(0:8) :: R2label, R2data

        nn = outputs%ntot_nodes
        ne = outputs%ntot_cells
        group = outputs%group
        call semname_xdmf(group, fnamef)

        R2label( 0) = "R2o0odX"; R2data( 0) = "R2_0_dX"
        R2label( 1) = "R2o0odY"; R2data( 1) = "R2_0_dY"
        R2label( 2) = "R2o0odZ"; R2data( 2) = "R2_0_dZ"
        R2label( 3) = "R2o1odX"; R2data( 3) = "R2_1_dX"
        R2label( 4) = "R2o1odY"; R2data( 4) = "R2_1_dY"
        R2label( 5) = "R2o1odZ"; R2data( 5) = "R2_1_dZ"
        R2label( 6) = "R2o2odX"; R2data( 6) = "R2_2_dX"
        R2label( 7) = "R2o2odY"; R2data( 7) = "R2_2_dY"
        R2label( 8) = "R2o2odZ"; R2data( 8) = "R2_2_dZ"

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

            if (out_variables(OUT_DEPLA)    == 1) call write_xdmf_attr_vector_nodes("Displ",      nn, i, group, "displ"     )
            if (out_variables(OUT_VITESSE)  == 1) call write_xdmf_attr_vector_nodes("Veloc",      nn, i, group, "veloc"     )
            if (out_variables(OUT_ACCEL)    == 1) call write_xdmf_attr_vector_nodes("Accel",      nn, i, group, "accel"     )
            if (out_variables(OUT_PRESSION) == 1) call write_xdmf_attr_scalar_nodes("Press_gll",  nn, i, group, "press_gll" )
            if (out_variables(OUT_PRESSION) == 1) call write_xdmf_attr_scalar_cells("Press_elem", ne, i, group, "press_elem")
            if (out_variables(OUT_EPS_VOL)  == 1) call write_xdmf_attr_scalar_cells("eps_vol",    ne, i, group, "eps_vol"   )
            !plus propre de faire passer DUDX en attribut vecteur ou equiv (mais
            !write_attr_vector pas adapté, dim forcé à 3?)
            if (out_variables(OUT_DUDX) == 1) then
                call write_xdmf_attr_scalar_nodes("dUxdx", nn, i, group, "dUxdx")
                call write_xdmf_attr_scalar_nodes("dUxdy", nn, i, group, "dUxdy")
                call write_xdmf_attr_scalar_nodes("dUxdz", nn, i, group, "dUxdz")
                call write_xdmf_attr_scalar_nodes("dUydx", nn, i, group, "dUydx")
                call write_xdmf_attr_scalar_nodes("dUydy", nn, i, group, "dUydy")
                call write_xdmf_attr_scalar_nodes("dUydz", nn, i, group, "dUydz")
                call write_xdmf_attr_scalar_nodes("dUzdx", nn, i, group, "dUzdx")
                call write_xdmf_attr_scalar_nodes("dUzdy", nn, i, group, "dUzdy")
                call write_xdmf_attr_scalar_nodes("dUzdz", nn, i, group, "dUzdz")
            end if
            if (out_variables(OUT_EPS_DEV) == 1) then
                call write_xdmf_attr_scalar_cells("eps_dev_xx", ne, i, group, "eps_dev_xx")
                call write_xdmf_attr_scalar_cells("eps_dev_yy", ne, i, group, "eps_dev_yy")
                call write_xdmf_attr_scalar_cells("eps_dev_zz", ne, i, group, "eps_dev_zz")
                call write_xdmf_attr_scalar_cells("eps_dev_xy", ne, i, group, "eps_dev_xy")
                call write_xdmf_attr_scalar_cells("eps_dev_xz", ne, i, group, "eps_dev_xz")
                call write_xdmf_attr_scalar_cells("eps_dev_yz", ne, i, group, "eps_dev_yz")
                if (Tdomain%nl_flag) then
                    call write_xdmf_attr_scalar_cells("eps_dev_pl_xx", ne, i, group, "eps_dev_pl_xx")
                    call write_xdmf_attr_scalar_cells("eps_dev_pl_yy", ne, i, group, "eps_dev_pl_yy")
                    call write_xdmf_attr_scalar_cells("eps_dev_pl_zz", ne, i, group, "eps_dev_pl_zz")
                    call write_xdmf_attr_scalar_cells("eps_dev_pl_xy", ne, i, group, "eps_dev_pl_xy")
                    call write_xdmf_attr_scalar_cells("eps_dev_pl_xz", ne, i, group, "eps_dev_pl_xz")
                    call write_xdmf_attr_scalar_cells("eps_dev_pl_yz", ne, i, group, "eps_dev_pl_yz")
                end if
            end if
            if (out_variables(OUT_STRESS_DEV) == 1) then
                call write_xdmf_attr_scalar_cells("sig_dev_xx", ne, i, group, "sig_dev_xx")
                call write_xdmf_attr_scalar_cells("sig_dev_yy", ne, i, group, "sig_dev_yy")
                call write_xdmf_attr_scalar_cells("sig_dev_zz", ne, i, group, "sig_dev_zz")
                call write_xdmf_attr_scalar_cells("sig_dev_xy", ne, i, group, "sig_dev_xy")
                call write_xdmf_attr_scalar_cells("sig_dev_xz", ne, i, group, "sig_dev_xz")
                call write_xdmf_attr_scalar_cells("sig_dev_yz", ne, i, group, "sig_dev_yz")
            end if
            if (out_variables(OUT_ENERGYP) == 1) call write_xdmf_attr_scalar_cells("P_energy", ne, i, group, "P_energy")
            if (out_variables(OUT_ENERGYK) == 1) call write_xdmf_attr_scalar_cells("K_energy", ne, i, group, "K_energy")
            if (out_variables(OUT_ENERGYD) == 1) then
                call write_xdmf_attr_scalar_cells("L_energy", ne, i, group, "L_energy")
                call write_xdmf_attr_scalar_cells("S_energy", ne, i, group, "S_energy")
                call write_xdmf_attr_scalar_cells("R_energy", ne, i, group, "R_energy")
            end if
            ! DOMAIN
            write(61,"(a)") '<Attribute Name="Domain" Center="Grid" AttributeType="Scalar">'
            write(61,"(a,I4,a)") '<DataItem Format="XML" NumberType="Int"  Dimensions="1">',group,'</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! MAT
!            call write_xdmf_attr_iscalar_cells_cst("Mat", ne, group, "Material")
            write(61,"(a)") '<Attribute Name="Mat" Center="Cell" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mat" Format="HDF" NumberType="Int"  Dimensions="',ne, &
                '">geometry',group,'.h5:/Material</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! PROC
!            call write_xdmf_attr_iscalar_cells_cst("Proc", ne, group, "Proc")
            write(61,"(a)") '<Attribute Name="Proc" Center="Cell" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Proc" Format="HDF" NumberType="Int"  Dimensions="',ne, &
                '">geometry',group,'.h5:/Proc</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! MASS
            write(61,"(a)") '<Attribute Name="Mass" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" NumberType="Float" Precision="4"  Dimensions="',nn, &
                '">geometry',group,'.h5:/Mass</DataItem>'
            write(61,"(a)") '</Attribute>'
#ifdef CPML
#ifdef DEBUG_CPML
            write(61,"(a)") '<Attribute Name="Alpha_PML" Center="Node" AttributeType="Vector">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Alpha_PML" Format="HDF" NumberType="Float" Precision="4" Dimensions="3 ',nn, &
               '">geometry',group,'.h5:/Alpha_PML</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Kappa_PML" Center="Node" AttributeType="Vector">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Kappa_PML" Format="HDF" NumberType="Float" Precision="4" Dimensions="3 ',nn, &
                '">geometry',group,'.h5:/Kappa_PML</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Dxi_K_PML" Center="Node" AttributeType="Vector">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Dxi_K_PML" Format="HDF" NumberType="Float" Precision="4" Dimensions="3 ',nn, &
                '">geometry',group,'.h5:/Dxi_K_PML</DataItem>'
            write(61,"(a)") '</Attribute>'

            call write_xdmf_attr_vector_nodes("R1o0", nn, i, group, "R1_0")
            call write_xdmf_attr_vector_nodes("R1o1", nn, i, group, "R1_1")
            call write_xdmf_attr_vector_nodes("R1o2", nn, i, group, "R1_2")

            do j = 0, 8
                call write_xdmf_attr_vector_nodes(trim(R2label(j)), nn, i, group, trim(R2data(j)))
            end do
#endif
!
#else
            ! ALPHA/DUMPSX
            write(61,"(a)") '<Attribute Name="Alpha" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Alpha" Format="HDF" NumberType="Float" Precision="4"  Dimensions="',nn, &
                '">geometry',group,'.h5:/Alpha</DataItem>'
            write(61,"(a)") '</Attribute>'
#endif
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
            if (Tdomain%out_var_snap(OUT_GRAD_LA) == 1) then
                ! GRAD LAMBDA
                write(61,"(a)") '<Attribute Name="GradLa_gll" Center="Node" AttributeType="Vector">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="GradLa_gll" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                    ' 3">geometry',group,'.h5:/GradLa_gll</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            if (Tdomain%out_var_snap(OUT_GRAD_MU) == 1) then
                ! GRAD MU
                write(61,"(a)") '<Attribute Name="GradMu_gll" Center="Node" AttributeType="Vector">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="GradMu_gll" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                    ' 3">geometry',group,'.h5:/GradMu_gll</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
            ! KAPPA
            write(61,"(a)") '<Attribute Name="Kappa" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="Kappa" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                '">geometry',group,'.h5:/Kappa</DataItem>'
            write(61,"(a)") '</Attribute>'
            !C11 MC
            if (Tdomain%aniso .EQV. .true.) then
                write(61,"(a)") '<Attribute Name="C11" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C11" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C11</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C12" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C12" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C12</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C13" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C13" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C13</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C14" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C14" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C14</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C15" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C15" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C15</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C16" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C16" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C16</DataItem>'
                write(61,"(a)") '</Attribute>'


                write(61,"(a)") '<Attribute Name="C22" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C22" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C22</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C23" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C23" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C23</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C24" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C24" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C24</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C25" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C25" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C25</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C26" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C26" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C26</DataItem>'
                write(61,"(a)") '</Attribute>'


                write(61,"(a)") '<Attribute Name="C33" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C33" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C33</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C34" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C34" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C34</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C35" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C35" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C35</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C36" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C36" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C36</DataItem>'
                write(61,"(a)") '</Attribute>'


                write(61,"(a)") '<Attribute Name="C44" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C44" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C44</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C45" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C45" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C45</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C46" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C46" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C46</DataItem>'
                write(61,"(a)") '</Attribute>'


                write(61,"(a)") '<Attribute Name="C55" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C55" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C55</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="C56" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C56" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C56</DataItem>'
                write(61,"(a)") '</Attribute>'


                write(61,"(a)") '<Attribute Name="C66" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="C66" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/C66</DataItem>'
                write(61,"(a)") '</Attribute>'


            end if
            ! anisotropic-density fluid: density tensor rho_ij + bulk modulus kappa
            if (Tdomain%any_fanisodom) then
                write(61,"(a)") '<Attribute Name="rho11" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="rho11" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/rho11</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="rho22" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="rho22" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/rho22</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="rho33" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="rho33" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/rho33</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="rho12" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="rho12" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/rho12</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="rho13" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="rho13" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/rho13</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="rho23" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="rho23" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/rho23</DataItem>'
                write(61,"(a)") '</Attribute>'

                write(61,"(a)") '<Attribute Name="kappa" Center="Node" AttributeType="Scalar">'
                write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="kappa" Format="HDF" NumberType="Float" Precision="4" Dimensions="',nn, &
                        '">geometry',group,'.h5:/Kappa</DataItem>'
                write(61,"(a)") '</Attribute>'
            end if
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

    subroutine write_constant_fields(Tdomain, fid, outputs)
        use dom_solid
        use dom_fluid_aniso
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        type (output_var_t), intent(inout) :: outputs
        !
        real(fpp), dimension(:),allocatable :: mass, jac
#ifdef CPML
        real(fpp), dimension(:), allocatable :: alpha_pml, kappa_pml, dxi_k_pml
        integer :: i1, i2, d0, d1
#else
        real(fpp), dimension(:),allocatable :: dumpsx
        real(fpp) :: dx, dy, dz, dt
#endif
        real(fpp), dimension(:),allocatable :: dens, lamb, mu, kappa
        real(fpp), dimension(:),allocatable :: C11, C12, C13, C14, C15, C16
        real(fpp), dimension(:),allocatable :: C22, C23, C24, C25, C26
        real(fpp), dimension(:),allocatable :: C33, C34, C35, C36
        real(fpp), dimension(:),allocatable :: C44, C45, C46
        real(fpp), dimension(:),allocatable :: C55, C56
        real(fpp), dimension(:),allocatable :: C66
        real(fpp), dimension(:),allocatable :: rho11, rho22, rho33, rho12, rho13, rho23
        real(fpp) :: rB11, rB22, rB33, rB12, rB13, rB23, rdet
        real(fpp), dimension(:,:,:,:), allocatable :: grad_La
        real(fpp), dimension(:,:,:,:), allocatable :: grad_Mu
        real(fpp), dimension(:,:),allocatable :: grad_La_n,grad_Mu_n 
        integer :: ngll, idx
        integer :: i, j, k, n, lnum, nnodes_tot, bnum, ee
        integer :: domain_type, imat
        integer :: nnodes, ncells

        nnodes = outputs%nnodes
        ncells = outputs%ncells
        allocate(mass(0:nnodes-1))
        allocate(jac(0:nnodes-1))
        allocate(dens(0:nnodes-1))
        allocate(lamb(0:nnodes-1))
        allocate(mu(0:nnodes-1))
        allocate(kappa(0:nnodes-1))
        
        if (Tdomain%aniso .EQV. .true.) then 
            allocate(C11(0:nnodes-1))
            allocate(C12(0:nnodes-1))
            allocate(C13(0:nnodes-1))
            allocate(C14(0:nnodes-1))
            allocate(C15(0:nnodes-1))
            allocate(C16(0:nnodes-1))

            allocate(C22(0:nnodes-1))
            allocate(C23(0:nnodes-1))
            allocate(C24(0:nnodes-1))
            allocate(C25(0:nnodes-1))
            allocate(C26(0:nnodes-1))

            allocate(C33(0:nnodes-1))
            allocate(C34(0:nnodes-1))
            allocate(C35(0:nnodes-1))
            allocate(C36(0:nnodes-1))

            allocate(C44(0:nnodes-1))
            allocate(C45(0:nnodes-1))
            allocate(C46(0:nnodes-1))

            allocate(C55(0:nnodes-1))
            allocate(C56(0:nnodes-1))

            allocate(C66(0:nnodes-1))

        end if
        if (Tdomain%any_fanisodom) then
            allocate(rho11(0:nnodes-1))
            allocate(rho22(0:nnodes-1))
            allocate(rho33(0:nnodes-1))
            allocate(rho12(0:nnodes-1))
            allocate(rho13(0:nnodes-1))
            allocate(rho23(0:nnodes-1))
        end if
        if (Tdomain%out_var_snap(OUT_GRAD_LA) == 1) then
            allocate(grad_La_n(0:2,0:nnodes-1))
            grad_La_n = 0d0
        end if 
        if (Tdomain%out_var_snap(OUT_GRAD_MU) == 1) then
            allocate(grad_Mu_n(0:2,0:nnodes-1))
            grad_Mu_n = 0d0
        end if

#ifdef CPML
        allocate(alpha_pml(0:3*nnodes-1))
        alpha_pml = -1.0
        allocate(kappa_pml(0:3*nnodes-1))
        kappa_pml = -1.0
        allocate(dxi_k_pml(0:3*nnodes-1))
        dxi_k_pml = -1.0
#else
        allocate(dumpsx(0:nnodes-1))
        dumpsx = 0d0
#endif

        mass = 0d0
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
            case (DM_SOLID_CG)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
                                mass(idx) = Tdomain%sdom%MassMat(Tdomain%sdom%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case (DM_SOLID_DG)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
                                mass(idx) = Tdomain%sdomdg%MassMat(Tdomain%sdomdg%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case (DM_SOLID_CG_PML)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
#ifdef CPML
                                d0 = Tdomain%spmldom%D0(ee, bnum)
                                alpha_pml(d0 + 3*idx) = Tdomain%spmldom%Alpha_0(ee, i, j, k, bnum)
                                kappa_pml(d0 + 3*idx) = Tdomain%spmldom%Kappa_0(ee, i, j, k, bnum)
                                dxi_k_pml(d0 + 3*idx) = Tdomain%spmldom%dxi_k_0(ee, i, j, k, bnum)

                                i1 = Tdomain%spmldom%I1(ee, bnum)
                                i2 = Tdomain%spmldom%I2(ee, bnum)
                                if (i1/=-1) then
                                    d1 = Tdomain%spmldom%D1(ee, bnum)
                                    alpha_pml(d1 + 3*idx) = Tdomain%spmldom%Alpha_1(i, j, k, i1)
                                    kappa_pml(d1 + 3*idx) = Tdomain%spmldom%Kappa_1(i, j, k, i1)
                                    dxi_k_pml(d1 + 3*idx) = Tdomain%spmldom%dxi_k_1(i, j, k, i1)
                                end if
                                if (i2/=-1) then
                                    alpha_pml(2 + 3*idx) = Tdomain%spmldom%Alpha_2(i, j, k, i2)
                                    kappa_pml(2 + 3*idx) = Tdomain%spmldom%Kappa_2(i, j, k, i2)
                                    dxi_k_pml(2 + 3*idx) = Tdomain%spmldom%dxi_k_2(i, j, k, i2)
                                end if

#else
                                mass(idx) = Tdomain%spmldom%MassMat(Tdomain%spmldom%Idom_(i,j,k,bnum,ee))
                                dt = 2_fpp*Tdomain%TimeD%dtmin
                                dx = ((1_fpp/Tdomain%spmldom%PMLDumpSx_(i,j,k,1,bnum,ee))-1.)/dt
                                dy = ((1_fpp/Tdomain%spmldom%PMLDumpSy_(i,j,k,1,bnum,ee))-1.)/dt
                                dz = ((1_fpp/Tdomain%spmldom%PMLDumpSz_(i,j,k,1,bnum,ee))-1.)/dt
                                dumpsx(idx) = dx+dy+dz
#endif
                            endif
                        end do
                    end do
                end do
            case (DM_FLUID_CG)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
                                mass(idx) = Tdomain%fdom%MassMat(Tdomain%fdom%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case (DM_FLUID_CG_ANISO)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
                                mass(idx) = Tdomain%fanisodom%MassMat(Tdomain%fanisodom%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case (DM_FLUID_CG_PML)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
#ifdef CPML
                                i1 = Tdomain%fpmldom%I1(ee, bnum)
                                i2 = Tdomain%fpmldom%I2(ee, bnum)
                                d0 = Tdomain%fpmldom%D0(ee, bnum)
                                alpha_pml(d0 + 3*idx) = Tdomain%fpmldom%Alpha_0(ee, i, j, k, bnum)
                                kappa_pml(d0 + 3*idx) = Tdomain%fpmldom%Kappa_0(ee, i, j, k, bnum)
                                dxi_k_pml(d0 + 3*idx) = Tdomain%fpmldom%dxi_k_0(ee, i, j, k, bnum)

                                if (i1/=-1) then
                                    d1 = Tdomain%fpmldom%D1(ee, bnum)
                                    alpha_pml(d1 + 3*idx) = Tdomain%fpmldom%Alpha_1(i, j, k, i1)
                                    kappa_pml(d1 + 3*idx) = Tdomain%fpmldom%Kappa_1(i, j, k, i1)
                                    dxi_k_pml(d1 + 3*idx) = Tdomain%fpmldom%dxi_k_1(i, j, k, i1)
                                endif
                                if (i2/=-1) then
                                    alpha_pml(2 + 3*idx) = Tdomain%fpmldom%Alpha_2(i, j, k, i2)
                                    kappa_pml(2 + 3*idx) = Tdomain%fpmldom%Kappa_2(i, j, k, i2)
                                    dxi_k_pml(2 + 3*idx) = Tdomain%fpmldom%dxi_k_2(i, j, k, i2)
                                end if
#else
                                mass(idx) = Tdomain%fpmldom%MassMat(Tdomain%fpmldom%Idom_(i,j,k,bnum,ee))
#endif
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
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_CG)
                                jac(idx) = Tdomain%sdom%Jacob_   (i,j,k,bnum,ee)
                            case (DM_SOLID_DG)
                                jac(idx) = Tdomain%sdomdg%Jacob_ (i,j,k,bnum,ee)
                            case (DM_SOLID_CG_PML)
                                jac(idx) = Tdomain%spmldom%Jacob_(i,j,k,bnum,ee)
                            case (DM_FLUID_CG)
                                jac(idx) = Tdomain%fdom%Jacob_      (i,j,k,bnum,ee)
                            case (DM_FLUID_CG_ANISO)
                                jac(idx) = Tdomain%fanisodom%Jacob_ (i,j,k,bnum,ee)
                            case (DM_FLUID_CG_PML)
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
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_CG)
                                dens(idx) = Tdomain%sdom%Density_        (i,j,k,bnum,ee)
                            case (DM_SOLID_DG)
                                dens(idx) = Tdomain%sdomdg%Density_      (i,j,k,bnum,ee)
                            case (DM_SOLID_CG_PML)
                                dens(idx) = Tdomain%spmldom%Density_     (i,j,k,bnum,ee)
                            case (DM_FLUID_CG)
                                dens(idx) = 1.0D0/Tdomain%fdom%IDensity_ (i,j,k,bnum,ee)
                            case (DM_FLUID_CG_ANISO)
                                ! isotropic-equivalent density = 3 / trace(rho^{-1})
                                dens(idx) = 3d0/(Tdomain%fanisodom%IDensTensor_(0,i,j,k,bnum,ee) + &
                                                 Tdomain%fanisodom%IDensTensor_(1,i,j,k,bnum,ee) + &
                                                 Tdomain%fanisodom%IDensTensor_(2,i,j,k,bnum,ee))
                            case (DM_FLUID_CG_PML)
#ifdef CPML
                                dens(idx) = 0. ! Tdomain%fpmldom%Density_(i,j,k,bnum,ee) ! TODO
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
            if (Tdomain%out_var_snap(OUT_GRAD_LA) == 1) then
                allocate(grad_La(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                call get_solid_dom_grad_lambda(Tdomain%sdom, Tdomain%specel(n)%lnum, grad_La)
            end if
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_CG)
                                lamb(idx) = Tdomain%sdom%Lambda_        (i,j,k,bnum,ee)
                                if (Tdomain%out_var_snap(OUT_GRAD_LA) == 1) grad_La_n(0:2,idx) = grad_La(i,j,k,0:2)
                            case (DM_SOLID_DG)
                                lamb(idx) = Tdomain%sdomdg%Lambda_      (i,j,k,bnum,ee)
                            case (DM_SOLID_CG_PML)
                                lamb(idx) = Tdomain%spmldom%Lambda_     (i,j,k,bnum,ee)
                            case (DM_FLUID_CG)
                                lamb(idx) = Tdomain%fdom%Lambda_        (i,j,k,bnum,ee)
                            case (DM_FLUID_CG_ANISO)
                                lamb(idx) = Tdomain%fanisodom%Lambda_(i,j,k,bnum,ee)
                            case (DM_FLUID_CG_PML)
                                lamb(idx) = Tdomain%fpmldom%Lambda_     (i,j,k,bnum,ee)
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
            if (allocated(grad_La)) deallocate(grad_La)
        end do

        ! mu
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
            if (Tdomain%out_var_snap(OUT_GRAD_MU) == 1) then
                allocate(grad_Mu(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                call get_solid_dom_grad_mu(Tdomain%sdom, Tdomain%specel(n)%lnum, grad_Mu)
            end if
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_CG)
                                mu(idx) = Tdomain%sdom%Mu_(i,j,k,bnum,ee)
                                if (Tdomain%out_var_snap(OUT_GRAD_MU) == 1) grad_Mu_n(0:2,idx) = grad_Mu(i,j,k,0:2)
                            case (DM_SOLID_DG)
                                mu(idx) = Tdomain%sdomdg%Mu_(i,j,k,bnum,ee)
                            case (DM_SOLID_CG_PML)
                                mu(idx) = Tdomain%spmldom%Mu_(i,j,k,bnum,ee)
                            case (DM_FLUID_CG)
                                mu(idx) = -1d0
                            case (DM_FLUID_CG_ANISO)
                                mu(idx) = -1d0
                            case (DM_FLUID_CG_PML)
                                mu(idx) = -1d0
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
            if (allocated(grad_Mu)) deallocate(grad_Mu)
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
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_CG)
                                kappa(idx) = Tdomain%sdom%Kappa_(i,j,k,bnum,ee)
                            case (DM_SOLID_DG)
                                kappa(idx) = -1d0
                            case (DM_SOLID_CG_PML)
                                kappa(idx) = -1d0
                            case (DM_FLUID_CG)
                                kappa(idx) = -1d0
                            case (DM_FLUID_CG_ANISO)
                                kappa(idx) = Tdomain%fanisodom%Lambda_(i,j,k,bnum,ee)
                            case (DM_FLUID_CG_PML)
                                kappa(idx) = -1d0
                            case default
                                stop "unknown domain"
                        end select
                    end do
                end do
            end do
        end do
        !!! C11 (faisable pour tout Cii)
        if (Tdomain%aniso .EQV. .true.) then
            do n = 0,Tdomain%n_elem-1
                if (.not. Tdomain%specel(n)%OUTPUT) cycle
                ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
                bnum = Tdomain%specel(n)%lnum/VCHUNK
                ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
                do k = 0,ngll-1
                        do j = 0,ngll-1
                                do i = 0,ngll-1
                                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                                        select case (Tdomain%specel(n)%domain)
                                           case (DM_SOLID_CG)
                                                C11(idx) = Tdomain%sdom%Cij_(0,i,j,k,bnum,ee)
                                                C12(idx) = Tdomain%sdom%Cij_(1,i,j,k,bnum,ee)
                                                C13(idx) = Tdomain%sdom%Cij_(2,i,j,k,bnum,ee)
                                                C14(idx) = Tdomain%sdom%Cij_(3,i,j,k,bnum,ee)
                                                C15(idx) = Tdomain%sdom%Cij_(4,i,j,k,bnum,ee)
                                                C16(idx) = Tdomain%sdom%Cij_(5,i,j,k,bnum,ee)

                                                C22(idx) = Tdomain%sdom%Cij_(6,i,j,k,bnum,ee)
                                                C23(idx) = Tdomain%sdom%Cij_(7,i,j,k,bnum,ee)
                                                C24(idx) = Tdomain%sdom%Cij_(8,i,j,k,bnum,ee)
                                                C25(idx) = Tdomain%sdom%Cij_(9,i,j,k,bnum,ee)
                                                C26(idx) = Tdomain%sdom%Cij_(10,i,j,k,bnum,ee)

                                                C33(idx) = Tdomain%sdom%Cij_(11,i,j,k,bnum,ee)
                                                C34(idx) = Tdomain%sdom%Cij_(12,i,j,k,bnum,ee)
                                                C35(idx) = Tdomain%sdom%Cij_(13,i,j,k,bnum,ee)
                                                C36(idx) = Tdomain%sdom%Cij_(14,i,j,k,bnum,ee)

                                                C44(idx) = Tdomain%sdom%Cij_(15,i,j,k,bnum,ee)
                                                C45(idx) = Tdomain%sdom%Cij_(16,i,j,k,bnum,ee)
                                                C46(idx) = Tdomain%sdom%Cij_(17,i,j,k,bnum,ee)
                                                
                                                C55(idx) = Tdomain%sdom%Cij_(18,i,j,k,bnum,ee)
                                                C56(idx) = Tdomain%sdom%Cij_(19,i,j,k,bnum,ee)
                                                
                                                C66(idx) = Tdomain%sdom%Cij_(20,i,j,k,bnum,ee)

                                           case (DM_SOLID_DG)
                                                C11(idx) = -1d0
                                           case (DM_SOLID_CG_PML)
                                                C11(idx) = -1d0
                                           case (DM_FLUID_CG)
                                                C11(idx) = -1d0
                                           case (DM_FLUID_CG_ANISO)
                                                C11(idx) = -1d0
                                           case (DM_FLUID_CG_PML)
                                                C11(idx) = -1d0
                                           case default
                                                stop "unknown domain"
                                         end select
                                  end do
                         end do
                end do
            end do
        end if
        !!! rho_ij for the anisotropic-density fluid (Capdeville & Cance 2015):
        !!! m_IDensTensor stores the inverse-density tensor rho^{-1}_ij; output the density
        !!! tensor rho_ij = (rho^{-1})^{-1} (the quantity itself, not its inverse). kappa is
        !!! filled above directly from m_Lambda. det~0 -> -1 (skip, do not abort).
        if (Tdomain%any_fanisodom) then
            do n = 0,Tdomain%n_elem-1
                if (.not. Tdomain%specel(n)%OUTPUT) cycle
                ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
                bnum = Tdomain%specel(n)%lnum/VCHUNK
                ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (Tdomain%specel(n)%domain == DM_FLUID_CG_ANISO) then
                                rB11 = Tdomain%fanisodom%IDensTensor_(0,i,j,k,bnum,ee)
                                rB22 = Tdomain%fanisodom%IDensTensor_(1,i,j,k,bnum,ee)
                                rB33 = Tdomain%fanisodom%IDensTensor_(2,i,j,k,bnum,ee)
                                rB12 = Tdomain%fanisodom%IDensTensor_(3,i,j,k,bnum,ee)
                                rB13 = Tdomain%fanisodom%IDensTensor_(4,i,j,k,bnum,ee)
                                rB23 = Tdomain%fanisodom%IDensTensor_(5,i,j,k,bnum,ee)
                                rdet = rB11*(rB22*rB33-rB23*rB23) &
                                     - rB12*(rB12*rB33-rB23*rB13) &
                                     + rB13*(rB12*rB23-rB22*rB13)
                                if (abs(rdet) <= tiny(1._fpp)) then
                                    rho11(idx) = -1d0; rho22(idx) = -1d0; rho33(idx) = -1d0
                                    rho12(idx) = -1d0; rho13(idx) = -1d0; rho23(idx) = -1d0
                                else
                                    rho11(idx) =  (rB22*rB33 - rB23*rB23)/rdet
                                    rho22(idx) =  (rB11*rB33 - rB13*rB13)/rdet
                                    rho33(idx) =  (rB11*rB22 - rB12*rB12)/rdet
                                    rho12(idx) = -(rB12*rB33 - rB13*rB23)/rdet
                                    rho13(idx) =  (rB12*rB23 - rB13*rB22)/rdet
                                    rho23(idx) = -(rB11*rB23 - rB12*rB13)/rdet
                                end if
                            else
                                rho11(idx) = -1d0; rho22(idx) = -1d0; rho33(idx) = -1d0
                                rho12(idx) = -1d0; rho13(idx) = -1d0; rho23(idx) = -1d0
                            end if
                        end do
                    end do
                end do
            end do
        end if
        call grp_write_real_1d(outputs, fid, "Mass", nnodes, mass, nnodes_tot)
#ifdef CPML
#ifdef DEBUG_CPML
        call grp_write_real_1d(outputs, fid, "Alpha_PML", 3*nnodes, alpha_pml, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Kappa_PML", 3*nnodes, kappa_pml, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Dxi_K_PML", 3*nnodes, dxi_k_pml, nnodes_tot)
#endif
#else
        call grp_write_real_1d(outputs, fid, "Alpha", nnodes, dumpsx, nnodes_tot)
#endif
        
        call grp_write_real_1d(outputs, fid, "Jac", nnodes, jac, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Dens", nnodes, dens, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Lamb", nnodes, lamb, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Mu", nnodes, mu, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Kappa", nnodes, kappa, nnodes_tot)
        
        if (Tdomain%aniso .EQV. .true.) then
            call grp_write_real_1d(outputs, fid, "C11", nnodes, C11, nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C12", nnodes, C12,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C13", nnodes, C13,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C14", nnodes, C14,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C15", nnodes, C15,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C16", nnodes, C16,  nnodes_tot)
        
            call grp_write_real_1d(outputs, fid, "C22", nnodes, C22,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C23", nnodes, C23,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C24", nnodes, C24,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C25", nnodes, C25,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C26", nnodes, C26,  nnodes_tot)

            call grp_write_real_1d(outputs, fid, "C33", nnodes, C33,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C34", nnodes, C34,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C35", nnodes, C35,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C36", nnodes, C36,  nnodes_tot)

            call grp_write_real_1d(outputs, fid, "C44", nnodes, C44,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C45", nnodes, C45,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C46", nnodes, C46,  nnodes_tot)

            call grp_write_real_1d(outputs, fid, "C55", nnodes, C55,  nnodes_tot)
            call grp_write_real_1d(outputs, fid, "C56", nnodes, C56,  nnodes_tot)

            call grp_write_real_1d(outputs, fid, "C66", nnodes, C66,  nnodes_tot)
        end if
        if (Tdomain%any_fanisodom) then
            ! density tensor rho_ij (kappa is written above as "Kappa")
            call grp_write_real_1d(outputs, fid, "rho11", nnodes, rho11, nnodes_tot)
            call grp_write_real_1d(outputs, fid, "rho22", nnodes, rho22, nnodes_tot)
            call grp_write_real_1d(outputs, fid, "rho33", nnodes, rho33, nnodes_tot)
            call grp_write_real_1d(outputs, fid, "rho12", nnodes, rho12, nnodes_tot)
            call grp_write_real_1d(outputs, fid, "rho13", nnodes, rho13, nnodes_tot)
            call grp_write_real_1d(outputs, fid, "rho23", nnodes, rho23, nnodes_tot)
        end if
        ! GRAD LAMBDA
        if (Tdomain%out_var_snap(OUT_GRAD_LA) == 1) then
            call grp_write_real_2d(outputs, fid, "GradLa_gll", 3, nnodes, grad_La_n, nnodes_tot)
        end if
        ! GRAD MU
        if (Tdomain%out_var_snap(OUT_GRAD_MU) == 1) then
            call grp_write_real_2d(outputs, fid, "GradMu_gll", 3, nnodes, grad_Mu_n, nnodes_tot)
        end if

        call grp_write_int_1d(outputs, fid, "Dom", nnodes, outputs%domains, nnodes_tot)
        deallocate(mass,jac)
        deallocate(dens, lamb, mu, kappa)
        if(allocated(grad_La)) deallocate(grad_La)
        if(allocated(grad_Mu)) deallocate(grad_Mu)
        if(allocated(grad_La_n)) deallocate(grad_La_n)
        if(allocated(grad_Mu_n)) deallocate(grad_Mu_n)
        
        if(allocated(C11)) deallocate(C11)
        if(allocated(C12)) deallocate(C12)
        if(allocated(C13)) deallocate(C13)
        if(allocated(C14)) deallocate(C14)
        if(allocated(C15)) deallocate(C15)
        if(allocated(C16)) deallocate(C16)

        if(allocated(C22)) deallocate(C22)
        if(allocated(C23)) deallocate(C23)
        if(allocated(C24)) deallocate(C24)
        if(allocated(C25)) deallocate(C25)
        if(allocated(C26)) deallocate(C26)

        if(allocated(C33)) deallocate(C33)
        if(allocated(C34)) deallocate(C34)
        if(allocated(C35)) deallocate(C35)
        if(allocated(C36)) deallocate(C36)

        if(allocated(C44)) deallocate(C44)
        if(allocated(C45)) deallocate(C45)
        if(allocated(C46)) deallocate(C46)

        if(allocated(C55)) deallocate(C55)
        if(allocated(C56)) deallocate(C56)

        if(allocated(C66)) deallocate(C66)

        if(allocated(rho11)) deallocate(rho11)
        if(allocated(rho22)) deallocate(rho22)
        if(allocated(rho33)) deallocate(rho33)
        if(allocated(rho12)) deallocate(rho12)
        if(allocated(rho13)) deallocate(rho13)
        if(allocated(rho23)) deallocate(rho23)

#ifdef CPML

        deallocate(alpha_pml)
        deallocate(kappa_pml)
        deallocate(dxi_k_pml)
#endif

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
