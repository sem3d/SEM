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

    subroutine grp_write_real_2d(outputs, parent_id, name, dim1, dim2, data, ntot_nodes)
        type (output_var_t), intent (INOUT):: outputs
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

        call MPI_Gatherv(data, dim1*dim2, MPI_DOUBLE_PRECISION, all_data, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, outputs%comm, ierr)
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
        type(output_var_t), intent(in) :: outputs
        integer(HID_T), intent(in) :: parent_id
        character(len=*), intent(in) :: fname
        real(fpp), dimension(:) :: field
        !
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: dset_id
        integer :: hdferr, ierr

        call MPI_Gatherv(field, outputs%ncells, MPI_DOUBLE_PRECISION, &
            outputs%all_data_1d_c, outputs%counts_c, outputs%displs_c, MPI_DOUBLE_PRECISION, 0,&
            outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = outputs%ntot_cells
            call create_dset(parent_id, trim(fname), H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, outputs%all_data_1d_c, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
    end subroutine write_1d_var_c

    subroutine write_1d_var_n(outputs, parent_id, fname, field)
        type(output_var_t), intent(in) :: outputs
        integer(HID_T), intent(in) :: parent_id
        character(len=*), intent(in) :: fname
        real(fpp), dimension(:) :: field
        !
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: dset_id
        integer :: hdferr, ierr

        call MPI_Gatherv(field, outputs%nnodes, MPI_DOUBLE_PRECISION, &
            outputs%all_data_1d_n, outputs%counts_n, outputs%displs_n, MPI_DOUBLE_PRECISION, 0,&
            outputs%comm, ierr)
        if (outputs%rank==0) then
            dims(1) = outputs%ntot_nodes
            call create_dset(parent_id, trim(fname), H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, outputs%all_data_1d_n, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
    end subroutine write_1d_var_n

    subroutine write_2d_var_vecn(outputs, parent_id, fname, field)
        type(output_var_t), intent(in) :: outputs
        integer(HID_T), intent(in) :: parent_id
        character(len=*), intent(in) :: fname
        real(fpp), dimension(:,:) :: field
        !
        integer(HSIZE_T), dimension(2) :: dims
        integer(HID_T) :: dset_id
        integer :: hdferr, ierr

        call MPI_Gatherv(field, 3*outputs%nnodes, MPI_DOUBLE_PRECISION, &
            outputs%all_data_2d_n, outputs%counts2d_n, outputs%displs2d_n, MPI_DOUBLE_PRECISION, 0,&
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
        type(output_var_t), intent(in) :: outputs
        integer(HID_T), intent(in) :: parent_id
        integer, dimension(0:), intent(in) :: out_variables
        !

        ! P_ENERGY
        if (out_variables(OUT_ENERGYP) == 1) then
            call write_1d_var_c(outputs, parent_id, "P_energy", outputs%P_energy)
        end if
        ! S_ENERGY
        if (out_variables(OUT_ENERGYS) == 1) then
            call write_1d_var_c(outputs, parent_id, "S_energy", outputs%S_energy)
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

        call MPI_Gatherv(data, dim1, MPI_DOUBLE_PRECISION, all_data, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, outputs%comm, ierr)
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

        call allocate_fields(outputs, Tdomain%out_variables, Tdomain%nl_flag)


        call write_constant_fields(Tdomain, fid, outputs)

        if (outputs%rank==0) then
            call h5fclose_f(fid, hdferr)
        endif

        if (rg==0) call write_master_xdmf(Tdomain)
    end subroutine write_snapshot_geom

    subroutine write_global_nodes(Tdomain, fid, outputs)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        type(output_var_t), intent(in) :: outputs
        !
        real, dimension(:,:), allocatable :: nodes
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
        if (out_flags(OUT_DEPLA     ) == 1) allocate(outputs%displ(0:2,0:nnodes-1))
        if (out_flags(OUT_VITESSE   ) == 1) allocate(outputs%veloc(0:2,0:nnodes-1))
        if (out_flags(OUT_ACCEL     ) == 1) allocate(outputs%accel(0:2,0:nnodes-1))
        if (out_flags(OUT_PRESSION  ) == 1) allocate(outputs%press_n(0:nnodes-1))
        ! sortie par element
        if (out_flags(OUT_ENERGYP   ) == 1) allocate(outputs%P_energy(0:ncells-1))
        if (out_flags(OUT_ENERGYS   ) == 1) allocate(outputs%S_energy(0:ncells-1))
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
        if (out_flags(OUT_ENERGYS   ) == 1) outputs%S_energy   = 0.
        if (out_flags(OUT_EPS_VOL   ) == 1) outputs%eps_vol    = 0.
        if (out_flags(OUT_PRESSION  ) == 1) outputs%press_c    = 0.
        if (out_flags(OUT_PRESSION  ) == 1) outputs%press_n    = 0.
        if (out_flags(OUT_EPS_DEV   ) == 1) outputs%eps_dev    = 0.
        if (out_flags(OUT_STRESS_DEV) == 1) outputs%sig_dev    = 0.
        if (out_flags(OUT_EPS_DEV_PL) == 1) outputs%eps_dev_pl = 0.
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
        use dom_fluid
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
        integer :: ngll
        integer :: i, j, k, n, m, ind
        integer :: nnodes
        type(Element), pointer :: el
        type(subdomain), pointer :: sub_dom_mat
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
        nnodes = outputs%nnodes
        out_variables(:) = Tdomain%out_variables(:)

        call create_dir_sorties(Tdomain, isort)
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
                        ind = outputs%irenum(el%Iglobnum(i,j,k))
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
                        if (out_variables(OUT_DEPLA)   == 1)  outputs%displ(0:2,ind) = fieldU(i,j,k,0:2)
                        if (out_variables(OUT_VITESSE) == 1)  outputs%veloc(0:2,ind) = fieldV(i,j,k,0:2)
                        if (out_variables(OUT_ACCEL)   == 1)  outputs%accel(0:2,ind) = fieldA(i,j,k,0:2)
                        if (out_variables(OUT_PRESSION) == 1) outputs%press_n(ind) = fieldP(i,j,k)
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
            if (out_variables(OUT_ENERGYS) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, S_energy, outputs%S_energy)
            endif
            if (out_variables(OUT_EPS_VOL) == 1) then
                call evaluate_cell_centers(ngll, GLLc, cell_start, eps_vol, outputs%eps_vol)
            endif
            if (out_variables(OUT_EPS_DEV) == 1) then
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

    subroutine write_xdmf(Tdomain, isort, outputs, out_variables)
        implicit none
        type (domain), intent (IN):: Tdomain
        type(output_var_t), intent(in) :: outputs
        integer, intent(in) :: isort
        integer, dimension(0:), intent(in) :: out_variables
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: i, nn, ne, group
        real :: time
        character(len=11) :: R2label(0:20)
#ifdef CPML
        integer :: j
#endif

        nn = outputs%ntot_nodes
        ne = outputs%ntot_cells
        group = outputs%group
        call semname_xdmf(group, fnamef)

        R2label( 0) = "R2_L120_uxx"
        R2label( 1) = "R2_L120_uxy"
        R2label( 2) = "R2_L120_uxz"
        R2label( 3) = "R2_L021_uyx"
        R2label( 4) = "R2_L021_uyy"
        R2label( 5) = "R2_L021_uyz"
        R2label( 6) = "R2_L012_uzx"
        R2label( 7) = "R2_L012_uzy"
        R2label( 8) = "R2_L012_uzz"
        R2label( 9) = "R2_L0_uyy"
        R2label(10) = "R2_L0_uyz"
        R2label(11) = "R2_L0_uzy"
        R2label(12) = "R2_L0_uzz"
        R2label(13) = "R2_L1_uxx"
        R2label(14) = "R2_L1_uxz"
        R2label(15) = "R2_L1_uzx"
        R2label(16) = "R2_L1_uzz"
        R2label(17) = "R2_L2_uxx"
        R2label(18) = "R2_L2_uxy"
        R2label(19) = "R2_L2_uyx"
        R2label(20) = "R2_L2_uyy"

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
            if (out_variables(OUT_ENERGYS) == 1) call write_xdmf_attr_scalar_cells("S_energy", ne, i, group, "S_energy")
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
#ifdef CPML
            write(61,"(a)") '<Attribute Name="GlobCoord_PML" Center="Node" AttributeType="Vector">'
            write(61,"(a,I9,a,I4.4,a)") '<DataItem Name="GlobCoord_PML" Format="HDF" NumberType="Float" Precision="4" Dimensions="3 ',nn, &
                '">geometry',group,'.h5:/GlobCoord_PML</DataItem>'
            write(61,"(a)") '</Attribute>'

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

            call write_xdmf_attr_vector_nodes("R1_x", nn, i, group, "R1_x")
            call write_xdmf_attr_vector_nodes("R1_y", nn, i, group, "R1_y")
            call write_xdmf_attr_vector_nodes("R1_z", nn, i, group, "R1_z")

            do j = 0, 20
                call write_xdmf_attr_scalar_nodes(trim(R2label(j)), nn, i, group, trim(R2label(j)))
            end do
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

    subroutine write_constant_fields(Tdomain, fid, outputs)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        type (output_var_t), intent(inout) :: outputs
        !
        real, dimension(:),allocatable :: mass, jac
#ifdef CPML
        real, dimension(:), allocatable :: globcoord_pml, alpha_pml, kappa_pml, dxi_k_pml
        real, dimension(:), allocatable :: R1_x, R1_y, R1_z
        real, dimension(:), allocatable :: R2_L120_uxx, R2_L120_uxy, R2_L120_uxz
        real, dimension(:), allocatable :: R2_L021_uyx, R2_L021_uyy, R2_L021_uyz
        real, dimension(:), allocatable :: R2_L012_uzx, R2_L012_uzy, R2_L012_uzz
        real, dimension(:), allocatable :: R2_L0_uyy, R2_L0_uyz, R2_L0_uzy, R2_L0_uzz
        real, dimension(:), allocatable :: R2_L1_uxx, R2_L1_uxz, R2_L1_uzx, R2_L1_uzz
        real, dimension(:), allocatable :: R2_L2_uxx, R2_L2_uxy, R2_L2_uyx, R2_L2_uyy
        integer :: dir
#else
        real, dimension(:),allocatable :: dumpsx
        real(fpp) :: dx, dy, dz, dt
#endif
        real, dimension(:),allocatable :: dens, lamb, mu, kappa
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
#ifdef CPML
        allocate(globcoord_pml(0:2+3*(nnodes-1))) ! 0-based: 2+ for idx 0
        globcoord_pml = 0.0

        allocate(alpha_pml(0:2+3*(nnodes-1))) ! 0-based: 2+ for idx 0
        alpha_pml = -1.0
        allocate(kappa_pml(0:2+3*(nnodes-1))) ! 0-based: 2+ for idx 0
        kappa_pml = -1.0
        allocate(dxi_k_pml(0:2+3*(nnodes-1))) ! 0-based: 2+ for idx 0
        dxi_k_pml = -1.0

        allocate(R1_x(0:2+3*(nnodes-1))) ! 0-based: 2+ for idx 0
        R1_x = 0.0
        allocate(R1_y(0:2+3*(nnodes-1))) ! 0-based: 2+ for idx 0
        R1_y = 0.0
        allocate(R1_z(0:2+3*(nnodes-1))) ! 0-based: 2+ for idx 0
        R1_z = 0.0

        allocate(R2_L120_uxx(0:nnodes-1), R2_L120_uxy(0:nnodes-1), R2_L120_uxz(0:nnodes-1)                       )
        allocate(R2_L021_uyx(0:nnodes-1), R2_L021_uyy(0:nnodes-1), R2_L021_uyz(0:nnodes-1)                       )
        allocate(R2_L012_uzx(0:nnodes-1), R2_L012_uzy(0:nnodes-1), R2_L012_uzz(0:nnodes-1)                       )
        allocate(R2_L0_uyy  (0:nnodes-1), R2_L0_uyz  (0:nnodes-1), R2_L0_uzy  (0:nnodes-1), R2_L0_uzz(0:nnodes-1))
        allocate(R2_L1_uxx  (0:nnodes-1), R2_L1_uxz  (0:nnodes-1), R2_L1_uzx  (0:nnodes-1), R2_L1_uzz(0:nnodes-1))
        allocate(R2_L2_uxx  (0:nnodes-1), R2_L2_uxy  (0:nnodes-1), R2_L2_uyx  (0:nnodes-1), R2_L2_uyy(0:nnodes-1))
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
            case (DM_SOLID)
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
            case (DM_SOLID_PML)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
#ifdef CPML
                                alpha_pml(Tdomain%spmldom%D0(ee, bnum) + 3*idx) = Tdomain%spmldom%Alpha_0(ee, i, j, k, bnum)
                                if(allocated(Tdomain%spmldom%Alpha_1)) & ! May NOT be allocated
                                alpha_pml(Tdomain%spmldom%D1(ee, bnum) + 3*idx) = Tdomain%spmldom%Alpha_1(i, j, k, 0)
                                if(allocated(Tdomain%spmldom%Alpha_2)) then ! May NOT be allocated
                                    dir = 0 + 1 + 2 ! All possible directions
                                    dir = dir - Tdomain%spmldom%I1(ee, bnum) - Tdomain%spmldom%I2(ee, bnum) ! Remove known directions
                                    alpha_pml(dir + 3*idx) = Tdomain%spmldom%Alpha_2(i, j, k, 0)
                                end if

                                kappa_pml(Tdomain%spmldom%D0(ee, bnum) + 3*idx) = Tdomain%spmldom%Kappa_0(ee, i, j, k, bnum)
                                if(allocated(Tdomain%spmldom%Kappa_1)) & ! May NOT be allocated
                                kappa_pml(Tdomain%spmldom%D1(ee, bnum) + 3*idx) = Tdomain%spmldom%Kappa_1(i, j, k, 0)
                                if(allocated(Tdomain%spmldom%Kappa_2)) then ! May NOT be allocated
                                    dir = 0 + 1 + 2 ! All possible directions
                                    dir = dir - Tdomain%spmldom%I1(ee, bnum) - Tdomain%spmldom%I2(ee, bnum) ! Remove known directions
                                    kappa_pml(dir + 3*idx) = Tdomain%spmldom%Kappa_2(i, j, k, 0)
                                end if

                                dxi_k_pml(Tdomain%spmldom%D0(ee, bnum) + 3*idx) = Tdomain%spmldom%dxi_k_0(ee, i, j, k, bnum)
                                if(allocated(Tdomain%spmldom%dxi_k_1)) & ! May NOT be allocated
                                dxi_k_pml(Tdomain%spmldom%D1(ee, bnum) + 3*idx) = Tdomain%spmldom%dxi_k_1(i, j, k, 0)
                                if(allocated(Tdomain%spmldom%dxi_k_2)) then ! May NOT be allocated
                                    dir = 0 + 1 + 2 ! All possible directions
                                    dir = dir - Tdomain%spmldom%I1(ee, bnum) - Tdomain%spmldom%I2(ee, bnum) ! Remove known directions
                                    dxi_k_pml(dir + 3*idx) = Tdomain%spmldom%dxi_k_2(i, j, k, 0)
                                end if

                                globcoord_pml(0 + 3*idx) = Tdomain%spmldom%GlobCoord(0, idx)
                                globcoord_pml(1 + 3*idx) = Tdomain%spmldom%GlobCoord(1, idx)
                                globcoord_pml(2 + 3*idx) = Tdomain%spmldom%GlobCoord(2, idx)

                                R2_L120_uxx(idx) = Tdomain%spmldom%R2_0(ee, 0, i, j, k, bnum)
                                R2_L120_uxy(idx) = Tdomain%spmldom%R2_0(ee, 1, i, j, k, bnum)
                                R2_L120_uxz(idx) = Tdomain%spmldom%R2_0(ee, 2, i, j, k, bnum)
                                R2_L021_uyx(idx) = Tdomain%spmldom%R2_0(ee, 3, i, j, k, bnum)
                                R2_L021_uyy(idx) = Tdomain%spmldom%R2_0(ee, 4, i, j, k, bnum)
                                R2_L021_uyz(idx) = Tdomain%spmldom%R2_0(ee, 5, i, j, k, bnum)
                                R2_L012_uzx(idx) = Tdomain%spmldom%R2_0(ee, 6, i, j, k, bnum)
                                R2_L012_uzy(idx) = Tdomain%spmldom%R2_0(ee, 7, i, j, k, bnum)
                                R2_L012_uzz(idx) = Tdomain%spmldom%R2_0(ee, 8, i, j, k, bnum)
                                select case(Tdomain%spmldom%D0(ee, bnum))
                                case(0)
                                    R1_x(0 + 3*idx) = Tdomain%spmldom%R1_0(ee, 0, i, j, k, bnum)
                                    R1_x(1 + 3*idx) = Tdomain%spmldom%R1_0(ee, 1, i, j, k, bnum)
                                    R1_x(2 + 3*idx) = Tdomain%spmldom%R1_0(ee, 2, i, j, k, bnum)

                                    R2_L0_uyy(idx) = Tdomain%spmldom%R2_0(ee, 4, i, j, k, bnum)
                                    R2_L0_uyz(idx) = Tdomain%spmldom%R2_0(ee, 5, i, j, k, bnum)
                                    R2_L0_uzy(idx) = Tdomain%spmldom%R2_0(ee, 7, i, j, k, bnum)
                                    R2_L0_uzz(idx) = Tdomain%spmldom%R2_0(ee, 8, i, j, k, bnum)
                                    R2_L1_uxx(idx) = 0.
                                    R2_L1_uxz(idx) = 0.
                                    R2_L1_uzx(idx) = 0.
                                    R2_L1_uzz(idx) = 0.
                                    R2_L2_uxx(idx) = 0.
                                    R2_L2_uxy(idx) = 0.
                                    R2_L2_uyx(idx) = 0.
                                    R2_L2_uyy(idx) = 0.
                                case(1)
                                    R1_y(0 + 3*idx) = Tdomain%spmldom%R1_0(ee, 0, i, j, k, bnum)
                                    R1_y(1 + 3*idx) = Tdomain%spmldom%R1_0(ee, 1, i, j, k, bnum)
                                    R1_y(2 + 3*idx) = Tdomain%spmldom%R1_0(ee, 2, i, j, k, bnum)

                                    R2_L0_uyy(idx) = 0.
                                    R2_L0_uyz(idx) = 0.
                                    R2_L0_uzy(idx) = 0.
                                    R2_L0_uzz(idx) = 0.
                                    R2_L1_uxx(idx) = Tdomain%spmldom%R2_0(ee, 0, i, j, k, bnum)
                                    R2_L1_uxz(idx) = Tdomain%spmldom%R2_0(ee, 2, i, j, k, bnum)
                                    R2_L1_uzx(idx) = Tdomain%spmldom%R2_0(ee, 6, i, j, k, bnum)
                                    R2_L1_uzz(idx) = Tdomain%spmldom%R2_0(ee, 8, i, j, k, bnum)
                                    R2_L2_uxx(idx) = 0.
                                    R2_L2_uxy(idx) = 0.
                                    R2_L2_uyx(idx) = 0.
                                    R2_L2_uyy(idx) = 0.
                                case(2)
                                    R1_z(0 + 3*idx) = Tdomain%spmldom%R1_0(ee, 0, i, j, k, bnum)
                                    R1_z(1 + 3*idx) = Tdomain%spmldom%R1_0(ee, 1, i, j, k, bnum)
                                    R1_z(2 + 3*idx) = Tdomain%spmldom%R1_0(ee, 2, i, j, k, bnum)

                                    R2_L0_uyy(idx) = 0.
                                    R2_L0_uyz(idx) = 0.
                                    R2_L0_uzy(idx) = 0.
                                    R2_L0_uzz(idx) = 0.
                                    R2_L1_uxx(idx) = 0.
                                    R2_L1_uxz(idx) = 0.
                                    R2_L1_uzx(idx) = 0.
                                    R2_L1_uzz(idx) = 0.
                                    R2_L2_uxx(idx) = Tdomain%spmldom%R2_0(ee, 0, i, j, k, bnum)
                                    R2_L2_uxy(idx) = Tdomain%spmldom%R2_0(ee, 1, i, j, k, bnum)
                                    R2_L2_uyx(idx) = Tdomain%spmldom%R2_0(ee, 3, i, j, k, bnum)
                                    R2_L2_uyy(idx) = Tdomain%spmldom%R2_0(ee, 4, i, j, k, bnum)
                                end select
#else
                                mass(idx) = Tdomain%spmldom%MassMat(Tdomain%spmldom%Idom_(i,j,k,bnum,ee))
                                dt = 2d0*Tdomain%TimeD%dtmin
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
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
                                mass(idx) = Tdomain%fdom%MassMat(Tdomain%fdom%Idom_(i,j,k,bnum,ee))
                            endif
                        end do
                    end do
                end do
            case (DM_FLUID_PML)
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (outputs%domains(idx)==domain_type) then
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
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
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
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                dens(idx) = Tdomain%sdom%Density_        (i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
                                dens(idx) = Tdomain%spmldom%Density_     (i,j,k,bnum,ee)
                            case (DM_FLUID)
                                dens(idx) = 1.0D0/Tdomain%fdom%IDensity_ (i,j,k,bnum,ee)
                            case (DM_FLUID_PML)
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
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                lamb(idx) = Tdomain%sdom%Lambda_        (i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
                                lamb(idx) = Tdomain%spmldom%Lambda_     (i,j,k,bnum,ee)
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
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                mu(idx) = Tdomain%sdom%Mu_(i,j,k,bnum,ee)
                            case (DM_SOLID_PML)
                                mu(idx) = Tdomain%spmldom%Mu_(i,j,k,bnum,ee)
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
                        idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
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

        call grp_write_real_1d(outputs, fid, "Mass", nnodes, mass, nnodes_tot)
#ifdef CPML
        call grp_write_real_1d(outputs, fid, "GlobCoord_PML", 3*nnodes, globcoord_pml, nnodes_tot)

        call grp_write_real_1d(outputs, fid, "Alpha_PML", 3*nnodes, alpha_pml, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Kappa_PML", 3*nnodes, kappa_pml, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Dxi_K_PML", 3*nnodes, dxi_k_pml, nnodes_tot)

        call grp_write_real_1d(outputs, fid, "R1_x", 3*nnodes, R1_x, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R1_y", 3*nnodes, R1_y, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R1_z", 3*nnodes, R1_z, nnodes_tot)

        call grp_write_real_1d(outputs, fid, "R2_L120_uxx", nnodes, R2_L120_uxx, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L120_uxy", nnodes, R2_L120_uxy, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L120_uxz", nnodes, R2_L120_uxz, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L021_uyx", nnodes, R2_L021_uyx, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L021_uyy", nnodes, R2_L021_uyy, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L021_uyz", nnodes, R2_L021_uyz, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L012_uzx", nnodes, R2_L012_uzx, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L012_uzy", nnodes, R2_L012_uzy, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L012_uzz", nnodes, R2_L012_uzz, nnodes_tot)

        call grp_write_real_1d(outputs, fid, "R2_L0_uyy", nnodes, R2_L0_uyy, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L0_uyz", nnodes, R2_L0_uyz, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L0_uzy", nnodes, R2_L0_uzy, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L0_uzz", nnodes, R2_L0_uzz, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L1_uxx", nnodes, R2_L1_uxx, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L1_uxz", nnodes, R2_L1_uxz, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L1_uzx", nnodes, R2_L1_uzx, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L1_uzz", nnodes, R2_L1_uzz, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L2_uxx", nnodes, R2_L2_uxx, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L2_uxy", nnodes, R2_L2_uxy, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L2_uyx", nnodes, R2_L2_uyx, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "R2_L2_uyy", nnodes, R2_L2_uyy, nnodes_tot)
#else
        call grp_write_real_1d(outputs, fid, "Alpha", nnodes, dumpsx, nnodes_tot)
#endif
        call grp_write_real_1d(outputs, fid, "Jac", nnodes, jac, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Dens", nnodes, dens, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Lamb", nnodes, lamb, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Mu", nnodes, mu, nnodes_tot)
        call grp_write_real_1d(outputs, fid, "Kappa", nnodes, kappa, nnodes_tot)

        call grp_write_int_1d(outputs, fid, "Dom", nnodes, outputs%domains, nnodes_tot)
        deallocate(mass,jac)
        deallocate(dens, lamb, mu, kappa)
#ifdef CPML
        deallocate(globcoord_pml)

        deallocate(alpha_pml)
        deallocate(kappa_pml)
        deallocate(dxi_k_pml)

        deallocate(R1_x)
        deallocate(R1_y)
        deallocate(R1_z)

        deallocate(R2_L120_uxx, R2_L120_uxy, R2_L120_uxz           )
        deallocate(R2_L021_uyx, R2_L021_uyy, R2_L021_uyz           )
        deallocate(R2_L012_uzx, R2_L012_uzy, R2_L012_uzz           )
        deallocate(R2_L0_uyy  , R2_L0_uyz  , R2_L0_uzy  , R2_L0_uzz)
        deallocate(R2_L1_uxx  , R2_L1_uxz  , R2_L1_uzx  , R2_L1_uzz)
        deallocate(R2_L2_uxx  , R2_L2_uxy  , R2_L2_uyx  , R2_L2_uyy)
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
