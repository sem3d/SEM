module msnapshots
    use sdomain
    use hdf5
    use sem_hdf5
    use semdatafiles
    use mpi
    use mfields
    use sem_c_config, only : sem_mkdir
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


    subroutine grp_write_fields(Tdomain, parent_id, dim2, displ, veloc, accel, press, ntot_nodes)
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: parent_id
        integer, intent(in) :: dim2
        real, dimension(0:2,0:dim2-1), intent(in) :: displ, veloc, accel
        real, dimension(0:dim2-1), intent(in) :: press
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

        ! On commence par recuperer les champs 1D
        ntot_nodes = 0
        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = ntot_nodes
                ntot_nodes = ntot_nodes+counts(n)
            end do
            allocate(all_data_1d(0:ntot_nodes-1))
        end if
        
        call MPI_Gatherv(press, dim2, MPI_DOUBLE_PRECISION, all_data_1d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = ntot_nodes
            call create_dset(parent_id, "press", H5T_IEEE_F32LE, dims(1), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_1d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
            deallocate(all_data_1d)
        end if

        if (Tdomain%output_rank==0) then
            do n=0, Tdomain%nb_output_procs-1
                displs(n) = displs(n)*3
                counts(n) = counts(n)*3
            end do
            allocate(all_data_2d(0:2,0:ntot_nodes-1))
        end if
        call MPI_Gatherv(veloc, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = 3
            dims(2) = ntot_nodes
            call create_dset_2d(parent_id, "veloc", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
        call MPI_Gatherv(displ, 3*dim2, MPI_DOUBLE_PRECISION, all_data_2d, counts, displs, &
            MPI_DOUBLE_PRECISION, 0, Tdomain%comm_output, ierr)
        if (Tdomain%output_rank==0) then
            dims(1) = 3
            dims(2) = ntot_nodes
            call create_dset_2d(parent_id, "displ", H5T_IEEE_F32LE, dims(1), dims(2), dset_id)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, all_data_2d, dims, hdferr)
            call h5dclose_f(dset_id, hdferr)
        end if
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
        integer :: status, ierr, domain_type
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
            domain_type = get_domain(Tdomain%specel(n))

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

    subroutine create_dir_sorties(Tdomain, rg, isort)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer,intent(in) :: isort, rg
        character(Len=MAX_FILE_SIZE) :: temp
        integer :: code, ierr

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
    subroutine write_snapshot_geom(Tdomain, rg)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: rg
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        integer :: hdferr, code, ierr
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes
        integer, dimension(0:Tdomain%n_proc-1) :: nodes_per_proc
        integer :: group
        integer, dimension(:), allocatable :: domains

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

        call write_global_nodes(Tdomain, rg, fid, irenum, nnodes)
        
        call write_elem_connectivity(Tdomain, rg, fid, irenum)

        call write_constant_fields(Tdomain, fid, irenum, nnodes, domains)

        if (Tdomain%output_rank==0) then
            call h5fclose_f(fid, hdferr)
        endif

        call mpi_gather(nnodes, 1, MPI_INTEGER, nodes_per_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (rg==0) call write_master_xdmf(Tdomain, nodes_per_proc)
    end subroutine write_snapshot_geom


    subroutine write_global_nodes(Tdomain, rg, fid, irenum, nnodes)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, intent(in) :: rg, nnodes
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

    subroutine write_elem_connectivity(Tdomain, rg, fid, irenum)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        integer, intent(in) :: rg
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
                        proc(count) = rg
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

    subroutine save_field_h5(Tdomain, rg, isort)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: rg, isort
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        real, dimension(:,:),allocatable :: displ, veloc, accel
        real, dimension(:), allocatable :: press
        real, dimension(:,:,:,:),allocatable :: field_displ, field_veloc, field_accel
        real, dimension(:,:,:),allocatable :: field_press
        real, dimension(:,:,:),allocatable :: field_phi, field_vphi
        integer :: domain_type
        integer, dimension(:), allocatable :: valence
        real, dimension(:,:,:,:), allocatable :: Depla
        integer :: hdferr
        integer :: ngllx, nglly, ngllz, idx
        integer :: i, j, k, n, i_dir
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes, group, nnodes_tot, mat
        integer, dimension(:), allocatable :: domains
        type(Element), pointer :: el
        type(subdomain), pointer :: sub_dom_mat

        call create_dir_sorties(Tdomain, rg, isort)

        call compute_saved_elements(Tdomain, irenum, nnodes, domains)

        allocate(displ(0:2,0:nnodes-1))
        allocate(veloc(0:2,0:nnodes-1))
        allocate(accel(0:2,0:nnodes-1))
        allocate(press(0:nnodes-1))
        allocate(valence(0:nnodes-1))

        ngllx = 0
        nglly = 0
        ngllz = 0
        valence(:) = 0
        veloc(:,:) = 0d0
        accel(:,:) = 0d0
        displ(:,:) = 0d0
        press(:) = 0d0

        do n = 0,Tdomain%n_elem-1
            el => Tdomain%specel(n)
            sub_dom_mat => Tdomain%sSubdomain(el%mat_index)
            if (.not. el%OUTPUT) cycle
            if (ngllx /= el%ngllx .or. &
                nglly /= el%nglly .or. &
                ngllz /= el%ngllz) then
                ngllx = el%ngllx
                nglly = el%nglly
                ngllz = el%ngllz
                if (allocated(field_veloc)) deallocate(field_veloc)
                if (allocated(field_accel)) deallocate(field_accel)
                if (allocated(field_press)) deallocate(field_press)
                if (allocated(field_displ)) deallocate(field_displ)
#if NEW_GLOBAL_METHOD
                if (allocated(field_phi)) deallocate(field_phi)
                if (allocated(field_vphi)) deallocate(field_vphi)
#endif
                allocate(field_displ(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                allocate(field_veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                allocate(field_accel(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                allocate(field_press(0:ngllx-1,0:nglly-1,0:ngllz-1))
#if NEW_GLOBAL_METHOD
                allocate(field_phi(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(field_vphi(0:ngllx-1,0:nglly-1,0:ngllz-1))
#endif
            endif
            
#if NEW_GLOBAL_METHOD
            domain_type = get_domain(el)
            select case(domain_type)
            case (DM_SOLID)
                call gather_field(el, field_displ, Tdomain%champs0%Depla, el%Isol)
                call gather_field(el, field_veloc, Tdomain%champs0%Veloc, el%Isol)
                call gather_field(el, field_accel, Tdomain%champs0%Forces, el%Isol)
                call pressure_solid(ngllx,nglly,ngllz,sub_dom_mat%htprimex,              &
                    sub_dom_mat%hprimey,sub_dom_mat%hprimez, &
                    el%InvGrad,field_displ, el%Lambda, el%Mu,field_press)
            case (DM_SOLID_PML)
                field_displ = 0
                call gather_field_pml(el, field_veloc, Tdomain%champs0%VelocPML, el%slpml%IsolPML)
                !call gather_field_pml(el, field_accel, Tdomain%champs1%ForcesPML, el%slpml%IsolPML)
                field_accel = 0
                field_press = 0
            case (DM_FLUID)
                field_displ = 0
                call gather_field_fluid(el, field_phi, Tdomain%champs0%Phi, el%IFlu)
                call fluid_velocity(ngllx,nglly,ngllz,sub_dom_mat%htprimex,              &
                            sub_dom_mat%hprimey,sub_dom_mat%hprimez, &
                            el%InvGrad,el%density,field_phi,field_veloc)
                call gather_field_fluid(el, field_vphi, Tdomain%champs0%VelPhi, el%IFlu)
                call fluid_velocity(ngllx,nglly,ngllz,sub_dom_mat%htprimex,              &
                            sub_dom_mat%hprimey,sub_dom_mat%hprimez, &
                            el%InvGrad,el%density,field_vphi,field_accel)
                field_press = -field_vphi
            case (DM_FLUID_PML)
                field_displ = 0
                field_veloc = 0
                field_accel = 0
                field_press = 0
            end select

            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        idx = irenum(el%Iglobnum(i,j,k))
                        if (domains(idx)==domain_type) then
                            valence(idx) = valence(idx)+1
                            displ(:,idx) = field_displ(i,j,k,:)
                            veloc(:,idx) = veloc(:,idx) + field_veloc(i,j,k,:)
                            accel(:,idx) = accel(:,idx) + field_accel(i,j,k,:)
                            press(idx) = field_press(i,j,k)
                        endif
                    enddo
                enddo
            enddo
#else
            call gather_elem_displ(Tdomain, n, field_displ)
            call gather_elem_veloc(Tdomain, n, field_veloc)
            call gather_elem_accel(Tdomain, n, field_accel)
            call gather_elem_press(Tdomain, n, field_press)

            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        idx = irenum(el%Iglobnum(i,j,k))
                        valence(idx) = valence(idx)+1
                        displ(:,idx) = field_displ(i,j,k,:)
                        veloc(:,idx) = veloc(:,idx) + field_veloc(i,j,k,:)
                        accel(:,idx) = accel(:,idx) + field_accel(i,j,k,:)
                        press(idx) = field_press(i,j,k)
                    enddo
                enddo
            enddo
#endif
        enddo

        ! normalization
        do i = 0,nnodes-1
            if (valence(i)/=0) then
                veloc(0:2,i) = veloc(0:2,i)/valence(i)
                accel(0:2,i) = accel(0:2,i)/valence(i)
            end if
        enddo

        if (Tdomain%output_rank==0) then
            group = rg/Tdomain%ngroup
            call semname_snap_result_file(group, isort, fnamef)
            call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        else
            fid = -1
        endif
        call grp_write_fields(Tdomain, fid, nnodes, displ, veloc, accel, press, nnodes_tot)

        if (Tdomain%output_rank==0) then
            call h5fclose_f(fid, hdferr)
            call write_xdmf(Tdomain, group, isort, nnodes_tot)
        endif
        deallocate(displ,veloc,accel,press,valence)
        if (allocated(field_displ)) deallocate(field_displ)
        if (allocated(field_veloc)) deallocate(field_veloc)
        if (allocated(field_accel)) deallocate(field_accel)
        if (allocated(field_press)) deallocate(field_press)
#if NEW_GLOBAL_METHOD
        if (allocated(field_phi)) deallocate(field_phi)
        if (allocated(field_vphi)) deallocate(field_vphi)
#endif
        call mpi_barrier(Tdomain%communicateur, hdferr)
    end subroutine save_field_h5

    subroutine write_master_xdmf(Tdomain, nodes_per_proc)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: n_procs, nelem, n_groups
        character (len=MAX_FILE_SIZE) :: fnamef
        integer, dimension(0:Tdomain%n_proc-1), intent(in) :: nodes_per_proc
        integer :: rg
        n_procs = Tdomain%n_proc
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
        do rg=0,n_groups-1
            write(61,"(a,I4.4,a)") '<xi:include href="mesh.',rg,'.xmf"/>'
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

    subroutine write_xdmf(Tdomain, rg, isort, nnodes)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer, intent(in) :: rg, isort, nnodes
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: i, nn, ne
        real :: time
        call semname_xdmf(rg, fnamef)

        nn = nnodes
        ne = Tdomain%n_hexa
        if (ne==0) then
            return
        end if
        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a,I4.4,a)") '<Grid CollectionType="Temporal" GridType="Collection" Name="space.',rg,'">'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mat" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',rg,'.h5:/Material</DataItem>'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Proc" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',rg,'.h5:/Proc</DataItem>'

        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Mass</DataItem>'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Jac" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Jac</DataItem>'
        time = 0
        do i=1,isort
            write(61,"(a,I4.4,a,I4.4,a)") '<Grid Name="mesh.',i,'.',rg,'">'
            write(61,"(a,F20.10,a)") '<Time Value="', time,'"/>'
            write(61,"(a,I8,a)") '<Topology Type="Hexahedron" NumberOfElements="',ne,'">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Int" Dimensions="',ne,' 8">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Elements'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Topology>'
            write(61,"(a)") '<Geometry Type="XYZ">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Nodes'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Geometry>'

            write(61,"(a,I8,a)") '<Attribute Name="Displ" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/displ'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a,I8,a)") '<Attribute Name="Veloc" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/veloc'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a,I8,a)") '<Attribute Name="Accel" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/accel'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a,I8,a)") '<Attribute Name="Pressure" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/press'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a)") '<Attribute Name="Domain" Center="Grid" AttributeType="Scalar" Dimensions="1">'
            write(61,"(a,I4,a)") '<DataItem Format="XML" Datatype="Int"  Dimensions="1">',rg,'</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I8,a)") '<Attribute Name="Mat" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',rg, &
                '"]/DataItem[@Name="Mat"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I8,a)") '<Attribute Name="Proc" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',rg, &
                '"]/DataItem[@Name="Proc"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I8,a)") '<Attribute Name="Mass" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',rg, &
                '"]/DataItem[@Name="Mass"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I8,a)") '<Attribute Name="Jac" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[@Name="space.',rg, &
                '"]/DataItem[@Name="Jac"]</DataItem>'
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
#if NEW_GLOBAL_METHOD
        integer :: domain_type
#endif
        

        allocate(mass(0:nnodes-1))
        allocate(jac(0:nnodes-1))
        mass = 0d0
#if ! NEW_GLOBAL_METHOD
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
#else
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            domain_type = get_domain(Tdomain%specel(n))
            select case(domain_type)
            case (DM_SOLID)
                do k = 0,ngllz-1
                    do j = 0,nglly-1
                        do i = 0,ngllx-1
                            idx = irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                            if (domains(idx)==domain_type) then
                                mass(idx) = Tdomain%MassMatSol(Tdomain%specel(n)%Isol(i,j,k))
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
                                mass(idx) = Tdomain%MassMatSolPML(Tdomain%specel(n)%slpml%IsolPML(i,j,k))
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
                                mass(idx) = Tdomain%MassMatFlu(Tdomain%specel(n)%IFlu(i,j,k))
                            endif
                        end do
                    end do
                end do
            case (DM_FLUID_PML)
                do k = 0,ngllz-1
                    do j = 0,nglly-1
                        do i = 0,ngllx-1
                            !TODO
                            stop "TODO : assemblage de la masse mat pour PML fluide !"
                        end do
                    end do
                end do
            end select
        enddo
        if (Tdomain%any_PML) deallocate(Tdomain%MassMatSolPml)
#endif
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

#if ! NEW_GLOBAL_METHOD
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
#endif
        call grp_write_real_1d(Tdomain, fid, "Mass", nnodes, mass, nnodes_tot)
        call grp_write_real_1d(Tdomain, fid, "Jac", nnodes, jac, nnodes_tot)
        deallocate(mass,jac)

    end subroutine write_constant_fields

end module msnapshots
