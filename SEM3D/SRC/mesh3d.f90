!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module mesh3d
    use sdomain
    use hdf5
    use sem_hdf5
    implicit none
contains

    subroutine  read_interface(inter, fid, pfx)
        use hdf5
        use sem_hdf5
        implicit none
        type(inter_num), intent(inout) :: inter
        integer(HID_T), intent(in) :: fid
        character(len=*), intent(in) :: pfx
        !
        integer, dimension(:,:), allocatable :: itemp2

        call read_attr_int(fid, "n_"//trim(pfx)//"_faces",    inter%surf0%n_faces)
        call read_attr_int(fid, "n_"//trim(pfx)//"_edges",    inter%surf0%n_edges)
        call read_attr_int(fid, "n_"//trim(pfx)//"_vertices", inter%surf0%n_vertices)

        inter%surf1%n_faces = inter%surf0%n_faces
        inter%surf1%n_edges = inter%surf0%n_edges
        inter%surf1%n_vertices = inter%surf0%n_vertices

        call allocate_interface(inter)

        if (inter%surf0%n_faces/=0) then
            call read_dataset(fid, trim(pfx)//"_faces", itemp2)
            inter%surf0%if_faces = itemp2(1,:)
            inter%surf1%if_faces = itemp2(2,:)
            deallocate(itemp2)
            call read_dataset(fid, trim(pfx)//"_orient", itemp2)
            inter%surf0%if_norm = itemp2(1,:)
            inter%surf1%if_norm = itemp2(2,:)
            deallocate(itemp2)
        end if
        if (inter%surf0%n_edges/=0) then
            call read_dataset(fid, trim(pfx)//"_edges", itemp2)
            inter%surf0%if_edges = itemp2(1,:)
            inter%surf1%if_edges = itemp2(2,:)
            deallocate(itemp2)
        end if
        if (inter%surf0%n_vertices/=0) then
            call read_dataset(fid, trim(pfx)//"_vertices", itemp2)
            inter%surf0%if_vertices = itemp2(1,:)
            inter%surf1%if_vertices = itemp2(2,:)
            deallocate(itemp2)
        end if
    end subroutine read_interface

    ! Reads the main attributes of a mesh file (dim, mesh sizes)
    subroutine read_mesh_attributes(Tdomain, fid)
        type(domain), intent(inout) :: Tdomain
        integer(HID_T) :: fid
        integer :: num_processors
        !
        call read_attr_int (fid, "ndims", Tdomain%n_dime)
        call read_attr_int(fid, "n_processors", num_processors)
        call read_attr_int(fid, "n_materials", Tdomain%n_mat)
        call read_attr_int(fid, "n_elements",  Tdomain%n_elem)
        call read_attr_int(fid, "n_faces",     Tdomain%n_face)
        call read_attr_int(fid, "n_edges",     Tdomain%n_edge)
        call read_attr_int(fid, "n_vertices",  Tdomain%n_vertex)
        !call read_attr_int(fid, "neumann_present_loc",  Tdomain%n_neumannfind)
        !
        if(Tdomain%n_dime /=3)   &
            stop "No general code for the time being: only 3D propagation"
        !
        if(Tdomain%nb_procs/=num_processors) then
            write(*,*) "Mesh file for rank", Tdomain%rank, "has incoherent number of processors"
            write(*,*) "you should check you meshfiles in sem/*"
            stop 1
        endif
    end subroutine read_mesh_attributes

    ! Read the definitions of the local mesh elements (elem, faces, edges, vertices)
    subroutine read_mesh_elements(Tdomain, fid)
        type(domain), intent(inout) :: Tdomain
        integer(HID_T) :: fid
        !
        integer,   allocatable, dimension(:,:) :: itemp2, itemp2b
        integer,   allocatable, dimension(:)   :: itemp
        real(fpp), allocatable, dimension(:,:) :: rtemp2
        integer :: i, mat
        !
        ! Global nodes' coordinates for each proc.
        call read_dataset(fid, "local_nodes", rtemp2)
        Tdomain%n_glob_nodes = size(rtemp2,2)
        allocate (Tdomain%Coord_nodes(0:Tdomain%n_dime-1,0:Tdomain%n_glob_nodes-1))
        Tdomain%Coord_nodes = rtemp2
        deallocate(rtemp2)

        ! Elements (material and solid or fluid and if fluid: Dirichlet boundary?)
        call read_dataset(fid, "material", itemp)
        if (Tdomain%n_elem /= size(itemp)) then
            write(*,*) "N_elem:", Tdomain%n_elem
            write(*,*) "itemp:", size(itemp)
            stop "Incoherent number of elements"
        end if

        allocate(Tdomain%specel(0:Tdomain%n_elem-1))
        allocate (Tdomain%subD_exist(0:Tdomain%n_mat-1))
        Tdomain%subD_exist(:) = .false.

        do i=0,Tdomain%n_elem-1
            call init_element(Tdomain%specel(i))
            Tdomain%specel(i)%mat_index = itemp(i+1)
            Tdomain%specel(i)%OUTPUT = .true.
            Tdomain%subD_exist(itemp(i+1)) = .true.
        end do

        deallocate(itemp)

        ! Read elements definitions
        ! n_nodes : number of control nodes (8 or 27)
        call read_dataset(fid, "elements", itemp2)
        Tdomain%n_nodes = size(itemp2,1)
        do i = 0, Tdomain%n_elem-1
            allocate(Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
            Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1) = itemp2(:,i+1)
            mat = Tdomain%specel(i)%mat_index

        end do
        deallocate(itemp2)

        ! Faces and elements properties related to faces
        call read_dataset(fid, "faces", itemp2)
        call read_dataset(fid, "faces_def", itemp2b)
        call read_dataset(fid, "faces_dom", itemp)

        allocate(Tdomain%sFace(0:Tdomain%n_face-1))
        do i=0,Tdomain%n_face-1
            call init_face(Tdomain%sFace(i))
            Tdomain%sFace(i)%inodes = itemp2b(:,i+1)
            Tdomain%sFace(i)%domain = itemp(i+1)
        enddo
        do i = 0, Tdomain%n_elem-1
            Tdomain%specel(i)%Near_Faces(0:5) = itemp2(:,i+1)
        enddo
        deallocate(itemp2, itemp2b, itemp)

        ! Edges
        call read_dataset(fid, "edges", itemp2)
        call read_dataset(fid, "edges_def", itemp2b)
        call read_dataset(fid, "edges_dom", itemp)

        allocate(Tdomain%sEdge(0:Tdomain%n_edge-1))
        do i=0,Tdomain%n_edge-1
            call init_edge(Tdomain%sEdge(i))
            Tdomain%sEdge(i)%inodes = itemp2b(:,i+1)
            Tdomain%sEdge(i)%domain = itemp(i+1)
        enddo
        do i = 0, Tdomain%n_elem-1
            Tdomain%specel(i)%Near_Edges(0:11) = itemp2(:,i+1)
        enddo
        deallocate(itemp2, itemp2b, itemp)

        ! Vertices
        call read_dataset(fid, "vertices", itemp2)
        call read_dataset(fid, "vertices_dom", itemp)

        allocate(Tdomain%sVertex(0:Tdomain%n_vertex-1))
        do i=0,Tdomain%n_vertex-1
            call init_vertex(Tdomain%sVertex(i))
            Tdomain%sVertex(i)%domain = itemp(i+1)
        enddo
        do i = 0,Tdomain%n_elem-1
            Tdomain%specel(i)%Near_Vertices(0:7) = itemp2(:,i+1)
        enddo
        deallocate(itemp2, itemp)

    end subroutine read_mesh_elements

    subroutine compute_material_boundaries(Tdomain)
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        !
        integer :: i, j, k, mat, nod
        logical :: pml
        real(fpp) :: Coord_local
        !
        do i = 0, Tdomain%n_mat-1
            Tdomain%sSubDomain(i)%MinBound_loc = MAX_DOUBLE
            Tdomain%sSubDomain(i)%MaxBound_loc = -MAX_DOUBLE
            Tdomain%sSubDomain(i)%present = .false.
        end do

        do i = 0, Tdomain%n_elem-1
            mat = Tdomain%specel(i)%mat_index
            Tdomain%sSubDomain(mat)%present = .true.
            pml = is_pml(Tdomain%sSubdomain(mat))

            do j = 0, Tdomain%n_nodes-1
                nod = Tdomain%specel(i)%Control_Nodes(j)
                do k = 0,2
                    Coord_local = Tdomain%Coord_nodes(k, nod)
                    if (pml) then
                         if (Tdomain%sSubDomain(mat)%pml_width(k) /= 0) then
                            Coord_local = Tdomain%sSubDomain(mat)%pml_pos(k)
                         end if
                    end if
                    !Max
                    if(Coord_local > Tdomain%sSubDomain(mat)%MaxBound_loc(k)) &
                        Tdomain%sSubDomain(mat)%MaxBound_loc(k) = Coord_local
                    !Min
                    if(Coord_local < Tdomain%sSubDomain(mat)%MinBound_loc(k)) &
                        Tdomain%sSubDomain(mat)%MinBound_loc(k) = Coord_local
                end do
            end do
        end do
        do i = 0, Tdomain%n_mat-1
            do k=0,2
                if (Tdomain%sSubDomain(i)%MinBound_loc(k)>Tdomain%sSubDomain(i)%MaxBound_loc(k)) then
                    ! Material is not present on this domain

                end if
            end do
        end do

    end subroutine compute_material_boundaries

    subroutine read_comm_proc_data(proc_id, scomm)
        integer(HID_T), intent(in) :: proc_id
        type(comm), intent(inout) :: scomm
        !
        integer, allocatable, dimension(:)   :: itemp
        integer :: j
        !
        call read_attr_int(proc_id, "proc_dest", scomm%dest)
        call read_attr_int(proc_id, "n_faces", scomm%nb_faces)
        call read_attr_int(proc_id, "n_edges", scomm%nb_edges)
        call read_attr_int(proc_id, "n_vertices", scomm%nb_vertices)

        if(scomm%nb_faces > 0)then
            allocate(scomm%faces(0:scomm%nb_faces-1))
            call read_dataset(proc_id, "faces", itemp)
            do j = 0,scomm%nb_faces-1
                scomm%faces(j)        = itemp(j+1)
            enddo
            deallocate(itemp)
        endif
        if(scomm%nb_edges > 0)then
            call read_dataset(proc_id, "edges", itemp)
            allocate(scomm%edges(0:scomm%nb_edges-1))
            do j = 0,scomm%nb_edges-1
                scomm%edges(j) = itemp(j+1)
            enddo
            deallocate(itemp)
        endif
        if(scomm%nb_vertices > 0)then
            call read_dataset(proc_id, "vertices", itemp)
            allocate(scomm%vertices(0:scomm%nb_vertices-1))
            do j = 0,scomm%nb_vertices-1
                scomm%vertices(j) = itemp(j+1)
            enddo
            deallocate(itemp)
        endif

    end subroutine read_comm_proc_data

    subroutine read_one_surface(surf, pfx, gid)
        implicit none
        type(surf_num), intent(inout) :: surf
        integer(HID_T), intent(in) :: gid
        character(len=*) :: pfx
        !
        integer, allocatable, dimension(:) :: itemp
        !
        call init_surface(surf)

        call read_attr_int(gid, "n_"//trim(pfx)//"_faces",    surf%n_faces)
        call read_attr_int(gid, "n_"//trim(pfx)//"_edges",    surf%n_edges)
        call read_attr_int(gid, "n_"//trim(pfx)//"_vertices", surf%n_vertices)

        call allocate_surface(surf)

        if (surf%n_faces/=0) then
            call read_dataset(gid, trim(pfx)//"_faces", itemp)
            surf%if_faces = itemp
            deallocate(itemp)
            call read_dataset(gid, trim(pfx)//"_orient", itemp)
            surf%if_norm = itemp
            deallocate(itemp)
        end if
        if (surf%n_edges/=0) then
            call read_dataset(gid, trim(pfx)//"_edges", itemp)
            surf%if_edges = itemp
            deallocate(itemp)
        end if
        if (surf%n_vertices/=0) then
            call read_dataset(gid, trim(pfx)//"_vertices", itemp)
            surf%if_vertices = itemp
            deallocate(itemp)
        end if
    end subroutine read_one_surface

    subroutine get_surface_domain(gid, domain)

        implicit none
        integer, intent(inout)             :: domain
        integer(HID_T), intent(in)         :: gid
        character(len=4), dimension(4)     :: pfx=(/"sl  ", &
            "fl  ", &
            "spml", &
            "fpml"/)
        integer, allocatable, dimension(:) :: itemp
        integer                            :: n_faces, i

        do i=1,4
            call read_attr_int(gid, "n_"//trim(pfx(i))//"_faces",n_faces)
            if (n_faces/=0) then
                call read_dataset(gid, trim(pfx(i))//"_faces_dom", itemp)
                if (MAXVAL(itemp).ne.MINVAL(itemp)) then
                    stop "Error : faces_dom uncorrectly denied in surfaces list"
                endif
                domain=MINVAL(itemp)
                deallocate(itemp)
            endif
        enddo

    end subroutine get_surface_domain

    subroutine read_surfaces(Tdomain, gid)
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer(HID_T), intent(in) :: gid
        !
        integer :: n_surfaces
        integer :: storage_type, nlinks, max_corder, ierr
        character(len=100) :: surfname
        integer(HSIZE_T) :: namesz, i
        integer(HID_T) :: surf_id
        call read_attr_int(gid, "n_surfaces", n_surfaces)
        !write(*,*) "NSURFACES=", n_surfaces
        allocate(Tdomain%sSurfaces(0:n_surfaces-1))
        ! Get group info
        call H5Gget_info_f(gid, storage_type, nlinks, max_corder, ierr)

        do i=0,nlinks-1
            call H5Lget_name_by_idx_f(gid, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, i, surfname, ierr, namesz)
            call H5Gopen_f(gid, trim(surfname), surf_id, ierr)
            !
            ! Read one surface
            !
            Tdomain%sSurfaces(i)%name = surfname
            Tdomain%sSurfaces(i)%cond_type = COND_NONE
            call read_one_surface(Tdomain%sSurfaces(i)%surf_sl  , "sl"  , surf_id)
            call read_one_surface(Tdomain%sSurfaces(i)%surf_sldg, "sldg", surf_id)
            call read_one_surface(Tdomain%sSurfaces(i)%surf_fl  , "fl"  , surf_id)
            call read_one_surface(Tdomain%sSurfaces(i)%surf_spml, "spml", surf_id)
            call read_one_surface(Tdomain%sSurfaces(i)%surf_fpml, "fpml", surf_id)
            call get_surface_domain(surf_id, Tdomain%sSurfaces(i)%domain)
        end do
    end subroutine read_surfaces

    subroutine read_mesh_file_h5(Tdomain)
        use mpi
        use sem_c_bindings
        use semdatafiles, only : MAX_FILE_SIZE, strrank
        use constants
        use surface_input

        implicit none
        !
        type(domain), intent(inout) :: Tdomain
        integer :: rg
        integer :: i
        logical :: neumann_log
        !
        integer(HID_T) :: fid, proc_id, surf_id
        integer :: hdferr, ierr
        character(len=10) :: proc_grp
        integer, allocatable, dimension(:)   :: nb_elems_per_proc
        character(Len=MAX_FILE_SIZE) :: fname
        !
        rg = Tdomain%rank
        !
        call init_hdf5()
        !
        fname = trim(adjustl(Tdomain%mesh_file))//".h5"
        if (sem_check_file(trim(fname))==0) then
            write (*,*) "Process ",rg, " can't open his mesh_file ", trim(fname)
            stop
        endif

        call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, fid, hdferr)
        if (hdferr .ne. 0) then
            write(*,*) "h5fopen : failed to open file : ",fname," rank:", rg
            stop 1
        end if
        !
        !-- Reading mesh properties
        neumann_log = .false.
        Tdomain%logicD%Neumann_local_present = .false.
        !
        ! solid - pml interface
        call init_interface(Tdomain%intSolPml)
        ! fluid - pml interface
        call init_interface(Tdomain%intFluPml)
        ! solid - fluid interface
        call init_interface(Tdomain%SF%intSolFlu)
        call init_interface(Tdomain%SF%intSolFluPml)

        call read_mesh_attributes(Tdomain, fid)

        !Subdomains allocation
        allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))
        !
        if (Tdomain%logicD%Neumann) then
            Tdomain%logicD%Neumann_local_present = .true.
            Tdomain%logicD%Neumann = .true.
            neumann_log = .true.
        endif

        if(neumann_log .neqv. Tdomain%logicD%Neumann) then
            write(*,*) rg,"neumann (input.spec)=",Tdomain%logicD%Neumann
            write(*,*) rg,"neumann (mesh file)=",neumann_log
            stop "Introduction of Neumann B.C.: mesh and input files not in coincidence."
        endif

        call read_mesh_elements(Tdomain, fid)

        call read_interface(Tdomain%intSolPml, fid, "spml")
        call read_interface(Tdomain%intFluPml, fid, "fpml")
        ! Solid-fluid properties, eventually
        call read_interface(Tdomain%SF%intSolFlu, fid, "sf")
        call read_interface(Tdomain%SF%intSolFluPml, fid, "sfpml")
        Tdomain%logicD%SF_local_present = ( &
            Tdomain%SF%intSolFlu%surf0%n_faces + &
            Tdomain%SF%intSolFlu%surf0%n_edges + &
            Tdomain%SF%intSolFlu%surf0%n_vertices + &
            Tdomain%SF%intSolFluPml%surf0%n_faces + &
            Tdomain%SF%intSolFluPml%surf0%n_edges + &
            Tdomain%SF%intSolFluPml%surf0%n_vertices) /= 0

        call h5gopen_f(fid, "Surfaces", surf_id, hdferr)
        call read_surfaces(Tdomain, surf_id)
        call h5gclose_f(surf_id, hdferr)
        ! Interproc communications
        Tdomain%tot_comm_proc = 0
        call read_attr_int(fid, "tot_comm_proc", Tdomain%tot_comm_proc)
        allocate (Tdomain%sComm(0:Tdomain%nb_procs-1))
        do i = 0,Tdomain%tot_comm_proc-1
            write(proc_grp,"(a,I4.4)") "Proc", i
            call h5gopen_f(fid, trim(adjustl(proc_grp)), proc_id, hdferr)
            call read_comm_proc_data(proc_id, Tdomain%sComm(i) )

        end do

        allocate(nb_elems_per_proc(0:8*((Tdomain%nb_procs-1)/8+1)))
        nb_elems_per_proc = 0
        call MPI_Gather(Tdomain%n_elem, 1, MPI_INTEGER, nb_elems_per_proc, 1, MPI_INTEGER, 0, Tdomain%communicateur, ierr)
        if (rg==0) then
            write(*,*) "Mesh read correctly, elements per proc:"
            do i=0,Tdomain%nb_procs-1,8
                write(*,'(I5.5,a,8I6)') i, ":", nb_elems_per_proc(i:i+7)
            end do
        end if
        deallocate(nb_elems_per_proc)

    end subroutine read_mesh_file_h5
    !!                                                                                                            !!
    !!
    !!
end module mesh3d

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
