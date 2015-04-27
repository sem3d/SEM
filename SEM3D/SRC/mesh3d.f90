!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module mesh3d

contains

subroutine read_mesh_file_h5(Tdomain)
    use sdomain
    use mpi
    use hdf5
    use sem_hdf5
    use sem_c_bindings
    use displayCarvalhol
    use semdatafiles, only : MAX_FILE_SIZE
    implicit none
    !
    type(domain), intent(inout) :: Tdomain
    integer :: rg, num_processors
    integer :: i,j
    logical :: neumann_log
    !
    integer(HID_T) :: fid, proc_id
    integer :: hdferr, ierr
    integer, allocatable, dimension(:,:) :: itemp2, itemp2b
    integer, allocatable, dimension(:)   :: itemp, itempb
    real,    allocatable, dimension(:,:) :: rtemp2
    !real,    allocatable, dimension(:)   :: rtemp
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
    !
    !-- Reading mesh properties
    call read_attr_int (fid, "ndim", Tdomain%n_dime)
    call read_attr_bool(fid, "solid_fluid", Tdomain%logicD%solid_fluid)
    call read_attr_bool(fid, "all_fluid", Tdomain%logicD%all_fluid)
    call read_attr_bool(fid, "neumann_present", neumann_log)
    call read_attr_bool(fid, "neumann_present_loc", Tdomain%logicD%Neumann_local_present)
    call read_attr_bool(fid, "curve", Tdomain%curve)
    call read_attr_bool(fid, "solid_fluid_loc", Tdomain%logicD%SF_local_present)
    call read_attr_int(fid, "n_processors", num_processors)
    call read_attr_int(fid, "n_materials", Tdomain%n_mat)
    call read_attr_int(fid, "n_elements",  Tdomain%n_elem)
    call read_attr_int(fid, "n_faces",     Tdomain%n_face)
    call read_attr_int(fid, "n_edges",     Tdomain%n_edge)
    call read_attr_int(fid, "n_vertices",  Tdomain%n_vertex)
    !
    if(Tdomain%n_dime /=3)   &
        stop "No general code for the time being: only 3D propagation"
    !
    if(Tdomain%nb_procs/=num_processors) then
        write(*,*) "Mesh file for rank",rg,"has incoherent number of processors"
        write(*,*) "you should check you meshfiles in sem/*"
        stop 1
    endif
    !
    if(neumann_log .neqv. Tdomain%logicD%Neumann) then
        write(*,*) rg,"neumann (input.spec)=",Tdomain%logicD%Neumann
        write(*,*) rg,"neumann (mesh file)=",neumann_log
        stop "Introduction of Neumann B.C.: mesh and input files not in coincidence."
    endif

    !Subdomains allocation
    allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))

    ! Global nodes' coordinates for each proc.
    call read_dataset(fid, "local_nodes", rtemp2)
    Tdomain%n_glob_nodes = size(rtemp2,2)
    allocate (Tdomain%Coord_nodes(0:Tdomain%n_dime-1,0:Tdomain%n_glob_nodes-1))
    Tdomain%Coord_nodes = rtemp2
    deallocate(rtemp2)

    !call dispCarvalhol(transpose(Tdomain%Coord_nodes), "transpose(Tdomain%Coord_nodes)")

    ! Elements (material and solid or fluid and if fluid: Dirichlet boundary?)
    call read_dataset(fid, "material", itemp2)
    if (Tdomain%n_elem /= size(itemp2,2)) then
        write(*,*) "N_elem:", Tdomain%n_elem
        write(*,*) "itemp:", size(itemp2,1), size(itemp2,2)
        stop "Incoherent number of elements"
    end if
    allocate(Tdomain%specel(0:Tdomain%n_elem-1))
    do i=0,Tdomain%n_elem-1
        call init_element(Tdomain%specel(i))
        Tdomain%specel(i)%mat_index = itemp2(1,i+1)
        Tdomain%specel(i)%solid = itemp2(2,i+1) .ne. 0
        Tdomain%specel(i)%fluid_dirich = itemp2(3,i+1) .ne. 0
        Tdomain%specel(i)%OUTPUT = .true.
    enddo
    deallocate(itemp2)

    ! Read elements definitions
    ! n_nodes : number of control nodes (8 or 27)
    call read_dataset(fid, "elements", itemp2)
    Tdomain%n_nodes = size(itemp2,1)
    do i = 0, Tdomain%n_elem-1
        allocate(Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
        Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1) = itemp2(:,i+1)
    enddo
    deallocate(itemp2)

    ! Faces and elements properties related to faces
    call read_dataset(fid, "faces", itemp2)
    call read_dataset(fid, "faces_map", itemp2b)

    allocate(Tdomain%sFace(0:Tdomain%n_face-1))
    do i=0,Tdomain%n_face-1
        call init_face(Tdomain%sFace(i))
    enddo
    do i = 0, Tdomain%n_elem-1
        Tdomain%specel(i)%Near_Faces(0:5) = itemp2(:,i+1)
        Tdomain%specel(i)%Orient_Faces(0:5) = itemp2b(:,i+1)
    enddo
    deallocate(itemp2, itemp2b)

    ! Edges
    call read_dataset(fid, "edges", itemp2)
    call read_dataset(fid, "edges_map", itemp2b)

    allocate(Tdomain%sEdge(0:Tdomain%n_edge-1))
    do i=0,Tdomain%n_edge-1
        call init_edge(Tdomain%sEdge(i))
    enddo
    do i = 0, Tdomain%n_elem-1
        Tdomain%specel(i)%Near_Edges(0:11) = itemp2(:,i+1)
        Tdomain%specel(i)%Orient_Edges(0:11) = itemp2b(:,i+1)
    enddo
    deallocate(itemp2, itemp2b)

    ! Vertices
    call read_dataset(fid, "vertices", itemp2)
    call read_dataset(fid, "vertices_to_global", itemp)

    allocate(Tdomain%sVertex(0:Tdomain%n_vertex-1))
    do i=0,Tdomain%n_vertex-1
        call init_vertex(Tdomain%sVertex(i))
        Tdomain%sVertex(i)%global_numbering = itemp(i+1)
    enddo
    do i = 0,Tdomain%n_elem-1
        Tdomain%specel(i)%Near_Vertices(0:7) = itemp2(:,i+1)
    enddo
    deallocate(itemp2, itemp)

    ! Solid-fluid properties, eventually
    if(Tdomain%logicD%SF_local_present)then
        ! Edges and their orientation, for each SF face
        call read_dataset(fid, "sf_face_near_edges", itemp2)
        call read_dataset(fid, "sf_face_near_edges_orient", itemp2b)
        Tdomain%SF%SF_n_faces = size(itemp2,2)
        allocate(Tdomain%SF%SF_face(0:Tdomain%SF%SF_n_faces-1))
        do i = 0, Tdomain%SF%SF_n_faces-1
            Tdomain%SF%SF_face(i)%Near_Edges(0:3) = itemp2(:,i+1)
            Tdomain%SF%SF_face(i)%Near_Edges_Orient(0:3) = itemp2b(:,i+1)
        end do
        deallocate(itemp2, itemp2b)
        ! Vertices for each SF face
        call read_dataset(fid, "sf_face_near_vertices", itemp2)
        do i = 0, Tdomain%SF%SF_n_faces-1
            Tdomain%SF%SF_face(i)%Near_Vertices(0:3) = itemp2(:,i+1)
        end do
        deallocate(itemp2)
        ! associated fluid (0) and solid (1) faces; orientation Solid/Fluid
        call read_dataset(fid, "sf_face_glob_interface", itemp2)
        call read_dataset(fid, "sf_face_orient", itemp)
        do i = 0, Tdomain%SF%SF_n_faces-1
            Tdomain%SF%SF_face(i)%Face(0:1) = itemp2(1:2,i+1)
            Tdomain%SF%SF_face(i)%Orient_Face = itemp(i+1)
            Tdomain%SF%SF_face(i)%PML = .false.
        end do
        deallocate(itemp, itemp2)
        ! SF edges
        call read_dataset(fid, "sf_edge_glob_interface", itemp2)
        call read_dataset(fid, "sf_edge_orient", itemp)
        Tdomain%SF%SF_n_edges = size(itemp2,2)
        allocate(Tdomain%SF%SF_edge(0:Tdomain%SF%SF_n_edges-1))
        do i = 0, Tdomain%SF%SF_n_edges-1
            Tdomain%SF%SF_edge(i)%Edge(0:1) = itemp2(1:2,i+1)
            Tdomain%SF%SF_edge(i)%Orient_Edge = itemp(i+1)
            Tdomain%SF%SF_edge(i)%PML = .false.
        end do
        deallocate(itemp, itemp2)
        ! SF vertices
        call read_dataset(fid, "sf_vertex_glob_interface", itemp2)
        Tdomain%SF%SF_n_vertices = size(itemp2, 2)
        allocate(Tdomain%SF%SF_vertex(0:Tdomain%SF%SF_n_vertices-1))
        do i = 0, Tdomain%SF%SF_n_vertices-1
            Tdomain%SF%SF_vertex(i)%Vertex(0:1) = itemp2(1:2,i+1)
            Tdomain%SF%SF_Vertex(i)%PML = .false.
        end do
        deallocate(itemp2)
    end if

    ! Neumann B.C. properties, eventually
    if(Tdomain%logicD%Neumann_local_present)then
        ! Neumann properties
        ! Neumann faces
        call read_dataset(fid, "neu_face_near_edges", itemp2)
        ! Edges and their orientation, for each Neumann face
        call read_dataset(fid, "neu_face_near_edges_orient", itemp2b)
        Tdomain%Neumann%Neu_n_faces = size(itemp2,2)
        allocate(Tdomain%Neumann%Neu_face(0:Tdomain%Neumann%Neu_n_faces-1))
        do i = 0, Tdomain%Neumann%Neu_n_faces-1
            Tdomain%Neumann%Neu_face(i)%Near_Edges(0:3) = itemp2(:,i+1)
            Tdomain%Neumann%Neu_face(i)%Near_Edges_Orient(0:3) = itemp2b(:,i+1)
        end do
        deallocate(itemp2, itemp2b)
        ! Vertices for each Neumann face
        call read_dataset(fid, "neu_face_near_vertices", itemp2)
        ! associated face
        call read_dataset(fid, "neu_face_glob_interface", itemp)
        do i = 0, Tdomain%Neumann%Neu_n_faces-1
            Tdomain%Neumann%Neu_face(i)%Near_Vertices(0:3) = itemp2(:,i+1)
            Tdomain%Neumann%Neu_face(i)%Face = itemp(i+1)
        end do
        deallocate(itemp, itemp2)
        ! Neumann edges
        call read_dataset(fid, "neu_edge_glob_interface", itemp)
        Tdomain%Neumann%Neu_n_edges = size(itemp,1)
        allocate(Tdomain%Neumann%Neu_edge(0:Tdomain%Neumann%Neu_n_edges-1))
        do i = 0, Tdomain%Neumann%Neu_n_edges-1
            Tdomain%Neumann%Neu_edge(i)%Edge = itemp(i+1)
        end do
        deallocate(itemp)
        ! Neumann vertices
        call read_dataset(fid,"neu_vertex_glob_interface", itemp)
        Tdomain%Neumann%Neu_n_vertices = size(itemp,1)
        allocate(Tdomain%Neumann%Neu_vertex(0:Tdomain%Neumann%Neu_n_vertices-1))
        do i = 0, Tdomain%Neumann%Neu_n_vertices-1
            Tdomain%Neumann%Neu_vertex(i)%Vertex = itemp(i+1)
        end do
        deallocate(itemp)
    end if

    ! Interproc communications
    Tdomain%tot_comm_proc = 0
    call read_attr_int(fid, "tot_comm_proc", Tdomain%tot_comm_proc)
    allocate (Tdomain%sComm(0:Tdomain%nb_procs-1))
    do i = 0,Tdomain%tot_comm_proc-1
        write(proc_grp,"(a,I4.4)") "Proc", i
        call h5gopen_f(fid, trim(adjustl(proc_grp)), proc_id, hdferr)
        call read_attr_int(proc_id, "proc_dest", Tdomain%sComm(i)%dest)
        call read_attr_int(proc_id, "n_faces", Tdomain%sComm(i)%nb_faces)
        call read_attr_int(proc_id, "n_edges", Tdomain%sComm(i)%nb_edges)
        call read_attr_int(proc_id, "n_vertices", Tdomain%sComm(i)%nb_vertices)
        if(Tdomain%logicD%SF_local_present)then
            call read_attr_int(proc_id, "n_sf_faces", Tdomain%sComm(i)%SF_nf_shared)
            call read_attr_int(proc_id, "n_sf_edges", Tdomain%sComm(i)%SF_ne_shared)
            call read_attr_int(proc_id, "n_sf_vertices", Tdomain%sComm(i)%SF_nv_shared)
        else
            Tdomain%sComm(i)%SF_nf_shared = 0
            Tdomain%sComm(i)%SF_ne_shared = 0
            Tdomain%sComm(i)%SF_nv_shared = 0
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            call read_attr_int(proc_id, "n_neu_edges", Tdomain%sComm(i)%Neu_ne_shared)
            call read_attr_int(proc_id, "n_neu_vertices", Tdomain%sComm(i)%Neu_nv_shared)
        else
            Tdomain%sComm(i)%Neu_ne_shared = 0
            Tdomain%sComm(i)%Neu_nv_shared = 0
        end if
        if(Tdomain%sComm(i)%nb_faces > 0)then
            allocate(Tdomain%sComm(i)%faces(0:Tdomain%sComm(i)%nb_faces-1))
            allocate(Tdomain%sComm(i)%orient_faces(0:Tdomain%sComm(i)%nb_faces-1))
            call read_dataset(proc_id, "faces", itemp)
            call read_dataset(proc_id, "faces_map", itempb)
            do j = 0,Tdomain%sComm(i)%nb_faces-1
                Tdomain%sComm(i)%faces(j)        = itemp(j+1)
                Tdomain%sComm(i)%orient_faces(j) = itempb(j+1)
            enddo
            deallocate(itemp, itempb)
        else
            nullify(Tdomain%sComm(i)%faces)
            nullify(Tdomain%sComm(i)%orient_faces)
        endif
        if(Tdomain%sComm(i)%nb_edges > 0)then
            call read_dataset(proc_id, "edges", itemp)
            call read_dataset(proc_id, "edges_map", itempb)
            allocate(Tdomain%sComm(i)%edges(0:Tdomain%sComm(i)%nb_edges-1))
            allocate(Tdomain%sComm(i)%orient_edges(0:Tdomain%sComm(i)%nb_edges-1))
            do j = 0,Tdomain%sComm(i)%nb_edges-1
                Tdomain%sComm(i)%edges(j) = itemp(j+1)
                Tdomain%sComm(i)%orient_edges(j) = itempb(j+1)
            enddo
            deallocate(itemp, itempb)
        else
            nullify(Tdomain%sComm(i)%edges)
            nullify(Tdomain%sComm(i)%orient_edges)
        endif
        if(Tdomain%sComm(i)%nb_vertices > 0)then
            call read_dataset(proc_id, "vertices", itemp)
            allocate(Tdomain%sComm(i)%vertices(0:Tdomain%sComm(i)%nb_vertices-1))
            do j = 0,Tdomain%sComm(i)%nb_vertices-1
                Tdomain%sComm(i)%vertices(j) = itemp(j+1)
            enddo
            deallocate(itemp)
        endif
        if(Tdomain%logicD%SF_local_present)then
            if(Tdomain%sComm(i)%SF_nf_shared > 0)then
                call read_dataset(proc_id, "sf_faces", itemp)
                allocate(Tdomain%sComm(i)%SF_faces_shared(0:Tdomain%sComm(i)%SF_nf_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_nf_shared-1
                    Tdomain%sComm(i)%SF_faces_shared(j) = itemp(j+1)
                enddo
                deallocate(itemp)
            else
                nullify(Tdomain%sComm(i)%SF_faces_shared)
            endif
            if(Tdomain%sComm(i)%SF_ne_shared > 0)then
                call read_dataset(proc_id,"sf_edges", itemp)
                call read_dataset(proc_id,"sf_edges_map", itempb)
                allocate(Tdomain%sComm(i)%SF_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                allocate(Tdomain%sComm(i)%SF_mapping_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_ne_shared-1
                    Tdomain%sComm(i)%SF_edges_shared(j) = itemp(j+1)
                    Tdomain%sComm(i)%SF_mapping_edges_shared(j) = itempb(j+1)
                enddo
                deallocate(itemp, itempb)
            else
                nullify(Tdomain%sComm(i)%SF_edges_shared)
                nullify(Tdomain%sComm(i)%SF_mapping_edges_shared)
            endif
            if(Tdomain%sComm(i)%SF_nv_shared > 0)then
                call read_dataset(proc_id, "sf_vertices", itemp)
                allocate(Tdomain%sComm(i)%SF_vertices_shared(0:Tdomain%sComm(i)%SF_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_nv_shared-1
                    Tdomain%sComm(i)%SF_vertices_shared(j) = itemp(j+1)
                enddo
                deallocate(itemp)
            else
                nullify(Tdomain%sComm(i)%SF_vertices_shared)
            endif
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            if(Tdomain%sComm(i)%Neu_ne_shared > 0)then
                call read_dataset(proc_id,"neu_edges", itemp)
                call read_dataset(proc_id,"neu_edges_map", itempb)
                allocate(Tdomain%sComm(i)%Neu_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                allocate(Tdomain%sComm(i)%Neu_mapping_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_ne_shared-1
                    Tdomain%sComm(i)%Neu_edges_shared(j) = itemp(j+1)
                    Tdomain%sComm(i)%Neu_mapping_edges_shared(j) = itempb(j+1)
                enddo
                deallocate(itemp, itempb)
            else
                nullify(Tdomain%sComm(i)%Neu_edges_shared)
                nullify(Tdomain%sComm(i)%Neu_mapping_edges_shared)
            endif
            if(Tdomain%sComm(i)%Neu_nv_shared > 0)then
                call read_dataset(proc_id,"neu_vertices", itemp)
                allocate(Tdomain%sComm(i)%Neu_vertices_shared(0:Tdomain%sComm(i)%Neu_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_nv_shared-1
                    Tdomain%sComm(i)%Neu_vertices_shared(j) = itemp(j+1)
                enddo
                deallocate(itemp)
            else
                nullify(Tdomain%sComm(i)%Neu_vertices_shared)
            endif
        else
            nullify(Tdomain%sComm(i)%Neu_edges_shared)
            nullify(Tdomain%sComm(i)%Neu_mapping_edges_shared)
            nullify(Tdomain%sComm(i)%Neu_vertices_shared)
        end if

    end do

    !write(*,*) "Mesh read correctly for proc #", rg
    if(rg == 0)then
        if(Tdomain%logicD%solid_fluid)then
            write(*,*) "  --> Propagation in solid-fluid media."
        else if(Tdomain%logicD%all_fluid)then
            write(*,*) "  --> Propagation in fluid media."
        else
            write(*,*) "  --> Propagation in solid media."
        end if
    end if
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
    !write(*,*) rg, "NFACES=", Tdomain%n_face

end subroutine read_mesh_file_h5


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
