!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module neumann

    ! anything that is related to interfaces on which Neumann BC are applied.

    use mesh_properties
    use sets

    implicit none

    type  :: Neu_face_type
       integer  :: ind_face
       integer, pointer  :: edges(:), vertex(:), orient_edge(:)
    end type Neu_face_type
    type  :: Neu_face
       integer  :: Elem,face, num
       integer  :: vertex(0:3), edge(0:3)
       !    integer  :: orient  ! orientation of the 1 face / in ref to the 0 one
    end type Neu_face
    type  :: Neu_edge
       integer  :: valence, coherency
       integer  :: vertex(0:1)
       integer, pointer  :: orient(:),faces(:),elem_ref(:,:)
       logical, pointer  :: local(:)
    end type Neu_edge
    type  :: Neu_vertex
       integer  :: node, valence
       integer, pointer  :: faces(:)
       logical, pointer  :: local(:)
    end type Neu_vertex
    type :: proc_obj
       integer, dimension(:,:), pointer  :: E,V
    end type proc_obj


contains




    !-----------------------------------------------------
    subroutine Neu_object_construct(n_elem,Ipointer,n_neu,Faces_Neumann,initnode,  &
        Neu_object,Neu_object_Face,Neu_object_n_faces)
        !- construction of Neumann objects: properties of Neumann interfaces
        implicit none
        integer, intent(in)   :: n_elem, n_neu
        integer, dimension(0:,0:), intent(in)  :: Ipointer, Faces_Neumann
        type(near_node), dimension(0:), intent(in)  ::  initnode
        logical, dimension(:), allocatable, intent(out)   ::  Neu_object
        type(Neu_Face_type), dimension(:,:), allocatable, intent(out)   ::  Neu_object_Face
        integer, dimension(:), allocatable, intent(out)   ::  Neu_object_n_faces

        integer   :: j,k,n,nf,nn,ok,num,node_ref
        integer, dimension(0:3)  :: corner
        type(near_entity), pointer  :: near_neighb => NULL()

        write(*,*) "  --> Construction of Neumann objects."
        allocate(Neu_object(0:n_elem-1),Neu_object_n_faces(0:n_elem-1))
        allocate(Neu_object_Face(0:n_elem-1,0:5))

        Neu_object(0:) = .false. ; Neu_object_n_faces(0:) = 0
        Neu_Object_face(0:,0:)%ind_face = -1

        ! loop on Neumann faces
        do n = 0,n_neu-1
            node_ref = Faces_Neumann(0,n)
            ! near elements, to shorten the search:
            near_neighb => initnode(node_ref)%ptr
            neu_search_0: do while(associated(near_neighb))
                nn = near_neighb%elem
                do nf = 0,5
                    call face2corner(Ipointer(0:,nn),nf,corner)
                    search_neu1 : do j = 0,3
                        num = Faces_Neumann(j,n)
                        ok = 0
                        search_neu2 : do k = 0,3
                            if(corner(k) == num)then
                                ok = 1
                                exit search_neu2
                            endif
                        end do search_neu2
                        if(ok == 0) exit search_neu1
                        if(j == 3) then  ! Ok, we found the element with the Neumann interface
                            Neu_Object(nn) = .true.
                            Neu_object_face(nn,Neu_object_n_faces(nn))%ind_face = nf
                            Neu_object_n_faces(nn) = Neu_object_n_faces(nn)+1
                            exit neu_search_0
                        endif
                    end do search_neu1
                end do
                near_neighb => near_neighb%pt   ! next element in the linked list
            end do neu_search_0

        end do

    end subroutine Neu_object_construct
    !------------------------------------------------------
    !------------------------------------------------------
    subroutine Neu_global_faces_construct(n_elem,Neu_n_global_faces,      &
        Neu_object,Neu_object_face,Neu_object_n_faces,  &
        Neu_global_faces)
        !- construction of GLOBAL Neumann faces
        implicit none
        integer, intent(in)       :: n_elem, Neu_n_global_faces
        logical, dimension(0:), intent(in)   ::  Neu_object
        type(Neu_Face_type), dimension(0:,0:), intent(in)   ::  Neu_object_Face
        integer, dimension(0:), intent(in)   ::  Neu_object_n_faces
        type(Neu_face), dimension(0:), intent(out) :: Neu_global_faces

        integer   :: n,nf,nnf

        write(*,*) "  --> Construction of global Neumann faces"
        write(*,*) "    --> Nb of global Neumann faces: ", Neu_n_global_faces
        nnf = 0
        do n = 0,n_elem-1
            if(Neu_object(n))then
                Neu_global_faces(nnf)%elem = -1
                do nf = 0, Neu_object_n_faces(n)-1
                    Neu_global_faces(nnf)%elem = n
                    Neu_global_faces(nnf)%face = Neu_object_Face(n,nf)%ind_face
                    ! absolute index in the element
                    Neu_global_faces(nnf)%num = nf
                    ! number of the face in the element
                    nnf = nnf+1    ! one more global face
                end do
            end if
        end do

    end subroutine Neu_global_faces_construct
    !------------------------------------------------------
    !------------------------------------------------------
    subroutine Neu_global_vertices_construct(Neu_n_global_faces,       &
        Ipointer,Neu_global_faces,Neu_global_node_to_vertex,    &
        Neu_n_global_vertices,Neu_global_vertices)

        implicit none
        integer, intent(in)   :: Neu_n_global_faces
        integer, dimension(0:,0:), intent(in) :: Ipointer
        type(Neu_face), dimension(0:Neu_n_global_faces-1), intent(inout) ::    &
            Neu_global_faces
        integer, dimension(0:), intent(out) :: Neu_global_node_to_vertex
        integer, intent(out)  :: Neu_n_global_vertices
        type(Neu_vertex), dimension(:), allocatable, intent(out)  :: Neu_global_vertices

        integer   :: i,j,k,nel,nf
        integer, dimension(0:3)   :: corner
        integer, dimension(0:4*Neu_n_global_faces-1)  :: node_list, node_valence
        integer, dimension(0:7,0:4*Neu_n_global_faces-1) :: Neu_vert_to_face

        write(*,*) "  --> Construction of global Neumann vertices"
        Neu_n_global_vertices = 0
        node_list(0:) = -1
        node_valence(0:) = 0
        Neu_global_node_to_vertex(0:) = -1

        do i = 0,Neu_n_global_faces-1
            nel = Neu_global_faces(i)%elem
            nf = Neu_global_faces(i)%face
            call face2corner(Ipointer(0:,nel),nf,corner)
            do j = 0,3
                k = already_in(corner(j),node_list,Neu_n_global_vertices)
                if(k < 0)then   ! node not seen
                    node_list(Neu_n_global_vertices) = corner(j)
                    node_valence(Neu_n_global_vertices) = 1
                    Neu_vert_to_face(0,Neu_n_global_vertices) = i
                    Neu_global_faces(i)%vertex(j) = Neu_n_global_vertices
                    Neu_n_global_vertices = Neu_n_global_vertices+1
                else   !  node already seen
                    Neu_vert_to_face(node_valence(k),k) = i
                    Neu_global_faces(i)%vertex(j) = k
                    node_valence(k) = node_valence(k)+1
                end if
            end do
        end do
        write(*,*) "    --> Nb of global Neumann vertices: ", Neu_n_global_vertices
        allocate(Neu_global_vertices(0:Neu_n_global_vertices-1))
        do i = 0,Neu_n_global_vertices-1
            Neu_global_vertices(i)%valence = node_valence(i)
            ! nb of faces seen by a global Neumann vertex
            Neu_global_vertices(i)%node = node_list(i)
            Neu_global_node_to_vertex(Neu_global_vertices(i)%node) = i
            !- now we link vertices to faces
            allocate(Neu_global_vertices(i)%faces(0:Neu_global_vertices(i)%valence-1))
            Neu_global_vertices(i)%faces(0:Neu_global_vertices(i)%valence-1) =    &
                Neu_vert_to_face(0:Neu_global_vertices(i)%valence-1,i)
        end do

    end subroutine Neu_global_vertices_construct
    !------------------------------------------------------
    !------------------------------------------------------
    subroutine Neu_global_edges_construct(Neu_n_global_faces,Ipointer,  &
        Neu_global_node_to_vertex, Neu_global_faces,          &
        Neu_n_global_edges,Neu_global_edges)

        implicit none
        integer, intent(in)      :: Neu_n_global_faces
        integer, dimension(0:), intent(in)   :: Neu_global_node_to_vertex
        integer, dimension(0:,0:), intent(in)   :: Ipointer
        type(Neu_face), dimension(0:Neu_n_global_faces-1)  :: Neu_global_faces
        integer, intent(out)   :: Neu_n_global_edges
        type(Neu_edge), dimension(:), allocatable, intent(out)  :: Neu_global_edges

        integer, dimension(0:4*Neu_n_global_faces-1)  :: edge_valence
        integer, dimension(0:7,0:4*Neu_n_global_faces-1) :: Neu_edge_to_face,  &
            Neu_edge_orient_tmp
        integer, dimension(0:4*Neu_n_global_faces-1,0:1) :: Neu_vertex_to_edge
        integer  :: i,j,n,nf,ne,nv0,nv1,ok
        integer, dimension(0:1)  :: corner_edge, e_corner, e_vertex
        integer, dimension(0:3)  :: corner, edge

        write(*,*) "  --> Construction of global Neumann edges"
        Neu_n_global_edges = 0

        do i = 0,Neu_n_global_faces-1
            n = Neu_global_faces(i)%elem
            nf = Neu_global_faces(i)%face
            call caract_face(nf,corner,edge)
            do j = 0,3
                corner_edge(0) = corner(j)
                if(j == 2) then
                    corner_edge(0) = corner(3)
                    corner_edge(1) = corner(2)
                else if (j == 3) then
                    corner_edge(0) = corner(0)
                    corner_edge(1) = corner(3)
                else
                    corner_edge(1) = corner(j+1)
                endif
                !- global node..
                e_corner(0) = Ipointer(corner_edge(0),n)
                e_corner(1) = Ipointer(corner_edge(1),n)
                !- .. and Neumann vertices associated:
                e_vertex(0) = Neu_global_node_to_vertex(e_corner(0))
                e_vertex(1) = Neu_global_node_to_vertex(e_corner(1))
                ok = 0
                findNeuedge: do ne = 0,Neu_n_global_edges-1  ! edge already seen?
                    nv0 = Neu_vertex_to_edge(ne,0)
                    nv1 = Neu_vertex_to_edge(ne,1)
                    if((nv0 == e_vertex(0) .and. nv1 == e_vertex(1)) .or. &
                        (nv0 == e_vertex(1) .and. nv1 == e_vertex(0)))then
                        ! yes, Neumann edge seen before
                        ok = 1
                        Neu_edge_to_face(edge_valence(ne),ne) = i
                        Neu_global_faces(i)%edge(j) = ne
                        edge_valence(ne) = edge_valence(ne)+1
                        exit findNeuedge
                    end if
                end do findNeuedge
                if(ok == 0)then    ! Neumann edge never seen
                    Neu_edge_to_face(0,Neu_n_global_edges) = i
                    Neu_vertex_to_edge(Neu_n_global_edges,0:1) = e_vertex(0:1)
                    Neu_global_faces(i)%edge(j) = Neu_n_global_edges
                    Neu_edge_orient_tmp(0,Neu_n_global_edges) = 0
                    edge_valence(Neu_n_global_edges) = 1
                    Neu_n_global_edges = Neu_n_global_edges+1
                end if
            end do
        end do
        write(*,*) "    --> Nb of global Neumann edges: ", Neu_n_global_edges
        allocate(Neu_global_edges(0:Neu_n_global_edges-1))
        do i = 0,Neu_n_global_edges-1
            Neu_global_edges(i)%vertex(0:1) = Neu_vertex_to_edge(i,0:1)
            Neu_global_edges(i)%valence = edge_valence(i)
            allocate(Neu_global_edges(i)%faces(0:edge_valence(i)-1),       &
                Neu_global_edges(i)%orient(0:edge_valence(i)-1))
            Neu_global_edges(i)%faces(0:edge_valence(i)-1) =        &
                Neu_edge_to_face(0:edge_valence(i)-1,i)
            Neu_global_edges(i)%orient(0:edge_valence(i)-1) =              &
                Neu_edge_orient_tmp(0:edge_valence(i)-1,i)
        end do

    end subroutine Neu_global_edges_construct
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Neu_vertices_proc_belong(n_proc,Neu_n_global_vertices,part,    &
        initnode,Neu_global_vertices)

        implicit none
        integer, intent(in)     :: n_proc,Neu_n_global_vertices
        integer, dimension(0:n_proc-1), intent(in)  :: part
        type(near_node), dimension(0:), intent(in)  ::  initnode
        type(Neu_vertex), dimension(0:Neu_n_global_vertices-1), intent(inout) ::  &
            Neu_global_vertices

        integer   :: i,nf
        type(near_entity), pointer  :: near_neighb => NULL()

        do i = 0, Neu_n_global_vertices-1
            allocate(Neu_global_vertices(i)%local(0:n_proc-1))
            Neu_global_vertices(i)%local(0:) = .false.
            nf = Neu_global_vertices(i)%node
            near_neighb => initnode(nf)%ptr
            do while(associated(near_neighb))
                Neu_global_vertices(i)%local(part(near_neighb%elem)) = .true.
                near_neighb => near_neighb%pt
            end do
        end do


    end subroutine Neu_vertices_proc_belong
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Neu_edges_proc_belong(n_proc,Neu_n_global_edges,     &
        Neu_global_vertices,Neu_global_edges)

        implicit none
        integer, intent(in)     :: n_proc,Neu_n_global_edges
        type(Neu_vertex), dimension(0:), intent(in)   :: Neu_global_vertices
        type(Neu_edge), dimension(0:Neu_n_global_edges-1), intent(inout) ::   &
            Neu_global_edges

        integer   :: i,n,nv0,nv1

        do i = 0,Neu_n_global_edges-1
            allocate(Neu_global_edges(i)%local(0:n_proc-1))
            Neu_global_edges(i)%local(0:) = .false.
            nv0 = Neu_global_edges(i)%vertex(0)
            nv1 = Neu_global_edges(i)%vertex(1)
            do n = 0,n_proc-1
                if(Neu_global_vertices(nv0)%local(n) .and.   &
                    Neu_global_vertices(nv1)%local(n))       &
                    Neu_global_edges(i)%local(n) = .true.
            end do
        end do


    end subroutine Neu_edges_proc_belong
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Neu_edges_reference_elem(nproc,Neu_n_global_edges,       &
        Elem_glob2loc,part,Ipointer,Neu_global_vertices,   &
        Neu_edges_near_elem,Neu_global_edges)

        implicit none
        integer, intent(in)   :: nproc,Neu_n_global_edges
        integer, dimension(0:), intent(in)   :: Elem_glob2loc, part
        integer, dimension(0:,0:), intent(in)   :: Ipointer
        type(Neu_vertex), dimension(0:), intent(in)   :: Neu_global_vertices
        type(near_elem), dimension(0:Neu_n_global_edges-1), intent(in)   ::   &
            Neu_edges_near_elem
        type(Neu_edge), dimension(0:Neu_n_global_edges-1), intent(inout)  ::   &
            Neu_global_edges

        logical, dimension(:), allocatable   :: Lproc
        integer   :: i,j,k,nv0,nv1,n,np,nel,num,ok
        integer, dimension(0:1)  :: npf,e_neighbor_corner
        type(near_entity), pointer   :: near_neighb => NULL()

        allocate(Lproc(0:nproc-1))
        do i = 0,Neu_n_global_edges-1
            Lproc = .true.
            allocate(Neu_global_edges(i)%elem_ref(0:nproc-1,0:1))
            Neu_global_edges(i)%elem_ref = -1

            nv0 = Neu_global_edges(i)%vertex(0)
            nv1 = Neu_global_edges(i)%vertex(1)
            npf(0) = Neu_global_vertices(nv0)%node
            npf(1) = Neu_global_vertices(nv1)%node
            near_neighb => Neu_edges_near_elem(i)%ptr

            do while(associated(near_neighb))
                nel = near_neighb%elem
                n = Elem_glob2loc(nel)
                np = part(nel)
                if(Lproc(np).and. Neu_global_edges(i)%local(np))then
                    do j = 0,1
                        num = npf(j)
                        ok = 0
                        do_k: do k = 0,7
                            if(num == Ipointer(k,nel))then
                                e_neighbor_corner(j) = k
                                ok = 1
                                exit do_k
                            end if
                        end do do_k
                    end do
                    if(j == 2 .and. ok == 1)then
                        Neu_global_edges(i)%elem_ref(np,0) = n
                        Neu_global_edges(i)%elem_ref(np,1) =         &
                            vertices2edge(e_neighbor_corner)
                        Lproc(np) = .false.
                    end if
                end if
                near_neighb => near_neighb%pt
            end do
        end do

        deallocate(Lproc)


    end subroutine Neu_edges_reference_elem
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Neu_local_faces_construct(proc,nproc,Neu_n_faces,          &
        Neu_n_global_faces,part,Elem_glob2loc,                 &
        Neu_global_face_present,Neu_global_faces,              &
        faces,Neu_faces,Neu_local_to_global_faces)

        implicit none
        integer, intent(in) :: proc,nproc,Neu_n_faces,Neu_n_global_faces
        logical, dimension(0:), intent(in)  :: Neu_global_face_present
        integer, dimension(0:), intent(in)  :: part,Elem_glob2loc
        type(Neu_face), dimension(0:), intent(in)  :: Neu_global_faces
        integer, dimension(0:,0:), intent(in)   :: faces
        integer, dimension(0:Neu_n_faces-1), intent(out)  :: Neu_faces
        integer, dimension(0:Neu_n_faces-1), intent(out) ::             &
            Neu_local_to_global_faces
        integer :: j,nf,nel,nnf,n

        j = 0    ! local Neumann faces counting
        do nf = 0,Neu_n_global_faces-1
            if(.not. Neu_global_face_present(nf)) cycle
            Neu_local_to_global_faces(j) = nf
            nel = Neu_global_faces(nf)%elem   ! global number of the element sharing the face
            nnf = Neu_global_faces(nf)%face   ! index of the face
            n = Elem_glob2loc(nel)
            Neu_faces(j) = faces(n,nnf)
            j = j+1    ! one more local Neumann face
        end do

    end subroutine Neu_local_faces_construct
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Neu_local_edges_construct(proc,nproc,Neu_n_edges,Neu_n_global_edges,    &
        part,glob2loc,Neu_global_edges,Neu_global_vertices,edges,              &
        Ipointer_local,Ipointer,which_elem_in_proc,Neu_global_to_local_edges,  &
        Neu_edges,Neu_edges_shared, Neu_mapping_edges_shared,Neu_ne_shared,    &
        MemoryNeu)

        implicit none
        integer, intent(in)  :: proc,nproc,Neu_n_edges,Neu_n_global_edges
        integer, dimension(0:), intent(in)  :: part,glob2loc
        type(Neu_edge), dimension(0:), intent(in)  :: Neu_global_edges
        type(Neu_vertex), dimension(0:), intent(in)  :: Neu_global_vertices
        integer, dimension(0:,0:), intent(in)   :: edges,Ipointer_local,   &
            Ipointer,which_elem_in_proc
        integer, dimension(0:), intent(out)  :: Neu_global_to_local_edges
        integer, dimension(0:Neu_n_edges-1), intent(out)  :: Neu_edges
        integer, dimension(0:,0:), intent(out)   :: Neu_edges_shared,   &
            Neu_mapping_edges_shared
        integer, dimension(0:nproc-1), intent(out)   :: Neu_ne_shared
        type(proc_obj), dimension(0:), intent(inout)  :: MemoryNeu
        integer :: i,j,nnf,num,n,ns,nes,ne,n0,nns,nv0,nels
        logical   ::  orient_proc,orient_o_proc

        Neu_edges(0:) = -1
        Neu_ne_shared(0:) = 0
        Neu_mapping_edges_shared(0:,0:) = -1

        j = 0   ! local Neumann edges counting
        do i = 0,Neu_n_global_edges-1
            if((.not. Neu_global_edges(i)%local(proc))) cycle
            !- yes, Neumann edge:
            orient_proc = .false.
            Neu_global_to_local_edges(i) = j
            n = Neu_global_edges(i)%elem_ref(proc,0)
            ne = Neu_global_edges(i)%elem_ref(proc,1)
            Neu_edges(j) = edges(n,ne)
            ! reference global vertex:
            nv0 = Neu_global_edges(i)%vertex(0)
            nnf = Neu_global_vertices(nv0)%node
            n0 = glob2loc(nnf)
            if(n0 == Ipointer_local(edge2vertex(ne),n)) orient_proc = .true.

            !- now we look for Neumann edges shared with other procs
            do num = 0,nproc-1
                if(num == proc) cycle
                if((.not. Neu_global_edges(i)%local(num))) cycle
                orient_o_proc = .false.
                ! now we have found a proc on which we have the same global Neumann edge
                ! interproc orientation
                ns = Neu_global_edges(i)%elem_ref(num,0)
                nels = which_elem_in_proc(num,ns)
                nes = Neu_global_edges(i)%elem_ref(num,1)
                nv0 = Neu_global_edges(i)%vertex(0)
                nns = Neu_global_vertices(nv0)%node
                if(nns == Ipointer(edge2vertex(nes),nels)) orient_o_proc = .true.

                if(num < proc)then   ! proc already seen
                    Neu_edges_shared(num,MemoryNeu(proc)%E(num,i)) = j
                    Neu_mapping_edges_shared(num,MemoryNeu(proc)%E(num,i)) =     &
                        merge(0,1,orient_proc .eqv. orient_o_proc)
                else   ! proc not seen yet
                    Neu_edges_shared(num,Neu_ne_shared(num)) = j
                    MemoryNeu(num)%E(proc,i) = Neu_ne_shared(num)
                    Neu_mapping_edges_shared(num,Neu_ne_shared(num)) =     &
                        merge(0,1,orient_proc .eqv. orient_o_proc)
                end if
                Neu_ne_shared(num) = Neu_ne_shared(num)+1
            end do
            !--
            j = j+1   ! one more local Neumann edge
        end do

    end subroutine Neu_local_edges_construct
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Neu_local_vertices_construct(proc,nproc,Neu_n_vertices,     &
        Neu_n_global_vertices,glob2loc,N_valid_Vertex,            &
        Neu_global_vertices,Neu_global_to_local_vertices,         &
        Neu_vertices,Neu_vertices_shared,Neu_nv_shared,MemoryNeu)

        implicit none
        integer, intent(in)  :: proc,nproc,Neu_n_vertices,Neu_n_global_vertices
        integer, dimension(0:), intent(in)  :: glob2loc,N_valid_Vertex
        type(Neu_vertex), dimension(0:), intent(in)  :: Neu_global_vertices
        integer, dimension(0:), intent(out)  :: Neu_global_to_local_vertices
        integer, dimension(0:Neu_n_vertices-1), intent(out) :: Neu_vertices
        integer, dimension(0:,0:), intent(out)   :: Neu_vertices_shared
        integer, dimension(0:nproc-1), intent(out)   :: Neu_nv_shared
        type(proc_obj), dimension(0:), intent(inout)  :: MemoryNeu
        integer :: i,j,nf,nnf,num

        Neu_vertices(0:) = -1
        Neu_nv_shared(0:) = 0

        j = 0   ! local Neu vertices counting
        do i = 0,Neu_n_global_vertices-1
            if((.not. Neu_global_vertices(i)%local(proc))) cycle
            Neu_global_to_local_vertices(i) = j
            nf = Neu_global_vertices(i)%node
            nnf = glob2loc(nf)
            Neu_vertices(j) = N_valid_Vertex(nnf)

            !- now we look for Neumann vertices shared with other procs
            do num = 0,nproc-1
                if(num == proc) cycle
                if((.not. Neu_global_vertices(i)%local(num))) cycle
                ! now we have found a proc on which we have the same global Neumann vertex
                if(num < proc)then   ! proc already seen
                    Neu_vertices_shared(num,MemoryNeu(proc)%V(num,i)) = j
                else   ! proc not seen yet
                    Neu_vertices_shared(num,Neu_nv_shared(num)) = j
                    MemoryNeu(num)%V(proc,i) = Neu_nv_shared(num)
                end if
                Neu_nv_shared(num) = Neu_nv_shared(num)+1
            end do

            j = j+1   ! one more local SF vertex
        end do


    end subroutine Neu_local_vertices_construct
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    subroutine Neu_faces_to_EV(proc,Neu_n_faces,Neu_local_to_global_faces,    &
        Neu_global_to_local_edges,Neu_global_to_local_vertices,         &
        Neu_global_faces,Neu_global_edges,Neu_faces,mapping_edges,      &
        Neu_Face_Near_Edges,Neu_Face_Near_Edges_Orient,Neu_Face_Near_Vertices)

        implicit none
        integer, intent(in)  :: proc,Neu_n_faces
        integer, dimension(0:), intent(in)  :: Neu_local_to_global_faces,Neu_Faces, &
            Neu_global_to_local_edges,Neu_global_to_local_vertices
        type(Neu_face), dimension(0:), intent(in)  :: Neu_global_faces
        type(Neu_edge), dimension(0:), intent(in)  :: Neu_global_edges
        integer, dimension(0:,0:), intent(in)  :: mapping_edges
        integer, dimension(0:,0:), intent(out)  :: Neu_Face_Near_Edges,   &
            Neu_Face_Near_Edges_Orient,Neu_Face_Near_Vertices
        integer   ::  i,j,k,n,ne,nel,nes,nv,nvl

        Neu_Face_Near_Edges(0:,0:) = -1
        Neu_Face_Near_Edges_Orient(0:,0) = -1
        Neu_Face_Near_Vertices(0:,0:) = -1

        do j = 0, Neu_n_faces-1
            i = Neu_local_to_global_faces(j)
            do k = 0,3
                ! edges
                ne = Neu_global_faces(i)%edge(k)
                nel = Neu_global_to_local_edges(ne)
                Neu_Face_Near_Edges(j,k) = nel
                n = Neu_global_edges(ne)%elem_ref(proc,0)
                nes = Neu_global_edges(ne)%elem_ref(proc,1)
                Neu_Face_Near_Edges_Orient(j,k) = mapping_edges(n,nes)
                ! vertices
                nv = Neu_global_faces(i)%vertex(k)
                nvl = Neu_global_to_local_vertices(nv)
                Neu_Face_Near_Vertices(j,k) = nvl
            end do

        end do

    end subroutine Neu_faces_to_EV
    !--------------------------------------------------------------
    !--------------------------------------------------------------

end module neumann

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
