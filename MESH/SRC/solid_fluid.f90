module solid_fluid

    ! all that is related to the solid-fluid interfaces

    use mesh_properties
    use sets

    implicit none


    type  :: Face_type
       integer  :: ind_face, neighbor_elem, face_neighbor, orient
    end type Face_type
    type  :: SF_face
       integer  :: Elem(0:1),face(0:1), num(0:1)   !- 0: fluid, 1: solid
       integer  :: vertex(0:3), edge(0:3)
       integer  :: orient  ! orientation of the 1 face / in ref to the 0 one
    end type SF_face
    type  :: SF_edge
       integer  :: valence, coherency
       integer  :: vertex(0:1)
       integer, pointer  :: orient(:), faces(:), elem_ref_fluid(:,:),   &
           elem_ref_solid(:,:)
       logical, pointer  :: local_fluid(:), local_solid(:)
    end type SF_edge
    type  :: SF_vertex
       integer  :: node(0:1), valence
       integer, pointer  :: faces(:)
       logical, pointer  :: local_fluid(:), local_solid(:)
    end type SF_vertex
    type :: process_obj
       integer, dimension(:,:), pointer  :: F,E,V
    end type process_obj


contains


    !-------------------------------------
    subroutine solid_fluid_init(n_nods, n_points, n_blocks, n_elem, tabmat,  &
        Ipointer, Material,any_fluid,all_fluid,solid_fluid,       &
        n_SF_nodes,nodes_nature,elem_solid,elem_contact)
        ! gives general pieces of info concerning solid/fluid properties
        implicit none
        integer, intent(in)  :: n_nods, n_points, n_blocks, n_elem
        character, dimension(0:n_blocks-1), intent(in)  :: tabmat
        integer, dimension(0:n_elem-1), intent(in)    :: Material
        integer, dimension(0:n_nods-1,0:n_elem-1), intent(in) :: Ipointer
        logical, intent(out)  :: any_fluid, all_fluid, solid_fluid
        integer, intent(out)  :: n_SF_nodes
        integer, dimension(:), allocatable, intent(out) :: nodes_nature
        logical, dimension(:), allocatable, intent(out) :: elem_solid, elem_contact

        write(*,*)
        write(*,*) "****************************************"
        write(*,*) "  Fluid case - if any : generalities"
        write(*,*) "****************************************"

        allocate(nodes_nature(0:n_points-1))
        allocate(elem_solid(0:n_elem-1),elem_contact(0:n_elem-1))
        elem_solid(0:) = .true.
        elem_contact(0:) = .false.

        call solid_fluid_bool(n_blocks,tabmat,any_fluid,all_fluid,solid_fluid)

        !- nature of elements (S,F,contact?) and of nodes (S,F,SF?)
        write(*,*) "  --> Nature of elements (S,F,contact?) and of nodes (S,F,SF?)"
        call nature_node_elem(tabmat,Ipointer,Material,n_nods,any_fluid,all_fluid,   &
            solid_fluid,nodes_nature,elem_solid,elem_contact,n_SF_nodes)

    end subroutine solid_fluid_init
    !------------------------------------------
    !------------------------------------------
    subroutine solid_fluid_bool(nmat,mattab,any_f,all_f,s_f)
        !- checking if fluid in the domain: any, or all, or nothing.
        implicit none
        integer, intent(in)   :: nmat
        character, dimension(0:), intent(in) :: mattab
        logical, intent(out)  :: any_f, all_f,s_f
        integer    :: i


        if (size(mattab) /= nmat)  &
            stop "In solid_fluid_bool: incompatibility for # of materials"
        ! init.
        any_f = .false. ; all_f = .true. ; s_f = .false.
        do i = 0, nmat-1
            if(mattab(i) == 'F' .or. mattab(i) == 'L') any_f = .true.
            if(mattab(i) == 'S' .or. mattab(i) == 'P') all_f = .false.
        end do

        !-
        if(any_f)then
            if(all_f)then
                print*,"  --> Domain entirely fluid."
            else
                s_f = .true.
                print*,"  --> Solid/fluid domain."
            end if
        else
            print*,"  --> Domain entirely solid."
        end if

        return
    end subroutine solid_fluid_bool
    !---------------------------------------
    subroutine nature_node_elem(mattab,Ipointer,mater,n_nodes,all_f,any_f,  &
        s_f,nods,elems,elemc,nbSF)
        implicit none
        integer   :: n,i,j,code,imat,nn
        character, intent(in) :: mattab(0:)
        integer, dimension(0:,0:), intent(in)  :: Ipointer
        integer, intent(in)   :: mater(0:), n_nodes
        logical, intent(in)   :: all_f,any_f,s_f
        integer, intent(out)  :: nods(0:), nbSF
        logical, intent(out)  :: elems(0:),elemc(0:)

        !- code for a node: 0 if fluid, 1 if solid, 2 if both
        if(.not. any_f) then  ! pure solid
            nods(0:) = 1
            elems(0:) = .true.
        end if
        if(all_f) then   ! pure fluid
            nods(0:) = 0
            elems(0:) = .false.
        end if
        if(s_f)then    ! solid-fluid case
            elems(:) = .false.
            nods(:) = -1
            !- Nature of nodes
            do n = 0, size(elems)-1
                imat = mater(n)
                if(mattab(imat) == 'S' .or. mattab(imat) == 'P')then
                    code = 1   ! solid
                else
                    code = 0   ! liquid
                end if
                !-
                do i = 0, n_nodes-1
                    j = Ipointer(i,n)
                    if(nods(j) == code .or. nods(j) == -1)then
                        nods(j) = code
                    else
                        nods(j) = 2   ! solid-fluid node
                    end if
                    if(i == 0)then
                        if(code == 1)then
                            elems(n) = .true.  ! solid element
                        else
                            elems(n) = .false.  ! fluid element
                        end if
                    else     ! checking for consistency
                        if((elems(n) .and. code == 0).or.(.not.elems(n) .and. code == 1))   &
                            stop " Pb in nature_node_elem: an element cannot be solid and fluid  &
                            &  at the same time"
                    end if
                end do
            end do
            ! elements in SF contact (or not)
            do n = 0, size(elems)-1
                do i = 0, n_nodes-1
                    j = Ipointer(i,n)
                    if(nods(j) == 2) elemc(n) = .true.
                end do
            end do
        end if

        ! solid-fluid nodes counting
        nn = 0
        do n = 0, size(nods)-1
            if(nods(n) == -1)   &
                stop "Pb in solid-fluid nodes attribution: one node is unknown"
            if(nods(n) == 2) nn = nn+1
        end do
        nbSF = nn
        write(*,*) "    --> Nb of solid-fluid nodes: ",nbSF

        return
    end subroutine nature_node_elem
    !-----------------------------------------------------
    !----------------------------------------
    subroutine SF_nodes_doubling(Nodes_SF,nods,npts)
        implicit none
        integer, intent(out)   :: Nodes_SF(0:,0:)
        integer, intent(in)    :: nods(0:),npts
        integer                :: i,j

        j = 0

        do i = 0,npts-1
            if(nods(i) == 2)then    ! solid-fluid node
                Nodes_SF(j,0) = i
                Nodes_SF(j,1) = npts+j    ! new node added
                j = j+1
            end if
        end do

        if(j /= size(Nodes_SF,1)) stop "Problem in SF_nodes_doubling."

        return
    end subroutine SF_nodes_doubling
    !----------------------------------------
    !-----------------------------------------------------
    subroutine SF_object_construct(n_elem,Ipointer,elem_contact, elem_solid,  &
        nodes_nature,dxadj,dxadjncy,SF_object,     &
        SF_object_Face,SF_object_n_faces,          &
        SF_n_global_faces)
        !- construction of SF objects: properties of solid-fluid interfaces
        implicit none
        integer, intent(in)   :: n_elem
        integer, dimension(0:,0:), intent(in)  :: Ipointer
        logical, dimension(0:n_elem-1), intent(in)  :: elem_contact,elem_solid
        integer, dimension(0:), intent(in)  :: nodes_nature, dxadj, dxadjncy
        logical, dimension(:,:), allocatable, intent(out)   ::  SF_object
        type(Face_type), dimension(:,:), allocatable, intent(out)   ::  SF_object_Face
        integer, dimension(:), allocatable, intent(out)   ::  SF_object_n_faces
        integer, intent(out)   :: SF_n_global_faces

        integer   :: i,j,k,l,n,nf,nn,nnn,ok,num,neighbor_edge,neighbor_face
        integer, dimension(0:1)  :: corner_edge
        integer, dimension(0:3)  :: corner, corner_fl, neighbor_corner,  &
            corner_s,edge_fl,edge_s

        SF_n_global_faces = 0
        write(*,*) "  --> Construction of Solid-Fluid objects."
        allocate(SF_object(0:1,0:n_elem-1),SF_object_n_faces(0:n_elem-1))
        allocate(SF_object_Face(0:n_elem-1,0:5))

        SF_object(0:,0:) = .false. ; SF_object_n_faces(0:) = 0
        SF_object_face(0:,0:)%orient = -1

        do n = 0, n_elem-1
            if(elem_contact(n) .and. elem_solid(n)) then   ! a solid element,
                ! we look for fluid neighbours
                do nf = 0,5
                    call face2corner(Ipointer(0:,n),nf,corner)
                    !--  fluid neighbours????
                    if(nodes_nature(corner(0)) == 2 .and. nodes_nature(corner(1)) == 2 .and. &
                        nodes_nature(corner(2)) == 2 .and. nodes_nature(corner(3)) == 2)then
                        ! Yes: SF face. Now we define it.
                        search0 : do nnn = dxadj(n),dxadj(n+1)-1
                            nn = dxadjncy(nnn)
                            if(elem_contact(nn) .and. .not.elem_solid(nn)) then
                                search1 : do j = 0,3   ! Does the face belong to this element?
                                    num = corner(j)
                                    ok = 0
                                    search2 : do k = 0,7
                                        if (Ipointer(k,nn) == num) then
                                            neighbor_corner(j) = k
                                            ok = 1
                                            exit search2
                                        endif
                                    enddo search2
                                    if(ok == 0) exit search1   ! No, so let's see another element
                                    if(j == 3) then   ! Yes
                                        corner_fl(0:3) = neighbor_corner(0:3)
                                        call sort(neighbor_corner,4)
                                        ! So which face of the neighbor is it ?
                                        neighbor_face = neighb_face(neighbor_corner)

                                        SF_object(0,nn) = .true.  ! fluid element
                                        SF_object(1,n) = .true.   ! solid element
                                        SF_object_face(n,SF_object_n_faces(n))%ind_face = nf
                                        SF_object_face(nn,SF_object_n_faces(nn))%ind_face = neighbor_face
                                        SF_object_face(n,SF_object_n_faces(n))%face_neighbor = neighbor_face
                                        SF_object_face(nn,SF_object_n_faces(nn))%face_neighbor = nf
                                        SF_object_face(n,SF_object_n_faces(n))%neighbor_elem = nn
                                        SF_object_face(nn,SF_object_n_faces(nn))%neighbor_elem = n
                                        !- orientation of a face with reference to the other one (symmetry)
                                        call face_orientation(corner_fl,4,neighbor_face,       &
                                            SF_object_face(nn,SF_object_n_faces(nn))%orient)
                                        SF_object_face(n,SF_object_n_faces(n))%orient  =   &
                                            SF_object_face(nn,SF_object_n_faces(nn))%orient
                                        !-- increment of the nb of SF faces for elements in (face)contact
                                        SF_object_n_faces(n) = SF_object_n_faces(n)+1
                                        SF_object_n_faces(nn) = SF_object_n_faces(nn)+1
                                        SF_n_global_faces = SF_n_global_faces+1
                                    end if
                                end do search1
                            end if
                        enddo search0
                    end if
                    !-------------------
                    !--------
                end do
            end if
        end do

    end subroutine SF_object_construct
    !------------------------------------------------------
    subroutine SF_global_faces_construct(n_elem, SF_n_global_faces,      &
        SF_object,SF_object_face,SF_object_n_faces,  &
        SF_global_faces)
        !- construction of GLOBAL solid-fluid faces
        implicit none
        integer, intent(in)       :: n_elem, SF_n_global_faces
        logical, dimension(0:,0:), intent(in)   ::  SF_object
        type(Face_type), dimension(0:,0:), intent(in)   ::  SF_object_Face
        integer, dimension(0:), intent(in)   ::  SF_object_n_faces
        type(SF_face), dimension(0:), intent(out) :: SF_global_faces

        integer   :: n,nf,nnf,i,nels

        write(*,*) "  --> Construction of global solid/fluid faces"
        write(*,*) "    --> Nb of global solid/fluid faces: ", SF_n_global_faces
        nnf = 0
        do n = 0,n_elem-1
            if(SF_object(0,n))then   ! look only at fluid elements:
                ! as we are in the global frame, it's ok
                SF_global_faces(nnf)%elem(:) = -1
                do nf = 0, SF_object_n_faces(n)-1
                    SF_global_faces(nnf)%elem(0) = n
                    SF_global_faces(nnf)%face(0) = SF_object_Face(n,nf)%ind_face
                    ! absolute index in the element
                    SF_global_faces(nnf)%num(0) = nf
                    ! number of the SF face in the fluid element
                    nels = SF_object_Face(n,nf)%neighbor_elem
                    SF_global_faces(nnf)%elem(1) = nels
                    SF_global_faces(nnf)%face(1) = SF_object_face(n,nf)%face_neighbor
                    do i = 0,SF_object_n_faces(nels)-1
                        if(SF_object_face(nels,i)%ind_face == SF_global_faces(nnf)%face(1)) exit
                    end do
                    SF_global_Faces(nnf)%num(1) = i
                    SF_global_faces(nnf)%orient = SF_object_face(n,nf)%orient
                    nnf = nnf+1    ! one more global face
                end do
            end if
        end do

    end subroutine SF_global_faces_construct
    !------------------------------------------------------
    !------------------------------------------------------
    subroutine SF_global_vertices_construct(n_SF_nodes,SF_n_global_faces,       &
        Ipointer,nodes_on_SF,SF_global_faces,SF_global_node_to_vertex,   &
        SF_n_global_vertices,SF_global_vertices)

        implicit none
        integer, intent(in)   :: n_SF_nodes, SF_n_global_faces
        integer, dimension(0:,0:), intent(in) :: Ipointer
        integer, dimension(0:n_SF_nodes-1,0:1), intent(in) :: nodes_on_SF
        type(SF_face), dimension(0:SF_n_global_faces-1), intent(inout) ::    &
            SF_global_faces
        integer, dimension(0:), intent(out) :: SF_global_node_to_vertex
        integer, intent(out)  :: SF_n_global_vertices
        type(SF_vertex), dimension(:), allocatable, intent(out)  :: SF_global_vertices

        integer   :: i,j,k,nel,nf
        integer, dimension(0:3)   :: corner
        integer, dimension(0:4*SF_n_global_faces-1)  :: node_list, node_valence
        integer, dimension(0:7,0:4*SF_n_global_faces-1) :: SF_vert_to_face

        write(*,*) "  --> Construction of global solid/fluid vertices"
        SF_n_global_vertices = 0
        node_list(0:) = -1
        node_valence(0:) = 0
        SF_global_node_to_vertex(0:) = -1

        do i = 0,SF_n_global_faces-1
            nel = SF_global_faces(i)%elem(0)
            nf = SF_global_faces(i)%face(0)
            call face2corner(Ipointer(0:,nel),nf,corner)
            do j = 0,3
                k = already_in(corner(j),node_list,SF_n_global_vertices)
                if(k < 0)then   ! node not seen
                    node_list(SF_n_global_vertices) = corner(j)
                    node_valence(SF_n_global_vertices) = 1
                    SF_vert_to_face(0,SF_n_global_vertices) = i
                    SF_global_faces(i)%vertex(j) = SF_n_global_vertices
                    SF_n_global_vertices = SF_n_global_vertices+1
                else   !  node already seen
                    SF_vert_to_face(node_valence(k),k) = i
                    SF_global_faces(i)%vertex(j) = k
                    node_valence(k) = node_valence(k)+1
                end if
            end do
        end do
        write(*,*) "    --> Nb of global solid/fluid vertices: ", SF_n_global_vertices
        allocate(SF_global_vertices(0:SF_n_global_vertices-1))
        do i = 0,SF_n_global_vertices-1
            SF_global_vertices(i)%valence = node_valence(i)
            ! nb of faces seen by a global SF vertex
            SF_global_vertices(i)%node(1) = node_list(i)
            ! old node is a solid one
            do j = 0,n_SF_nodes-1
                if(node_list(i) == Nodes_on_SF(j,0)) SF_global_vertices(i)%node(0) =   &
                    Nodes_on_SF(j,1)
            end do
            SF_global_node_to_vertex(SF_global_vertices(i)%node(0)) = i
            SF_global_node_to_vertex(SF_global_vertices(i)%node(1)) = i
            !- now we link vertices to faces
            allocate(SF_global_vertices(i)%faces(0:SF_global_vertices(i)%valence-1))
            SF_global_vertices(i)%faces(0:SF_global_vertices(i)%valence-1) =    &
                SF_vert_to_face(0:SF_global_vertices(i)%valence-1,i)
        end do

    end subroutine SF_global_vertices_construct
    !------------------------------------------------------
    !------------------------------------------------------
    subroutine SF_global_edges_construct(SF_n_global_faces,Ipointer,  &
        SF_global_node_to_vertex, SF_global_faces,          &
        SF_n_global_edges,SF_global_edges)

        implicit none
        integer, intent(in)      :: SF_n_global_faces
        integer, dimension(0:), intent(in)   :: SF_global_node_to_vertex
        integer, dimension(0:,0:), intent(in)   :: Ipointer
        type(SF_face), dimension(0:SF_n_global_faces-1)  :: SF_global_faces
        integer, intent(out)   :: SF_n_global_edges
        type(SF_edge), dimension(:), allocatable, intent(out)  :: SF_global_edges

        integer, dimension(0:4*SF_n_global_faces-1)  :: edge_valence
        integer, dimension(0:7,0:4*SF_n_global_faces-1) :: SF_edge_to_face,  &
            SF_edge_orient_tmp
        integer, dimension(0:4*SF_n_global_faces-1,0:1) :: SF_vertex_to_edge
        integer  :: i,j,n,nf,ne,nv0,nv1,ok
        integer, dimension(0:1)  :: corner_edge, e_corner, e_vertex
        integer, dimension(0:3)  :: corner_fl, edge_fl

        write(*,*) "  --> Construction of global solid/fluid edges"
        SF_n_global_edges = 0

        do i = 0,SF_n_global_faces-1
            n = SF_global_faces(i)%elem(0)
            nf = SF_global_faces(i)%face(0)
            call caract_face(nf,corner_fl,edge_fl)
            do j = 0,3
                corner_edge(0) = corner_fl(j)
                if(j == 2) then
                    corner_edge(0) = corner_fl(3)
                    corner_edge(1) = corner_fl(2)
                else if (j == 3) then
                    corner_edge(0) = corner_fl(0)
                    corner_edge(1) = corner_fl(3)
                else
                    corner_edge(1) = corner_fl(j+1)
                endif
                !- global nodes..
                e_corner(0) = Ipointer(corner_edge(0),n)
                e_corner(1) = Ipointer(corner_edge(1),n)
                !- .. and SF vertices associated:
                e_vertex(0) = SF_global_node_to_vertex(e_corner(0))
                e_vertex(1) = SF_global_node_to_vertex(e_corner(1))
                ok = 0
                findSFedge: do ne = 0,SF_n_global_edges-1  ! edge already seen?
                    nv0 = SF_vertex_to_edge(ne,0)
                    nv1 = SF_vertex_to_edge(ne,1)
                    if((nv0 == e_vertex(0) .and. nv1 == e_vertex(1)) .or. &
                        (nv0 == e_vertex(1) .and. nv1 == e_vertex(0)))then
                        ! yes, SF edge seen before
                        ok = 1
                        SF_edge_to_face(edge_valence(ne),ne) = i
                        SF_global_faces(i)%edge(j) = ne
                        if(nv0 == e_vertex(0))then
                            SF_edge_orient_tmp(edge_valence(ne),ne) = 0
                        else
                            SF_edge_orient_tmp(edge_valence(ne),ne) = 1
                        end if
                        edge_valence(ne) = edge_valence(ne)+1
                        exit findSFedge
                    end if
                end do findSFedge
                if(ok == 0)then    ! SF edge never seen
                    SF_edge_to_face(0,SF_n_global_edges) = i
                    SF_vertex_to_edge(SF_n_global_edges,0:1) = e_vertex(0:1)
                    SF_global_faces(i)%edge(j) = SF_n_global_edges
                    SF_edge_orient_tmp(0,SF_n_global_edges) = 0
                    edge_valence(SF_n_global_edges) = 1
                    SF_n_global_edges = SF_n_global_edges+1
                end if
            end do
        end do
        write(*,*) "    --> Nb of global solid/fluid edges: ", SF_n_global_edges
        allocate(SF_global_edges(0:SF_n_global_edges-1))
        do i = 0,SF_n_global_edges-1
            SF_global_edges(i)%vertex(0:1) = SF_vertex_to_edge(i,0:1)
            SF_global_edges(i)%valence = edge_valence(i)
            allocate(SF_global_edges(i)%faces(0:edge_valence(i)-1),       &
                SF_global_edges(i)%orient(0:edge_valence(i)-1))
            SF_global_edges(i)%faces(0:edge_valence(i)-1) =        &
                SF_edge_to_face(0:edge_valence(i)-1,i)
            SF_global_edges(i)%orient(0:edge_valence(i)-1) =              &
                SF_edge_orient_tmp(0:edge_valence(i)-1,i)
        end do

    end subroutine SF_global_edges_construct
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine SF_vertices_proc_belong(n_proc,SF_n_global_vertices,part,    &
        elem_solid,initnode,SF_global_vertices)

        implicit none
        integer, intent(in)     :: n_proc,SF_n_global_vertices
        integer, dimension(0:n_proc-1), intent(in)  :: part
        logical, dimension(0:), intent(in)    :: elem_solid
        type(near_node), dimension(0:), intent(in)  ::  initnode
        type(SF_vertex), dimension(0:SF_n_global_vertices-1), intent(inout) ::  &
            SF_global_vertices

        integer   :: i,nf,ns
        type(near_entity), pointer  :: near_neighb => NULL()

        do i = 0, SF_n_global_vertices-1
            allocate(SF_global_vertices(i)%local_fluid(0:n_proc-1))
            SF_global_vertices(i)%local_fluid(0:) = .false.
            allocate(SF_global_vertices(i)%local_solid(0:n_proc-1))
            SF_global_vertices(i)%local_solid(0:) = .false.
            nf = SF_global_vertices(i)%node(0)
            ns = SF_global_vertices(i)%node(1)
            near_neighb => initnode(ns)%ptr
            do while(associated(near_neighb))
                if(elem_solid(near_neighb%elem))then
                    SF_global_vertices(i)%local_solid(part(near_neighb%elem)) = .true.
                else
                    SF_global_vertices(i)%local_fluid(part(near_neighb%elem)) = .true.
                end if
                near_neighb => near_neighb%pt
            end do
        end do


    end subroutine SF_vertices_proc_belong
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine SF_edges_proc_belong(n_proc,SF_n_global_edges,     &
        SF_global_vertices,SF_global_edges)

        implicit none
        integer, intent(in)     :: n_proc,SF_n_global_edges
        type(SF_vertex), dimension(0:), intent(in)   :: SF_global_vertices
        type(SF_edge), dimension(0:SF_n_global_edges-1), intent(inout) ::   &
            SF_global_edges

        integer   :: i,n,nv0,nv1

        do i = 0,SF_n_global_edges-1
            allocate(SF_global_edges(i)%local_fluid(0:n_proc-1))
            SF_global_edges(i)%local_fluid(0:) = .false.
            allocate(SF_global_edges(i)%local_solid(0:n_proc-1))
            SF_global_edges(i)%local_solid(0:) = .false.
            nv0 = SF_global_edges(i)%vertex(0)
            nv1 = SF_global_edges(i)%vertex(1)
            do n = 0,n_proc-1
                if(SF_global_vertices(nv0)%local_fluid(n) .and.   &
                    SF_global_vertices(nv1)%local_fluid(n))       &
                    SF_global_edges(i)%local_fluid(n) = .true.
                if(SF_global_vertices(nv0)%local_solid(n) .and.   &
                    SF_global_vertices(nv1)%local_solid(n))       &
                    SF_global_edges(i)%local_solid(n) = .true.
            end do
        end do


    end subroutine SF_edges_proc_belong
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine SF_edges_reference_elem(nproc,SF_n_global_edges,Elem_glob2loc,  &
        part,Ipointer,elem_solid,SF_global_vertices,                 &
        SF_edges_near_elem,SF_global_edges)

        implicit none
        integer, intent(in)   :: nproc,SF_n_global_edges
        integer, dimension(0:), intent(in)   :: Elem_glob2loc, part
        integer, dimension(0:,0:), intent(in)   :: Ipointer
        logical, dimension(0:), intent(in)     :: elem_solid
        type(SF_vertex), dimension(0:), intent(in)   :: SF_global_vertices
        type(near_elem), dimension(0:SF_n_global_edges-1), intent(in)   ::   &
            SF_edges_near_elem
        type(SF_edge), dimension(0:SF_n_global_edges-1), intent(inout)  ::   &
            SF_global_edges

        logical, dimension(:), allocatable   :: Lproc_solid,Lproc_fluid
        integer   :: i,j,k,nv0,nv1,n,np,nel,num,ok
        integer, dimension(0:1)  :: npf,nps,e_neighbor_corner
        type(near_entity), pointer   :: near_neighb => NULL()

        allocate(Lproc_solid(0:nproc-1),Lproc_fluid(0:nproc-1))
        do i = 0,SF_n_global_edges-1
            Lproc_solid(:) = .true. ; Lproc_fluid(:) = .true.
            allocate(SF_global_edges(i)%elem_ref_fluid(0:nproc-1,0:1))
            allocate(SF_global_edges(i)%elem_ref_solid(0:nproc-1,0:1))
            SF_global_edges(i)%elem_ref_fluid = -1
            SF_global_edges(i)%elem_ref_solid = -1

            nv0 = SF_global_edges(i)%vertex(0)
            nv1 = SF_global_edges(i)%vertex(1)
            npf(0) = SF_global_vertices(nv0)%node(0)
            nps(0) = SF_global_vertices(nv0)%node(1)
            npf(1) = SF_global_vertices(nv1)%node(0)
            nps(1) = SF_global_vertices(nv1)%node(1)
            near_neighb => SF_edges_near_elem(i)%ptr

            do while(associated(near_neighb))
                nel = near_neighb%elem
                n = Elem_glob2loc(nel)
                np = part(nel)
                ! fluid part
                if((.not. elem_solid(nel)) .and. Lproc_fluid(np)      &
                    .and. SF_global_edges(i)%local_fluid(np))then
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
                        SF_global_edges(i)%elem_ref_fluid(np,0) = n
                        SF_global_edges(i)%elem_ref_fluid(np,1) =         &
                            vertices2edge(e_neighbor_corner)
                        Lproc_fluid(np) = .false.
                    end if
                end if
                ! solid part
                if((elem_solid(nel)) .and. Lproc_solid(np)            &
                    .and. SF_global_edges(i)%local_solid(np))then
                    do j = 0,1
                        num = nps(j)
                        ok = 0
                        do_k2: do k = 0,7
                            if(num == Ipointer(k,nel))then
                                e_neighbor_corner(j) = k
                                ok = 1
                                exit do_k2
                            end if
                        end do do_k2
                    end do
                    if(j == 2 .and. ok == 1)then
                        SF_global_edges(i)%elem_ref_solid(np,0) = n
                        SF_global_edges(i)%elem_ref_solid(np,1) =     &
                            vertices2edge(e_neighbor_corner)
                        Lproc_solid(np) = .false.
                    end if
                end if
                near_neighb => near_neighb%pt
            end do
        end do

        deallocate(Lproc_fluid,Lproc_solid)


    end subroutine SF_edges_reference_elem
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine SF_local_faces_construct(proc,nproc,SF_n_faces,SF_n_global_faces,  &
        part,Elem_glob2loc,SF_global_face_present,SF_global_faces,        &
        SF_object_n_faces,SF_object_face,faces,SF_faces,SF_faces_shared,  &
        SF_nf_shared,SF_face_orient,SF_local_to_global_faces,MemorySF)

        implicit none
        integer, intent(in)   :: proc,nproc,SF_n_faces,SF_n_global_faces
        logical, dimension(0:), intent(in)  :: SF_global_face_present
        integer, dimension(0:), intent(in)  :: part,Elem_glob2loc
        type(SF_face), dimension(0:), intent(in)  :: SF_global_faces
        integer, dimension(0:), intent(in)  :: SF_object_n_faces
        type(Face_type), dimension(0:,0:), intent(in)   :: SF_object_face
        integer, dimension(0:,0:), intent(in)   :: faces
        integer, dimension(0:SF_n_faces-1,0:1), intent(out)  :: SF_faces
        integer, dimension(0:,0:), intent(out)   :: SF_faces_shared
        integer, dimension(0:nproc-1), intent(out)   :: SF_nf_shared
        integer, dimension(0:SF_n_faces-1), intent(out) :: SF_face_orient, &
            SF_local_to_global_faces
        type(process_obj), dimension(0:), intent(inout)  :: MemorySF
        integer :: i,j,k,nf,nel,nnf,num,n

        SF_nf_shared(0:) = 0
        SF_face_orient(0:) = -1
        j = 0    ! local SF faces counting
        do nf = 0,SF_n_global_faces-1
            if(.not. SF_global_face_present(nf)) cycle
            SF_local_to_global_faces(j) = nf
            !- fluid side
            nel = SF_global_faces(nf)%elem(0)   ! global number of the fluid element sharing the face
            nnf = SF_global_faces(nf)%face(0)   ! number of the face
            if(part(nel) == proc)then   ! the fluid face is on the proc
                n = Elem_glob2loc(nel)
                SF_faces(j,0) = faces(n,nnf)
                do k = 0,SF_object_n_faces(nel)-1
                    if(SF_object_face(nel,k)%ind_face == nnf) exit
                end do
                SF_face_orient(j) = SF_Object_face(nel,k)%orient
            else   ! the fluid face is not on this proc, but the solid side is!!
                SF_faces(j,0) = -1
                num = part(nel)
                if(num < proc)then   ! proc already seen: fluid face is already
                    !   identified as a shared SF face
                    SF_faces_shared(num,MemorySF(proc)%F(num,nf)) = j
                else     ! proc not seen yet
                    SF_faces_shared(num,SF_nf_shared(num)) = j
                    MemorySF(num)%F(proc,nf) = SF_nf_shared(num)
                end if
                SF_nf_shared(num) = SF_nf_shared(num)+1
            end if
            !- solid side
            nel = SF_global_faces(nf)%elem(1)   ! global number of the solid element sharing the face
            nnf = SF_global_faces(nf)%face(1)   ! number of the face
            if(part(nel) == proc)then   ! the solid face is on the proc
                n = Elem_glob2loc(nel)
                SF_faces(j,1) = faces(n,nnf)
                do k = 0,SF_object_n_faces(nel)-1
                    if(SF_object_face(nel,k)%ind_face == nnf) exit
                end do
                SF_face_orient(j) = SF_Object_face(nel,k)%orient
            else   ! the solid face is not on this proc, but the fluid side is!!
                SF_faces(j,1) = -1
                num = part(nel)
                if(num < proc)then   ! proc already seen: solid face is already
                    !    identified as a shared SF face
                    SF_faces_shared(num,MemorySF(proc)%F(num,nf)) = j
                else     ! proc not seen yet
                    SF_faces_shared(num,SF_nf_shared(num)) = j
                    MemorySF(num)%F(proc,nf) = SF_nf_shared(num)
                end if
                SF_nf_shared(num) = SF_nf_shared(num)+1
            end if
            j = j+1    ! one more local SF face
        end do


    end subroutine SF_local_faces_construct
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine SF_local_edges_construct(proc,nproc,SF_n_edges,SF_n_global_edges,    &
        part,glob2loc,SF_global_edges,SF_global_vertices,edges,              &
        Ipointer_local,Ipointer,which_elem_in_proc,SF_global_to_local_edges, &
        SF_edges,SF_edges_shared, SF_mapping_edges_shared,SF_ne_shared,      &
        SF_mapping_edges,MemorySF)

        implicit none
        integer, intent(in)  :: proc,nproc,SF_n_edges,SF_n_global_edges
        integer, dimension(0:), intent(in)  :: part,glob2loc
        type(SF_edge), dimension(0:), intent(in)  :: SF_global_edges
        type(SF_vertex), dimension(0:), intent(in)  :: SF_global_vertices
        integer, dimension(0:,0:), intent(in)   :: edges,Ipointer_local,   &
            Ipointer,which_elem_in_proc
        integer, dimension(0:), intent(out)  :: SF_global_to_local_edges
        integer, dimension(0:SF_n_edges-1,0:1), intent(out)  :: SF_edges
        integer, dimension(0:,0:), intent(out)   :: SF_edges_shared,   &
            SF_mapping_edges_shared
        integer, dimension(0:nproc-1), intent(out)   :: SF_ne_shared
        integer, dimension(0:SF_n_edges-1), intent(out) :: SF_mapping_edges
        type(process_obj), dimension(0:), intent(inout)  :: MemorySF
        integer :: i,j,k,nf,nel,nnf,num,n,ns,nes,ne,n0,n1,nns,nv0,nels
        logical   ::  orient_fluid,orient_solid,orient_fluid_loc,orient_solid_loc,  &
            orient_o_proc_fluid,orient_o_proc_solid


        SF_edges(0:,0:) = -1
        SF_mapping_edges(0:) = -1
        SF_ne_shared(0:) = 0
        SF_mapping_edges_shared(0:,0:) = -1

        j = 0   ! local SF edges counting
        do i = 0,SF_n_global_edges-1
            if((.not. SF_global_edges(i)%local_fluid(proc)) .and.    &
                      (.not. SF_global_edges(i)%local_solid(proc))) cycle
            orient_fluid_loc = .false. ; orient_solid_loc = .false.
           !- yes, SF edge:
            SF_global_to_local_edges(i) = j   
           !- fluid side
            if(SF_global_edges(i)%local_fluid(proc))then
                n = SF_global_edges(i)%elem_ref_fluid(proc,0)
                ne = SF_global_edges(i)%elem_ref_fluid(proc,1)
                SF_edges(j,0) = edges(n,ne)
                nv0 = SF_global_edges(i)%vertex(0)
                nnf = SF_global_vertices(nv0)%node(0)
                n0 = glob2loc(nnf)
                if(n0 == Ipointer_local(edge2vertex(ne),n)) orient_fluid_loc = .true.
            end if
           !- solid side
            if(SF_global_edges(i)%local_solid(proc))then
                ns = SF_global_edges(i)%elem_ref_solid(proc,0)
                nes = SF_global_edges(i)%elem_ref_solid(proc,1)
                SF_edges(j,1) = edges(ns,nes)
                nv0 = SF_global_edges(i)%vertex(0)
                nns = SF_global_vertices(nv0)%node(1)
                n1 = glob2loc(nns)
                if(n1 == Ipointer_local(edge2vertex(nes),ns)) orient_solid_loc = .true.
            end if  

           ! intraproc orientation:
            SF_mapping_edges(j) = merge(0,1,orient_fluid_loc .eqv. orient_solid_loc)

           !- now we look for SF edges shared with other procs
            do num = 0,nproc-1
                if(num == proc) cycle
                if((.not. SF_global_edges(i)%local_fluid(num)) .and.    &
                           (.not. SF_global_edges(i)%local_solid(num))) cycle
                if(((.not. SF_global_edges(i)%local_fluid(proc)) .and.    &
                   (.not. SF_global_edges(i)%local_fluid(num))).or.    &
                   ((.not. SF_global_edges(i)%local_solid(proc)) .and. &
                   (.not. SF_global_edges(i)%local_solid(num)))) cycle
              ! elimination of edges which do not exchange Solid/Fluid information (only S/S or F/F)
              ! now we have found another proc on which we have the same global SF edge
                orient_solid = .false. ; orient_fluid = .false.
                orient_o_proc_solid = .true. ; orient_o_proc_fluid = .true.
              ! interproc orientation
                if(SF_global_edges(i)%local_fluid(proc) .and. SF_global_edges(i)%local_solid(num))then
                    ns = SF_global_edges(i)%elem_ref_solid(num,0)
                    nels = which_elem_in_proc(num,ns)
                    nes = SF_global_edges(i)%elem_ref_solid(num,1)
                    nv0 = SF_global_edges(i)%vertex(0)
                    nns = SF_global_vertices(nv0)%node(1)
                    if(nns == Ipointer(edge2vertex(nes),nels)) orient_solid = .true.
                    orient_o_proc_fluid = merge(.true.,.false.,orient_fluid_loc .eqv. orient_solid)
                end if 
                if(SF_global_edges(i)%local_solid(proc) .and. SF_global_edges(i)%local_fluid(num))then
                    ns = SF_global_edges(i)%elem_ref_fluid(num,0)
                    nels = which_elem_in_proc(num,ns)
                    nes = SF_global_edges(i)%elem_ref_fluid(num,1)
                    nv0 = SF_global_edges(i)%vertex(0)
                    nns = SF_global_vertices(nv0)%node(0)
                    if(nns == Ipointer(edge2vertex(nes),nels)) orient_fluid = .true.
                    orient_o_proc_solid = merge(.true.,.false.,orient_fluid .eqv. orient_solid_loc)
                end if 
                if(orient_o_proc_solid .neqv. orient_o_proc_fluid) stop "Problem in SF edge orientation"
                if(num < proc)then   ! proc already seen
                    SF_edges_shared(num,MemorySF(proc)%E(num,i)) = j
                else   ! proc not seen yet
                    SF_edges_shared(num,SF_ne_shared(num)) = j
                    MemorySF(num)%E(proc,i) = SF_ne_shared(num)
                end if
                SF_mapping_edges_shared(num,SF_ne_shared(num)) = merge(0,1,orient_o_proc_fluid)
                SF_ne_shared(num) = SF_ne_shared(num)+1
            end do 
 !--  
            j = j+1   ! one more local SF edge
        end do


    end subroutine SF_local_edges_construct
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine SF_local_vertices_construct(proc,nproc,SF_n_vertices,              &
        SF_n_global_vertices,glob2loc,N_valid_Vertex,SF_global_vertices, &
        SF_global_to_local_vertices,SF_vertices,SF_vertices_shared,      &
        SF_nv_shared,MemorySF)

        implicit none
        integer, intent(in)  :: proc,nproc,SF_n_vertices,SF_n_global_vertices
        integer, dimension(0:), intent(in)  :: glob2loc,N_valid_Vertex
        type(SF_vertex), dimension(0:), intent(in)  :: SF_global_vertices
        integer, dimension(0:), intent(out)  :: SF_global_to_local_vertices
        integer, dimension(0:SF_n_vertices-1,0:1), intent(out) :: SF_vertices
        integer, dimension(0:,0:), intent(out)   :: SF_vertices_shared
        integer, dimension(0:nproc-1), intent(out)   :: SF_nv_shared
        type(process_obj), dimension(0:), intent(inout)  :: MemorySF
        integer :: i,j,k,nf,nel,nnf,num,n,ns,nes,ne,n0,n1,nns,nv0,nels

        SF_vertices(0:,0:) = -1
        SF_nv_shared(0:) = 0

        j = 0   ! local SF vertices counting
        do i = 0,SF_n_global_vertices-1
            if((.not. SF_global_vertices(i)%local_fluid(proc)) .and.    &
                (.not. SF_global_vertices(i)%local_solid(proc))) cycle
            SF_global_to_local_vertices(i) = j
            !- fluid side
            if(SF_global_vertices(i)%local_fluid(proc))then
                nf = SF_global_vertices(i)%node(0)
                nnf = glob2loc(nf)
                SF_vertices(j,0) = N_valid_Vertex(nnf)
            end if
            !- solid side
            if(SF_global_vertices(i)%local_solid(proc))then
                ns = SF_global_vertices(i)%node(1)
                nns = glob2loc(ns)
                SF_vertices(j,1) = N_valid_Vertex(nns)
            end if

            !- now we look for SF vertices shared with other procs
            do num = 0,nproc-1
                if(num == proc) cycle
                if((.not. SF_global_vertices(i)%local_fluid(num)) .and.   &
                    (.not. SF_global_vertices(i)%local_solid(num))) cycle
                if(((.not. SF_global_vertices(i)%local_fluid(proc)) .and.     &
                    (.not. SF_global_vertices(i)%local_fluid(num))).or.    &
                    ((.not. SF_global_vertices(i)%local_solid(proc)) .and. &
                    (.not. SF_global_vertices(i)%local_solid(num)))) cycle
                ! elimination of vertices which do not exchange Solid/Fluid information (only S/S or F/F)
                ! now we have found a proc on which we have the same global SF vertex
                if(num < proc)then   ! proc already seen
                    SF_vertices_shared(num,MemorySF(proc)%V(num,i)) = j
                else   ! proc not seen yet
                    SF_vertices_shared(num,SF_nv_shared(num)) = j
                    MemorySF(num)%V(proc,i) = SF_nv_shared(num)
                end if
                SF_nv_shared(num) = SF_nv_shared(num)+1
            end do

            j = j+1   ! one more local SF vertex
        end do


    end subroutine SF_local_vertices_construct
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine SF_faces_to_EV(proc,SF_n_faces,SF_local_to_global_faces,       &
        SF_global_to_local_edges,SF_global_to_local_vertices,           &
        SF_global_faces,SF_global_edges,SF_faces,mapping_edges,         &
        SF_Face_Near_Edges,SF_Face_Near_Edges_Orient,SF_Face_Near_Vertices)

        implicit none
        integer, intent(in)  :: proc,SF_n_faces
        integer, dimension(0:), intent(in)  :: SF_local_to_global_faces, &
            SF_global_to_local_edges,SF_global_to_local_vertices
        type(SF_face), dimension(0:), intent(in)  :: SF_global_faces
        type(SF_edge), dimension(0:), intent(in)  :: SF_global_edges
        integer, dimension(0:,0:), intent(in)  :: SF_faces,mapping_edges
        integer, dimension(0:,0:), intent(out)  :: SF_Face_Near_Edges,   &
            SF_Face_Near_Edges_Orient,SF_Face_Near_Vertices
        integer   ::  i,j,k,n,ne,nel,nes,nv,nvl

        SF_Face_Near_Edges(0:,0:) = -1
        SF_Face_Near_Edges_Orient(0:,0) = -1
        SF_Face_Near_Vertices(0:,0:) = -1
        do j = 0, SF_n_faces-1
            i = SF_local_to_global_faces(j)
            do k = 0,3
                ! edges
                ne = SF_global_faces(i)%edge(k)
                nel = SF_global_to_local_edges(ne)
                SF_Face_Near_Edges(j,k) = nel
                if(SF_faces(j,0) > -1)then
                    n = SF_global_edges(ne)%elem_ref_fluid(proc,0)
                    nes = SF_global_edges(ne)%elem_ref_fluid(proc,1)
                else
                    n = SF_global_edges(ne)%elem_ref_solid(proc,0)
                    nes = SF_global_edges(ne)%elem_ref_solid(proc,1)
                end if
                SF_Face_Near_Edges_Orient(j,k) = mapping_edges(n,nes)
                ! vertices
                nv = SF_global_faces(i)%vertex(k)
                nvl = SF_global_to_local_vertices(nv)
                SF_Face_Near_Vertices(j,k) = nvl
            end do

        end do

    end subroutine SF_faces_to_EV
    !--------------------------------------------------------------
    !--------------------------------------------------------------
end module solid_fluid
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
