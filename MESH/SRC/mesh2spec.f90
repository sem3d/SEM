module mesh2spec

    use mesh_properties
    use sets
    use partition_mesh
    use solid_fluid
    use local_mesh_properties
    use neumann

contains

    !---------------------
    subroutine gen_mesh(tabmat,nproc,pml_b,pml_t,pml_bottom,nmatref)

        implicit none

        !- general --
        integer   :: iunit,choice,i_ex, n,nn,nnn,i,j,k,l,idummy,icount,ok,num, &
            nel,nf,ns,ne,nv,np,i_count,nx,ny,nz,nnf,ind_f,nelf,nels,  &
            nns,kf,ks,ii,jj,ind_e,nv0,nv1,nvl,n0,n1,nes,nfile,nface
        character(len=13) :: meshfilename

        !- mesh --
        integer   :: n_nods, n_elem, n_blocks, n_points ! n_blocks corresponds
        ! to the number of different materials (cf Cubit)
        integer   :: corner(0:3), neighbor_corner(0:3), corner_edge(0:1),      &
            orient_corner(0:3), e_corner(0:1), e_vertex(0:1),         &
            e_neighbor_corner(0:1), e_orient_corner(0:1)
        integer   :: neighbor, neighbor_face, neighbor_edge
        integer, allocatable :: elmnts(:), material(:), Ipointer(:,:)
        real, allocatable  :: xco(:),yco(:),zco(:), Gcoord(:,:)
        integer  :: n_old_points, n_dim, vert
        logical  :: pml_bool, curve

        !- nearest elements
        type(near_entity), pointer                  :: near_neighb => NULL(),  &
            o_proc_near_neighb => NULL()
        type(near_node), dimension(:), allocatable  :: initnode,prevnode,currnode
        type(near_elem), dimension(:), allocatable  :: near_elem_set,     &
            SF_edges_near_elem, Neu_edges_near_elem

        !- materials --
        integer      :: n_SF_nodes
        integer, intent(in), optional  :: nmatref
        character, dimension(0:), intent(in)  :: tabmat
        integer, intent(in), optional   :: pml_b, pml_t, pml_bottom
        logical    :: any_fluid, all_fluid, solid_fluid
        integer, allocatable :: nodes_nature(:)
        logical, allocatable :: elem_contact(:), elem_solid(:)
        !- solid/fluid --
        integer, allocatable :: Nodes_on_SF(:,:), SF_object_n_faces(:),             &
            Elem_Edge_Ref(:), node_list(:), node_valence(:),    &
            SF_vert_to_face(:,:), SF_global_node_to_vertex(:),  &
            edge_valence(:), SF_edge_to_face(:,:),              &
            SF_edge_orient_tmp(:,:), SF_vertex_to_edge(:,:),    &
            SF_global_to_local_edges(:),                        &
            SF_global_to_local_vertices(:),                     &
            SF_local_to_global_faces(:)
        logical, allocatable :: SF_object(:,:), SF_global_face_present(:),          &
            SF_global_edge_present(:), ok_tmp(:)
        type(Face_type), allocatable :: SF_object_Face(:,:)
        type(SF_face), allocatable :: SF_global_faces(:)
        type(SF_edge), allocatable :: SF_global_edges(:)
        type(SF_vertex), allocatable :: SF_global_vertices(:)
        integer   :: corner_fl(0:3), corner_s(0:3), edge_fl(0:3),              &
            edge_s(0:3), npf(0:1),nps(0:1)
        integer   :: SF_n_global_faces, SF_n_global_vertices, SF_n_global_edges
        !- interfaces: Neumann BC
        integer   :: Neu_n_global_faces, Neu_n_global_vertices, Neu_n_global_edges
        logical, dimension(:), allocatable  :: Neu_object
        integer, dimension(:), allocatable :: Neu_object_n_faces,              &
            Neu_global_node_to_vertex,                    &
            Neu_global_to_local_edges(:),                   &
            Neu_global_to_local_vertices(:),                &
            Neu_local_to_global_faces(:)
        integer, dimension(:,:), allocatable :: Faces_Neumann
        logical   :: Neumann_present = .false.
        type(Neu_face_type), dimension(:,:), allocatable   :: Neu_object_face
        type(Neu_face), dimension(:), allocatable  :: Neu_global_faces
        type(Neu_vertex), dimension(:), allocatable  :: Neu_global_vertices
        type(Neu_edge), dimension(:), allocatable  :: Neu_global_edges
        logical, dimension(:), allocatable  :: Neu_global_face_present,   &
            Neu_global_edge_present

        !- interfaces: Plane Wave
        integer   :: n_PW
        integer, dimension(:,:), allocatable :: Faces_PW
        !- partitioning - see Metis user doc for parameter value
        integer, intent(in)  :: nproc
        integer, allocatable :: dxadj(:), dxadjncy(:), part(:)
        !- local meshes
        type(souvenir), dimension(:), allocatable :: memory
        type(process_obj), dimension(:), allocatable :: MemorySF
        type(proc_obj), dimension(:), allocatable :: MemoryNeu
        logical, allocatable :: L_Proc(:), Lproc_fluid(:), Lproc_solid(:),        &
            which_points_inside(:), L_vertex(:)
        integer  :: nparts,proc, n_vertices, n_faces, n_edges, n_points_local,    &
            NmaxVertex, local_n
        integer, dimension(:), allocatable :: counter, elem_glob2loc, glob2loc,   &
            N_valid_Vertex, vertex_to_glob,        &
            nelem_in_proc, which_nodes, nf_shared, &
            ne_shared, nv_shared, node_loc2glob
        integer, dimension(:,:), allocatable :: elmnts_local, which_elem_in_proc,  &
            Ipointer_local, vertices_shared, &
            vertices
        integer, dimension(:,:), allocatable :: faces, faces_shared, mapping_faces, &
            mapping_faces_shared
        integer, dimension(:,:), allocatable :: edges, edges_shared, mapping_edges, &
            mapping_edges_shared
        ! local Neumann objects
        logical  :: Neu_present_local
        integer   :: Neu_n_faces,Neu_n_edges,Neu_n_vertices
        integer, dimension(:), allocatable  :: Neu_faces, Neu_edges,         &
            Neu_ne_shared, Neu_vertices, Neu_nv_shared
        integer, allocatable, dimension(:,:) :: Neu_face_Near_edges,               &
            Neu_Face_Near_Edges_Orient,Neu_Face_Near_Vertices,  &
            Neu_edges_shared, Neu_vertices_shared,              &
            Neu_mapping_edges_shared

        !- local solid-fluid objects
        logical   :: SF_present_local, orient_fluid, orient_solid
        integer   :: n_solid, ne_solid
        integer   :: SF_n_faces, SF_n_edges, SF_n_vertices
        logical, allocatable   :: SF_log_edge(:), SF_log_vertex(:)
        integer, allocatable, dimension(:) :: SF_face_orient, SF_edge_orient,     &
            SF_nf_shared, SF_ne_shared,        &
            SF_nv_shared, SF_mapping_edges
        integer, allocatable, dimension(:,:) :: SF_faces, SF_edges, SF_vertices,  &
            SF_face_Near_edges, SF_Face_Near_Edges_Orient, &
            SF_Face_Near_Vertices, SF_faces_shared,        &
            SF_edges_shared, SF_vertices_shared,           &
            SF_mapping_faces_shared, SF_mapping_edges_shared


        !--------------------------------------------------------------
        !- INITIALIZATIONS OF DIFFERENT MESH AND MATERIAL PROPERTIES
        call mesh_init_3D(n_nods,n_points,n_elem,n_blocks,                  &
            xco,yco,zco,Ipointer,Material,tabmat,             &
            Neu_n_global_faces,n_PW,Faces_Neumann,Faces_PW,   &
            pml_b,pml_t,pml_bottom,nmatref)
        if(Neu_n_global_faces > 0) Neumann_present = .true.

        !------------------------------------------------------
        !- LOOKING FOR NEAREST NEIGHBOURS FOR EACH ELEMENT
        write(*,*)
        write(*,*) "*********************************************"
        write(*,*) "  Finding nearest elements.."
        write(*,*) "*********************************************"

        !- looking for elements containing each node: loop on nodes,
        !    and addition of each new element containing the node (linked list)
        allocate(initnode(0:n_points-1),currnode(0:n_points-1),prevnode(0:n_points-1))
        do n = 0,n_elem-1
            do i = 0,7
                j = Ipointer(i,n)
                if(.not. associated(initnode(j)%ptr))then
                    allocate(initnode(j)%ptr) ; currnode(j)%ptr => initnode(j)%ptr
                else
                    allocate(currnode(j)%ptr) ; prevnode(j)%ptr%pt => currnode(j)%ptr
                end if
                currnode(j)%ptr%elem = n ; prevnode(j)%ptr => currnode(j)%ptr
            end do
        end do
        !- union of sets
        allocate(near_elem_set(0:n_elem-1))
        call find_near_elem(Ipointer,initnode,near_elem_set)
        !- sorting sets
        do n = 0,n_elem-1
            call entity_sort(near_elem_set(n)%ptr)
        end do


        !--------------------------------------------------
        !-   FLUID PART: INITIALIZATIONS FOR GLOBAL PROPERTIES
        call solid_fluid_init(n_nods,n_points,n_blocks,n_elem,tabmat,Ipointer,   &
            Material,any_fluid,all_fluid,solid_fluid,n_SF_nodes,    &
            nodes_nature,elem_solid,elem_contact)


        !--------------------------------------------
        !-  New number of control points (nodes), if fluid present
        n_old_points = n_points ; n_points = n_points+n_SF_nodes
        allocate(Gcoord(0:n_points-1,0:2))
        do n = 0,n_old_points-1
            Gcoord(n,0) = xco(n)
            Gcoord(n,1) = yco(n)
            Gcoord(n,2) = zco(n)
        end do
        deallocate(xco,yco,zco)


        !-----------------------------------------------------------------------
        !-  PARTITIONING:
        !--    here we use the METIS library (perhaps try Scotch later), v.4.0.3
        !--    Eventually: upgrade to Metis 5.0, but the syntax has changed --

        call part_mesh_3D(n_elem,n_points,Ipointer,nproc,dxadj,dxadjncy,part)
        write(*,*) "  --> Partition done."


        !-----------------------------------------------------
        !-   NEUMANN INTERFACE : GLOBAL PROPERTIES
        if(Neumann_present)then
            write(*,*)
            write(*,*) "*************************************************"
            write(*,*) "  --> Neumann case : global properties"
            write(*,*) "*************************************************"
            !- construction of Neumann Objects: pieces of info
            !      relative to Neumann interfaces

            call Neu_object_construct(n_elem,Ipointer,Neu_n_global_faces,   &
                Faces_Neumann,initnode,Neu_object,            &
                Neu_object_Face,Neu_object_n_faces)

            !do n = 0,n_elem-1
            !   if(Neu_object(n))then
            !    print*,n,(Neu_object_face(n,j)%ind_face,j=0,Neu_object_n_faces(n)-1)
            !   end if
            !end do
            !read*

            !- now we construct GLOBAL Neumann faces; it shall be
            !-   of great help when constructing local Neumann faces
            allocate(Neu_global_faces(0:Neu_n_global_faces-1))
            allocate(Neu_global_face_present(0:Neu_n_global_faces-1))
            call Neu_global_faces_construct(n_elem,Neu_n_global_faces,      &
                Neu_object,Neu_object_face,Neu_object_n_faces,  &
                Neu_global_faces)

            !do n = 0,Neu_n_global_faces-1
            !    print*,"NEU",n,Neu_global_faces(n)%elem,Neu_global_faces(n)%face,Neu_global_faces(n)%num
            !end do
            !read*

            !- now we construct GLOBAL Neumann vertices
            allocate(Neu_global_node_to_vertex(0:n_points-1))
            call Neu_global_vertices_construct(Neu_n_global_faces,         &
                Ipointer,Neu_global_faces,Neu_global_node_to_vertex,    &
                Neu_n_global_vertices,Neu_global_vertices)
            !do n = 0,Neu_n_global_vertices-1
            !    print*,"NEU",n,Neu_global_vertices(n)%node
            !    print*,"NEU2",Neu_global_vertices(n)%faces(0:)
            !end do
            !read*
            !do n = 0,n_points-1
            !    if(Neu_global_node_to_vertex(n) < 0) cycle
            !    print*,"NEU3",n,Neu_global_node_to_vertex(n)
            !end do
            !read*

            !- now we construct GLOBAL Neumann edges
            call Neu_global_edges_construct(Neu_n_global_faces,Ipointer,  &
                Neu_global_node_to_vertex, Neu_global_faces,         &
                Neu_n_global_edges,Neu_global_edges)

            allocate(Neu_global_edge_present(0:Neu_n_global_edges-1))
            !do n = 0,Neu_n_global_edges-1
            !    print*,"NEU",n,Neu_global_edges(n)%vertex
            !    print*,"NEU2",Neu_global_edges(n)%faces(0:)
            !end do
            !read*


            !- search for elements sharing the same Neumann edge:
            !     intersection of sets
            allocate(Neu_edges_near_elem(0:Neu_n_global_edges-1))
            do i = 0, Neu_n_global_edges-1
                nv0 = Neu_global_edges(i)%vertex(0)
                nv1 = Neu_global_edges(i)%vertex(1)
                n0 = Neu_global_vertices(nv0)%node
                n1 = Neu_global_vertices(nv1)%node
                call entity_intersect(initnode(n0)%ptr,initnode(n1)%ptr,Neu_edges_near_elem(i)%ptr)
                call entity_sort(Neu_edges_near_elem(i)%ptr)
                ! call list_entity(Neu_edges_near_elem(i)%ptr)
            end do

            !- which procs do Neumann vertices belong to?
            call Neu_vertices_proc_belong(nproc,Neu_n_global_vertices,part,    &
                initnode,Neu_global_vertices)
            !do n = 0,Neu_n_global_vertices-1
            !   print*,"GANCEL",n,Neu_global_vertices(n)%local
            !end do
            !read*

            !- which procs do Neumann edges belong to?
            call Neu_edges_proc_belong(nproc,Neu_n_global_edges,     &
                Neu_global_vertices,Neu_global_edges)

            !do n = 0,Neu_n_global_edges-1
            !   print*,"GANCEL2",n,Neu_global_edges(n)%local
            !end do
            !read*

        end if     !- temporary end of Neumann case


        !--------------------------------------------------------------
        !-   BACK TO THE SOLID/FLUID GLOBAL PROPERTIES (if any SF interfaces)
        if(solid_fluid)then
            write(*,*)
            write(*,*) "*************************************************"
            write(*,*) "  --> Solid-Fluid case : global properties"
            write(*,*) "*************************************************"
            write(*,*) "  --> Number of control nodes if fluid present: ", n_points
            !- doubling of SF nodes
            write(*,*) "  --> Doubling of SF nodes."
            allocate(Nodes_on_SF(0:n_SF_nodes-1,0:1))
            call SF_nodes_doubling(Nodes_on_SF,nodes_nature,n_old_points)
            !- coordinates of the "new" fluid nodes
            do i = 0,n_SF_nodes-1
                j = Nodes_on_SF(i,0)
                Gcoord(n_old_points+i,0) = Gcoord(j,0)
                Gcoord(n_old_points+i,1) = Gcoord(j,1)
                Gcoord(n_old_points+i,2) = Gcoord(j,2)
            end do
            !- construction of Solid-Fluid Objects: pieces of info
            !      relative to SF interfaces
            call SF_object_construct(n_elem,Ipointer,elem_contact, elem_solid,    &
                nodes_nature,dxadj,dxadjncy,SF_object,     &
                SF_object_Face,SF_object_n_faces,          &
                SF_n_global_faces)

            !- now we construct GLOBAL solid/fluid faces; it shall be
            !-   of great help when constructing local SF faces
            allocate(SF_global_faces(0:SF_n_global_faces-1))
            allocate(SF_global_face_present(0:SF_n_global_faces-1))
            call SF_global_faces_construct(n_elem, SF_n_global_faces,        &
                SF_object,SF_object_face,SF_object_n_faces,  &
                SF_global_faces)


            !- now we construct GLOBAL solid/fluid vertices
            allocate(SF_global_node_to_vertex(0:n_points-1))
            call SF_global_vertices_construct(n_SF_nodes,SF_n_global_faces,       &
                Ipointer,nodes_on_SF,SF_global_faces,SF_global_node_to_vertex, &
                SF_n_global_vertices,SF_global_vertices)

            !- now we construct GLOBAL solid/fluid edges
            call SF_global_edges_construct(SF_n_global_faces,Ipointer,    &
                SF_global_node_to_vertex, SF_global_faces,          &
                SF_n_global_edges,SF_global_edges)
            allocate(SF_global_edge_present(0:SF_n_global_edges-1))


            !- search for elements sharing the same SF edge:
            !     intersection of sets
            allocate(SF_edges_near_elem(0:SF_n_global_edges-1))
            do i = 0, SF_n_global_edges-1
                nv0 = SF_global_edges(i)%vertex(0)
                nv1 = SF_global_edges(i)%vertex(1)
                n0 = SF_global_vertices(nv0)%node(1)
                n1 = SF_global_vertices(nv1)%node(1)
                call entity_intersect(initnode(n0)%ptr,initnode(n1)%ptr,SF_edges_near_elem(i)%ptr)
                call entity_sort(SF_edges_near_elem(i)%ptr)
                !   call list_entity(SF_edges_near_elem(i)%ptr)
            end do

            !- if Solid-Fluid: doubled nodes -> change of indices
            call indices_modif(n_nods,n_SF_nodes,SF_n_global_faces,Nodes_on_SF,    &
                SF_global_faces,Ipointer)


            !- which procs do SF vertices belong to?
            call SF_vertices_proc_belong(nproc,SF_n_global_vertices,part,    &
                elem_solid,initnode,SF_global_vertices)

            !- which procs do SF edges belong to?
            call SF_edges_proc_belong(nproc,SF_n_global_edges,     &
                SF_global_vertices,SF_global_edges)


        end if     !- temporary end of solid/fluid case



        !-------------------------------------------------------
        !-------------------------------------------------------
        !-------------------------------------------------------
        !-  CREATION OF SEM-style files
        !-------------------------------------------------------
        !-------------------------------------------------------
        !-------------------------------------------------------
        !-------------------------------------------------------
        write(*,*)
        write(*,*) "****************************************"
        write(*,*) " CREATION OF SEM-style DATA STUCTURES"
        write(*,*) "****************************************"

        !-------------------------------------------------
        !-  CREATION OF SEM-style files : GENERALITIES

        !- defining the number of elements per processor
        write(*,*) "  --> Number of elements per proc."
        allocate(nelem_in_proc(0:nproc-1))
        nelem_in_proc = 0
        do nel = 0,n_elem-1
            nelem_in_proc(part(nel)) = nelem_in_proc(part(nel))+1
        end do

        !- defining for each processor which elements are inside
        ! "Counter" refers to the local numberings and "nel" to the global one
        write(*,*) "  --> Defining global/local correspondences for elements."
        allocate(counter(0:nproc-1))
        counter = 0
        allocate(which_elem_in_proc(0:nproc-1,0:maxval(nelem_in_proc)-1))
        allocate(Elem_glob2loc(0:n_elem-1))  ! gives the local number from the global one
        do nel = 0,n_elem-1
            num = part(nel)
            which_elem_in_proc(num,counter(num)) = nel
            Elem_glob2loc(nel) = counter(num)
            counter(num) = counter(num) + 1
        enddo
        deallocate(counter)
        !- allocating "memory":
        !   -> establishes correspondence between objects shared by different processors
        allocate(memory(0:n_elem-1))
        do nel = 0,n_elem-1
            if(part(nel) /= nproc-1)then
                allocate(memory(nel)%rank(part(nel)+1:nproc-1))
                do proc = part(nel)+1,nproc-1
                    allocate(memory(nel)%rank(proc)%E(0:25))
                enddo
            endif
        enddo

        if(Neumann_present)then    !-- Neumann local part
            !- looking for reference elements for Neumann edges
            call Neu_edges_reference_elem(nproc,Neu_n_global_edges,       &
                Elem_glob2loc,part,Ipointer,Neu_global_vertices,   &
                Neu_edges_near_elem,Neu_global_edges)

            !do n = 0, Neu_n_global_edges-1
            !    print*,"REF. EDGE: ",n
            !    print*,Neu_global_edges(n)%elem_ref(0:,0)
            !    print*,Neu_global_edges(n)%elem_ref(0:,1)
            !    print*
            !end do
            !read*


            !- local <-> global correspondence
            allocate(Neu_global_to_local_edges(0:Neu_n_global_edges-1))
            allocate(Neu_global_to_local_vertices(0:Neu_n_global_vertices-1))

            !- Neumann memory
            allocate(MemoryNeu(0:nproc-1))
            do i = 0,nproc-1
                allocate(MemoryNeu(i)%E(0:nproc-1,0:Neu_n_global_edges-1))
                allocate(MemoryNeu(i)%V(0:nproc-1,0:Neu_n_global_vertices-1))
            end do
        end if   !-  temporary end of local Neumann properties


        if(solid_fluid)then    !--  Solid/fluid local part
            !- looking for reference elements for SF edges
            call SF_edges_reference_elem(nproc,SF_n_global_edges,Elem_glob2loc,    &
                part,Ipointer,elem_solid,SF_global_vertices,                 &
                SF_edges_near_elem,SF_global_edges)

            !- local <-> global correspondence
            allocate(SF_global_to_local_edges(0:SF_n_global_edges-1))
            allocate(SF_global_to_local_vertices(0:SF_n_global_vertices-1))

            !- SF memory
            allocate(MemorySF(0:nproc-1))
            do i = 0,nproc-1
                allocate(MemorySF(i)%F(0:nproc-1,0:SF_n_global_faces-1))
                allocate(MemorySF(i)%E(0:nproc-1,0:SF_n_global_edges-1))
                allocate(MemorySF(i)%V(0:nproc-1,0:SF_n_global_vertices-1))
            end do
        end if   !-  temporary end of local solid/fluid properties


        !- global nodes
        allocate(glob2loc(0:n_points-1))
        allocate(which_points_inside(0:n_points-1))

        !- we consider each processor, and create all structures for the meshfiles
        call system("rm -f mesh4spec.???")
        meshfilename(1:10) = "mesh4spec."


        !-------------------------------------------------
        !-  CREATION OF SEM-style files : LOOP ON PROCS
        write(*,*) "  --> Loop on procs!"
        do proc = 0,nproc-1
            write(*,*) "*---------------------*"
            write(*,*) "   --> Proc. # --  ", proc

            SF_present_local = .false.   !- solid_fluid properties
            Neu_present_local = .false.   !- Neumann properties

            !- which nodes belong to the processor?
            !- these nodes will be sorted according to the global numbering
            !- Creation of local nodes. Correspondences local <-> global nodes.
            allocate(Ipointer_local(0:n_nods-1,0:nelem_in_proc(proc)-1))
            call local_nodes_definition(proc,n_nods,n_points,nelem_in_proc,    &
                which_elem_in_proc,Ipointer,n_points_local,          &
                which_points_inside,glob2loc,Ipointer_local,node_loc2glob)
            !- vertices construction
            allocate(vertices(0:nelem_in_proc(proc)-1,0:7))
            call local_vertex_construct(proc,n_points_local,nelem_in_proc,        &
                Ipointer_local,n_vertices,vertices,vertex_to_glob,N_valid_Vertex)

            !- Neumann objects locally?
            if(Neumann_present)then
                Neu_global_face_present(:) = .false.
                Neu_global_edge_present(:) = .false.
                do i = 0,Neu_n_global_vertices-1
                    if(Neu_global_vertices(i)%local(proc))then
                        Neu_present_local = .true.
                        exit
                    end if
                end do
                if(Neu_present_local) write(*,*) "    --> Neumann objects present in proc. ",proc
                Neu_global_to_local_edges(0:) = -1
                Neu_global_to_local_vertices(0:) = -1
            end if

            !- SF objects locally? => only a SF node is enough, that means: an element in contact
            if(solid_fluid)then
                SF_global_face_present(:) = .false.
                SF_global_edge_present(:) = .false.
                do i = 0,nelem_in_proc(proc)-1
                    nel = which_elem_in_proc(proc,i)
                    if(elem_contact(nel))then
                        SF_present_local = .true.
                        exit
                    end if
                end do
                if(SF_present_local) write(*,*) "    --> Solid/fluid objects present in proc. ",proc
                SF_global_to_local_edges(0:) = -1
                SF_global_to_local_vertices(0:) = -1
            end if


            !- Now let's turn to faces: 6 associated  to each element
            !  Are also determined: the faces shared between different procs
            allocate(faces(0:nelem_in_proc(proc)-1,0:5))
            allocate(mapping_faces(0:nelem_in_proc(proc)-1,0:5))
            allocate(faces_shared(0:nproc-1,0:6*nelem_in_proc(proc)-1))
            allocate(mapping_faces_shared(0:nproc-1,0:6*nelem_in_proc(proc)-1))
            allocate(nf_shared(0:nproc-1))
            call local_faces_construct(proc,part,nelem_in_proc,which_elem_in_proc,   &
                Ipointer,dxadj,dxadjncy,n_faces,faces,mapping_faces,   &
                faces_shared,mapping_faces_shared,nf_shared,memory)


            !- associating to each element: 12 edges
            !- defining the edges shared with other processors
            allocate(edges(0:nelem_in_proc(proc)-1,0:11))
            allocate(mapping_edges(0:nelem_in_proc(proc)-1,0:11))
            allocate(Elem_Edge_Ref(0:12*nelem_in_proc(proc)-1))
            allocate(edges_shared(0:nproc-1,0:12*nelem_in_proc(proc)-1))
            allocate(mapping_edges_shared(0:nproc-1,0:12*nelem_in_proc(proc)-1))
            allocate(ne_shared(0:nproc-1))
            call local_edges_construct(nproc,proc,part,nelem_in_proc,Elem_glob2loc,   &
                which_elem_in_proc,Ipointer,near_elem_set,n_edges,           &
                edges,mapping_edges,edges_shared,mapping_edges_shared,       &
                ne_shared,Elem_edge_ref,memory)

            !- defining vertices shared with other processors
            allocate(vertices_shared(0:nproc-1,0:7*nelem_in_proc(proc)-1))
            allocate(nv_shared(0:nproc-1))
            call local_vertices_comm(n_vertices,proc,nproc,nelem_in_proc,part,   &
                vertex_to_glob,node_loc2glob,which_elem_in_proc,Ipointer,  &
                vertices,near_elem_set,vertices_shared,nv_shared,memory)


            !- now we include eventual Neumann objects
            if(Neu_present_local)then
                write(*,*) "    --> Neumann local meshing properties."
                Neu_n_faces = 0 ; Neu_n_edges = 0 ; Neu_n_vertices = 0
                !- local Neumann faces ---------------------
                ! - how many?
                do nf = 0, Neu_n_global_faces-1
                    if(part(Neu_global_faces(nf)%elem) == proc)then
                        Neu_global_face_present(nf) = .true.
                        Neu_n_faces = Neu_n_faces+1
                    end if
                end do

                allocate(Neu_faces(0:Neu_n_faces-1))
                allocate(Neu_local_to_global_faces(0:Neu_n_faces-1))

                call Neu_local_faces_construct(proc,nproc,Neu_n_faces,      &
                    Neu_n_global_faces,part,Elem_glob2loc,               &
                    Neu_global_face_present,Neu_global_faces,            &
                    faces,Neu_faces,Neu_local_to_global_faces)

                !       do nf = 0,Neu_n_faces-1
                !           print*,"GANCELLITA",proc,nf,Neu_faces(nf)
                !       end do
                !       read*

                !- local Neumann edges  ---------------------
                ! how many?
                do i = 0, Neu_n_global_edges-1
                    if((.not. Neu_global_edges(i)%local(proc))) cycle
                    Neu_n_edges = Neu_n_edges+1
                end do

                allocate(Neu_edges(0:Neu_n_edges-1))
                allocate(Neu_ne_shared(0:nproc-1))
                allocate(Neu_edges_shared(0:nproc-1,0:Neu_n_edges-1))
                allocate(Neu_mapping_edges_shared(0:nproc-1,0:Neu_n_edges-1))

                call Neu_local_edges_construct(proc,nproc,Neu_n_edges,Neu_n_global_edges, &
                    part,glob2loc,Neu_global_edges,Neu_global_vertices,edges,              &
                    Ipointer_local,Ipointer,which_elem_in_proc,Neu_global_to_local_edges,  &
                    Neu_edges,Neu_edges_shared, Neu_mapping_edges_shared,Neu_ne_shared,    &
                    MemoryNeu)

                !- local Neumann vertices
                ! how many?
                do i = 0, Neu_n_global_vertices-1
                    if((.not. Neu_global_vertices(i)%local(proc))) cycle
                    Neu_n_vertices = Neu_n_vertices+1
                end do

                allocate(Neu_vertices(0:Neu_n_vertices-1))
                allocate(Neu_nv_shared(0:nproc-1))
                allocate(Neu_vertices_shared(0:nproc-1,0:Neu_n_vertices-1))

                call Neu_local_vertices_construct(proc,nproc,Neu_n_vertices,     &
                    Neu_n_global_vertices,glob2loc,N_valid_Vertex,              &
                    Neu_global_vertices,Neu_global_to_local_vertices,           &
                    Neu_vertices,Neu_vertices_shared,Neu_nv_shared,MemoryNeu)

                !- finally: Neumann edges and vertices associated to Neumann faces
                !   -> important for the calculation of properties of normals to faces.
                allocate(Neu_Face_Near_Edges(0:Neu_n_faces-1,0:3))
                allocate(Neu_Face_Near_Edges_Orient(0:Neu_n_faces-1,0:3))
                allocate(Neu_Face_Near_Vertices(0:Neu_n_faces-1,0:3))

                call Neu_faces_to_EV(proc,Neu_n_faces,Neu_local_to_global_faces,  &
                    Neu_global_to_local_edges,Neu_global_to_local_vertices,         &
                    Neu_global_faces,Neu_global_edges,Neu_faces,mapping_edges,      &
                    Neu_Face_Near_Edges,Neu_Face_Near_Edges_Orient,Neu_Face_Near_Vertices)


            end if   !-   end of local Neumann properties



            !- now we include eventual SF objects
            if(SF_present_local)then
                write(*,*) "    --> Solid/fluid local meshing properties."
                SF_n_faces = 0 ; SF_n_edges = 0 ; SF_n_vertices = 0
                !- local SF faces ---------------------
                ! - how many?
                do nf = 0, SF_n_global_faces-1
                    if(part(SF_global_faces(nf)%elem(0)) == proc .or.   &
                        part(SF_global_faces(nf)%elem(1)) == proc)then
                        SF_global_face_present(nf) = .true.
                        SF_n_faces = SF_n_faces+1
                    end if
                end do

                allocate(SF_faces(0:SF_n_faces-1,0:1))
                allocate(SF_nf_shared(0:nproc-1))
                allocate(SF_faces_shared(0:nproc-1,0:SF_n_faces-1))
                allocate(SF_face_orient(0:SF_n_faces-1))
                allocate(SF_local_to_global_faces(0:SF_n_faces-1))

                call SF_local_faces_construct(proc,nproc,SF_n_faces,SF_n_global_faces,  &
                    part,Elem_glob2loc,SF_global_face_present,SF_global_faces,        &
                    SF_object_n_faces,SF_object_face,faces,SF_faces,SF_faces_shared,  &
                    SF_nf_shared,SF_face_orient,SF_local_to_global_faces,MemorySF)

                !- local SF edges  ---------------------
                ! how many?
                do i = 0, SF_n_global_edges-1
                    if((.not. SF_global_edges(i)%local_fluid(proc)) .and.   &
                        (.not. SF_global_edges(i)%local_solid(proc))) cycle
                    SF_n_edges = SF_n_edges+1
                end do

                allocate(SF_edges(0:SF_n_edges-1,0:1))
                allocate(SF_mapping_edges(0:SF_n_edges-1))
                allocate(SF_ne_shared(0:nproc-1))
                allocate(SF_edges_shared(0:nproc-1,0:SF_n_edges-1))
                allocate(SF_mapping_edges_shared(0:nproc-1,0:SF_n_edges-1))

                call SF_local_edges_construct(proc,nproc,SF_n_edges,SF_n_global_edges, &
                    part,glob2loc,SF_global_edges,SF_global_vertices,edges,             &
                    Ipointer_local,Ipointer,which_elem_in_proc,SF_global_to_local_edges,&
                    SF_edges,SF_edges_shared, SF_mapping_edges_shared,SF_ne_shared,     &
                    SF_mapping_edges,MemorySF)

                !- local SF vertices
                ! how many?
                do i = 0, SF_n_global_vertices-1
                    if((.not. SF_global_vertices(i)%local_fluid(proc)) .and.    &
                        (.not. SF_global_vertices(i)%local_solid(proc))) cycle
                    SF_n_vertices = SF_n_vertices+1
                end do

                allocate(SF_vertices(0:SF_n_vertices-1,0:1))
                allocate(SF_nv_shared(0:nproc-1))
                allocate(SF_vertices_shared(0:nproc-1,0:SF_n_vertices-1))

                call SF_local_vertices_construct(proc,nproc,SF_n_vertices,            &
                    SF_n_global_vertices,glob2loc,N_valid_Vertex,SF_global_vertices, &
                    SF_global_to_local_vertices,SF_vertices,SF_vertices_shared,      &
                    SF_nv_shared,MemorySF)


                !- finally: SF edges and vertices associated to SF faces -> important for
                !       the calculation of properties of normals to faces.
                allocate(SF_Face_Near_Edges(0:SF_n_faces-1,0:3))
                allocate(SF_Face_Near_Edges_Orient(0:SF_n_faces-1,0:3))
                allocate(SF_Face_Near_Vertices(0:SF_n_faces-1,0:3))

                call SF_faces_to_EV(proc,SF_n_faces,SF_local_to_global_faces,    &
                    SF_global_to_local_edges,SF_global_to_local_vertices,          &
                    SF_global_faces,SF_global_edges,SF_faces,mapping_edges,        &
                    SF_Face_Near_Edges,SF_Face_Near_Edges_Orient,SF_Face_Near_Vertices)


            end if   !-   end of local SF properties


            !- writing the meshfile
            write(meshfilename(11:13), '(i3.3)') proc
            open(11, file = trim(meshfilename))
            n_dim = 3 ; curve = .false.
            write(11,"(1i6,a)") n_dim
            write(11,*) solid_fluid
            write(11,*) all_fluid
            write(11,*) Neumann_present
            write(11,*) "Local nodes"
            write(11,"(1i6,a)") n_points_local
            write(11,*) curve
            do n = 0,n_points_local-1
                i_count = node_loc2glob(n)
                write (11,*) Gcoord(i_count,0:2)
            enddo
            write(11,*) "Nb elements"
            write(11,"(1i6,a)") nelem_in_proc(proc)
            write(11,*) "Materials"
            write(11,"(1i6,a)") n_blocks
            do n = 0,nelem_in_proc(proc)-1
                nel = which_elem_in_proc(proc,n)
                write(11,"(i6,1x,l7)") Material(nel),elem_solid(nel)
            enddo
            write(11,*) "Global nodes for elements"
            write(11,"(1i6,a)") n_nods
            do n = 0,nelem_in_proc(proc)-1
                write(11,"(8i6)") (Ipointer_local(i,n),i=0,n_nods-1)
            enddo
            write(11,*) "Faces"
            write(11,"(1i6,a)") n_faces
            do n = 0,nelem_in_proc(proc)-1
                write(11,"(6i6)") (faces(n,i),i=0,5)
                write(11,"(6i6)") (mapping_faces(n,i),i=0,5)
            enddo
            write(11,*) "Edges"
            write(11,"(1i6,a)") n_edges
            do n = 0,nelem_in_proc(proc)-1
                write(11,"(12i6)") (edges(n,i),i=0,11)
                write(11,"(12i6)") (mapping_edges(n,i),i=0,11)
            enddo
            write(11,*) "Vertices"
            write(11,"(1i6,a)") n_vertices
            do n = 0,nelem_in_proc(proc)-1
                write(11,"(8i6)") (vertices(n,i),i=0,7)
            enddo
            write(11,*) "Vertex <-> global nodes"
            do n = 0,n_vertices-1
                write(11,*) vertex_to_glob(n)
            end do
            if(solid_fluid) then
                write(11,*)
                write(11,*) SF_present_local
                if(SF_present_local)then
                    write(11,*) "Solid Fluid properties"
                    write(11,*) "Solid-Fluid Faces"
                    write(11,*) SF_n_faces
                    write(11,*) "4 Edges for each face"
                    do n = 0,SF_n_faces-1
                        write(11,"(4i6)") (SF_Face_Near_Edges(n,i),i=0,3)
                        write(11,"(4i6)") (SF_Face_Near_Edges_Orient(n,i),i=0,3)
                    enddo
                    write(11,*) "4 Vertices for each face"
                    do n = 0,SF_n_faces-1
                        write(11,"(4i6)") (SF_Face_Near_Vertices(n,i),i=0,3)
                    enddo
                    write(11,*) "Glob number of fluid and solid faces and orientation of solid compared to fluid"
                    do n = 0,SF_n_faces-1
                        write(11,"(3i6)") SF_Faces(n,0),SF_Faces(n,1),SF_Face_Orient(n)
                    enddo
                    write(11,*) "Glob number of fluid and solid edges and orientation of solid compared to fluid"
                    write(11,*) SF_n_edges
                    do n = 0,SF_n_edges-1
                        write(11,"(3i6)") SF_Edges(n,0),SF_Edges(n,1),SF_mapping_edges(n)
                    enddo
                    write(11,*) "Glob number SF vertices"
                    write(11,*) SF_n_vertices
                    do n = 0,SF_n_vertices-1
                        write(11,"(2i6)") SF_Vertices(n,0),SF_Vertices(n,1)
                    enddo
                end if
            endif
            if(Neumann_present) then
                write(11,*)
                write(11,*) Neu_present_local
                if(Neu_present_local)then
                    write(11,*) "Neumann interface properties"
                    write(11,*) "Neumann Faces"
                    write(11,*) Neu_n_faces
                    write(11,*) "4 Edges for each face"
                    do n = 0,Neu_n_faces-1
                        write(11,"(4i6)") (Neu_Face_Near_Edges(n,i),i=0,3)
                        write(11,"(4i6)") (Neu_Face_Near_Edges_Orient(n,i),i=0,3)
                    enddo
                    write(11,*) "4 Vertices for each face"
                    do n = 0,Neu_n_faces-1
                        write(11,"(4i6)") (Neu_Face_Near_Vertices(n,i),i=0,3)
                    enddo
                    write(11,*) "Glob number Neumann faces"
                    do n = 0,Neu_n_faces-1
                        write(11,"(i6)") Neu_Faces(n)
                    enddo
                    write(11,*) "Glob number Neumann edges"
                    write(11,*) Neu_n_edges
                    do n = 0,Neu_n_edges-1
                        write(11,"(i6)") Neu_Edges(n)
                    enddo
                    write(11,*) "Glob number Neumann vertices"
                    write(11,*) Neu_n_vertices
                    do n = 0,Neu_n_vertices-1
                        write(11,"(i6)") Neu_Vertices(n)
                    enddo
                end if
            endif
            write(11,*)

            write(11,*) "Communications between procs."
            write(11,"(1i6,a)") nproc, "  #  Nb of procs"
            do n = 0,nproc-1
                write(11,"(3i6)") nf_shared(n),ne_shared(n),nv_shared(n)
                if(SF_present_local)then
                    write(11,"(3i6)") SF_nf_shared(n), SF_ne_shared(n),SF_nv_shared(n)
                end if
                if(Neu_present_local)then
                    write(11,"(3i6)")  Neu_ne_shared(n),Neu_nv_shared(n)
                end if
                do nf = 0,nf_shared(n)-1
                    write(11,"(2i6)") faces_shared(n,nf),mapping_faces_shared(n,nf)
                enddo
                do ne = 0,ne_shared(n)-1
                    write(11,"(2i6)") edges_shared(n,ne),mapping_edges_shared(n,ne)
                enddo
                do nv = 0,nv_shared(n)-1
                    write(11,"(2i6)") vertices_shared(n,nv)
                enddo
                if(SF_present_local)then
                    do nf = 0,SF_nf_shared(n)-1
                        write(11,"(1i6)") SF_faces_shared(n,nf)
                    enddo
                    do ne = 0,SF_ne_shared(n)-1
                        write(11,"(2i6)") SF_edges_shared(n,ne),SF_mapping_edges_shared(n,ne)
                    enddo
                    do nv = 0,SF_nv_shared(n)-1
                        write(11,"(2i6)") SF_vertices_shared(n,nv)
                    enddo
                end if
                if(Neu_present_local)then
                    do ne = 0,Neu_ne_shared(n)-1
                        write(11,"(2i6)") Neu_edges_shared(n,ne),Neu_mapping_edges_shared(n,ne)
                    enddo
                    do nv = 0,Neu_nv_shared(n)-1
                        write(11,"(i6)") Neu_vertices_shared(n,nv)
                    enddo
                end if

            enddo
            close(11)


            ! - deallocations on each proc
            deallocate(Ipointer_local)
            deallocate(node_loc2glob)
            deallocate(N_valid_Vertex)
            deallocate(vertices)
            deallocate(vertex_to_glob)
            deallocate(faces)
            deallocate(mapping_faces)
            deallocate(faces_shared)
            deallocate(mapping_faces_shared)
            deallocate(nf_shared)
            deallocate(edges)
            deallocate(mapping_edges)
            deallocate(Elem_Edge_ref)
            deallocate(edges_shared)
            deallocate(mapping_edges_shared)
            deallocate(ne_shared)
            deallocate(vertices_shared)
            deallocate(nv_shared)

            ! deallocations if SF objects present locally
            if(SF_present_local)then
                deallocate(SF_faces)
                deallocate(SF_nf_shared)
                deallocate(SF_face_orient)
                deallocate(SF_faces_shared)
                deallocate(SF_edges)
                deallocate(SF_mapping_edges)
                deallocate(SF_ne_shared)
                deallocate(SF_edges_shared)
                deallocate(SF_mapping_edges_shared)
                deallocate(SF_vertices)
                deallocate(SF_nv_shared)
                deallocate(SF_vertices_shared)
                deallocate(SF_local_to_global_faces)
                deallocate(SF_Face_Near_Edges)
                deallocate(SF_Face_Near_Edges_Orient)
                deallocate(SF_Face_Near_Vertices)
            end if
            if(Neu_present_local)then
                deallocate(Neu_faces)
                deallocate(Neu_edges)
                deallocate(Neu_ne_shared)
                deallocate(Neu_edges_shared)
                deallocate(Neu_mapping_edges_shared)
                deallocate(Neu_vertices)
                deallocate(Neu_nv_shared)
                deallocate(Neu_vertices_shared)
                deallocate(Neu_local_to_global_faces)
                deallocate(Neu_Face_Near_Edges)
                deallocate(Neu_Face_Near_Edges_Orient)
                deallocate(Neu_Face_Near_Vertices)
            end if


            !- end of the loop on processors -!
        end do

        deallocate(which_points_inside)
        !---------------------------------------
        deallocate(elem_solid)
        deallocate(elem_contact)
        deallocate(nodes_nature)
        deallocate(Material)
        deallocate(Ipointer)
        if(solid_fluid)then
            deallocate(SF_global_to_local_edges)
            deallocate(SF_global_to_local_vertices)
        end if
        if(Neumann_present)then
            deallocate(Neu_global_to_local_edges)
            deallocate(Neu_global_to_local_vertices)
        end if


        write(*,*)
        write(*,*) "****************************************"
        write(*,*) "      --- NOW SEM CAN BE USED.. ---"
        write(*,*) "****************************************"


    contains

        !---------------------------------------------------------------------
        !-------------------------------------------------------------------
        subroutine indices_modif(n_nods,n_SF_nodes,SF_n_global_faces,Nodes_on_SF,    &
            SF_global_faces,Ipointer)

            implicit none
            integer, intent(in)   :: n_nods,n_SF_nodes,SF_n_global_faces
            integer, dimension(0:n_SF_nodes-1,0:1), intent(in)  :: Nodes_on_SF
            type(SF_face), dimension(0:SF_n_global_faces-1), intent(in) ::   &
                SF_global_faces
            integer, dimension(0:,0:), intent(inout)   :: Ipointer

            integer   :: i,j,n,nf

            do nf = 0, SF_n_global_faces-1
                n = SF_global_faces(nf)%elem(0)
                do i = 0,n_nods-1
                    do j = 0, n_SF_nodes-1
                        if(Ipointer(i,n) == Nodes_on_SF(j,0))then
                            Ipointer(i,n) = Nodes_on_SF(j,1)
                            exit
                        end if
                    end do
                end do
            end do

        end subroutine indices_modif
        !-----------------------------------------

    end subroutine gen_mesh

    !-------------------------------------------------

end module mesh2spec
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
