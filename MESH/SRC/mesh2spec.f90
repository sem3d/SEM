module mesh2spec

    use mesh_properties
    use sets
    use partition_mesh
    use solid_fluid
    use local_mesh_properties
    use neumann
    use free_surf_fluid

    type :: local_sf_info
        logical :: present_local
        integer :: n_faces
        integer :: n_edges
        integer :: n_vertices

        integer, allocatable, dimension(:)   :: ne_shared
        integer, allocatable, dimension(:)   :: nv_shared
        integer, allocatable, dimension(:)   :: nf_shared
        integer, allocatable, dimension(:)   :: face_orient
        integer, allocatable, dimension(:)   :: mapping_edges

        integer, allocatable, dimension(:,:) :: Faces
        integer, allocatable, dimension(:,:) :: Edges
        integer, allocatable, dimension(:,:) :: Vertices
        integer, allocatable, dimension(:,:) :: Face_Near_Edges
        integer, allocatable, dimension(:,:) :: Face_Near_Edges_Orient
        integer, allocatable, dimension(:,:) :: Face_Near_Vertices
        integer, allocatable, dimension(:,:) :: faces_shared
        integer, allocatable, dimension(:,:) :: edges_shared
        integer, allocatable, dimension(:,:) :: vertices_shared
        integer, allocatable, dimension(:,:) :: mapping_edges_shared

        integer, allocatable, dimension(:)   :: local_to_global_faces

    end type local_sf_info

    type :: local_neu_info
        logical  :: present_local
        integer  :: n_faces
        integer  :: n_edges
        integer  :: n_vertices

        integer, allocatable, dimension(:)   :: ne_shared
        integer, allocatable, dimension(:)   :: nv_shared
        integer, allocatable, dimension(:)   :: faces
        integer, allocatable, dimension(:)   :: edges
        integer, allocatable, dimension(:)   :: vertices

        integer, allocatable, dimension(:,:) :: Face_Near_edges
        integer, allocatable, dimension(:,:) :: Face_Near_Edges_Orient
        integer, allocatable, dimension(:,:) :: Face_Near_Vertices
        integer, allocatable, dimension(:,:) :: edges_shared
        integer, allocatable, dimension(:,:) :: vertices_shared
        integer, allocatable, dimension(:,:) :: mapping_edges_shared

        integer, allocatable, dimension(:)   :: local_to_global_faces

    end type local_neu_info

    type :: shared_info
        integer, dimension(:), allocatable :: nf
        integer, dimension(:), allocatable :: ne
        integer, dimension(:), allocatable :: nv
        integer, dimension(:,:), allocatable :: vertices
        integer, dimension(:,:), allocatable :: faces
        integer, dimension(:,:), allocatable :: edges
        integer, dimension(:,:), allocatable :: mapping_edges
        integer, dimension(:,:), allocatable :: mapping_faces
     end type shared_info
contains

    !---------------------
    subroutine gen_mesh(nproc)

        implicit none

        !- general --
        integer   :: n,i,j,num, &
            nel,nf,  &
            nv0,nv1,n0,n1
        character(len=14) :: meshfilename

        !- mesh --
        integer   :: n_nods, n_elem, n_blocks, n_points ! n_blocks corresponds
        ! to the number of different materials (cf Cubit)
        integer, allocatable :: material(:), Ipointer(:,:)
        real, allocatable  :: xco(:),yco(:),zco(:), Gcoord(:,:)
        integer  :: n_old_points

        !- nearest elements
        type(near_node), dimension(:), allocatable  :: initnode,prevnode,currnode
        type(near_elem), dimension(:), allocatable  :: near_elem_set,     &
            SF_edges_near_elem, Neu_edges_near_elem
        type(near_proc), dimension(:), allocatable  :: elem_near_proc

        !- materials --
        character, dimension(:), allocatable  ::  tabmat
        integer      :: n_SF_nodes
        logical    :: any_fluid, all_fluid, solid_fluid
        integer, allocatable :: nodes_nature(:)
        logical, allocatable :: elem_contact(:), elem_solid(:)
        logical, allocatable :: elem_fluid_dirich(:)
        !- solid/fluid --
        integer, allocatable :: Nodes_on_SF(:,:), SF_object_n_faces(:),             &
            Elem_Edge_Ref(:),     &
            SF_global_node_to_vertex(:),  &
            SF_global_to_local_edges(:),                        &
            SF_global_to_local_vertices(:)
        logical, allocatable :: SF_object(:,:), SF_global_face_present(:),          &
            SF_global_edge_present(:)
        type(Face_type), allocatable :: SF_object_Face(:,:)
        type(SF_face), allocatable :: SF_global_faces(:)
        type(SF_edge), allocatable :: SF_global_edges(:)
        type(SF_vertex), allocatable :: SF_global_vertices(:)
        integer   :: SF_n_global_faces, SF_n_global_vertices, SF_n_global_edges
        !- interfaces: Neumann BC
        integer   :: Neu_n_global_faces, Neu_n_global_vertices, Neu_n_global_edges
        logical, dimension(:), allocatable  :: Neu_object
        integer, dimension(:), allocatable :: Neu_object_n_faces,              &
            Neu_global_node_to_vertex,                    &
            Neu_global_to_local_edges(:),                   &
            Neu_global_to_local_vertices(:)
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
        integer, pointer, dimension(:) :: dxadj, dxadjncy
        integer, allocatable, dimension(:) :: part
        !- local meshes
        type(souvenir), dimension(:), allocatable :: memory
        type(process_obj), dimension(:), allocatable :: MemorySF_F,MemorySF_E,MemorySF_V
        type(proc_obj), dimension(:), allocatable :: MemoryNeu
        logical, allocatable :: which_points_inside(:)
        integer  :: proc, n_vertices, n_faces, n_edges, n_points_local
        integer, dimension(:), allocatable :: counter, elem_glob2loc, glob2loc,   &
            N_valid_Vertex, vertex_to_glob,        &
            nelem_in_proc,  node_loc2glob
        integer, dimension(:,:), allocatable :: which_elem_in_proc,  &
            Ipointer_local, vertices
        integer, dimension(:,:), allocatable :: faces,  mapping_faces
        integer, dimension(:,:), allocatable :: edges, mapping_edges

        ! Stuff shared between procs
        type(shared_info) :: shared
        ! local Neumann objects
        type(local_neu_info) :: Neu

        !- local solid-fluid objects
        type(local_sf_info) :: SF




        !--------------------------------------------------------------
        !- INITIALIZATIONS OF DIFFERENT MESH AND MATERIAL PROPERTIES
        call mesh_init_3D(n_nods,n_points,n_elem,n_blocks,                  &
            xco,yco,zco,Ipointer,Material,tabmat,             &
            Neu_n_global_faces,n_PW,Faces_Neumann,Faces_PW)
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

        call part_mesh_3D(n_nods,n_elem,n_points,Ipointer,nproc,dxadj,dxadjncy,part)
        write(*,*) "  --> Partition done."

        allocate(elem_near_proc(0:n_elem-1))
        call find_near_proc(n_elem,near_elem_set,part,elem_near_proc)
         

        !-----------------------------------------------------------
        !-  BOUNDARY CONDITIONS
        !- Dirichlet boundary condition in fluid: free surface
        if(any_fluid)then
            allocate(elem_fluid_dirich(0:n_elem-1))
            call find_fluid_elem_freesurf(n_elem,dxadj,dxadjncy,   &
                Ipointer,elem_solid,elem_fluid_dirich)
        end if


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

            !- now we construct GLOBAL Neumann faces; it shall be
            !-   of great help when constructing local Neumann faces
            allocate(Neu_global_faces(0:Neu_n_global_faces-1))
            allocate(Neu_global_face_present(0:Neu_n_global_faces-1))
            call Neu_global_faces_construct(n_elem,Neu_n_global_faces,      &
                Neu_object,Neu_object_face,Neu_object_n_faces,  &
                Neu_global_faces)


            !- now we construct GLOBAL Neumann vertices
            allocate(Neu_global_node_to_vertex(0:n_points-1))
            call Neu_global_vertices_construct(Neu_n_global_faces,         &
                Ipointer,Neu_global_faces,Neu_global_node_to_vertex,    &
                Neu_n_global_vertices,Neu_global_vertices)

            !- now we construct GLOBAL Neumann edges
            call Neu_global_edges_construct(Neu_n_global_faces,Ipointer,  &
                Neu_global_node_to_vertex, Neu_global_faces,         &
                Neu_n_global_edges,Neu_global_edges)

            allocate(Neu_global_edge_present(0:Neu_n_global_edges-1))


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
            end do

            !- which procs do Neumann vertices belong to?
            call Neu_vertices_proc_belong(nproc,Neu_n_global_vertices,part,    &
                initnode,Neu_global_vertices)

            !- which procs do Neumann edges belong to?
            call Neu_edges_proc_belong(nproc,Neu_n_global_edges,     &
                Neu_global_vertices,Neu_global_edges)


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
                 !  call list_entity(SF_edges_near_elem(i)%ptr) ; read*
            end do

            !- if Solid-Fluid: doubled nodes -> change of indices
            call indices_modif(n_nods,n_SF_nodes,SF_n_global_faces,Nodes_on_SF,    &
                SF_global_faces,Ipointer)


            !- which procs do SF vertices belong to?
            call SF_vertices_proc_belong(nproc,SF_n_global_vertices,part,    &
                elem_solid,initnode,SF_global_vertices)

            !- which procs do SF edges belong to?
            call SF_edges_proc_belong(nproc,part,elem_solid,SF_n_global_edges, &
                 SF_global_vertices,SF_edges_near_elem,SF_global_edges)


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
            nullify(memory(nel)%rank)
            if(part(nel) /= nproc-1)then
                if(elem_near_proc(nel)%nb > 0)  &
                    allocate(memory(nel)%rank(0:elem_near_proc(nel)%nb-1))
                do proc = 0,elem_near_proc(nel)%nb-1
                    allocate(memory(nel)%rank(proc)%E(0:25))
                enddo
            endif
        enddo

        if(Neumann_present)then    !-- Neumann local part
            !- looking for reference elements for Neumann edges
            call Neu_edges_reference_elem(nproc,Neu_n_global_edges,       &
                Elem_glob2loc,part,Ipointer,Neu_global_vertices,   &
                Neu_edges_near_elem,Neu_global_edges)


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
            allocate(MemorySF_F(0:SF_n_global_faces-1))
            allocate(MemorySF_E(0:SF_n_global_edges-1))
            allocate(MemorySF_V(0:SF_n_global_vertices-1))

            do i = 0,SF_n_global_faces-1
                allocate(MemorySF_F(i)%procs(0:0,0:0))
                MemorySF_F(i)%procs(0:0,0:0) = -1
            end do
            do i = 0,SF_n_global_edges-1
                if(SF_global_edges(i)%local_proc%nr > 0)then
                    allocate(MemorySF_E(i)%procs(0:SF_global_edges(i)%local_proc%nr-1,0:SF_global_edges(i)%local_proc%nr-1))
                    MemorySF_E(i)%procs(0:,0:) = -1
                end if
            end do
            do i = 0,SF_n_global_vertices-1
                if(SF_global_vertices(i)%local_proc%nr > 0)then
                    allocate(MemorySF_V(i)%procs(0:SF_global_vertices(i)%local_proc%nr-1,0:SF_global_vertices(i)%local_proc%nr-1))
                    MemorySF_V(i)%procs(0:,0:) = -1
                end if
            end do

        end if   !-  temporary end of local solid/fluid properties


        !- global nodes
        allocate(glob2loc(0:n_points-1))
        allocate(which_points_inside(0:n_points-1))

        !- we consider each processor, and create all structures for the meshfiles
        call system("rm -f mesh4spec.????.h5")


        meshfilename(1:12) = "mesh4spec.h5"
        call write_global_mesh_file_h5(meshfilename, n_points, n_elem, n_nods, Gcoord, Ipointer, material );

        meshfilename(1:10) = "mesh4spec."

        !-------------------------------------------------
        !-  CREATION OF SEM-style files : LOOP ON PROCS
        write(*,*) "  --> Loop on procs!"
        do proc = 0,nproc-1
            write(*,*) "*---------------------*"
            write(*,*) "   --> Proc. # --  ", proc

            SF%present_local = .false.   !- solid_fluid properties
            Neu%present_local = .false.   !- Neumann properties

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
                        Neu%present_local = .true.
                        exit
                    end if
                end do
                if(Neu%present_local) write(*,*) "    --> Neumann objects present in proc. ",proc
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
                        SF%present_local = .true.
                        exit
                    end if
                end do
                if(SF%present_local) write(*,*) "    --> Solid/fluid objects present in proc. ",proc
                SF_global_to_local_edges(0:) = -1
                SF_global_to_local_vertices(0:) = -1
            end if


            !- Now let's turn to faces: 6 associated  to each element
            !  Are also determined: the faces shared between different procs
            allocate(faces(0:nelem_in_proc(proc)-1,0:5))
            allocate(mapping_faces(0:nelem_in_proc(proc)-1,0:5))
            allocate(shared%faces(0:nproc-1,0:6*nelem_in_proc(proc)-1))
            allocate(shared%mapping_faces(0:nproc-1,0:6*nelem_in_proc(proc)-1))
            allocate(shared%nf(0:nproc-1))
            call local_faces_construct(proc,part,nelem_in_proc,which_elem_in_proc,Elem_glob2loc,   &
                Ipointer,dxadj,dxadjncy,elem_near_proc,n_faces,faces,mapping_faces,  &
                shared%faces,shared%mapping_faces,shared%nf,memory)


            !- associating to each element: 12 edges
            !- defining the edges shared with other processors
            allocate(edges(0:nelem_in_proc(proc)-1,0:11))
            allocate(mapping_edges(0:nelem_in_proc(proc)-1,0:11))
            allocate(Elem_Edge_Ref(0:12*nelem_in_proc(proc)-1))
            allocate(shared%edges(0:nproc-1,0:12*nelem_in_proc(proc)-1))
            allocate(shared%mapping_edges(0:nproc-1,0:12*nelem_in_proc(proc)-1))
            allocate(shared%ne(0:nproc-1))
            call local_edges_construct(nproc,proc,part,nelem_in_proc,Elem_glob2loc,   &
                which_elem_in_proc,Ipointer,near_elem_set,elem_near_proc,n_edges,     &
                edges,mapping_edges,shared%edges,shared%mapping_edges,       &
                shared%ne,Elem_edge_ref,memory)

            !- defining vertices shared with other processors
            allocate(shared%vertices(0:nproc-1,0:7*nelem_in_proc(proc)-1))
            allocate(shared%nv(0:nproc-1))
            call local_vertices_comm(n_vertices,proc,nproc,nelem_in_proc,part,   &
                vertex_to_glob,node_loc2glob,which_elem_in_proc,Ipointer,  &
                vertices,near_elem_set,elem_near_proc,shared%vertices,shared%nv,memory)


            !- now we include eventual Neumann objects
            if(Neu%present_local)then
                write(*,*) "    --> Neumann local meshing properties."
                Neu%n_faces = 0 ; Neu%n_edges = 0 ; Neu%n_vertices = 0
                !- local Neumann faces ---------------------
                ! - how many?
                do nf = 0, Neu_n_global_faces-1
                    if(part(Neu_global_faces(nf)%elem) == proc)then
                        Neu_global_face_present(nf) = .true.
                        Neu%n_faces = Neu%n_faces+1
                    end if
                end do

                allocate(Neu%faces(0:Neu%n_faces-1))
                allocate(Neu%local_to_global_faces(0:Neu%n_faces-1))

                call Neu_local_faces_construct(proc,nproc,Neu%n_faces,      &
                    Neu_n_global_faces,part,Elem_glob2loc,               &
                    Neu_global_face_present,Neu_global_faces,            &
                    faces,Neu%faces,Neu%local_to_global_faces)

                !- local Neumann edges  ---------------------
                ! how many?
                do i = 0, Neu_n_global_edges-1
                    if((.not. Neu_global_edges(i)%local(proc))) cycle
                    Neu%n_edges = Neu%n_edges+1
                end do

                allocate(Neu%edges(0:Neu%n_edges-1))
                allocate(Neu%ne_shared(0:nproc-1))
                allocate(Neu%edges_shared(0:nproc-1,0:Neu%n_edges-1))
                allocate(Neu%mapping_edges_shared(0:nproc-1,0:Neu%n_edges-1))

                call Neu_local_edges_construct(proc,nproc,Neu%n_edges,Neu_n_global_edges, &
                    part,glob2loc,Neu_global_edges,Neu_global_vertices,edges,              &
                    Ipointer_local,Ipointer,which_elem_in_proc,Neu_global_to_local_edges,  &
                    Neu%edges,Neu%edges_shared, Neu%mapping_edges_shared,Neu%ne_shared,    &
                    MemoryNeu)

                !- local Neumann vertices
                ! how many?
                do i = 0, Neu_n_global_vertices-1
                    if((.not. Neu_global_vertices(i)%local(proc))) cycle
                    Neu%n_vertices = Neu%n_vertices+1
                end do

                allocate(Neu%vertices(0:Neu%n_vertices-1))
                allocate(Neu%nv_shared(0:nproc-1))
                allocate(Neu%vertices_shared(0:nproc-1,0:Neu%n_vertices-1))

                call Neu_local_vertices_construct(proc,nproc,Neu%n_vertices,     &
                    Neu_n_global_vertices,glob2loc,N_valid_Vertex,              &
                    Neu_global_vertices,Neu_global_to_local_vertices,           &
                    Neu%vertices,Neu%vertices_shared,Neu%nv_shared,MemoryNeu)

                !- finally: Neumann edges and vertices associated to Neumann faces
                !   -> important for the calculation of properties of normals to faces.
                allocate(Neu%Face_Near_Edges(0:Neu%n_faces-1,0:3))
                allocate(Neu%Face_Near_Edges_Orient(0:Neu%n_faces-1,0:3))
                allocate(Neu%Face_Near_Vertices(0:Neu%n_faces-1,0:3))

                call Neu_faces_to_EV(proc,Neu%n_faces,Neu%local_to_global_faces,  &
                    Neu_global_to_local_edges,Neu_global_to_local_vertices,         &
                    Neu_global_faces,Neu_global_edges,Neu%faces,mapping_edges,      &
                    Neu%Face_Near_Edges,Neu%Face_Near_Edges_Orient,Neu%Face_Near_Vertices)


            end if   !-   end of local Neumann properties



            !- now we include eventual SF objects
            if(SF%present_local)then
                SF%n_faces = 0 ; SF%n_edges = 0 ; SF%n_vertices = 0
                write(*,*) "    --> Solid/fluid local meshing properties."
                !- local SF faces ---------------------
                ! - how many?
                do nf = 0, SF_n_global_faces-1
                    if(part(SF_global_faces(nf)%elem(0)) == proc .or.   &
                        part(SF_global_faces(nf)%elem(1)) == proc)then
                        SF_global_face_present(nf) = .true.
                        SF%n_faces = SF%n_faces+1
                    end if
                end do

                allocate(SF%faces(0:SF%n_faces-1,0:1))
                allocate(SF%nf_shared(0:nproc-1))
                allocate(SF%faces_shared(0:nproc-1,0:SF%n_faces-1))
                allocate(SF%face_orient(0:SF%n_faces-1))
                allocate(SF%local_to_global_faces(0:SF%n_faces-1))

                call SF_local_faces_construct(proc,nproc,SF%n_faces,SF_n_global_faces,  &
                    part,Elem_glob2loc,SF_global_face_present,SF_global_faces,        &
                    SF_object_n_faces,SF_object_face,faces,SF%faces,SF%faces_shared,  &
                    SF%nf_shared,SF%face_orient,SF%local_to_global_faces,MemorySF_F)

                !- local SF edges  ---------------------
                ! how many?
                do i = 0, SF_n_global_edges-1
                    if(already_in(proc,SF_global_edges(i)%local_proc%elem,    &
                       SF_global_edges(i)%local_proc%nr)>= 0)                 &
                           SF%n_edges = SF%n_edges+1
                end do

                allocate(SF%edges(0:SF%n_edges-1,0:1))
                allocate(SF%mapping_edges(0:SF%n_edges-1))
                allocate(SF%ne_shared(0:nproc-1))
                allocate(SF%edges_shared(0:nproc-1,0:SF%n_edges-1))
                allocate(SF%mapping_edges_shared(0:nproc-1,0:SF%n_edges-1))

                call SF_local_edges_construct(proc,nproc,SF%n_edges,SF_n_global_edges, &
                    part,glob2loc,SF_global_edges,SF_global_vertices,edges,             &
                    Ipointer_local,Ipointer,which_elem_in_proc,SF_global_to_local_edges,&
                    SF%edges,SF%edges_shared, SF%mapping_edges_shared,SF%ne_shared,     &
                    SF%mapping_edges,MemorySF_E)

                !- local SF vertices
                ! how many?
                do i = 0, SF_n_global_vertices-1
                    if(already_in(proc,SF_global_vertices(i)%local_proc%elem,    &
                       SF_global_vertices(i)%local_proc%nr) >= 0)                &
                        SF%n_vertices = SF%n_vertices+1
                end do

                allocate(SF%vertices(0:SF%n_vertices-1,0:1))
                allocate(SF%nv_shared(0:nproc-1))
                allocate(SF%vertices_shared(0:nproc-1,0:SF%n_vertices-1))

                call SF_local_vertices_construct(proc,nproc,SF%n_vertices,            &
                    SF_n_global_vertices,glob2loc,N_valid_Vertex,SF_global_vertices, &
                    SF_global_to_local_vertices,SF%vertices,SF%vertices_shared,      &
                    SF%nv_shared,MemorySF_V)


                !- finally: SF edges and vertices associated to SF faces -> important for
                !       the calculation of properties of normals to faces.
                allocate(SF%Face_Near_Edges(0:SF%n_faces-1,0:3))
                allocate(SF%Face_Near_Edges_Orient(0:SF%n_faces-1,0:3))
                allocate(SF%Face_Near_Vertices(0:SF%n_faces-1,0:3))

                call SF_faces_to_EV(proc,SF%n_faces,SF%local_to_global_faces,    &
                    SF_global_to_local_edges,SF_global_to_local_vertices,          &
                    SF_global_faces,SF_global_edges,SF%faces,mapping_edges,        &
                    SF%Face_Near_Edges,SF%Face_Near_Edges_Orient,SF%Face_Near_Vertices)


            end if   !-   end of local SF properties


            !- writing the meshfile
            write(meshfilename(11:14), '(i4.4)') proc

!            call write_mesh_file_txt(meshfilename,solid_fluid,all_fluid,Neumann_present, &
!                n_elem,n_points,n_points_local,n_blocks,n_edges,n_faces,n_nods,n_vertices, &
!                SF, Neu, shared, &
!                material, elem_solid, Ipointer_local, &
!                faces, mapping_faces, edges, mapping_edges, vertices, &
!                node_loc2glob,Gcoord,vertex_to_glob, &
!                which_elem_in_proc,nelem_in_proc,proc,nproc)

            call write_mesh_file_h5(meshfilename//".h5",solid_fluid,all_fluid,Neumann_present, &
                n_elem,n_points,n_points_local,n_blocks,n_edges,n_faces,n_nods,n_vertices, &
                SF, Neu, shared, &
                material, elem_solid,elem_fluid_dirich, Ipointer_local, &
                faces, mapping_faces, edges, mapping_edges, vertices, &
                node_loc2glob,Gcoord,vertex_to_glob, &
                which_elem_in_proc,nelem_in_proc,proc,nproc)


            ! - deallocations on each proc
            deallocate(Ipointer_local)
            deallocate(node_loc2glob)
            deallocate(N_valid_Vertex)
            deallocate(vertices)
            deallocate(vertex_to_glob)
            deallocate(faces)
            deallocate(mapping_faces)
            deallocate(shared%faces)
            deallocate(shared%mapping_faces)
            deallocate(shared%nf)
            deallocate(edges)
            deallocate(mapping_edges)
            deallocate(Elem_Edge_ref)
            deallocate(shared%edges)
            deallocate(shared%mapping_edges)
            deallocate(shared%ne)
            deallocate(shared%vertices)
            deallocate(shared%nv)

            ! deallocations if SF objects present locally
            if(SF%present_local)then
                deallocate(SF%faces)
                deallocate(SF%nf_shared)
                deallocate(SF%face_orient)
                deallocate(SF%faces_shared)
                deallocate(SF%edges)
                deallocate(SF%mapping_edges)
                deallocate(SF%ne_shared)
                deallocate(SF%edges_shared)
                deallocate(SF%mapping_edges_shared)
                deallocate(SF%vertices)
                deallocate(SF%nv_shared)
                deallocate(SF%vertices_shared)
                deallocate(SF%local_to_global_faces)
                deallocate(SF%Face_Near_Edges)
                deallocate(SF%Face_Near_Edges_Orient)
                deallocate(SF%Face_Near_Vertices)
            end if
            if(Neu%present_local)then
                deallocate(Neu%faces)
                deallocate(Neu%edges)
                deallocate(Neu%ne_shared)
                deallocate(Neu%edges_shared)
                deallocate(Neu%mapping_edges_shared)
                deallocate(Neu%vertices)
                deallocate(Neu%nv_shared)
                deallocate(Neu%vertices_shared)
                deallocate(Neu%local_to_global_faces)
                deallocate(Neu%Face_Near_Edges)
                deallocate(Neu%Face_Near_Edges_Orient)
                deallocate(Neu%Face_Near_Vertices)
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
        if(any_fluid) deallocate(elem_fluid_dirich)
        if(solid_fluid)then
            deallocate(SF_global_to_local_edges)
            deallocate(SF_global_to_local_vertices)
        end if
        if(Neumann_present)then
            deallocate(Neu_global_to_local_edges)
            deallocate(Neu_global_to_local_vertices)
        end if

        do nel = 0,n_elem-1
            if(part(nel) /= nproc-1)then
                if (associated(memory(nel)%rank)) then
                    do proc = 0,elem_near_proc(nel)%nb-1
                        deallocate(memory(nel)%rank(proc)%E)
                    enddo
                    deallocate(memory(nel)%rank)
                endif
            endif
        enddo
        deallocate(memory)



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

    subroutine write_mesh_file_txt(meshfilename,solid_fluid,all_fluid,Neumann_present, &
        n_elem,n_points,n_points_local,n_blocks,n_edges,n_faces,n_nods,n_vertices, &
        SF, Neu, shared, &
        material, elem_solid, Ipointer_local, &
        faces, mapping_faces, edges, mapping_edges, vertices, &
        node_loc2glob,Gcoord,vertex_to_glob, &
        which_elem_in_proc,nelem_in_proc,proc,nproc)
        implicit none
        character(len=14),intent(in) :: meshfilename
        integer, intent(in)          :: proc, nproc
        logical, intent(in)          :: all_fluid, solid_fluid
        logical, intent(in)          :: Neumann_present
        integer, intent(in)          :: n_elem,n_points, n_points_local, n_blocks
        integer, intent(in)          :: n_edges, n_faces, n_nods, n_vertices
        real, intent(in), dimension(0:n_points-1,0:2)      :: Gcoord
        integer, intent(in), dimension(0:n_points_local-1) :: node_loc2glob
        integer, intent(in), dimension(0:nproc-1)          :: nelem_in_proc
        integer, intent(in), dimension(0:nproc-1,0:maxval(nelem_in_proc)-1) :: which_elem_in_proc
        integer, intent(in), dimension(0:n_elem-1) :: material
        logical, intent(in), dimension(0:n_elem-1) :: elem_solid
        integer, intent(in), dimension(0:n_nods-1,0:nelem_in_proc(proc)-1) :: Ipointer_local
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:5)  :: faces
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:5)  :: mapping_faces
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:11) :: edges
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:11) :: mapping_edges
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:7)  :: vertices
        integer, intent(in), dimension(0:8*nelem_in_proc(proc)-1)    :: vertex_to_glob
        type(local_sf_info), intent(in) :: SF
        type(local_neu_info), intent(in) :: Neu
        type(shared_info), intent(in) :: shared
        !
        logical :: curve
        integer :: n_dim, i, i_count, n, ne, nv, nf, nel
        !
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
            write(11,*) SF%present_local
            if(SF%present_local)then
                write(11,*) "Solid Fluid properties"
                write(11,*) "Solid-Fluid Faces"
                write(11,*) SF%n_faces
                write(11,*) "4 Edges for each face"
                do n = 0,SF%n_faces-1
                    write(11,"(4i6)") (SF%Face_Near_Edges(n,i),i=0,3)
                    write(11,"(4i6)") (SF%Face_Near_Edges_Orient(n,i),i=0,3)
                enddo
                write(11,*) "4 Vertices for each face"
                do n = 0,SF%n_faces-1
                    write(11,"(4i6)") (SF%Face_Near_Vertices(n,i),i=0,3)
                enddo
                write(11,*) "Glob number of fluid and solid faces and orientation of solid compared to fluid"
                do n = 0,SF%n_faces-1
                    write(11,"(3i6)") SF%Faces(n,0),SF%Faces(n,1),SF%Face_Orient(n)
                enddo
                write(11,*) "Glob number of fluid and solid edges and orientation of solid compared to fluid"
                write(11,*) SF%n_edges
                do n = 0,SF%n_edges-1
                    write(11,"(3i6)") SF%Edges(n,0),SF%Edges(n,1),SF%mapping_edges(n)
                enddo
                write(11,*) "Glob number SF vertices"
                write(11,*) SF%n_vertices
                do n = 0,SF%n_vertices-1
                    write(11,"(2i6)") SF%Vertices(n,0),SF%Vertices(n,1)
                enddo
            end if
        endif
        if(Neumann_present) then
            write(11,*)
            write(11,*) Neu%present_local
            if(Neu%present_local)then
                write(11,*) "Neumann interface properties"
                write(11,*) "Neumann Faces"
                write(11,*) Neu%n_faces
                write(11,*) "4 Edges for each face"
                do n = 0,Neu%n_faces-1
                    write(11,"(4i6)") (Neu%Face_Near_Edges(n,i),i=0,3)
                    write(11,"(4i6)") (Neu%Face_Near_Edges_Orient(n,i),i=0,3)
                enddo
                write(11,*) "4 Vertices for each face"
                do n = 0,Neu%n_faces-1
                    write(11,"(4i6)") (Neu%Face_Near_Vertices(n,i),i=0,3)
                enddo
                write(11,*) "Glob number Neumann faces"
                do n = 0,Neu%n_faces-1
                    write(11,"(i6)") Neu%Faces(n)
                enddo
                write(11,*) "Glob number Neumann edges"
                write(11,*) Neu%n_edges
                do n = 0,Neu%n_edges-1
                    write(11,"(i6)") Neu%Edges(n)
                enddo
                write(11,*) "Glob number Neumann vertices"
                write(11,*) Neu%n_vertices
                do n = 0,Neu%n_vertices-1
                    write(11,"(i6)") Neu%Vertices(n)
                enddo
            end if
        endif
        write(11,*)

        write(11,*) "Communications between procs."
        write(11,"(1i6,a)") nproc, "  #  Nb of procs"
        do n = 0,nproc-1
            write(11,"(3i6)") shared%nf(n),shared%ne(n),shared%nv(n)
            if(SF%present_local)then
                write(11,"(3i6)") SF%nf_shared(n), SF%ne_shared(n),SF%nv_shared(n)
            end if
            if(Neu%present_local)then
                write(11,"(3i6)")  Neu%ne_shared(n),Neu%nv_shared(n)
            end if
            do nf = 0,shared%nf(n)-1
                write(11,"(2i6)") shared%faces(n,nf),shared%mapping_faces(n,nf)
            enddo
            do ne = 0,shared%ne(n)-1
                write(11,"(2i6)") shared%edges(n,ne),shared%mapping_edges(n,ne)
            enddo
            do nv = 0,shared%nv(n)-1
                write(11,"(2i6)") shared%vertices(n,nv)
            enddo
            if(SF%present_local)then
                do nf = 0,SF%nf_shared(n)-1
                    write(11,"(1i6)") SF%faces_shared(n,nf)
                enddo
                do ne = 0,SF%ne_shared(n)-1
                    write(11,"(2i6)") SF%edges_shared(n,ne),SF%mapping_edges_shared(n,ne)
                enddo
                do nv = 0,SF%nv_shared(n)-1
                    write(11,"(2i6)") SF%vertices_shared(n,nv)
                enddo
            end if
            if(Neu%present_local)then
                do ne = 0,Neu%ne_shared(n)-1
                    write(11,"(2i6)") Neu%edges_shared(n,ne),Neu%mapping_edges_shared(n,ne)
                enddo
                do nv = 0,Neu%nv_shared(n)-1
                    write(11,"(i6)") Neu%vertices_shared(n,nv)
                enddo
            end if

        enddo
        close(11)
    end subroutine write_mesh_file_txt


   subroutine write_global_mesh_file_h5(meshfilename, n_points, n_elem, &
        n_nods, Gcoord, Ipointer, material)
        use hdf5
        use sem_hdf5
        implicit none
        character(len=*),intent(in) :: meshfilename
        integer, intent(in)          :: n_nods, n_points, n_elem
        real, intent(in), dimension(0:n_points-1,0:2)      :: Gcoord
        integer, intent(in), dimension(0:n_nods-1,0:n_elem-1) :: Ipointer
        integer, intent(in), dimension(0:n_elem-1) :: material

        integer :: hdferr, i
        integer(HID_T) :: fid
        integer, dimension(0:7,0:n_elem-1) :: Ipointer_vertices


        call init_hdf5()
        call h5fcreate_f(meshfilename, H5F_ACC_TRUNC_F, fid, hdferr)

        call write_attr_int(fid, "n_vertices", n_points)
        call write_attr_int(fid, "n_elements", n_elem)

        call write_dataset(fid, "local_nodes", transpose(Gcoord), hdferr)

        do i=0,n_elem-1
            Ipointer_vertices(0:7, i) = Ipointer(0:7, i)
        enddo

        call write_dataset(fid, "elements", Ipointer_vertices, hdferr)

        call write_dataset(fid, "material", material, hdferr)


        call h5fclose_f(fid, hdferr)


    end subroutine write_global_mesh_file_h5



    subroutine write_mesh_file_h5(meshfilename,solid_fluid,all_fluid,Neumann_present, &
        n_elem,n_points,n_points_local,n_blocks,n_edges,n_faces,n_nods,n_vertices, &
        SF, Neu, shared, &
        material, elem_solid,elem_fluid_dirich, Ipointer_local, &
        faces, mapping_faces, edges, mapping_edges, vertices, &
        node_loc2glob,Gcoord,vertex_to_glob, &
        which_elem_in_proc,nelem_in_proc,proc,nproc)
        use hdf5
        use sem_hdf5
        implicit none
        character(len=17),intent(in) :: meshfilename
        integer, intent(in)          :: proc, nproc
        logical, intent(in)          :: all_fluid, solid_fluid
        logical, intent(in)          :: Neumann_present
        integer, intent(in)          :: n_elem,n_points, n_points_local, n_blocks
        integer, intent(in)          :: n_edges, n_faces, n_nods, n_vertices
        real, intent(in), dimension(0:n_points-1,0:2)      :: Gcoord
        integer, intent(in), dimension(0:n_points_local-1) :: node_loc2glob
        integer, intent(in), dimension(0:nproc-1)          :: nelem_in_proc
        integer, intent(in), dimension(0:nproc-1,0:maxval(nelem_in_proc)-1) :: which_elem_in_proc
        integer, intent(in), dimension(0:n_elem-1) :: material
        logical, intent(in), dimension(0:n_elem-1) :: elem_solid, elem_fluid_dirich
        integer, intent(in), dimension(0:n_nods-1,0:nelem_in_proc(proc)-1) :: Ipointer_local
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:5)  :: faces
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:5)  :: mapping_faces
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:11) :: edges
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:11) :: mapping_edges
        integer, intent(in), dimension(0:nelem_in_proc(proc)-1,0:7)  :: vertices
        integer, intent(in), dimension(0:8*nelem_in_proc(proc)-1)    :: vertex_to_glob
        type(local_sf_info), intent(in) :: SF
        type(local_neu_info), intent(in) :: Neu
        type(shared_info), intent(in) :: shared
        !
        logical :: curve
        integer :: n_dim, i, i_count, n, ne, nv, nf, nel
        integer(HSIZE_T), dimension(2) :: dims
        real, allocatable, dimension(:,:) :: rtemp2
        integer, allocatable, dimension(:,:) :: itemp2
        integer(HID_T) :: dset_id, fid, proc_id
        integer :: hdferr
        character(Len=20) :: proc_grp
        !
        curve = .false.
        !
        call init_hdf5()
        !
        call h5fcreate_f(meshfilename, H5F_ACC_TRUNC_F, fid, hdferr)
        !
        !! Attributes
        call write_attr_int(fid, "ndim", 3)
        call write_attr_int(fid, "n_processors", nproc)
        call write_attr_int(fid, "n_materials", n_blocks)
        call write_attr_int(fid, "n_elements", nelem_in_proc(proc))
        call write_attr_int(fid, "n_faces", n_faces)
        call write_attr_int(fid, "n_edges", n_edges)
        call write_attr_int(fid, "n_vertices", n_vertices)
        call write_attr_bool(fid, "solid_fluid", solid_fluid)
        call write_attr_bool(fid, "solid_fluid_loc", SF%present_local)
        call write_attr_bool(fid, "all_fluid", all_fluid)
        if(solid_fluid .and. SF%present_local)then
            call write_attr_int(fid, "n_SF_faces", SF%n_faces)
            call write_attr_int(fid, "n_SF_edges", SF%n_edges)
            call write_attr_int(fid, "n_SF_vertices", SF%n_vertices)
        end if
        call write_attr_bool(fid, "neumann_present", neumann_present)
        call write_attr_bool(fid, "neumann_present_loc", Neu%present_local)
        call write_attr_bool(fid, "curve", curve)



        dims(1) = 3
        dims(2) = n_points_local

        ! Local nodes

        allocate(rtemp2(0:2,0:n_points_local-1))
        do n = 0,n_points_local-1
            i_count = node_loc2glob(n)
            rtemp2(0:2,n) = Gcoord(i_count,0:2)
        enddo
        call write_dataset(fid, "local_nodes", rtemp2, hdferr)
        deallocate(rtemp2)

        ! Material table and solid fluid flag
        dims(1) = 2
        dims(2) = nelem_in_proc(proc)

        allocate(itemp2(0:2,0:dims(2)-1))
        do i=0,nelem_in_proc(proc)-1
            nel = which_elem_in_proc(proc,i)
            itemp2(0,i) = Material(nel)
            if (elem_solid(nel)) then
                itemp2(1,i) = 1
            else
                itemp2(1,i) = 0
            end if
            if (all_fluid .OR. solid_fluid) then
                if (elem_fluid_dirich(nel)) then
                    itemp2(2,i) = 1
                else
                    itemp2(2,i) = 0
                end if
            else
                itemp2(2,i) = 0
            end if
        end do
        call write_dataset(fid, "material", itemp2, hdferr)
        deallocate(itemp2)

        ! Global node indexes per element
        call write_dataset(fid, "elements", Ipointer_local, hdferr)

        !Faces
        call write_dataset(fid, "faces", transpose(faces), hdferr)
        call write_dataset(fid, "faces_map", transpose(mapping_faces), hdferr)

        !Edges
        call write_dataset(fid, "edges", transpose(edges), hdferr)
        call write_dataset(fid, "edges_map", transpose(mapping_edges), hdferr)

        ! Vertices
        call write_dataset(fid, "vertices", transpose(vertices), hdferr)
        call write_dataset(fid, "vertices_to_global", vertex_to_glob(0:n_vertices-1), hdferr)

        if (solid_fluid .and. SF%present_local) then
            call write_dataset(fid, "sf_face_near_edges", transpose(SF%Face_Near_Edges), hdferr)
            call write_dataset(fid, "sf_face_near_edges_orient", transpose(SF%Face_Near_Edges_Orient), hdferr)
            call write_dataset(fid, "sf_face_near_vertices", transpose(SF%Face_Near_Vertices), hdferr)
            call write_dataset(fid, "sf_face_glob_interface", transpose(SF%faces), hdferr)
            call write_dataset(fid, "sf_face_orient", SF%face_orient, hdferr)
            call write_dataset(fid, "sf_edge_glob_interface", transpose(SF%edges), hdferr)
            call write_dataset(fid, "sf_edge_orient", SF%mapping_edges, hdferr)
            call write_dataset(fid, "sf_vertex_glob_interface", transpose(SF%vertices), hdferr)
        end if

        if (Neumann_present .and. Neu%present_local) then
            call write_dataset(fid, "neu_face_near_edges", transpose(Neu%Face_Near_Edges), hdferr)
            call write_dataset(fid, "neu_face_near_edges_orient", transpose(Neu%Face_Near_Edges_Orient), hdferr)
            call write_dataset(fid, "neu_face_near_vertices", transpose(Neu%Face_Near_Vertices), hdferr)
            call write_dataset(fid, "neu_face_glob_interface", Neu%faces, hdferr)
            call write_dataset(fid, "neu_edge_glob_interface", Neu%edges, hdferr)
            call write_dataset(fid, "neu_vertex_glob_interface", Neu%vertices, hdferr)

        end if

        ! Communications
        do n=0,nproc-1
            write(proc_grp,"(a,I4.4)") "Proc", n
            call h5gcreate_f(fid, trim(adjustl(proc_grp)), proc_id, hdferr, 0_SIZE_T)

            call write_attr_int(proc_id, "n_faces", shared%nf(n))
            call write_attr_int(proc_id, "n_edges", shared%ne(n))
            call write_attr_int(proc_id, "n_vertices", shared%nv(n))
            if (SF%present_local) then
                call write_attr_int(proc_id, "n_sf_faces", SF%nf_shared(n))
                call write_attr_int(proc_id, "n_sf_edges", SF%ne_shared(n))
                call write_attr_int(proc_id, "n_sf_vertices", SF%nv_shared(n))
            else
                call write_attr_int(proc_id, "n_sf_faces", 0)
                call write_attr_int(proc_id, "n_sf_edges", 0)
                call write_attr_int(proc_id, "n_sf_vertices", 0)
            endif
            if (Neu%present_local) then
                call write_attr_int(proc_id, "n_neu_edges", Neu%ne_shared(n))
                call write_attr_int(proc_id, "n_neu_vertices", Neu%nv_shared(n))
            else
                call write_attr_int(proc_id, "n_neu_edges", 0)
                call write_attr_int(proc_id, "n_neu_vertices", 0)
            endif


            call write_dataset(proc_id, "faces", shared%faces(n,0:shared%nf(n)-1), hdferr)
            call write_dataset(proc_id, "faces_map", shared%mapping_faces(n,0:shared%nf(n)-1), hdferr)
            call write_dataset(proc_id, "edges", shared%edges(n,0:shared%ne(n)-1), hdferr)
            call write_dataset(proc_id, "edges_map", shared%mapping_edges(n,0:shared%ne(n)-1), hdferr)
            call write_dataset(proc_id, "vertices", shared%vertices(n,0:shared%nv(n)-1), hdferr)

            if (SF%present_local) then
                call write_dataset(proc_id, "sf_faces", SF%faces_shared(n,0:SF%nf_shared(n)-1), hdferr)
                call write_dataset(proc_id, "sf_edges", SF%edges_shared(n,0:SF%ne_shared(n)-1), hdferr)
                call write_dataset(proc_id, "sf_edges_map", SF%mapping_edges_shared(n,0:SF%ne_shared(n)-1), hdferr)
                call write_dataset(proc_id, "sf_vertices", SF%vertices_shared(n,0:SF%nv_shared(n)-1), hdferr)
            end if
            if(Neu%present_local)then
                call write_dataset(proc_id, "neu_edges", Neu%edges_shared(n,0:Neu%ne_shared(n)-1), hdferr)
                call write_dataset(proc_id, "neu_edges_map", Neu%mapping_edges_shared(n,0:Neu%ne_shared(n)-1), hdferr)
                call write_dataset(proc_id, "neu_vertices", Neu%vertices_shared(n,0:Neu%nv_shared(n)-1), hdferr)
            end if
            call h5gclose_f(proc_id, hdferr)
        end do


        call h5fclose_f(fid, hdferr)

    end subroutine write_mesh_file_h5

end module mesh2spec
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
