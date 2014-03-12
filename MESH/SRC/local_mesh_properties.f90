module local_mesh_properties

    use sets
    use mesh_properties

    implicit none

    !- types associated to memory of faces edges and vertices,
    !   already seen in the process of constructing local entities
    type :: processor
       integer, dimension(:), pointer :: E
    end type processor

    type :: souvenir
       type(processor), dimension(:), pointer :: rank
    end type souvenir


contains


    !------------------------------------------------------------
    !------------------------------------------------------------
    subroutine local_nodes_definition(proc,n_nods,n_points,nelem_in_proc,  &
        which_elem_in_proc,Ipointer,n_points_local,           &
        which_points_inside,glob2loc,Ipointer_local,node_loc2glob)

        implicit none
        integer, intent(in)    :: proc,n_nods,n_points
        integer, dimension(0:), intent(in) :: nelem_in_proc
        integer, dimension(0:,0:), intent(in)   :: which_elem_in_proc,Ipointer
        integer, intent(out)     :: n_points_local
        logical, dimension(0:n_points-1), intent(out)  :: which_points_inside
        integer, dimension(0:n_points-1), intent(out)  :: glob2loc
        integer, dimension(0:,0:), intent(out)    :: Ipointer_local
        integer, dimension(:), allocatable, intent(out)  :: node_loc2glob

        integer     :: n,nel,i,j,num,idummy

        write(*,*) "    --> Nodes on proc. ",proc
        which_points_inside(0:) = .false.
        n_points_local = 0 ; glob2loc(0:) = -1
        do n = 0,nelem_in_proc(proc)-1
            nel = which_elem_in_proc(proc,n)
            do i = 0,n_nods-1
                num = Ipointer(i,nel)
                which_points_inside(num) = .true.
            end do
        end do
        !- correspondence global <-> local nodes
        !-
        do n = 0, n_points-1
            if(which_points_inside(n))then
                glob2loc(n) = n_points_local
                n_points_local = n_points_local+1
            end if
        end do
        !-
        allocate(node_loc2glob(0:n_points_local-1))
        j = 0
        do n = 0, n_points-1
            if(which_points_inside(n))then
                node_loc2glob(j) = n
                j = j+1
            end if
        end do

        !- now we construct the set of local nodes
        do n = 0, nelem_in_proc(proc)-1
            nel = which_elem_in_proc(proc,n)
            do i = 0, n_nods-1
                idummy = Ipointer(i,nel)
                Ipointer_local(i,n) = glob2loc(idummy)
            end do
        end do

    end subroutine local_nodes_definition
    !------------------------------------------------------------
    !------------------------------------------------------------
    subroutine local_vertex_construct(proc,n_points_local,nelem_in_proc,        &
        Ipointer_local,n_vertices,vertices,vertex_to_glob,N_valid_Vertex)

        implicit none

        integer, intent(in)     :: proc, n_points_local
        integer, dimension(0:), intent(in)   :: nelem_in_proc
        integer, dimension(0:,0:), intent(in)   :: Ipointer_local
        integer, intent(out)     :: n_vertices
        integer, dimension(0:,0:), intent(out)   :: vertices
        integer, dimension(:), allocatable, intent(out)  :: vertex_to_glob,   &
            N_valid_Vertex

        integer   :: nMaxVertex, n,i,j,idummy


        write(*,*) "    --> Construction of local vertices."
        nMaxVertex = 8*nelem_in_proc(proc)
        print*, "      --> Max.number of vertices: ", nMaxVertex
        allocate(vertex_to_glob(0:nMaxVertex-1))
        allocate(N_valid_Vertex(0:n_points_local-1))
        N_valid_vertex(0:) = -1
        n_vertices = 0
        do n = 0,nelem_in_proc(proc)-1
            do j = 0,7
                idummy = Ipointer_local(j,n)
                if(N_valid_vertex(idummy) < 0)then   ! vertex never seen
                    N_valid_Vertex(idummy) = n_vertices
                    vertices(n,j) = n_vertices
                    vertex_to_glob(n_vertices) = idummy
                    n_vertices = n_vertices+1
                else
                    vertices(n,j) = N_valid_Vertex(idummy)
                end if
            end do
        end do


    end subroutine local_vertex_construct
    !------------------------------------------------------------
    !------------------------------------------------------------
    subroutine local_faces_construct(proc,part,nelem_in_proc,which_elem_in_proc,Elem_glob2loc,   &
        Ipointer,dxadj,dxadjncy,elem_near_proc,n_faces,faces,mapping_faces,     &
        faces_shared,mapping_faces_shared,nf_shared,memory)

        implicit none
        integer, intent(in)   :: proc
        integer, dimension(0:), intent(in)   :: part,nelem_in_proc,dxadj,dxadjncy,Elem_glob2loc
        integer, dimension(0:,0:), intent(in) :: which_elem_in_proc,Ipointer
        type(near_proc), dimension(0:), intent(in) :: elem_near_proc 
        integer, intent(out)    :: n_faces
        integer, dimension(0:,0:), intent(out)  :: faces,mapping_faces
        integer, dimension(0:,0:), intent(out)  :: faces_shared,mapping_faces_shared
        integer, dimension(0:), intent(out)     :: nf_shared
        type(souvenir), dimension(0:), intent(inout)  :: memory
        integer    :: i,j,k,n,nel,nf,neighbor,num,ok,neighbor_face,i_count,procrank
        integer, dimension(0:3)  ::  corner,neighbor_corner,orient_corner

        faces = -1
        mapping_faces = -1
        mapping_faces_shared = -1
        nf_shared = 0
        n_faces = 0
        do n = 0,nelem_in_proc(proc)-1
            nel = which_elem_in_proc(proc,n)
            do nf = 0,5
                call face2corner(Ipointer(0:,nel),nf,corner)

                find0 : do i = dxadj(nel), dxadj(nel+1)-1
                    ! we look at a neighbor of the element
                    neighbor = dxadjncy(i)
                    find1 : do j = 0,3   ! does the face belong to this neighbour?
                        num = corner(j)
                        ok = 0
                        find2 : do k = 0,7
                            if(Ipointer(k,neighbor) == num)then
                                neighbor_corner(j) = k
                                ok = 1
                                exit find2
                            endif
                        enddo find2
                        if(ok == 0) exit find1   ! No, so let's see another neighbor
                        if(j == 3) then   ! Yes
                            if(part(neighbor) == proc) then  ! The neighbor is on the same processor
                                if(neighbor > nel) then   ! The neighbor is an element we've never seen
                                    faces(n,nf) = n_faces
                                    mapping_faces(n,nf) = 0
                                    n_faces = n_faces + 1
                                else   ! The neighbor is an element already seen
                                    orient_corner(0:3) = neighbor_corner(0:3)
                                    call sort(neighbor_corner,4)
                                    neighbor_face = neighb_face(neighbor_corner)
                               !     g2l : do i_count = 0,n-1
                               !         if(which_elem_in_proc(proc,i_count) == neighbor) then
                                            i_count = Elem_glob2loc(neighbor)
                                            faces(n,nf) = faces(i_count,neighbor_face)
                               !             exit g2l
                               !         endif
                               !     enddo g2l
                                    !! Coherency
                                    call face_orientation(orient_corner,4,neighbor_face,mapping_faces(n,nf))
                                endif
                            else   ! The neighbor is not on the same processor than the element
                                faces(n,nf) = n_faces
                                mapping_faces(n,nf) = 0
                                num = part(neighbor)
                                ! Ensuring the correspondence between the faces shared
                                if(num < proc) then   ! we've already seen the processor of the neighbor
                                    !! Coherency
                                    orient_corner(0:3) = neighbor_corner(0:3)
                                    call sort(neighbor_corner,4)
                                    neighbor_face = neighb_face(neighbor_corner)
                                    procrank = already_in(proc,elem_near_proc(neighbor)%list,elem_near_proc(neighbor)%nb)
                                    if(procrank < 0) stop "Gros probleme avec les procs"
                                    faces_shared(num,memory(neighbor)%rank(procrank)%E(neighbor_face))= n_faces
                                    call face_orientation(orient_corner,4,neighbor_face,    &
                                        mapping_faces_shared(num,memory(neighbor)%rank(procrank)%E(neighbor_face)))
                                else   ! We've never seen the processor of the neighbor
                                    procrank = already_in(num,elem_near_proc(nel)%list,elem_near_proc(nel)%nb)
                                    if(procrank < 0) stop "Gros probleme avec les procs"
                                    faces_shared(num,nf_shared(num)) = n_faces
                                    memory(nel)%rank(procrank)%E(nf) = nf_shared(num)
                                    !! Coherency
                                    orient_corner(0:3) = neighbor_corner(0:3)
                                    call sort(neighbor_corner,4)
                                    neighbor_face = neighb_face(neighbor_corner)
                                    call face_orientation(orient_corner,4,neighbor_face,&
                                        mapping_faces_shared(num,nf_shared(num)))
                                endif
                                nf_shared(num) = nf_shared(num) + 1
                                n_faces = n_faces + 1
                            endif
                            exit find0
                        endif
                    end do find1
                end do find0
                if(ok == 0) then   ! The face is not shared by a neighbor
                    faces(n,nf) = n_faces
                    mapping_faces(n,nf) = 0
                    n_faces = n_faces + 1
                endif
            enddo
        enddo

    end subroutine local_faces_construct
    !------------------------------------------------------------
    !------------------------------------------------------------
    subroutine local_edges_construct(nproc,proc,part,nelem_in_proc,Elem_glob2loc, &
        which_elem_in_proc,Ipointer,near_elem_set,elem_near_proc,n_edges,         &
        edges,mapping_edges,edges_shared,mapping_edges_shared,       &
        ne_shared,Elem_edge_ref,memory)

        implicit none
        integer, intent(in)  :: nproc,proc
        integer, dimension(0:), intent(in) :: part,nelem_in_proc,Elem_glob2loc
        integer, dimension(0:,0:), intent(in)    :: which_elem_in_proc,Ipointer
        type(near_elem), dimension(0:), intent(in)   :: near_elem_set
        type(near_proc), dimension(0:), intent(in)   :: elem_near_proc
        integer, intent(out)  :: n_edges
        integer, dimension(0:,0:), intent(out) :: edges,mapping_edges
        integer, dimension(0:,0:), intent(out) :: edges_shared,mapping_edges_shared
        integer, dimension(0:), intent(out)  :: ne_shared,Elem_Edge_Ref
        type(souvenir), dimension(0:), intent(inout)  :: memory
        integer   :: i,j,k,n,ne,nel,ok,num,neighbor,neighbor_edge,procrank
        integer, dimension(0:1)  :: e_corner,e_neighbor_corner,e_orient_corner
        logical, dimension(0:nproc-1)  :: L_Proc
        type(near_entity), pointer  :: near_neighb => NULL(),   &
            o_proc_near_neighb => NULL()

        mapping_edges = -1
        mapping_edges_shared = -1
        ne_shared = 0
        n_edges = 0
        do n = 0,nelem_in_proc(proc)-1
            nel = which_elem_in_proc(proc,n)
            do ne = 0,11
                call edge2corner(Ipointer(0:,nel),ne,e_corner)
                near_neighb => near_elem_set(nel)%ptr  ! pointer on the linked list of nearest elements
                findbis0 : do while(associated(near_neighb)) ! loop on nearest neighbors only: already seen on this proc?
                    ok = 0
                    neighbor = near_neighb%elem
                    if(part(neighbor) /= proc .or. Elem_glob2loc(neighbor) >= n)then
                        near_neighb => near_neighb%pt
                        cycle findbis0
                    end if
                    findbis1 : do j = 0,1
                        num = e_corner(j)
                        ok = 0
                        findbis2 : do k = 0,7
                            if (Ipointer(k,neighbor) == num) then
                                e_neighbor_corner(j) = k
                                ok = 1
                                exit findbis2
                            endif
                        enddo findbis2
                        if (ok == 0) exit findbis1
                        if (j == 1) then   ! YES
                            e_orient_corner(0:1) = e_neighbor_corner(0:1)
                            call sort(e_neighbor_corner,2)
                            neighbor_edge = neighb_edge(e_neighbor_corner)
                            edges(n,ne) = edges(Elem_glob2loc(neighbor),neighbor_edge)
                            !! Coherency
                            call edge_orientation(e_orient_corner,2,neighbor_edge,mapping_edges(n,ne))
                            exit findbis0
                        endif
                    enddo findbis1
                    near_neighb => near_neighb%pt   ! next element in the linked list
                enddo findbis0
                if(ok == 0 .or. n == 0) then   ! no: edge never seen on the same proc before
                    edges(n,ne) = n_edges
                    mapping_edges(n,ne) = 0
                    Elem_Edge_Ref(n_edges) = n
                    n_edges = n_edges + 1

                    ! now look to other procs
                    L_Proc = .true.
                    L_Proc(proc) = .false.
                    o_proc_near_neighb => near_elem_set(nel)%ptr
                    do while(associated(o_proc_near_neighb))
                        i = o_proc_near_neighb%elem ! edge shared with another on another proc?
                        if (L_Proc(part(i))) then
                            findter1 : do j = 0,1
                                num = e_corner(j)
                                ok = 0
                                findter2 : do k = 0,7
                                    if (Ipointer(k,i) == num) then
                                        e_neighbor_corner(j) = k
                                        ok = 1
                                        exit findter2
                                    endif
                                enddo findter2
                                if (ok == 0) exit findter1
                                if (j == 1) then   ! YES
                                    num = part(i)
!!! Ensuring the correspondence between the edges shared !!!
                                    if (num<proc) then   ! It deals with a processor we've ever seen
                                        e_orient_corner(0:1) = e_neighbor_corner(0:1)
                                        call sort(e_neighbor_corner,2)
                                        neighbor_edge = neighb_edge(e_neighbor_corner)
                                        procrank = already_in(proc,elem_near_proc(i)%list,elem_near_proc(i)%nb)
                                        if(procrank < 0) stop "Gros probleme avec les procs"
                                        edges_shared(num,memory(i)%rank(procrank)%E(neighbor_edge+6)) = edges(n,ne)
                                        !! Coherency
                                        call edge_orientation(e_orient_corner,2,neighbor_edge,&
                                            mapping_edges_shared(num,memory(i)%rank(procrank)%E(neighbor_edge+6)))
                                    else   ! It deals with a processor we've never seen
                                        procrank = already_in(num,elem_near_proc(nel)%list,elem_near_proc(nel)%nb)
                                        if(procrank < 0) stop "Gros probleme avec les procs"
                                        edges_shared(num,ne_shared(num)) = edges(n,ne)
                                        memory(nel)%rank(procrank)%E(ne+6) = ne_shared(num)
                                        ! Coherency
                                        e_orient_corner(0:1) = e_neighbor_corner(0:1)
                                        call sort(e_neighbor_corner,2)
                                        neighbor_edge = neighb_edge(e_neighbor_corner)

                                        call edge_orientation(e_orient_corner,2,neighbor_edge,&
                                            mapping_edges_shared(num,ne_shared(num)))
                                    endif
                                    ne_shared(num) = ne_shared(num) + 1
                                    L_Proc(num) = .false.
                                endif
                            enddo findter1
                        endif
                        o_proc_near_neighb => o_proc_near_neighb%pt  ! next elem in the linked list
                    enddo
                endif
            enddo
        enddo

    end subroutine local_edges_construct
    !------------------------------------------------------------
    !------------------------------------------------------------
    subroutine local_vertices_comm(n_vertices,proc,n_proc,nelem_in_proc,part,  &
        vertex_to_glob,node_loc2glob,which_elem_in_proc,Ipointer,    &
        vertices,near_elem_set,elem_near_proc,vertices_shared,nv_shared,memory)

        implicit none
        integer, intent(in)  :: n_vertices,proc,n_proc
        integer, dimension(0:), intent(in)  :: nelem_in_proc,part,    &
            vertex_to_glob,node_loc2glob
        integer, dimension(0:,0:), intent(in)  :: which_elem_in_proc, Ipointer,  &
            vertices
        type(near_elem), dimension(0:), intent(in)   :: near_elem_set
        type(near_proc), dimension(0:), intent(in)   :: elem_near_proc
        integer, dimension(0:,0:), intent(out) :: vertices_shared
        integer, dimension(0:), intent(out)   :: nv_shared
        type(souvenir), dimension(0:), intent(inout)  :: memory
        integer   :: n,nv,vert,local_n,icount,nel,i,k,num,procrank
        logical, dimension(0:n_vertices-1)  :: L_vertex
        logical, dimension(0:n_proc-1)  :: L_Proc
        type(near_entity), pointer  :: near_neighb => NULL(),   &
            o_proc_near_neighb => NULL()


        nv_shared = 0
        L_vertex(0:) = .false.

        do n = 0,nelem_in_proc(proc)-1
            vert_loop: do nv = 0,7
                vert = vertices(n,nv)
                if(L_vertex(vert)) cycle vert_loop
                L_vertex(vert) = .true.
                local_n = vertex_to_glob(vert)
                icount = node_loc2glob(local_n)
                nel = which_elem_in_proc(proc,n)
                L_proc(0:) = .true. ; L_proc(proc) = .false.
                o_proc_near_neighb => near_elem_set(nel)%ptr
                do while(associated(o_proc_near_neighb))
                    i = o_proc_near_neighb%elem
                    if(L_proc(part(i)))then
                        do k = 0,7
                            if(Ipointer(k,i) == icount)then
                                num = part(i)
                                if(num < proc)then   ! it deals with a processor we've already seen
                                    procrank = already_in(proc,elem_near_proc(i)%list,elem_near_proc(i)%nb)
                                    if(procrank < 0) stop "Gros probleme avec les procs"
                                    vertices_shared(num,memory(i)%rank(procrank)%E(18+k)) = vert
                                else   ! it deals with a processor we've never seen
                                    procrank = already_in(num,elem_near_proc(nel)%list,elem_near_proc(nel)%nb)
                                    if(procrank < 0) stop "Gros probleme avec les procs"
                                    vertices_shared(num,nv_shared(num)) = vert
                                    memory(nel)%rank(procrank)%E(18+nv) = nv_shared(num)
                                end if
                                nv_shared(num) = nv_shared(num)+1
                                L_proc(num) = .false.
                            end if
                        end do
                    end if
                    o_proc_near_neighb => o_proc_near_neighb%pt
                end do
            end do vert_loop
        end do


    end subroutine local_vertices_comm
    !------------------------------------------------------------
    !------------------------------------------------------------
end module local_mesh_properties

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
