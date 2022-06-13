!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!
!  This module handles orientation of elements
!
!    Z
!    ^
!    |
!       7------6
!      /|     /|   Y
!     / |    / |  /
!    4------5  | /
!    |  3---|--2
!    | /    | /
!    |/     |/
!    0------1   --> X(xi)
!
!
module orientation
    use selement
    use sfaces
    use sedges
    use svertices
    use sdomain
    use mindex
    implicit none


contains

    !! Apply material properties set on elements to its adjacent faces, edges, vertices
    subroutine apply_mat_to_faces(Tdomain)
        type(domain), intent(inout) :: Tdomain
        !
        integer :: i, j, k
        integer :: ngll
        integer :: nf, mat, dom
        integer :: node, dir
        integer, dimension(0:3) :: elface

        do i = 0,Tdomain%n_face-1
            Tdomain%sFace(i)%orphan = .true.
            Tdomain%sFace(i)%elem_0 = -1
            Tdomain%sFace(i)%elem_1 = -1
        end do
        do i = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(i)%mat_index
            dom = Tdomain%sSubDomain(mat)%dom
            ngll = domain_ngll(Tdomain, Tdomain%specel(i)%domain)

            do j = 0,5
                nf = Tdomain%specel(i)%Near_Faces(j)
                Tdomain%sFace(nf)%orphan = .false.
                if (Tdomain%sFace(nf)%domain /= dom) then
                    write(*,*) "Error: inconsistency detected in apply_mat_to_faces"
                    stop 1
                end if
                do k=0,3
                    elface(k) = Tdomain%specel(i)%Control_nodes(face_def(k,j))
                end do

                call rel_orient(Tdomain%sFace(nf)%inodes, elface, node, dir)
                if (node==-1) then
                    write(*,*) "Error: faces orientation inconsistency in apply_mat_to_faces"
                    stop 1
                end if
                if (dir==1) then
                    Tdomain%sFace(nf)%elem_0 = i
                else
                    Tdomain%sFace(nf)%elem_1 = i
                endif
            end do
        end do
    end subroutine apply_mat_to_faces

    subroutine apply_mat_to_edges(Tdomain)
        type(domain), intent(inout) :: Tdomain
        !
        integer :: i, j, k
        integer :: ngll
        integer :: ne, mat, dom
        integer, dimension(0:1) :: eledge

        do i = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(i)%mat_index
            dom = Tdomain%sSubDomain(mat)%dom
            ngll = domain_ngll(Tdomain, Tdomain%specel(i)%domain)

            do j = 0,11
                ne = Tdomain%specel(i)%Near_Edges(j)
                do k=0,1
                    eledge(k) = Tdomain%specel(i)%Control_nodes(edge_def(k,j))
                end do
                if (Tdomain%sEdge(ne)%domain /= dom) then
                    write(*,*) "Error: inconsistency detected in apply_mat_to_edges"
                    stop 1
                end if
            end do
        end do

    end subroutine apply_mat_to_edges

    subroutine apply_mat_to_vertices(Tdomain)
        type(domain), intent(inout) :: Tdomain
        !
        integer :: i, j
        integer :: nv, mat, dom

        do i = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(i)%mat_index
            dom = Tdomain%sSubDomain(mat)%dom
            do j = 0,7
                nv = Tdomain%specel(i)%Near_Vertices(j)
                if (Tdomain%sVertex(nv)%domain /= dom) then
                    write(*,*) "Error: inconsistency detected in apply_mat_to_vertices"
                    stop 1
                end if
            end do
        end do
    end subroutine apply_mat_to_vertices

    ! make sure face,edge,vertices domains are set
    ! and that ngll on both sides match
    !
    ! single surface element without a 3d hex associated can happen
    ! at processor boundary. It's the role of this function to make
    ! sure those are correctly initialised
    subroutine apply_interface(Tdomain, inter, d0, d1, check)
        type(domain), intent(inout) :: Tdomain
        type(inter_num), intent(in) :: inter
        integer, intent(in) :: d0, d1
        logical, intent(in) :: check
        !
        integer :: k
        integer :: i0, i1
        integer :: interface_ok

        interface_ok = 0
        do k=0,inter%surf0%n_faces-1
            i0 = inter%surf0%if_faces(k)
            i1 = inter%surf1%if_faces(k)
            if (.not.check) cycle
            if ((Tdomain%sFace(i0)%domain /= d0).or.(Tdomain%sFace(i1)%domain /= d1)) then
                write(*,*) "Inconsistency detected, unhandled interface (face)"
                write(*,*) "Face0:", i0, "domain=", Tdomain%sFace(i0)%domain, "expected:", d0
                write(*,*) "Face1:", i1, "domain=", Tdomain%sFace(i1)%domain, "expected:", d1
                stop 1
            end if
        end do
        do k=0,inter%surf0%n_edges-1
            i0 = inter%surf0%if_edges(k)
            i1 = inter%surf1%if_edges(k)
            if (.not.check) cycle
            if ((Tdomain%sEdge(i0)%domain /= d0).or.(Tdomain%sEdge(i1)%domain /= d1)) then
                write(*,*) "Inconsistency detected, unhandled interface (edge)"
                write(*,*) "Edge0:", i0, "domain=", Tdomain%sEdge(i0)%domain, "expected:", d0
                write(*,*) "Edge1:", i1, "domain=", Tdomain%sEdge(i1)%domain, "expected:", d1
                stop 1
            end if
        end do
        do k=0,inter%surf0%n_vertices-1
            i0 = inter%surf0%if_vertices(k)
            i1 = inter%surf1%if_vertices(k)
            if (.not.check) cycle
            if ((Tdomain%sVertex(i0)%domain /= d0).or.(Tdomain%sVertex(i1)%domain /= d1)) then
                write(*,*) "Inconsistency detected, unhandled interface (vertex)"
                write(*,*) "Vertex0:", i0, "domain=", Tdomain%sVertex(i0)%domain, "expected:", d0
                write(*,*) "Vertex1:", i1, "domain=", Tdomain%sVertex(i1)%domain, "expected:", d1
                stop 1
            end if
        end do
    end subroutine apply_interface
end module orientation

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
