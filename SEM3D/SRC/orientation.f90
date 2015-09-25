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
        integer, dimension(0:2) :: ngll
        integer :: nf, mat, dom
        integer :: node, dir
        integer, dimension(0:3) :: elface

        do i = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(i)%mat_index
            dom = get_domain(Tdomain%sSubDomain(mat))
            ngll(0) = Tdomain%specel(i)%ngllx
            ngll(1) = Tdomain%specel(i)%nglly
            ngll(2) = Tdomain%specel(i)%ngllz

            do j = 0,5
                nf = Tdomain%specel(i)%Near_Faces(j)
                do k=0,3
                    !write(*,*) 'k=', k, 'f(k,j)=', face_def(k,j), size(Tdomain%specel(i)%Control_nodes)
                    elface(k) = Tdomain%specel(i)%Control_nodes(face_def(k,j))
                end do
                call rel_orient(Tdomain%sFace(nf)%inodes, elface, node, dir)
                if ((node==0) .or. (node==2)) then
                    Tdomain%sFace(nf)%ngll1 = ngll(face_dir(0,j))
                    Tdomain%sFace(nf)%ngll2 = ngll(face_dir(1,j))
                else
                    Tdomain%sFace(nf)%ngll1 = ngll(face_dir(1,j))
                    Tdomain%sFace(nf)%ngll2 = ngll(face_dir(0,j))
                end if
                Tdomain%sFace(nf)%domain = dom
            end do
        end do

    end subroutine apply_mat_to_faces

    subroutine apply_mat_to_edges(Tdomain)
        type(domain), intent(inout) :: Tdomain
        !
        integer :: i, j, k
        integer, dimension(0:2) :: ngll
        integer :: ne, mat, dom
        integer, dimension(0:1) :: eledge

        do i = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(i)%mat_index
            dom = get_domain(Tdomain%sSubDomain(mat))
            ngll(0) = Tdomain%specel(i)%ngllx
            ngll(1) = Tdomain%specel(i)%nglly
            ngll(2) = Tdomain%specel(i)%ngllz

            do j = 0,11
                ne = Tdomain%specel(i)%Near_Edges(j)
                do k=0,1
                    eledge(k) = Tdomain%specel(i)%Control_nodes(edge_def(k,j))
                end do
                Tdomain%sEdge(ne)%ngll = ngll(edge_axis(j))
                Tdomain%sEdge(ne)%domain = dom
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
            dom = get_domain(Tdomain%sSubDomain(mat))
            do j = 0,7
                nv = Tdomain%specel(i)%Near_Vertices(j)
                Tdomain%sVertex(nv)%domain = dom
            end do
        end do
    end subroutine apply_mat_to_vertices
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
