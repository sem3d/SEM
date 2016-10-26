!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file indexation.f90
!!\brief calcul de la numerotation des ddl
!!
!<
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

module mindex
    implicit none

    integer, parameter, dimension(0:2,0:7) :: vertex_def = reshape( (/ &
        0,0,0, &
        1,0,0, &
        1,1,0, &
        0,1,0, &
        0,0,1, &
        1,0,1, &
        1,1,1, &
        0,1,1 /), (/ 3, 8 /) )

    ! This needs to match RefFace in meshpart.cpp
    integer, parameter, dimension(0:3,0:5) :: face_def = reshape( (/ &
        0, 1, 2, 3, &
        0, 4, 5, 1, &
        1, 5, 6, 2, &
        3, 2, 6, 7, &
        0, 3, 7, 4, &
        4, 7, 6, 5 /), (/ 4, 6 /) )

    ! axes of the element corresponding to the face according to the reference def of face_def
    ! ie if face_def(a,b,c,d)  face_dir(0)=axis of elem to go from a to b
    ! face_dir(1)=axis of elem to go from a to d
    integer, parameter, dimension(0:1,0:5) :: face_dir = reshape( (/ &
        0, 1, &  ! x,y
        2, 0, &  ! z,x
        2, 1, &  ! z,y
        0, 2, &  ! x,z
        1, 2, &  ! y,z
        1, 0 /), (/ 2, 6 /))

    integer, parameter, dimension(0:1,0:11) :: edge_def = reshape( (/ &
        0, 1, &  ! 0
        1, 2, &  ! 1
        3, 2, &  ! 2
        0, 3, &  ! 3
        4, 5, &  ! 4
        5, 6, &  ! 5
        7, 6, &  ! 6
        4, 7, &  ! 7
        0, 4, &  ! 8
        1, 5, &  ! 9
        2, 6, &  ! 10
        3, 7 /), (/ 2, 12 /) )

    integer, parameter, dimension(0:11) :: edge_axis = (/ &
        0, 1, 0, 1, &
        0, 1, 0, 1, &
        2, 2, 2, 2 /)

contains

    ! From two face definitions (4 nodes each), return in node
    ! -1 if they are not the same
    ! 0..3 : the position of refface[0] in face
    ! dir : 1 if refface[1] == face[mod(node+1,4)]
    !     : -1 if refface[1] == face[mod(node-1,4)]
    !
    ! XXX: refface and face are supposed to be the same, we don't check all
    ! the nodes for equality
    subroutine rel_orient(refface, face, node, dir)
        integer, dimension(0:3), intent(in) :: refface
        integer, dimension(0:3), intent(in) :: face
        integer, intent(out) :: node
        integer, intent(out) :: dir
        !
        integer :: i
        do i = 0, 3
            if (refface(0)==face(i)) then
                node = i
                if (refface(1)==face(mod(i+1,4))) then
                    dir = 1
                    return
                else if (refface(1)==face(mod(i+3,4))) then
                    dir = -1
                    return
                end if
                node = -1
                dir = -1
                return
            end if
        end do
        node = -1
        dir = -1
    end subroutine rel_orient

    ! Return for an element with ngll(0:2) glls and its face nf, the arrays idxi, idxj, idxk
    ! such that face_ngll(i,j) = elem_ngll(idx0(i,j), idx1(i,j), idx2(i,j) )
    ! with idx0(i,j) = i0(0)+i*di(0)+j*dj(0) ...
    ! refface describes the nodes of the face
    ! nf is the number (0..5) of the face of the element
    ! face describes the nodes of the face nf of the element
    subroutine ind_elem_face(ngll, nf, refface, face, i0, di, dj)
        integer, intent(in) :: ngll
        integer, intent(in) :: nf
        integer, intent(in), dimension(0:3) :: refface, face
        integer, intent(out), dimension(0:2) :: i0, di, dj
        !
        integer :: node, dir, dir1, dir2, dir3
        integer :: ngll1
        call rel_orient(refface, face, node, dir)
        ! dir1..dir3 are the axis numbers (0,1,2) of respectively
        ! i index of face, j index of face, fixed index of face within
        ! the element
        dir1 = face_dir(0, nf)
        dir2 = face_dir(1, nf)
        dir3 = 3-dir1-dir2
        di(dir3) = 0
        dj(dir3) = 0
        ngll1 = ngll-1
        select case (nf)
        case (0,1,4)
            i0(dir3) = 0
            di(dir3) = 0
            dj(dir3) = 0
        case (2,3,5)
            i0(dir3) = ngll1
            di(dir3) = 0
            dj(dir3) = 0
        end select
        if (dir>0) then
            select case (node)
            case (0)
                ! Face    Elem
                ! d c     d c
                ! a b     a b
                i0(dir1) = 0
                i0(dir2) = 0
                di(dir1) = 1
                di(dir2) = 0
                dj(dir1) = 0
                dj(dir2) = 1
            case (1)
                ! Face    Elem
                ! d c     c b
                ! a b     d a
                i0(dir1) = ngll1
                i0(dir2) = 0
                di(dir1) = 0
                di(dir2) = 1
                dj(dir1) = -1
                dj(dir2) = 0
            case (2)
                ! Face    Elem
                ! d c     b a
                ! a b     c d
                i0(dir1) = ngll1
                i0(dir2) = ngll1
                di(dir1) = -1
                di(dir2) = 0
                dj(dir1) = 0
                dj(dir2) = -1
            case (3)
                ! Face    Elem
                ! d c     a d
                ! a b     b c
                i0(dir1) = 0
                i0(dir2) = ngll1
                di(dir1) = 0
                di(dir2) = -1
                dj(dir1) = 1
                dj(dir2) = 0
            end select
        else
            select case (node)
            case (0)
                ! Face    Elem
                ! d c     b c
                ! a b     a d
                i0(dir1) = 0
                i0(dir2) = 0
                di(dir1) = 0
                di(dir2) = 1
                dj(dir1) = 1
                dj(dir2) = 0
            case (1)
                ! Face    Elem
                ! d c     c d
                ! a b     b a
                i0(dir1) = ngll1
                i0(dir2) = 0
                di(dir1) = -1
                di(dir2) = 0
                dj(dir1) = 0
                dj(dir2) = 1
            case (2)
                ! Face    Elem
                ! d c     d a
                ! a b     c b
                i0(dir1) = ngll1
                i0(dir2) = ngll1
                di(dir1) = 0
                di(dir2) = -1
                dj(dir1) = -1
                dj(dir2) = 0
            case (3)
                ! Face    Elem
                ! d c     a b
                ! a b     d c
                i0(dir1) = 0
                i0(dir2) = ngll1
                di(dir1) = 1
                di(dir2) = 0
                dj(dir1) = 0
                dj(dir2) = -1
            end select
        endif
    end subroutine ind_elem_face
    !------------------------------------------------------------
    !------------------------------------------------------------
    subroutine ind_elem_edge(ngll, ne, refedge, edge, i0, di)
        integer, intent(in) :: ngll
        integer, intent(in) :: ne
        integer, intent(in), dimension(0:1) :: refedge, edge
        integer, intent(out), dimension(0:2) :: i0, di
        !
        integer :: axis
        !
        di = 0
        i0 = 0
        axis = edge_axis(ne)
        ! Edges _| to X
        select case(ne)
        case (3,11,7,8)
            i0(0) = 0
        case (1,10,5,9)
            i0(0) = ngll-1
        end select
        ! Edges _| to Y
        select case(ne)
        case (0,9,4,8)
            i0(1) = 0
        case (2,10,6,11)
            i0(1) = ngll-1
        end select
        select case(ne)
        case (4,5,6,7)
            i0(2) = ngll-1
        case (0,1,2,3)
            i0(2) = 0
        end select
        if (refedge(0)==edge(0)) then
            ! direct order
            di(axis) = 1
        else
            i0(axis) = ngll-1
            di(axis) = -1
        end if
    end subroutine ind_elem_edge
end module mindex
!---------------------------------------------------------------------------
!-------------------------------------------------------------------------

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
