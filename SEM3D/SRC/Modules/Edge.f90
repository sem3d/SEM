!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file Edge.f90
!! \brief
!!
!<

module sedges
    implicit none
    type :: edge

        integer :: ngll
        integer :: domain
        integer, dimension (:), allocatable :: Iglobnum_Edge
        integer, dimension (:), allocatable :: Idom
        integer, dimension(0:1) :: inodes

        logical :: Abs
    end type edge

contains

    ! ###########################################################
    subroutine init_edge(ed)
        type(Edge), intent(inout) :: ed
        !
        ed%ngll = 0
        ed%domain = -1
        ed%Abs = .false.
    end subroutine init_edge

end module sedges

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
