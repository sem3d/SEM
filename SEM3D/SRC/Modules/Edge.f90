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
        real, dimension (:,:), allocatable :: Forces, Forces3
        logical                            :: PML
    end type edge

contains

    ! ###########################################################
    subroutine init_edge(ed)
        type(Edge), intent(inout) :: ed
        !
        ed%ngll = 0
        ed%domain = -1
    end subroutine init_edge

    
    subroutine allocate_edge_force(ed)
        
      implicit none
      type(Edge), intent(inout) :: ed
             
      allocate(ed%Forces(1:ed%ngll-2,0:2))
      allocate(ed%Forces3(1:ed%ngll-2,0:2))
      ed%Forces = 0.0
      ed%Forces3= 0.0

    end subroutine allocate_edge_force 

    subroutine free_edge_force(ed)
     
      implicit none
      type(edge), intent(inout) :: ed
     
      deallocate(ed%Forces)
      deallocate(ed%Forces3)
    end subroutine free_edge_force

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
