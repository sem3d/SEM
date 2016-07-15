!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Vertex.f90
!!\brief Assure la gestion des Vertex.
!!
!<

module svertices

    implicit none
    type :: vertex
       integer :: domain
       integer :: Iglobnum_Vertex
       integer :: Idom
       integer :: inode
       real, dimension (:), allocatable :: Forces, Forces3
       logical                          :: PML
    end type vertex

contains
    subroutine init_vertex(ve)
        type(Vertex), intent(inout) :: ve
        !
        ve%domain = -1
        ve%Iglobnum_Vertex = -1
    end subroutine init_vertex


   subroutine allocate_vertex_force(ve)
    
        implicit none
        type(vertex), intent(inout) :: ve
    
        allocate(ve%Forces(0:2))
        allocate(ve%Forces3(0:2))
        ve%Forces = 0.0
        ve%Forces3= 0.0

   end subroutine allocate_vertex_force

   subroutine free_vertex_force(ve)
      implicit none
      type(vertex), intent(inout) :: ve
         
       deallocate(ve%Forces)
       deallocate(ve%Forces3)
   end subroutine free_vertex_force

end module svertices

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
