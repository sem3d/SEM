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

    type :: vertex_pml
       real, dimension (:), allocatable :: Iveloc1, Iveloc2, Iveloc3
       real, dimension (:), allocatable :: Ivx, Ivy, Ivz
    end type vertex_pml

    type :: vertex
       integer  :: mat_index
       logical :: PML, Abs, FPML
       integer :: Iglobnum_Vertex, global_numbering
       ! Lien entre ngll et numérotation des champs globaux
       integer :: Renum

       ! solid-fluid
       logical :: solid, fluid_dirich
       type(vertex_pml), pointer :: spml

       !! Couplage Externe
       real, dimension (:), allocatable :: ForcesExt
       real :: tsurfsem

    end type vertex

contains
    subroutine init_vertex(ve)
        type(Vertex), intent(inout) :: ve

        ve%PML = .false.
        ve%Abs = .false.
        ve%FPML = .false.
        ve%solid = .true.
        ve%global_numbering = -1
        ve%Iglobnum_Vertex = -1
    end subroutine init_vertex

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
