!>
!!\file Vertex.f90
!!\brief Assure la gestion des Vertex.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module svertices

    ! Modified by Gaetano 31/01/2005
    ! Modified by Paul 06/11/2005

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
#ifdef COUPLAGE
       real, dimension (:), allocatable :: ForcesMka
       real :: tsurfsem
#endif

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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
