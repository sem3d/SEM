!>
!! \file Edge.f90
!! \brief
!!
!<

module sedges

    type :: edge_pml
       real, dimension (:,:), allocatable :: IVeloc1, Iveloc2, Iveloc3
       real, dimension (:), allocatable :: Ivx, Ivy, Ivz

    end type edge_pml
    type :: edge

       logical :: PML, Abs, FPML

       integer :: ngll,mat_index
       integer, dimension (:), allocatable :: Iglobnum_Edge

       ! Lien entre ngll et numérotation des champs globaux
       integer, dimension (:), allocatable :: Renum

       ! solid-fluid
       logical  :: solid, fluid_dirich

       type(edge_pml), pointer :: spml
#ifdef COUPLAGE
       real, dimension (:,:), allocatable :: ForcesMka
       !     integer, dimension (:,:), allocatable :: FlagMka
       real, dimension (:), allocatable :: tsurfsem
#endif


    end type edge

contains

    ! ###########################################################
    subroutine init_edge(ed)
        type(Edge), intent(inout) :: ed

        ed%PML = .false.
        ed%Abs = .false.
        ed%FPML = .false.
        ed%ngll = 0
        ed%solid = .true.
    end subroutine init_edge

end module sedges
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
