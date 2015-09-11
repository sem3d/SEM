!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Face.f90
!!\brief Gère les faces des éléments.
!!
!<

module sfaces


    type :: face
       logical :: PML, Abs, FPML
       integer :: ngll1, ngll2, dir, Which_Elem, mat_index
       integer, dimension (:,:), allocatable :: Iglobnum_Face
       ! Lien entre ngll et numérotation des champs globaux
       integer, dimension (:,:), allocatable :: Renum

       ! solid-fluid
       logical :: solid, fluid_dirich

       !! Couplage Externe
!       real, dimension (:,:,:), allocatable :: ForcesExt
!       real, dimension (:,:), allocatable :: tsurfsem
    end type face

contains
    ! ############################################################

    subroutine init_face(fc)
        type(Face), intent(inout) :: fc

        fc%PML = .false.
        fc%Abs = .false.
        fc%FPML = .false.
        fc%ngll1 = 0
        fc%ngll2 = 0
        fc%dir = -1
        fc%Which_Elem = -1
        fc%solid = .true.
    end subroutine init_face

end module sfaces

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
