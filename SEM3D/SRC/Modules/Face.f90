!>
!!\file Face.f90
!!\brief Gère les faces des éléments.
!!
!<

module sfaces

    type :: face_pml
       real, dimension (:,:,:), allocatable :: IVeloc1, IVeloc2, IVeloc3
       real, dimension (:,:), allocatable :: Ivx, Ivy, Ivz
    end type face_pml

    type :: face
       logical :: PML, Abs, FPML
       integer :: ngll1, ngll2, dir, Which_Elem, mat_index
       integer, dimension (:,:), allocatable :: Iglobnum_Face
       ! Lien entre ngll et num�rotation des champs globaux
       integer, dimension (:,:), allocatable :: Renum

       ! solid-fluid
       logical :: solid, fluid_dirich
       ! pml
       type(face_pml), pointer :: spml

       !! Couplage Externe
       real, dimension (:,:,:), allocatable :: ForcesExt
       real, dimension (:,:), allocatable :: tsurfsem
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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
