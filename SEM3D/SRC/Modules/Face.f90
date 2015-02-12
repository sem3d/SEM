!>
!!\file Face.f90
!!\brief Gère les faces des éléments.
!!
!<

module sfaces

    type :: face_pml
       real, dimension (:,:,:), allocatable :: Forces1, Forces2, Forces3, Veloc1, Veloc2, Veloc3
       real, dimension (:,:,:), allocatable :: DumpVx, DumpVy, DumpVz, DumpMass
       real, dimension (:,:,:), allocatable :: IVeloc1, IVeloc2, IVeloc3
       real, dimension (:,:), allocatable :: Ivx, Ivy, Ivz
       real, dimension(:,:), allocatable :: ForcesFl1, ForcesFl2, ForcesFl3, VelPhi1, VelPhi2, VelPhi3
    end type face_pml

    type :: face
       logical :: PML, Abs, FPML
       integer :: ngll1, ngll2, dir, Which_Elem, mat_index
       integer, dimension (:), allocatable :: FaceNum
       integer, dimension (:,:), allocatable :: Iglobnum_Face
       real, dimension (:,:), allocatable  :: MassMat
       real, dimension (:,:,:), allocatable :: Forces, Displ, Veloc, Accel, V0

       !! solid-fluid
       logical :: solid, fluid_dirich
       real, dimension(:,:), allocatable :: ForcesFl, Phi, VelPhi, AccelPhi, VelPhi0

       !! pml
       type(face_pml), pointer :: spml

       !! Couplage Externe
       real, dimension (:,:,:), allocatable :: ForcesExt
       real, dimension (:,:), allocatable :: tsurfsem
    end type face

contains

    ! ############################################################
    !>
    !! \brief
    !!
    !! \param type (Face), intent (INOUT) F
    !<
    subroutine Prediction_Face_Veloc (F, dt)
        implicit none
        type (Face), intent (INOUT) :: F
        real, intent(in) :: dt

        F%Forces(:,:,:) = F%Displ(:,:,:)
        F%V0(:,:,:) = F%Veloc(:,:,:)
        return
    end subroutine Prediction_Face_Veloc

    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Prediction_Face_VelPhi(F,dt)
        implicit none
        type(Face), intent(inout) :: F
        real, intent(in) :: dt

        F%VelPhi0(:,:) = F%VelPhi(:,:)
        F%ForcesFl(:,:) = F%Phi(:,:)
        return
    end subroutine Prediction_Face_VelPhi

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Face), intent (INOUT) F
    !! \param real, intent (IN) dt
    !<
    subroutine Correction_Face_Veloc (F, dt)
        implicit none

        type (Face), intent (INOUT) :: F
        !real, intent (IN) :: bega, gam1
        real, intent (IN) :: dt
        integer :: i

        do i = 0,2
            F%Forces(:,:,i) = F%MassMat(:,:) * F%Forces(:,:,i)
        enddo

        F%Veloc(:,:,:) = F%v0(:,:,:) + dt * F%Forces(:,:,:)
        F%Accel(:,:,:) = (F%Veloc(:,:,:)-F%V0(:,:,:))/dt
        F%Displ(:,:,:) = F%Displ(:,:,:) +  dt * F%Veloc(:,:,:)

        return
    end subroutine Correction_Face_Veloc

    subroutine Correction_Face_VelPhi(F,dt)
        implicit none

        type(Face), intent(inout) :: F
        real, intent(in) :: dt

        F%ForcesFl(:,:) = F%MassMat(:,:) * F%ForcesFl(:,:)
        F%VelPhi(:,:) = F%VelPhi0(:,:) + dt * F%ForcesFl(:,:)
        if(F%fluid_dirich) F%VelPhi = 0d0
        F%AccelPhi(:,:) = (F%VelPhi(:,:)-F%VelPhi0(:,:))/dt
        F%Phi(:,:) = F%Phi(:,:) + dt * F%VelPhi(:,:)
        return
    end subroutine Correction_Face_VelPhi


    ! ##########################################################
    !>
    !! \brief
    !!
    !! \param type (Face), intent (INOUT) F
    !! \param real, intent (IN) dt
    !<
    subroutine Correction_Face_PML_Veloc (F, dt)

        implicit none

        type (Face), intent (INOUT) :: F
        real, intent (IN) :: dt

        integer :: i


        do i = 0,2
            F%spml%Veloc1(:,:,i) = F%spml%DumpVx(:,:,0) * F%spml%Veloc1(:,:,i) + dt * F%spml%DumpVx(:,:,1) * F%spml%Forces1(:,:,i)
            F%spml%Veloc2(:,:,i) = F%spml%DumpVy(:,:,0) * F%spml%Veloc2(:,:,i) + dt * F%spml%DumpVy(:,:,1) * F%spml%Forces2(:,:,i)
            F%spml%Veloc3(:,:,i) = F%spml%DumpVz(:,:,0) * F%spml%Veloc3(:,:,i) + dt * F%spml%DumpVz(:,:,1) * F%spml%Forces3(:,:,i)
        enddo

        F%Veloc = F%spml%Veloc1 + F%spml%Veloc2 + F%spml%Veloc3

        F%Displ(:,:,:) = F%Displ(:,:,:) +  dt * F%Veloc(:,:,:)

        if (F%Abs) then
            F%Veloc = 0
        endif

        F%V0 = F%Veloc
        
        return
    end subroutine Correction_Face_PML_Veloc

    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Correction_Face_PML_VelPhi(F,dt)

        implicit none

        type(Face), intent(inout) :: F
        real, intent(in) :: dt


        F%spml%VelPhi1(:,:) = F%spml%DumpVx(:,:,0) * F%spml%VelPhi1(:,:) +     &
            dt * F%spml%DumpVx(:,:,1) * F%spml%ForcesFl1(:,:)
        F%spml%VelPhi2(:,:) = F%spml%DumpVy(:,:,0) * F%spml%VelPhi2(:,:) +     &
            dt * F%spml%DumpVy(:,:,1) * F%spml%ForcesFl2(:,:)
        F%spml%VelPhi3(:,:) = F%spml%DumpVz(:,:,0) * F%spml%VelPhi3(:,:) +     &
            dt * F%spml%DumpVz(:,:,1) * F%spml%ForcesFl3(:,:)

        F%VelPhi = F%spml%VelPhi1 + F%spml%VelPhi2 + F%spml%VelPhi3

        if(F%Abs .or. F%fluid_dirich)then
            F%VelPhi = 0
        endif

        F%Phi = F%Phi+dt*F%Velphi

        return
    end subroutine Correction_Face_PML_VelPhi

    ! ############################################################
    !>
    !! \fn subroutine Correction_Face_FPML_Veloc (F, dt, fil)
    !! \brief
    !!
    !! \param type (Face toto), intent (INOUT) F
    !! \param real, intent (IN toto) dt
    !! \param real, intent (IN toto) fil
    !<


    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Correction_Face_FPML_Veloc (F, dt, fil)
        implicit none

        type (Face), intent (INOUT) :: F
        real, intent (IN) :: dt, fil

        integer :: i
        real :: fil2
        real, dimension (1:F%ngll1-2,1:F%ngll2-2) :: Ausiliar_velocity

        fil2 = fil**2

        do i = 0,2
            Ausiliar_Velocity = F%spml%Veloc1(:,:,i)
            F%spml%Veloc1(:,:,i) = F%spml%DumpVx(:,:,0) * F%spml%Veloc1(:,:,i) + dt * F%spml%DumpVx(:,:,1) * F%spml%Forces1(:,:,i) + F%spml%Ivx * F%spml%Iveloc1(:,:,i)
            F%spml%Iveloc1(:,:,i) = Fil2 * F%spml%Iveloc1(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%spml%Veloc1(:,:,i))

            Ausiliar_Velocity = F%spml%Veloc2(:,:,i)
            F%spml%Veloc2(:,:,i) = F%spml%DumpVy(:,:,0) * F%spml%Veloc2(:,:,i) + dt * F%spml%DumpVy(:,:,1) * F%spml%Forces2(:,:,i) + F%spml%Ivy * F%spml%Iveloc2(:,:,i)
            F%spml%Iveloc2(:,:,i) = Fil2 * F%spml%Iveloc2(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%spml%Veloc2(:,:,i))

            Ausiliar_Velocity = F%spml%Veloc3(:,:,i)
            F%spml%Veloc3(:,:,i) = F%spml%DumpVz(:,:,0) * F%spml%Veloc3(:,:,i) + dt * F%spml%DumpVz(:,:,1) * F%spml%Forces3(:,:,i) + F%spml%Ivz * F%spml%Iveloc3(:,:,i)
            F%spml%Iveloc3(:,:,i) = Fil2 * F%spml%Iveloc3(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%spml%Veloc3(:,:,i))
        enddo

        F%Veloc = F%spml%Veloc1 + F%spml%Veloc2 + F%spml%Veloc3

        if (F%Abs) then
            F%Veloc = 0
        endif

        return
    end subroutine Correction_Face_FPML_Veloc

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
