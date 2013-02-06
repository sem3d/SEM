!>
!!\file Face.f90
!!\brief Gère les faces des éléments.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sfaces

    ! Modified by Gaetano Festa 24/2/2005
    ! Modified by Paul Cupillard 06/11/2005

    type :: face
       logical :: PML, Abs, FPML
       integer :: ngll1, ngll2, dir, Which_Elem
       integer, dimension (:), pointer :: FaceNum
       integer, dimension (:,:), pointer :: Iglobnum_Face
       real, dimension (:,:), pointer  :: MassMat
       real, dimension (:,:,:), pointer :: Forces, Displ, Veloc, Accel, V0
       real, dimension (:,:,:), pointer :: Forces1, Forces2, Forces3, Veloc1, Veloc2, Veloc3
       real, dimension (:,:,:), pointer :: DumpVx, DumpVy, DumpVz, DumpMass
       real, dimension (:,:,:), pointer :: IVeloc1, IVeloc2, IVeloc3
       real, dimension (:,:), pointer :: Ivx, Ivy, Ivz
       logical :: solid
       ! solid-fluid
       real, dimension(:,:), pointer :: ForcesFl, Phi, VelPhi, AccelPhi, VelPhi0
       real, dimension(:,:), pointer :: ForcesFl1, ForcesFl2, ForcesFl3, VelPhi1, VelPhi2, VelPhi3
#ifdef MKA3D
       real, dimension (:,:,:), pointer :: ForcesMka
       real, dimension (:,:), pointer :: tsurfsem
#endif
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
        real :: xmas

#ifdef MKA3D
        xmas = 0.
        if (  F%tsurfsem(1,1) > 0. ) then
            xmas = 1.
        endif
#endif

        do i = 0,2
#ifdef MKA3D
            F%Forces(:,:,i) = F%MassMat(:,:) * F%Forces(:,:,i)/(1.+xmas)
#else
            F%Forces(:,:,i) = F%MassMat(:,:) * F%Forces(:,:,i)
#endif
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
            F%Veloc1(:,:,i) = F%DumpVx(:,:,0) * F%Veloc1(:,:,i) + dt * F%DumpVx(:,:,1) * F%Forces1(:,:,i)
            F%Veloc2(:,:,i) = F%DumpVy(:,:,0) * F%Veloc2(:,:,i) + dt * F%DumpVy(:,:,1) * F%Forces2(:,:,i)
            F%Veloc3(:,:,i) = F%DumpVz(:,:,0) * F%Veloc3(:,:,i) + dt * F%DumpVz(:,:,1) * F%Forces3(:,:,i)
        enddo

        F%Veloc = F%Veloc1 + F%Veloc2 + F%Veloc3

        F%Displ(:,:,:) = F%Displ(:,:,:) +  dt * F%Veloc(:,:,:)

        if (F%Abs) then
            F%Veloc = 0
        endif

        return
    end subroutine Correction_Face_PML_Veloc

    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine Correction_Face_PML_VelPhi(F,dt)

        implicit none

        type(Face), intent(inout) :: F
        real, intent(in) :: dt


        F%VelPhi1(:,:) = F%DumpVx(:,:,0) * F%VelPhi1(:,:) +     &
            dt * F%DumpVx(:,:,1) * F%ForcesFl1(:,:)
        F%VelPhi2(:,:) = F%DumpVy(:,:,0) * F%VelPhi2(:,:) +     &
            dt * F%DumpVy(:,:,1) * F%ForcesFl2(:,:)
        F%VelPhi3(:,:) = F%DumpVz(:,:,0) * F%VelPhi3(:,:) +     &
            dt * F%DumpVz(:,:,1) * F%ForcesFl3(:,:)

        F%VelPhi = F%VelPhi1 + F%VelPhi2 + F%VelPhi3

        if(F%Abs)then
            F%VelPhi = 0
        endif

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
            Ausiliar_Velocity = F%Veloc1(:,:,i)
            F%Veloc1(:,:,i) = F%DumpVx(:,:,0) * F%Veloc1(:,:,i) + dt * F%DumpVx(:,:,1) * F%Forces1(:,:,i) + F%Ivx * F%Iveloc1(:,:,i)
            F%Iveloc1(:,:,i) = Fil2 * F%Iveloc1(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%Veloc1(:,:,i))

            Ausiliar_Velocity = F%Veloc2(:,:,i)
            F%Veloc2(:,:,i) = F%DumpVy(:,:,0) * F%Veloc2(:,:,i) + dt * F%DumpVy(:,:,1) * F%Forces2(:,:,i) + F%Ivy * F%Iveloc2(:,:,i)
            F%Iveloc2(:,:,i) = Fil2 * F%Iveloc2(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%Veloc2(:,:,i))

            Ausiliar_Velocity = F%Veloc3(:,:,i)
            F%Veloc3(:,:,i) = F%DumpVz(:,:,0) * F%Veloc3(:,:,i) + dt * F%DumpVz(:,:,1) * F%Forces3(:,:,i) + F%Ivz * F%Iveloc3(:,:,i)
            F%Iveloc3(:,:,i) = Fil2 * F%Iveloc3(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%Veloc3(:,:,i))
        enddo

        F%Veloc = F%Veloc1 + F%Veloc2 + F%Veloc3

        if (F%Abs) then
            F%Veloc = 0
        endif

        return
    end subroutine Correction_Face_FPML_Veloc


    ! ############################################################
    subroutine get_vel_face (F, Vfree, ngll1, ngll2, dt, logic, orient)
        implicit none

        integer, intent (IN) :: ngll1, ngll2, orient
        real, intent (IN) :: dt
        type (Face), intent (IN) :: F
        real, dimension (1:ngll1-2,1:ngll2-2,0:2), intent (INOUT) :: Vfree
        logical, intent (IN) :: logic

        integer :: i,j

        if (.not. F%PML) then

            if (logic) then
                select case (orient)
                case (0)
                    do i = 0,2
                        Vfree(1:ngll1-2,1:ngll2-2,i) = Vfree(1:ngll1-2,1:ngll2-2,i) - ( F%V0(1:ngll1-2,1:ngll2-2,i) + &
                            dt*F%MassMat(1:ngll1-2,1:ngll2-2)*F%Forces(1:ngll1-2,1:ngll2-2,i) )
                    enddo
                case (1)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-i,j,0:2) + dt*F%MassMat(ngll1-1-i,j)*F%Forces(ngll1-1-i,j,0:2) )
                        enddo
                    enddo
                case (2)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(i,ngll2-1-j,0:2) + dt*F%MassMat(i,ngll2-1-j)*F%Forces(i,ngll2-1-j,0:2) )
                        enddo
                    enddo
                case (3)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-i,ngll2-1-j,0:2) + dt*F%MassMat(ngll1-1-i,ngll2-1-j)*F%Forces(ngll1-1-i,ngll2-1-j,0:2) )
                        enddo
                    enddo
                case (4)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(j,i,0:2) + dt*F%MassMat(j,i)*F%Forces(j,i,0:2) )
                        enddo
                    enddo
                case (5)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-j,i,0:2) + dt*F%MassMat(ngll1-1-j,i)*F%Forces(ngll1-1-j,i,0:2) )
                        enddo
                    enddo
                case (6)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(j,ngll2-1-i,0:2) + dt*F%MassMat(j,ngll2-1-i)*F%Forces(j,ngll2-1-i,0:2) )
                        enddo
                    enddo
                case (7)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-j,ngll2-1-i,0:2) + dt*F%MassMat(ngll1-1-j,ngll2-1-i)*F%Forces(ngll1-1-j,ngll2-1-i,0:2) )
                        enddo
                    enddo
                end select
            else
                do i = 0,2
                    Vfree(1:ngll1-2,1:ngll2-2,i) =  F%V0(1:ngll1-2,1:ngll2-2,i) + dt*F%MassMat(1:ngll1-2,1:ngll2-2)*F%Forces(1:ngll1-2,1:ngll2-2,i)
                enddo
            endif

        else

            if (logic) then
                select case (orient)
                case (0)
                    do i = 0,2
                        Vfree(1:ngll1-2,1:ngll2-2,i) = Vfree(1:ngll1-2,1:ngll2-2,i) - ( F%DumpVz(1:ngll1-2,1:ngll2-2,0) * F%Veloc3(1:ngll1-2,1:ngll2-2,i) + &
                            dt * F%DumpVz(1:ngll1-2,1:ngll2-2,1) * F%Forces3(1:ngll1-2,1:ngll2-2,i) )
                    enddo
                case (1)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll1-1-i,j,0) * F%Veloc3(ngll1-1-i,j,0:2) + &
                                dt * F%DumpVz(ngll1-1-i,j,1) * F%Forces3(ngll1-1-i,j,0:2) )
                        enddo
                    enddo
                case (2)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(i,ngll2-1-j,0) * F%Veloc3(i,ngll2-1-j,0:2) + &
                                dt * F%DumpVz(i,ngll2-1-j,1) * F%Forces3(i,ngll2-1-j,0:2) )
                        enddo
                    enddo
                case (3)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll1-1-i,ngll2-1-j,0) * F%Veloc3(ngll1-1-i,ngll2-1-j,0:2) + &
                                dt * F%DumpVz(ngll1-1-i,ngll2-1-j,1) * F%Forces3(ngll1-1-i,ngll2-1-j,0:2) )
                        enddo
                    enddo
                case (4)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(j,i,0) * F%Veloc3(j,i,0:2) + &
                                dt * F%DumpVz(j,i,1) * F%Forces3(j,i,0:2) )
                        enddo
                    enddo
                case (5)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll2-1-j,i,0) * F%Veloc3(ngll2-1-j,i,0:2) + &
                                dt * F%DumpVz(ngll2-1-j,i,1) * F%Forces3(ngll2-1-j,i,0:2) )
                        enddo
                    enddo
                case (6)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(j,ngll1-1-i,0) * F%Veloc3(j,ngll1-1-i,0:2) + &
                                dt * F%DumpVz(j,ngll1-1-i,1) * F%Forces3(j,ngll1-1-i,0:2) )
                        enddo
                    enddo
                case (7)
                    do j = 1, ngll2-2
                        do i = 1, ngll1-2
                            Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll2-1-j,ngll1-1-i,0) * F%Veloc3(ngll2-1-j,ngll1-1-i,0:2) + &
                                dt * F%DumpVz(ngll2-1-j,ngll1-1-i,1) * F%Forces3(ngll2-1-j,ngll1-1-i,0:2) )
                        enddo
                    enddo
                end select
            else
                do i = 0,2
                    Vfree(1:ngll1-2,1:ngll2-2,i) =  F%DumpVz(1:ngll1-2,1:ngll2-2,0) * F%Veloc3(1:ngll1-2,1:ngll2-2,i) + &
                        dt * F%DumpVz(1:ngll1-2,1:ngll2-2,1) * F%Forces3(1:ngll1-2,1:ngll2-2,i)
                enddo
            endif

        endif

        return
    end subroutine get_vel_face

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
        nullify(fc%FaceNum)
        nullify(fc%Iglobnum_Face)
        nullify(fc%MassMat)
        nullify(fc%Forces)
        nullify(fc%Displ)
        nullify(fc%Veloc)
        nullify(fc%Accel)
        nullify(fc%V0)
        nullify(fc%Forces1)
        nullify(fc%Forces2)
        nullify(fc%Forces3)
        nullify(fc%Veloc1)
        nullify(fc%Veloc2)
        nullify(fc%Veloc3)
        nullify(fc%DumpVx)
        nullify(fc%DumpVy)
        nullify(fc%DumpVz)
        nullify(fc%DumpMass)
        nullify(fc%IVeloc1)
        nullify(fc%IVeloc2)
        nullify(fc%IVeloc3)
        nullify(fc%Ivx)
        nullify(fc%Ivy)
        nullify(fc%Ivz)
        nullify(fc%ForcesFl)
        nullify(fc%Phi)
        nullify(fc%VelPhi)
        nullify(fc%AccelPhi)
        nullify(fc%VelPhi0)
        nullify(fc%ForcesFl1)
        nullify(fc%ForcesFl2)
        nullify(fc%ForcesFl3)
        nullify(fc%VelPhi1)
        nullify(fc%VelPhi2)
        nullify(fc%VelPhi3)
#ifdef MKA3D
        nullify(fc%ForcesMka)
        nullify(fc%tsurfsem)
#endif
    end subroutine init_face

end module sfaces
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
