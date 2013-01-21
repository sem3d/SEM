!>
!! \file Edge.f90
!! \brief
!!
!<

module sedges

    type :: edge

       logical :: PML, Abs, FPML

       integer :: ngll
       integer, dimension (:), pointer :: Iglobnum_Edge,EdgeNum,Which_Elem,Which_EdgeinElem

       real, dimension (:), pointer  :: MassMat
       real, dimension (:,:), pointer :: Forces, Displ, Veloc, Accel, V0
       real, dimension (:,:), pointer :: DumpMass, DumpVx, DumpVy, DumpVz
       real, dimension (:,:), pointer :: Veloc1, Veloc2, Veloc3
       real, dimension (:,:), pointer :: Forces1, Forces2, Forces3
       real, dimension (:,:), pointer :: IVeloc1, Iveloc2, Iveloc3
       real, dimension (:), pointer :: Ivx, Ivy, Ivz
       logical  :: solid
       ! solid-fluid
       real, dimension(:), pointer :: ForcesFl, Phi, VelPhi, AccelPhi, VelPhi0
       real, dimension(:), pointer :: ForcesFl1, ForcesFl2, ForcesFl3, VelPhi1, VelPhi2, VelPhi3

#ifdef MKA3D
       real, dimension (:,:), pointer :: ForcesMka
       !     integer, dimension (:,:), pointer :: FlagMka
       real, dimension (:), pointer :: tsurfsem
#endif


    end type edge

contains

    ! ############################################################
    !subroutine Prediction_Edge_Veloc (E, alpha, bega, dt)
    subroutine Prediction_Edge_Veloc (E,dt)
        implicit none

        type (Edge), intent (INOUT) :: E
        real, intent(in) :: dt

        E%Forces(:,:) = E%Displ(:,:)
        E%V0(:,:) = E%Veloc(:,:)

        return
    end subroutine Prediction_Edge_Veloc

    !------------------------------------------------------------------
    !------------------------------------------------------------------
    subroutine Prediction_Edge_VelPhi(E,dt)
        implicit none

        type(Edge), intent(inout) :: E
        real, intent(in) :: dt

        E%VelPhi0(:) = E%VelPhi(:)
        E%ForcesFl(:) = E%Phi(:)

        return
    end subroutine Prediction_Edge_VelPhi

    ! ###########################################################
    !subroutine Correction_Edge_Veloc (E, bega, gam1, dt)
    subroutine Correction_Edge_Veloc (E, dt)
        implicit none

        type (Edge), intent (INOUT) :: E
        real, intent (IN) :: dt
        integer :: i, ngll, j

        real :: xmas

#ifdef MKA3D
        xmas = 0.
        if (  E%tsurfsem(1) > 0. ) then
            xmas = 1.
            ! ??? Masse forfaitaire = 1/2 pour faces de couplage ?
            !    test volI2b
            !              xmas = 0.
        endif
#endif

        ngll = E%ngll

        do i = 0,2
            do j=1,ngll-2
#ifdef MKA3D
                !   cas inter1
                !          E%Forces(j,i) = E%MassMat(j) * E%Forces(j,i)/(1.+ E%FlagMka(j,i))
                E%Forces(j,i) = E%MassMat(j) * E%Forces(j,i)/(1.+ xmas)
#else
                E%Forces(j,i) =  E%MassMat(j) * E%Forces(j,i)
#endif
            enddo
        enddo

        !  modification schema en temps mariotti pour couplage
        !E%Veloc(:,:) = E%v0(:,:) + dt * E%Forces(:,:)
        !E%Accel(:,:) = E%Accel(:,:) + gam1 / dt * (E%Veloc(:,:)-E%V0(:,:))
        !E%Displ(:,:) = E%Displ(:,:) + bega * dt * (E%Veloc(:,:)+E%V0(:,:))
        E%Veloc(:,:) = E%v0(:,:) + dt * E%Forces(:,:)
        E%Accel(:,:) = (E%Veloc(:,:)-E%V0(:,:))/dt
        E%Displ(:,:) = E%Displ(:,:) + dt * E%Veloc(:,:)

        return
    end subroutine Correction_Edge_Veloc

    !------------------------------------------------------------------
    !------------------------------------------------------------------
    subroutine Correction_Edge_VelPhi(E,dt)
        implicit none

        type(Edge), intent(inout) :: E
        real, intent(in) :: dt

        E%ForcesFl(:) = E%MassMat(:) * E%ForcesFl(:)

        E%VelPhi(:) = E%VelPhi0(:) + dt * E%ForcesFl(:)
        E%AccelPhi(:) = (E%VelPhi(:)-E%VelPhi0(:))/dt
        E%Phi(:) = E%Phi(:) + dt * E%VelPhi(:)


        return
    end subroutine Correction_Edge_VelPhi

    ! ############################################################
    subroutine Correction_Edge_PML_Veloc (E, dt)

        implicit none

        type (Edge), intent (INOUT) :: E
        real, intent (IN) :: dt

        integer :: i

        do i = 0,2
            E%Veloc1(:,i) = E%DumpVx(:,0) * E%Veloc1(:,i) + dt * E%DumpVx(:,1) * E%Forces1(:,i)
            E%Veloc2(:,i) = E%DumpVy(:,0) * E%Veloc2(:,i) + dt * E%DumpVy(:,1) * E%Forces2(:,i)
            E%Veloc3(:,i) = E%DumpVz(:,0) * E%Veloc3(:,i) + dt * E%DumpVz(:,1) * E%Forces3(:,i)
        enddo

        E%Veloc = E%Veloc1 + E%Veloc2 + E%Veloc3

        if (E%Abs) then
            E%Veloc = 0
        endif

        return
    end subroutine Correction_Edge_PML_Veloc

    !------------------------------------------------------------------
    !------------------------------------------------------------------
    subroutine Correction_Edge_PML_VelPhi(E,dt)

        implicit none

        type(Edge), intent(inout) :: E
        real, intent(in) :: dt

        E%VelPhi1(:) = E%DumpVx(:,0) * E%VelPhi1(:) + dt * E%DumpVx(:,1) * E%ForcesFl1(:)
        E%VelPhi2(:) = E%DumpVy(:,0) * E%VelPhi2(:) + dt * E%DumpVy(:,1) * E%ForcesFl2(:)
        E%VelPhi3(:) = E%DumpVz(:,0) * E%VelPhi3(:) + dt * E%DumpVz(:,1) * E%ForcesFl3(:)

        E%VelPhi = E%VelPhi1 + E%VelPhi2 + E%VelPhi3

        if(E%Abs)then
            E%VelPhi = 0
        endif

        return
    end subroutine Correction_Edge_PML_VelPhi

    ! ############################################
    subroutine Correction_Edge_FPML_Veloc (E, dt)

        implicit none

        type (Edge), intent (INOUT) :: E
        real, intent (IN) :: dt

        integer :: i
        real :: fil2
        real, dimension (1:E%ngll-1) :: Ausiliar_velocity

        do i = 0,2
            Ausiliar_velocity = E%Veloc1(:,i)
            E%Veloc1(:,i) = E%DumpVx(:,0) * E%Veloc1(:,i) + dt * E%DumpVx(:,1) * E%Forces1(:,i) + E%Ivx * E%Iveloc1(:,i)
            E%Iveloc1(:,i) = Fil2 * E%Iveloc1(:,i) + 0.5 * (1-Fil2) *  (Ausiliar_velocity + E%Veloc1(:,i))

            Ausiliar_velocity = E%Veloc2(:,i)
            E%Veloc2(:,i) = E%DumpVy(:,0) * E%Veloc2(:,i) + dt * E%DumpVy(:,1) * E%Forces2(:,i) + E%Ivy * E%Iveloc2(:,i)
            E%Iveloc2(:,i) = Fil2 * E%Iveloc2(:,i) + 0.5 * (1-Fil2) *  (Ausiliar_velocity + E%Veloc2(:,i))

            Ausiliar_velocity = E%Veloc3(:,i)
            E%Veloc3(:,i) = E%DumpVz(:,0) * E%Veloc3(:,i) + dt * E%DumpVz(:,1) * E%Forces3(:,i) + E%Ivz * E%Iveloc3(:,i)
            E%Iveloc3(:,i) = Fil2 * E%Iveloc3(:,i) + 0.5 * (1-Fil2) *  (Ausiliar_velocity + E%Veloc3(:,i))
        enddo

        E%Veloc = E%Veloc1 + E%Veloc2 + E%Veloc3

        if (E%Abs) then
            E%Veloc = 0
        endif

        return
    end subroutine Correction_Edge_FPML_Veloc

    ! ############################################

    subroutine get_vel_edge (E, Vfree, ngll, dt, logic, orient)
        implicit none

        integer, intent (IN) :: ngll, orient
        real, intent (IN) :: dt
        type (Edge), intent (IN) :: E
        real, dimension (1:ngll-2,0:2), intent (INOUT) :: Vfree
        logical, intent (IN) :: logic

        integer :: i,j

        if (.not. E%PML) then

            if (logic) then
                if ( orient ==0 ) then
                    do i = 0,2
                        Vfree(1:ngll-2,i) = Vfree(1:ngll-2,i) - ( E%V0(1:ngll-2,i) + dt*E%MassMat(1:ngll-2)*E%Forces(1:ngll-2,i) )
                    enddo
                else
                    do i = 0,2
                        do j = 1, ngll-2
                            Vfree(j,i) = Vfree(j,i) - ( E%V0(ngll-1-j,i) + dt*E%MassMat(ngll-1-j)*E%Forces(ngll-1-j,i) )
                        enddo
                    enddo
                endif
            else
                do i = 0,2
                    Vfree(1:ngll-2,i) =  E%V0(1:ngll-2,i) + dt*E%MassMat(1:ngll-2)*E%Forces(1:ngll-2,i)
                enddo
            endif

        else

            if (logic) then
                if ( orient ==0 ) then
                    do i = 0,2
                        Vfree(1:ngll-2,i) = Vfree(1:ngll-2,i) - ( E%DumpVz(1:ngll-2,0) * E%Veloc3(1:ngll-2,i) &
                            + dt * E%DumpVz(1:ngll-2,1) * E%Forces3(1:ngll-2,i) )
                    enddo
                else
                    do i = 0,2
                        do j = 1, ngll-2
                            Vfree(j,i) = Vfree(j,i) - ( E%DumpVz(ngll-1-j,0) * E%Veloc3(ngll-1-j,i) + dt * E%DumpVz(ngll-1-j,1) * E%Forces3(ngll-1-j,i) )
                        enddo
                    enddo
                endif
            else
                do i = 0,2
                    Vfree(1:ngll-2,i) =  E%DumpVz(1:ngll-2,0) * E%Veloc3(1:ngll-2,i) + dt * E%DumpVz(1:ngll-2,1) * E%Forces3(1:ngll-2,i)
                enddo
            endif

        endif

        return
    end subroutine get_vel_edge

    ! ###########################################################
    subroutine init_edge(ed)
        type(Edge), intent(inout) :: ed

        ed%PML = .false.
        ed%Abs = .false.
        ed%FPML = .false.
        ed%ngll = 0
        ed%solid = .true.
        nullify(ed%Which_Elem)
        nullify(ed%Which_EdgeinElem)
        nullify(ed%EdgeNum)
        nullify(ed%Iglobnum_Edge)
        nullify(ed%MassMat)
        nullify(ed%Forces)
        nullify(ed%Displ)
        nullify(ed%Veloc)
        nullify(ed%Accel)
        nullify(ed%V0)
        nullify(ed%Forces1)
        nullify(ed%Forces2)
        nullify(ed%Forces3)
        nullify(ed%Veloc1)
        nullify(ed%Veloc2)
        nullify(ed%Veloc3)
        nullify(ed%DumpVx)
        nullify(ed%DumpVy)
        nullify(ed%DumpVz)
        nullify(ed%DumpMass)
        nullify(ed%IVeloc1)
        nullify(ed%IVeloc2)
        nullify(ed%IVeloc3)
        nullify(ed%Ivx)
        nullify(ed%Ivy)
        nullify(ed%Ivz)
        nullify(ed%ForcesFl)
        nullify(ed%Phi)
        nullify(ed%VelPhi)
        nullify(ed%AccelPhi)
        nullify(ed%VelPhi0)
        nullify(ed%ForcesFl1)
        nullify(ed%ForcesFl2)
        nullify(ed%ForcesFl3)
        nullify(ed%VelPhi1)
        nullify(ed%VelPhi2)
        nullify(ed%VelPhi3)
#ifdef MKA3D
        nullify(ed%ForcesMka)
        nullify(ed%tsurfsem)
#endif
    end subroutine init_edge

end module sedges
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
