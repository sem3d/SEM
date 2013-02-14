!>
!! \file Edge.f90
!! \brief
!!
!<

module sedges

    type :: edge_pml
       real, dimension (:,:), allocatable :: DumpMass, DumpVx, DumpVy, DumpVz
       real, dimension (:,:), allocatable :: Veloc1, Veloc2, Veloc3
       real, dimension (:,:), allocatable :: Forces1, Forces2, Forces3
       real, dimension (:,:), allocatable :: IVeloc1, Iveloc2, Iveloc3
       real, dimension (:), allocatable :: Ivx, Ivy, Ivz
       real, dimension(:), allocatable :: ForcesFl1, ForcesFl2, ForcesFl3, VelPhi1, VelPhi2, VelPhi3
    end type edge_pml
    type :: edge

       logical :: PML, Abs, FPML

       integer :: ngll
       integer, dimension (:), allocatable :: Iglobnum_Edge,EdgeNum
       !integer, dimension (:), allocatable :: Which_Elem,Which_EdgeinElem

       real, dimension (:), allocatable  :: MassMat
       real, dimension (:,:), allocatable :: Forces, Displ, Veloc, Accel, V0
       logical  :: solid
       ! solid-fluid
       real, dimension(:), allocatable :: ForcesFl, Phi, VelPhi, AccelPhi, VelPhi0

       type(edge_pml), allocatable :: spml
#ifdef MKA3D
       real, dimension (:,:), allocatable :: ForcesMka
       !     integer, dimension (:,:), allocatable :: FlagMka
       real, dimension (:), allocatable :: tsurfsem
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
    subroutine Correction_Edge_Veloc (E, dt)
        implicit none

        type(Edge), intent(inout) :: E
        real, intent (in) :: dt
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
            E%spml%Veloc1(:,i) = E%spml%DumpVx(:,0) * E%spml%Veloc1(:,i) + dt * E%spml%DumpVx(:,1) * E%spml%Forces1(:,i)
            E%spml%Veloc2(:,i) = E%spml%DumpVy(:,0) * E%spml%Veloc2(:,i) + dt * E%spml%DumpVy(:,1) * E%spml%Forces2(:,i)
            E%spml%Veloc3(:,i) = E%spml%DumpVz(:,0) * E%spml%Veloc3(:,i) + dt * E%spml%DumpVz(:,1) * E%spml%Forces3(:,i)
        enddo

        E%Veloc = E%spml%Veloc1 + E%spml%Veloc2 + E%spml%Veloc3
        E%Displ(:,:) = E%Displ(:,:) + dt * E%Veloc(:,:)

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

        E%spml%VelPhi1(:) = E%spml%DumpVx(:,0) * E%spml%VelPhi1(:) + dt * E%spml%DumpVx(:,1) * E%spml%ForcesFl1(:)
        E%spml%VelPhi2(:) = E%spml%DumpVy(:,0) * E%spml%VelPhi2(:) + dt * E%spml%DumpVy(:,1) * E%spml%ForcesFl2(:)
        E%spml%VelPhi3(:) = E%spml%DumpVz(:,0) * E%spml%VelPhi3(:) + dt * E%spml%DumpVz(:,1) * E%spml%ForcesFl3(:)

        E%VelPhi = E%spml%VelPhi1 + E%spml%VelPhi2 + E%spml%VelPhi3

        if(E%Abs)then
            E%VelPhi = 0
        endif

        return
    end subroutine Correction_Edge_PML_VelPhi

    ! ############################################
    subroutine Correction_Edge_FPML_Veloc (E, dt, fil)

        implicit none

        type (Edge), intent (INOUT) :: E
        real, intent (IN) :: dt, fil

        integer :: i
        real :: fil2
        real, dimension (1:E%ngll-1) :: Ausiliar_velocity

        fil2 = fil**2
        do i = 0,2
            Ausiliar_velocity = E%spml%Veloc1(:,i)
            E%spml%Veloc1(:,i) = E%spml%DumpVx(:,0) * E%spml%Veloc1(:,i) + dt * E%spml%DumpVx(:,1) * E%spml%Forces1(:,i) + E%spml%Ivx * E%spml%Iveloc1(:,i)
            E%spml%Iveloc1(:,i) = Fil2 * E%spml%Iveloc1(:,i) + 0.5 * (1-Fil2) *  (Ausiliar_velocity + E%spml%Veloc1(:,i))

            Ausiliar_velocity = E%spml%Veloc2(:,i)
            E%spml%Veloc2(:,i) = E%spml%DumpVy(:,0) * E%spml%Veloc2(:,i) + dt * E%spml%DumpVy(:,1) * E%spml%Forces2(:,i) + E%spml%Ivy * E%spml%Iveloc2(:,i)
            E%spml%Iveloc2(:,i) = Fil2 * E%spml%Iveloc2(:,i) + 0.5 * (1-Fil2) *  (Ausiliar_velocity + E%spml%Veloc2(:,i))

            Ausiliar_velocity = E%spml%Veloc3(:,i)
            E%spml%Veloc3(:,i) = E%spml%DumpVz(:,0) * E%spml%Veloc3(:,i) + dt * E%spml%DumpVz(:,1) * E%spml%Forces3(:,i) + E%spml%Ivz * E%spml%Iveloc3(:,i)
            E%spml%Iveloc3(:,i) = Fil2 * E%spml%Iveloc3(:,i) + 0.5 * (1-Fil2) *  (Ausiliar_velocity + E%spml%Veloc3(:,i))
        enddo

        E%Veloc = E%spml%Veloc1 + E%spml%Veloc2 + E%spml%Veloc3

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
                        Vfree(1:ngll-2,i) = Vfree(1:ngll-2,i) - ( E%spml%DumpVz(1:ngll-2,0) * E%spml%Veloc3(1:ngll-2,i) &
                            + dt * E%spml%DumpVz(1:ngll-2,1) * E%spml%Forces3(1:ngll-2,i) )
                    enddo
                else
                    do i = 0,2
                        do j = 1, ngll-2
                            Vfree(j,i) = Vfree(j,i) - ( E%spml%DumpVz(ngll-1-j,0) * E%spml%Veloc3(ngll-1-j,i) + dt * E%spml%DumpVz(ngll-1-j,1) * E%spml%Forces3(ngll-1-j,i) )
                        enddo
                    enddo
                endif
            else
                do i = 0,2
                    Vfree(1:ngll-2,i) =  E%spml%DumpVz(1:ngll-2,0) * E%spml%Veloc3(1:ngll-2,i) + dt * E%spml%DumpVz(1:ngll-2,1) * E%spml%Forces3(1:ngll-2,i)
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
    end subroutine init_edge

end module sedges
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
