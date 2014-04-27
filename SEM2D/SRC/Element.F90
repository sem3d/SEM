!>
!!\file Element.F90
!!\brief contient les méthodes qui assure la gestion du type Element.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module selement
    implicit none

    ! Last Modification 8/11/2004

    type :: element

       integer :: mat_index,ngllx,ngllz
       integer, dimension (:), allocatable :: Control_nodes
       integer, dimension (0:3) :: Near_Face   !
       integer, dimension (0:3) :: Near_Vertex !
       integer, dimension (:,:), allocatable :: Iglobnum

       real, dimension (:,:), allocatable :: Jacob,Density, Lambda, Mu,MassMat
       real, dimension(:,:,:), allocatable :: Forces,Stress,Veloc,Displ,Accel,V0
       real, dimension(:,:,:), allocatable :: ACoeff
       real, dimension(:,:,:,:), allocatable :: InvGrad
       logical :: OUTPUT

       ! PML allocation
       logical :: PML
       real, dimension (:,:,:), allocatable :: Stress1,Stress2
       real, dimension (:,:,:), allocatable :: Veloc1, Veloc2,Forces1,Forces2
       real, dimension (:,:,:), allocatable :: DumpMass
       real, dimension (:,:,:), allocatable :: DumpSx,DumpSz,DumpVx,DumpVz

       ! FPML allocation
       logical :: FPML
       real, dimension (:,:), allocatable ::  Isx, Isz, Ivx, Ivz
       real, dimension (:,:,:), allocatable :: Istress1, IStress2, Iveloc1, Iveloc2

       ! CPML allocation
       logical :: CPML
       real, dimension (:,:), allocatable ::  Ax, Bx, Az, Bz

       real dist_max !!Ajout Gsa 03/10 - taille caracteristique de l'element
    end type element

contains

    subroutine init_element(el)
        implicit none
        type(Element), intent(INOUT) :: el

        el%mat_index=-1
        el%PML = .false.
        el%OUTPUT = .true.
        el%FPML = .false.
        el%dist_max = 0.0

    end subroutine init_element

    ! ############################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<


    !subroutine Prediction_Elem_Veloc (Elem,alpha,bega, dt)
    subroutine Prediction_Elem_Veloc (Elem)
        implicit none

        type (Element), intent (INOUT) :: Elem
        !real, intent (IN) :: bega, alpha, dt

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz
        !  test mariotti
        !     Elem%Forces(1:ngllx-2,1:ngllz-2,0:1) = Elem%Displ + dt * Elem%Veloc + dt**2 * (0.5 - bega) * Elem%Accel
        !     Elem%V0 = Elem%Veloc
        !     Elem%Forces(1:ngllx-2,1:ngllz-2,0:1) = alpha * Elem%Forces(1:ngllx-2,1:ngllz-2,0:1) + (1-alpha) * Elem%Displ
        Elem%Forces(1:ngllx-2,1:ngllz-2,0:1) = Elem%Displ
        Elem%V0 = Elem%Veloc

        return
    end subroutine Prediction_Elem_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !<


    !  subroutine Correction_Elem_Veloc (Elem, bega, gam1,dt)
    subroutine Correction_Elem_Veloc (Elem, dt)
        implicit none

        type (Element), intent (INOUT) :: Elem
        !real, intent (IN) :: bega, gam1
        real, intent (IN) :: dt

        integer :: ngllx, ngllz,i

        ! test

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz
        do i = 0,1
            Elem%Forces(1:ngllx-2,1:ngllz-2,i) = Elem%MassMat(:,:)  * Elem%Forces(1:ngllx-2,1:ngllz-2,i)
        enddo
        !       print*,' valeur Mass apres ',Elem%Forces(1,1,0)
        !  test mariotti
        !      Elem%Veloc(:,:,:)  = Elem%v0(:,:,:)+ dt * Elem%Forces(1:ngllx-2,1:ngllz-2,:)
        !      Elem%Accel  = Elem%Accel + gam1 /dt * (Elem%Veloc-Elem%V0)
        !      Elem%Displ = Elem%Displ + bega *dt * (Elem%Veloc+Elem%V0)
        Elem%Veloc(:,:,:)  = Elem%v0(:,:,:)+ dt * Elem%Forces(1:ngllx-2,1:ngllz-2,:)
        Elem%Accel  = (Elem%Veloc-Elem%V0)/dt
        Elem%Displ = Elem%Displ + dt * Elem%Veloc

        return
    end subroutine Correction_Elem_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) dt
    !<


    subroutine Correction_Elem_PML_Veloc (Elem, dt)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) ::  dt
        integer :: ngllx, ngllz,i

        ngllx = Elem%ngllx; ngllz=Elem%ngllz

        do i = 0,1
            Elem%Veloc1(:,:,i) = Elem%DumpVx(:,:,0) * Elem%Veloc1(:,:,i) + dt * &
                Elem%DumpVx(:,:,1)*Elem%Forces1(1:ngllx-2,1:ngllz-2,i)
            Elem%Veloc2(:,:,i) =Elem%DumpVz(:,:,0) * Elem%Veloc2(:,:,i) + Dt * &
                Elem%DumpVz(:,:,1)*Elem%Forces2(1:ngllx-2,1:ngllz-2,i)
        enddo
        Elem%V0 = Elem%Veloc
        Elem%Veloc = Elem%Veloc1 + Elem%Veloc2

        return
    end subroutine Correction_Elem_PML_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) fil
    !<


    subroutine Correction_Elem_FPML_Veloc (Elem, dt, fil)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) ::  dt, fil

        integer :: ngllx, ngllz,i
        real :: fil2
        real, dimension (1:Elem%ngllx-2,1:Elem%ngllz-2) :: Ausiliar_Velocity

        ngllx = Elem%ngllx; ngllz=Elem%ngllz
        fil2 = fil**2

        do i = 0,1
            Ausiliar_Velocity = Elem%Veloc1(:,:,i)
            Elem%Veloc1(:,:,i) = Elem%DumpVx(:,:,0) * Elem%Veloc1(:,:,i) + dt * &
                Elem%DumpVx(:,:,1)*Elem%Forces1(1:ngllx-2,1:ngllz-2,i) + Elem%Ivx * Elem%Iveloc1(:,:,i)
            Elem%Iveloc1(:,:,i) = Fil2 * Elem%Iveloc1(:,:,i) + 0.5 * (1.-Fil2) *  &
                (Ausiliar_Velocity + Elem%Veloc1(:,:,i) )

            Ausiliar_Velocity = Elem%Veloc2(:,:,i)
            Elem%Veloc2(:,:,i) =Elem%DumpVz(:,:,0) * Elem%Veloc2(:,:,i) + Dt * &
                Elem%DumpVz(:,:,1)*Elem%Forces2(1:ngllx-2,1:ngllz-2,i) + Elem%Ivz * Elem%IVeloc2(:,:,i)
            Elem%Iveloc2(:,:,i) = Fil2 * Elem%Iveloc2(:,:,i) + 0.5 * (1.-Fil2) * &
                (Ausiliar_Velocity + Elem%Veloc2(:,:,i) )
        enddo

        Elem%V0 = Elem%Veloc
        Elem%Veloc = Elem%Veloc1 + Elem%Veloc2

        return
    end subroutine Correction_Elem_FPML_Veloc
    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTmat
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hmatz
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vxloc
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vzloc
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<


    subroutine Prediction_Elem_PML_Veloc (Elem,alpha, bega, dt,Vxloc,Vzloc,Hmatz, HTmat)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) ::  hTmat
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hmatz
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT)  ::Vxloc, Vzloc
        real, intent (IN) :: bega, dt, alpha

        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: s0,s1,s2,s3

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx; ngllz = Elem%ngllz

        !Elem%Stress0 = Elem%Stress

        VxLoc(1:ngllx-2,1:ngllz-2)  = (0.5+alpha) * Elem%Veloc(:,:,0) + dt *(0.5-bega) *Elem%Accel(:,:,0) + (0.5-alpha) * Elem%V0(:,:,0)
        VzLoc(1:ngllx-2,1:ngllz-2)  =(0.5+alpha) * Elem%Veloc(:,:,1) + dt *(0.5-bega) *Elem%Accel(:,:,1)+ (0.5-alpha) * Elem%V0(:,:,1)


        s0 = MATMUL (HTmat,VxLoc)
        s2 = MATMUL (HTmat,VzLoc)
        s1 = MATMUL (VxLoc,Hmatz)
        s3 = MATMUL (VzLoc,Hmatz)

        Elem%Stress1(:,:,0) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,0) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,0) * s0 + &
            Elem%Acoeff(:,:,2) * s2 )
        Elem%Stress2(:,:,0) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,0) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,3) * s3 + &
            Elem%Acoeff(:,:,1)* s1 )

        Elem%Stress1(:,:,1) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,1) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,4) * s0 + &
            Elem%Acoeff(:,:,6) * s2 )
        Elem%Stress2(:,:,1) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,1) + Elem%DumpSz(:,:,1) * Dt *( Elem%Acoeff(:,:,7) * s3 + &
            Elem%Acoeff(:,:,5) * s1 )

        Elem%Stress1(:,:,2) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,2) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,10) * s2 + &
            Elem%Acoeff(:,:,8) * s0)
        Elem%Stress2(:,:,2) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,2) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,9) * s1 + &
            Elem%Acoeff(:,:,11) * s3 )

        Elem%Stress = Elem%Stress1 + Elem%Stress2
        !Elem%Stress = alpha  * Elem%Stress + (1-alpha) * Elem%Stress0
        return
    end subroutine Prediction_Elem_PML_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTmat
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hmatz
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vxloc
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vzloc
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !! \param real, intent (IN) fil
    !<


    subroutine Prediction_Elem_FPML_Veloc (Elem,alpha, bega, dt,Vxloc,Vzloc,Hmatz, HTmat,fil)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) ::  hTmat
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hmatz
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT)  ::Vxloc, Vzloc
        real, intent (IN) :: bega, dt, alpha, fil

        real :: fil2
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: s0,s1,s2,s3, Stress_ausiliar

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx; ngllz = Elem%ngllz
        fil2 = fil**2

        !Elem%Stress0 = Elem%Stress

        VxLoc(1:ngllx-2,1:ngllz-2)  = (0.5+alpha) * Elem%Veloc(:,:,0) + dt *(0.5-bega) *Elem%Accel(:,:,0) + (0.5-alpha) * Elem%V0(:,:,0)
        VzLoc(1:ngllx-2,1:ngllz-2)  =(0.5+alpha) * Elem%Veloc(:,:,1) + dt *(0.5-bega) *Elem%Accel(:,:,1)+ (0.5-alpha) * Elem%V0(:,:,1)


        s0 = MATMUL (HTmat,VxLoc)
        s2 = MATMUL (HTmat,VzLoc)
        s1 = MATMUL (VxLoc,Hmatz)
        s3 = MATMUL (VzLoc,Hmatz)

        Stress_Ausiliar = Elem%Stress1(:,:,0)
        Elem%Stress1(:,:,0) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,0) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,0) * s0 + &
            Elem%Acoeff(:,:,2) * s2 ) + Elem%Isx * Elem%Istress1(:,:,0)
        Elem%Istress1 (:,:,0) =  Fil2 * Elem%Istress1 (:,:,0) + 0.5 * (1.-Fil2) * (Elem%Stress1(:,:,0) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress2(:,:,0)
        Elem%Stress2(:,:,0) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,0) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,3) * s3 + &
            Elem%Acoeff(:,:,1)* s1 )  + Elem%Isz * Elem%Istress2(:,:,0)
        Elem%Istress2 (:,:,0) =  Fil2 * Elem%Istress2 (:,:,0) + 0.5 * (1.-Fil2) * (Elem%Stress2(:,:,0) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress1(:,:,1)
        Elem%Stress1(:,:,1) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,1) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,4) * s0 + &
            Elem%Acoeff(:,:,6) * s2 ) + Elem%Isx * Elem%Istress1(:,:,1)
        Elem%Istress1 (:,:,1) =  Fil2 * Elem%Istress1 (:,:,1) + 0.5 * (1.-Fil2) * (Elem%Stress1(:,:,1) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress2(:,:,1)
        Elem%Stress2(:,:,1) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,1) + Elem%DumpSz(:,:,1) * Dt *( Elem%Acoeff(:,:,7) * s3 + &
            Elem%Acoeff(:,:,5) * s1 ) + Elem%Isz * Elem%Istress2(:,:,1)
        Elem%Istress2 (:,:,1) =  Fil2 * Elem%Istress2 (:,:,1) + 0.5 * (1.-Fil2) * (Elem%Stress2(:,:,1) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress1(:,:,2)
        Elem%Stress1(:,:,2) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,2) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,10) * s2 + &
            Elem%Acoeff(:,:,8) * s0) + Elem%Isx * Elem%Istress1(:,:,2)
        Elem%Istress1 (:,:,2) =  Fil2 * Elem%Istress1 (:,:,2) + 0.5 * (1.-Fil2) * (Elem%Stress1(:,:,2) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress2(:,:,2)
        Elem%Stress2(:,:,2) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,2) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,9) * s1 + &
            Elem%Acoeff(:,:,11) * s3 ) + Elem%Isz * Elem%Istress2(:,:,2)
        Elem%Istress2 (:,:,2) =  Fil2 * Elem%Istress2 (:,:,2) + 0.5 * (1.-Fil2) * (Elem%Stress2(:,:,2) + Stress_Ausiliar)

        Elem%Stress = Elem%Stress1 + Elem%Stress2
        return
    end subroutine Prediction_Elem_FPML_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTmat
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hmatz
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vxloc
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vzloc
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !! \param real, intent (IN) fil
    !<
    subroutine Prediction_Elem_CPML_Veloc (Elem, alpha, bega, dt, Vxloc, Vzloc, Hmatz, HTmat, fil)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) ::  hTmat
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hmatz
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) ::Vxloc, Vzloc
        real, intent (IN) :: bega, dt, alpha, fil

        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: s0,s1,s2,s3, Stress_ausiliar

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx; ngllz = Elem%ngllz

        VxLoc(1:ngllx-2,1:ngllz-2) = (0.5+alpha) * Elem%Veloc(:,:,0) &
                                   + dt *(0.5-bega)*Elem%Accel(:,:,0) + (0.5-alpha)*Elem%V0(:,:,0)
        VzLoc(1:ngllx-2,1:ngllz-2) = (0.5+alpha) * Elem%Veloc(:,:,1) &
                                   + dt *(0.5-bega)*Elem%Accel(:,:,1) + (0.5-alpha)*Elem%V0(:,:,1)
        s0 = MATMUL (HTmat,VxLoc)
        s2 = MATMUL (HTmat,VzLoc)
        s1 = MATMUL (VxLoc,Hmatz)
        s3 = MATMUL (VzLoc,Hmatz)

        Elem%Stress(:,:,0) = Elem%Stress(:,:,0) +  Dt * (Elem%Acoeff(:,:,0) * s0 + Elem%Acoeff(:,:,2) * s2 )

######################### following is OLD ########################
        Elem%Stress1(:,:,0) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,0) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,0) * s0 + &
            Elem%Acoeff(:,:,2) * s2 )
        Elem%Stress2(:,:,0) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,0) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,3) * s3 + &
            Elem%Acoeff(:,:,1)* s1 )

        Elem%Stress1(:,:,1) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,1) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,4) * s0 + &
            Elem%Acoeff(:,:,6) * s2 )
        Elem%Stress2(:,:,1) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,1) + Elem%DumpSz(:,:,1) * Dt *( Elem%Acoeff(:,:,7) * s3 + &
            Elem%Acoeff(:,:,5) * s1 )

        Elem%Stress1(:,:,2) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,2) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,10) * s2 + &
            Elem%Acoeff(:,:,8) * s0)
        Elem%Stress2(:,:,2) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,2) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,9) * s1 + &
            Elem%Acoeff(:,:,11) * s3 )

        return
    end subroutine Prediction_Elem_CPML_Veloc

    ! ###########################################################

    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<


    subroutine  compute_InternalForces_Elem (Elem,hprime, hTprime, hprimez,hTprimez)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime, hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: Uxloc, Uzloc, dUx_dxi, dUx_deta,  dUz_dxi, dUz_deta, s0

        Uxloc =Elem%Forces (:,:,0)
        Uzloc = Elem%Forces (:,:,1)

        dUx_dxi = MATMUL ( hTprime, Uxloc)
        dUz_dxi = MATMUL ( hTprime, Uzloc )

        dUx_deta = MATMUL ( Uxloc , hprimez )
        dUz_deta = MATMUL ( Uzloc , hprimez )

        s0 = Elem%Acoeff(:,:,0)*dUx_dxi + Elem%Acoeff(:,:,1)*dUx_deta + Elem%Acoeff(:,:,2)*dUz_dxi + Elem%Acoeff(:,:,3)*dUz_deta
        Uxloc = MATMUL ( hprime, s0 )

        s0 = Elem%Acoeff(:,:,1)*dUx_dxi + Elem%Acoeff(:,:,4)*dUx_deta + Elem%Acoeff(:,:,5)*dUz_dxi + Elem%Acoeff(:,:,6)*dUz_deta
        Uxloc = Uxloc +  MATMUL( s0 , hTprimez )

        s0 = Elem%Acoeff(:,:,2)*dUx_dxi + Elem%Acoeff(:,:,5)*dUx_deta + Elem%Acoeff(:,:,7)*dUz_dxi + Elem%Acoeff(:,:,8)*dUz_deta
        Uzloc = MATMUL ( hprime, s0  )

        s0= Elem%Acoeff(:,:,3)*dUx_dxi + Elem%Acoeff(:,:,6)*dUx_deta + Elem%Acoeff(:,:,8)*dUz_dxi + Elem%Acoeff(:,:,9)*dUz_deta
        Uzloc = Uzloc + MATMUL( s0 , hTprimez )

        Elem%Forces(:,:,0) = Uxloc
        Elem%Forces(:,:,1) = Uzloc

        return
    end subroutine compute_InternalForces_Elem

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<


    subroutine  compute_InternalForces_PML_Elem (Elem,hprime, hTprimez)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: s0,s1

        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,0) + Elem%Acoeff(:,:,14) * Elem%Stress(:,:,2)
        s1= MATMUL (hprime,s0)
        Elem%Forces1(:,:,0) = s1
        s0 = Elem%Acoeff(:,:,15) * Elem%Stress (:,:,2)  + Elem%Acoeff(:,:,13) * Elem%Stress (:,:,0)
        s1  =   MATMUL (s0,hTprimez)

        Elem%Forces2(:,:,0) = s1
        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,2) + Elem%Acoeff(:,:,14)* Elem%Stress(:,:,1)
        s1 = MATMUL (hprime,s0)
        Elem%Forces1(:,:,1) = s1
        s0 = Elem%Acoeff(:,:,15) * Elem%Stress(:,:,1) + Elem%Acoeff(:,:,13) * Elem%Stress (:,:,2)
        s1  =   MATMUL (s0,hTprimez)
        Elem%Forces2(:,:,1) = s1

        Elem%Forces = Elem%Forces1 + Elem%Forces2

        return
    end subroutine compute_InternalForces_PML_Elem


    ! ###########################################################
    !>
    !! \brief This subroutine computes the energy of the elastic deformation of the element
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !! \param real, intent (INOUT) E_elas
    !<

    subroutine  compute_Elastic_Energy (Elem, hTprime, hprimez, E_elas)
        implicit none

        type (Element), intent (IN) :: Elem
        real, intent (INOUT) :: E_elas
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
        real, dimension ( 0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: dUx_dxi, dUx_deta,  dUz_dxi, dUz_deta
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: s0, EMat, Uxloc, Uzloc

        ! This subroutine is called outside the Newmark scheme, so the local
        ! displacements are store in Forces
        Uxloc =Elem%Forces (:,:,0)
        Uzloc = Elem%Forces (:,:,1)

        dUx_dxi = MATMUL ( hTprime, Uxloc)
        dUz_dxi = MATMUL ( hTprime, Uzloc )
        dUx_deta = MATMUL ( Uxloc , hprimez )
        dUz_deta = MATMUL ( Uzloc , hprimez )

         s0 =  Elem%Acoeff(:,:,0)*dUx_dxi + Elem%Acoeff(:,:,1)*dUx_deta &
             + Elem%Acoeff(:,:,2)*dUz_dxi + Elem%Acoeff(:,:,3)*dUz_deta
         EMat = dUx_dxi * s0

         s0 =  Elem%Acoeff(:,:,2)*dUx_dxi + Elem%Acoeff(:,:,5)*dUx_deta &
             + Elem%Acoeff(:,:,7)*dUz_dxi + Elem%Acoeff(:,:,8)*dUz_deta
         EMat = EMat + dUz_dxi * s0

         s0 =  Elem%Acoeff(:,:,1)*dUx_dxi + Elem%Acoeff(:,:,4)*dUx_deta &
             + Elem%Acoeff(:,:,5)*dUz_dxi + Elem%Acoeff(:,:,6)*dUz_deta
         EMat = EMat + dUx_deta * s0

         s0 =  Elem%Acoeff(:,:,3)*dUx_dxi + Elem%Acoeff(:,:,6)*dUx_deta &
             + Elem%Acoeff(:,:,8)*dUz_dxi + Elem%Acoeff(:,:,9)*dUz_deta
         EMat = EMat + dUz_deta * s0

         E_elas = -0.5 * sum(EMat)

    end subroutine compute_Elastic_Energy


    ! ###########################################################
    !>
    !! \brief This subroutine computes the kinetic energy of the inner nodes
    !!  of an element
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (INOUT) E_kin
    !<

    subroutine  compute_Kinetic_Energy (Elem, Dt, E_kin)
        implicit none

        type (Element), intent (IN) :: Elem
        real, intent (IN)    :: Dt
        real, intent (INOUT) :: E_kin
        real, dimension (1:Elem%ngllx-2, 1:Elem%ngllz-2)      :: Ener_Mat
        real, dimension (1:Elem%ngllx-2, 1:Elem%ngllz-2, 0:1) :: Vel_half
        integer :: ngllx, ngllz

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz

        Vel_half(:,:,:) = Elem%Veloc(:,:,:) + 0.5 * dt * Elem%Forces(1:ngllx-2,1:ngllz-2,:)
        Ener_Mat (:,:)  = 1./Elem%MassMat(:,:) * ( Vel_half(:,:,0)*Vel_half(:,:,0) &
                                                  +Vel_half(:,:,1)*Vel_half(:,:,1))
        E_kin = 0.5 * sum(Ener_Mat)

    end subroutine compute_Kinetic_Energy

end module selement
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
