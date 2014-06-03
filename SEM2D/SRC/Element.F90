!>
!!\file Element.F90
!!\brief contient les méthodes qui assure la gestion du type Element.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module selement
    use constants
    implicit none


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
       real, dimension (:,:), allocatable :: Axi, Bxi, Aeta, Beta, Axi_prime, Aeta_prime
       real, dimension (:,:), allocatable :: PsiVxxi, PsiVxeta, PsiVzxi, PsiVzeta, PsiSxxxi
       real, dimension (:,:), allocatable :: PsiSxxeta, PsiSzzxi, PsiSzzeta, PsiSxzxi, PsiSxzeta

       real dist_max !!Ajout Gsa 03/10 - taille caracteristique de l'element

       ! DG
       integer :: type_DG
       logical :: acoustic
       real, dimension (:),    allocatable :: Coeff_Integr_Faces
       real, dimension(:,:,:), allocatable :: Strain
       real, dimension(:,:,:), allocatable :: Vect_RK
       ! HDG
       real, dimension(:,:), allocatable :: Normal_nodes
       real, dimension(:,:), allocatable :: MatPen, TracFace, Vhat

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
        el%type_DG = GALERKIN_CONT

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
    subroutine Prediction_Elem_CPML_Veloc (Elem, alpha, bega, dt, Vxloc, Vzloc, Hmatz, HTmat)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) ::  hTmat
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hmatz
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) ::Vxloc, Vzloc
        real, intent (IN) :: bega, dt, alpha

        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: s0,s1,s2,s3

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

        ! Updating convolution :
        Elem%PsiVxxi (:,:) = Elem%Bxi (:,:) * Elem%PsiVxxi (:,:) + Elem%Axi (:,:) * s0
        Elem%PsiVxeta(:,:) = Elem%Beta(:,:) * Elem%PsiVxeta(:,:) + Elem%Aeta(:,:) * s1
        Elem%PsiVzxi (:,:) = Elem%Bxi (:,:) * Elem%PsiVzxi (:,:) + Elem%Axi (:,:) * s2
        Elem%PsiVzeta(:,:) = Elem%Beta(:,:) * Elem%PsiVzeta(:,:) + Elem%Aeta(:,:) * s3

        ! Updating PML Stresses in the PML
        Elem%Stress(:,:,0)= Elem%Stress(:,:,0) + Dt*(Elem%Acoeff(:,:,0)*(s0 + Elem%PsiVxxi (:,:)) &
                                                    +Elem%Acoeff(:,:,2)*(s2 + Elem%PsiVzxi (:,:)) &
                                                    +Elem%Acoeff(:,:,3)*(s3 + Elem%PsiVzeta(:,:)) &
                                                    +Elem%Acoeff(:,:,1)*(s1 + Elem%PsiVxeta(:,:)))
        Elem%Stress(:,:,1)= Elem%Stress(:,:,1) + Dt*(Elem%Acoeff(:,:,4)*(s0 + Elem%PsiVxxi (:,:)) &
                                                    +Elem%Acoeff(:,:,6)*(s2 + Elem%PsiVzxi (:,:)) &
                                                    +Elem%Acoeff(:,:,7)*(s3 + Elem%PsiVzeta(:,:)) &
                                                    +Elem%Acoeff(:,:,5)*(s1 + Elem%PsiVxeta(:,:)))
        Elem%Stress(:,:,2)= Elem%Stress(:,:,2) + Dt*(Elem%Acoeff(:,:,8)*(s0 + Elem%PsiVxxi (:,:)) &
                                                    +Elem%Acoeff(:,:,10)*(s2+ Elem%PsiVzxi (:,:)) &
                                                    +Elem%Acoeff(:,:,11)*(s3+ Elem%PsiVzeta(:,:)) &
                                                    +Elem%Acoeff(:,:,9)*(s1 + Elem%PsiVxeta(:,:)))
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
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1) :: Uxloc, Uzloc, dUx_dxi, dUx_deta,  dUz_dxi, dUz_deta, s0

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


    subroutine  compute_InternalForces_DG_Weak (Elem,hprime,hTprimez)
      implicit none

      type (Element), intent (INOUT) :: Elem
      real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime
      real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
      real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: aux1, aux2

      aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,0)
      aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,0)
      Elem%Forces(:,:,0) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      aux1 = Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,1)
      aux2 = Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,1)
      Elem%Forces(:,:,1) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,0)
      aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,0)
      Elem%Forces(:,:,2) = 0.5 * (MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez))

      aux1 = (Elem%Acoeff(:,:,4) + Elem%Acoeff(:,:,5))*Elem%Strain(:,:,0) + Elem%Acoeff(:,:,4)*Elem%Strain(:,:,1) &
           + Elem%Acoeff(:,:,7)*Elem%Strain(:,:,2)
      aux2 = (Elem%Acoeff(:,:,8) + Elem%Acoeff(:,:,9))*Elem%Strain(:,:,0) + Elem%Acoeff(:,:,8)*Elem%Strain(:,:,1) &
           + Elem%Acoeff(:,:,11)*Elem%Strain(:,:,2)
      Elem%Forces(:,:,3) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      aux1 = (Elem%Acoeff(:,:,7) + Elem%Acoeff(:,:,6))*Elem%Strain(:,:,1) + Elem%Acoeff(:,:,6)*Elem%Strain(:,:,0) &
           + Elem%Acoeff(:,:,5)*Elem%Strain(:,:,2)
      aux2 = (Elem%Acoeff(:,:,11) + Elem%Acoeff(:,:,10))*Elem%Strain(:,:,1) + Elem%Acoeff(:,:,10)*Elem%Strain(:,:,0) &
           + Elem%Acoeff(:,:,9)*Elem%Strain(:,:,2)
      Elem%Forces(:,:,4) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      Elem%Forces(:,:,:) = - Elem%Forces(:,:,:)

      return
    end subroutine compute_InternalForces_DG_Weak


    ! ###########################################################

    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !<

    subroutine  compute_InternalForces_DG_Strong (Elem,hTprime,hprimez)
      implicit none

      type (Element), intent (INOUT) :: Elem
      real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
      real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
      real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: aux1, aux2

      aux1 = MATMUL(hTprime,Elem%Veloc(:,:,0))
      aux2 = MATMUL(Elem%Veloc(:,:,0),hprimez)
      Elem%Forces(:,:,0) = Elem%Acoeff(:,:,0)*aux1 + Elem%Acoeff(:,:,1)*aux2

      aux1 = MATMUL(hTprime,Elem%Veloc(:,:,1))
      aux2 = MATMUL(Elem%Veloc(:,:,1),hprimez)
      Elem%Forces(:,:,1) = Elem%Acoeff(:,:,2)*aux1 + Elem%Acoeff(:,:,3)*aux2

      aux1 = Elem%Acoeff(:,:,0) * MATMUL(hTprime,Elem%Veloc(:,:,1)) +  Elem%Acoeff(:,:,1) * MATMUL(Elem%Veloc(:,:,1),hprimez)
      aux2 = Elem%Acoeff(:,:,2) * MATMUL(hTprime,Elem%Veloc(:,:,0)) +  Elem%Acoeff(:,:,3) * MATMUL(Elem%Veloc(:,:,0),hprimez)
      Elem%Forces(:,:,2) = 0.5 * (aux1 + aux2)

      aux1 = (Elem%Acoeff(:,:,4) + Elem%Acoeff(:,:,5)) * MATMUL(hTprime,Elem%Strain(:,:,0)) &
           + Elem%Acoeff(:,:,4) * MATMUL(hTprime,Elem%Strain(:,:,1)) &
           + Elem%Acoeff(:,:,7) * MATMUL(hTprime,Elem%Strain(:,:,2))
      aux2 = (Elem%Acoeff(:,:,8) + Elem%Acoeff(:,:,9)) * MATMUL(Elem%Strain(:,:,0),hprimez) &
           + Elem%Acoeff(:,:,8) * MATMUL(Elem%Strain(:,:,1),hprimez) &
           + Elem%Acoeff(:,:,11)* MATMUL(Elem%Strain(:,:,2),hprimez)
      Elem%Forces(:,:,3) = aux1 + aux2

      aux1 = (Elem%Acoeff(:,:,7) + Elem%Acoeff(:,:,6)) * MATMUL(hTprime,Elem%Strain(:,:,1)) &
           + Elem%Acoeff(:,:,6) * MATMUL(hTprime,Elem%Strain(:,:,0)) &
           + Elem%Acoeff(:,:,5) * MATMUL(hTprime,Elem%Strain(:,:,2))
      aux2 = (Elem%Acoeff(:,:,11) + Elem%Acoeff(:,:,10)) * MATMUL(Elem%Strain(:,:,1),hprimez) &
           + Elem%Acoeff(:,:,10)* MATMUL(Elem%Strain(:,:,0),hprimez) &
           + Elem%Acoeff(:,:,9) * MATMUL(Elem%Strain(:,:,2),hprimez)
      Elem%Forces(:,:,4) = aux1 + aux2

      return
    end subroutine compute_InternalForces_DG_Strong


    ! ###########################################################
    !>
    !! \brief Thus subroutine adds the contribution of the memory PML terms
    !! for the ADE-PML in case of DG or HDG.
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine  add_Psi4PML (Elem)
        implicit none

        type (Element), intent (INOUT) :: Elem

        Elem%Forces(:,:,0) = Elem%Forces(:,:,0) + Elem%PsiVxxi (:,:) * Elem%Acoeff(:,:,0) &
                                                + Elem%PsiVxeta(:,:) * Elem%Acoeff(:,:,1) &
        Elem%Forces(:,:,1) = Elem%Forces(:,:,1) + Elem%PsiVzxi (:,:) * Elem%Acoeff(:,:,2) &
                                                + Elem%PsiVzeta(:,:) * Elem%Acoeff(:,:,3)
        Elem%Forces(:,:,2) = Elem%Forces(:,:,2) + 0.5*(Elem%PsiVxxi (:,:) * Elem%Acoeff(:,:,2) &
                                                      +Elem%PsiVxeta(:,:) * Elem%Acoeff(:,:,3) &
                                                      +Elem%PsiVzxi (:,:) * Elem%Acoeff(:,:,0) &
                                                      +Elem%PsiVzeta(:,:) * Elem%Acoeff(:,:,1))
        Elem%Forces(:,:,3) = Elem%Forces(:,:,3) + Elem%PsiSxxxi (:,:) * Elem%Acoeff(:,:,0) &
                                                + Elem%PsiSxxeta(:,:) * Elem%Acoeff(:,:,1) &
                                                + Elem%PsiSxzxi (:,:) * Elem%Acoeff(:,:,2) &
                                                + Elem%PsiSxzeta(:,:) * Elem%Acoeff(:,:,3)
        Elem%Forces(:,:,4) = Elem%Forces(:,:,4) + Elem%PsiSxzxi (:,:) * Elem%Acoeff(:,:,0) &
                                                + Elem%PsiSxzeta(:,:) * Elem%Acoeff(:,:,1) &
                                                + Elem%PsiSzzxi (:,:) * Elem%Acoeff(:,:,2) &
                                                + Elem%PsiSzzeta(:,:) * Elem%Acoeff(:,:,3)
        return
    end subroutine add_Psi4PML


    ! ###########################################################
    !>
    !! \brief Thus subroutine is used in a ADE-PML framework :
    !! the memory terms PsiS*** and PsiV*** are updated for the next iteration.
    !! In a RK4 framework, the ADE for the Psi variables is solved using RK4,
    !! like the others variables.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) coeff1
    !! \param real, intent (IN) coeff2
    !<
    subroutine  update_Psi_RK4 (Elem,coeff1,coeff2,Dt)
        implicit none

        type(Element), intent (INOUT) :: Elem
        real, intent(IN) :: coeff1
        real, intent(IN) :: coeff2
        real, intent(IN) :: Dt
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1) :: PsiTrace
        real, dimension(0:ngllx+1,0:ngllz+1,0:9) :: smbr
        integer          :: ngx, ngz

        ! Defining Second Member of time evolution equation for the Psi
        Saux = Elem%Stain(:,:,0)
        smbr(:,:,0) = -Elem%Bx(:,:) * Elem%PsiSxxxi(:,:) - MATMUL(Htprimex,Saux)

        ! UPDATING in time using usual LSERK4
        Elem%Psi_RK(:,:,:)  = coeff1 * Elem%Psi_RK(:,:,:) + Dt*smbr(:,:,:)
        Elem%PsiSxxxi (:,:) = Elem%PsiSxxxi (:,:) + coeff2 * Elem%Psi_RK(:,:,0)
        Elem%PsiSxxeta(:,:) = Elem%PsiSxxeta(:,:) + coeff2 * Elem%Psi_RK(:,:,1)
        Elem%PsiSzzxi (:,:) = Elem%PsiSzzxi (:,:) + coeff2 * Elem%Psi_RK(:,:,2)
        Elem%PsiSzzeta(:,:) = Elem%PsiSzzeta(:,:) + coeff2 * Elem%Psi_RK(:,:,3)
        Elem%PsiSxzxi (:,:) = Elem%PsiSxzxi (:,:) + coeff2 * Elem%Psi_RK(:,:,4)
        Elem%PsiSxzeta(:,:) = Elem%PsiSxzeta(:,:) + coeff2 * Elem%Psi_RK(:,:,5)
        Elem%PsiVxxi (:,:)  = Elem%PsiVxxi (:,:)  + coeff2 * Elem%Psi_RK(:,:,6)
        Elem%PsiVxeta(:,:)  = Elem%PsiVxeta(:,:)  + coeff2 * Elem%Psi_RK(:,:,7)
        Elem%PsiVzxi (:,:)  = Elem%PsiVxxi (:,:)  + coeff2 * Elem%Psi_RK(:,:,8)
        Elem%PsiVzeta(:,:)  = Elem%PsiVxeta(:,:)  + coeff2 * Elem%Psi_RK(:,:,9)

    end subroutine update_Psi_RK4

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
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<

    subroutine  compute_InternalForces_CPML_Elem (Elem, hprime, hTprime, hprimez, hTprimez)
        implicit none

        type (Element), intent (INOUT) :: Elem

        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: s0,s1
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime, hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        s0 = MATMUL(hTprime,Elem%Axi)
        s1 = MATMUL(Elem%Aeta,hprimez)

        ! Updating convolution :
        Elem%PsiSxxxi (:,:) = Elem%Bxi (:,:) * Elem%PsiSxxxi (:,:) - Elem%Axi_prime (:,:) * Elem%Stress(:,:,0) & ! - Elem%Stress(:,:,0) * s0 &
            - Elem%Acoeff(:,:,17) * MATMUL(hprime, Elem%Stress(:,:,0)*Elem%Axi(:,:)*Elem%Acoeff(:,:,16))
        Elem%PsiSxxeta(:,:) = Elem%Beta(:,:) * Elem%PsiSxxeta(:,:) - Elem%Aeta_prime(:,:) * Elem%Stress(:,:,0) & ! - Elem%Stress(:,:,0) * s1 &
            - Elem%Acoeff(:,:,17) * MATMUL(Elem%Stress(:,:,0)*Elem%Aeta(:,:)*Elem%Acoeff(:,:,16), hTprimez)

        Elem%PsiSzzxi (:,:) = Elem%Bxi (:,:) * Elem%PsiSzzxi (:,:) - Elem%Axi_prime (:,:) * Elem%Stress(:,:,1) & ! - Elem%Stress(:,:,1) * s0 &
            - Elem%Acoeff(:,:,17) * MATMUL(hprime, Elem%Stress(:,:,1)*Elem%Axi(:,:)*Elem%Acoeff(:,:,16))
        Elem%PsiSzzeta(:,:) = Elem%Beta(:,:) * Elem%PsiSzzeta(:,:) - Elem%Aeta_prime(:,:) * Elem%Stress(:,:,1) & ! - Elem%Stress(:,:,1) * s1 &
            - Elem%Acoeff(:,:,17) * MATMUL(Elem%Stress(:,:,1)*Elem%Aeta(:,:)*Elem%Acoeff(:,:,16), hTprimez)

        Elem%PsiSxzxi (:,:) = Elem%Bxi (:,:) * Elem%PsiSxzxi (:,:) - Elem%Axi_prime (:,:) * Elem%Stress(:,:,2) & ! - Elem%Stress(:,:,2) * s0 &
            - Elem%Acoeff(:,:,17) * MATMUL(hprime, Elem%Stress(:,:,2)*Elem%Axi(:,:)*Elem%Acoeff(:,:,16))
        Elem%PsiSxzeta(:,:) = Elem%Beta(:,:) * Elem%PsiSxzeta(:,:) - Elem%Aeta_prime(:,:) * Elem%Stress(:,:,2) & ! - Elem%Stress(:,:,2) * s1 &
            - Elem%Acoeff(:,:,17) * MATMUL(Elem%Stress(:,:,2)*Elem%Aeta(:,:)*Elem%Acoeff(:,:,16), hTprimez)

        ! Updating Forces :
        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,0) + Elem%Acoeff(:,:,14) * Elem%Stress(:,:,2)
        s1 = Elem%Acoeff(:,:,13) * Elem%Stress(:,:,0) + Elem%Acoeff(:,:,15) * Elem%Stress(:,:,2)
        Elem%Forces(:,:,0) = MATMUL(hprime,s0) + MATMUL(s1,hTprimez) &
                           - Elem%Acoeff(:,:,12)*Elem%PsiSxxxi(:,:) - Elem%Acoeff(:,:,13)*Elem%PsiSxxeta(:,:) &
                           - Elem%Acoeff(:,:,14)*Elem%PsiSxzxi(:,:) - Elem%Acoeff(:,:,15)*Elem%PsiSxzeta(:,:)

        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,2) + Elem%Acoeff(:,:,14) * Elem%Stress(:,:,1)
        s1 = Elem%Acoeff(:,:,13) * Elem%Stress(:,:,2) + Elem%Acoeff(:,:,15) * Elem%Stress(:,:,1)
        Elem%Forces(:,:,1) = MATMUL(hprime,s0) + MATMUL(s1,hTprimez) &
                           - Elem%Acoeff(:,:,12)*Elem%PsiSxzxi(:,:) - Elem%Acoeff(:,:,13)*Elem%PsiSxzeta(:,:) &
                           - Elem%Acoeff(:,:,14)*Elem%PsiSzzxi(:,:) - Elem%Acoeff(:,:,15)*Elem%PsiSzzeta(:,:)

        return
    end subroutine compute_InternalForces_CPML_Elem

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine inversion_massmat(Elem)

      implicit none
      type (Element), intent (INOUT) :: Elem
      integer  :: ngllx,ngllz,i

      if(Elem%type_DG==GALERKIN_CONT) then
         ngllx = Elem%ngllx ; ngllz = Elem%ngllz
         do i = 0,1
            Elem%Forces(1:ngllx-2,1:ngllz-2,i) = Elem%MassMat(:,:)  * Elem%Forces(1:ngllx-2,1:ngllz-2,i)
         enddo
      else
         do i = 0,2
            Elem%Forces(:,:,i) = Elem%Acoeff(:,:,12) * Elem%Forces(:,:,i)
         enddo
         do i = 3,4
            Elem%Forces(:,:,i) = Elem%MassMat(:,:) * Elem%Forces(:,:,i)
         enddo
      endif

    end subroutine inversion_massmat


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
    !! \brief This subroutine computes the energy of the elastic deformation of the element.
    !!  Suitable for DG elements only.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (INOUT) E_elas
    !<

    subroutine  compute_Elastic_Energy_DG (Elem, E_elas)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (INOUT) :: E_elas
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: sigma11, sigma22, sigma12, EMat

        sigma11 = (Elem%Lambda + 2*Elem%Mu) * Elem%Strain(:,:,0) + Elem%Lambda * Elem%Strain(:,:,1)
        sigma22 = (Elem%Lambda + 2*Elem%Mu) * Elem%Strain(:,:,1) + Elem%Lambda * Elem%Strain(:,:,0)
        sigma12 = 2*Elem%Mu * Elem%Strain(:,:,2)

        EMat = Elem%Strain(:,:,0)*sigma11 + 2.*Elem%Strain(:,:,2)*sigma12 &
             + Elem%Strain(:,:,1)*sigma22
        EMat = EMat * 1. / Elem%Acoeff(:,:,12)  ! Pour multiplier par Jac*wx*wz
        E_elas = 0.5 * sum(EMat)

  end subroutine compute_Elastic_Energy_DG


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


    ! ###########################################################
    !>
    !! \brief This subroutine computes the kinetic energy using all the nodes
    !!  of an element. It suitable for Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (INOUT) E_kin
    !<

    subroutine  compute_Kinetic_Energy_DG (Elem, E_kin)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (INOUT) :: E_kin
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: Ener_Mat

        Ener_Mat = 1./Elem%MassMat(:,:) * ( Elem%Veloc(:,:,0)*Elem%Veloc(:,:,0) &
                                           +Elem%Veloc(:,:,1)*Elem%Veloc(:,:,1))
        E_kin = 0.5 * sum(Ener_Mat)

    end subroutine compute_Kinetic_Energy_DG


! ###########################################################
    !>
    !! \brief This subroutine computes the Penalization Matrix MatPen
    !! of an element on each external node of the element.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_MatPen (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        real, dimension (0:Elem%ngllx-1) :: Zp_x, Zs_x
        real, dimension (0:Elem%ngllz-1) :: Zp_z, Zs_z
        integer    :: imin, imax, ngllx, ngllz

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz

        ! Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Zp_x(:) = sqrt(Elem%Density(0:ngllx-1,0) *(Elem%Lambda(0:ngllx-1,0)+2.*Elem%Mu(0:ngllx-1,0)))
        Zs_x(:) = sqrt(Elem%Density(0:ngllx-1,0) * Elem%Mu(0:ngllx-1,0))
        Elem%MatPen(imin:imax,0) = Zp_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_x(:)-Zs_x(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

        ! Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Zp_z(:) = sqrt(Elem%Density(ngllx-1,0:ngllz-1) * (Elem%Lambda(ngllx-1,0:ngllz-1) &
                                                         + 2.*Elem%Mu(ngllx-1,0:ngllz-1)))
        Zs_z(:) = sqrt(Elem%Density(ngllx-1,0:ngllz-1) *  Elem%Mu(ngllx-1,0:ngllz-1))
        Elem%MatPen(imin:imax,0) = Zp_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_z(:)-Zs_z(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

        ! Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Zp_x(:) = sqrt(Elem%Density(0:ngllx-1,ngllz-1) * (Elem%Lambda(0:ngllx-1,ngllz-1) &
                                                         + 2.*Elem%Mu(0:ngllx-1,ngllz-1)))
        Zs_x(:) = sqrt(Elem%Density(0:ngllx-1,ngllz-1) *  Elem%Mu(0:ngllx-1,ngllz-1))
        Elem%MatPen(imin:imax,0) = Zp_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_x(:)-Zs_x(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

        ! Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Zp_z(:) = sqrt(Elem%Density(0,0:ngllz-1) *(Elem%Lambda(0,0:ngllz-1)+2.*Elem%Mu(0,0:ngllz-1)))
        Zs_z(:) = sqrt(Elem%Density(0,0:ngllz-1) * Elem%Mu(0,0:ngllz-1))
        Elem%MatPen(imin:imax,0) = Zp_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_z(:)-Zs_z(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

    end subroutine compute_MatPen


! ###########################################################
    !>
    !! \brief This subroutine computes the numerical Tractions at the
    !! face of an element on each external node of the element.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_TracFace (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        integer    :: imin, imax, ngx, ngz
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! For the Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(0:ngx-1,0) + 2*Elem%Mu(0:ngx-1,0)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(0:ngx-1,0,0) &
                                    + Elem%Lambda(0:ngx-1,0) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,0,1) &
                                    + 2*Elem%Mu(0:ngx-1,0)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,0,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(0:ngx-1,0) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,0,0) &
                                    +(Elem%Lambda(0:ngx-1,0) + 2*Elem%Mu(0:ngx-1,0)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(0:ngx-1,0,1) &
                                    + 2*Elem%Mu(0:ngx-1,0)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,0,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(0:ngx-1,0,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,0,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,0,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(0:ngx-1,0,1)
        ! For the Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(ngx-1,0:ngz-1) + 2*Elem%Mu(ngx-1,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(ngx-1,0:ngz-1,0) &
                                    + Elem%Lambda(ngx-1,0:ngz-1) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(ngx-1,0:ngz-1,1) &
                                    + 2*Elem%Mu(ngx-1,0:ngz-1)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(ngx-1,0:ngz-1,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(ngx-1,0:ngz-1) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(ngx-1,0:ngz-1,0) &
                                    +(Elem%Lambda(ngx-1,0:ngz-1) + 2*Elem%Mu(ngx-1,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(ngx-1,0:ngz-1,1) &
                                    + 2*Elem%Mu(ngx-1,0:ngz-1)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(ngx-1,0:ngz-1,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(ngx-1,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(ngx-1,0:ngz-1,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(ngx-1,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(ngx-1,0:ngz-1,1)
        ! For the Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(0:ngx-1,ngz-1) + 2*Elem%Mu(0:ngx-1,ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(0:ngx-1,ngz-1,0) &
                                    + Elem%Lambda(0:ngx-1,ngz-1) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,ngz-1,1) &
                                    + 2*Elem%Mu(0:ngx-1,ngz-1)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,ngz-1,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(0:ngx-1,ngz-1) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,ngz-1,0) &
                                    +(Elem%Lambda(0:ngx-1,ngz-1) + 2*Elem%Mu(0:ngx-1,ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(0:ngx-1,ngz-1,1) &
                                    + 2*Elem%Mu(0:ngx-1,ngz-1)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,ngz-1,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(0:ngx-1,ngz-1,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,ngz-1,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,ngz-1,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(0:ngx-1,ngz-1,1)
        ! For the Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(0,0:ngz-1) + 2*Elem%Mu(0,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(0,0:ngz-1,0) &
                                    + Elem%Lambda(0,0:ngz-1) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0,0:ngz-1,1) &
                                    + 2*Elem%Mu(0,0:ngz-1)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0,0:ngz-1,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(0,0:ngz-1) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0,0:ngz-1,0) &
                                    +(Elem%Lambda(0,0:ngz-1) + 2*Elem%Mu(0,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(0,0:ngz-1,1) &
                                    + 2*Elem%Mu(0,0:ngz-1)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0,0:ngz-1,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(0,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0,0:ngz-1,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(0,0:ngz-1,1)
    end subroutine compute_TracFace

! ###########################################################
    !>
    !! \brief This subroutine computes the strain traces at the
    !! on each external node of the element, and add these strains to the
    !! Elem%Forces after weighting for integrals computation.
    !! It is suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_Traces (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1) :: Vxz
        integer    :: imin, imax, ngx, ngz
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! Adding contribution of Vhat on the trace Traction :
        Elem%TracFace(:,0) = Elem%TracFace(:,0) + Elem%MatPen(:,0)*Elem%Vhat(:,0) &
                                                + Elem%MatPen(:,2)*Elem%Vhat(:,1)
        Elem%TracFace(:,1) = Elem%TracFace(:,1) + Elem%MatPen(:,2)*Elem%Vhat(:,0) &
                                                + Elem%MatPen(:,1)*Elem%Vhat(:,1)
        ! Weights for lineic integrals computation :
        Elem%TracFace(:,0) = Elem%TracFace(:,0) * Elem%Coeff_Integr_Faces(:)
        Elem%TracFace(:,1) = Elem%TracFace(:,1) * Elem%Coeff_Integr_Faces(:)

        Vxz(:) = 0.5 * ( Elem%Vhat(:,0) * Elem%Normal_Nodes(:,1) &
                       + Elem%Vhat(:,1) * Elem%Normal_Nodes(:,0))* Elem%Coeff_Integr_Faces(:)
        Elem%Vhat(:,0) = Elem%Vhat(:,0) * Elem%Normal_Nodes(:,0) * Elem%Coeff_Integr_Faces(:)
        Elem%Vhat(:,1) = Elem%Vhat(:,1) * Elem%Normal_Nodes(:,1) * Elem%Coeff_Integr_Faces(:)

        ! Adding Strain and velocities traces to Elem%Forces :
        ! For the Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Elem%Forces(0:ngx-1,0,0) = Elem%Forces(0:ngx-1,0,0) + Elem%Vhat(imin:imax,0)
        Elem%Forces(0:ngx-1,0,1) = Elem%Forces(0:ngx-1,0,1) + Elem%Vhat(imin:imax,1)
        Elem%Forces(0:ngx-1,0,2) = Elem%Forces(0:ngx-1,0,2) + Vxz(imin:imax)
        Elem%Forces(0:ngx-1,0,3) = Elem%Forces(0:ngx-1,0,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(0:ngx-1,0,4) = Elem%Forces(0:ngx-1,0,4) + Elem%TracFace(imin:imax,1)
        ! For the Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Elem%Forces(ngx-1,0:ngz-1,0) = Elem%Forces(ngx-1,0:ngz-1,0) + Elem%Vhat(imin:imax,0)
        Elem%Forces(ngx-1,0:ngz-1,1) = Elem%Forces(ngx-1,0:ngz-1,1) + Elem%Vhat(imin:imax,1)
        Elem%Forces(ngx-1,0:ngz-1,2) = Elem%Forces(ngx-1,0:ngz-1,2) + Vxz(imin:imax)
        Elem%Forces(ngx-1,0:ngz-1,3) = Elem%Forces(ngx-1,0:ngz-1,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(ngx-1,0:ngz-1,4) = Elem%Forces(ngx-1,0:ngz-1,4) + Elem%TracFace(imin:imax,1)
        ! For the Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Elem%Forces(0:ngx-1,ngz-1,0) = Elem%Forces(0:ngx-1,ngz-1,0) + Elem%Vhat(imin:imax,0)
        Elem%Forces(0:ngx-1,ngz-1,1) = Elem%Forces(0:ngx-1,ngz-1,1) + Elem%Vhat(imin:imax,1)
        Elem%Forces(0:ngx-1,ngz-1,2) = Elem%Forces(0:ngx-1,ngz-1,2) + Vxz(imin:imax)
        Elem%Forces(0:ngx-1,ngz-1,3) = Elem%Forces(0:ngx-1,ngz-1,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(0:ngx-1,ngz-1,4) = Elem%Forces(0:ngx-1,ngz-1,4) + Elem%TracFace(imin:imax,1)
        ! For the Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Elem%Forces(0,0:ngz-1,0) = Elem%Forces(0,0:ngz-1,0) + Elem%Vhat(imin:imax,0)
        Elem%Forces(0,0:ngz-1,1) = Elem%Forces(0,0:ngz-1,1) + Elem%Vhat(imin:imax,1)
        Elem%Forces(0,0:ngz-1,2) = Elem%Forces(0,0:ngz-1,2) + Vxz(imin:imax)
        Elem%Forces(0,0:ngz-1,3) = Elem%Forces(0,0:ngz-1,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(0,0:ngz-1,4) = Elem%Forces(0,0:ngz-1,4) + Elem%TracFace(imin:imax,1)

    end subroutine compute_Traces


! ###########################################################
!>
!!\brief Subroutine that computes the "imin" and "imax" indexes  which corresponds to the
!! begining and the end of the arrays that contain data for the local face "w_face" of element "Elem".
!!\version 1.0
!!\date 03/04/2014
!! This subroutine is used with DG and HDG elements.
!<
subroutine get_iminimax (Elem,w_face,imin,imax)

    implicit none

    type(element), intent(IN) :: Elem
    integer, intent (IN)      :: w_face
    integer, intent(INOUT)    :: imin
    integer, intent(INOUT)    :: imax

    select case (w_face)
    case(0)
        imin = 0
        imax = Elem%ngllx-1
    case(1)
        imin = Elem%ngllx
        imax = Elem%ngllx   + Elem%ngllz-1
    case(2)
        imin = Elem%ngllx   + Elem%ngllz
        imax = 2*Elem%ngllx + Elem%ngllz-1
    case(3)
        imin = 2*Elem%ngllx + Elem%ngllz
        imax = 2*Elem%ngllx + 2*Elem%ngllz -1
    end select

end subroutine get_iminimax

! ###########################################################

end module selement
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
