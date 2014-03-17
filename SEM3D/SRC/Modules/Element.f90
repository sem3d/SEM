!>
!!\file Element.f90
!!\brief contient les méthodes qui assure la gestion du type Element.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module selement
    use deriv3d
    implicit none

    type :: element_solid
       real, dimension(:,:,:,:), allocatable :: Cij

       real, dimension(:,:,:,:), allocatable :: ACoeff
       real, dimension(:,:,:,:), allocatable :: Forces,Veloc,Displ,Accel,V0

        ! Attenuation
       real, dimension (:,:,:), allocatable :: Q, Qs, Qp, onemSbeta, onemPbeta, &
           epsilonvol_, &
           epsilondev_xx_,epsilondev_yy_,epsilondev_xy_,epsilondev_xz_,epsilondev_yz_

       real, dimension(:,:,:,:), allocatable :: &
           factor_common_3, alphaval_3,betaval_3,gammaval_3, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_, &
           factor_common_P, alphaval_P,betaval_P,gammaval_P, R_vol_

    end type element_solid

    type :: element_fluid
       real, dimension(:,:,:,:), allocatable :: ACoeff
       ! fluid part
       real, dimension(:,:,:), allocatable:: Phi,VelPhi0,VelPhi,AccelPhi
       ! fluid part
       real, dimension(:,:,:), allocatable:: ForcesFl

    end type element_fluid

    type :: element_solid_pml
       ! TODO move pml related data here
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: Residual_Stress1, Residual_Stress2, Residual_Stress3
       real, dimension(:,:,:,:), allocatable :: DumpSx,DumpSy,DumpSz
       real, dimension(:,:,:,:), allocatable :: Forces1,Forces2,Forces3
       real, dimension(:,:,:,:), allocatable :: Veloc1,Veloc2,Veloc3
       ! FPML
       real, dimension(:,:,:), allocatable :: Isx, Isy, Isz
       real, dimension(:,:,:), allocatable :: Ivx, Ivy, Ivz
       real, dimension(:,:,:,:), allocatable :: Iveloc1, Iveloc2, Iveloc3
       real, dimension(:,:,:,:), allocatable :: I_Diagonal_Stress1, I_Diagonal_Stress2, I_Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: I_Residual_Stress1, I_Residual_Stress2
       real, dimension(:,:,:,:), allocatable :: DumpVx,DumpVy,DumpVz, DumpMass
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress, Residual_Stress
       real, dimension(:,:), allocatable :: Normales, Inv_Normales
    end type element_solid_pml

    type :: element_fluid_pml
        real, dimension(:,:,:,:), allocatable :: DumpSx,DumpSy,DumpSz
        real, dimension(:,:,:,:), allocatable :: DumpVx,DumpVy,DumpVz, DumpMass
        real, dimension(:,:,:,:), allocatable :: Veloc,Veloc1,Veloc2,Veloc3
        real, dimension(:,:,:), allocatable :: ForcesFl1,ForcesFl2,ForcesFl3,VelPhi1,VelPhi2,VelPhi3
    end type element_fluid_pml

    type :: element

       integer :: mat_index, ngllx, nglly, ngllz
       integer, dimension (:), allocatable :: Control_nodes
       integer, dimension (0:5) :: Near_Faces, Orient_Faces
       integer, dimension (0:11) :: Near_Edges, Orient_Edges
       integer, dimension (0:7) :: Near_Vertices

       integer, dimension (:,:,:), allocatable :: Iglobnum,Num
       real, dimension (:,:,:), allocatable :: Jacob
       real, dimension (:), allocatable :: wgtx, wgty, wgtz
       real, dimension(:,:,:,:,:), allocatable :: InvGrad
       real, dimension (:,:,:), allocatable :: Density, MassMat

       real, dimension (:,:,:), allocatable :: Lambda, Mu, Kappa
       ! flag for PML allocation
       logical :: PML, FPML
       logical  :: solid, fluid_dirich
       ! Whether this element will be part of snapshot outputs
       logical :: OUTPUT

       type(element_solid), allocatable :: sl
       type(element_fluid), allocatable :: fl
       type(element_solid_pml), allocatable :: slpml
       type(element_fluid_pml), allocatable :: flpml

       real :: dist_max !! taille caracteristique de l'element
    end type element


    interface
       subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
           CHARACTER*1        TRANSA, TRANSB
           INTEGER            M, N, K, LDA, LDB, LDC
           DOUBLE PRECISION   ALPHA, BETA
           DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
       end subroutine DGEMM
    end interface

contains


    !--------------------------------------------------------------
    !>
    !! \fn subroutine Prediction_Elem_Veloc (Elem,alpha,bega, gam1,dt)
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<
    subroutine Prediction_Elem_Veloc(Elem)

        implicit none

        type(Element), intent(inout) :: Elem
        integer :: ngllx, nglly, ngllz


        ngllx = Elem%ngllx ; nglly = Elem%nglly ; ngllz = Elem%ngllz

        Elem%sl%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%sl%Displ
        Elem%sl%V0 = Elem%sl%Veloc

        return
    end subroutine Prediction_Elem_Veloc
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    subroutine Prediction_Elem_VelPhi(Elem,dt)
        implicit none
        type(Element), intent(inout) :: Elem
        real, intent(in) :: dt
        integer :: ngllx, nglly, ngllz

        ngllx = Elem%ngllx ;  nglly = Elem%nglly ; ngllz = Elem%ngllz
        Elem%fl%VelPhi0(:,:,:) = Elem%fl%VelPhi(:,:,:)
        Elem%fl%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) = Elem%fl%Phi(:,:,:)
        return
    end subroutine Prediction_Elem_VelPhi
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !>
    !! \fn subroutine Correction_Elem_Veloc (Elem, bega, gam1,dt)
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) dt
    !<
    subroutine Correction_Elem_Veloc(Elem,dt)
        implicit none
        type(Element), intent(inout) :: Elem
        real, intent(in) :: dt
        integer :: ngllx, nglly, ngllz, i


        ngllx = Elem%ngllx ;  nglly = Elem%nglly ; ngllz = Elem%ngllz
        do i = 0,2
            Elem%sl%Forces(1:ngllx-2,1:nglly-2, 1:ngllz-2,i) = Elem%MassMat(:,:,:)*    &
                Elem%sl%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
        enddo
        Elem%sl%Veloc(:,:,:,:) = Elem%sl%v0(:,:,:,:) + dt * Elem%sl%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,:)
        Elem%sl%Accel  =  (Elem%sl%Veloc-Elem%sl%V0)/dt
        Elem%sl%Displ = Elem%sl%Displ + dt * Elem%sl%Veloc

        return
    end subroutine Correction_Elem_Veloc
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    subroutine Correction_Elem_VelPhi(Elem,dt)
        implicit none
        type(Element), intent(inout) :: Elem
        real, intent(in) :: dt
        integer :: ngllx, nglly, ngllz

        ngllx = Elem%ngllx ;  nglly = Elem%nglly ; ngllz = Elem%ngllz

        Elem%fl%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) = Elem%MassMat(:,:,:)*    &
            Elem%fl%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%fl%VelPhi(:,:,:) = Elem%fl%VelPhi0(:,:,:)+ dt * Elem%fl%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%fl%AccelPhi  = (Elem%fl%VelPhi-Elem%fl%VelPhi0)/dt
        Elem%fl%Phi = Elem%fl%Phi + dt*Elem%fl%VelPhi

        return
    end subroutine Correction_Elem_VelPhi
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !>
    !! \fn subroutine compute_InternalForces_Elem (Elem,hprime, hTprime, hprimez,hTprimez)
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<
    subroutine compute_InternalForces_Elem(Elem,hprimex,htprimex,hprimey,htprimey,hprimez,htprimez)

        implicit none

        type(Element), intent(inout) :: Elem
        real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hprimex,htprimex
        real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hprimey,htprimey
        real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hprimez,hTprimez

        integer :: n_z, m1, m2, m3

        real, dimension(0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1) :: dUx_dxi, dUx_deta, dUx_dzeta, &
            dUy_dxi, dUy_deta, dUy_dzeta, &
            dUz_dxi, dUz_deta, dUz_dzeta, &
            t1, s0, Uxloc, Uyloc, Uzloc

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz


        !- gradient at GLL points
        ! dUx_(dxi,deta,dzeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,0),dUx_dxi,dUx_deta,dUx_dzeta)

        ! dUy_(dxi,deta,dzeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,1),dUy_dxi,dUy_deta,dUy_dzeta)

        ! dUz_(dxi,deta,dzeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,2),dUz_dxi,dUz_deta,dUz_dzeta)


        !- Internal forces
        t1 = Elem%sl%Acoeff(:,:,:,0)*dUx_dxi + Elem%sl%Acoeff(:,:,:,1)*dUx_deta + Elem%sl%Acoeff(:,:,:,2)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,3)*dUy_dxi + Elem%sl%Acoeff(:,:,:,4)*dUy_deta + Elem%sl%Acoeff(:,:,:,5)*dUy_dzeta +  &
            Elem%sl%Acoeff(:,:,:,6)*dUz_dxi + Elem%sl%Acoeff(:,:,:,7)*dUz_deta + Elem%sl%Acoeff(:,:,:,8)*dUz_dzeta

        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uxloc,m1)

        t1 = Elem%sl%Acoeff(:,:,:,1)*dUx_dxi + Elem%sl%Acoeff(:,:,:,9)*dUx_deta + Elem%sl%Acoeff(:,:,:,10)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,11)*dUy_dxi + Elem%sl%Acoeff(:,:,:,12)*dUy_deta + Elem%sl%Acoeff(:,:,:,13)*dUy_dzeta + &
            Elem%sl%Acoeff(:,:,:,14)*dUz_dxi + Elem%sl%Acoeff(:,:,:,15)*dUz_deta + Elem%sl%Acoeff(:,:,:,16)*dUz_dzeta

        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Uxloc = s0 + Uxloc

        t1 = Elem%sl%Acoeff(:,:,:,2)*dUx_dxi + Elem%sl%Acoeff(:,:,:,10)*dUx_deta + Elem%sl%Acoeff(:,:,:,17)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,18)*dUy_dxi + Elem%sl%Acoeff(:,:,:,19)*dUy_deta + Elem%sl%Acoeff(:,:,:,20)*dUy_dzeta +&
            Elem%sl%Acoeff(:,:,:,21)*dUz_dxi + Elem%sl%ACoeff(:,:,:,22)*dUz_deta + Elem%sl%Acoeff(:,:,:,23)*dUz_dzeta

        call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Uxloc = s0 + Uxloc

        t1 = Elem%sl%Acoeff(:,:,:,3)*dUx_dxi + Elem%sl%Acoeff(:,:,:,11)*dUx_deta + Elem%sl%Acoeff(:,:,:,18)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,24)*dUy_dxi + Elem%sl%Acoeff(:,:,:,25)*dUy_deta + Elem%sl%Acoeff(:,:,:,26)*dUy_dzeta + &
            Elem%sl%Acoeff(:,:,:,27)*dUz_dxi + Elem%sl%Acoeff(:,:,:,28)*dUz_deta + Elem%sl%Acoeff(:,:,:,29)*dUz_dzeta

        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uyloc,m1)

        t1 = Elem%sl%Acoeff(:,:,:,4)*dUx_dxi + Elem%sl%Acoeff(:,:,:,12)*dUx_deta + Elem%sl%Acoeff(:,:,:,19)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,25)*dUy_dxi + Elem%sl%Acoeff(:,:,:,30)*dUy_deta + Elem%sl%Acoeff(:,:,:,31)*dUy_dzeta + &
            Elem%sl%Acoeff(:,:,:,32)*dUz_dxi + Elem%sl%Acoeff(:,:,:,33)*dUz_deta + Elem%sl%Acoeff(:,:,:,34)*dUz_dzeta

        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Uyloc = s0 + Uyloc

        t1 = Elem%sl%Acoeff(:,:,:,5)*dUx_dxi + Elem%sl%Acoeff(:,:,:,13)*dUx_deta + Elem%sl%Acoeff(:,:,:,20)* dUx_dzeta +&
            Elem%sl%Acoeff(:,:,:,26)*dUy_dxi + Elem%sl%Acoeff(:,:,:,31)*dUy_deta + Elem%sl%Acoeff(:,:,:,35)*dUy_dzeta + &
            Elem%sl%Acoeff(:,:,:,36)*dUz_dxi + Elem%sl%Acoeff(:,:,:,37)*dUz_deta + Elem%sl%Acoeff(:,:,:,38)*dUz_dzeta

        call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Uyloc = s0 + Uyloc

        t1 = Elem%sl%Acoeff(:,:,:,6)*dUx_dxi + Elem%sl%Acoeff(:,:,:,14)*dUx_deta + Elem%sl%Acoeff(:,:,:,21)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,27)*dUy_dxi + Elem%sl%Acoeff(:,:,:,32)*dUy_deta + Elem%sl%Acoeff(:,:,:,36)*dUy_dzeta + &
            Elem%sl%Acoeff(:,:,:,39)*dUz_dxi + Elem%sl%Acoeff(:,:,:,40)*dUz_deta + Elem%sl%Acoeff(:,:,:,41)*dUz_dzeta

        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uzloc,m1)

        t1 = Elem%sl%Acoeff(:,:,:,7)*dUx_dxi + Elem%sl%Acoeff(:,:,:,15)*dUx_deta + Elem%sl%Acoeff(:,:,:,22)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,28)*dUy_dxi + Elem%sl%Acoeff(:,:,:,33)*dUy_deta + Elem%sl%Acoeff(:,:,:,37)*dUy_dzeta + &
            Elem%sl%Acoeff(:,:,:,40)*dUz_dxi + Elem%sl%Acoeff(:,:,:,42)*dUz_deta + Elem%sl%Acoeff(:,:,:,43)*dUz_dzeta

        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Uzloc = s0 + Uzloc

        t1 = Elem%sl%Acoeff(:,:,:,8)*dUx_dxi + Elem%sl%Acoeff(:,:,:,16)*dUx_deta + Elem%sl%Acoeff(:,:,:,23)*dUx_dzeta + &
            Elem%sl%Acoeff(:,:,:,29)*dUy_dxi + Elem%sl%Acoeff(:,:,:,34)*dUy_deta + Elem%sl%Acoeff(:,:,:,38)*dUy_dzeta + &
            Elem%sl%Acoeff(:,:,:,41)*dUz_dxi + Elem%sl%Acoeff(:,:,:,43)*dUz_deta + Elem%sl%Acoeff(:,:,:,44)*dUz_dzeta

        call DGEMM ('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Uzloc = Uzloc + s0

        Elem%sl%Forces(:,:,:,0) = Uxloc
        Elem%sl%Forces(:,:,:,1) = Uyloc
        Elem%sl%Forces(:,:,:,2) = Uzloc

        return
    end subroutine compute_InternalForces_Elem
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    subroutine compute_InternalForces_Elem_Fluid(Elem,hprimex,htprimex,hprimey,htprimey,hprimez,htprimez)

        implicit none

        type(Element), intent(inout) :: Elem
        real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hprimex,htprimex
        real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hprimey,htprimey
        real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hprimez,hTprimez
        integer :: n_z,m1,m2,m3
        real, dimension(0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1) ::    &
            dPhi_dxi,dPhi_deta,dPhi_dzeta,t1,s0,Floc

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz


        !- gradients at GLLs points
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%fl%ForcesFl(:,:,:), &
            dPhi_dxi,dPhi_deta,dPhi_dzeta)

        !- Internal Forces
        t1 = Elem%fl%Acoeff(:,:,:,0)*dPhi_dxi + Elem%fl%Acoeff(:,:,:,1)*dPhi_deta + Elem%fl%Acoeff(:,:,:,2)*dPhi_dzeta
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Floc,m1)

        t1 = Elem%fl%Acoeff(:,:,:,1)*dPhi_dxi + Elem%fl%Acoeff(:,:,:,3)*dPhi_deta + Elem%fl%Acoeff(:,:,:,4)*dPhi_dzeta
        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Floc = s0 + Floc

        t1 = Elem%fl%Acoeff(:,:,:,2)*dPhi_dxi + Elem%fl%Acoeff(:,:,:,4)*dPhi_deta + Elem%fl%Acoeff(:,:,:,5)*dPhi_dzeta
        call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Floc = s0 + Floc

        !-
        Elem%fl%ForcesFl(:,:,:) = Floc

        return
    end subroutine compute_InternalForces_Elem_Fluid


    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !>
    !! \fn subroutine Prediction_Elem_PML_Veloc (Elem,alpha, bega, dt,Vxloc,Vzloc,Hmatz, HTmat)
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
    subroutine Prediction_Elem_PML_Veloc(Elem,bega,dt,hTprimex,Hprimey,Hprimez,rg,n)

        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) :: bega, dt
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        integer, intent (IN) :: rg, n

        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVx_dxi,dVx_deta,dVx_dzeta, &
            dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta

        integer :: m1, m2,m3


        m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz
  ! useless as bega = 1/2 in general
        Elem%sl%Forces(1:m1-2,1:m2-2, 1:m3-2,0:2)  = Elem%sl%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%sl%Accel(:,:,:,:)

        ! partial of velocity components with respect to xi,eta,zeta
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,0),dVx_dxi,dVx_deta,dVx_dzeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,1),dVy_dxi,dVy_deta,dVy_dzeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,2),dVz_dxi,dVz_deta,dVz_dzeta)


  ! Stress_xx
   ! (Stress_xx)^x
        Elem%slpml%Diagonal_Stress1 (:,:,:,0) = Elem%slpml%DumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,0) + &
            Elem%slpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,0) * dVx_dxi + Elem%sl%Acoeff(:,:,:,1) * dVx_deta + Elem%sl%Acoeff(:,:,:,2) * dVx_dzeta)
   ! (Stress_xx)^y
        Elem%slpml%Diagonal_Stress2 (:,:,:,0) = Elem%slpml%DumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,0) + &
            Elem%slpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,3) * dVy_dxi + Elem%sl%Acoeff(:,:,:,4) * dVy_deta + Elem%sl%Acoeff(:,:,:,5) * dVy_dzeta)
   ! (Stress_xx)^z
        Elem%slpml%Diagonal_Stress3 (:,:,:,0) = Elem%slpml%DumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,0) + &
            Elem%slpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,6) * dVz_dxi + Elem%sl%Acoeff(:,:,:,7) * dVz_deta + Elem%sl%Acoeff(:,:,:,8) * dVz_dzeta)

  ! Stress_yy
        Elem%slpml%Diagonal_Stress1 (:,:,:,1) = Elem%slpml%DumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,1) + &
            Elem%slpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,9) * dVx_dxi + Elem%sl%Acoeff(:,:,:,10) * dVx_deta + Elem%sl%Acoeff(:,:,:,11) * dVx_dzeta)
        Elem%slpml%Diagonal_Stress2 (:,:,:,1) = Elem%slpml%DumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,1) + &
            Elem%slpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,12) * dVy_dxi + Elem%sl%Acoeff(:,:,:,13) * dVy_deta + Elem%sl%Acoeff(:,:,:,14) * dVy_dzeta)
        Elem%slpml%Diagonal_Stress3 (:,:,:,1) = Elem%slpml%DumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,1) + &
            Elem%slpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,6) * dVz_dxi + Elem%sl%Acoeff(:,:,:,7) * dVz_deta + Elem%sl%Acoeff(:,:,:,8) * dVz_dzeta)

  ! Stress_zz
        Elem%slpml%Diagonal_Stress1 (:,:,:,2) = Elem%slpml%DumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,2) + &
            Elem%slpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,9) * dVx_dxi + Elem%sl%Acoeff(:,:,:,10) * dVx_deta + Elem%sl%Acoeff(:,:,:,11) * dVx_dzeta)
        Elem%slpml%Diagonal_Stress2 (:,:,:,2) = Elem%slpml%DumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,2) + &
            Elem%slpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,3) * dVy_dxi + Elem%sl%Acoeff(:,:,:,4) * dVy_deta + Elem%sl%Acoeff(:,:,:,5) * dVy_dzeta)
        Elem%slpml%Diagonal_Stress3 (:,:,:,2) = Elem%slpml%DumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,2) + &
            Elem%slpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,15) * dVz_dxi + Elem%sl%Acoeff(:,:,:,16) * dVz_deta + Elem%sl%Acoeff(:,:,:,17) * dVz_dzeta)

        Elem%slpml%Diagonal_Stress = Elem%slpml%Diagonal_Stress1 + Elem%slpml%Diagonal_Stress2 + Elem%slpml%Diagonal_Stress3

  ! Stress_xy
        Elem%slpml%Residual_Stress1 (:,:,:,0) = Elem%slpml%DumpSx(:,:,:,0) * Elem%slpml%Residual_Stress1 (:,:,:,0) + &
            Elem%slpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,18) * dVy_dxi + Elem%sl%Acoeff(:,:,:,19) * dVy_deta + Elem%sl%Acoeff(:,:,:,20) * dVy_dzeta)
        Elem%slpml%Residual_Stress2 (:,:,:,0) = Elem%slpml%DumpSy(:,:,:,0) * Elem%slpml%Residual_Stress2 (:,:,:,0) + &
            Elem%slpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,21) * dVx_dxi + Elem%sl%Acoeff(:,:,:,22) * dVx_deta + Elem%sl%Acoeff(:,:,:,23) * dVx_dzeta)
        Elem%slpml%Residual_Stress3 (:,:,:,0) = Elem%slpml%DumpSz(:,:,:,0) * Elem%slpml%Residual_Stress3 (:,:,:,0)

  ! Stress_xz
        Elem%slpml%Residual_Stress1 (:,:,:,1) = Elem%slpml%DumpSx(:,:,:,0) * Elem%slpml%Residual_Stress1 (:,:,:,1) + &
            Elem%slpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,18) * dVz_dxi + Elem%sl%Acoeff(:,:,:,19) * dVz_deta + Elem%sl%Acoeff(:,:,:,20) * dVz_dzeta)
        Elem%slpml%Residual_Stress2 (:,:,:,1) = Elem%slpml%DumpSy(:,:,:,0) * Elem%slpml%Residual_Stress2 (:,:,:,1)
        Elem%slpml%Residual_Stress3 (:,:,:,1) = Elem%slpml%DumpSz(:,:,:,0) * Elem%slpml%Residual_Stress3 (:,:,:,1) + &
            Elem%slpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,24) * dVx_dxi + Elem%sl%Acoeff(:,:,:,25) * dVx_deta + Elem%sl%Acoeff(:,:,:,26) * dVx_dzeta)

  ! Stress_yz
        Elem%slpml%Residual_Stress1 (:,:,:,2) = Elem%slpml%DumpSx(:,:,:,0) * Elem%slpml%Residual_Stress1 (:,:,:,2)
        Elem%slpml%Residual_Stress2 (:,:,:,2) = Elem%slpml%DumpSy(:,:,:,0) * Elem%slpml%Residual_Stress2 (:,:,:,2) + &
            Elem%slpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,21) * dVz_dxi + Elem%sl%Acoeff(:,:,:,22) * dVz_deta + Elem%sl%Acoeff(:,:,:,23) * dVz_dzeta)
        Elem%slpml%Residual_Stress3 (:,:,:,2) = Elem%slpml%DumpSz(:,:,:,0) * Elem%slpml%Residual_Stress3 (:,:,:,2) + &
            Elem%slpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,24) * dVy_dxi + Elem%sl%Acoeff(:,:,:,25) * dVy_deta + Elem%sl%Acoeff(:,:,:,26) * dVy_dzeta)

        Elem%slpml%Residual_Stress = Elem%slpml%Residual_Stress1 + Elem%slpml%Residual_Stress2 + Elem%slpml%Residual_Stress3

        return
    end subroutine Prediction_Elem_PML_Veloc
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    subroutine Prediction_Elem_PML_VelPhi(Elem,bega,dt,hTprimex,Hprimey,Hprimez)
        ! same as previously, but for fluid part
        implicit none

        type(Element), intent(inout) :: Elem
        real, intent(in) :: bega, dt
        real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hTprimex
        real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hprimey
        real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hprimez

        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVelPhi_dxi,dVelPhi_deta,dVelPhi_dzeta
        integer :: m1,m2,m3


        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz


        ! prediction in the element
        Elem%fl%ForcesFl(1:m1-2,1:m2-2,1:m3-2) = Elem%fl%VelPhi(:,:,:) +dt*(0.5-bega)*Elem%fl%AccelPhi(:,:,:)

        ! d(rho*Phi)_d(xi,eta,zeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%fl%ForcesFl(:,:,:), &
            dVelPhi_dxi,dVelPhi_deta,dVelPhi_dzeta)

        ! prediction for (physical) velocity (which is the equivalent of a stress, here)
        ! V_x^x
        Elem%flpml%Veloc1(:,:,:,0) = Elem%flpml%dumpSx(:,:,:,0) * Elem%flpml%Veloc1(:,:,:,0) + Elem%flpml%dumpSx(:,:,:,1) * Dt *  &
            (Elem%fl%Acoeff(:,:,:,0) * dVelPhi_dxi + Elem%fl%Acoeff(:,:,:,1) * dVelPhi_deta + Elem%fl%Acoeff(:,:,:,2) * dVelPhi_dzeta)
        ! V_x^y
        Elem%flpml%Veloc2(:,:,:,0) = Elem%flpml%dumpSy(:,:,:,0) * Elem%flpml%Veloc2(:,:,:,0)
        ! V_x^z
        Elem%flpml%Veloc3(:,:,:,0) = Elem%flpml%dumpSz(:,:,:,0) * Elem%flpml%Veloc3(:,:,:,0)
        ! V_y^x
        Elem%flpml%Veloc1(:,:,:,1) = Elem%flpml%dumpSx(:,:,:,0) * Elem%flpml%Veloc1(:,:,:,1)
        ! V_y^y
        Elem%flpml%Veloc2(:,:,:,1) = Elem%flpml%dumpSy(:,:,:,0) * Elem%flpml%Veloc2(:,:,:,1) + Elem%flpml%dumpSy(:,:,:,1) * Dt *  &
            (Elem%fl%Acoeff(:,:,:,3) * dVelPhi_dxi + Elem%fl%Acoeff(:,:,:,4) * dVelPhi_deta + Elem%fl%Acoeff(:,:,:,5) * dVelPhi_dzeta)
        ! V_y^z
        Elem%flpml%Veloc3(:,:,:,1) = Elem%flpml%dumpSz(:,:,:,0) * Elem%flpml%Veloc3(:,:,:,1)
        ! V_z^x
        Elem%flpml%Veloc1(:,:,:,2) = Elem%flpml%dumpSx(:,:,:,0) * Elem%flpml%Veloc1(:,:,:,2)
        ! V_z^y
        Elem%flpml%Veloc2(:,:,:,2) = Elem%flpml%dumpSy(:,:,:,0) * Elem%flpml%Veloc2(:,:,:,2)
        ! V_z^z
        Elem%flpml%Veloc3(:,:,:,2) = Elem%flpml%dumpSz(:,:,:,0) * Elem%flpml%Veloc3(:,:,:,2) + Elem%flpml%dumpSz(:,:,:,1) * Dt *  &
            (Elem%fl%Acoeff(:,:,:,6) * dVelPhi_dxi + Elem%fl%Acoeff(:,:,:,7) * dVelPhi_deta + Elem%fl%Acoeff(:,:,:,8) * dVelPhi_dzeta)

        ! total velocity vector after dumping = sum of splitted parts
        Elem%flpml%Veloc(:,:,:,:) = Elem%flpml%Veloc1(:,:,:,:) + Elem%flpml%Veloc2(:,:,:,:) + Elem%flpml%Veloc3(:,:,:,:)

        return
    end subroutine Prediction_Elem_PML_VelPhi
    !------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------
    !>
    !! \fn subroutine Prediction_Elem_FPML_Veloc (Elem,alpha, bega, dt,Vxloc,Vzloc,Hmatz, HTmat,fil)
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
    subroutine Prediction_Elem_FPML_Veloc(Elem,bega,dt,hTprimex,Hprimey,Hprimez,rg,n,fil)

        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) :: bega, dt, fil
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        integer, intent (IN) :: rg, n

        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVx_dxi,dVx_deta,dVx_dzeta, &
            dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta, Stress_ausiliar

        integer :: m1, m2,m3 ,n_z
        real :: fil2


        m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz
        fil2 = fil**2

        Elem%sl%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2) = Elem%sl%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%sl%Accel(:,:,:,:)

        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,0),dVx_dxi,dVx_deta,dVx_dzeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,1),dVy_dxi,dVy_deta,dVy_dzeta)
        call elem_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%sl%Forces(:,:,:,2),dVz_dxi,dVz_deta,dVz_dzeta)

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress1 (:,:,:,0)
        Elem%slpml%Diagonal_Stress1 (:,:,:,0) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,0) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,0) * dVx_dxi + Elem%sl%Acoeff(:,:,:,1) * dVx_deta + &
            Elem%sl%Acoeff(:,:,:,2) * dVx_dzeta) + Elem%slpml%Isx * Elem%slpml%I_Diagonal_stress1 (:,:,:,0)
        Elem%slpml%I_Diagonal_Stress1  (:,:,:,0)= Fil2* Elem%slpml%I_Diagonal_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress1 (:,:,:,0) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress2 (:,:,:,0)
        Elem%slpml%Diagonal_Stress2 (:,:,:,0) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,0) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,3) * dVy_dxi + Elem%sl%Acoeff(:,:,:,4) * dVy_deta + &
            Elem%sl%Acoeff(:,:,:,5) * dVy_dzeta) + Elem%slpml%Isy * Elem%slpml%I_Diagonal_stress2 (:,:,:,0)
        Elem%slpml%I_Diagonal_Stress2  (:,:,:,0)= Fil2* Elem%slpml%I_Diagonal_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress2 (:,:,:,0) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress3 (:,:,:,0)
        Elem%slpml%Diagonal_Stress3 (:,:,:,0) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,0) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,6) * dVz_dxi + Elem%sl%Acoeff(:,:,:,7) * dVz_deta + &
            Elem%sl%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%slpml%Isz * Elem%slpml%I_Diagonal_stress3 (:,:,:,0)
        Elem%slpml%I_Diagonal_Stress3  (:,:,:,0)= Fil2* Elem%slpml%I_Diagonal_Stress3(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress3 (:,:,:,0) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress1 (:,:,:,1)
        Elem%slpml%Diagonal_Stress1 (:,:,:,1) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,1) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,9) * dVx_dxi + Elem%sl%Acoeff(:,:,:,10) * dVx_deta + Elem%sl%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%slpml%Isx * Elem%slpml%I_Diagonal_stress1 (:,:,:,1)
        Elem%slpml%I_Diagonal_Stress1  (:,:,:,1)= Fil2* Elem%slpml%I_Diagonal_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress1 (:,:,:,1) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress2 (:,:,:,1)
        Elem%slpml%Diagonal_Stress2 (:,:,:,1) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,1) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,12) * dVy_dxi + Elem%sl%Acoeff(:,:,:,13) * dVy_deta + Elem%sl%Acoeff(:,:,:,14) * dVy_dzeta) + Elem%slpml%Isy * Elem%slpml%I_Diagonal_stress2 (:,:,:,1)
        Elem%slpml%I_Diagonal_Stress2  (:,:,:,1)= Fil2* Elem%slpml%I_Diagonal_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress2 (:,:,:,1) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress3 (:,:,:,1)
        Elem%slpml%Diagonal_Stress3 (:,:,:,1) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,1) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,6) * dVz_dxi + Elem%sl%Acoeff(:,:,:,7) * dVz_deta + Elem%sl%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%slpml%Isz * Elem%slpml%I_Diagonal_stress3 (:,:,:,1)
        Elem%slpml%I_Diagonal_Stress3  (:,:,:,1)= Fil2* Elem%slpml%I_Diagonal_Stress3(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress3 (:,:,:,1) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress1 (:,:,:,2)
        Elem%slpml%Diagonal_Stress1 (:,:,:,2) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,2) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,9) * dVx_dxi + Elem%sl%Acoeff(:,:,:,10) * dVx_deta + Elem%sl%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%slpml%Isx * Elem%slpml%I_Diagonal_stress1 (:,:,:,2)
        Elem%slpml%I_Diagonal_Stress1  (:,:,:,2)= Fil2* Elem%slpml%I_Diagonal_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress1 (:,:,:,2) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress2 (:,:,:,2)
        Elem%slpml%Diagonal_Stress2 (:,:,:,2) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,2) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,3) * dVy_dxi + Elem%sl%Acoeff(:,:,:,4) * dVy_deta + Elem%sl%Acoeff(:,:,:,5) * dVy_dzeta)+ Elem%slpml%Isy * Elem%slpml%I_Diagonal_stress2 (:,:,:,2)
        Elem%slpml%I_Diagonal_Stress2  (:,:,:,2)= Fil2* Elem%slpml%I_Diagonal_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress2 (:,:,:,2) )

        Stress_Ausiliar = Elem%slpml%Diagonal_Stress3 (:,:,:,2)
        Elem%slpml%Diagonal_Stress3 (:,:,:,2) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,2) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,15) * dVz_dxi + Elem%sl%Acoeff(:,:,:,16) * dVz_deta + Elem%sl%Acoeff(:,:,:,17) * dVz_dzeta) + Elem%slpml%Isz * Elem%slpml%I_Diagonal_stress3 (:,:,:,2)
        Elem%slpml%I_Diagonal_Stress3  (:,:,:,2)= Fil2* Elem%slpml%I_Diagonal_Stress3(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%Diagonal_Stress3 (:,:,:,2) )


        Elem%slpml%Diagonal_Stress = Elem%slpml%Diagonal_Stress1 + Elem%slpml%Diagonal_Stress2 + Elem%slpml%Diagonal_Stress3


        Stress_Ausiliar = Elem%slpml%residual_Stress1 (:,:,:,0)
        Elem%slpml%residual_Stress1 (:,:,:,0) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%residual_Stress1 (:,:,:,0) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,18) * dVy_dxi + Elem%sl%Acoeff(:,:,:,19) * dVy_deta + Elem%sl%Acoeff(:,:,:,20) * dVy_dzeta) + Elem%slpml%Isx * Elem%slpml%I_Residual_stress1(:,:,:,0)
        Elem%slpml%I_Residual_Stress1  (:,:,:,0)= Fil2* Elem%slpml%I_Residual_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%residual_Stress1 (:,:,:,0) )

        Stress_Ausiliar = Elem%slpml%residual_Stress2 (:,:,:,0)
        Elem%slpml%residual_Stress2 (:,:,:,0) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%residual_Stress2 (:,:,:,0) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,21) * dVx_dxi + Elem%sl%Acoeff(:,:,:,22) * dVx_deta + Elem%sl%Acoeff(:,:,:,23) * dVx_dzeta) + Elem%slpml%Isy * Elem%slpml%I_Residual_stress2(:,:,:,0)
        Elem%slpml%I_Residual_Stress2  (:,:,:,0) = Fil2* Elem%slpml%I_Residual_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%residual_Stress2 (:,:,:,0) )

        Stress_Ausiliar = Elem%slpml%residual_Stress1 (:,:,:,1)
        Elem%slpml%residual_Stress1 (:,:,:,1) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%residual_Stress1 (:,:,:,1) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,18) * dVz_dxi + Elem%sl%Acoeff(:,:,:,19) * dVz_deta + Elem%sl%Acoeff(:,:,:,20) * dVz_dzeta) + Elem%slpml%Isx * Elem%slpml%I_Residual_stress1(:,:,:,1)
        Elem%slpml%I_Residual_Stress1  (:,:,:,1)= Fil2* Elem%slpml%I_Residual_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%residual_Stress1 (:,:,:,1) )

        Stress_Ausiliar = Elem%slpml%residual_Stress2 (:,:,:,1)
        Elem%slpml%residual_Stress2 (:,:,:,1) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%residual_Stress2 (:,:,:,1) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,24) * dVx_dxi + Elem%sl%Acoeff(:,:,:,25) * dVx_deta + Elem%sl%Acoeff(:,:,:,26) * dVx_dzeta) + Elem%slpml%Isz * Elem%slpml%I_Residual_stress2(:,:,:,1)
        Elem%slpml%I_Residual_Stress2  (:,:,:,1) = Fil2* Elem%slpml%I_Residual_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%residual_Stress2 (:,:,:,1) )

        Stress_Ausiliar = Elem%slpml%residual_Stress1 (:,:,:,2)
        Elem%slpml%residual_Stress1 (:,:,:,2) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%residual_Stress1 (:,:,:,2) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,21) * dVz_dxi + Elem%sl%Acoeff(:,:,:,22) * dVz_deta + Elem%sl%Acoeff(:,:,:,23) * dVz_dzeta) + Elem%slpml%Isy * Elem%slpml%I_Residual_stress1(:,:,:,2)
        Elem%slpml%I_Residual_Stress1  (:,:,:,2)= Fil2* Elem%slpml%I_Residual_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%residual_Stress1 (:,:,:,2) )

        Stress_Ausiliar = Elem%slpml%residual_Stress2 (:,:,:,2)
        Elem%slpml%residual_Stress2 (:,:,:,2) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%residual_Stress2 (:,:,:,2) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,24) * dVy_dxi + Elem%sl%Acoeff(:,:,:,25) * dVy_deta + Elem%sl%Acoeff(:,:,:,26) * dVy_dzeta)+ Elem%slpml%Isz * Elem%slpml%I_Residual_stress2(:,:,:,2)
        Elem%slpml%I_Residual_Stress2  (:,:,:,2) = Fil2* Elem%slpml%I_Residual_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%slpml%residual_Stress2 (:,:,:,2) )

        Elem%slpml%residual_Stress = Elem%slpml%residual_Stress1 + Elem%slpml%residual_Stress2

        return
    end subroutine Prediction_Elem_FPML_Veloc
    !------------------------------------------------------------------
    !------------------------------------------------------------------
    !>
    !! \fn subroutine Correction_Elem_PML_Veloc (Elem, dt)
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) dt
    !<
    subroutine Correction_Elem_PML_Veloc(Elem,dt)

        implicit none

        type(Element), intent(inout) :: Elem
        real, intent(in) :: dt
        integer :: ngllx, nglly, ngllz, i

        ngllx = Elem%ngllx; nglly = Elem%nglly; ngllz = Elem%ngllz


        do i = 0,2
            Elem%slpml%Veloc1(:,:,:,i) = Elem%slpml%dumpVx(:,:,:,0) * Elem%slpml%Veloc1(:,:,:,i) +    &
                dt * Elem%slpml%dumpVx(:,:,:,1)*Elem%slpml%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
            Elem%slpml%Veloc2(:,:,:,i) = Elem%slpml%dumpVy(:,:,:,0) * Elem%slpml%Veloc2(:,:,:,i) +    &
                dt * Elem%slpml%dumpVy(:,:,:,1)*Elem%slpml%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
            Elem%slpml%Veloc3(:,:,:,i) = Elem%slpml%dumpVz(:,:,:,0) * Elem%slpml%Veloc3(:,:,:,i) +    &
                dt * Elem%slpml%dumpVz(:,:,:,1)*Elem%slpml%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
        enddo

        Elem%sl%Veloc = Elem%slpml%Veloc1 + Elem%slpml%Veloc2 + Elem%slpml%Veloc3
        ! Usefull only for traces and debug
        Elem%sl%Displ = Elem%sl%Displ + dt * Elem%sl%Veloc

        return
    end subroutine Correction_Elem_PML_Veloc

    !------------------------------------------------------------------
    !------------------------------------------------------------------
    subroutine Correction_Elem_PML_VelPhi(Elem,dt)

        implicit none

        type(Element), intent(inout) :: Elem
        real, intent(in) :: dt
        integer :: ngllx, nglly, ngllz


        ngllx = Elem%ngllx ; nglly = Elem%nglly ; ngllz = Elem%ngllz

        Elem%flpml%VelPhi1(:,:,:) = Elem%slpml%dumpVx(:,:,:,0) * Elem%flpml%VelPhi1(:,:,:) +    &
            dt * Elem%slpml%dumpVx(:,:,:,1)*Elem%flpml%ForcesFl1(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%flpml%VelPhi2(:,:,:) = Elem%slpml%dumpVy(:,:,:,0) * Elem%flpml%VelPhi2(:,:,:) +    &
            dt * Elem%slpml%dumpVy(:,:,:,1)*Elem%flpml%ForcesFl2(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%flpml%VelPhi3(:,:,:) = Elem%slpml%dumpVz(:,:,:,0) * Elem%flpml%VelPhi3(:,:,:) +    &
            dt * Elem%slpml%dumpVz(:,:,:,1)*Elem%flpml%ForcesFl3(1:ngllx-2,1:nglly-2,1:ngllz-2)

        Elem%fl%VelPhi = Elem%flpml%VelPhi1 + Elem%flpml%VelPhi2 + Elem%flpml%VelPhi3

        Elem%fl%Phi = Elem%fl%Phi + dt * Elem%fl%VelPhi

        return
    end subroutine Correction_Elem_PML_VelPhi
    !------------------------------------------------------------------
    !------------------------------------------------------------------
    !>
    !! \fn subroutine Correction_Elem_FPML_Veloc (Elem, dt, fil)
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) fil
    !<
    subroutine Correction_Elem_FPML_Veloc (Elem, dt, fil)

        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) :: dt, fil

        integer :: ngllx, nglly, ngllz, i
        real :: fil2
        real, dimension (1:Elem%ngllx-2,1:Elem%nglly-2,1:Elem%ngllz-2) :: Ausiliar_velocity

        ngllx = Elem%ngllx; nglly = Elem%nglly; ngllz = Elem%ngllz
        fil2 = fil**2

        do i = 0,2
            Ausiliar_Velocity =  Elem%slpml%Veloc1(:,:,:,i)
            Elem%slpml%Veloc1(:,:,:,i) = Elem%slpml%dumpVx(:,:,:,0) * Elem%slpml%Veloc1(:,:,:,i) + dt * Elem%slpml%dumpVx(:,:,:,1)*Elem%slpml%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%slpml%Ivx * Elem%slpml%Iveloc1 (:,:,:,i)
            Elem%slpml%Iveloc1 (:,:,:,i)  = Fil2*Elem%slpml%Iveloc1 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%slpml%Veloc1(:,:,:,i))

            Ausiliar_Velocity =  Elem%slpml%Veloc2(:,:,:,i)
            Elem%slpml%Veloc2(:,:,:,i) = Elem%slpml%dumpVy(:,:,:,0) * Elem%slpml%Veloc2(:,:,:,i) + dt * Elem%slpml%dumpVy(:,:,:,1)*Elem%slpml%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%slpml%Ivy * Elem%slpml%Iveloc2 (:,:,:,i)
            Elem%slpml%Iveloc2 (:,:,:,i)  = Fil2*Elem%slpml%Iveloc2 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%slpml%Veloc2(:,:,:,i))

            Ausiliar_Velocity =  Elem%slpml%Veloc3(:,:,:,i)
            Elem%slpml%Veloc3(:,:,:,i) = Elem%slpml%dumpVz(:,:,:,0) * Elem%slpml%Veloc3(:,:,:,i) + dt * Elem%slpml%dumpVz(:,:,:,1)*Elem%slpml%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%slpml%Ivz * Elem%slpml%Iveloc3 (:,:,:,i)
            Elem%slpml%Iveloc3 (:,:,:,i)  = Fil2*Elem%slpml%Iveloc3 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%slpml%Veloc3(:,:,:,i))
        enddo

        Elem%sl%Veloc = Elem%slpml%Veloc1 + Elem%slpml%Veloc2 + Elem%slpml%Veloc3

        return
    end subroutine Correction_Elem_FPML_Veloc

    !------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------
    !>
    !! \fn subroutine compute_InternalForces_PML_Elem (Elem,hprime, hTprimez)
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<
    subroutine  compute_InternalForces_PML_Elem (Elem, hprimex, hTprimey, htprimez)

        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hTprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez

        integer :: m1, m2, m3, n_z
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0,s1


        m1 = Elem%ngllx;  m2 = Elem%nglly;  m3 = Elem%ngllz

        s0 = Elem%sl%Acoeff(:,:,:,27) * Elem%slpml%Diagonal_Stress(:,:,:,0) + Elem%sl%Acoeff(:,:,:,28) * Elem%slpml%residual_Stress(:,:,:,0) + &
            Elem%sl%Acoeff(:,:,:,29) * Elem%slpml%residual_Stress(:,:,:,1)
        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1, s0(:,:,:), m1, 0., s1, m1 )
        Elem%slpml%Forces1(:,:,:,0) = s1

        s0 = Elem%sl%Acoeff(:,:,:,30) * Elem%slpml%Diagonal_Stress(:,:,:,0) + Elem%sl%Acoeff(:,:,:,31) * Elem%slpml%residual_Stress(:,:,:,0) + &
            Elem%sl%Acoeff(:,:,:,32) * Elem%slpml%residual_Stress(:,:,:,1)
        do n_z = 0,m3-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
        enddo
        Elem%slpml%Forces2(:,:,:,0) = s1

        s0 = Elem%sl%Acoeff(:,:,:,33) * Elem%slpml%Diagonal_Stress(:,:,:,0) + Elem%sl%Acoeff(:,:,:,34) * Elem%slpml%residual_Stress(:,:,:,0) + &
            Elem%sl%Acoeff(:,:,:,35) * Elem%slpml%residual_Stress(:,:,:,1)

        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(:,:,:), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
        Elem%slpml%Forces3(:,:,:,0) = s1

        s0 = Elem%sl%Acoeff(:,:,:,27) * Elem%slpml%residual_Stress(:,:,:,0) + Elem%sl%Acoeff(:,:,:,28) * Elem%slpml%Diagonal_Stress(:,:,:,1) + &
            Elem%sl%Acoeff(:,:,:,29) * Elem%slpml%residual_Stress(:,:,:,2)
        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(:,:,:) ,m1, 0., s1, m1 )
        Elem%slpml%Forces1(:,:,:,1) = s1

        s0 = Elem%sl%Acoeff(:,:,:,30) * Elem%slpml%residual_Stress(:,:,:,0) + Elem%sl%Acoeff(:,:,:,31) * Elem%slpml%Diagonal_Stress(:,:,:,1) + &
            Elem%sl%Acoeff(:,:,:,32) * Elem%slpml%residual_Stress(:,:,:,2)
        do n_z = 0,m3-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
        enddo
        Elem%slpml%Forces2(:,:,:,1) = s1

        s0 = Elem%sl%Acoeff(:,:,:,33) * Elem%slpml%residual_Stress(:,:,:,0) + Elem%sl%Acoeff(:,:,:,34) * Elem%slpml%Diagonal_Stress(:,:,:,1) + &
            Elem%sl%Acoeff(:,:,:,35) * Elem%slpml%residual_Stress(:,:,:,2)

        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(:,:,:), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
        Elem%slpml%Forces3(:,:,:,1) = s1


        s0 = Elem%sl%Acoeff(:,:,:,27) * Elem%slpml%residual_Stress(:,:,:,1) + Elem%sl%Acoeff(:,:,:,28) * Elem%slpml%residual_Stress(:,:,:,2) + &
            Elem%sl%Acoeff(:,:,:,29) * Elem%slpml%Diagonal_Stress(:,:,:,2)
        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(:,:,:) ,m1, 0., s1, m1 )
        Elem%slpml%Forces1(:,:,:,2) = s1

        s0 = Elem%sl%Acoeff(:,:,:,30) * Elem%slpml%residual_Stress(:,:,:,1) + Elem%sl%Acoeff(:,:,:,31) * Elem%slpml%residual_Stress(:,:,:,2) + &
            Elem%sl%Acoeff(:,:,:,32) * Elem%slpml%Diagonal_Stress(:,:,:,2)
        do n_z = 0,m3-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
        enddo
        Elem%slpml%Forces2(:,:,:,2) = s1

        s0 = Elem%sl%Acoeff(:,:,:,33) * Elem%slpml%residual_Stress(:,:,:,1) + Elem%sl%Acoeff(:,:,:,34) * Elem%slpml%residual_Stress(:,:,:,2) + &
            Elem%sl%Acoeff(:,:,:,35) * Elem%slpml%Diagonal_Stress(:,:,:,2)

        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(:,:,:), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
        Elem%slpml%Forces3(:,:,:,2) = s1

        Elem%sl%Forces = Elem%slpml%Forces1 + Elem%slpml%Forces2 + Elem%slpml%Forces3

        return
    end subroutine compute_InternalForces_PML_Elem
    !------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------
    subroutine compute_InternalForces_PML_Elem_Fl(Elem,hprimex,hTprimey,htprimez)

        implicit none

        type(Element), intent(inout) :: Elem
        real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hprimex
        real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hTprimey
        real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hTprimez

        integer :: m1, m2, m3, n_z
        real, dimension(0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1)  :: s0,s1


        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        Elem%flpml%ForcesFl1(:,:,:) = 0d0
        Elem%flpml%ForcesFl2(:,:,:) = 0d0
        Elem%flpml%ForcesFl3(:,:,:) = 0d0
        s0 = 0d0 ; s1 = 0d0


        ! forces associated to V_x
        s0 = Elem%fl%Acoeff(:,:,:,9) * Elem%flpml%Veloc(:,:,:,0)
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(:,:,:),m1,0.,s1,m1)
        Elem%flpml%ForcesFl1(:,:,:) = s1

        s0 = Elem%fl%Acoeff(:,:,:,10) * Elem%flpml%Veloc(:,:,:,0)
        do n_z = 0,m3-1
            call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
        enddo
        Elem%flpml%ForcesFl1(:,:,:) = s1+Elem%flpml%ForcesFl1(:,:,:)

        s0 = Elem%fl%Acoeff(:,:,:,11) * Elem%flpml%Veloc(:,:,:,0)
        call DGEMM('N','N',m1*m2,m3,m3,1.,s0(:,:,:),m1*m2,htprimez,m3,0.,s1,m1*m2)
        Elem%flpml%ForcesFl1(:,:,:) = s1+Elem%flpml%ForcesFl1(:,:,:)

        ! forces associated to V_y
        s0 = Elem%fl%Acoeff(:,:,:,12) * Elem%flpml%Veloc(:,:,:,1)
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(:,:,:),m1,0.,s1,m1)
        Elem%flpml%ForcesFl2(:,:,:) = s1

        s0 = Elem%fl%Acoeff(:,:,:,13) * Elem%flpml%Veloc(:,:,:,1)
        do n_z = 0,m3-1
            call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
        enddo
        Elem%flpml%ForcesFl2(:,:,:) = s1+Elem%flpml%ForcesFl2(:,:,:)

        s0 = Elem%fl%Acoeff(:,:,:,14) * Elem%flpml%Veloc(:,:,:,1)
        call DGEMM('N','N',m1*m2,m3,m3,1.,s0(:,:,:),m1*m2,htprimez,m3,0.,s1,m1*m2)
        Elem%flpml%ForcesFl2(:,:,:) = s1+Elem%flpml%ForcesFl2(:,:,:)

        ! forces associated to V_z
        s0 = Elem%fl%Acoeff(:,:,:,15) * Elem%flpml%Veloc(:,:,:,2)
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(:,:,:),m1,0.,s1,m1)
        Elem%flpml%ForcesFl3(:,:,:) = s1

        s0 = Elem%fl%Acoeff(:,:,:,16) * Elem%flpml%Veloc(:,:,:,2)
        do n_z = 0,m3-1
            call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
        enddo
        Elem%flpml%ForcesFl3(:,:,:) = s1+Elem%flpml%ForcesFl3(:,:,:)

        s0 = Elem%fl%Acoeff(:,:,:,17) * Elem%flpml%Veloc(:,:,:,2)
        call DGEMM('N','N',m1*m2,m3,m3,1.,s0(:,:,:),m1*m2,htprimez,m3,0.,s1,m1*m2)
        Elem%flpml%ForcesFl3(:,:,:) = s1+Elem%flpml%ForcesFl3(:,:,:)


        Elem%fl%ForcesFl(:,:,:) = Elem%flpml%ForcesFl1(:,:,:) + Elem%flpml%ForcesFl2(:,:,:) + Elem%flpml%ForcesFl3(:,:,:)

        return
    end subroutine compute_InternalForces_PML_Elem_Fl
    !------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------

    ! ###########################################################
    ! ###########################################################
    subroutine Prediction_Elem_PML_Veloc_curve (Elem, bega, dt, hTprimex, Hprimey, Hprimez)

        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) :: bega, dt
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez

        integer :: m1,m2,m3, n_x,n_y,n_z
        real, dimension(0:2,0:2) :: inv_norm, tmp, res, norm
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVx_dxi,dVx_deta,dVx_dzeta, &
            dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1, 0:2, 0:2) :: VxMat,VyMat,VzMat, VxTerms,VyTerms,VzTerms


        ! On calcule les derivees de la vitesse par rapport a xi, eta et zeta

        m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz

        Elem%sl%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2) = Elem%sl%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%sl%Accel(:,:,:,:)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%sl%Forces(:,:,:,0), m1, 0., dVx_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%sl%Forces(:,:,n_z,0), m1, hprimey, m2, 0., dVx_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%sl%Forces(:,:,:,0), m1*m2, hprimez, m3, 0., dVx_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%sl%Forces(:,:,:,1), m1, 0., dVy_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%sl%Forces(:,:,n_z,1), m1, hprimey, m2, 0., dVy_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%sl%Forces(:,:,:,1), m1*m2, hprimez, m3, 0., dVy_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%sl%Forces(:,:,:,2), m1, 0., dVz_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%sl%Forces(:,:,n_z,2), m1, hprimey, m2, 0., dVz_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%sl%Forces(:,:,:,2), m1*m2, hprimez, m3, 0., dVz_dzeta, m1*m2)


        ! On cree VxMat, VyMat et VzMat:

        VxMat(:,:,:,0,0) = Elem%sl%Acoeff(:,:,:,0) * dVx_dxi + Elem%sl%Acoeff(:,:,:,1) * dVx_deta + Elem%sl%Acoeff(:,:,:,2) * dVx_dzeta
        VxMat(:,:,:,0,1) = Elem%sl%Acoeff(:,:,:,12) * dVx_dxi + Elem%sl%Acoeff(:,:,:,13) * dVx_deta + Elem%sl%Acoeff(:,:,:,14) * dVx_dzeta
        VxMat(:,:,:,0,2) = Elem%sl%Acoeff(:,:,:,15) * dVx_dxi + Elem%sl%Acoeff(:,:,:,16) * dVx_deta + Elem%sl%Acoeff(:,:,:,17) * dVx_dzeta
        VxMat(:,:,:,1,0) = Elem%sl%Acoeff(:,:,:,9) * dVx_dxi + Elem%sl%Acoeff(:,:,:,10) * dVx_deta + Elem%sl%Acoeff(:,:,:,11) * dVx_dzeta
        VxMat(:,:,:,1,1) = Elem%sl%Acoeff(:,:,:,3) * dVx_dxi + Elem%sl%Acoeff(:,:,:,4) * dVx_deta + Elem%sl%Acoeff(:,:,:,5) * dVx_dzeta
        VxMat(:,:,:,1,2) = Elem%sl%Acoeff(:,:,:,6) * dVx_dxi + Elem%sl%Acoeff(:,:,:,7) * dVx_deta + Elem%sl%Acoeff(:,:,:,8) * dVx_dzeta
        VxMat(:,:,:,2,0) = Elem%sl%Acoeff(:,:,:,18) * dVx_dxi + Elem%sl%Acoeff(:,:,:,19) * dVx_deta + Elem%sl%Acoeff(:,:,:,20) * dVx_dzeta
        VxMat(:,:,:,2,1) = Elem%sl%Acoeff(:,:,:,21) * dVx_dxi + Elem%sl%Acoeff(:,:,:,22) * dVx_deta + Elem%sl%Acoeff(:,:,:,23) * dVx_dzeta
        VxMat(:,:,:,2,2) = Elem%sl%Acoeff(:,:,:,24) * dVx_dxi + Elem%sl%Acoeff(:,:,:,25) * dVx_deta + Elem%sl%Acoeff(:,:,:,26) * dVx_dzeta

        VyMat(:,:,:,0,0) = Elem%sl%Acoeff(:,:,:,0) * dVy_dxi + Elem%sl%Acoeff(:,:,:,1) * dVy_deta + Elem%sl%Acoeff(:,:,:,2) * dVy_dzeta
        VyMat(:,:,:,0,1) = Elem%sl%Acoeff(:,:,:,12) * dVy_dxi + Elem%sl%Acoeff(:,:,:,13) * dVy_deta + Elem%sl%Acoeff(:,:,:,14) * dVy_dzeta
        VyMat(:,:,:,0,2) = Elem%sl%Acoeff(:,:,:,15) * dVy_dxi + Elem%sl%Acoeff(:,:,:,16) * dVy_deta + Elem%sl%Acoeff(:,:,:,17) * dVy_dzeta
        VyMat(:,:,:,1,0) = Elem%sl%Acoeff(:,:,:,9) * dVy_dxi + Elem%sl%Acoeff(:,:,:,10) * dVy_deta + Elem%sl%Acoeff(:,:,:,11) * dVy_dzeta
        VyMat(:,:,:,1,1) = Elem%sl%Acoeff(:,:,:,3) * dVy_dxi + Elem%sl%Acoeff(:,:,:,4) * dVy_deta + Elem%sl%Acoeff(:,:,:,5) * dVy_dzeta
        VyMat(:,:,:,1,2) = Elem%sl%Acoeff(:,:,:,6) * dVy_dxi + Elem%sl%Acoeff(:,:,:,7) * dVy_deta + Elem%sl%Acoeff(:,:,:,8) * dVy_dzeta
        VyMat(:,:,:,2,0) = Elem%sl%Acoeff(:,:,:,18) * dVy_dxi + Elem%sl%Acoeff(:,:,:,19) * dVy_deta + Elem%sl%Acoeff(:,:,:,20) * dVy_dzeta
        VyMat(:,:,:,2,1) = Elem%sl%Acoeff(:,:,:,21) * dVy_dxi + Elem%sl%Acoeff(:,:,:,22) * dVy_deta + Elem%sl%Acoeff(:,:,:,23) * dVy_dzeta
        VyMat(:,:,:,2,2) = Elem%sl%Acoeff(:,:,:,24) * dVy_dxi + Elem%sl%Acoeff(:,:,:,25) * dVy_deta + Elem%sl%Acoeff(:,:,:,26) * dVy_dzeta

        VzMat(:,:,:,0,0) = Elem%sl%Acoeff(:,:,:,0) * dVz_dxi + Elem%sl%Acoeff(:,:,:,1) * dVz_deta + Elem%sl%Acoeff(:,:,:,2) * dVz_dzeta
        VzMat(:,:,:,0,1) = Elem%sl%Acoeff(:,:,:,12) * dVz_dxi + Elem%sl%Acoeff(:,:,:,13) * dVz_deta + Elem%sl%Acoeff(:,:,:,14) * dVz_dzeta
        VzMat(:,:,:,0,2) = Elem%sl%Acoeff(:,:,:,15) * dVz_dxi + Elem%sl%Acoeff(:,:,:,16) * dVz_deta + Elem%sl%Acoeff(:,:,:,17) * dVz_dzeta
        VzMat(:,:,:,1,0) = Elem%sl%Acoeff(:,:,:,9) * dVz_dxi + Elem%sl%Acoeff(:,:,:,10) * dVz_deta + Elem%sl%Acoeff(:,:,:,11) * dVz_dzeta
        VzMat(:,:,:,1,1) = Elem%sl%Acoeff(:,:,:,3) * dVz_dxi + Elem%sl%Acoeff(:,:,:,4) * dVz_deta + Elem%sl%Acoeff(:,:,:,5) * dVz_dzeta
        VzMat(:,:,:,1,2) = Elem%sl%Acoeff(:,:,:,6) * dVz_dxi + Elem%sl%Acoeff(:,:,:,7) * dVz_deta + Elem%sl%Acoeff(:,:,:,8) * dVz_dzeta
        VzMat(:,:,:,2,0) = Elem%sl%Acoeff(:,:,:,18) * dVz_dxi + Elem%sl%Acoeff(:,:,:,19) * dVz_deta + Elem%sl%Acoeff(:,:,:,20) * dVz_dzeta
        VzMat(:,:,:,2,1) = Elem%sl%Acoeff(:,:,:,21) * dVz_dxi + Elem%sl%Acoeff(:,:,:,22) * dVz_deta + Elem%sl%Acoeff(:,:,:,23) * dVz_dzeta
        VzMat(:,:,:,2,2) = Elem%sl%Acoeff(:,:,:,24) * dVz_dxi + Elem%sl%Acoeff(:,:,:,25) * dVz_deta + Elem%sl%Acoeff(:,:,:,26) * dVz_dzeta


        ! On cree VxTerms, VyTerms et VzTerms:

        inv_norm = transpose(Elem%slpml%Inv_Normales)
        do n_x = 0,m1-1
            do n_y = 0,m2-1
                do n_z = 0,m3-1
                    tmp(0:2,0:2) = VxMat(n_x,n_y,n_z,0:2,0:2)
                    call DGEMM ('N', 'N', 3, 3, 3, 1., tmp, 3, inv_norm, 3, 0., res, 3)
                    VxTerms(n_x,n_y,n_z,0:2,0:2) = res(0:2,0:2)

                    tmp(0:2,0:2) = VyMat(n_x,n_y,n_z,0:2,0:2)
                    call DGEMM ('N', 'N', 3, 3, 3, 1., tmp, 3, inv_norm, 3, 0., res, 3)
                    VyTerms(n_x,n_y,n_z,0:2,0:2) = res(0:2,0:2)

                    tmp(0:2,0:2) = VzMat(n_x,n_y,n_z,0:2,0:2)
                    call DGEMM ('N', 'N', 3, 3, 3, 1., tmp, 3, inv_norm, 3, 0., res, 3)
                    VzTerms(n_x,n_y,n_z,0:2,0:2) = res(0:2,0:2)
                enddo
            enddo
        enddo


        ! On calcule Diagonal_Stress:

        norm = Elem%slpml%Normales

        Elem%slpml%Diagonal_Stress1(:,:,:,0) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1(:,:,:,0) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,0,0) + norm(1,0)*VyTerms(:,:,:,1,0) + norm(2,0)*VzTerms(:,:,:,1,0))
        Elem%slpml%Diagonal_Stress1(:,:,:,1) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1(:,:,:,1) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,1,0) + norm(1,0)*VyTerms(:,:,:,0,0) + norm(2,0)*VzTerms(:,:,:,1,0))
        Elem%slpml%Diagonal_Stress1(:,:,:,2) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1(:,:,:,2) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,1,0) + norm(1,0)*VyTerms(:,:,:,1,0) + norm(2,0)*VzTerms(:,:,:,0,0))

        Elem%slpml%Diagonal_Stress2(:,:,:,0) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2(:,:,:,0) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,0,1) + norm(1,1)*VyTerms(:,:,:,1,1) + norm(2,1)*VzTerms(:,:,:,1,1))
        Elem%slpml%Diagonal_Stress2(:,:,:,1) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2(:,:,:,1) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,1,1) + norm(1,1)*VyTerms(:,:,:,0,1) + norm(2,1)*VzTerms(:,:,:,1,1))
        Elem%slpml%Diagonal_Stress2(:,:,:,2) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2(:,:,:,2) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,1,1) + norm(1,1)*VyTerms(:,:,:,1,1) + norm(2,1)*VzTerms(:,:,:,0,1))

        Elem%slpml%Diagonal_Stress3(:,:,:,0) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3(:,:,:,0) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,0,2) + norm(1,2)*VyTerms(:,:,:,1,2) + norm(2,2)*VzTerms(:,:,:,1,2))
        Elem%slpml%Diagonal_Stress3(:,:,:,1) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3(:,:,:,1) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,1,2) + norm(1,2)*VyTerms(:,:,:,0,2) + norm(2,2)*VzTerms(:,:,:,1,2))
        Elem%slpml%Diagonal_Stress3(:,:,:,2) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3(:,:,:,2) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,1,2) + norm(1,2)*VyTerms(:,:,:,1,2) + norm(2,2)*VzTerms(:,:,:,0,2))

        Elem%slpml%Diagonal_Stress = Elem%slpml%Diagonal_Stress1 + Elem%slpml%Diagonal_Stress2 + Elem%slpml%Diagonal_Stress3


        ! On calcule Residual_Stress:

        Elem%slpml%residual_Stress1(:,:,:,0) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%residual_Stress1(:,:,:,0) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VyTerms(:,:,:,2,0) + norm(1,0)*VxTerms(:,:,:,2,0))
        Elem%slpml%residual_Stress1(:,:,:,1) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%residual_Stress1(:,:,:,1) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VzTerms(:,:,:,2,0) + norm(2,0)*VxTerms(:,:,:,2,0))
        Elem%slpml%residual_Stress1(:,:,:,2) = Elem%slpml%dumpSx(:,:,:,0) * Elem%slpml%residual_Stress1(:,:,:,2) + &
            Elem%slpml%dumpSx(:,:,:,1) * Dt * (norm(1,0)*VzTerms(:,:,:,2,0) + norm(2,0)*VyTerms(:,:,:,2,0))

        Elem%slpml%residual_Stress2(:,:,:,0) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%residual_Stress2(:,:,:,0) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VyTerms(:,:,:,2,1) + norm(1,1)*VxTerms(:,:,:,2,1))
        Elem%slpml%residual_Stress2(:,:,:,1) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%residual_Stress2(:,:,:,1) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VzTerms(:,:,:,2,1) + norm(2,1)*VxTerms(:,:,:,2,1))
        Elem%slpml%residual_Stress2(:,:,:,2) = Elem%slpml%dumpSy(:,:,:,0) * Elem%slpml%residual_Stress2(:,:,:,2) + &
            Elem%slpml%dumpSy(:,:,:,1) * Dt * (norm(1,1)*VzTerms(:,:,:,2,1) + norm(2,1)*VyTerms(:,:,:,2,1))

        Elem%slpml%residual_Stress3(:,:,:,0) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%residual_Stress3(:,:,:,0) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VyTerms(:,:,:,2,2) + norm(1,2)*VxTerms(:,:,:,2,2))
        Elem%slpml%residual_Stress3(:,:,:,1) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%residual_Stress3(:,:,:,1) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VzTerms(:,:,:,2,2) + norm(2,2)*VxTerms(:,:,:,2,2))
        Elem%slpml%residual_Stress3(:,:,:,2) = Elem%slpml%dumpSz(:,:,:,0) * Elem%slpml%residual_Stress3(:,:,:,2) + &
            Elem%slpml%dumpSz(:,:,:,1) * Dt * (norm(1,2)*VzTerms(:,:,:,2,2) + norm(2,2)*VyTerms(:,:,:,2,2))

        Elem%slpml%residual_Stress = Elem%slpml%residual_Stress1 + Elem%slpml%residual_Stress2 + Elem%slpml%residual_Stress3


        return
    end subroutine Prediction_Elem_PML_Veloc_curve


    subroutine init_element(el)
        type(element), intent(inout) :: el

        el%mat_index=-1
        el%ngllx=0
        el%nglly=0
        el%ngllz=0
        el%PML = .false.
        el%solid = .true.

    end subroutine init_element


!    subroutine compute_local_coordinates(el, coord, mat, xs,ys,zs, xi, eta, zeta)
!        use ssubdomains
!        use shape_geom_3d
!        type(element), intent(in) :: el
!        type(Subdomain), intent(in) :: mat
!        real, intent(in) :: xs,ys,zs
!        real, dimension(0:,0:), intent(in) :: coord
!        real, intent(out) :: xi, eta, zeta
!        real :: x,y,z
!        integer :: i,j,k, ipt
!        integer :: ic,jc,kc
!        real :: dmin, d, derr
!        real, dimension(0:2,0:2) :: LocInvGrad, Jac
!        real, dimension(0:7) :: xco, yco, zco
!        real, dimension(0:2) :: vx,dx
!
!        !
!        ! On cherche le pt de gauss le plus proche
!        !
!        dmin = 1e30
!        do i = 0,mat%ngllx-1
!            do j = 0,mat%nglly-1
!                do k = 0,mat%ngllz-1
!                    ipt = el%Iglobnum(i,j,k)
!                    d = sqrt ((coord(0,ipt)-xs)**2 + (coord(1,ipt)-ys)**2 + (coord(2,ipt)-zs)**2)
!                    if (d <= dmin) then
!                        dmin = d
!                        ic = i
!                        jc = j
!                        kc = k
!                    endif
!                enddo
!            enddo
!        enddo
!        do i=0,7
!            xco(i) = coord(0, el%Control_Nodes(i))
!            yco(i) = coord(1, el%Control_Nodes(i))
!            zco(i) = coord(2, el%Control_Nodes(i))
!        end do
!        do
!            x = f_p(xco, xi, eta, zeta)
!            y = f_p(yco, xi, eta, zeta)
!            z = f_p(zco, xi, eta, zeta)
!
!            LocInvGrad(0,0) = der_dx_dxi(xco,eta,zeta)
!            LocInvGrad(1,0) = der_dx_deta(xco,xi,zeta)
!            LocInvGrad(2,0) = der_dx_dzeta(xco,xi,eta)
!            LocInvGrad(0,1) = der_dy_dxi(yco,eta,zeta)
!            LocInvGrad(1,1) = der_dy_deta(yco,xi,zeta)
!            LocInvGrad(2,1) = der_dy_dzeta(yco,xi,eta)
!            LocInvGrad(0,2) = der_dz_dxi(zco,eta,zeta)
!            LocInvGrad(1,2) = der_dz_deta(zco,xi,zeta)
!            LocInvGrad(2,2) = der_dz_dzeta(zco,xi,eta)
!
!            call invert_3d(LocInvGrad,Jac)
!
!            derr = (x-xs)**2+(y-ys)**2+(z-zs)**2
!            if (derr<1e-14) exit
!        end do
!
!    end subroutine compute_local_coordinates

end module selement

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
