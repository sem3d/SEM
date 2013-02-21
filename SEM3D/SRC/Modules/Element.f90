!>
!!\file Element.f90
!!\brief contient les méthodes qui assure la gestion du type Element.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module selement

    implicit none

    type :: element_pml
       ! TODO move pml related data here
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: Residual_Stress1, Residual_Stress2, Residual_Stress3
       real, dimension(:,:,:,:), allocatable :: DumpSx,DumpSy,DumpSz
       real, dimension(:,:,:,:), allocatable :: Forces1,Forces2,Forces3
       real, dimension(:,:,:,:), allocatable :: Veloc1,Veloc2,Veloc3
       ! FPML
       real, dimension(:,:,:), allocatable :: Isx, Isy, Isz, Ivx, Ivy, Ivz
       real, dimension(:,:,:,:), allocatable :: Iveloc1, Iveloc2, Iveloc3
       real, dimension(:,:,:,:), allocatable :: I_Diagonal_Stress1, I_Diagonal_Stress2, I_Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: I_Residual_Stress1, I_Residual_Stress2
       real, dimension(:,:,:,:), allocatable :: DumpVx,DumpVy,DumpVz, DumpMass
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress, Residual_Stress
       real, dimension(:,:), allocatable :: Normales, Inv_Normales
       real, dimension(:,:,:), allocatable :: ForcesFl1,ForcesFl2,ForcesFl3,VelPhi1,VelPhi2,VelPhi3
    end type element_pml

    type :: element

       integer :: mat_index, ngllx, nglly, ngllz
       integer, dimension (:), allocatable :: Control_nodes

       integer, dimension (0:5) :: Near_Faces, Orient_Faces
       integer, dimension (0:11) :: Near_Edges, Orient_Edges

       integer, dimension (0:7) :: Near_Vertices
       integer, dimension (:,:,:), allocatable :: Iglobnum,Num
       real, dimension (:,:,:), allocatable :: Jacob, Density, Lambda, Mu, MassMat
       real, dimension (:,:,:), allocatable :: Kappa,Q, Qs, Qp, onemSbeta, onemPbeta, &
           epsilonvol_, &
           epsilondev_xx_,epsilondev_yy_,epsilondev_xy_,epsilondev_xz_,epsilondev_yz_
       real, dimension (:), allocatable :: wgtx, wgty, wgtz
       real, dimension(:,:,:,:), allocatable :: ACoeff, Forces,Veloc,Displ,Accel,V0
       real, dimension(:,:,:,:), allocatable :: Cij, &
           factor_common_3, alphaval_3,betaval_3,gammaval_3, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_, &
           factor_common_P, alphaval_P,betaval_P,gammaval_P, R_vol_
       real, dimension(:,:,:,:,:), allocatable :: InvGrad
       ! fluid part
       real, dimension(:,:,:), allocatable:: Phi,VelPhi0,VelPhi,AccelPhi
       ! PML allocation
       logical :: PML, FPML
       type(element_pml), allocatable :: spml
       ! fluid part
       real, dimension(:,:,:), allocatable:: ForcesFl

       ! solid-fluid
       logical  :: solid
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

        Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%Displ
        Elem%V0 = Elem%Veloc

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
        Elem%VelPhi0(:,:,:) = Elem%VelPhi(:,:,:)
        Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) = Elem%Phi(:,:,:)
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
            Elem%Forces(1:ngllx-2,1:nglly-2, 1:ngllz-2,i) = Elem%MassMat(:,:,:)*    &
                Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
        enddo
        Elem%Veloc(:,:,:,:) = Elem%v0(:,:,:,:) + dt * Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,:)
        Elem%Accel  =  (Elem%Veloc-Elem%V0)/dt
        Elem%Displ = Elem%Displ + dt * Elem%Veloc

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

        Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) = Elem%MassMat(:,:,:)*    &
            Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%VelPhi(:,:,:) = Elem%VelPhi0(:,:,:)+ dt * Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%AccelPhi  = (Elem%VelPhi-Elem%VelPhi0)/dt
        Elem%Phi = Elem%Phi + dt*Elem%VelPhi

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
        call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%Forces(:,:,:,0),m1,0.,dUx_dxi,m1)
        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,Elem%Forces(:,:,n_z,0),m1,hprimey,m2,0.,dUx_deta(:,:,n_z),m1)
        enddo
        call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%Forces(:,:,:,0),m1*m2,hprimez,m3,0.,dUx_dzeta,m1*m2)
        ! dUy_(dxi,deta,dzeta)
        call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%Forces(:,:,:,1),m1,0.,dUy_dxi,m1)
        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,Elem%Forces(:,:,n_z,1),m1,hprimey,m2,0.,dUy_deta(:,:,n_z),m1)
        enddo
        call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%Forces(:,:,:,1),m1*m2,hprimez,m3,0.,dUy_dzeta,m1*m2)
        ! dUz_(dxi,deta,dzeta)
        call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%Forces(:,:,:,2),m1,0.,dUz_dxi,m1)
        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,Elem%Forces(:,:,n_z,2),m1,hprimey,m2,0.,dUz_deta(:,:,n_z),m1)
        enddo
        call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%Forces(:,:,:,2),m1*m2,hprimez,m3,0.,dUz_dzeta,m1*m2)


        !- Internal forces
        t1 = Elem%Acoeff(:,:,:,0)*dUx_dxi + Elem%Acoeff(:,:,:,1)*dUx_deta + Elem%Acoeff(:,:,:,2)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,3)*dUy_dxi + Elem%Acoeff(:,:,:,4)*dUy_deta + Elem%Acoeff(:,:,:,5)*dUy_dzeta +  &
            Elem%Acoeff(:,:,:,6)*dUz_dxi + Elem%Acoeff(:,:,:,7)*dUz_deta + Elem%Acoeff(:,:,:,8)*dUz_dzeta

        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uxloc,m1)

        t1 = Elem%Acoeff(:,:,:,1)*dUx_dxi + Elem%Acoeff(:,:,:,9)*dUx_deta + Elem%Acoeff(:,:,:,10)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,11)*dUy_dxi + Elem%Acoeff(:,:,:,12)*dUy_deta + Elem%Acoeff(:,:,:,13)*dUy_dzeta + &
            Elem%Acoeff(:,:,:,14)*dUz_dxi + Elem%Acoeff(:,:,:,15)*dUz_deta + Elem%Acoeff(:,:,:,16)*dUz_dzeta

        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Uxloc = s0 + Uxloc

        t1 = Elem%Acoeff(:,:,:,2)*dUx_dxi + Elem%Acoeff(:,:,:,10)*dUx_deta + Elem%Acoeff(:,:,:,17)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,18)*dUy_dxi + Elem%Acoeff(:,:,:,19)*dUy_deta + Elem%Acoeff(:,:,:,20)*dUy_dzeta +&
            Elem%Acoeff(:,:,:,21)*dUz_dxi + Elem%ACoeff(:,:,:,22)*dUz_deta + Elem%Acoeff(:,:,:,23)*dUz_dzeta

        call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Uxloc = s0 + Uxloc

        t1 = Elem%Acoeff(:,:,:,3)*dUx_dxi + Elem%Acoeff(:,:,:,11)*dUx_deta + Elem%Acoeff(:,:,:,18)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,24)*dUy_dxi + Elem%Acoeff(:,:,:,25)*dUy_deta + Elem%Acoeff(:,:,:,26)*dUy_dzeta + &
            Elem%Acoeff(:,:,:,27)*dUz_dxi + Elem%Acoeff(:,:,:,28)*dUz_deta + Elem%Acoeff(:,:,:,29)*dUz_dzeta

        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uyloc,m1)

        t1 = Elem%Acoeff(:,:,:,4)*dUx_dxi + Elem%Acoeff(:,:,:,12)*dUx_deta + Elem%Acoeff(:,:,:,19)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,25)*dUy_dxi + Elem%Acoeff(:,:,:,30)*dUy_deta + Elem%Acoeff(:,:,:,31)*dUy_dzeta + &
            Elem%Acoeff(:,:,:,32)*dUz_dxi + Elem%Acoeff(:,:,:,33)*dUz_deta + Elem%Acoeff(:,:,:,34)*dUz_dzeta

        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Uyloc = s0 + Uyloc

        t1 = Elem%Acoeff(:,:,:,5)*dUx_dxi + Elem%Acoeff(:,:,:,13)*dUx_deta + Elem%Acoeff(:,:,:,20)* dUx_dzeta +&
            Elem%Acoeff(:,:,:,26)*dUy_dxi + Elem%Acoeff(:,:,:,31)*dUy_deta + Elem%Acoeff(:,:,:,35)*dUy_dzeta + &
            Elem%Acoeff(:,:,:,36)*dUz_dxi + Elem%Acoeff(:,:,:,37)*dUz_deta + Elem%Acoeff(:,:,:,38)*dUz_dzeta

        call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Uyloc = s0 + Uyloc

        t1 = Elem%Acoeff(:,:,:,6)*dUx_dxi + Elem%Acoeff(:,:,:,14)*dUx_deta + Elem%Acoeff(:,:,:,21)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,27)*dUy_dxi + Elem%Acoeff(:,:,:,32)*dUy_deta + Elem%Acoeff(:,:,:,36)*dUy_dzeta + &
            Elem%Acoeff(:,:,:,39)*dUz_dxi + Elem%Acoeff(:,:,:,40)*dUz_deta + Elem%Acoeff(:,:,:,41)*dUz_dzeta

        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uzloc,m1)

        t1 = Elem%Acoeff(:,:,:,7)*dUx_dxi + Elem%Acoeff(:,:,:,15)*dUx_deta + Elem%Acoeff(:,:,:,22)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,28)*dUy_dxi + Elem%Acoeff(:,:,:,33)*dUy_deta + Elem%Acoeff(:,:,:,37)*dUy_dzeta + &
            Elem%Acoeff(:,:,:,40)*dUz_dxi + Elem%Acoeff(:,:,:,42)*dUz_deta + Elem%Acoeff(:,:,:,43)*dUz_dzeta

        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Uzloc = s0 + Uzloc

        t1 = Elem%Acoeff(:,:,:,8)*dUx_dxi + Elem%Acoeff(:,:,:,16)*dUx_deta + Elem%Acoeff(:,:,:,23)*dUx_dzeta + &
            Elem%Acoeff(:,:,:,29)*dUy_dxi + Elem%Acoeff(:,:,:,34)*dUy_deta + Elem%Acoeff(:,:,:,38)*dUy_dzeta + &
            Elem%Acoeff(:,:,:,41)*dUz_dxi + Elem%Acoeff(:,:,:,43)*dUz_deta + Elem%Acoeff(:,:,:,44)*dUz_dzeta

        call DGEMM ('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Uzloc = Uzloc + s0

        Elem%Forces(:,:,:,0) = Uxloc
        Elem%Forces(:,:,:,1) = Uyloc
        Elem%Forces(:,:,:,2) = Uzloc

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

        !- Modification: potential -> density*potentiel
        Elem%ForcesFl(:,:,:) = Elem%density(:,:,:)*Elem%ForcesFl(:,:,:)

        !- gradients at GLLs points
        ! d(rho*Phi)_dxi
        call DGEMM('N','N',m1,m2*m3,m1,1d0,htprimex,m1,Elem%ForcesFl(:,:,:),m1,0d0,dPhi_dxi,m1)
        ! d(rho*Phi)_deta
        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1d0,Elem%ForcesFl(:,:,n_z),m1,hprimey,m2,0d0,dPhi_deta(:,:,n_z),m1)
        enddo
        ! d(rho*Phi)_dzeta
        call DGEMM('N','N',m1*m2,m3,m3,1d0,Elem%ForcesFl(:,:,:),m1*m2,hprimez,m3,0d0,dPhi_dzeta,m1*m2)


        !- Internal Forces
        t1 = Elem%Acoeff(:,:,:,0)*dPhi_dxi + Elem%Acoeff(:,:,:,1)*dPhi_deta + Elem%Acoeff(:,:,:,2)*dPhi_dzeta
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Floc,m1)

        t1 = Elem%Acoeff(:,:,:,1)*dPhi_dxi + Elem%Acoeff(:,:,:,3)*dPhi_deta + Elem%Acoeff(:,:,:,4)*dPhi_dzeta
        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
        enddo
        Floc = s0 + Floc

        t1 = Elem%Acoeff(:,:,:,2)*dPhi_dxi + Elem%Acoeff(:,:,:,4)*dPhi_deta + Elem%Acoeff(:,:,:,5)*dPhi_dzeta
        call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
        Floc = s0 + Floc

        !-
        Elem%ForcesFl(:,:,:) = Floc

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

        integer :: m1, m2,m3 ,n_z


        m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz
  ! useless as bega = 1/2 in general
        Elem%Forces(1:m1-2,1:m2-2, 1:m3-2,0:2)  = Elem%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%Accel(:,:,:,:)

  ! partial of velocity components with respect to xi,eta,zeta
        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,0) ,m1, 0., dVx_dxi, m1 )
        do n_z = 0,Elem%ngllz-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,0), m1, hprimey ,m2, 0., dVx_deta(:,:,n_z),m1 )
        enddo
        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,0), m1*m2, hprimez ,m3, 0., dVx_dzeta, m1*m2 )

        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,1) ,m1, 0., dVy_dxi, m1 )
        do n_z = 0,Elem%ngllz-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,1), m1, hprimey ,m2, 0., dVy_deta(:,:,n_z),m1 )
        enddo
        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,1), m1*m2, hprimez ,m3, 0., dVy_dzeta, m1*m2 )

        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,2) ,m1, 0., dVz_dxi, m1 )
        do n_z = 0,Elem%ngllz-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,2), m1, hprimey ,m2, 0., dVz_deta(:,:,n_z),m1)
        enddo
        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,2), m1*m2, hprimez ,m3, 0., dVz_dzeta, m1*m2 )


  ! Stress_xx
   ! (Stress_xx)^x
        Elem%spml%Diagonal_Stress1 (:,:,:,0) = Elem%spml%DumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1 (:,:,:,0) + &
            Elem%spml%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + Elem%Acoeff(:,:,:,2) * dVx_dzeta)
   ! (Stress_xx)^y
        Elem%spml%Diagonal_Stress2 (:,:,:,0) = Elem%spml%DumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2 (:,:,:,0) + &
            Elem%spml%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
   ! (Stress_xx)^z
        Elem%spml%Diagonal_Stress3 (:,:,:,0) = Elem%spml%DumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3 (:,:,:,0) + &
            Elem%spml%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

  ! Stress_yy
        Elem%spml%Diagonal_Stress1 (:,:,:,1) = Elem%spml%DumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1 (:,:,:,1) + &
            Elem%spml%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
        Elem%spml%Diagonal_Stress2 (:,:,:,1) = Elem%spml%DumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2 (:,:,:,1) + &
            Elem%spml%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta + Elem%Acoeff(:,:,:,14) * dVy_dzeta)
        Elem%spml%Diagonal_Stress3 (:,:,:,1) = Elem%spml%DumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3 (:,:,:,1) + &
            Elem%spml%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

  ! Stress_zz
        Elem%spml%Diagonal_Stress1 (:,:,:,2) = Elem%spml%DumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1 (:,:,:,2) + &
            Elem%spml%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
        Elem%spml%Diagonal_Stress2 (:,:,:,2) = Elem%spml%DumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2 (:,:,:,2) + &
            Elem%spml%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
        Elem%spml%Diagonal_Stress3 (:,:,:,2) = Elem%spml%DumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3 (:,:,:,2) + &
            Elem%spml%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta + Elem%Acoeff(:,:,:,17) * dVz_dzeta)

        Elem%spml%Diagonal_Stress = Elem%spml%Diagonal_Stress1 + Elem%spml%Diagonal_Stress2 + Elem%spml%Diagonal_Stress3

  ! Stress_xy
        Elem%spml%Residual_Stress1 (:,:,:,0) = Elem%spml%DumpSx(:,:,:,0) * Elem%spml%Residual_Stress1 (:,:,:,0) + &
            Elem%spml%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta + Elem%Acoeff(:,:,:,20) * dVy_dzeta)
        Elem%spml%Residual_Stress2 (:,:,:,0) = Elem%spml%DumpSy(:,:,:,0) * Elem%spml%Residual_Stress2 (:,:,:,0) + &
            Elem%spml%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta + Elem%Acoeff(:,:,:,23) * dVx_dzeta)
        Elem%spml%Residual_Stress3 (:,:,:,0) = Elem%spml%DumpSz(:,:,:,0) * Elem%spml%Residual_Stress3 (:,:,:,0)

  ! Stress_xz
        Elem%spml%Residual_Stress1 (:,:,:,1) = Elem%spml%DumpSx(:,:,:,0) * Elem%spml%Residual_Stress1 (:,:,:,1) + &
            Elem%spml%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta + Elem%Acoeff(:,:,:,20) * dVz_dzeta)
        Elem%spml%Residual_Stress2 (:,:,:,1) = Elem%spml%DumpSy(:,:,:,0) * Elem%spml%Residual_Stress2 (:,:,:,1)
        Elem%spml%Residual_Stress3 (:,:,:,1) = Elem%spml%DumpSz(:,:,:,0) * Elem%spml%Residual_Stress3 (:,:,:,1) + &
            Elem%spml%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta + Elem%Acoeff(:,:,:,26) * dVx_dzeta)

  ! Stress_yz
        Elem%spml%Residual_Stress1 (:,:,:,2) = Elem%spml%DumpSx(:,:,:,0) * Elem%spml%Residual_Stress1 (:,:,:,2)
        Elem%spml%Residual_Stress2 (:,:,:,2) = Elem%spml%DumpSy(:,:,:,0) * Elem%spml%Residual_Stress2 (:,:,:,2) + &
            Elem%spml%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta + Elem%Acoeff(:,:,:,23) * dVz_dzeta)
        Elem%spml%Residual_Stress3 (:,:,:,2) = Elem%spml%DumpSz(:,:,:,0) * Elem%spml%Residual_Stress3 (:,:,:,2) + &
            Elem%spml%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta + Elem%Acoeff(:,:,:,26) * dVy_dzeta)

        Elem%spml%Residual_Stress = Elem%spml%Residual_Stress1 + Elem%spml%Residual_Stress2 + Elem%spml%Residual_Stress3


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
        integer :: m1,m2,m3,n_z


        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        ! prediction in the element
        Elem%ForcesFl(1:m1-2,1:m2-2,1:m3-2) = Elem%VelPhi(:,:,:) +dt*(0.5-bega)*Elem%AccelPhi(:,:,:)
        ! potential -> -pressure
        Elem%ForcesFl(:,:,:) = Elem%Density(:,:,:)*Elem%VelPhi(:,:,:)

        ! d(rho*Phi)_d(xi,eta,zeta)
        call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%ForcesFl(:,:,:),m1,0.,dVelPhi_dxi,m1)
        do n_z = 0,Elem%ngllz-1
            call DGEMM('N','N',m1,m2,m2,1.,Elem%ForcesFl(:,:,n_z),m1,hprimey,m2,0.,dVelPhi_deta(:,:,n_z),m1)
        enddo
        call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%ForcesFl(:,:,:),m1*m2,hprimez,m3,0.,dVelPhi_dzeta,m1*m2)

        ! prediction for (physical) velocity (which is the equivalent of a stress, here)
        ! V_x^x
        Elem%spml%Veloc1(:,:,:,0) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Veloc1(:,:,:,0) + Elem%spml%dumpSx(:,:,:,1) * Dt *  &
            (Elem%Acoeff(:,:,:,0) * dVelPhi_dxi + Elem%Acoeff(:,:,:,1) * dVelPhi_deta + Elem%Acoeff(:,:,:,2) * dVelPhi_dzeta)
        ! V_x^y
        Elem%spml%Veloc2(:,:,:,0) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Veloc2(:,:,:,0)
        ! V_x^z
        Elem%spml%Veloc3(:,:,:,0) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Veloc3(:,:,:,0)
        ! V_y^x
        Elem%spml%Veloc1(:,:,:,1) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Veloc1(:,:,:,1)
        ! V_y^y
        Elem%spml%Veloc2(:,:,:,1) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Veloc2(:,:,:,1) + Elem%spml%dumpSy(:,:,:,1) * Dt *  &
            (Elem%Acoeff(:,:,:,3) * dVelPhi_dxi + Elem%Acoeff(:,:,:,4) * dVelPhi_deta + Elem%Acoeff(:,:,:,5) * dVelPhi_dzeta)
        ! V_y^z
        Elem%spml%Veloc3(:,:,:,1) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Veloc3(:,:,:,1)
        ! V_z^x
        Elem%spml%Veloc1(:,:,:,2) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Veloc1(:,:,:,2)
        ! V_z^y
        Elem%spml%Veloc2(:,:,:,2) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Veloc2(:,:,:,2)
        ! V_z^z
        Elem%spml%Veloc3(:,:,:,2) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Veloc3(:,:,:,2) + Elem%spml%dumpSz(:,:,:,1) * Dt *  &
            (Elem%Acoeff(:,:,:,6) * dVelPhi_dxi + Elem%Acoeff(:,:,:,7) * dVelPhi_deta + Elem%Acoeff(:,:,:,8) * dVelPhi_dzeta)

        ! total velocity vector after dumping = sum of splitted parts
        Elem%Veloc(:,:,:,:) = Elem%spml%Veloc1(:,:,:,:) + Elem%spml%Veloc2(:,:,:,:) + Elem%spml%Veloc3(:,:,:,:)

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

        Elem%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2) = Elem%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%Accel(:,:,:,:)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,0), m1, 0., dVx_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,0), m1, hprimey, m2, 0., dVx_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,0), m1*m2, hprimez, m3, 0., dVx_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,1), m1, 0., dVy_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,1), m1, hprimey, m2, 0., dVy_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,1), m1*m2, hprimez, m3, 0., dVy_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,2), m1, 0., dVz_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,2), m1, hprimey, m2, 0., dVz_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,2), m1*m2, hprimez, m3, 0., dVz_dzeta, m1*m2)


        Stress_Ausiliar = Elem%spml%Diagonal_Stress1 (:,:,:,0)
        Elem%spml%Diagonal_Stress1 (:,:,:,0) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1 (:,:,:,0) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + &
            Elem%Acoeff(:,:,:,2) * dVx_dzeta) + Elem%spml%Isx * Elem%spml%I_Diagonal_stress1 (:,:,:,0)
        Elem%spml%I_Diagonal_Stress1  (:,:,:,0)= Fil2* Elem%spml%I_Diagonal_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress1 (:,:,:,0) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress2 (:,:,:,0)
        Elem%spml%Diagonal_Stress2 (:,:,:,0) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2 (:,:,:,0) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + &
            Elem%Acoeff(:,:,:,5) * dVy_dzeta) + Elem%spml%Isy * Elem%spml%I_Diagonal_stress2 (:,:,:,0)
        Elem%spml%I_Diagonal_Stress2  (:,:,:,0)= Fil2* Elem%spml%I_Diagonal_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress2 (:,:,:,0) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress3 (:,:,:,0)
        Elem%spml%Diagonal_Stress3 (:,:,:,0) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3 (:,:,:,0) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + &
            Elem%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%spml%Isz * Elem%spml%I_Diagonal_stress3 (:,:,:,0)
        Elem%spml%I_Diagonal_Stress3  (:,:,:,0)= Fil2* Elem%spml%I_Diagonal_Stress3(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress3 (:,:,:,0) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress1 (:,:,:,1)
        Elem%spml%Diagonal_Stress1 (:,:,:,1) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1 (:,:,:,1) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%spml%Isx * Elem%spml%I_Diagonal_stress1 (:,:,:,1)
        Elem%spml%I_Diagonal_Stress1  (:,:,:,1)= Fil2* Elem%spml%I_Diagonal_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress1 (:,:,:,1) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress2 (:,:,:,1)
        Elem%spml%Diagonal_Stress2 (:,:,:,1) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2 (:,:,:,1) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta + Elem%Acoeff(:,:,:,14) * dVy_dzeta) + Elem%spml%Isy * Elem%spml%I_Diagonal_stress2 (:,:,:,1)
        Elem%spml%I_Diagonal_Stress2  (:,:,:,1)= Fil2* Elem%spml%I_Diagonal_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress2 (:,:,:,1) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress3 (:,:,:,1)
        Elem%spml%Diagonal_Stress3 (:,:,:,1) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3 (:,:,:,1) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%spml%Isz * Elem%spml%I_Diagonal_stress3 (:,:,:,1)
        Elem%spml%I_Diagonal_Stress3  (:,:,:,1)= Fil2* Elem%spml%I_Diagonal_Stress3(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress3 (:,:,:,1) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress1 (:,:,:,2)
        Elem%spml%Diagonal_Stress1 (:,:,:,2) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1 (:,:,:,2) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%spml%Isx * Elem%spml%I_Diagonal_stress1 (:,:,:,2)
        Elem%spml%I_Diagonal_Stress1  (:,:,:,2)= Fil2* Elem%spml%I_Diagonal_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress1 (:,:,:,2) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress2 (:,:,:,2)
        Elem%spml%Diagonal_Stress2 (:,:,:,2) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2 (:,:,:,2) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)+ Elem%spml%Isy * Elem%spml%I_Diagonal_stress2 (:,:,:,2)
        Elem%spml%I_Diagonal_Stress2  (:,:,:,2)= Fil2* Elem%spml%I_Diagonal_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress2 (:,:,:,2) )

        Stress_Ausiliar = Elem%spml%Diagonal_Stress3 (:,:,:,2)
        Elem%spml%Diagonal_Stress3 (:,:,:,2) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3 (:,:,:,2) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta + Elem%Acoeff(:,:,:,17) * dVz_dzeta) + Elem%spml%Isz * Elem%spml%I_Diagonal_stress3 (:,:,:,2)
        Elem%spml%I_Diagonal_Stress3  (:,:,:,2)= Fil2* Elem%spml%I_Diagonal_Stress3(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%Diagonal_Stress3 (:,:,:,2) )


        Elem%spml%Diagonal_Stress = Elem%spml%Diagonal_Stress1 + Elem%spml%Diagonal_Stress2 + Elem%spml%Diagonal_Stress3


        Stress_Ausiliar = Elem%spml%residual_Stress1 (:,:,:,0)
        Elem%spml%residual_Stress1 (:,:,:,0) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%residual_Stress1 (:,:,:,0) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta + Elem%Acoeff(:,:,:,20) * dVy_dzeta) + Elem%spml%Isx * Elem%spml%I_Residual_stress1(:,:,:,0)
        Elem%spml%I_Residual_Stress1  (:,:,:,0)= Fil2* Elem%spml%I_Residual_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%residual_Stress1 (:,:,:,0) )

        Stress_Ausiliar = Elem%spml%residual_Stress2 (:,:,:,0)
        Elem%spml%residual_Stress2 (:,:,:,0) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%residual_Stress2 (:,:,:,0) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta + Elem%Acoeff(:,:,:,23) * dVx_dzeta) + Elem%spml%Isy * Elem%spml%I_Residual_stress2(:,:,:,0)
        Elem%spml%I_Residual_Stress2  (:,:,:,0) = Fil2* Elem%spml%I_Residual_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%residual_Stress2 (:,:,:,0) )

        Stress_Ausiliar = Elem%spml%residual_Stress1 (:,:,:,1)
        Elem%spml%residual_Stress1 (:,:,:,1) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%residual_Stress1 (:,:,:,1) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta + Elem%Acoeff(:,:,:,20) * dVz_dzeta) + Elem%spml%Isx * Elem%spml%I_Residual_stress1(:,:,:,1)
        Elem%spml%I_Residual_Stress1  (:,:,:,1)= Fil2* Elem%spml%I_Residual_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%residual_Stress1 (:,:,:,1) )

        Stress_Ausiliar = Elem%spml%residual_Stress2 (:,:,:,1)
        Elem%spml%residual_Stress2 (:,:,:,1) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%residual_Stress2 (:,:,:,1) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta + Elem%Acoeff(:,:,:,26) * dVx_dzeta) + Elem%spml%Isz * Elem%spml%I_Residual_stress2(:,:,:,1)
        Elem%spml%I_Residual_Stress2  (:,:,:,1) = Fil2* Elem%spml%I_Residual_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%residual_Stress2 (:,:,:,1) )

        Stress_Ausiliar = Elem%spml%residual_Stress1 (:,:,:,2)
        Elem%spml%residual_Stress1 (:,:,:,2) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%residual_Stress1 (:,:,:,2) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta + Elem%Acoeff(:,:,:,23) * dVz_dzeta) + Elem%spml%Isy * Elem%spml%I_Residual_stress1(:,:,:,2)
        Elem%spml%I_Residual_Stress1  (:,:,:,2)= Fil2* Elem%spml%I_Residual_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%residual_Stress1 (:,:,:,2) )

        Stress_Ausiliar = Elem%spml%residual_Stress2 (:,:,:,2)
        Elem%spml%residual_Stress2 (:,:,:,2) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%residual_Stress2 (:,:,:,2) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta + Elem%Acoeff(:,:,:,26) * dVy_dzeta)+ Elem%spml%Isz * Elem%spml%I_Residual_stress2(:,:,:,2)
        Elem%spml%I_Residual_Stress2  (:,:,:,2) = Fil2* Elem%spml%I_Residual_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * &
            (Stress_Ausiliar +Elem%spml%residual_Stress2 (:,:,:,2) )

        Elem%spml%residual_Stress = Elem%spml%residual_Stress1 + Elem%spml%residual_Stress2

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
            Elem%spml%Veloc1(:,:,:,i) = Elem%spml%dumpVx(:,:,:,0) * Elem%spml%Veloc1(:,:,:,i) +    &
                dt * Elem%spml%dumpVx(:,:,:,1)*Elem%spml%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
            Elem%spml%Veloc2(:,:,:,i) = Elem%spml%dumpVy(:,:,:,0) * Elem%spml%Veloc2(:,:,:,i) +    &
                dt * Elem%spml%dumpVy(:,:,:,1)*Elem%spml%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
            Elem%spml%Veloc3(:,:,:,i) = Elem%spml%dumpVz(:,:,:,0) * Elem%spml%Veloc3(:,:,:,i) +    &
                dt * Elem%spml%dumpVz(:,:,:,1)*Elem%spml%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
        enddo

        Elem%Veloc = Elem%spml%Veloc1 + Elem%spml%Veloc2 + Elem%spml%Veloc3
        ! Usefull only for traces and debug
        Elem%Displ = Elem%Displ + dt * Elem%Veloc

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

        Elem%spml%VelPhi1(:,:,:) = Elem%spml%dumpVx(:,:,:,0) * Elem%spml%VelPhi1(:,:,:) +    &
            dt * Elem%spml%dumpVx(:,:,:,1)*Elem%spml%ForcesFl1(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%spml%VelPhi2(:,:,:) = Elem%spml%dumpVy(:,:,:,0) * Elem%spml%VelPhi2(:,:,:) +    &
            dt * Elem%spml%dumpVy(:,:,:,1)*Elem%spml%ForcesFl2(1:ngllx-2,1:nglly-2,1:ngllz-2)
        Elem%spml%VelPhi3(:,:,:) = Elem%spml%dumpVz(:,:,:,0) * Elem%spml%VelPhi3(:,:,:) +    &
            dt * Elem%spml%dumpVz(:,:,:,1)*Elem%spml%ForcesFl3(1:ngllx-2,1:nglly-2,1:ngllz-2)

        Elem%VelPhi = Elem%spml%VelPhi1 + Elem%spml%VelPhi2 + Elem%spml%VelPhi3

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
            Ausiliar_Velocity =  Elem%spml%Veloc1(:,:,:,i)
            Elem%spml%Veloc1(:,:,:,i) = Elem%spml%dumpVx(:,:,:,0) * Elem%spml%Veloc1(:,:,:,i) + dt * Elem%spml%dumpVx(:,:,:,1)*Elem%spml%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%spml%Ivx * Elem%spml%Iveloc1 (:,:,:,i)
            Elem%spml%Iveloc1 (:,:,:,i)  = Fil2*Elem%spml%Iveloc1 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%spml%Veloc1(:,:,:,i))

            Ausiliar_Velocity =  Elem%spml%Veloc2(:,:,:,i)
            Elem%spml%Veloc2(:,:,:,i) = Elem%spml%dumpVy(:,:,:,0) * Elem%spml%Veloc2(:,:,:,i) + dt * Elem%spml%dumpVy(:,:,:,1)*Elem%spml%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%spml%Ivy * Elem%spml%Iveloc2 (:,:,:,i)
            Elem%spml%Iveloc2 (:,:,:,i)  = Fil2*Elem%spml%Iveloc2 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%spml%Veloc2(:,:,:,i))

            Ausiliar_Velocity =  Elem%spml%Veloc3(:,:,:,i)
            Elem%spml%Veloc3(:,:,:,i) = Elem%spml%dumpVz(:,:,:,0) * Elem%spml%Veloc3(:,:,:,i) + dt * Elem%spml%dumpVz(:,:,:,1)*Elem%spml%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%spml%Ivz * Elem%spml%Iveloc3 (:,:,:,i)
            Elem%spml%Iveloc3 (:,:,:,i)  = Fil2*Elem%spml%Iveloc3 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%spml%Veloc3(:,:,:,i))
        enddo

        Elem%Veloc = Elem%spml%Veloc1 + Elem%spml%Veloc2 + Elem%spml%Veloc3

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

        s0 = Elem%Acoeff(:,:,:,27) * Elem%spml%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%spml%residual_Stress(:,:,:,0) + &
            Elem%Acoeff(:,:,:,29) * Elem%spml%residual_Stress(:,:,:,1)
        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1, s0(:,:,:), m1, 0., s1, m1 )
        Elem%spml%Forces1(:,:,:,0) = s1

        s0 = Elem%Acoeff(:,:,:,30) * Elem%spml%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%spml%residual_Stress(:,:,:,0) + &
            Elem%Acoeff(:,:,:,32) * Elem%spml%residual_Stress(:,:,:,1)
        do n_z = 0,m3-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
        enddo
        Elem%spml%Forces2(:,:,:,0) = s1

        s0 = Elem%Acoeff(:,:,:,33) * Elem%spml%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%spml%residual_Stress(:,:,:,0) + &
            Elem%Acoeff(:,:,:,35) * Elem%spml%residual_Stress(:,:,:,1)

        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(:,:,:), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
        Elem%spml%Forces3(:,:,:,0) = s1

        s0 = Elem%Acoeff(:,:,:,27) * Elem%spml%residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%spml%Diagonal_Stress(:,:,:,1) + &
            Elem%Acoeff(:,:,:,29) * Elem%spml%residual_Stress(:,:,:,2)
        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(:,:,:) ,m1, 0., s1, m1 )
        Elem%spml%Forces1(:,:,:,1) = s1

        s0 = Elem%Acoeff(:,:,:,30) * Elem%spml%residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%spml%Diagonal_Stress(:,:,:,1) + &
            Elem%Acoeff(:,:,:,32) * Elem%spml%residual_Stress(:,:,:,2)
        do n_z = 0,m3-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
        enddo
        Elem%spml%Forces2(:,:,:,1) = s1

        s0 = Elem%Acoeff(:,:,:,33) * Elem%spml%residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%spml%Diagonal_Stress(:,:,:,1) + &
            Elem%Acoeff(:,:,:,35) * Elem%spml%residual_Stress(:,:,:,2)

        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(:,:,:), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
        Elem%spml%Forces3(:,:,:,1) = s1


        s0 = Elem%Acoeff(:,:,:,27) * Elem%spml%residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,28) * Elem%spml%residual_Stress(:,:,:,2) + &
            Elem%Acoeff(:,:,:,29) * Elem%spml%Diagonal_Stress(:,:,:,2)
        call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(:,:,:) ,m1, 0., s1, m1 )
        Elem%spml%Forces1(:,:,:,2) = s1

        s0 = Elem%Acoeff(:,:,:,30) * Elem%spml%residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,31) * Elem%spml%residual_Stress(:,:,:,2) + &
            Elem%Acoeff(:,:,:,32) * Elem%spml%Diagonal_Stress(:,:,:,2)
        do n_z = 0,m3-1
            call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
        enddo
        Elem%spml%Forces2(:,:,:,2) = s1

        s0 = Elem%Acoeff(:,:,:,33) * Elem%spml%residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,34) * Elem%spml%residual_Stress(:,:,:,2) + &
            Elem%Acoeff(:,:,:,35) * Elem%spml%Diagonal_Stress(:,:,:,2)

        call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(:,:,:), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
        Elem%spml%Forces3(:,:,:,2) = s1

        Elem%Forces = Elem%spml%Forces1 + Elem%spml%Forces2 + Elem%spml%Forces3

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

        ! forces associated to V_x
        s0 = Elem%Acoeff(:,:,:,9) * Elem%Veloc(:,:,:,0)
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(:,:,:),m1,0.,s1,m1)
        Elem%spml%ForcesFl1(:,:,:) = s1

        s0 = Elem%Acoeff(:,:,:,10) * Elem%Veloc(:,:,:,0)
        do n_z = 0,m3-1
            call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
        enddo
        Elem%spml%ForcesFl1(:,:,:) = s1+Elem%spml%ForcesFl1(:,:,:)

        s0 = Elem%Acoeff(:,:,:,11) * Elem%Veloc(:,:,:,0)
        call DGEMM('N','N',m1*m2,m3,m3,1.,s0(:,:,:),m1*m2,htprimez,m3,0.,s1,m1*m2)
        Elem%spml%ForcesFl1(:,:,:) = s1+Elem%spml%ForcesFl1(:,:,:)

        ! forces associated to V_y
        s0 = Elem%Acoeff(:,:,:,12) * Elem%Veloc(:,:,:,1)
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(:,:,:),m1,0.,s1,m1)
        Elem%spml%ForcesFl2(:,:,:) = s1

        s0 = Elem%Acoeff(:,:,:,13) * Elem%Veloc(:,:,:,1)
        do n_z = 0,m3-1
            call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
        enddo
        Elem%spml%ForcesFl2(:,:,:) = s1+Elem%spml%ForcesFl2(:,:,:)

        s0 = Elem%Acoeff(:,:,:,14) * Elem%Veloc(:,:,:,1)
        call DGEMM('N','N',m1*m2,m3,m3,1.,s0(:,:,:),m1*m2,htprimez,m3,0.,s1,m1*m2)
        Elem%spml%ForcesFl2(:,:,:) = s1+Elem%spml%ForcesFl2(:,:,:)

        ! forces associated to V_z
        s0 = Elem%Acoeff(:,:,:,15) * Elem%Veloc(:,:,:,2)
        call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(:,:,:),m1,0.,s1,m1)
        Elem%spml%ForcesFl3(:,:,:) = s1

        s0 = Elem%Acoeff(:,:,:,16) * Elem%Veloc(:,:,:,2)
        do n_z = 0,m3-1
            call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
        enddo
        Elem%spml%ForcesFl3(:,:,:) = s1+Elem%spml%ForcesFl3(:,:,:)

        s0 = Elem%Acoeff(:,:,:,17) * Elem%Veloc(:,:,:,2)
        call DGEMM('N','N',m1*m2,m3,m3,1.,s0(:,:,:),m1*m2,htprimez,m3,0.,s1,m1*m2)
        Elem%spml%ForcesFl3(:,:,:) = s1+Elem%spml%ForcesFl3(:,:,:)


        Elem%ForcesFl(:,:,:) = Elem%spml%ForcesFl1(:,:,:) + Elem%spml%ForcesFl2(:,:,:) + Elem%spml%ForcesFl3(:,:,:)

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

        Elem%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2) = Elem%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%Accel(:,:,:,:)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,0), m1, 0., dVx_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,0), m1, hprimey, m2, 0., dVx_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,0), m1*m2, hprimez, m3, 0., dVx_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,1), m1, 0., dVy_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,1), m1, hprimey, m2, 0., dVy_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,1), m1*m2, hprimez, m3, 0., dVy_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(:,:,:,2), m1, 0., dVz_dxi, m1)
        do n_z = 0,m3-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,2), m1, hprimey, m2, 0., dVz_deta(:,:,n_z),m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,2), m1*m2, hprimez, m3, 0., dVz_dzeta, m1*m2)


        ! On cree VxMat, VyMat et VzMat:

        VxMat(:,:,:,0,0) = Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + Elem%Acoeff(:,:,:,2) * dVx_dzeta
        VxMat(:,:,:,0,1) = Elem%Acoeff(:,:,:,12) * dVx_dxi + Elem%Acoeff(:,:,:,13) * dVx_deta + Elem%Acoeff(:,:,:,14) * dVx_dzeta
        VxMat(:,:,:,0,2) = Elem%Acoeff(:,:,:,15) * dVx_dxi + Elem%Acoeff(:,:,:,16) * dVx_deta + Elem%Acoeff(:,:,:,17) * dVx_dzeta
        VxMat(:,:,:,1,0) = Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta
        VxMat(:,:,:,1,1) = Elem%Acoeff(:,:,:,3) * dVx_dxi + Elem%Acoeff(:,:,:,4) * dVx_deta + Elem%Acoeff(:,:,:,5) * dVx_dzeta
        VxMat(:,:,:,1,2) = Elem%Acoeff(:,:,:,6) * dVx_dxi + Elem%Acoeff(:,:,:,7) * dVx_deta + Elem%Acoeff(:,:,:,8) * dVx_dzeta
        VxMat(:,:,:,2,0) = Elem%Acoeff(:,:,:,18) * dVx_dxi + Elem%Acoeff(:,:,:,19) * dVx_deta + Elem%Acoeff(:,:,:,20) * dVx_dzeta
        VxMat(:,:,:,2,1) = Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta + Elem%Acoeff(:,:,:,23) * dVx_dzeta
        VxMat(:,:,:,2,2) = Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta + Elem%Acoeff(:,:,:,26) * dVx_dzeta

        VyMat(:,:,:,0,0) = Elem%Acoeff(:,:,:,0) * dVy_dxi + Elem%Acoeff(:,:,:,1) * dVy_deta + Elem%Acoeff(:,:,:,2) * dVy_dzeta
        VyMat(:,:,:,0,1) = Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta + Elem%Acoeff(:,:,:,14) * dVy_dzeta
        VyMat(:,:,:,0,2) = Elem%Acoeff(:,:,:,15) * dVy_dxi + Elem%Acoeff(:,:,:,16) * dVy_deta + Elem%Acoeff(:,:,:,17) * dVy_dzeta
        VyMat(:,:,:,1,0) = Elem%Acoeff(:,:,:,9) * dVy_dxi + Elem%Acoeff(:,:,:,10) * dVy_deta + Elem%Acoeff(:,:,:,11) * dVy_dzeta
        VyMat(:,:,:,1,1) = Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta
        VyMat(:,:,:,1,2) = Elem%Acoeff(:,:,:,6) * dVy_dxi + Elem%Acoeff(:,:,:,7) * dVy_deta + Elem%Acoeff(:,:,:,8) * dVy_dzeta
        VyMat(:,:,:,2,0) = Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta + Elem%Acoeff(:,:,:,20) * dVy_dzeta
        VyMat(:,:,:,2,1) = Elem%Acoeff(:,:,:,21) * dVy_dxi + Elem%Acoeff(:,:,:,22) * dVy_deta + Elem%Acoeff(:,:,:,23) * dVy_dzeta
        VyMat(:,:,:,2,2) = Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta + Elem%Acoeff(:,:,:,26) * dVy_dzeta

        VzMat(:,:,:,0,0) = Elem%Acoeff(:,:,:,0) * dVz_dxi + Elem%Acoeff(:,:,:,1) * dVz_deta + Elem%Acoeff(:,:,:,2) * dVz_dzeta
        VzMat(:,:,:,0,1) = Elem%Acoeff(:,:,:,12) * dVz_dxi + Elem%Acoeff(:,:,:,13) * dVz_deta + Elem%Acoeff(:,:,:,14) * dVz_dzeta
        VzMat(:,:,:,0,2) = Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta + Elem%Acoeff(:,:,:,17) * dVz_dzeta
        VzMat(:,:,:,1,0) = Elem%Acoeff(:,:,:,9) * dVz_dxi + Elem%Acoeff(:,:,:,10) * dVz_deta + Elem%Acoeff(:,:,:,11) * dVz_dzeta
        VzMat(:,:,:,1,1) = Elem%Acoeff(:,:,:,3) * dVz_dxi + Elem%Acoeff(:,:,:,4) * dVz_deta + Elem%Acoeff(:,:,:,5) * dVz_dzeta
        VzMat(:,:,:,1,2) = Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta
        VzMat(:,:,:,2,0) = Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta + Elem%Acoeff(:,:,:,20) * dVz_dzeta
        VzMat(:,:,:,2,1) = Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta + Elem%Acoeff(:,:,:,23) * dVz_dzeta
        VzMat(:,:,:,2,2) = Elem%Acoeff(:,:,:,24) * dVz_dxi + Elem%Acoeff(:,:,:,25) * dVz_deta + Elem%Acoeff(:,:,:,26) * dVz_dzeta


        ! On cree VxTerms, VyTerms et VzTerms:

        inv_norm = transpose(Elem%spml%Inv_Normales)
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

        norm = Elem%spml%Normales

        Elem%spml%Diagonal_Stress1(:,:,:,0) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1(:,:,:,0) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,0,0) + norm(1,0)*VyTerms(:,:,:,1,0) + norm(2,0)*VzTerms(:,:,:,1,0))
        Elem%spml%Diagonal_Stress1(:,:,:,1) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1(:,:,:,1) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,1,0) + norm(1,0)*VyTerms(:,:,:,0,0) + norm(2,0)*VzTerms(:,:,:,1,0))
        Elem%spml%Diagonal_Stress1(:,:,:,2) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%Diagonal_Stress1(:,:,:,2) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VxTerms(:,:,:,1,0) + norm(1,0)*VyTerms(:,:,:,1,0) + norm(2,0)*VzTerms(:,:,:,0,0))

        Elem%spml%Diagonal_Stress2(:,:,:,0) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2(:,:,:,0) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,0,1) + norm(1,1)*VyTerms(:,:,:,1,1) + norm(2,1)*VzTerms(:,:,:,1,1))
        Elem%spml%Diagonal_Stress2(:,:,:,1) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2(:,:,:,1) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,1,1) + norm(1,1)*VyTerms(:,:,:,0,1) + norm(2,1)*VzTerms(:,:,:,1,1))
        Elem%spml%Diagonal_Stress2(:,:,:,2) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%Diagonal_Stress2(:,:,:,2) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VxTerms(:,:,:,1,1) + norm(1,1)*VyTerms(:,:,:,1,1) + norm(2,1)*VzTerms(:,:,:,0,1))

        Elem%spml%Diagonal_Stress3(:,:,:,0) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3(:,:,:,0) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,0,2) + norm(1,2)*VyTerms(:,:,:,1,2) + norm(2,2)*VzTerms(:,:,:,1,2))
        Elem%spml%Diagonal_Stress3(:,:,:,1) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3(:,:,:,1) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,1,2) + norm(1,2)*VyTerms(:,:,:,0,2) + norm(2,2)*VzTerms(:,:,:,1,2))
        Elem%spml%Diagonal_Stress3(:,:,:,2) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%Diagonal_Stress3(:,:,:,2) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VxTerms(:,:,:,1,2) + norm(1,2)*VyTerms(:,:,:,1,2) + norm(2,2)*VzTerms(:,:,:,0,2))

        Elem%spml%Diagonal_Stress = Elem%spml%Diagonal_Stress1 + Elem%spml%Diagonal_Stress2 + Elem%spml%Diagonal_Stress3


        ! On calcule Residual_Stress:

        Elem%spml%residual_Stress1(:,:,:,0) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%residual_Stress1(:,:,:,0) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VyTerms(:,:,:,2,0) + norm(1,0)*VxTerms(:,:,:,2,0))
        Elem%spml%residual_Stress1(:,:,:,1) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%residual_Stress1(:,:,:,1) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (norm(0,0)*VzTerms(:,:,:,2,0) + norm(2,0)*VxTerms(:,:,:,2,0))
        Elem%spml%residual_Stress1(:,:,:,2) = Elem%spml%dumpSx(:,:,:,0) * Elem%spml%residual_Stress1(:,:,:,2) + &
            Elem%spml%dumpSx(:,:,:,1) * Dt * (norm(1,0)*VzTerms(:,:,:,2,0) + norm(2,0)*VyTerms(:,:,:,2,0))

        Elem%spml%residual_Stress2(:,:,:,0) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%residual_Stress2(:,:,:,0) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VyTerms(:,:,:,2,1) + norm(1,1)*VxTerms(:,:,:,2,1))
        Elem%spml%residual_Stress2(:,:,:,1) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%residual_Stress2(:,:,:,1) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (norm(0,1)*VzTerms(:,:,:,2,1) + norm(2,1)*VxTerms(:,:,:,2,1))
        Elem%spml%residual_Stress2(:,:,:,2) = Elem%spml%dumpSy(:,:,:,0) * Elem%spml%residual_Stress2(:,:,:,2) + &
            Elem%spml%dumpSy(:,:,:,1) * Dt * (norm(1,1)*VzTerms(:,:,:,2,1) + norm(2,1)*VyTerms(:,:,:,2,1))

        Elem%spml%residual_Stress3(:,:,:,0) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%residual_Stress3(:,:,:,0) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VyTerms(:,:,:,2,2) + norm(1,2)*VxTerms(:,:,:,2,2))
        Elem%spml%residual_Stress3(:,:,:,1) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%residual_Stress3(:,:,:,1) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (norm(0,2)*VzTerms(:,:,:,2,2) + norm(2,2)*VxTerms(:,:,:,2,2))
        Elem%spml%residual_Stress3(:,:,:,2) = Elem%spml%dumpSz(:,:,:,0) * Elem%spml%residual_Stress3(:,:,:,2) + &
            Elem%spml%dumpSz(:,:,:,1) * Dt * (norm(1,2)*VzTerms(:,:,:,2,2) + norm(2,2)*VyTerms(:,:,:,2,2))

        Elem%spml%residual_Stress = Elem%spml%residual_Stress1 + Elem%spml%residual_Stress2 + Elem%spml%residual_Stress3


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
