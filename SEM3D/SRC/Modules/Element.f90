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
    use blas
    implicit none

    type :: element_solid
       real, dimension(:,:,:,:), allocatable :: Cij

       real, dimension(:,:,:,:), allocatable :: ACoeff

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
    end type element_fluid

    type :: element_pml
        real, dimension(:,:,:,:), allocatable :: DumpSx,DumpSy,DumpSz
        real, dimension(:,:,:,:), allocatable :: DumpMass
    end type element_pml

    type :: element_solid_pml
       ! TODO move pml related data here
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: Residual_Stress1, Residual_Stress2, Residual_Stress3
       real, dimension(:,:,:,:), allocatable :: Forces1,Forces2,Forces3
       ! FPML
       real, dimension(:,:,:), allocatable :: Isx, Isy, Isz
       real, dimension(:,:,:), allocatable :: Ivx, Ivy, Ivz
       real, dimension(:,:,:,:), allocatable :: Iveloc1, Iveloc2, Iveloc3
       real, dimension(:,:,:,:), allocatable :: I_Diagonal_Stress1, I_Diagonal_Stress2, I_Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: I_Residual_Stress1, I_Residual_Stress2
       real, dimension(:,:,:,:), allocatable :: DumpMass
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress, Residual_Stress
       real, dimension(:,:), allocatable :: Normales, Inv_Normales
       integer, dimension (:,:,:), allocatable :: ISolPml
    end type element_solid_pml

    type :: element_fluid_pml
        integer, dimension (:,:,:), allocatable :: IFluPml
        real, dimension(:,:,:,:), allocatable :: Veloc
        real, dimension(:,:,:,:), allocatable :: Veloc1,Veloc2,Veloc3
        real, dimension(:,:,:), allocatable :: ForcesFl1,ForcesFl2,ForcesFl3
        real, dimension(:,:,:), allocatable :: VelPhi1,VelPhi2,VelPhi3
    end type element_fluid_pml

    type :: element

       integer :: mat_index, ngllx, nglly, ngllz
       integer, dimension (:), allocatable :: Control_nodes
       integer, dimension (0:5) :: Near_Faces, Orient_Faces
       integer, dimension (0:11) :: Near_Edges, Orient_Edges
       integer, dimension (0:7) :: Near_Vertices

       integer, dimension (:,:,:), allocatable :: Iglobnum,Num
       real, dimension (:,:,:), allocatable :: Jacob
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
       type(element_pml), allocatable :: xpml
       type(element_solid_pml), allocatable :: slpml
       type(element_fluid_pml), allocatable :: flpml
       integer, dimension (:,:,:), allocatable :: ISol, IFlu
       real :: dist_max !! taille caracteristique de l'element
    end type element


!    interface
!       subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
!           CHARACTER*1        TRANSA, TRANSB
!           INTEGER            M, N, K, LDA, LDB, LDC
!           DOUBLE PRECISION   ALPHA, BETA
!           DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!       end subroutine DGEMM
!    end interface

contains

    integer function get_domain(el)
        use constants
        implicit none
        type(Element), intent(INOUT) :: el

        if (el%solid) then
            if (el%PML) then
                get_domain = DM_SOLID_PML
            else
                get_domain = DM_SOLID
            endif
        else ! Fluid
            if (el%PML) then
                get_domain = DM_FLUID_PML
            else
                get_domain = DM_FLUID
            endif
        endif
        return
    end function get_domain

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
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    subroutine Prediction_Elem_PML_Veloc(Elem,bega,dt,hprimex,Hprimey,Hprimez, ngll_pmls, Vitesses, Forces)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) :: bega, dt
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        integer, intent(in) :: ngll_pmls
        real, dimension (0:ngll_pmls-1,0:2), intent (IN) :: Vitesses, Forces

        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVx_dxi,dVx_deta,dVx_dzeta, &
            dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta
        real, dimension (:,:,:,:), allocatable :: Veloc
        integer :: m1, m2, m3, ind, i, j, k, i_dir

        m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz

        allocate(Veloc(0:m1-1,0:m2-1,0:m3-1,0:2))
        do i_dir = 0,2
            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                        ind = Elem%slpml%ISolPml(i,j,k)
                        Veloc(i,j,k,i_dir) = Vitesses(ind,i_dir) + Vitesses(ind+1,i_dir) + Vitesses(ind+2,i_dir) + &
                                             (dt * (0.5-bega) * (Forces(ind,i_dir)+Forces(ind+1,i_dir)+Forces(ind+2,i_dir)))
                    enddo
                enddo
            enddo
        enddo

        ! partial of velocity components with respect to xi,eta,zeta
        call elem_part_deriv(m1,m2,m3,hprimex,hprimey,hprimez,Veloc(:,:,:,0),dVx_dxi,dVx_deta,dVx_dzeta)
        call elem_part_deriv(m1,m2,m3,hprimex,hprimey,hprimez,Veloc(:,:,:,1),dVy_dxi,dVy_deta,dVy_dzeta)
        call elem_part_deriv(m1,m2,m3,hprimex,hprimey,hprimez,Veloc(:,:,:,2),dVz_dxi,dVz_deta,dVz_dzeta)

        deallocate(Veloc)

  ! Stress_xx
   ! (Stress_xx)^x
        Elem%slpml%Diagonal_Stress1 (:,:,:,0) = Elem%xpml%DumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,0) + &
            Elem%xpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,0) * dVx_dxi + Elem%sl%Acoeff(:,:,:,1) * dVx_deta + Elem%sl%Acoeff(:,:,:,2) * dVx_dzeta)
   ! (Stress_xx)^y
        Elem%slpml%Diagonal_Stress2 (:,:,:,0) = Elem%xpml%DumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,0) + &
            Elem%xpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,3) * dVy_dxi + Elem%sl%Acoeff(:,:,:,4) * dVy_deta + Elem%sl%Acoeff(:,:,:,5) * dVy_dzeta)
   ! (Stress_xx)^z
        Elem%slpml%Diagonal_Stress3 (:,:,:,0) = Elem%xpml%DumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,0) + &
            Elem%xpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,6) * dVz_dxi + Elem%sl%Acoeff(:,:,:,7) * dVz_deta + Elem%sl%Acoeff(:,:,:,8) * dVz_dzeta)

  ! Stress_yy
        Elem%slpml%Diagonal_Stress1 (:,:,:,1) = Elem%xpml%DumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,1) + &
            Elem%xpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,9) * dVx_dxi + Elem%sl%Acoeff(:,:,:,10) * dVx_deta + Elem%sl%Acoeff(:,:,:,11) * dVx_dzeta)
        Elem%slpml%Diagonal_Stress2 (:,:,:,1) = Elem%xpml%DumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,1) + &
            Elem%xpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,12) * dVy_dxi + Elem%sl%Acoeff(:,:,:,13) * dVy_deta + Elem%sl%Acoeff(:,:,:,14) * dVy_dzeta)
        Elem%slpml%Diagonal_Stress3 (:,:,:,1) = Elem%xpml%DumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,1) + &
            Elem%xpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,6) * dVz_dxi + Elem%sl%Acoeff(:,:,:,7) * dVz_deta + Elem%sl%Acoeff(:,:,:,8) * dVz_dzeta)

  ! Stress_zz
        Elem%slpml%Diagonal_Stress1 (:,:,:,2) = Elem%xpml%DumpSx(:,:,:,0) * Elem%slpml%Diagonal_Stress1 (:,:,:,2) + &
            Elem%xpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,9) * dVx_dxi + Elem%sl%Acoeff(:,:,:,10) * dVx_deta + Elem%sl%Acoeff(:,:,:,11) * dVx_dzeta)
        Elem%slpml%Diagonal_Stress2 (:,:,:,2) = Elem%xpml%DumpSy(:,:,:,0) * Elem%slpml%Diagonal_Stress2 (:,:,:,2) + &
            Elem%xpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,3) * dVy_dxi + Elem%sl%Acoeff(:,:,:,4) * dVy_deta + Elem%sl%Acoeff(:,:,:,5) * dVy_dzeta)
        Elem%slpml%Diagonal_Stress3 (:,:,:,2) = Elem%xpml%DumpSz(:,:,:,0) * Elem%slpml%Diagonal_Stress3 (:,:,:,2) + &
            Elem%xpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,15) * dVz_dxi + Elem%sl%Acoeff(:,:,:,16) * dVz_deta + Elem%sl%Acoeff(:,:,:,17) * dVz_dzeta)

        Elem%slpml%Diagonal_Stress = Elem%slpml%Diagonal_Stress1 + Elem%slpml%Diagonal_Stress2 + Elem%slpml%Diagonal_Stress3

  ! Stress_xy
        Elem%slpml%Residual_Stress1 (:,:,:,0) = Elem%xpml%DumpSx(:,:,:,0) * Elem%slpml%Residual_Stress1 (:,:,:,0) + &
            Elem%xpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,18) * dVy_dxi + Elem%sl%Acoeff(:,:,:,19) * dVy_deta + Elem%sl%Acoeff(:,:,:,20) * dVy_dzeta)
        Elem%slpml%Residual_Stress2 (:,:,:,0) = Elem%xpml%DumpSy(:,:,:,0) * Elem%slpml%Residual_Stress2 (:,:,:,0) + &
            Elem%xpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,21) * dVx_dxi + Elem%sl%Acoeff(:,:,:,22) * dVx_deta + Elem%sl%Acoeff(:,:,:,23) * dVx_dzeta)
        Elem%slpml%Residual_Stress3 (:,:,:,0) = Elem%xpml%DumpSz(:,:,:,0) * Elem%slpml%Residual_Stress3 (:,:,:,0)

  ! Stress_xz
        Elem%slpml%Residual_Stress1 (:,:,:,1) = Elem%xpml%DumpSx(:,:,:,0) * Elem%slpml%Residual_Stress1 (:,:,:,1) + &
            Elem%xpml%DumpSx(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,18) * dVz_dxi + Elem%sl%Acoeff(:,:,:,19) * dVz_deta + Elem%sl%Acoeff(:,:,:,20) * dVz_dzeta)
        Elem%slpml%Residual_Stress2 (:,:,:,1) = Elem%xpml%DumpSy(:,:,:,0) * Elem%slpml%Residual_Stress2 (:,:,:,1)
        Elem%slpml%Residual_Stress3 (:,:,:,1) = Elem%xpml%DumpSz(:,:,:,0) * Elem%slpml%Residual_Stress3 (:,:,:,1) + &
            Elem%xpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,24) * dVx_dxi + Elem%sl%Acoeff(:,:,:,25) * dVx_deta + Elem%sl%Acoeff(:,:,:,26) * dVx_dzeta)

  ! Stress_yz
        Elem%slpml%Residual_Stress1 (:,:,:,2) = Elem%xpml%DumpSx(:,:,:,0) * Elem%slpml%Residual_Stress1 (:,:,:,2)
        Elem%slpml%Residual_Stress2 (:,:,:,2) = Elem%xpml%DumpSy(:,:,:,0) * Elem%slpml%Residual_Stress2 (:,:,:,2) + &
            Elem%xpml%DumpSy(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,21) * dVz_dxi + Elem%sl%Acoeff(:,:,:,22) * dVz_deta + Elem%sl%Acoeff(:,:,:,23) * dVz_dzeta)
        Elem%slpml%Residual_Stress3 (:,:,:,2) = Elem%xpml%DumpSz(:,:,:,0) * Elem%slpml%Residual_Stress3 (:,:,:,2) + &
            Elem%xpml%DumpSz(:,:,:,1) * Dt * (Elem%sl%Acoeff(:,:,:,24) * dVy_dxi + Elem%sl%Acoeff(:,:,:,25) * dVy_deta + Elem%sl%Acoeff(:,:,:,26) * dVy_dzeta)

        Elem%slpml%Residual_Stress = Elem%slpml%Residual_Stress1 + Elem%slpml%Residual_Stress2 + Elem%slpml%Residual_Stress3

        return
    end subroutine Prediction_Elem_PML_Veloc
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------
    subroutine  compute_InternalForces_PML_Elem (Elem, hprimex, hTprimey, htprimez, ngll, Forces)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hTprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
        integer, intent(in) :: ngll
        real, dimension (0:ngll-1,0:2), intent (INOUT) :: Forces

        integer :: m1, m2, m3, n_z, ind, i, j, k
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

        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    ind = Elem%slpml%ISolPml(i,j,k)
                    Forces(ind,:) = Forces(ind,:)+Elem%slpml%Forces1(i,j,k,:)
                    Forces(ind+1,:) = Forces(ind+1,:)+Elem%slpml%Forces2(i,j,k,:)
                    Forces(ind+2,:) = Forces(ind+2,:)+Elem%slpml%Forces3(i,j,k,:)
                enddo
            enddo
        enddo

        return
    end subroutine compute_InternalForces_PML_Elem
    !------------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------------


    subroutine init_element(el)
        type(element), intent(inout) :: el

        el%mat_index=-1
        el%ngllx=0
        el%nglly=0
        el%ngllz=0
        el%PML = .false.
        el%solid = .true.

    end subroutine init_element

end module selement

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
