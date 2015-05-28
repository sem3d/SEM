!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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

        ! Attenuation
       real, dimension (:,:,:), allocatable :: Q, Qs, Qp, onemSbeta, onemPbeta, &
           epsilonvol_, &
           epsilondev_xx_,epsilondev_yy_,epsilondev_xy_,epsilondev_xz_,epsilondev_yz_

       real, dimension(:,:,:,:), allocatable :: &
           factor_common_3, alphaval_3,betaval_3,gammaval_3, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_, &
           factor_common_P, alphaval_P,betaval_P,gammaval_P, R_vol_

    end type element_solid

    type :: element_fluid
    end type element_fluid

    type :: element_pml
        real, dimension(:,:,:,:), allocatable :: DumpSx,DumpSy,DumpSz
        real, dimension(:,:,:,:), allocatable :: DumpMass
    end type element_pml

    type :: element_solid_pml
       integer, dimension (:,:,:), allocatable :: ISolPml
       ! TODO move pml related data here
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: Residual_Stress1, Residual_Stress2, Residual_Stress3
       ! FPML
       real, dimension(:,:,:), allocatable :: Isx, Isy, Isz
       real, dimension(:,:,:), allocatable :: Ivx, Ivy, Ivz
       real, dimension(:,:,:,:), allocatable :: Iveloc1, Iveloc2, Iveloc3
       real, dimension(:,:,:,:), allocatable :: I_Diagonal_Stress1, I_Diagonal_Stress2, I_Diagonal_Stress3
       real, dimension(:,:,:,:), allocatable :: I_Residual_Stress1, I_Residual_Stress2
       real, dimension(:,:,:,:), allocatable :: Diagonal_Stress, Residual_Stress
       real, dimension(:,:), allocatable :: Normales, Inv_Normales
    end type element_solid_pml

    type :: element_fluid_pml
        integer, dimension (:,:,:), allocatable :: IFluPml
        real, dimension(:,:,:,:), allocatable :: Veloc
        real, dimension(:,:,:,:), allocatable :: Veloc1,Veloc2,Veloc3
    end type element_fluid_pml

    type :: element

       integer :: mat_index, ngllx, nglly, ngllz
       integer, dimension (:), allocatable :: Control_nodes
       integer, dimension (0:5) :: Near_Faces, Orient_Faces
       integer, dimension (0:11) :: Near_Edges, Orient_Edges
       integer, dimension (0:7) :: Near_Vertices

       integer, dimension (:,:,:), allocatable :: Iglobnum
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
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
