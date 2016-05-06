!>
!!\file champs_solid.f90
!!\brief Contient la définition du type champs pour un domaine solide
!!
!<

module champs_solid

    use constants
    implicit none

    type :: champssolid

        !! Solide
        real(fpp), dimension(:,:), allocatable :: Forces
        real(fpp), dimension(:,:), allocatable :: Depla
        real(fpp), dimension(:,:), allocatable :: Veloc

    end type champssolid

    type :: lmc_param ! LAMAITRE & CHABOCHE NONLINEAR PARAMETERS
        real, dimension (:,:,:,:), allocatable :: m_biso
        real, dimension (:,:,:,:), allocatable :: m_rinf
        real, dimension (:,:,:,:), allocatable :: m_Ckin
        real, dimension (:,:,:,:), allocatable :: m_kkin
        real, dimension (:,:,:,:), allocatable :: m_syld
    end type lmc_param

    type :: nl_parameters  ! STRUCTURE CONTAINING NL PARAMETER SETS
        type(lmc_param), allocatable :: LMC
    end type nl_parameters
    
    !! ATTENTION: voir index.h en ce qui concerne les champs dont les noms commencent par m_
    type domain_solid
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre de gll dans chaque element du domaine
        integer :: ngll

        ! Nombre total de gll du domaine (assembles)
        integer :: nglltot

        ! Nombre d'elements dans le domaine
        integer :: nbelem

        ! Nombre d'elements alloues dans le domaine (>=nbelem)
        integer :: nbelem_alloc
        integer :: nb_chunks ! nbelem_alloc == nb_chunks*CHUNK

        ! Points, poids de gauss et derivees
        real(fpp), dimension (:), allocatable :: GLLc
        real(fpp), dimension (:), allocatable :: GLLw
        real(fpp), dimension (:,:), allocatable :: hprime
        real(fpp), dimension (:,:), allocatable :: hTprime

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real(fpp), dimension(:), allocatable :: MassMat

        logical :: aniso
        real(fpp), dimension (:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Kappa, m_Density
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_Cij

        real(fpp), dimension(:,:,:,:),     allocatable :: m_Jacob
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_InvGrad

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! Index of a gll node within the physical domain
        integer, dimension (:,:,:,:), allocatable :: m_Idom ! Idom copied from element

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolid) :: champs0
        type(champssolid) :: champs1
        ! Attenuation
        integer :: n_sls
        real(fpp) :: dt
        real(fpp), dimension(:,:,:,:),   allocatable :: m_Q
        real(fpp), dimension(:,:,:,:),   allocatable :: m_Qs
        real(fpp), dimension(:,:,:,:),   allocatable :: m_Qp
        real(fpp), dimension(:,:,:,:),   allocatable :: m_epsilonvol
        real(fpp), dimension(:,:,:,:),   allocatable :: m_epsilondev_xx
        real(fpp), dimension(:,:,:,:),   allocatable :: m_epsilondev_yy
        real(fpp), dimension(:,:,:,:),   allocatable :: m_epsilondev_xy
        real(fpp), dimension(:,:,:,:),   allocatable :: m_epsilondev_xz
        real(fpp), dimension(:,:,:,:),   allocatable :: m_epsilondev_yz
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_R_xx
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_R_yy
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_R_xy
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_R_xz
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_R_yz
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_R_vol
        ! Coef RK4
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_omega_tau_s
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_agamma_mu
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_agamma_kappa
        !
        real(fpp), dimension(:,:,:,:),   allocatable :: m_onemSbeta
        real(fpp), dimension(:,:,:,:),   allocatable :: m_onemPbeta
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_factor_common_P
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_alphaval_P
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_betaval_P
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_gammaval_P
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_factor_common_3
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_alphaval_3
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_betaval_3
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_gammaval_3
        
        ! Lemaitre-Chaboche non linear model
        integer :: nl_law
        type(nl_parameters),allocatable :: nl_param
        real(fpp), dimension(:,:,:,:),   allocatable :: m_radius
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_stress 
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_center
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_eps_ep
        real(fpp), dimension(:,:,:,:,:), allocatable :: m_eps_pl

    end type domain_solid

    contains

end module champs_solid

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
