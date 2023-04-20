!>
!!\file champs_solid.f90
!!\brief Contient la définition du type champs pour un domaine solide
!!
!<

module champs_solid

    use constants
    use mdombase
    implicit none

    type PWfield
        real(kind=8), dimension(:,:), allocatable :: displ
        real(kind=8), dimension(:,:), allocatable :: veloc
        real(kind=8), dimension(:,:), allocatable :: accel
        logical                                   :: Exist
    end type PWfield

    type :: champssolid
        !! Solide
        real(fpp), dimension(:,:), allocatable :: Depla
        real(fpp), dimension(:,:), allocatable :: Veloc
    end type champssolid

    type :: lmc_param ! LEMAITRE & CHABOCHE NONLINEAR PARAMETERS
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_biso
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_rinf
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_ckin
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_kkin
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_syld
    end type lmc_param

    type :: nl_parameters  ! STRUCTURE CONTAINING NL PARAMETER SETS
        type(lmc_param), allocatable :: LMC
    end type nl_parameters

    ! Mirror
    type :: time_mirror_sl
        integer :: n_glltot, n_gll, n_t
        real(fpp) :: d_t
        integer, dimension(:,:,:,:), allocatable :: map
        real(fpp), dimension(:,:), allocatable :: coords
        real(fpp), dimension(:,:), allocatable :: fields
        real(fpp), dimension(:), allocatable :: winfunc
    end type time_mirror_sl

    !! ATTENTION: voir index.h en ce qui concerne les champs dont les noms commencent par m_
    type, extends(dombase) :: domain_solid
        ! D'abord, les données membres qui ne sont pas modifiées
        logical :: aniso
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Kappa, m_Density
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_Cij
        ! Mirror
        !!! GB logical :: use_mirror
        integer :: mirror_type
        type(time_mirror_sl) :: mirror_sl

        ! pre-allocated temporary mem
        real(fpp), dimension(:,:,:,:,:,:),   allocatable :: Forces
        real(fpp), dimension(:,:,:,:,:,:),   allocatable :: Depla
        real(fpp), dimension(:,:,:,:,:,:),   allocatable :: Veloc
        real(fpp), dimension(:,:,:,:,:,:),   allocatable :: Sigma

        ! Champs
        type(champssolid), dimension(:), allocatable :: champs
        ! Attenuation
        integer :: n_sls

        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_Qs
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_Qp
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_epsilonvol
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_epsilondev_xx
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_epsilondev_yy
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_epsilondev_xy
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_epsilondev_xz
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_epsilondev_yz
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R_xx
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R_yy
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R_xy
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R_xz
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R_yz
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R_vol
        ! Coef RK4
        real(fpp), dimension(:), allocatable :: omega_tau_s
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_agamma_mu
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_agamma_kappa
        !
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_onemSbeta
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_onemPbeta
        ! for plane waves
        type(PWfield)                                  :: PlaneW
        ! Lemaitre-Chaboche non linear model
        integer :: nl_law
        type(nl_parameters),allocatable :: nl_param
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_radius
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_stress
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_center
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_strain
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_plstrain
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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
