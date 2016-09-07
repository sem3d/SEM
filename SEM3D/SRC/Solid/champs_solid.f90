!>
!!\file champs_solid.f90
!!\brief Contient la définition du type champs pour un domaine solide
!!
!<

module champs_solid

    use constants
    use mdombase
    implicit none

    type :: champssolid

        !! Solide
        real(fpp), dimension(:,:), allocatable :: Forces
        real(fpp), dimension(:,:), allocatable :: Depla
        real(fpp), dimension(:,:), allocatable :: Veloc

    end type champssolid

    !! ATTENTION: voir index.h en ce qui concerne les champs dont les noms commencent par m_
    type, extends(dombase) :: domain_solid
        ! D'abord, les données membres qui ne sont pas modifiées
        logical :: aniso
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Kappa, m_Density
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_Cij


        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolid) :: champs0
        type(champssolid) :: champs1
        ! Attenuation
        integer :: n_sls
        real(fpp) :: dt
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_Q
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
