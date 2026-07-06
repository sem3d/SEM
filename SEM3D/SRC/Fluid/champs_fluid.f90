!>
!!\file champs_fluid.f90
!!\brief Contient la définition du type champs pour un domaine fluide
!!
!! Anisotropic-density materials (Capdeville & Cance 2015) are handled in this
!! same domain, mirroring dom_solid: m_IDensTensor generalises the scalar
!! m_IDensity = 1/rho to a symmetric tensor rho^{-1}_ij (6 indep comps: 11, 22,
!! 33, 12, 13, 23), selected by the `aniso` flag (allocated only if .true.).
!<

module champs_fluid

    use constants
    use mdombase
    implicit none

    type champsfluid

        !! Fluide
        real(fpp), dimension(:), allocatable :: ForcesFl
        real(fpp), dimension(:), allocatable :: Phi
        real(fpp), dimension(:), allocatable :: VelPhi

    end type champsfluid

    ! Mirror
    type :: time_mirror_fl
        integer :: n_glltot, n_gll
        integer, dimension(:,:,:,:), allocatable :: map
        real(fpp), dimension(:,:), allocatable :: coords
        real(fpp), dimension(:,:), allocatable :: fields
        real(fpp), dimension(:), allocatable :: winfunc
    end type time_mirror_fl

    type, extends(dombase) :: domain_fluid
        ! D'abord, les données membres qui ne sont pas modifiées
        logical :: aniso
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_IDensity ! Inverse of density
        ! Anisotropic-density materials only (allocated iff aniso): rho^{-1}_ij,
        ! index 0=11, 1=22, 2=33, 3=12, 4=13, 5=23
        real(fpp), dimension (:,:,:,:,:,:), allocatable :: m_IDensTensor
        ! Mirror
        !!! GB logical :: use_mirror
        integer :: mirror_type
        type(time_mirror_fl) :: mirror_fl

        ! Champs
        type(champsfluid), dimension(:), allocatable :: champs
    end type domain_fluid

    contains

end module champs_fluid

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
