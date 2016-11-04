!>
!!\file champs_fluidpml.f90
!!\brief Contient la définition du type champs pour un domaine fluide
!!
!<

module champs_fluidpml

    use constants
    use ssubdomains
    use mdombase
    implicit none

    type champsfluidpml

        !! Fluide
        real(fpp), dimension(:), allocatable :: ForcesFl
        real(fpp), dimension(:), allocatable :: Phi
        real(fpp), dimension(:), allocatable :: VelPhi

    end type champsfluidpml

    type, extends(dombase) :: domain_fluidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_IDensity ! Inverse of density

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluidpml) :: champs0
        type(champsfluidpml) :: champs1

        ! Element materials
        type(subdomain), dimension (:), pointer :: sSubDomain ! Point to Tdomain%sSubDomain

        ! CPML parameters
        real(fpp) :: cpml_c
        real(fpp) :: cpml_n
        real(fpp) :: cpml_rc
        real(fpp) :: cpml_kappa_0, cpml_kappa_1
        real(fpp) :: alphamax

        ! Integration Rxx
        real(fpp) :: dt

    end type domain_fluidpml

    contains

end module champs_fluidpml

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
