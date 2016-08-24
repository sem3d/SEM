!>
!!\file champs_fluid.f90
!!\brief Contient la définition du type champs pour un domaine fluide
!!
!<

module champs_fluid

    use constants
    use sdombase
    implicit none

    type champsfluid

        !! Fluide
        real(fpp), dimension(:), allocatable :: ForcesFl
        real(fpp), dimension(:), allocatable :: Phi
        real(fpp), dimension(:), allocatable :: VelPhi

    end type champsfluid

    type, extends(dombase) :: domain_fluid
        ! D'abord, les données membres qui ne sont pas modifiées

        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_IDensity ! Inverse of density

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluid) :: champs0
        type(champsfluid) :: champs1
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
