!>
!!\file champs_fluidpml.f90
!!\brief Contient la définition du type champs pour un domaine fluide PML
!!
!<

module champs_fluidpml

    use constants
    use mdombase
    implicit none

    type :: champsfluidpml

        !! Fluide PML
        real(fpp), dimension(:,:), allocatable   :: fpml_Forces
        real(fpp), dimension(:,:), allocatable   :: fpml_Phi
        real(fpp), dimension(:,:), allocatable   :: fpml_VelPhi
        real(fpp), dimension(:,:,:), allocatable :: fpml_DumpV

    end type champsfluidpml

    type, extends(dombase) :: domain_fluidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        real(fpp), dimension(:,:), allocatable :: DumpMass
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda, m_Density

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluidpml) :: champs0
        type(champsfluidpml) :: champs1
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_PMLVeloc
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_PMLDumpSx
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_PMLDumpSy
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_PMLDumpSz
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
