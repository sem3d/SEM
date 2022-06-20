!>
!!\file champs_solidpml.f90
!!\brief Contient la définition du type champs pour un domaine solide PML
!!
!<

module champs_solidpml

    use constants
    use mdombase
    implicit none

    type :: champssolidpml

        !! Solide PML
        real(fpp), dimension(:,:,:), allocatable :: ForcesPML
        real(fpp), dimension(:,:,:), allocatable :: VelocPML

    end type champssolidpml

    type, extends(dombase) :: domain_solidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        real(fpp), dimension(:,:), allocatable :: DumpMass
        real(fpp), dimension(:,:,:), allocatable :: DumpV
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Density


        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolidpml), dimension(:), allocatable :: champs

        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_Stress1
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_Stress2
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_Stress3
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_PMLDumpSx
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_PMLDumpSy
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_PMLDumpSz
    end type domain_solidpml

    contains

end module champs_solidpml

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
