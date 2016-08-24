!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Subdomains.f90
!!\brief Calcul des coefficients de Lame.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssubdomains
    use constants
    implicit none

    type Subdomain
        integer          :: dom ! The computation domain SOLID/FLUID/SPML/FPML
        integer          :: material_definition

        !! Numerotation gll
        integer :: NGLL

        !! Definition materiau solide, isotrope
        real(fpp) :: Pspeed, Sspeed, Ddensity
        real(fpp) :: DLambda, DMu
        real(fpp) :: DKappa

        !! Definition materiau solide anisotrope
        ! TODO

        !! ATTENUATION
        real(fpp) :: Qmu, Qpression

        !! PML
        real(fpp), dimension(0:2) :: pml_pos, pml_width
        integer :: npow
        real(fpp) :: Apow

        !! Boundaries for material initialisation from file
        real(fpp), dimension(0:2) :: MinBound, MaxBound, MinBound_Loc, MaxBound_Loc
        character(len=1024), dimension(0:2) :: propFilePath

    end type Subdomain

contains

    !>
    !! \fn subroutine Lame_coefficients (S)
    !! \brief
    !!
    !! \param type (Subdomain) S
    !<
    subroutine Lame_coefficients (S)

        type (Subdomain) :: S

        S%DMu = S%Sspeed**2 * S%Ddensity
        S%DLambda = (S%Pspeed**2 - 2 * S%Sspeed **2 ) * S%Ddensity
        S%DKappa = S%DLambda + 2.*S%DMu /3.

    end subroutine Lame_coefficients

    logical function is_pml(mat)
        type(Subdomain), intent(in) :: mat

        is_pml = .false.
        select case (mat%dom)
        case(DM_SOLID)
            is_pml = .false.
        case(DM_SOLID_PML)
            is_pml = .true.
        case(DM_FLUID)
            is_pml = .false.
        case(DM_FLUID_PML)
            is_pml = .true.
        end select
    end function is_pml

end module ssubdomains

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
