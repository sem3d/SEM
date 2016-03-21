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
        character(len=1) :: material_type
        character(len=1) :: initial_material_type
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

        !! RANDOM
        integer :: corrMod
        integer :: assocMat = -1
        integer :: seedStart
        integer  , dimension(:), allocatable :: margiFirst
        integer  , dimension(:), allocatable :: chosenSeed
        real(fpp), dimension(:), allocatable :: varProp
        real(fpp), dimension(:), allocatable :: corrL
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
        select case (mat%material_type)
        case('S')
            is_pml = .false.
        case('P')
            is_pml = .true.
        case('F')
            is_pml = .false.
        case('L')
            is_pml = .true.
        end select
    end function is_pml

    integer function get_domain(mat)
        use constants
        implicit none
        type(Subdomain), intent(in) :: mat

        get_domain = DM_SOLID
        select case (mat%material_type)
        case('S')
            get_domain = DM_SOLID
        case('P')
            get_domain = DM_SOLID_PML
        case('F')
            get_domain = DM_FLUID
        case('L')
            get_domain = DM_FLUID_PML
        case default
            stop "unknown domain"
        end select
        return
    end function get_domain

    logical function is_rand(mat)
        type(Subdomain), intent(in) :: mat

        is_rand = .false.
        select case (mat%initial_material_type)
        case('R')
            is_rand = .true.
        case('S')
            is_rand = .false.
        case('P')
            is_rand = .false.
        case('F')
            is_rand = .false.
        case('L')
            is_rand = .false.
        end select
    end function is_rand

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
