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
    implicit none
    !
    type LMC_properties

        ! variables d'écrouissage kinematic et isotrope de Lamaitre et Chaboche
        real :: sigma_yld   ! first yielding limit
        real :: C_kin       ! variable for kinematic hardening
        real :: kapa_kin    ! variable for kinematic hardening
        real :: b_iso       ! variable for isotropic hardening
        real :: Rinf_iso    ! variable for isotropic hardening

    end type LMC_properties
    !
    type nl_properties
        type(LMC_properties) :: LMC_prop
    end type nl_properties
    !
    type Subdomain

        logical :: Filtering, Px, Py, Pz, Left, Forward, Down

        integer :: NGLLx, NGLLy, NGLLz, npow

        !  modif mariotti fevrier 2007 cea
        ! Qmu en plus
        real :: Pspeed, Sspeed, Ddensity, Dt, Apow, freq, DLambda, DMu, Qmu, Q
        real :: DKappa, Qpression
        real, dimension (:), pointer :: GLLcx, GLLwx
        real, dimension (:,:), pointer :: hprimex, hTprimex
        real, dimension (:), pointer :: GLLcy, GLLwy
        real, dimension (:,:), pointer :: hprimey, hTprimey
        real, dimension (:), pointer :: GLLcz, GLLwz
        real, dimension (:,:), pointer :: hprimez, hTprimez

        character(len=1) :: material_type
        integer          :: material_definition
        !Modification to accept random media
        character(len = 15) :: corrMod
        integer             :: assocMat = -1
        integer             :: seedStart
        !integer             :: nElem = 0 !number of elements in each subdomain (by proc) - mesh3d.f90(362)
        !integer            , dimension(:)   , allocatable :: elemList !List of elements in "Tdomain%specel(:)" that belong to this subdomain (by proc)
        !logical            , dimension(:,:) , allocatable :: globCoordMask
        character(len = 30), dimension(:)   , allocatable :: margiFirst
        real               , dimension(:)   , allocatable :: varProp
        integer            , dimension(:)   , allocatable :: chosenSeed
        real               , dimension(:)   , allocatable :: corrL
        real               , dimension(0:2) :: MinBound, MaxBound
        type(nl_properties) :: nl_prop

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
        end select
        return
    end function get_domain

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
