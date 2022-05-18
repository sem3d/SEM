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
    type LMC_properties

        ! variables d'écrouissage kinematic et isotrope de Lamaitre et Chaboche
        real(fpp) :: sigma_yld   ! first yielding limit
        real(fpp) :: C_kin       ! variable for kinematic hardening
        real(fpp) :: kapa_kin    ! variable for kinematic hardening
        real(fpp) :: b_iso       ! variable for isotropic hardening
        real(fpp) :: Rinf_iso    ! variable for isotropic hardening

    end type LMC_properties
    !
    type nl_properties
        type(LMC_properties) :: LMC_prop
    end type nl_properties
    !

    type PropertyField
        character(len=1024) :: propFilePath
        character(len=100) :: propName ! name of property and HDF5 group
        real(fpp), dimension(0:2) :: MinBound, MaxBound, step
        ! XXX: need to handle fields of different sizes...
        integer, dimension(0:2) :: NN ! dimension of the grid for this property
        integer, dimension(0:2) :: imin, imax
        real(fpp), dimension(:,:,:), allocatable :: var
        type(nl_properties) :: nl_prop
    end type PropertyField

    type Subdomain
        integer          :: dom ! The computation domain SOLID/FLUID/SPML/FPML
        integer          :: material_definition
        integer          :: deftype
        logical          :: present ! true if an element with this mat exists on this cpu
        !! Numerotation gll
        integer :: NGLL

        !! Definition materiau solide, isotrope
        real(fpp) :: Pspeed, Sspeed, Ddensity
        real(fpp) :: DLambda, DMu
        real(fpp) :: DE, DNu
        real(fpp) :: DKappa
        real(kind=8) :: dt
        !! Definition materiau solide anisotrope
        ! TODO

        !! NONLINEAR LEMAITRE-CHABOCHE
        ! real(fpp) :: DSyld,DCkin,DKkin
        real(fpp) :: DRinf,DBiso
        real(fpp) :: DNlkp
        
        !! ATTENUATION
        real(fpp) :: Qmu, Qpression
        real(fpp) :: Qs, Qp
        !! PML
        real(fpp), dimension(0:2) :: pml_pos, pml_width
        integer :: npow
        real(fpp) :: Apow

        !! Boundaries for material initialisation from file
        real(fpp), dimension(0:2) :: MinBound_Loc, MaxBound_Loc

        ! three variables on a 3D grid used for intializing from fields in files
        ! depending on material_definition we can have
        ! Vp(v1) Vs(v2) Rho(v3)
        ! Lambda(v1) Mu(v2) Rho(v3) ...
        type(PropertyField), dimension(5) :: pf

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
        case(DM_SOLID_CG)
            is_pml = .false.
        case(DM_SOLID_CG_PML)
            is_pml = .true.
        case(DM_FLUID_CG)
            is_pml = .false.
        case(DM_FLUID_CG_PML)
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
