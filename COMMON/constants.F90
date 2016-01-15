!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file constants.F90
!! \brief Ce module contient les definitions des constantes (math, physique, systeme, parametres) du code
!!
!<

MODULE constants
    IMPLICIT none
    ! Precision
    integer, parameter :: FPP=kind(0D0)
    ! Constantes mathematique
    ! Les valeurs suivantes et leurs noms sont tirees de math.h
    real(KIND=8), parameter :: M_E = 2.7182818284590452354D0         ! e
    real(KIND=8), parameter :: M_LOG2E = 1.4426950408889634074D0     ! log_2 e
    real(KIND=8), parameter :: M_LOG10E = 0.43429448190325182765D0   ! log_10 e
    real(KIND=8), parameter :: M_LN2 = 0.69314718055994530942D0      ! log_e 2
    real(KIND=8), parameter :: M_LN10 = 2.30258509299404568402D0     ! log_e 10
    real(KIND=8), parameter :: M_PI = 3.14159265358979323846D0       ! pi
    real(KIND=8), parameter :: M_PI_2 = 1.57079632679489661923D0     ! pi/2
    real(KIND=8), parameter :: M_PI_4 = 0.78539816339744830962D0     ! pi/4
    real(KIND=8), parameter :: M_1_PI = 0.31830988618379067154D0     ! 1/pi
    real(KIND=8), parameter :: M_2_PI = 0.63661977236758134308D0     ! 2/pi
    real(KIND=8), parameter :: M_2_SQRTPI = 1.12837916709551257390D0 ! 2/sqrt(pi)
    real(KIND=8), parameter :: M_SQRT2 = 1.41421356237309504880D0    ! sqrt(2)
    real(KIND=8), parameter :: M_SQRT1_2 = 0.70710678118654752440D0  ! 1/sqrt(2)
    real(KIND=8), parameter :: M_1_3  = 0.33333333333333333333D0 ! 1/3
    real(KIND=8), parameter :: tol_nl = 0.0100000000000000000D0 ! tolerance for nonlinear analysis

    ! Constantes physiques

    ! Parametres systemes
    integer     , parameter :: MAX_INT = 2147483647
    real(KIND=8), parameter :: MAX_DOUBLE = 1.79769313486231570e+307

    ! Parametres algorithmes & code
    integer, parameter :: nProp = 3

    ! Les fichiers capteurs sont ecrits toutes les NCAPT_CACHE sorties
    integer, parameter :: NCAPT_CACHE=100


    ! Materials
    integer, parameter :: MATERIAL_CONSTANT   = 1
    integer, parameter :: MATERIAL_GRADIENT   = 2
    integer, parameter :: MATERIAL_EARTHCHUNK = 3
    integer, parameter :: MATERIAL_PREM       = 4
    integer, parameter :: MATERIAL_MULTIPLE   = 5

    ! DOMAINS (par ordre de priorite pour les sauvegardes)
    integer, parameter :: DM_SOLID = 4
    integer, parameter :: DM_SOLID_PML = 2
    integer, parameter :: DM_FLUID = 3
    integer, parameter :: DM_FLUID_PML = 1

    ! VARIABLES DE SORTIES
    integer, parameter :: OUT_ENERGYP    = 0
    integer, parameter :: OUT_ENERGYS    = 1
    integer, parameter :: OUT_EPS_VOL    = 2
    integer, parameter :: OUT_DEPLA      = 3
    integer, parameter :: OUT_VITESSE    = 4
    integer, parameter :: OUT_ACCEL      = 5
    integer, parameter :: OUT_PRESSION   = 6
    integer, parameter :: OUT_EPS_DEV    = 7
    integer, parameter :: OUT_STRESS_DEV = 8

    ! TYPE DE CONDITION pour les surfaces
    integer, parameter :: COND_NONE     = 0  ! not assigned/uninitialized
    integer, parameter :: COND_DIRICH   = 1
    integer, parameter :: COND_NEUMANN  = 2

    integer, parameter, dimension(0:8) :: OUT_VAR_DIMS_3D = (/ 1, 1, 1, 3, 3, 3, 1, 6, 6 /)
    real(KIND=8), dimension(0:5), parameter :: Miso = M_1_3*(/1, 1, 1, 0, 0, 0/) ! projection vector to get isotropic stress
CONTAINS


    subroutine unused(var)
        real :: var
        var = var
    end subroutine unused
END MODULE constants

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
