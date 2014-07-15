!>
!! \file constants.F90
!! \brief Ce module contient les definitions des constantes (math, physique, systeme, parametres) du code
!!
!<

MODULE constants
    IMPLICIT none
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

    ! Constantes physiques

    ! Parametres systemes

    ! Parametres algorithmes & code
    integer, parameter :: TIME_INTEG_NEWMARK=0
    integer, parameter :: TIME_INTEG_RK4=1
    integer, parameter :: TIME_INTEG_NEWMARK_DG=2
    integer, parameter :: GALERKIN_CONT=0
    integer, parameter :: GALERKIN_DG_STRONG=1
    integer, parameter :: GALERKIN_DG_WEAK=2
    integer, parameter :: GALERKIN_HDG_RP=3
    integer, parameter :: FLUX_NONE=0
    integer, parameter :: FLUX_CENTERED=1
    integer, parameter :: FLUX_GODUNOV=2
    integer, parameter :: FLUX_CUSTOM_LG=3
    integer, parameter :: FLUX_HDG=4
    integer, parameter :: DG_BC_FREE=0
    integer, parameter :: DG_BC_ABS=1
    integer, parameter :: DG_BC_REFL=2

    ! Les fichiers capteurs sont ecrits toutes les NCAPT_CACHE sorties
    integer, parameter :: NCAPT_CACHE=100


    ! Materials
    integer, parameter :: MATERIAL_CONSTANT   = 1
    integer, parameter :: MATERIAL_GRADIENT   = 2
    integer, parameter :: MATERIAL_EARTHCHUNK = 3


CONTAINS

END MODULE constants
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
