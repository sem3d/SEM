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
#ifdef SINGLEPRECISION
    integer, parameter :: FPP=kind(0E0)
    real(fpp), parameter :: MAX_DOUBLE = 1E15
    real(fpp), parameter :: SMALLFPP=1e-3
#else
    integer, parameter :: FPP=kind(0D0)
    real(fpp), parameter :: MAX_DOUBLE = 1.79769313486231570d+307
    real(fpp), parameter :: SMALLFPP=1e-10
#endif
    ! Constantes mathematique
    ! Les valeurs suivantes et leurs noms sont tirees de math.h
    real(fpp), parameter :: M_E        = 2.7182818284590452354_fpp  ! e
    real(fpp), parameter :: M_LOG2E    = 1.4426950408889634074_fpp  ! log_2 e
    real(fpp), parameter :: M_LOG10E   = 0.43429448190325182765_fpp ! log_10 e
    real(fpp), parameter :: M_LN2      = 0.69314718055994530942_fpp ! log_e 2
    real(fpp), parameter :: M_LN10     = 2.30258509299404568402_fpp ! log_e 10
    real(fpp), parameter :: M_PI       = 3.14159265358979323846_fpp ! pi
    real(fpp), parameter :: M_PI_2     = 1.57079632679489661923_fpp ! pi/2
    real(fpp), parameter :: M_PI_4     = 0.78539816339744830962_fpp ! pi/4
    real(fpp), parameter :: M_1_PI     = 0.31830988618379067154_fpp ! 1/pi
    real(fpp), parameter :: M_2_PI     = 0.63661977236758134308_fpp ! 2/pi
    real(fpp), parameter :: M_2_SQRTPI = 1.12837916709551257390_fpp ! 2/sqrt(pi)
    real(fpp), parameter :: M_SQRT2    = 1.41421356237309504880_fpp ! sqrt(2)
    real(fpp), parameter :: M_SQRT1_2  = 0.70710678118654752440_fpp ! 1/sqrt(2)
    real(fpp), parameter :: M_1_3      = 0.33333333333333333333_fpp ! 1/3
    real(fpp), parameter :: zero=0._fpp,one=1.0_fpp
    real(fpp), parameter :: half=0.5_fpp,two=2.0_fpp,three=3.0_fpp
    real(fpp), parameter :: deps=1e-12_fpp

    ! Constantes physiques

    ! Parametres systemes
    integer     , parameter :: MAX_INT = 2147483647

    ! Parametres algorithmes & code
    integer, parameter :: TIME_INTEG_NEWMARK=0
    integer, parameter :: TIME_INTEG_RK4=1
    integer, parameter :: TIME_INTEG_LDDRK64=2
    integer, parameter :: TIME_INTEG_MIDPOINT=3
    integer, parameter :: TIME_INTEG_MIDPOINT_ITER=4
    integer, parameter :: TIME_INTEG_EXPLICIT=0
    integer, parameter :: TIME_INTEG_SEMI_IMPLICIT=1
    integer, parameter :: TIME_INTEG_IMPLICIT=2
    integer, parameter :: GALERKIN_CONT=0
    integer, parameter :: GALERKIN_DG_STRONG=1
    integer, parameter :: GALERKIN_DG_WEAK=2
    integer, parameter :: GALERKIN_HDG_RP=3
    integer, parameter :: COUPLE_CG_HDG=4
    integer, parameter :: FLUX_NONE=0
    integer, parameter :: FLUX_CENTERED=1
    integer, parameter :: FLUX_GODUNOV=2
    integer, parameter :: FLUX_CUSTOM_LG=3
    integer, parameter :: FLUX_HDG=4
    integer, parameter :: DG_BC_FREE=0
    integer, parameter :: DG_BC_ABS=1
    integer, parameter :: DG_BC_REFL=2
    logical, parameter :: COMPUTE_VHAT=.true.
    logical, parameter :: NOT_COMPUTE_VHAT=.false.

    ! Les fichiers capteurs sont ecrits toutes les NCAPT_CACHE sorties
    integer, parameter :: NCAPT_CACHE=100
    integer, parameter :: CAPT_INTERPOLATED=0
    integer, parameter :: CAPT_NEAREST_NODE=1

    ! Materials
    integer, parameter :: MATERIAL_CONSTANT   = 1
    integer, parameter :: MATERIAL_GRADIENT   = 2 ! DONT USE
    integer, parameter :: MATERIAL_EARTHCHUNK = 3 ! DONT USE
    integer, parameter :: MATERIAL_PREM       = 4 ! DONT USE
    integer, parameter :: MATERIAL_RANDOM     = 5 ! DONT USE
    integer, parameter :: MATERIAL_FILE       = 6

    ! Material definition
    integer, parameter :: MATDEF_VP_VS_RHO     = 0
    integer, parameter :: MATDEF_E_NU_RHO      = 1
    integer, parameter :: MATDEF_LAMBDA_MU_RHO = 2
    integer, parameter :: MATDEF_KAPPA_MU_RHO  = 3
    integer, parameter :: MATDEF_HOOKE_RHO     = 4
    integer, parameter :: MATDEF_NLKP_VS_RHO   = 5
    integer, parameter :: MATDEF_NU_VS_RHO     = 6
    integer, parameter :: MATDEF_VTI_ANISO     = 14
    ! Heterogeneous damping
    integer, parameter :: MATDEF_VP_VS_RHO_D     = 7
    integer, parameter :: MATDEF_E_NU_RHO_D      = 8
    integer, parameter :: MATDEF_LAMBDA_MU_RHO_D = 9
    integer, parameter :: MATDEF_KAPPA_MU_RHO_D  = 10
    integer, parameter :: MATDEF_HOOKE_RHO_D     = 11
    integer, parameter :: MATDEF_NLKP_VS_RHO_D   = 12
    integer, parameter :: MATDEF_NU_VS_RHO_D     = 13

    ! DOMAINS (par ordre de priorite pour les sauvegardes)
    integer, parameter :: DM_FLUID_DG = 6
    integer, parameter :: DM_SOLID_DG = 5
    integer, parameter :: DM_SOLID_CG = 4
    integer, parameter :: DM_SOLID_CG_PML = 2
    integer, parameter :: DM_FLUID_CG = 3
    integer, parameter :: DM_FLUID_CG_PML = 1
    integer, parameter :: DM_MAX = 6

    ! VARIABLES DE SORTIES
    integer, parameter :: OUT_ENERGYP      = 0
    integer, parameter :: OUT_ENERGYS      = 1
    integer, parameter :: OUT_EPS_VOL      = 2
    integer, parameter :: OUT_DEPLA        = 3
    integer, parameter :: OUT_VITESSE      = 4
    integer, parameter :: OUT_ACCEL        = 5
    integer, parameter :: OUT_PRESSION     = 6
    integer, parameter :: OUT_EPS_DEV      = 7
    integer, parameter :: OUT_STRESS_DEV   = 8
    integer, parameter :: OUT_TOTAL_ENERGY = 9
    integer, parameter :: OUT_EPS_DEV_PL   = 10
    integer, parameter :: OUT_DUDX         = 11
    integer, parameter :: OUT_GRAD_LA      = 12
    integer, parameter :: OUT_GRAD_MU      = 13
    integer, parameter :: OUT_LAST=13  ! Numero de la derniere variable
    character(len=10), dimension(0:OUT_LAST) :: OUT_VAR_NAMES = (/ &
        "EnergyP   ", &
        "EnergyS   ", &
        "Eps Vol   ", &
        "Displ     ", &
        "Veloc     ", &
        "Accel     ", &
        "Pressure  ", &
        "Eps Dev   ", &
        "Stress Dev", &
        "Tot_Energy", &
        "Eps Dev Pl", &
        "DUDX      ", &
        "GradLambda", &
        "GradMu    " /)
    integer, parameter, dimension(0:OUT_LAST) :: OUT_VAR_DIMS_3D = (/ 1, 1, 1, 3, 3, 3, 1, 6, 6, 5, 6, 9, 3, 3/)

    integer, parameter :: CPT_INTERP = 0
    integer, parameter :: CPT_ENERGY = 1

    ! TYPE DE CONDITION pour les surfaces
    integer, parameter :: COND_NONE     = 0  ! not assigned/uninitialized
    integer, parameter :: COND_DIRICH   = 1
    integer, parameter :: COND_NEUMANN  = 2

    !METHOD
    integer, parameter :: ISOTROPIC = 1, &
                          SHINOZUKA = 2, &
                          RANDOMIZATION = 3, &
                          FFT = 4
    !Correlation Model
    integer, parameter :: cm_GAUSSIAN = 1, &
                          cm_COS = 2,      &
                          cm_KARMAN = 3
    !First-order Marginal Density
    integer, parameter :: fom_GAUSSIAN = 1, &
                          fom_LOGNORMAL = 2
    !Mesh Mode
    integer, parameter :: msh_AUTO = 1, msh_UNV = 2

    integer, parameter :: SCREEN=6
    integer, parameter :: buf_RF=1024 !Buffer for text
    ! Constants for referencing arrays in CPML code
    integer, parameter :: CPML_ORDER2    = 0
    integer, parameter :: CPML_MIDPOINT1 = 1
    integer, parameter :: CPML_MIDPOINT2 = 2
    integer, parameter :: CPML_ORDER1    = 3
    integer, parameter :: k012 = 0, k021 = 1, k120 = 2
    integer, parameter :: CPML_DIR_X  =0
    integer, parameter :: CPML_DIR_Y  =1
    integer, parameter :: CPML_DIR_Z  =2
    integer, parameter :: CPML_DIR_XY =3
    integer, parameter :: CPML_DIR_XZ =4
    integer, parameter :: CPML_DIR_YZ =5
    integer, parameter :: CPML_DIR_XYZ=6
    integer, parameter :: CMP_ABC=0
    integer, parameter :: CMP_AAC=1
    integer, parameter :: CMP_ABA=2
    integer, parameter :: CMP_ABB=3
    integer, parameter :: CMP_AAA=4

    ! 3D INVARIANT CONSTANTS
    ! projection vector to get isotropic stress
    real(fpp), dimension(0:2,0:2), parameter   :: Mmatrix(0:2,0:2) = one
    real(fpp), dimension(0:5),     parameter   :: Mvector  = (/one,one,one,zero,zero,zero/)
    real(fpp), dimension(0:5),     parameter   :: Avector  = (/one,one,one,two,two,two/)
    real(fpp), dimension(0:5),     parameter   :: A1vector = (/one,one,one,half,half,half/)
    real(fpp), parameter, dimension(0:5,0:5) :: &
        Amatrix  = reshape((/ &
        one , zero, zero, zero, zero, zero,&
        zero, one , zero, zero, zero, zero, &
        zero, zero, one , zero, zero, zero, &
        zero, zero, zero, two , zero, zero, &
        zero, zero, zero, zero, two , zero, &
        zero, zero, zero, zero, zero, two   &
        /), (/6,6/))
    real(fpp), parameter, dimension(0:5,0:5) :: &
        A1matrix = reshape((/ &
        one , zero, zero, zero, zero, zero,&
        zero, one , zero, zero, zero, zero, &
        zero, zero, one , zero, zero, zero, &
        zero, zero, zero, half, zero, zero, &
        zero, zero, zero, zero, half, zero, &
        zero, zero, zero, zero, zero, half  &
        /), (/6,6/))

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
!! vim: set sw=4 ts=4 et tw=80 smartindent :
