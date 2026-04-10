!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
! wrap subroutine in module to get arg type check at build time

module m_calcul_forces_fluid
implicit none
contains
#include "index.h"
#include "optims.h"
#include "loops.h"
#include "gllopt.h"

#define PROCNAMEBASE calcul_forces_fl_

#ifndef TEST_FLUID_ACC
#define TEST_FLUID_ACC USE_ACC_FOR_CPU
#endif

#if GENGLL4
#undef NGLLVAL
#define NGLLVAL 4
#include "calcul_fluid_main.inc"
#endif

#if GENGLL5
#undef NGLLVAL
#define NGLLVAL 5
#include "calcul_fluid_main.inc"
#endif

#if GENGLL6
#undef NGLLVAL
#define NGLLVAL 6
#include "calcul_fluid_main.inc"
#endif

#if GENGLL7
#undef NGLLVAL
#define NGLLVAL 7
#include "calcul_fluid_main.inc"
#endif

#if GENGLL8
#undef NGLLVAL
#define NGLLVAL 8
#include "calcul_fluid_main.inc"
#endif

#if GENGLL9
#undef NGLLVAL
#define NGLLVAL 9
#include "calcul_fluid_main.inc"
#endif

#if GENGLLN
#undef NGLLVAL
#include "calcul_fluid_main.inc"
#endif


end module m_calcul_forces_fluid

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
