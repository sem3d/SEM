!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module m_calcul_forces_fluid ! wrap subroutine in module to get arg type check at build time
contains
#include "index.h"
#include "gllopt.h"

#define PROCNAMEBASE() calcul_forces_fl_

#ifndef TEST_FLUID_ACC
#define TEST_FLUID_ACC 0
#endif

#if defined(OPENACC) || TEST_FLUID_ACC==1
#if GENGLL4
#undef NGLLVAL
#define NGLLVAL 4
#include "calcul_forces_fluid_acc.inc"
#endif

#if GENGLL5
#undef NGLLVAL
#define NGLLVAL 5
#include "calcul_forces_fluid_acc.inc"
#endif

#if GENGLL6
#undef NGLLVAL
#define NGLLVAL 6
#include "calcul_forces_fluid_acc.inc"
#endif

#if GENGLL7
#undef NGLLVAL
#define NGLLVAL 7
#include "calcul_forces_fluid_acc.inc"
#endif

#if GENGLL8
#undef NGLLVAL
#define NGLLVAL 8
#include "calcul_forces_fluid_acc.inc"
#endif

#if GENGLL9
#undef NGLLVAL
#define NGLLVAL 9
#include "calcul_forces_fluid_acc.inc"
#endif

#if GENGLLN
#undef NGLLVAL
#include "calcul_forces_fluid_acc.inc"
#endif

#else
#if GENGLL4
#undef NGLLVAL
#define NGLLVAL 4
#include "calcul_forces_fluid.inc"
#endif

#if GENGLL5
#undef NGLLVAL
#define NGLLVAL 5
#include "calcul_forces_fluid.inc"
#endif

#if GENGLL6
#undef NGLLVAL
#define NGLLVAL 6
#include "calcul_forces_fluid.inc"
#endif

#if GENGLL7
#undef NGLLVAL
#define NGLLVAL 7
#include "calcul_forces_fluid.inc"
#endif

#if GENGLL8
#undef NGLLVAL
#define NGLLVAL 8
#include "calcul_forces_fluid.inc"
#endif

#if GENGLL9
#undef NGLLVAL
#define NGLLVAL 9
#include "calcul_forces_fluid.inc"
#endif

#if GENGLLN
#undef NGLLVAL
#include "calcul_forces_fluid.inc"
#endif

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
