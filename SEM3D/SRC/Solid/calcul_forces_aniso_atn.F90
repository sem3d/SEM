!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"
#include "gllopt.h"

module m_calcul_forces_aniso_atn ! wrap subroutine in module to get arg type check at build time
    use constants
    implicit none
contains

#define ANISO
#define ATTENUATION
#define PROCNAMEBASE() calcul_forces_aniso_atn_
#define PROCNAMEBASE_ATN() attenuation_aniso_update_

#if GENGLL4
#define NGLLVAL 4
#include "calcul_forces_solid.inc"
#endif

#if GENGLL5
#undef NGLLVAL
#define NGLLVAL 5
#include "calcul_forces_solid.inc"
#endif

#if GENGLL6
#undef NGLLVAL
#define NGLLVAL 6
#include "calcul_forces_solid.inc"
#endif

#if GENGLL7
#undef NGLLVAL
#define NGLLVAL 7
#include "calcul_forces_solid.inc"
#endif

#if GENGLL8
#undef NGLLVAL
#define NGLLVAL 8
#include "calcul_forces_solid.inc"
#endif

#if GENGLL9
#undef NGLLVAL
#define NGLLVAL 9
#include "calcul_forces_solid.inc"
#endif

#if GENGLLN
#undef NGLLVAL
#define NGLL_GEN
#include "calcul_forces_solid.inc"
#endif

end module m_calcul_forces_aniso_atn

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
