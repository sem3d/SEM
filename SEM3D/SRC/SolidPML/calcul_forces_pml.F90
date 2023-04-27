!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP, CNRS
!!

!========================================================================
! Ce fichier contient la routine principale de calcul des PML solides
!
! L'implémentation de ces routines se trouve dans les fichiers .inc
! inclus par les directives #include
!
! Le but de ceci est de permettre au compilateur de "voir" le nombre
! de ngll utilisé réellement, on spécialise les routines pour les ngll
! 4,5,6,7,8,9. Une routine générique est fournie pour les autres
! valeurs.
!
! On gagne 15 à 20% de performance ainsi.
!
!========================================================================

#include "index.h"
#include "gllopt.h"

module m_calcul_forces_pml ! wrap subroutine in module to get arg type check at build time
    use constants
    implicit none
contains

#ifndef TEST_FORCE_PML
#define TEST_FORCE_PML 0
#endif

#define PROCNAMEBASE() calcul_forces_pml_


#if defined(OPENACC) || TEST_FORCE_PML==1
#if GENGLL4
#undef NGLLVAL
#define NGLLVAL 4
#include "calcul_forces_int_pml_acc.inc"
#endif

#if GENGLL5
#undef NGLLVAL
#define NGLLVAL 5
#include "calcul_forces_int_pml_acc.inc"
#endif

#if GENGLL6
#undef NGLLVAL
#define NGLLVAL 6
#include "calcul_forces_int_pml_acc.inc"
#endif

#if GENGLL7
#undef NGLLVAL
#define NGLLVAL 7
#include "calcul_forces_int_pml_acc.inc"
#endif

#if GENGLL8
#undef NGLLVAL
#define NGLLVAL 8
#include "calcul_forces_int_pml_acc.inc"
#endif

#if GENGLL9
#undef NGLLVAL
#define NGLLVAL 9
#include "calcul_forces_int_pml_acc.inc"
#endif

#if GENGLLN
#undef NGLLVAL
#include "calcul_forces_int_pml_acc.inc"
#endif

#else

#if GENGLL4
#undef NGLLVAL
#define NGLLVAL 4
#include "calcul_forces_int_pml.inc"
#endif

#if GENGLL5
#undef NGLLVAL
#define NGLLVAL 5
#include "calcul_forces_int_pml.inc"
#endif

#if GENGLL6
#undef NGLLVAL
#define NGLLVAL 6
#include "calcul_forces_int_pml.inc"
#endif

#if GENGLL7
#undef NGLLVAL
#define NGLLVAL 7
#include "calcul_forces_int_pml.inc"
#endif

#if GENGLL8
#undef NGLLVAL
#define NGLLVAL 8
#include "calcul_forces_int_pml.inc"
#endif

#if GENGLL9
#undef NGLLVAL
#define NGLLVAL 9
#include "calcul_forces_int_pml.inc"
#endif

#if GENGLLN
#undef NGLLVAL
#include "calcul_forces_int_pml.inc"
#endif
#endif

end module m_calcul_forces_pml

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
