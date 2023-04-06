!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

!========================================================================
! Ce fichier contient la routine principale de calcul des forces
! solides pour les 4 cas isotrope/anisotrope avec ou sans atténuation.
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

module m_calcul_forces_iso ! wrap subroutine in module to get arg type check at build time
    use constants
    implicit none
contains

#if 1
    subroutine calcul_forces_iso(dom,bnum,Fox,Foy,Foz,Depla,Sigma)
        !$acc routine worker
        use champs_solid
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: bnum
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5), intent(inout) :: Sigma

        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_iso,(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma))
            NGLLDISPATCHCALL_5(calcul_forces_iso,(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma))
            NGLLDISPATCHCALL_6(calcul_forces_iso,(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma))
            NGLLDISPATCHCALL_7(calcul_forces_iso,(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma))
            NGLLDISPATCHCALL_8(calcul_forces_iso,(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma))
            NGLLDISPATCHCALL_9(calcul_forces_iso,(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma))
            NGLLDISPATCHCALL_N(calcul_forces_iso,(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma))
        end select
    end subroutine calcul_forces_iso
#endif

#ifndef TEST_FORCE
#define TEST_FORCE 0
#endif

#undef ATTENUATION
#define PROCNAMEBASE() calcul_forces_iso_

#if TEST_FORCE==1
#if GENGLL4
#undef NGLLVAL
#define NGLLVAL 4
#include "calcul_forces_solid_iso.inc"
#endif

#if GENGLL5
#undef NGLLVAL
#define NGLLVAL 5
#include "calcul_forces_solid_iso.inc"
#endif

#if GENGLL6
#undef NGLLVAL
#define NGLLVAL 6
#include "calcul_forces_solid_iso.inc"
#endif

#if GENGLL7
#undef NGLLVAL
#define NGLLVAL 7
#include "calcul_forces_solid_iso.inc"
#endif

#if GENGLL8
#undef NGLLVAL
#define NGLLVAL 8
#include "calcul_forces_solid_iso.inc"
#endif

#if GENGLL9
#undef NGLLVAL
#define NGLLVAL 9
#include "calcul_forces_solid_iso.inc"
#endif

#if GENGLLN
#undef NGLLVAL
#include "calcul_forces_solid_iso.inc"
#endif

#else

#if GENGLL4
#undef NGLLVAL
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
#include "calcul_forces_solid.inc"
#endif
#endif

end module m_calcul_forces_iso

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
