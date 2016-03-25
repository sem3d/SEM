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

module m_calcul_forces ! wrap subroutine in module to get arg type check at build time
    use constants
    implicit none
contains

    subroutine calcul_forces_iso(dom,bnum,Fox,Foy,Foz,Depla)
        use champs_solid
        use deriv3d
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: bnum
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla

        select case(dom%ngll)
        case(4)
            call calcul_forces_iso_4(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case(5)
            call calcul_forces_iso_5(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (6)
            call calcul_forces_iso_6(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (7)
            call calcul_forces_iso_7(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (8)
            call calcul_forces_iso_8(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (9)
            call calcul_forces_iso_9(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case default
            call calcul_forces_iso_n(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        end select
    end subroutine calcul_forces_iso

#undef ATTENUATION
#define NGLLVAL 4
#define PROCNAME calcul_forces_iso_4
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 5
#define PROCNAME calcul_forces_iso_5
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 6
#define PROCNAME calcul_forces_iso_6
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 7
#define PROCNAME calcul_forces_iso_7
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 8
#define PROCNAME calcul_forces_iso_8
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 9
#define PROCNAME calcul_forces_iso_9
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define PROCNAME calcul_forces_iso_n
#include "calcul_forces_solid.inc"

    subroutine calcul_forces_aniso(dom,bnum,Fox,Foy,Foz,Depla)
        use champs_solid
        use deriv3d
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: bnum
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla

        select case(dom%ngll)
        case(4)
            call calcul_forces_aniso_4(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case(5)
            call calcul_forces_aniso_5(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (6)
            call calcul_forces_aniso_6(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (7)
            call calcul_forces_aniso_7(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (8)
            call calcul_forces_aniso_8(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (9)
            call calcul_forces_aniso_9(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case default
            call calcul_forces_aniso_n(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        end select
    end subroutine calcul_forces_aniso

#undef ATTENUATION
#define ANISO
#define NGLLVAL 4
#define PROCNAME calcul_forces_aniso_4
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 5
#define PROCNAME calcul_forces_aniso_5
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 6
#define PROCNAME calcul_forces_aniso_6
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 7
#define PROCNAME calcul_forces_aniso_7
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 8
#define PROCNAME calcul_forces_aniso_8
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 9
#define PROCNAME calcul_forces_aniso_9
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define PROCNAME calcul_forces_aniso_n
#include "calcul_forces_solid.inc"

end module m_calcul_forces

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
