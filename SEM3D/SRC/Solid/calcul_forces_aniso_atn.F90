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
    subroutine calcul_forces_aniso_atn(dom,bnum,Fox,Foy,Foz,Depla)
        !$acc routine worker
        use champs_solid
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: bnum
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: Sigma

        select case(dom%ngll)
#if GENGLL4
        case(4)
            call calcul_forces_aniso_atn_4(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma)
#endif
#if GENGLL5
        case(5)
            call calcul_forces_aniso_atn_5(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma)
#endif
#if GENGLL6
        case (6)
            call calcul_forces_aniso_atn_6(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma)
#endif
#if GENGLL7
        case (7)
            call calcul_forces_aniso_atn_7(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma)
#endif
#if GENGLL8
        case (8)
            call calcul_forces_aniso_atn_8(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma)
#endif
#if GENGLL9
        case (9)
            call calcul_forces_aniso_atn_9(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma)
#endif
#if GENGLLN
        case default
            call calcul_forces_aniso_atn_n(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla,Sigma)
#endif
        end select
    end subroutine calcul_forces_aniso_atn

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
