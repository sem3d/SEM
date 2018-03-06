!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"

module m_calcul_forces_atn ! wrap subroutine in module to get arg type check at build time
    use constants
    implicit none
contains

#define RK4_attenu_coefs(dt,omega_tau_s,alphaval,betaval,gammaval) \
dt_tau = -dt*omega_tau_s;\
alphaval = 1d0 + dt_tau + 0.5d0 * dt_tau**2 + dt_tau**3 *(1d0/6d0) + dt_tau**4 *(1d0/24.d0); \
betaval  = dt*(0.5d0 + dt_tau * (1d0/3.d0) + dt_tau**2 *(1d0/8d0) + dt_tau**3 *(1d0/24.d0)); \
gammaval = dt*(0.5d0 + dt_tau * (1d0/6.d0) + dt_tau**2 *(1d0/24d0))

    subroutine calcul_forces_iso_atn(dom,bnum,Fox,Foy,Foz,Depla)
        use champs_solid
        use deriv3d
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: bnum
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla

        select case(dom%ngll)
        case(4)
            call calcul_forces_iso_atn_4(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case(5)
            call calcul_forces_iso_atn_5(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (6)
            call calcul_forces_iso_atn_6(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (7)
            call calcul_forces_iso_atn_7(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (8)
            call calcul_forces_iso_atn_8(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (9)
            call calcul_forces_iso_atn_9(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case default
            call calcul_forces_iso_atn_n(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        end select
    end subroutine calcul_forces_iso_atn

#define ATTENUATION
#define NGLLVAL 4
#undef PROCNAME
#define PROCNAME calcul_forces_iso_atn_4
#define PROCNAME_ATN attenuation_update_4
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 5
#define PROCNAME calcul_forces_iso_atn_5
#define PROCNAME_ATN attenuation_update_5
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 6
#define PROCNAME calcul_forces_iso_atn_6
#define PROCNAME_ATN attenuation_update_6
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 7
#define PROCNAME calcul_forces_iso_atn_7
#define PROCNAME_ATN attenuation_update_7
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 8
#define PROCNAME calcul_forces_iso_atn_8
#define PROCNAME_ATN attenuation_update_8
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 9
#define PROCNAME calcul_forces_iso_atn_9
#define PROCNAME_ATN attenuation_update_9
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLL_GEN
#define PROCNAME calcul_forces_iso_atn_n
#define PROCNAME_ATN attenuation_update_n
#include "calcul_forces_solid.inc"

    subroutine calcul_forces_aniso_atn(dom,bnum,Fox,Foy,Foz,Depla)
        use champs_solid
        use deriv3d
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: bnum
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla

        select case(dom%ngll)
        case(4)
            call calcul_forces_aniso_atn_4(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case(5)
            call calcul_forces_aniso_atn_5(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (6)
            call calcul_forces_aniso_atn_6(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (7)
            call calcul_forces_aniso_atn_7(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (8)
            call calcul_forces_aniso_atn_8(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (9)
            call calcul_forces_aniso_atn_9(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case default
            call calcul_forces_aniso_atn_n(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        end select
    end subroutine calcul_forces_aniso_atn

#define ANISO
#define ATTENUATION
#define NGLLVAL 4
#undef PROCNAME
#undef PROCNAME_ATN
#define PROCNAME calcul_forces_aniso_atn_4
#define PROCNAME_ATN attenuation_aniso_update_4
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 5
#define PROCNAME calcul_forces_aniso_atn_5
#define PROCNAME_ATN attenuation_aniso_update_5
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 6
#define PROCNAME calcul_forces_aniso_atn_6
#define PROCNAME_ATN attenuation_aniso_update_6
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 7
#define PROCNAME calcul_forces_aniso_atn_7
#define PROCNAME_ATN attenuation_aniso_update_7
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 8
#define PROCNAME calcul_forces_aniso_atn_8
#define PROCNAME_ATN attenuation_aniso_update_8
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 9
#define PROCNAME calcul_forces_aniso_atn_9
#define PROCNAME_ATN attenuation_aniso_update_9
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLL_GEN
#define PROCNAME calcul_forces_aniso_atn_n
#define PROCNAME_ATN attenuation_aniso_update_n
#include "calcul_forces_solid.inc"


end module m_calcul_forces_atn

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
