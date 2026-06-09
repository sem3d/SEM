!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module m_calcul_forces_fluid_aniso
contains
#include "index.h"

    subroutine calcul_forces_fluid_aniso(dom,ngll,bnum,FFl,P)
        use sdomain
        use deriv3d
        implicit none
        type(domain_fluid_aniso), intent(INOUT) :: dom
        integer, intent(in) :: ngll
        integer, intent(in) :: bnum
        !
        integer :: nblocks
        real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: FFl
        real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(in)  :: P
        nblocks = dom%nblocks
        select case(ngll)
        case(4)
            call calcul_forces_fluid_aniso_4(ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_Kij,FFl,P)
        case(5)
            call calcul_forces_fluid_aniso_5(ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_Kij,FFl,P)
        case(6)
            call calcul_forces_fluid_aniso_6(ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_Kij,FFl,P)
        case(7)
            call calcul_forces_fluid_aniso_7(ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_Kij,FFl,P)
        case(8)
            call calcul_forces_fluid_aniso_8(ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_Kij,FFl,P)
        case(9)
            call calcul_forces_fluid_aniso_9(ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_Kij,FFl,P)
        case default
            call calcul_forces_fluid_aniso_n(ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_Kij,FFl,P)
        end select
    end subroutine calcul_forces_fluid_aniso

#define NGLLVAL 4
#define PROCNAME calcul_forces_fluid_aniso_4
#include "calcul_forces_fluid_aniso.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 5
#define PROCNAME calcul_forces_fluid_aniso_5
#include "calcul_forces_fluid_aniso.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 6
#define PROCNAME calcul_forces_fluid_aniso_6
#include "calcul_forces_fluid_aniso.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 7
#define PROCNAME calcul_forces_fluid_aniso_7
#include "calcul_forces_fluid_aniso.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 8
#define PROCNAME calcul_forces_fluid_aniso_8
#include "calcul_forces_fluid_aniso.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 9
#define PROCNAME calcul_forces_fluid_aniso_9
#include "calcul_forces_fluid_aniso.inc"
#undef NGLLVAL
#undef PROCNAME
#define PROCNAME calcul_forces_fluid_aniso_n
#include "calcul_forces_fluid_aniso.inc"

end module m_calcul_forces_fluid_aniso

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
