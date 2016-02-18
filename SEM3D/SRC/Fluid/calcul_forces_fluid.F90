!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module m_calcul_forces_fluid ! wrap subroutine in module to get arg type check at build time
contains
#include "index.h"

    subroutine calcul_forces_fluid(dom,ngll,lnum,FFl,Phi)
        use sdomain
        use deriv3d
        implicit none
        type(domain_fluid), intent (INOUT) :: dom
        integer, intent(in) :: ngll
        integer, intent(in) :: lnum
        real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: FFl
        real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Phi

        select case(ngll)
        case(4)
            call calcul_forces_fluid_4(dom,ngll,lnum,FFl,Phi)
        case(5)
            call calcul_forces_fluid_5(dom,ngll,lnum,FFl,Phi)
        case (6)
            call calcul_forces_fluid_6(dom,ngll,lnum,FFl,Phi)
        case (7)
            call calcul_forces_fluid_7(dom,ngll,lnum,FFl,Phi)
        case (8)
            call calcul_forces_fluid_8(dom,ngll,lnum,FFl,Phi)
        case (9)
            call calcul_forces_fluid_9(dom,ngll,lnum,FFl,Phi)
        case default
            call calcul_forces_fluid_n(dom,ngll,lnum,FFl,Phi)
        end select
    end subroutine calcul_forces_fluid

#define NGLLVAL 4
#define PROCNAME calcul_forces_fluid_4
#include "calcul_forces_fluid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 5
#define PROCNAME calcul_forces_fluid_5
#include "calcul_forces_fluid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 6
#define PROCNAME calcul_forces_fluid_6
#include "calcul_forces_fluid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 7
#define PROCNAME calcul_forces_fluid_7
#include "calcul_forces_fluid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 8
#define PROCNAME calcul_forces_fluid_8
#include "calcul_forces_fluid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 9
#define PROCNAME calcul_forces_fluid_9
#include "calcul_forces_fluid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL n
#define PROCNAME calcul_forces_fluid_n
#define NGLL_GEN
#include "calcul_forces_fluid.inc"

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
