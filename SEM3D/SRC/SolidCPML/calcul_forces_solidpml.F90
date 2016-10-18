!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"

module m_calcul_forces_solidpml
    use constants
    implicit none
contains

    subroutine calcul_forces_solidpml(dom,bnum,Fox,Foy,Foz,Depla)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: bnum
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1),     intent(inout) :: Fox,Foz,Foy
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in)  :: Depla

        select case(dom%ngll)
        case(4)
            call calcul_forces_solidpml_4(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case(5)
            call calcul_forces_solidpml_5(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (6)
            call calcul_forces_solidpml_6(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (7)
            call calcul_forces_solidpml_7(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (8)
            call calcul_forces_solidpml_8(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case (9)
            call calcul_forces_solidpml_9(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        case default
            call calcul_forces_solidpml_n(dom,dom%ngll,bnum,Fox,Foy,Foz,Depla)
        end select
    end subroutine calcul_forces_solidpml

    subroutine compute_a3a4a5_bar(dom, i, j, k, bnum, ee, a3b, a4b, a5b)
        use sdomain
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent (OUT) :: a3b, a4b, a5b
!        !
!        real(fpp) :: k012,a0,a1,a2,b0,b1,b2
!        k012 = dom%Kappa(ee, 0, i, j, k, bnum)*dom%Kappa(ee, 1, i, j, k, bnum)*dom%Kappa(ee, 2, i, j, k, bnum)
!        a0 = dom%Alpha(ee, 0, i, j, k, bnum)
!        a1 = dom%Alpha(ee, 1, i, j, k, bnum)
!        a2 = dom%Alpha(ee, 2, i, j, k, bnum)
!        b0 = a0 + dom%dxi_k(ee, 0, i, j, k, bnum)
!        b1 = a1 + dom%dxi_k(ee, 1, i, j, k, bnum)
!        b2 = a2 + dom%dxi_k(ee, 2, i, j, k, bnum)
!        a3b = k012*(b0-a0)*(b1-a0)*(b2-a0)/((a1-a0)*(a2-a0)) ! a_012
!        a4b = k012 ! a_
!        a5b = k012
    end subroutine compute_a3a4a5_bar
    !
    subroutine compute_L_convolution_terms(dom, i, j, k, bnum, ee, Rx, Ry, Rz)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(out) :: Rx, Ry, Rz
        !
        real(fpp) :: a3b, a4b, a5b
        call compute_a3a4a5_bar(dom, i, j, k, bnum, ee, a3b, a4b, a5b)
!        Rx = wijk*( &
!            a3b*dom%R(ee,0,0,i,j,k,bnum) + &
!            a4b*dom%R(ee,1,0,i,j,k,bnum) + &
!            a5b*dom%R(ee,2,0,i,j,k,bnum))
!        Ry = wijk*( &
!            a3b*dom%R(ee,0,1,i,j,k,bnum) + &
!            a4b*dom%R(ee,1,1,i,j,k,bnum) + &
!            a5b*dom%R(ee,2,1,i,j,k,bnum))
!        Rz = wijk*( &
!            a3b*dom%R(ee,0,2,i,j,k,bnum) + &
!            a4b*dom%R(ee,1,2,i,j,k,bnum) + &
!            a5b*dom%R(ee,2,2,i,j,k,bnum))
    end subroutine compute_L_convolution_terms
    subroutine compute_convolution_terms(dom, i, j, k, bnum, ee, L0, L1, L2, L012, L021, L120)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(out) :: L0, L1, L2, L012, L021, L120
!        !
!        integer :: ind
!        real(fpp) :: b1b, b2b, b3b
!
!        ind = dom%Idom_(i,j,k,bnum,ee)
!        L0 = dom%Kappa(ee,0,i,j,k,bnum) + dom%Kappa(ee,0,i,j,k,bnum)*dom%dxi_k(ee,0,i,j,k,bnum)*dom%R(0,ind,0) ! 1, x
!        L1 = dom%Kappa(ee,1,i,j,k,bnum) + dom%Kappa(ee,1,i,j,k,bnum)*dom%dxi_k(ee,1,i,j,k,bnum)*dom%R(0,ind,1) ! 1, y
!        L2 = dom%Kappa(ee,2,i,j,k,bnum) + dom%Kappa(ee,2,i,j,k,bnum)*dom%dxi_k(ee,2,i,j,k,bnum)*dom%R(0,ind,2) ! 1, z
!        b1b = 0.
!        b2b = 0.
!        b3b = 0.
!        L120 = L1*L2/L0 + b1b*1. + b2b*1. + b3b*1.
!        b1b = 0.
!        b2b = 0.
!        b3b = 0.
!        L021 = L0*L2/L1 + b1b*1. + b2b*1. + b3b*1.
!        b1b = 0.
!        b2b = 0.
!        b3b = 0.
!        L012 = L0*L1/L2 + b1b*1. + b2b*1. + b3b*1.
    end subroutine compute_convolution_terms

#define NGLLVAL 4
#define PROCNAME calcul_forces_solidpml_4
#include "calcul_forces_solidpml.inc"
#undef NGLLVAL
#undef PROCNAME

#define NGLLVAL 5
#define PROCNAME calcul_forces_solidpml_5
#include "calcul_forces_solidpml.inc"
#undef NGLLVAL
#undef PROCNAME

#define NGLLVAL 6
#define PROCNAME calcul_forces_solidpml_6
#include "calcul_forces_solidpml.inc"
#undef NGLLVAL
#undef PROCNAME

#define NGLLVAL 7
#define PROCNAME calcul_forces_solidpml_7
#include "calcul_forces_solidpml.inc"
#undef NGLLVAL
#undef PROCNAME

#define NGLLVAL 8
#define PROCNAME calcul_forces_solidpml_8
#include "calcul_forces_solidpml.inc"
#undef NGLLVAL
#undef PROCNAME

#define NGLLVAL 9
#define PROCNAME calcul_forces_solidpml_9
#include "calcul_forces_solidpml.inc"
#undef NGLLVAL
#undef PROCNAME

#undef NGLLVAL
#define PROCNAME calcul_forces_solidpml_n
#include "calcul_forces_solidpml.inc"
#undef NGLLVAL
#undef PROCNAME

end module m_calcul_forces_solidpml

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
