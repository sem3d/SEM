!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module m_calcul_forces_fluidpml ! wrap subroutine in module to get arg type check at build time
    use constants
    use pml
    integer, parameter :: kB012=0, kB021=1, kB120=2
    integer, parameter :: CPML_INTEG = CPML_ORDER2
contains

#include "index.h"

    subroutine calcul_forces_fluidpml(dom,ngll,bnum,FFl,Phi)
        use sdomain
        use deriv3d
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: ngll
        integer, intent(in) :: bnum
        !
        integer :: nblocks
        real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: FFl
        real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Phi
        nblocks = dom%nblocks
        select case(ngll)
        case(4)
          call calcul_forces_fluidpml_4(dom,ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_IDensity,FFl,Phi)
        case(5)
            call calcul_forces_fluidpml_5(dom,ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_IDensity,FFl,Phi)
        case (6)
            call calcul_forces_fluidpml_6(dom,ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_IDensity,FFl,Phi)
        case (7)
            call calcul_forces_fluidpml_7(dom,ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_IDensity,FFl,Phi)
        case (8)
            call calcul_forces_fluidpml_8(dom,ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_IDensity,FFl,Phi)
        case (9)
            call calcul_forces_fluidpml_9(dom,ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_IDensity,FFl,Phi)
        case default
            call calcul_forces_fluidpml_n(dom,ngll,nblocks,bnum,dom%hprime,dom%htprime,dom%gllw, &
                dom%m_InvGrad,dom%m_Jacob,dom%m_IDensity,FFl,Phi)
        end select
    end subroutine calcul_forces_fluidpml

    subroutine compute_L_convolution_terms(dom, i, j, k, bnum, ee, PhiNew, R)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: PhiNew
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(out) :: R
        !
        real(fpp) :: a3b
        real(fpp) :: PhiOld
        real(fpp) :: k0, d0, a0, dt, cf0,cf1, cf2
        integer :: n1, n2

        n1 = dom%I1(ee,bnum)
        n2 = dom%I2(ee,bnum)

        PhiOld = dom%PhiOld(ee,i,j,k,bnum)
        ! XXX valable pour ndir=1
        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        if (n1==-1 .and. n2==-1) then
            ! Update convolution term (implicit midpoint)
            call cpml_compute_coefs(CPML_INTEG, a0, dt, cf0, cf1, cf2)
            dom%R1_0(ee,i,j,k,bnum) = cf0*dom%R1_0(ee,i,j,k,bnum) + cf1*PhiNew + cf2*PhiOld
            a3b = k0*a0*a0*d0
            R = a3b*dom%R1_0(ee,i,j,k,bnum)
        else
            R = 0_fpp
        end if
        ! Save PhiOld
        dom%PhiOld(ee,i,j,k,bnum) = PhiNew
    end subroutine compute_L_convolution_terms

    subroutine compute_convolution_terms(dom, i, j, k, bnum, ee, DPhiNew, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:2) :: DPhiNew
        real(fpp), intent(out), dimension(0:2) :: LC
        !
        real(fpp), dimension(0:2) :: DPhiOld
        integer :: n1, n2
        n1 = dom%I1(ee,bnum)
        n2 = dom%I2(ee,bnum)

        DPhiOld = dom%DPhiOld(ee,:,i,j,k,bnum)
        if (n2==-1) then
            if (n1==-1) then
                call compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DPhiNew, DPhiOld, LC)
            else
                call compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DPhiNew, DPhiOld, LC)
            end if
        else
            call compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DPhiNew, DPhiOld, LC)
        end if
        ! Save current value as previous
        dom%DPhiOld(ee,:,i,j,k,bnum) = DPhiNew
    end subroutine compute_convolution_terms

    subroutine compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DPhiNew, DPhi, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:2) :: DPhi, DPhiNew
        real(fpp), intent(out), dimension(0:2) :: LC
        !
        integer   :: dim0, r
        real(fpp) :: k0, d0, a0, dt, cf0, cf1, cf2
        real(fpp), dimension(0:2) :: b0, b1
        real(fpp), dimension(0:2) :: cf

        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        if (k0<0) then
            write(*,*) "Erreur:", bnum, ee, i,j,k, k0
            stop 1
        endif
        b0(:) = 1d0
        dim0 = dom%D0(ee,bnum)
        cf(:) = a0
        ! Watch out here: b1(kB012) is one of b1 b2 b3 depending on which one is not zero
        !
        select case(dim0)
        case(0)
            b0(kB012) = k0
            b0(kB021) = k0
            b0(kB120) = 1./k0
            b1(kB012) = k0*d0  !b1
            b1(kB021) = k0*d0  !b1
            b1(kB120) = -d0/k0 !b3
            cf(0) = a0+d0
        case(1)
            b0(kB012) = k0
            b0(kB021) = 1./k0
            b0(kB120) = k0
            b1(kB012) = k0*d0
            b1(kB021) = -d0/k0
            b1(kB120) = k0*d0
            cf(1) = a0+d0
        case(2)
            b0(kB012) = 1./k0
            b0(kB021) = k0
            b0(kB120) = k0
            b1(kB012) = -d0/k0
            b1(kB021) = k0*d0
            b1(kB120) = k0*d0
            cf(2) = a0+d0
        case default
            stop 1
        end select

        ! update convolution terms
        do r=0,2
            call cpml_compute_coefs(CPML_INTEG, cf(r), dt, cf0, cf1, cf2)
            dom%R2_0(ee,r,i,j,k,bnum) = cf0*dom%R2_0(ee,r,i,j,k,bnum)+cf1*DPhiNew(r)+cf2*DPhi(r)
        end do

        ! We add the terms in Dirac with the (only) convolution term
        LC(0) = b0(kB120)*DPhiNew(0) + b1(kB120)*dom%R2_0(ee,0,i,j,k,bnum)
        LC(1) = b0(kB021)*DPhiNew(1) + b1(kB021)*dom%R2_0(ee,1,i,j,k,bnum)
        LC(2) = b0(kB012)*DPhiNew(2) + b1(kB012)*dom%R2_0(ee,2,i,j,k,bnum)
    end subroutine compute_convolution_terms_1d

    ! Compute convolution terms with atn in 2 directions
    subroutine compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DPhiNew, DPhi, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:2) :: DPhi, DPhiNew
        real(fpp), intent(out), dimension(0:2) :: LC
        !
        stop 1
    end subroutine compute_convolution_terms_2d

    ! Compute convolution terms with atn in 3 directions
    subroutine compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DPhiNew, DPhi, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:2) :: DPhi, DPhiNew
        real(fpp), intent(out), dimension(0:2) :: LC
        !
        !
        stop 1
    end subroutine compute_convolution_terms_3d

#define NGLLVAL 4
#define PROCNAME calcul_forces_fluidpml_4
#include "calcul_forces_fluidpml.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 5
#define PROCNAME calcul_forces_fluidpml_5
#include "calcul_forces_fluidpml.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 6
#define PROCNAME calcul_forces_fluidpml_6
#include "calcul_forces_fluidpml.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 7
#define PROCNAME calcul_forces_fluidpml_7
#include "calcul_forces_fluidpml.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 8
#define PROCNAME calcul_forces_fluidpml_8
#include "calcul_forces_fluidpml.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 9
#define PROCNAME calcul_forces_fluidpml_9
#include "calcul_forces_fluidpml.inc"
#undef NGLLVAL
#undef PROCNAME
#define PROCNAME calcul_forces_fluidpml_n
#include "calcul_forces_fluidpml.inc"

end module m_calcul_forces_fluidpml

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
