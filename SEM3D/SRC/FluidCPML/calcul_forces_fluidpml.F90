!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module m_calcul_forces_fluidpml ! wrap subroutine in module to get arg type check at build time
    use constants
    use pml
    implicit none
    integer, parameter :: kB012=2, kB021=1, kB120=0
    integer, parameter :: CPML_INTEG = CPML_ORDER2
contains

#include "index.h"
    function isclose(a,b)
        real(fpp) :: a,b
        logical :: isclose
        isclose = abs(a-b)<SMALLFPP
    end function isclose

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
        real(fpp) :: a3b, a4b, a5b
        real(fpp) :: PhiOld, R0, R1
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt, cf0,cf1, cf2

        integer :: n1, n2

        n1 = dom%I1(ee,bnum)
        n2 = dom%I2(ee,bnum)

        PhiOld = dom%PhiOld(ee,i,j,k,bnum)
        ! XXX valable pour ndir=1
        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        ! Update convolution term (implicit midpoint)
        call cpml_compute_coefs(CPML_INTEG, a0, dt, cf0, cf1, cf2)
        R0 = cf0*dom%R1_0(ee,i,j,k,bnum) + cf1*PhiNew + cf2*PhiOld
        dom%R1_0(ee,i,j,k,bnum) = R0
        if (n1==-1 .and. n2==-1) then
            a3b = k0*a0*a0*d0
            R = a3b*R0
        else
            k1 = dom%Kappa_1(i,j,k,n1)
            a1 = dom%Alpha_1(i,j,k,n1)
            d1 = dom%dxi_k_1(i,j,k,n1)
            call cpml_compute_coefs(CPML_INTEG, a1, dt, cf0, cf1, cf2)
            if (.not. isclose(a0,a1)) then
                R1 = cf0*dom%R1_1(i,j,k,n1) + cf1*PhiNew + cf2*PhiOld
            else
                R1 = cf0*dom%R1_1(i,j,k,n1) + cf1*(R0+PhiNew) + cf2*(R0+PhiOld)
            end if
            dom%R1_1(i,j,k,n1) = R1
            if (n2==-1) then
                if (.not. isclose(a0,a1)) then
                    a3b = k0*k1*a0*a0*d0*(d1+a1-a0)/(a1-a0)
                    a4b = k0*k1*a1*a1*d1*(d0+a0-a1)/(a0-a1)
                else
                    a3b = k0*k1*a0*(a0*d0+a0*d1-2*d0*d1)
                    a4b = k0*k1*a0*a0*d0*(a1-a0+d1)
                end if
                R = a3b*R0 + a4b*R1
            else
                call cpml_compute_coefs(CPML_INTEG, a2, dt, cf0, cf1, cf2)
                dom%R1_2(i,j,k,n2) = cf0*dom%R1_2(i,j,k,n2) + cf1*PhiNew + cf2*PhiOld
                k2 = dom%Kappa_2(i,j,k,n2)
                a2 = dom%Alpha_2(i,j,k,n2)
                d2 = dom%dxi_k_2(i,j,k,n2)
                a3b = k0*k1*k2*a0*a0*d0*(d1+a1-a0)*(d2+a2-a0)/((a1-a0)*(a2-a0))
                a4b = k0*k1*k2*a1*a1*d1*(d0+a0-a1)*(d2+a2-a1)/((a0-a1)*(a2-a1))
                a5b = k0*k1*k2*a2*a2*d2*(d0+a0-a2)*(d1+a1-a2)/((a0-a2)*(a1-a2))
                R = a3b*dom%R1_0(ee,i,j,k,bnum) + a4b*dom%R1_1(i,j,k,n1) + a5b*dom%R1_2(i,j,k,n2)
            end if
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
        do r=0,2
            LC(r) = b0(r)*DPhiNew(r) + b1(r)*dom%R2_0(ee,r,i,j,k,bnum)
        end do
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
        integer   :: dim0, dim1, r, i1
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: dt, cf0, cf1, cf2, r0
        real(fpp) :: AA, BB, CC, DD ! common subexpressions
        real(fpp), dimension(0:2) :: b0, b1, b2
        real(fpp), dimension(0:2) :: e0, e1

        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        i1 = dom%I1(ee,bnum)
        k1 = dom%Kappa_1(i,j,k,i1)
        a1 = dom%Alpha_1(i,j,k,i1)
        d1 = dom%dxi_k_1(i,j,k,i1)
        if (k0<0.or.k1<0) then
            write(*,*) "Erreur:", bnum, ee, i,j,k, k0
            stop 1
        endif
        dim0 = dom%D0(ee,bnum)
        dim1 = dom%D1(ee,bnum)
        AA = (a0+d0 - (a1+d1))
        BB = (a0+d0 -  a1)
        CC = (a0-a1)
        DD = (a0 - (a1+d1))
        if (dim0==0.and.dim1==1) then
            ! X=0,Y=1
            ! R2_0 : (e-b0 * dXdx, e-a0 * dXdy, e-a0 * DXdz)
            ! R2_1 : (e-a1 * dXdx, e-b1 * DXdy, e-a1 * dXdz)
            b0(kB120) = k1/k0
            b0(kB021) = k0/k1
            b0(kB012) = k0*k1

            e0(kB120) = a0+d0
            e1(kB120) = a1
            if (.not. isclose(e0(kB120),e1(kB120))) then
                b1(kB120) = -b0(kB120)*d0*AA/BB !b3_120
                b2(kB120) =  b0(kB120)*d1*CC/BB !b1_120
            else
                b1(kB120) = b0(kB120)*(d1+a0-a1)
                b2(kB120) = b0(kB120)*(a0-a1)*d1
            end if

            e0(kB021) = a0
            e1(kB021) = a1+d1
            if (.not. isclose(e0(kB021),e1(kB021))) then
                b1(kB021) = b0(kB021)*d0*CC/DD  !b1_021
                b2(kB021) = b0(kB021)*d1*AA/DD  !b3_021
            else
                b1(kB021) = b0(kB021)*(d0+a1-a0)
                b2(kB021) = b0(kB021)*(a1-a0)*d0
            end if
            e0(kB012) = a0
            e1(kB012) = a1
            if (.not. isclose(e0(kB012),e1(kB012))) then
                b1(kB012) = b0(kB012)*d0*DD/CC  !b1_012
                b2(kB012) = b0(kB012)*d1*BB/CC  !b2_012
            else
                b1(kB012) = b0(kB012)*(d0+d1)   !b1_012
                b2(kB012) = b0(kB012)*d1*d0     !b2_012
            end if
        else if (dim0==0.and.dim1==2) then
            ! X=0,Z=1
            b0(kB120) = k1/k0
            b0(kB021) = k0*k1
            b0(kB012) = k0/k1

            e0(kB120) = a0+d0
            e1(kB120) = a1
            b1(kB120) = -b0(kB120)*d0*AA/BB !b3_120
            b2(kB120) =  b0(kB120)*d1*CC/BB !b2_120

            e0(kB021) = a0
            e1(kB021) = a1
            b1(kB021) = b0(kB021)*d0*DD/CC  !b1_021
            b2(kB021) = b0(kB021)*d1*BB/CC  !b2_021

            e0(kB012) = a0
            e1(kB012) = a1+d1
            b1(kB012) = b0(kB012)*d0*CC/DD  !b1_012
            b2(kB012) = b0(kB012)*d1*AA/DD  !b3_012
        else if (dim0==1.and.dim1==2) then
            ! Y=0,Z=1
            b0(kB120) = k0*k1
            b0(kB021) = k1/k0
            b0(kB012) = k0/k1

            e0(kB120) =  a0
            e1(kB120) =  a1
            b1(kB120) =  b0(kB120)*d0*DD/CC !b1_120
            b2(kB120) =  b0(kB120)*d1*BB/CC !b2_120

            e0(kB021) =  a0+d0
            e1(kB021) =  a1
            b1(kB021) = -b0(kB021)*d0*AA/BB  !b3_021
            b2(kB021) =  b0(kB021)*d1*CC/BB  !b2_021

            e0(kB012) =  a0
            e1(kB012) =  a1+d1
            b1(kB012) =  b0(kB012)*d0*CC/DD  !b2_012
            b2(kB012) =  b0(kB012)*d1*AA/DD  !b3_012
        else
            stop 1
        end if

        ! update convolution terms
        do r=0,2  ! ie 120, 021, 012
            call cpml_compute_coefs(CPML_INTEG, e0(r), dt, cf0, cf1, cf2)
            dom%R2_0(ee,r,i,j,k,bnum) = cf0*dom%R2_0(ee,r,i,j,k,bnum)+cf1*DPhiNew(r)+cf2*DPhi(r)

            call cpml_compute_coefs(CPML_INTEG, e1(r), dt, cf0, cf1, cf2)
            if (.not. isclose(e0(r),e1(r))) then
                dom%R2_1(r,i,j,k,i1) = cf0*dom%R2_1(r,i,j,k,i1)+cf1*DPhiNew(r)+cf2*DPhi(r)
            else
                R0 = dom%R2_0(ee,r,i,j,k,bnum)
                dom%R2_1(r,i,j,k,i1) = cf0*dom%R2_1(r,i,j,k,i1)+cf1*(R0+DPhiNew(r))+cf2*(R0+DPhi(r))
            end if
            !write(*,*) i,j,k, r, e0(r), e1(r), b1(r), b2(r), dom%R2_0(ee,r,i,j,k,bnum), dom%R2_1(r,i,j,k,i1)
        end do
        b1 = 0.0
        b2 = 0.0
        ! We add the terms in Dirac with the (only) convolution term
        do r=0,2 ! ie 120, 021, 012 
            LC(r)=b0(r)*DPhiNew(r)+b1(r)*dom%R2_0(ee,r,i,j,k,bnum)+b2(r)*dom%R2_1(r,i,j,k,i1)
        end do

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
