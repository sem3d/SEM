!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module m_calcul_forces_fluidpml ! wrap subroutine in module to get arg type check at build time
    use constants
    use pml
    implicit none
    integer, parameter :: kB012=2, kB021=1, kB120=0
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
        real(fpp) :: a3b, a4b, a5b
        real(fpp) :: PhiOld, R0, R1
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt, cf0,cf1, cf2

        integer :: i1, i2, sel

        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        PhiOld = dom%PhiOld(ee,i,j,k,bnum)
        ! XXX valable pour ndir=1
        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        ! Update convolution term (implicit midpoint)
        call cpml_compute_coefs(dom%cpml_integ, a0, dt, cf0, cf1, cf2)
        R0 = cf0*dom%R1_0(ee,i,j,k,bnum) + cf1*PhiNew + cf2*PhiOld
        dom%R1_0(ee,i,j,k,bnum) = R0
        if (i1==-1 .and. i2==-1) then
            a3b = k0*a0*a0*d0
            R = a3b*R0
        else
            k1 = dom%Kappa_1(i,j,k,i1)
            a1 = dom%Alpha_1(i,j,k,i1)
            d1 = dom%dxi_k_1(i,j,k,i1)
            call cpml_compute_coefs(dom%cpml_integ, a1, dt, cf0, cf1, cf2)
            if (.not. isclose(a0,a1)) then
                R1 = cf0*dom%R1_1(i,j,k,i1) + cf1*PhiNew + cf2*PhiOld
            else
                R1 = cf0*dom%R1_1(i,j,k,i1) + (cf1+cf2)*R0
            end if
            dom%R1_1(i,j,k,i1) = R1
            if (i2==-1) then
                if (.not. isclose(a0,a1)) then
                    a3b = k0*k1*a0*a0*d0*(d1+a1-a0)/(a1-a0)
                    a4b = k0*k1*a1*a1*d1*(d0+a0-a1)/(a0-a1)
                else
                    a3b = k0*k1*a0*(a0*d0+a0*d1-2*d0*d1)
                    a4b = k0*k1*a0*a0*d0*(a1-a0+d1)
                end if
                R = a3b*R0 + a4b*R1
            else
                k2 = dom%Kappa_2(i,j,k,i2)
                a2 = dom%Alpha_2(i,j,k,i2)
                d2 = dom%dxi_k_2(i,j,k,i2)
                call cpml_compute_coefs(dom%cpml_integ, a2, dt, cf0, cf1, cf2)
                sel = compare_roots(a0,a1,a2)
                select case(sel)
                case (CMP_ABC)
                    dom%R1_2(i,j,k,i2) = cf0*dom%R1_2(i,j,k,i2) + cf1*PhiNew + cf2*PhiOld
                    call get_coefs_L_abc(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABA)
                    dom%R1_2(i,j,k,i2) = cf0*dom%R1_2(i,j,k,i2) + (cf1+cf2)*R0
                    call get_coefs_L_aba(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAC)
                    dom%R1_2(i,j,k,i2) = cf0*dom%R1_2(i,j,k,i2) + cf1*PhiNew + cf2*PhiOld
                    call get_coefs_L_aac(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABB)
                    dom%R1_2(i,j,k,i2) = cf0*dom%R1_2(i,j,k,i2) + (cf1+cf2)*R1
                    call get_coefs_L_abb(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAA)
                    dom%R1_2(i,j,k,i2) = cf0*dom%R1_2(i,j,k,i2) + 2*(cf1+cf2)*R1
                    call get_coefs_L_aaa(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case default
                    stop 1
                end select
                R = a3b*R0 + a4b*R1 + a5b*dom%R1_2(i,j,k,i2)
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
        integer :: i1, i2
        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        DPhiOld = dom%DPhiOld(ee,:,i,j,k,bnum)
        if (i2==-1) then
            if (i1==-1) then
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
            call cpml_compute_coefs(dom%cpml_integ, cf(r), dt, cf0, cf1, cf2)
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
        integer   :: r, i1
        real(fpp) :: dt, cf0, cf1, cf2, r0
        real(fpp), dimension(0:2) :: b0, b1, b2
        real(fpp), dimension(0:2) :: e0, e1

        call get_coefs_Lijk_2d(dom, ee, bnum, i, j, k, kB120, kB021, kB012, b0, b1, b2, e0, e1)

        ! update convolution terms
        i1 = dom%I1(ee,bnum)
        dt = dom%dt
        do r=0,2  ! ie 120, 021, 012
            call cpml_compute_coefs(dom%cpml_integ, e0(r), dt, cf0, cf1, cf2)
            dom%R2_0(ee,r,i,j,k,bnum) = cf0*dom%R2_0(ee,r,i,j,k,bnum)+cf1*DPhiNew(r)+cf2*DPhi(r)

            call cpml_compute_coefs(dom%cpml_integ, e1(r), dt, cf0, cf1, cf2)
            if (.not. isclose(e0(r),e1(r))) then
                dom%R2_1(r,i,j,k,i1) = cf0*dom%R2_1(r,i,j,k,i1)+cf1*DPhiNew(r)+cf2*DPhi(r)
            else
                R0 = dom%R2_0(ee,r,i,j,k,bnum)
                dom%R2_1(r,i,j,k,i1) = cf0*dom%R2_1(r,i,j,k,i1)+cf1*(R0)+cf2*(R0)
            end if
            !write(*,"(4I2,A3,2F12.3,A3,2F12.3,A3,2E16.8)") i,j,k, r, " E=",e0(r), e1(r)," B=", b1(r), b2(r), " R=", dom%R2_0(ee,r,i,j,k,bnum), dom%R2_1(r,i,j,k,i1)
        end do
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
        integer   :: r, i1, i2
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt
        real(fpp) :: cf00, cf01, cf02
        real(fpp) :: cf10, cf11, cf12
        real(fpp) :: cf20, cf21, cf22
        real(fpp), dimension(0:2) :: b0, b1, b2, b3
        real(fpp), dimension(0:2,0:2) :: e
        real(fpp) :: r0, r1, r2
        integer, dimension(0:2) :: sel

        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        i1 = dom%I1(ee,bnum)
        k1 = dom%Kappa_1(i,j,k,i1)
        a1 = dom%Alpha_1(i,j,k,i1)
        d1 = dom%dxi_k_1(i,j,k,i1)
        i2 = dom%I2(ee,bnum)
        k2 = dom%Kappa_2(i,j,k,i2)
        a2 = dom%Alpha_2(i,j,k,i2)
        d2 = dom%dxi_k_2(i,j,k,i2)
        if (k0<0.or.k1<0.or.k2<0) then
            write(*,*) "Erreur:", bnum, ee, i,j,k, k0
            stop 1
        endif
        !
        e(0:2,kB012) = (/ a0, a1, a2+d2 /)
        e(0:2,kB120) = (/ a1, a2, a0+d0 /)
        e(0:2,kB021) = (/ a0, a2, a1+d1 /)
        call get_coefs_Lijk_3d(k1,k2,k0,a1,a2,a0,d1,d2,d0,b0(kB120),b1(kB120),b2(kB120),b3(kB120),sel(0))
        call get_coefs_Lijk_3d(k0,k2,k1,a0,a2,a1,d0,d2,d1,b0(kB021),b1(kB021),b2(kB021),b3(kB021),sel(1))
        call get_coefs_Lijk_3d(k0,k1,k2,a0,a1,a2,d0,d1,d2,b0(kB012),b1(kB012),b2(kB012),b3(kB012),sel(2))

        do r=0,2  ! ie 120, 021, 012
            R0 = dom%R2_0(ee,r,i,j,k,bnum)
            R1 = dom%R2_1(r,i,j,k,i1)
            R2 = dom%R2_2(r,i,j,k,i2)
            call cpml_compute_coefs(dom%cpml_integ, e(0,r), dt, cf00, cf01, cf02)
            call cpml_compute_coefs(dom%cpml_integ, e(1,r), dt, cf10, cf11, cf12)
            call cpml_compute_coefs(dom%cpml_integ, e(2,r), dt, cf20, cf21, cf22)
            select case(sel(r))
            case (CMP_ABC)
                R0 = cf00*R0 + cf01*DPhiNew(r) + cf02*DPhi(r)
                R1 = cf10*R1 + cf11*DPhiNew(r) + cf12*DPhi(r)
                R2 = cf20*R2 + cf21*DPhiNew(r) + cf22*DPhi(r)
            case (CMP_AAC)
                R0 = cf00*R0 + cf01*DPhiNew(r) + cf02*DPhi(r)
                R1 = cf10*R1 + cf11*R0         + cf12*R0
                R2 = cf20*R2 + cf21*DPhiNew(r) + cf22*DPhi(r)
            case (CMP_ABA)
                R0 = cf00*R0 + cf01*DPhiNew(r) + cf02*DPhi(r)
                R1 = cf10*R1 + cf11*DPhiNew(r) + cf12*DPhi(r)
                R2 = cf20*R2 + cf21*R0         + cf22*R0
            case (CMP_ABB)
                R0 = cf00*R0 + cf01*DPhiNew(r) + cf02*DPhi(r)
                R1 = cf10*R1 + cf11*DPhiNew(r) + cf12*DPhi(r)
                R2 = cf20*R2 + cf21*R1         + cf22*R1
            case (CMP_AAA)
                R0 = cf00*R0 + cf01*DPhiNew(r) + cf02*DPhi(r)
                R1 = cf10*R1 + cf11*R0         + cf12*R0
                R2 = cf20*R2 + cf21*R1         + cf22*R1
            case default
                stop 1
            end select
            dom%R2_0(ee,r,i,j,k,bnum) = R0
            dom%R2_1(r,i,j,k,i1) = R1
            dom%R2_2(r,i,j,k,i2) = R2
        end do

        do r=0,2 ! ie 120, 021, 012
            LC(r) = b0(r)*DPhiNew(r) + &
                    b1(r)*dom%R2_0(ee,r,i,j,k,bnum) + &
                    b2(r)*dom%R2_1(r,i,j,k,i1) + &
                    b3(r)*dom%R2_2(r,i,j,k,i2)
        end do

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
