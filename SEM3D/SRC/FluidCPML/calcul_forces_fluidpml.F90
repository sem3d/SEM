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

    subroutine compute_L_convolution_terms(dom, i, j, k, bnum, ee, Phi, R)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: Phi
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(out) :: R
        !
        real(fpp) :: a3b, a4b, a5b
        real(fpp) :: R0n, R1n, R2n, dR0, dR1, dR2
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt, cf0,cf1

        integer :: i1, i2, sel

        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        ! XXX valable pour ndir=1
        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        ! Update convolution term (implicit midpoint)
        call cpml_coefs_midpoint2(a0, dt, cf0, cf1)
        dR0 = cf0*dom%R1_0(ee,i,j,k,bnum) + cf1*Phi
        R0n = dom%R1_0(ee,i,j,k,bnum) + 0.5*dR0
        dom%R1_0(ee,i,j,k,bnum) = dom%R1_0(ee,i,j,k,bnum) + dR0

        if (i1==-1 .and. i2==-1) then
            a3b = k0*a0*a0*d0
            R = a3b*R0n
        else
            k1 = dom%Kappa_1(i,j,k,i1)
            a1 = dom%Alpha_1(i,j,k,i1)
            d1 = dom%dxi_k_1(i,j,k,i1)
            call cpml_coefs_midpoint2(a1, dt, cf0, cf1)
            if (.not. isclose(a0,a1)) then
                dR1 = cf0*dom%R1_1(i,j,k,i1) + cf1*Phi
            else
                dR1 = cf0*dom%R1_1(i,j,k,i1) + cf1*R0n
            end if
            R1n = dom%R1_1(i,j,k,i1) + 0.5*dR1
            dom%R1_1(i,j,k,i1) = dom%R1_1(i,j,k,i1) + dR1
            if (i2==-1) then
                if (.not. isclose(a0,a1)) then
                    a3b = k0*k1*a0*a0*d0*(d1+a1-a0)/(a1-a0)
                    a4b = k0*k1*a1*a1*d1*(d0+a0-a1)/(a0-a1)
                else
                    a3b = k0*k1*a0*(a0*d0+a0*d1-2*d0*d1)
                    a4b = k0*k1*a0*a0*d0*(a1-a0+d1)
                end if
                R = a3b*R0n + a4b*R1n
            else
                k2 = dom%Kappa_2(i,j,k,i2)
                a2 = dom%Alpha_2(i,j,k,i2)
                d2 = dom%dxi_k_2(i,j,k,i2)
                call cpml_coefs_midpoint2(a2, dt, cf0, cf1)
                sel = compare_roots(a0,a1,a2)
                select case(sel)
                case (CMP_ABC)
                    dR2 = cf0*dom%R1_2(i,j,k,i2) + cf1*Phi
                    call get_coefs_L_abc(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABA)
                    dR2 = cf0*dom%R1_2(i,j,k,i2) + cf1*R0n
                    call get_coefs_L_aba(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAC)
                    dR2 = cf0*dom%R1_2(i,j,k,i2) + cf1*Phi
                    call get_coefs_L_aac(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABB)
                    dR2 = cf0*dom%R1_2(i,j,k,i2) + cf1*R1n
                    call get_coefs_L_abb(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAA)
                    dR2 = cf0*dom%R1_2(i,j,k,i2) + 2*cf1*R1n
                    call get_coefs_L_aaa(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case default
                    stop 1
                end select
                R2n = dom%R1_2(i,j,k,i2) + 0.5*dR2
                dom%R1_2(i,j,k,i2) = dom%R1_2(i,j,k,i2) + dR2
                R = a3b*R0n + a4b*R1n + a5b*R2n
            end if
        end if
    end subroutine compute_L_convolution_terms

    subroutine compute_convolution_terms(dom, i, j, k, bnum, ee, DPhi, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1, 0:2) :: DPhi
        real(fpp), intent(out), dimension(0:VCHUNK-1, 0:2) :: LC
        !
        integer :: i1, i2
        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        if (i2==-1) then
            if (i1==-1) then
                call compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DPhi, LC)
            else
                call compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DPhi, LC)
            end if
        else
            call compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DPhi, LC)
        end if
    end subroutine compute_convolution_terms

    subroutine compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DPhi, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1, 0:2) :: DPhi
        real(fpp), intent(out), dimension(0:VCHUNK-1, 0:2) :: LC
        !
        integer   :: dim0, r
        real(fpp) :: k0, d0, a0, dt, cf0, cf1, R0, R0n, dR0
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
            call cpml_coefs_midpoint2(cf(r), dt, cf0, cf1)
            R0 = dom%R2_0(ee,r,i,j,k,bnum)
            dR0 = cf0*R0+cf1*DPhi(ee,r)
            R0n = R0 + 0.5*dR0
            dom%R2_0(ee,r,i,j,k,bnum) = R0 + dR0
            LC(ee,r) = b0(r)*DPhi(ee,r) + b1(r)*R0n
        end do

    end subroutine compute_convolution_terms_1d

    ! Compute convolution terms with atn in 2 directions
    subroutine compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DPhi, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1,0:2) :: DPhi
        real(fpp), intent(out), dimension(0:VCHUNK-1,0:2) :: LC
        !
        integer   :: r, i1
        real(fpp) :: dt, cf0, cf1, R0, R0n, dR0, R1, R1n, dR1
        real(fpp), dimension(0:2) :: b0, b1, b2
        real(fpp), dimension(0:2) :: e0, e1

        call get_coefs_Lijk_2d(dom, ee, bnum, i, j, k, kB120, kB021, kB012, b0, b1, b2, e0, e1)

        ! update convolution terms
        i1 = dom%I1(ee,bnum)
        dt = dom%dt
        do r=0,2  ! ie 120, 021, 012
            call cpml_coefs_midpoint2(e0(r), dt, cf0, cf1)
            R0 = dom%R2_0(ee,r,i,j,k,bnum)
            dR0 = cf0*R0 + cf1*DPhi(ee,r)
            dom%R2_0(ee,r,i,j,k,bnum) = R0 + dR0
            R0n = R0 + 0.5*dR0
            call cpml_coefs_midpoint2(e1(r), dt, cf0, cf1)
            R1 = dom%R2_1(r,i,j,k,i1)
            if (.not. isclose(e0(r),e1(r))) then
                dR1 = cf0*R1 + cf1*DPhi(ee,r)
            else
                dR1 = cf0*R1 + cf1*R0n
                dom%R2_1(r,i,j,k,i1) = R1 + dR1
            end if
            R1n = R1 + 0.5*dR1
            dom%R2_1(r,i,j,k,i1) = R1 + dR1
            LC(ee,r) = b0(r)*DPhi(ee,r)+b1(r)*R0n+b2(r)*R1n
        end do
    end subroutine compute_convolution_terms_2d

    ! Compute convolution terms with atn in 3 directions
    subroutine compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DPhi, LC)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1,0:2) :: DPhi
        real(fpp), intent(out), dimension(0:VCHUNK-1,0:2) :: LC
        !
        integer   :: r, i1, i2
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt
        real(fpp) :: cf00, cf01
        real(fpp) :: cf10, cf11
        real(fpp) :: cf20, cf21
        real(fpp), dimension(0:2) :: b0, b1, b2, b3
        real(fpp), dimension(0:2,0:2) :: e
        real(fpp) :: R0, R1, R2, dR0, dR1, dR2
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
            call cpml_coefs_midpoint2(e(0,r), dt, cf00, cf01)
            call cpml_coefs_midpoint2(e(1,r), dt, cf10, cf11)
            call cpml_coefs_midpoint2(e(2,r), dt, cf20, cf21)
            select case(sel(r))
            case (CMP_ABC)
                dR0 = cf00*R0 + cf01*DPhi(ee,r)
                dR1 = cf10*R1 + cf11*DPhi(ee,r)
                dR2 = cf20*R2 + cf21*DPhi(ee,r)
            case (CMP_AAC)
                dR0 = cf00*R0 + cf01*DPhi(ee,r)
                dR1 = cf10*R1 + cf11*(R0+0.5*dR0)
                dR2 = cf20*R2 + cf21*DPhi(ee,r)
            case (CMP_ABA)
                dR0 = cf00*R0 + cf01*DPhi(ee,r)
                dR1 = cf10*R1 + cf11*DPhi(ee,r)
                dR2 = cf20*R2 + cf21*(R0+0.5*dR0)
            case (CMP_ABB)
                dR0 = cf00*R0 + cf01*DPhi(ee,r)
                dR1 = cf10*R1 + cf11*DPhi(ee,r)
                dR2 = cf20*R2 + cf21*(R1+0.5*dR1)
            case (CMP_AAA)
                dR0 = cf00*R0 + cf01*DPhi(ee,r)
                dR1 = cf10*R1 + cf11*(R0+0.5*dR0)
                dR2 = cf20*R2 + 2*cf21*(R1+0.5*dR1)
            case default
                stop 1
            end select
            dom%R2_0(ee,r,i,j,k,bnum) = R0 + dR0
            dom%R2_1(r,i,j,k,i1)      = R1 + dR1
            dom%R2_2(r,i,j,k,i2)      = R2 + dR2
            R0 = R0 + 0.5*dR0
            R1 = R1 + 0.5*dR1
            R2 = R2 + 0.5*dR2
            LC(ee,r) = b0(r)*DPhi(ee,r) + b1(r)*R0 + b2(r)*R1 + b3(r)*R2
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
