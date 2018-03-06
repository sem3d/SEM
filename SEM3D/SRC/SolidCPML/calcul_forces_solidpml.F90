!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"

module m_calcul_forces_solidpml
    use constants
    use pml
    implicit none

    integer, parameter :: DXX=0
    integer, parameter :: DXY=1
    integer, parameter :: DXZ=2
    integer, parameter :: DYX=3
    integer, parameter :: DYY=4
    integer, parameter :: DYZ=5
    integer, parameter :: DZX=6
    integer, parameter :: DZY=7
    integer, parameter :: DZZ=8
    integer, parameter :: L120_DXX = DXX
    integer, parameter :: L120_DYX = DYX
    integer, parameter :: L120_DZX = DZX
    integer, parameter :: L021_DXY = DXY
    integer, parameter :: L021_DYY = DYY
    integer, parameter :: L021_DZY = DZY
    integer, parameter :: L012_DXZ = DXZ
    integer, parameter :: L012_DYZ = DYZ
    integer, parameter :: L012_DZZ = DZZ
    integer, parameter :: L0_DZZ   =  9
    integer, parameter :: L0_DZY   = 10
    integer, parameter :: L0_DYY   = 11
    integer, parameter :: L0_DYZ   = 12
    integer, parameter :: L1_DZZ   = 13
    integer, parameter :: L1_DZX   = 14
    integer, parameter :: L1_DXZ   = 15
    integer, parameter :: L1_DXX   = 16
    integer, parameter :: L2_DYY   = 17
    integer, parameter :: L2_DYX   = 18
    integer, parameter :: L2_DXY   = 19
    integer, parameter :: L2_DXX   = 20
    integer, parameter :: kB120=0, kB021=1, kB012=2, kB0=3, kB1=4, kB2=5
    integer, parameter :: kB0YY=3,  kB0YZ=4,  kB0ZY=5,  kB0ZZ=6
    integer, parameter :: kB1XX=7,  kB1XZ=8,  kB1ZX=9,  kB1ZZ=10
    integer, parameter :: kB2XX=11, kB2XY=12, kB2YX=13, kB2YY=14
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
    !
    subroutine compute_L_convolution_terms(dom, i, j, k, bnum, ee, U, R)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        real(fpp), intent(in), dimension(0:2) :: U
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), dimension(0:2), intent(out) :: R
        !
        real(fpp) :: a3b, a4b, a5b
        ! dRi/dt
        real(fpp), dimension(0:2) :: dR0, dR1, dR2, R0n, R1n, R2n
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
        call cpml_coefs_midpoint2(a0, dt, cf0, cf1)
        dR0 = cf0*dom%R1_0(ee,:,i,j,k,bnum) + cf1*U
        R0n = dom%R1_0(ee,:,i,j,k,bnum) + 0.5*dR0
        dom%R1_0(ee,:,i,j,k,bnum) = dom%R1_0(ee,:,i,j,k,bnum) + dR0

        if (i1==-1 .and. i2==-1) then
            a3b = k0*a0*a0*d0
            R = a3b*R0n
        else
            k1 = dom%Kappa_1(i,j,k,i1)
            a1 = dom%Alpha_1(i,j,k,i1)
            d1 = dom%dxi_k_1(i,j,k,i1)
            call cpml_coefs_midpoint2(a1, dt, cf0, cf1)
            if (.not. isclose(a0,a1)) then
                dR1 = cf0*dom%R1_1(:,i,j,k,i1) + cf1*U
            else
                ! integration avec R0+0.5dR0/dt
                dR1 = cf0*dom%R1_1(:,i,j,k,i1) + cf1*R0n
            end if
            R1n = dom%R1_1(:,i,j,k,i1) + 0.5*dR1
            dom%R1_1(:,i,j,k,i1) = dom%R1_1(:,i,j,k,i1) + dR1
            if (i2==-1) then
                if (.not. isclose(a0,a1)) then
                    a3b = k0*k1*a0*a0*d0*(d1+a1-a0)/(a1-a0)
                    a4b = k0*k1*a1*a1*d1*(d0+a0-a1)/(a0-a1)
                else
                    a3b = k0*k1*a0*(a0*d0+a0*d1-2*d0*d1)
                    a4b = k0*k1*a0*a0*d0*d1
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
                    dR2 = cf0*dom%R1_2(:,i,j,k,i2) + cf1*U
                    call get_coefs_L_abc(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABA)
                    dR2 = cf0*dom%R1_2(:,i,j,k,i2) + cf1*R0n
                    call get_coefs_L_aba(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAC)
                    dR2 = cf0*dom%R1_2(:,i,j,k,i2) + cf1*U
                    call get_coefs_L_aac(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABB)
                    dR2 = cf0*dom%R1_2(:,i,j,k,i2) + cf1*R1n
                    call get_coefs_L_abb(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAA)
                    dR2 = cf0*dom%R1_2(:,i,j,k,i2) + 2*cf1*R1n
                    call get_coefs_L_aaa(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case default
                    stop 1
                end select
                R2n = dom%R1_2(:,i,j,k,i2) + 0.5*dR2
                dom%R1_2(:,i,j,k,i2) = dom%R1_2(:,i,j,k,i2) + dR2
                R = a3b*R0n + a4b*R1n + a5b*R2n
            end if
        end if
    end subroutine compute_L_convolution_terms

    subroutine compute_convolution_terms(dom, i, j, k, bnum, ee, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1, 0:8) :: DUDV
        real(fpp), intent(out), dimension(0:VCHUNK-1, 0:20) :: LC
        !
        integer :: i1, i2
        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        if (i2==-1) then
            if (i1==-1) then
                call compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DUDV, LC)
            else
                call compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DUDV, LC)
            end if
        else
            call compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DUDV, LC)
        end if
    end subroutine compute_convolution_terms

    subroutine get_coefs_Li_1d(dom, ee, bnum, i, j, k, b0, b1)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: ee, bnum, i, j, k
        real(fpp), dimension(0:5), intent (OUT) :: b0, b1
        !
        real(fpp) :: k0, d0
        integer :: dim0

        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)

        b0(kB0:kB0+2) = 1.
        b1(kB0:kB0+2) = 0.

        dim0 = dom%D0(ee,bnum)
        ! assumes kB1 = kB0+1 and kB2 = kB0+2
        b0(kB0+dim0) = k0
        b1(kB0+dim0) = d0*k0
    end subroutine get_coefs_Li_1d

    subroutine get_coefs_Li_2d(dom, ee, bnum, i, j, k, b0, b1, b2)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: ee, bnum, i, j, k
        real(fpp), dimension(0:2), intent (OUT) :: b0, b1, b2
        !
        real(fpp) :: k0, d0, k1, d1
        integer :: dim0, dim1, i1

        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)

        i1 = dom%I1(ee,bnum)
        k1 = dom%Kappa_1(i,j,k,i1)
        d1 = dom%dxi_k_1(i,j,k,i1)

        b0 = 1.
        b1 = 0.
        b2 = 0.

        dim0 = dom%D0(ee,bnum)
        ! assumes kB1 = kB0+1 and kB2 = kB0+2
        b0(dim0) = k0
        b1(dim0) = d0*k0

        dim1 = dom%D1(ee,bnum)
        ! assumes kB1 = kB0+1 and kB2 = kB0+2
        b0(dim1) = k1
        b2(dim1) = d1*k1
    end subroutine get_coefs_Li_2d

    subroutine get_coefs_Li_3d(dom, ee, bnum, i, j, k, b0, b1, b2, b3)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: ee, bnum, i, j, k
        real(fpp), dimension(0:2), intent (OUT) :: b0, b1, b2, b3
        !
        real(fpp) :: k0, d0, k1, d1, k2, d2
        integer :: i1, i2

        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)

        i1 = dom%I1(ee,bnum)
        k1 = dom%Kappa_1(i,j,k,i1)
        d1 = dom%dxi_k_1(i,j,k,i1)

        i2 = dom%I2(ee,bnum)
        k2 = dom%Kappa_2(i,j,k,i2)
        d2 = dom%dxi_k_2(i,j,k,i2)

        b0 = (/ k0,    k1,    k2    /)
        b1 = (/ d0*k0, 0d0,   0d0   /)
        b2 = (/ 0d0,   d1*k1, 0d0   /)
        b3 = (/ 0d0,   0d0,   d2*k2 /)
    end subroutine get_coefs_Li_3d

    subroutine compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1, 0:8) :: DUDV
        real(fpp), intent(out), dimension(0:VCHUNK-1, 0:20) :: LC
        !
        integer   :: dim0, r
        real(fpp) :: k0, d0, a0, dt, dR, cf0a, cf1a, cf0b, cf1b
        real(fpp), dimension(0:5) :: b0, b1
        real(fpp), dimension(0:8) :: cf0, cf1, R2_0n

        ! Initialize
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        if (k0<0) then
            write(*,*) "Erreur:", bnum, ee, i,j,k, k0
            stop 1
        endif
        b0(:) = 1d0
        b1(:) = 0d0
        ! only two values possible for alpha_0, so we precompute the
        ! integration coefficients
        call cpml_coefs_midpoint2(a0, dt, cf0a, cf1a)
        call cpml_coefs_midpoint2(a0+d0, dt, cf0b, cf1b)
        cf0(:) = cf0a
        cf1(:) = cf1a

        ! Li
        call get_coefs_Li_1d(dom, ee, bnum, i, j, k, b0, b1)

        ! Lijk: Watch out here, b1(kB012) is one of b1 b2 b3 depending on which one is not zero
        dim0 = dom%D0(ee,bnum)
        select case(dim0)
        case(0)
            b0(kB012) = k0
            b0(kB021) = k0
            b0(kB120) = 1./k0
            b1(kB012) = k0*d0  !b1
            b1(kB021) = k0*d0  !b1
            b1(kB120) = -d0/k0 !b3
            cf0(DXX) = cf0b
            cf0(DYX) = cf0b
            cf0(DZX) = cf0b
            cf1(DXX) = cf1b
            cf1(DYX) = cf1b
            cf1(DZX) = cf1b
        case(1)
            b0(kB012) = k0
            b0(kB021) = 1./k0
            b0(kB120) = k0
            b1(kB012) = k0*d0
            b1(kB021) = -d0/k0
            b1(kB120) = k0*d0
            cf0(DXY) = cf0b
            cf0(DYY) = cf0b
            cf0(DZY) = cf0b
            cf1(DXY) = cf1b
            cf1(DYY) = cf1b
            cf1(DZY) = cf1b
        case(2)
            b0(kB012) = 1./k0
            b0(kB021) = k0
            b0(kB120) = k0
            b1(kB012) = -d0/k0
            b1(kB021) = k0*d0
            b1(kB120) = k0*d0
            cf0(DXZ) = cf0b
            cf0(DYZ) = cf0b
            cf0(DZZ) = cf0b
            cf1(DXZ) = cf1b
            cf1(DYZ) = cf1b
            cf1(DZZ) = cf1b
        case default
            stop 1
        end select

        ! update convolution terms
        dt = dom%dt
        do r=0,8
            ! This loop relies on the fact that r=DXX or r=L120_DXX coincide for the right equations
            dR = cf0(r)*dom%R2_0(ee,r,i,j,k,bnum)+cf1(r)*DUDV(ee, r)
            R2_0n(r) = dom%R2_0(ee,r,i,j,k,bnum) + 0.5*dR
            dom%R2_0(ee,r,i,j,k,bnum) = dom%R2_0(ee,r,i,j,k,bnum)+dR
        end do

        ! We add the terms in Dirac with the (only) convolution term
        LC(ee, L120_DXX) = b0(kB120)*DUDV(ee,DXX) + b1(kB120)*R2_0n(DXX)
        LC(ee, L2_DYY  ) = b0(kB2  )*DUDV(ee,DYY) + b1(kB2  )*R2_0n(DYY)
        LC(ee, L1_DZZ  ) = b0(kB1  )*DUDV(ee,DZZ) + b1(kB1  )*R2_0n(DZZ)
        LC(ee, L120_DYX) = b0(kB120)*DUDV(ee,DYX) + b1(kB120)*R2_0n(DYX)
        LC(ee, L2_DYX  ) = b0(kB2  )*DUDV(ee,DYX) + b1(kB2  )*R2_0n(DYX)
        LC(ee, L120_DZX) = b0(kB120)*DUDV(ee,DZX) + b1(kB120)*R2_0n(DZX)
        LC(ee, L1_DZX  ) = b0(kB1  )*DUDV(ee,DZX) + b1(kB1  )*R2_0n(DZX)
        LC(ee, L021_DXY) = b0(kB021)*DUDV(ee,DXY) + b1(kB021)*R2_0n(DXY)
        LC(ee, L2_DXY  ) = b0(kB2  )*DUDV(ee,DXY) + b1(kB2  )*R2_0n(DXY)
        LC(ee, L2_DXX  ) = b0(kB2  )*DUDV(ee,DXX) + b1(kB2  )*R2_0n(DXX)
        LC(ee, L021_DYY) = b0(kB021)*DUDV(ee,DYY) + b1(kB021)*R2_0n(DYY)
        LC(ee, L0_DZZ  ) = b0(kB0  )*DUDV(ee,DZZ) + b1(kB0  )*R2_0n(DZZ)
        LC(ee, L021_DZY) = b0(kB021)*DUDV(ee,DZY) + b1(kB021)*R2_0n(DZY)
        LC(ee, L0_DZY  ) = b0(kB0  )*DUDV(ee,DZY) + b1(kB0  )*R2_0n(DZY)
        LC(ee, L012_DXZ) = b0(kB012)*DUDV(ee,DXZ) + b1(kB012)*R2_0n(DXZ)
        LC(ee, L1_DXZ  ) = b0(kB1  )*DUDV(ee,DXZ) + b1(kB1  )*R2_0n(DXZ)
        LC(ee, L012_DYZ) = b0(kB012)*DUDV(ee,DYZ) + b1(kB012)*R2_0n(DYZ)
        LC(ee, L0_DYZ  ) = b0(kB0  )*DUDV(ee,DYZ) + b1(kB0  )*R2_0n(DYZ)
        LC(ee, L1_DXX  ) = b0(kB1  )*DUDV(ee,DXX) + b1(kB1  )*R2_0n(DXX)
        LC(ee, L0_DYY  ) = b0(kB0  )*DUDV(ee,DYY) + b1(kB0  )*R2_0n(DYY)
        LC(ee, L012_DZZ) = b0(kB012)*DUDV(ee,DZZ) + b1(kB012)*R2_0n(DZZ)
    end subroutine compute_convolution_terms_1d

    ! Compute convolution terms with atn in 2 directions
    subroutine compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1, 0:8) :: DUDV
        real(fpp), intent(out), dimension(0:VCHUNK-1, 0:20) :: LC
        !
        integer :: i1, r
        real(fpp), dimension(0:14) :: b0, b1, b2
        real(fpp), dimension(0:2) :: c0, c1, c2, rr0, rr1
        real(fpp), dimension(0:8) :: e0, e1, R2_0n, R2_1n
        real(fpp) :: dt, cf0, cf1, dR0, dR1
        integer :: dim0, dim1, dimsel
        ! Initialize
        b0 = 1.
        b1 = 0.
        b2 = 0.
        e0 = dom%Alpha_0(ee,i,j,k,bnum)
        i1 = dom%I1(ee,bnum)
        e1 = dom%Alpha_1(i,j,k,i1)
        dim0 = dom%D0(ee,bnum)
        dim1 = dom%D1(ee,bnum)
        ! A single flag for selection on the only 3 possible cases
        ! for dim0/dim1 : (0,1) (0,2), (1,2)
        dimsel=0
        if (dim0==0.and.dim1==2) dimsel=1
        if (dim0==1.and.dim1==2) dimsel=2

        ! Lijk
        call get_coefs_Lijk_2d(dom, ee, bnum, i, j, k, kB120, kB021, kB012, b0, b1, b2, rr0, rr1)
        e0(DXX) = rr0(kB120); e1(DXX) = rr1(kB120)
        e0(DYX) = rr0(kB120); e1(DYX) = rr1(kB120)
        e0(DZX) = rr0(kB120); e1(DZX) = rr1(kB120)
        e0(DXY) = rr0(kB021); e1(DXY) = rr1(kB021)
        e0(DYY) = rr0(kB021); e1(DYY) = rr1(kB021)
        e0(DZY) = rr0(kB021); e1(DZY) = rr1(kB021)
        e0(DXZ) = rr0(kB012); e1(DXZ) = rr1(kB012)
        e0(DYZ) = rr0(kB012); e1(DYZ) = rr1(kB012)
        e0(DZZ) = rr0(kB012); e1(DZZ) = rr1(kB012)

        ! Li
        call get_coefs_Li_2d(dom, ee, bnum, i, j, k, c0, c1, c2)
        b0(kB0YY) = c0(0); b0(kB0YZ) = c0(0); b0(kB0ZY) = c0(0); b0(kB0ZZ) = c0(0)
        b0(kB1XX) = c0(1); b0(kB1XZ) = c0(1); b0(kB1ZX) = c0(1); b0(kB1ZZ) = c0(1)
        b0(kB2XX) = c0(2); b0(kB2XY) = c0(2); b0(kB2YX) = c0(2); b0(kB2YY) = c0(2)
        b1(kB0YY) = c1(0); b1(kB0YZ) = c1(0); b1(kB0ZY) = c1(0); b1(kB0ZZ) = c1(0)
        b1(kB1XX) = c1(1); b1(kB1XZ) = c1(1); b1(kB1ZX) = c1(1); b1(kB1ZZ) = c1(1)
        b1(kB2XX) = c1(2); b1(kB2XY) = c1(2); b1(kB2YX) = c1(2); b1(kB2YY) = c1(2)
        b2(kB0YY) = c2(0); b2(kB0YZ) = c2(0); b2(kB0ZY) = c2(0); b2(kB0ZZ) = c2(0)
        b2(kB1XX) = c2(1); b2(kB1XZ) = c2(1); b2(kB1ZX) = c2(1); b2(kB1ZZ) = c2(1)
        b2(kB2XX) = c2(2); b2(kB2XY) = c2(2); b2(kB2YX) = c2(2); b2(kB2YY) = c2(2)
        ! If double root, need to convole exp(-alpha*t) stored in R2_0 (instead of R2_1 which contains t*exp(-alpha*t))
        ! L1 : handle cases when you have double roots
        if (isclose(e0(DXX), e1(DXX))) then
            if (dimsel==0) then
              b1(kB1XX) = b2(kB1XX); b2(kB1XX) = 0. ! R2_1 -> R2_0
            endif
            b1(kB2XX) = b2(kB2XX); b2(kB2XX) = 0. ! R2_1 -> R2_0
        end if
        if (isclose(e0(DXY), e1(DXY))) then
            b1(kB2XY) = b2(kB2XY); b2(kB2XY) = 0. ! R2_1 -> R2_0
        end if
        if (isclose(e0(DXZ), e1(DXZ))) then
            if (dimsel==0) then
                b1(kB1XZ) = b2(kB1XZ); b2(kB1XZ) = 0. ! R2_1 -> R2_0
            endif
        end if
        if (isclose(e0(DYX), e1(DYX))) then
            b1(kB2YX) = b2(kB2YX); b2(kB2YX) = 0. ! R2_1 -> R2_0
        end if
        if (isclose(e0(DYY), e1(DYY))) then
            b1(kB2YY) = b2(kB2YY); b2(kB2YY) = 0. ! R2_1 -> R2_0
        end if
        !if (isclose(e0(DYZ), e1(DYZ))) then
        !end if
        if (isclose(e0(DZX), e1(DZX))) then
            if (dimsel==0) then
                b1(kB1ZX) = b2(kB1ZX); b2(kB1ZX) = 0. ! R2_1 -> R2_0
            endif
        end if
        !if (isclose(e0(DZY), e1(DZY))) then
        !end if
        if (isclose(e0(DZZ), e1(DZZ))) then
            if (dimsel==0) then
                b1(kB1ZZ) = b2(kB1ZZ); b2(kB1ZZ) = 0. ! R2_1 -> R2_0
            endif
        end if

        ! Convolution term update
        dt = dom%dt
        do r=0,8
            call cpml_coefs_midpoint2(e0(r), dt, cf0, cf1)
            dR0 = cf0*dom%R2_0(ee,r,i,j,k,bnum)+cf1*DUDV(ee, r)
            R2_0n(r) = dom%R2_0(ee,r,i,j,k,bnum) + 0.5*dR0
            dom%R2_0(ee,r,i,j,k,bnum) = dom%R2_0(ee,r,i,j,k,bnum)+dR0

            call cpml_coefs_midpoint2(e1(r), dt, cf0, cf1)
            if (.not. isclose(e0(r),e1(r))) then
                dR1 = cf0*dom%R2_1(r,i,j,k,i1)+cf1*DUDV(ee, r)
            else
                dR1 = cf0*dom%R2_1(r,i,j,k,i1)+cf1*R2_0n(r)
            end if
            R2_1n(r) = dom%R2_1(r,i,j,k,i1) + 0.5*dR1
            dom%R2_1(r,i,j,k,i1) = dom%R2_1(r,i,j,k,i1)+dR1
        end do

        ! Sum of convolution terms at T=nDt
        LC(ee, L120_DXX) = b0(kB120)*DUDV(ee, DXX) + b1(kB120)*R2_0n(DXX) + b2(kB120)*R2_1n(DXX)
        LC(ee, L120_DYX) = b0(kB120)*DUDV(ee, DYX) + b1(kB120)*R2_0n(DYX) + b2(kB120)*R2_1n(DYX)
        LC(ee, L120_DZX) = b0(kB120)*DUDV(ee, DZX) + b1(kB120)*R2_0n(DZX) + b2(kB120)*R2_1n(DZX)
        LC(ee, L021_DXY) = b0(kB021)*DUDV(ee, DXY) + b1(kB021)*R2_0n(DXY) + b2(kB021)*R2_1n(DXY)
        LC(ee, L021_DYY) = b0(kB021)*DUDV(ee, DYY) + b1(kB021)*R2_0n(DYY) + b2(kB021)*R2_1n(DYY)
        LC(ee, L021_DZY) = b0(kB021)*DUDV(ee, DZY) + b1(kB021)*R2_0n(DZY) + b2(kB021)*R2_1n(DZY)
        LC(ee, L012_DXZ) = b0(kB012)*DUDV(ee, DXZ) + b1(kB012)*R2_0n(DXZ) + b2(kB012)*R2_1n(DXZ)
        LC(ee, L012_DYZ) = b0(kB012)*DUDV(ee, DYZ) + b1(kB012)*R2_0n(DYZ) + b2(kB012)*R2_1n(DYZ)
        LC(ee, L012_DZZ) = b0(kB012)*DUDV(ee, DZZ) + b1(kB012)*R2_0n(DZZ) + b2(kB012)*R2_1n(DZZ)
        LC(ee, L0_DYY  ) = b0(kB0YY)*DUDV(ee, DYY) + b1(kB0YY)*R2_0n(DYY) + b2(kB0YY)*R2_1n(DYY)
        LC(ee, L0_DYZ  ) = b0(kB0YZ)*DUDV(ee, DYZ) + b1(kB0YZ)*R2_0n(DYZ) + b2(kB0YZ)*R2_1n(DYZ)
        LC(ee, L0_DZY  ) = b0(kB0ZY)*DUDV(ee, DZY) + b1(kB0ZY)*R2_0n(DZY) + b2(kB0ZY)*R2_1n(DZY)
        LC(ee, L0_DZZ  ) = b0(kB0ZZ)*DUDV(ee, DZZ) + b1(kB0ZZ)*R2_0n(DZZ) + b2(kB0ZZ)*R2_1n(DZZ)
        LC(ee, L1_DXX  ) = b0(kB1XX)*DUDV(ee, DXX) + b1(kB1XX)*R2_0n(DXX) + b2(kB1XX)*R2_1n(DXX)
        LC(ee, L1_DXZ  ) = b0(kB1XZ)*DUDV(ee, DXZ) + b1(kB1XZ)*R2_0n(DXZ) + b2(kB1XZ)*R2_1n(DXZ)
        LC(ee, L1_DZX  ) = b0(kB1ZX)*DUDV(ee, DZX) + b1(kB1ZX)*R2_0n(DZX) + b2(kB1ZX)*R2_1n(DZX)
        LC(ee, L1_DZZ  ) = b0(kB1ZZ)*DUDV(ee, DZZ) + b1(kB1ZZ)*R2_0n(DZZ) + b2(kB1ZZ)*R2_1n(DZZ)
        LC(ee, L2_DXX  ) = b0(kB2XX)*DUDV(ee, DXX) + b1(kB2XX)*R2_0n(DXX) + b2(kB2XX)*R2_1n(DXX)
        LC(ee, L2_DXY  ) = b0(kB2XY)*DUDV(ee, DXY) + b1(kB2XY)*R2_0n(DXY) + b2(kB2XY)*R2_1n(DXY)
        LC(ee, L2_DYX  ) = b0(kB2YX)*DUDV(ee, DYX) + b1(kB2YX)*R2_0n(DYX) + b2(kB2YX)*R2_1n(DYX)
        LC(ee, L2_DYY  ) = b0(kB2YY)*DUDV(ee, DYY) + b1(kB2YY)*R2_0n(DYY) + b2(kB2YY)*R2_1n(DYY)
    end subroutine compute_convolution_terms_2d

    ! Compute convolution terms with atn in 3 directions
    subroutine compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:VCHUNK-1, 0:8) :: DUDV
        real(fpp), intent(out), dimension(0:VCHUNK-1, 0:20) :: LC
        !
        integer   :: r, i1, i2
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt
        real(fpp) :: cf00, cf01, cf10, cf11, cf20, cf21
        real(fpp) :: R0, R1, R2, dR0, dR1, dR2
        real(fpp), dimension(0:14) :: b0, b1, b2, b3
        real(fpp), dimension(0:2) :: c0, c1, c2, c3
        real(fpp), dimension(0:2, 0:8) :: e
        integer, dimension(0:8) :: sel
        real(fpp), dimension(0:8) :: R0n, R1n, R2n
        integer :: sel120, sel021, sel012

        ! Initialize
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
        b0 = 1.
        b1 = 0.
        b2 = 0.
        b3 = 0.
        e = 0.
        sel = CMP_ABC

        ! Lijk
        e(0:2,L120_DXX) = (/ a0+d0, a1, a2 /)
        e(0:2,L120_DYX) = (/ a0+d0, a1, a2 /)
        e(0:2,L120_DZX) = (/ a0+d0, a1, a2 /)
        e(0:2,L021_DXY) = (/ a0, a1+d1, a2 /)
        e(0:2,L021_DYY) = (/ a0, a1+d1, a2 /)
        e(0:2,L021_DZY) = (/ a0, a1+d1, a2 /)
        e(0:2,L012_DXZ) = (/ a0, a1, a2+d2 /)
        e(0:2,L012_DYZ) = (/ a0, a1, a2+d2 /)
        e(0:2,L012_DZZ) = (/ a0, a1, a2+d2 /)
        call get_coefs_Lijk_3d(k1,k2,k0,a1,a2,a0,d1,d2,d0,b0(kB120),b2(kB120),b3(kB120),b1(kB120),sel120)
        call get_coefs_Lijk_3d(k0,k2,k1,a0,a2,a1,d0,d2,d1,b0(kB021),b1(kB021),b3(kB021),b2(kB021),sel021)
        call get_coefs_Lijk_3d(k0,k1,k2,a0,a1,a2,d0,d1,d2,b0(kB012),b1(kB012),b2(kB012),b3(kB012),sel012)
        do r=0,8
            sel(r) = compare_roots(e(0,r), e(1,r), e(2,r))
        enddo

        ! Li
        call get_coefs_Li_3d(dom, ee, bnum, i, j, k, c0, c1, c2, c3)
        b0(kB0YY) = c0(0); b0(kB0YZ) = c0(0); b0(kB0ZY) = c0(0); b0(kB0ZZ) = c0(0)
        b1(kB0YY) = c1(0); b1(kB0YZ) = c1(0); b1(kB0ZY) = c1(0); b1(kB0ZZ) = c1(0)
        b2(kB0YY) = c2(0); b2(kB0YZ) = c2(0); b2(kB0ZY) = c2(0); b2(kB0ZZ) = c2(0)
        b3(kB0YY) = c3(0); b3(kB0YZ) = c3(0); b3(kB0ZY) = c3(0); b3(kB0ZZ) = c3(0)
        b0(kB1XX) = c0(1); b0(kB1XZ) = c0(1); b0(kB1ZX) = c0(1); b0(kB1ZZ) = c0(1)
        b1(kB1XX) = c1(1); b1(kB1XZ) = c1(1); b1(kB1ZX) = c1(1); b1(kB1ZZ) = c1(1)
        b2(kB1XX) = c2(1); b2(kB1XZ) = c2(1); b2(kB1ZX) = c2(1); b2(kB1ZZ) = c2(1)
        b3(kB1XX) = c3(1); b3(kB1XZ) = c3(1); b3(kB1ZX) = c3(1); b3(kB1ZZ) = c3(1)
        b0(kB2XX) = c0(2); b0(kB2XY) = c0(2); b0(kB2YX) = c0(2); b0(kB2YY) = c0(2)
        b1(kB2XX) = c1(2); b1(kB2XY) = c1(2); b1(kB2YX) = c1(2); b1(kB2YY) = c1(2)
        b2(kB2XX) = c2(2); b2(kB2XY) = c2(2); b2(kB2YX) = c2(2); b2(kB2YY) = c2(2)
        b3(kB2XX) = c3(2); b3(kB2XY) = c3(2); b3(kB2YX) = c3(2); b3(kB2YY) = c3(2)
        call update_roots_L1_3d(sel(DXX), kB1XX, b1, b2, b3)
        call update_roots_L2_3d(sel(DXX), kB2XX, b1, b2, b3)
        call update_roots_L2_3d(sel(DXY), kB2XY, b1, b2, b3)
        call update_roots_L1_3d(sel(DXZ), kB1XZ, b1, b2, b3)
        call update_roots_L2_3d(sel(DYX), kB2YX, b1, b2, b3)
        call update_roots_L2_3d(sel(DYY), kB2YY, b1, b2, b3)
        call update_roots_L1_3d(sel(DZX), kB1ZX, b1, b2, b3)
        call update_roots_L1_3d(sel(DZZ), kB1ZZ, b1, b2, b3)

        ! Convolution term update
        dt = dom%dt
        do r=0,8
            R0 = dom%R2_0(ee,r,i,j,k,bnum)
            R1 = dom%R2_1(r,i,j,k,i1)
            R2 = dom%R2_2(r,i,j,k,i2)
            call cpml_coefs_midpoint2(e(0,r), dt, cf00, cf01)
            call cpml_coefs_midpoint2(e(1,r), dt, cf10, cf11)
            call cpml_coefs_midpoint2(e(2,r), dt, cf20, cf21)
            select case(sel(r))
            case (CMP_ABC)
                dR0 = cf00*R0 + cf01*DUDV(ee, r)
                dR1 = cf10*R1 + cf11*DUDV(ee, r)
                dR2 = cf20*R2 + cf21*DUDV(ee, r)
            case (CMP_AAC)
                dR0 = cf00*R0 + cf01*DUDV(ee, r)
                dR1 = cf10*R1 + cf11*(R0+0.5*dR0)
                dR2 = cf20*R2 + cf21*DUDV(ee, r)
            case (CMP_ABA)
                dR0 = cf00*R0 + cf01*DUDV(ee, r)
                dR1 = cf10*R1 + cf11*DUDV(ee, r)
                dR2 = cf20*R2 + cf21*(R0 + 0.5*dR0)
            case (CMP_ABB)
                dR0 = cf00*R0 + cf01*DUDV(ee, r)
                dR1 = cf10*R1 + cf11*DUDV(ee, r)
                dR2 = cf20*R2 + cf21*(R1+0.5*dR1)
            case (CMP_AAA)
                dR0 = cf00*R0 + cf01*DUDV(ee, r)
                dR1 = cf10*R1 + cf11*(R0 + 0.5*dR0)
                dR2 = cf20*R2 + 2*cf21*(R1 + 0.5*dR1)
            case default
                stop 1
            end select
            R0n(r) = R0 + 0.5*dR0
            R1n(r) = R1 + 0.5*dR1
            R2n(r) = R2 + 0.5*dR2
            dom%R2_0(ee,r,i,j,k,bnum) = R0 + dR0
            dom%R2_1(r,i,j,k,i1)      = R1 + dR1
            dom%R2_2(r,i,j,k,i2)      = R2 + dR2
        end do

        ! Convolve
        LC(ee, L120_DXX)=b0(kB120)*DUDV(ee, DXX)+b1(kB120)*R0n(DXX)+b2(kB120)*R1n(DXX)+b3(kB120)*R2n(DXX)
        LC(ee, L120_DYX)=b0(kB120)*DUDV(ee, DYX)+b1(kB120)*R0n(DYX)+b2(kB120)*R1n(DYX)+b3(kB120)*R2n(DYX)
        LC(ee, L120_DZX)=b0(kB120)*DUDV(ee, DZX)+b1(kB120)*R0n(DZX)+b2(kB120)*R1n(DZX)+b3(kB120)*R2n(DZX)
        LC(ee, L021_DXY)=b0(kB021)*DUDV(ee, DXY)+b1(kB021)*R0n(DXY)+b2(kB021)*R1n(DXY)+b3(kB021)*R2n(DXY)
        LC(ee, L021_DYY)=b0(kB021)*DUDV(ee, DYY)+b1(kB021)*R0n(DYY)+b2(kB021)*R1n(DYY)+b3(kB021)*R2n(DYY)
        LC(ee, L021_DZY)=b0(kB021)*DUDV(ee, DZY)+b1(kB021)*R0n(DZY)+b2(kB021)*R1n(DZY)+b3(kB021)*R2n(DZY)
        LC(ee, L012_DXZ)=b0(kB012)*DUDV(ee, DXZ)+b1(kB012)*R0n(DXZ)+b2(kB012)*R1n(DXZ)+b3(kB012)*R2n(DXZ)
        LC(ee, L012_DYZ)=b0(kB012)*DUDV(ee, DYZ)+b1(kB012)*R0n(DYZ)+b2(kB012)*R1n(DYZ)+b3(kB012)*R2n(DYZ)
        LC(ee, L012_DZZ)=b0(kB012)*DUDV(ee, DZZ)+b1(kB012)*R0n(DZZ)+b2(kB012)*R1n(DZZ)+b3(kB012)*R2n(DZZ)
        LC(ee, L0_DYY  )=b0(kB0YY)*DUDV(ee, DYY)+b1(kB0YY)*R0n(DYY)+b2(kB0YY)*R1n(DYY)+b3(kB0YY)*R2n(DYY)
        LC(ee, L0_DYZ  )=b0(kB0YZ)*DUDV(ee, DYZ)+b1(kB0YZ)*R0n(DYZ)+b2(kB0YZ)*R1n(DYZ)+b3(kB0YZ)*R2n(DYZ)
        LC(ee, L0_DZY  )=b0(kB0ZY)*DUDV(ee, DZY)+b1(kB0ZY)*R0n(DZY)+b2(kB0ZY)*R1n(DZY)+b3(kB0ZY)*R2n(DZY)
        LC(ee, L0_DZZ  )=b0(kB0ZZ)*DUDV(ee, DZZ)+b1(kB0ZZ)*R0n(DZZ)+b2(kB0ZZ)*R1n(DZZ)+b3(kB0ZZ)*R2n(DZZ)
        LC(ee, L1_DXX  )=b0(kB1XX)*DUDV(ee, DXX)+b1(kB1XX)*R0n(DXX)+b2(kB1XX)*R1n(DXX)+b3(kB1XX)*R2n(DXX)
        LC(ee, L1_DXZ  )=b0(kB1XZ)*DUDV(ee, DXZ)+b1(kB1XZ)*R0n(DXZ)+b2(kB1XZ)*R1n(DXZ)+b3(kB1XZ)*R2n(DXZ)
        LC(ee, L1_DZX  )=b0(kB1ZX)*DUDV(ee, DZX)+b1(kB1ZX)*R0n(DZX)+b2(kB1ZX)*R1n(DZX)+b3(kB1ZX)*R2n(DZX)
        LC(ee, L1_DZZ  )=b0(kB1ZZ)*DUDV(ee, DZZ)+b1(kB1ZZ)*R0n(DZZ)+b2(kB1ZZ)*R1n(DZZ)+b3(kB1ZZ)*R2n(DZZ)
        LC(ee, L2_DXX  )=b0(kB2XX)*DUDV(ee, DXX)+b1(kB2XX)*R0n(DXX)+b2(kB2XX)*R1n(DXX)+b3(kB2XX)*R2n(DXX)
        LC(ee, L2_DXY  )=b0(kB2XY)*DUDV(ee, DXY)+b1(kB2XY)*R0n(DXY)+b2(kB2XY)*R1n(DXY)+b3(kB2XY)*R2n(DXY)
        LC(ee, L2_DYX  )=b0(kB2YX)*DUDV(ee, DYX)+b1(kB2YX)*R0n(DYX)+b2(kB2YX)*R1n(DYX)+b3(kB2YX)*R2n(DYX)
        LC(ee, L2_DYY  )=b0(kB2YY)*DUDV(ee, DYY)+b1(kB2YY)*R0n(DYY)+b2(kB2YY)*R1n(DYY)+b3(kB2YY)*R2n(DYY)
    end subroutine compute_convolution_terms_3d

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
