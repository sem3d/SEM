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
    subroutine compute_L_convolution_terms(dom, i, j, k, bnum, ee, Unew, R)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        real(fpp), intent(in), dimension(0:2) :: Unew
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), dimension(0:2), intent(out) :: R
        !
        real(fpp) :: a3b, a4b, a5b
        real(fpp), dimension(0:2) :: Uold, R0, R1
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt, cf0,cf1,cf2
        integer :: i1, i2, sel

        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        Uold = dom%Uold(ee,:,i,j,k,bnum)
        ! XXX valable pour ndir=1
        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        call cpml_compute_coefs(dom%cpml_integ, a0, dt, cf0, cf1, cf2)
        R0 = cf0*dom%R1_0(ee,:,i,j,k,bnum) + cf1*Unew(:) + cf2*Uold(:)
        dom%R1_0(ee,:,i,j,k,bnum) = R0

        if (i1==-1 .and. i2==-1) then
            a3b = k0*a0*a0*d0
            R = a3b*R0
        else
            k1 = dom%Kappa_1(i,j,k,i1)
            a1 = dom%Alpha_1(i,j,k,i1)
            d1 = dom%dxi_k_1(i,j,k,i1)
            call cpml_compute_coefs(dom%cpml_integ, a1, dt, cf0, cf1, cf2)
            if (.not. isclose(a0,a1)) then
                R1 = cf0*dom%R1_1(:,i,j,k,i1) + cf1*UNew + cf2*UOld
            else
                R1 = cf0*dom%R1_1(:,i,j,k,i1) + (cf1+cf2)*R0
            end if
            dom%R1_1(:,i,j,k,i1) = R1
            if (i2==-1) then
                if (.not. isclose(a0,a1)) then
                    a3b = k0*k1*a0*a0*d0*(d1+a1-a0)/(a1-a0)
                    a4b = k0*k1*a1*a1*d1*(d0+a0-a1)/(a0-a1)
                else
                    a3b = k0*k1*a0*(a0*d0+a0*d1-2*d0*d1)
                    a4b = k0*k1*a0*a0*d0*d1
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
                    dom%R1_2(:,i,j,k,i2) = cf0*dom%R1_2(:,i,j,k,i2) + cf1*UNew + cf2*UOld
                    call get_coefs_L_abc(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABA)
                    dom%R1_2(:,i,j,k,i2) = cf0*dom%R1_2(:,i,j,k,i2) + (cf1+cf2)*R0
                    call get_coefs_L_aba(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAC)
                    dom%R1_2(:,i,j,k,i2) = cf0*dom%R1_2(:,i,j,k,i2) + cf1*UNew + cf2*UOld
                    call get_coefs_L_aac(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_ABB)
                    dom%R1_2(:,i,j,k,i2) = cf0*dom%R1_2(:,i,j,k,i2) + (cf1+cf2)*R1
                    call get_coefs_L_abb(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case (CMP_AAA)
                    dom%R1_2(:,i,j,k,i2) = cf0*dom%R1_2(:,i,j,k,i2) + 2*(cf1+cf2)*R1
                    call get_coefs_L_aaa(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
                case default
                    stop 1
                end select
                R = a3b*R0 + a4b*R1 + a5b*dom%R1_2(:,i,j,k,i2)
            end if
        end if

        ! Save Uold
        dom%Uold(ee,:,i,j,k,bnum) = Unew
    end subroutine compute_L_convolution_terms

    subroutine compute_convolution_terms(dom, i, j, k, bnum, ee, DUDVnew, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:8) :: DUDVnew
        real(fpp), intent(out), dimension(0:20) :: LC
        !
        real(fpp), dimension(0:8) :: DUDV
        integer :: i1, i2
        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        DUDV = dom%DUDVold(ee,:,i,j,k,bnum) !0.5d0*DUDVnew+0.5d0*dom%DUDVold(ee,:,i,j,k,bnum)
        if (i2==-1) then
            if (i1==-1) then
                call compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DUDVnew, DUDV, LC)
            else
                call compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DUDVnew, DUDV, LC)
            end if
        else
            call compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DUDVnew, DUDV, LC)
        end if
        ! Save current value as previous
        dom%DUDVold(ee,:,i,j,k,bnum) = DUDVnew
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

    subroutine compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DUDVn, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:8) :: DUDV, DUDVn
        real(fpp), intent(out), dimension(0:20) :: LC
        !
        integer   :: dim0, r
        real(fpp) :: k0, d0, a0, dt, cf0, cf1, cf2
        real(fpp), dimension(0:5) :: b0, b1
        real(fpp), dimension(0:8) :: cf

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
        cf(:) = a0

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
            cf(DXX) = a0+d0
            cf(DYX) = a0+d0
            cf(DZX) = a0+d0
        case(1)
            b0(kB012) = k0
            b0(kB021) = 1./k0
            b0(kB120) = k0
            b1(kB012) = k0*d0
            b1(kB021) = -d0/k0
            b1(kB120) = k0*d0
            cf(DXY) = a0+d0
            cf(DYY) = a0+d0
            cf(DZY) = a0+d0
        case(2)
            b0(kB012) = 1./k0
            b0(kB021) = k0
            b0(kB120) = k0
            b1(kB012) = -d0/k0
            b1(kB021) = k0*d0
            b1(kB120) = k0*d0
            cf(DXZ) = a0+d0
            cf(DYZ) = a0+d0
            cf(DZZ) = a0+d0
        case default
            stop 1
        end select

        ! update convolution terms
        dt = dom%dt
        do r=0,8
            ! This loop relies on the fact that r=DXX or r=L120_DXX coincide for the right equations
            call cpml_compute_coefs(dom%cpml_integ, cf(r), dt, cf0, cf1, cf2)
            dom%R2_0(ee,r,i,j,k,bnum) = cf0*dom%R2_0(ee,r,i,j,k,bnum)+cf1*DUDVn(r)+cf2*DUDV(r)
        end do

        ! We add the terms in Dirac with the (only) convolution term
        LC(L120_DXX) = b0(kB120)*DUDVn(DXX) + b1(kB120)*dom%R2_0(ee,DXX,i,j,k,bnum)
        LC(L2_DYY  ) = b0(kB2  )*DUDVn(DYY) + b1(kB2  )*dom%R2_0(ee,DYY,i,j,k,bnum)
        LC(L1_DZZ  ) = b0(kB1  )*DUDVn(DZZ) + b1(kB1  )*dom%R2_0(ee,DZZ,i,j,k,bnum)
        LC(L120_DYX) = b0(kB120)*DUDVn(DYX) + b1(kB120)*dom%R2_0(ee,DYX,i,j,k,bnum)
        LC(L2_DYX  ) = b0(kB2  )*DUDVn(DYX) + b1(kB2  )*dom%R2_0(ee,DYX,i,j,k,bnum)
        LC(L120_DZX) = b0(kB120)*DUDVn(DZX) + b1(kB120)*dom%R2_0(ee,DZX,i,j,k,bnum)
        LC(L1_DZX  ) = b0(kB1  )*DUDVn(DZX) + b1(kB1  )*dom%R2_0(ee,DZX,i,j,k,bnum)
        LC(L021_DXY) = b0(kB021)*DUDVn(DXY) + b1(kB021)*dom%R2_0(ee,DXY,i,j,k,bnum)
        LC(L2_DXY  ) = b0(kB2  )*DUDVn(DXY) + b1(kB2  )*dom%R2_0(ee,DXY,i,j,k,bnum)
        LC(L2_DXX  ) = b0(kB2  )*DUDVn(DXX) + b1(kB2  )*dom%R2_0(ee,DXX,i,j,k,bnum)
        LC(L021_DYY) = b0(kB021)*DUDVn(DYY) + b1(kB021)*dom%R2_0(ee,DYY,i,j,k,bnum)
        LC(L0_DZZ  ) = b0(kB0  )*DUDVn(DZZ) + b1(kB0  )*dom%R2_0(ee,DZZ,i,j,k,bnum)
        LC(L021_DZY) = b0(kB021)*DUDVn(DZY) + b1(kB021)*dom%R2_0(ee,DZY,i,j,k,bnum)
        LC(L0_DZY  ) = b0(kB0  )*DUDVn(DZY) + b1(kB0  )*dom%R2_0(ee,DZY,i,j,k,bnum)
        LC(L012_DXZ) = b0(kB012)*DUDVn(DXZ) + b1(kB012)*dom%R2_0(ee,DXZ,i,j,k,bnum)
        LC(L1_DXZ  ) = b0(kB1  )*DUDVn(DXZ) + b1(kB1  )*dom%R2_0(ee,DXZ,i,j,k,bnum)
        LC(L012_DYZ) = b0(kB012)*DUDVn(DYZ) + b1(kB012)*dom%R2_0(ee,DYZ,i,j,k,bnum)
        LC(L0_DYZ  ) = b0(kB0  )*DUDVn(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum)
        LC(L1_DXX  ) = b0(kB1  )*DUDVn(DXX) + b1(kB1  )*dom%R2_0(ee,DXX,i,j,k,bnum)
        LC(L0_DYY  ) = b0(kB0  )*DUDVn(DYY) + b1(kB0  )*dom%R2_0(ee,DYY,i,j,k,bnum)
        LC(L012_DZZ) = b0(kB012)*DUDVn(DZZ) + b1(kB012)*dom%R2_0(ee,DZZ,i,j,k,bnum)
    end subroutine compute_convolution_terms_1d

    ! Compute convolution terms with atn in 2 directions
    subroutine compute_convolution_terms_2d(dom, i, j, k, bnum, ee, DUDVn, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:8) :: DUDV, DUDVn
        real(fpp), intent(out), dimension(0:20) :: LC
        !
        integer :: i1, r
        real(fpp), dimension(0:14) :: b0, b1, b2
        real(fpp), dimension(0:2) :: c0, c1, c2, rr0, rr1
        real(fpp), dimension(0:8) :: e0, e1
        real(fpp) :: dt, cf0, cf1, cf2, R0
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
            call cpml_compute_coefs(dom%cpml_integ, e0(r), dt, cf0, cf1, cf2)
            dom%R2_0(ee,r,i,j,k,bnum) = cf0*dom%R2_0(ee,r,i,j,k,bnum)+cf1*DUDVn(r)+cf2*DUDV(r)

            call cpml_compute_coefs(dom%cpml_integ, e1(r), dt, cf0, cf1, cf2)
            if (.not. isclose(e0(r),e1(r))) then
                dom%R2_1(r,i,j,k,i1) = cf0*dom%R2_1(r,i,j,k,i1)+cf1*DUDVn(r)+cf2*DUDV(r)
            else
                R0 = dom%R2_0(ee,r,i,j,k,bnum)
                dom%R2_1(r,i,j,k,i1) = cf0*dom%R2_1(r,i,j,k,i1)+cf1*(R0)+cf2*(R0)
            end if
        end do

        ! Convolve
        LC(L120_DXX) = b0(kB120)*DUDVn(DXX) + b1(kB120)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DXX,i,j,k,i1)
        LC(L120_DYX) = b0(kB120)*DUDVn(DYX) + b1(kB120)*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DYX,i,j,k,i1)
        LC(L120_DZX) = b0(kB120)*DUDVn(DZX) + b1(kB120)*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DZX,i,j,k,i1)
        LC(L021_DXY) = b0(kB021)*DUDVn(DXY) + b1(kB021)*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DXY,i,j,k,i1)
        LC(L021_DYY) = b0(kB021)*DUDVn(DYY) + b1(kB021)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DYY,i,j,k,i1)
        LC(L021_DZY) = b0(kB021)*DUDVn(DZY) + b1(kB021)*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DZY,i,j,k,i1)
        LC(L012_DXZ) = b0(kB012)*DUDVn(DXZ) + b1(kB012)*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DXZ,i,j,k,i1)
        LC(L012_DYZ) = b0(kB012)*DUDVn(DYZ) + b1(kB012)*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DYZ,i,j,k,i1)
        LC(L012_DZZ) = b0(kB012)*DUDVn(DZZ) + b1(kB012)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DZZ,i,j,k,i1)
        LC(L0_DYY  ) = b0(kB0YY)*DUDVn(DYY) + b1(kB0YY)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB0YY)*dom%R2_1(DYY,i,j,k,i1)
        LC(L0_DYZ  ) = b0(kB0YZ)*DUDVn(DYZ) + b1(kB0YZ)*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB0YZ)*dom%R2_1(DYZ,i,j,k,i1)
        LC(L0_DZY  ) = b0(kB0ZY)*DUDVn(DZY) + b1(kB0ZY)*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB0ZY)*dom%R2_1(DZY,i,j,k,i1)
        LC(L0_DZZ  ) = b0(kB0ZZ)*DUDVn(DZZ) + b1(kB0ZZ)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB0ZZ)*dom%R2_1(DZZ,i,j,k,i1)
        LC(L1_DXX  ) = b0(kB1XX)*DUDVn(DXX) + b1(kB1XX)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB1XX)*dom%R2_1(DXX,i,j,k,i1)
        LC(L1_DXZ  ) = b0(kB1XZ)*DUDVn(DXZ) + b1(kB1XZ)*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB1XZ)*dom%R2_1(DXZ,i,j,k,i1)
        LC(L1_DZX  ) = b0(kB1ZX)*DUDVn(DZX) + b1(kB1ZX)*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB1ZX)*dom%R2_1(DZX,i,j,k,i1)
        LC(L1_DZZ  ) = b0(kB1ZZ)*DUDVn(DZZ) + b1(kB1ZZ)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB1ZZ)*dom%R2_1(DZZ,i,j,k,i1)
        LC(L2_DXX  ) = b0(kB2XX)*DUDVn(DXX) + b1(kB2XX)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB2XX)*dom%R2_1(DXX,i,j,k,i1)
        LC(L2_DXY  ) = b0(kB2XY)*DUDVn(DXY) + b1(kB2XY)*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB2XY)*dom%R2_1(DXY,i,j,k,i1)
        LC(L2_DYX  ) = b0(kB2YX)*DUDVn(DYX) + b1(kB2YX)*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB2YX)*dom%R2_1(DYX,i,j,k,i1)
        LC(L2_DYY  ) = b0(kB2YY)*DUDVn(DYY) + b1(kB2YY)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB2YY)*dom%R2_1(DYY,i,j,k,i1)
    end subroutine compute_convolution_terms_2d

    ! Compute convolution terms with atn in 3 directions
    subroutine compute_convolution_terms_3d(dom, i, j, k, bnum, ee, DUDVn, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:8) :: DUDV, DUDVn
        real(fpp), intent(out), dimension(0:20) :: LC
        !
        integer   :: r, i1, i2
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        real(fpp) :: k2, d2, a2
        real(fpp) :: dt
        real(fpp) :: cf00, cf01, cf02, cf10, cf11, cf12, cf20, cf21, cf22
        real(fpp) :: R0, R1, R2
        real(fpp), dimension(0:14) :: b0, b1, b2, b3
        real(fpp), dimension(0:2) :: c0, c1, c2, c3
        real(fpp), dimension(0:2, 0:8) :: e
        integer, dimension(0:8) :: sel
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
        e(0:2,kB012) = (/ a0, a1, a2+d2 /)
        e(0:2,kB120) = (/ a0+d0, a1, a2 /)
        e(0:2,kB021) = (/ a0, a1+d1, a2 /)
        call get_coefs_Lijk_3d(k1,k2,k0,a1,a2,a0,d1,d2,d0,b0(kB120),b1(kB120),b2(kB120),b3(kB120),sel120)
        call get_coefs_Lijk_3d(k0,k2,k1,a0,a2,a1,d0,d2,d1,b0(kB021),b1(kB021),b2(kB021),b3(kB021),sel021)
        call get_coefs_Lijk_3d(k0,k1,k2,a0,a1,a2,d0,d1,d2,b0(kB012),b1(kB012),b2(kB012),b3(kB012),sel012)
        sel(DXX) = sel120; sel(DXY) = sel021; sel(DXZ) = sel012;
        sel(DYX) = sel120; sel(DYY) = sel021; sel(DYZ) = sel012;
        sel(DZX) = sel120; sel(DZY) = sel021; sel(DZZ) = sel012;

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
            call cpml_compute_coefs(dom%cpml_integ, e(0,r), dt, cf00, cf01, cf02)
            call cpml_compute_coefs(dom%cpml_integ, e(1,r), dt, cf10, cf11, cf12)
            call cpml_compute_coefs(dom%cpml_integ, e(2,r), dt, cf20, cf21, cf22)
            select case(sel(r))
            case (CMP_ABC)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*DUDVn(r) + cf12*DUDV(r)
                R2 = cf20*R2 + cf21*DUDVn(r) + cf22*DUDV(r)
            case (CMP_AAC)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*R0       + cf12*R0
                R2 = cf20*R2 + cf21*DUDVn(r) + cf22*DUDV(r)
            case (CMP_ABA)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*DUDVn(r) + cf12*DUDV(r)
                R2 = cf20*R2 + cf21*R0       + cf22*R0
            case (CMP_ABB)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*DUDVn(r) + cf12*DUDV(r)
                R2 = cf20*R2 + cf21*R1       + cf22*R1
            case (CMP_AAA)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*R0       + cf12*R0
                R2 = cf20*R2 + 2*(cf21*R1    + cf22*R1)
            case default
                stop 1
            end select
            dom%R2_0(ee,r,i,j,k,bnum) = R0
            dom%R2_1(r,i,j,k,i1) = R1
            dom%R2_2(r,i,j,k,i2) = R2
        end do

        ! Convolve
        LC(L120_DXX) = b0(kB120)*DUDVn(DXX) + b1(kB120)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DXX,i,j,k,i1) + b3(kB120)*dom%R2_2(DXX,i,j,k,i2)
        LC(L120_DYX) = b0(kB120)*DUDVn(DYX) + b1(kB120)*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DYX,i,j,k,i1) + b3(kB120)*dom%R2_2(DYX,i,j,k,i2)
        LC(L120_DZX) = b0(kB120)*DUDVn(DZX) + b1(kB120)*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DZX,i,j,k,i1) + b3(kB120)*dom%R2_2(DZX,i,j,k,i2)
        LC(L021_DXY) = b0(kB021)*DUDVn(DXY) + b1(kB021)*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DXY,i,j,k,i1) + b3(kB021)*dom%R2_2(DXY,i,j,k,i2)
        LC(L021_DYY) = b0(kB021)*DUDVn(DYY) + b1(kB021)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DYY,i,j,k,i1) + b3(kB021)*dom%R2_2(DYY,i,j,k,i2)
        LC(L021_DZY) = b0(kB021)*DUDVn(DZY) + b1(kB021)*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DZY,i,j,k,i1) + b3(kB021)*dom%R2_2(DZY,i,j,k,i2)
        LC(L012_DXZ) = b0(kB012)*DUDVn(DXZ) + b1(kB012)*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DXZ,i,j,k,i1) + b3(kB012)*dom%R2_2(DXZ,i,j,k,i2)
        LC(L012_DYZ) = b0(kB012)*DUDVn(DYZ) + b1(kB012)*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DYZ,i,j,k,i1) + b3(kB012)*dom%R2_2(DYZ,i,j,k,i2)
        LC(L012_DZZ) = b0(kB012)*DUDVn(DZZ) + b1(kB012)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DZZ,i,j,k,i1) + b3(kB012)*dom%R2_2(DZZ,i,j,k,i2)
        LC(L0_DYY  ) = b0(kB0YY)*DUDVn(DYY) + b1(kB0YY)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB0YY)*dom%R2_1(DYY,i,j,k,i1) + b3(kB0YY)*dom%R2_2(DYY,i,j,k,i2)
        LC(L0_DYZ  ) = b0(kB0YZ)*DUDVn(DYZ) + b1(kB0YZ)*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB0YZ)*dom%R2_1(DYZ,i,j,k,i1) + b3(kB0YZ)*dom%R2_2(DYZ,i,j,k,i2)
        LC(L0_DZY  ) = b0(kB0ZY)*DUDVn(DZY) + b1(kB0ZY)*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB0ZY)*dom%R2_1(DZY,i,j,k,i1) + b3(kB0ZY)*dom%R2_2(DZY,i,j,k,i2)
        LC(L0_DZZ  ) = b0(kB0ZZ)*DUDVn(DZZ) + b1(kB0ZZ)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB0ZZ)*dom%R2_1(DZZ,i,j,k,i1) + b3(kB0ZZ)*dom%R2_2(DZZ,i,j,k,i2)
        LC(L1_DXX  ) = b0(kB1XX)*DUDVn(DXX) + b1(kB1XX)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB1XX)*dom%R2_1(DXX,i,j,k,i1) + b3(kB1XX)*dom%R2_2(DXX,i,j,k,i2)
        LC(L1_DXZ  ) = b0(kB1XZ)*DUDVn(DXZ) + b1(kB1XZ)*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB1XZ)*dom%R2_1(DXZ,i,j,k,i1) + b3(kB1XZ)*dom%R2_2(DXZ,i,j,k,i2)
        LC(L1_DZX  ) = b0(kB1ZX)*DUDVn(DZX) + b1(kB1ZX)*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB1ZX)*dom%R2_1(DZX,i,j,k,i1) + b3(kB1ZX)*dom%R2_2(DZX,i,j,k,i2)
        LC(L1_DZZ  ) = b0(kB1ZZ)*DUDVn(DZZ) + b1(kB1ZZ)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB1ZZ)*dom%R2_1(DZZ,i,j,k,i1) + b3(kB1ZZ)*dom%R2_2(DZZ,i,j,k,i2)
        LC(L2_DXX  ) = b0(kB2XX)*DUDVn(DXX) + b1(kB2XX)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB2XX)*dom%R2_1(DXX,i,j,k,i1) + b3(kB2XX)*dom%R2_2(DXX,i,j,k,i2)
        LC(L2_DXY  ) = b0(kB2XY)*DUDVn(DXY) + b1(kB2XY)*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB2XY)*dom%R2_1(DXY,i,j,k,i1) + b3(kB2XY)*dom%R2_2(DXY,i,j,k,i2)
        LC(L2_DYX  ) = b0(kB2YX)*DUDVn(DYX) + b1(kB2YX)*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB2YX)*dom%R2_1(DYX,i,j,k,i1) + b3(kB2YX)*dom%R2_2(DYX,i,j,k,i2)
        LC(L2_DYY  ) = b0(kB2YY)*DUDVn(DYY) + b1(kB2YY)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB2YY)*dom%R2_1(DYY,i,j,k,i1) + b3(kB2YY)*dom%R2_2(DYY,i,j,k,i2)
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
