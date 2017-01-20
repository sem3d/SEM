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
    integer, parameter :: kB012=0, kB021=1, kB120=2, kB0=3, kB1=4, kB2=5
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
    subroutine compute_L_convolution_terms(dom, i, j, k, bnum, ee, Unew, Rx, Ry, Rz)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        real(fpp), intent(in), dimension(0:2) :: Unew
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(out) :: Rx, Ry, Rz
        !
        real(fpp) :: a3b
        real(fpp), dimension(0:2) :: Uold
        real(fpp) :: k0, d0, a0, dt, cf0,cf1,cf2
        integer :: n1, n2

        n1 = dom%I1(ee,bnum)
        n2 = dom%I2(ee,bnum)

        Uold = dom%Uold(ee,:,i,j,k,bnum)
        ! XXX valable pour ndir=1
        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        call cpml_compute_coefs(dom%cpml_integ, a0, dt, cf0, cf1, cf2)
        dom%R1_0(ee,0,i,j,k,bnum) = cf0*dom%R1_0(ee,0,i,j,k,bnum) + cf1*Unew(0) + cf2*Uold(0)
        dom%R1_0(ee,1,i,j,k,bnum) = cf0*dom%R1_0(ee,1,i,j,k,bnum) + cf1*Unew(1) + cf2*Uold(1)
        dom%R1_0(ee,2,i,j,k,bnum) = cf0*dom%R1_0(ee,2,i,j,k,bnum) + cf1*Unew(2) + cf2*Uold(2)

        if (n1==-1 .and. n2==-1) then
            a3b = k0*a0*a0*d0
            Rx = a3b*dom%R1_0(ee,0,i,j,k,bnum)
            Ry = a3b*dom%R1_0(ee,1,i,j,k,bnum)
            Rz = a3b*dom%R1_0(ee,2,i,j,k,bnum)
        else
            Rx = 0d0
            Ry = 0d0
            Rz = 0d0
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
        integer :: n1, n2
        n1 = dom%I1(ee,bnum)
        n2 = dom%I2(ee,bnum)

        DUDV = dom%DUDVold(ee,:,i,j,k,bnum) !0.5d0*DUDVnew+0.5d0*dom%DUDVold(ee,:,i,j,k,bnum)
        if (n2==-1) then
            if (n1==-1) then
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

    subroutine get_coefs_Li(dom, ee, bnum, i, j, k, b0, b1)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: ee, bnum, i, j, k
        real(fpp), dimension(0:5), intent (OUT) :: b0, b1
        !
        real(fpp) :: a0, k0, d0
        integer :: dim0

        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        if (k0<0) then
            write(*,*) "Erreur:", bnum, ee, i,j,k, k0
            stop 1
        endif

        dim0 = dom%D0(ee,bnum)
        select case(dim0)
        case(0)
            b0(kB0) = k0
            b0(kB1) = 1.
            b0(kB2) = 1.
            b1(kB0) = d0*k0
            b1(kB1) = 0.
            b1(kB2) = 0.
        case(1)
            b0(kB0) = 1.
            b0(kB1) = k0
            b0(kB2) = 1.
            b1(kB0) = 0.
            b1(kB1) = d0*k0
            b1(kB2) = 0.
        case(2)
            b0(kB0) = 1.
            b0(kB1) = 1.
            b0(kB2) = k0
            b1(kB0) = 0.
            b1(kB1) = 0.
            b1(kB2) = d0*k0
        case default
            stop 1
        end select
    end subroutine get_coefs_Li

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
        call get_coefs_Li(dom, ee, bnum, i, j, k, b0, b1)

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
        real(fpp), dimension(0:5) :: b0, b1, b2
        real(fpp), dimension(0:8) :: e0, e1
        real(fpp) :: dt, cf0, cf1, cf2, R0

        ! Initialize
        b0 = 1.
        b1 = 0.
        b2 = 0.
        e0 = dom%Alpha_0(ee,i,j,k,bnum)
        i1 = dom%I1(ee,bnum)
        e1 = dom%Alpha_1(i,j,k,i1)

        ! Lijk
        call get_coefs_Lijk_2d(dom, ee, bnum, i, j, k, kB120, kB021, kB012, b0, b1, b2, e0, e1)

        ! Li
        call get_coefs_Li(dom, ee, bnum, i, j, k, b0, b1)

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
        LC(L2_DYY  ) = b0(kB2  )*DUDVn(DYY) + b1(kB2  )*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DYY,i,j,k,i1)
        LC(L1_DZZ  ) = b0(kB1  )*DUDVn(DZZ) + b1(kB1  )*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DZZ,i,j,k,i1)
        LC(L120_DYX) = b0(kB120)*DUDVn(DYX) + b1(kB120)*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DYX,i,j,k,i1)
        LC(L2_DYX  ) = b0(kB2  )*DUDVn(DYX) + b1(kB2  )*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DYX,i,j,k,i1)
        LC(L120_DZX) = b0(kB120)*DUDVn(DZX) + b1(kB120)*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DZX,i,j,k,i1)
        LC(L1_DZX  ) = b0(kB1  )*DUDVn(DZX) + b1(kB1  )*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DZX,i,j,k,i1)
        LC(L021_DXY) = b0(kB021)*DUDVn(DXY) + b1(kB021)*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DXY,i,j,k,i1)
        LC(L2_DXY  ) = b0(kB2  )*DUDVn(DXY) + b1(kB2  )*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DXY,i,j,k,i1)
        LC(L2_DXX  ) = b0(kB2  )*DUDVn(DXX) + b1(kB2  )*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DXX,i,j,k,i1)
        LC(L021_DYY) = b0(kB021)*DUDVn(DYY) + b1(kB021)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DYY,i,j,k,i1)
        LC(L0_DZZ  ) = b0(kB0  )*DUDVn(DZZ) + b1(kB0  )*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DZZ,i,j,k,i1)
        LC(L021_DZY) = b0(kB021)*DUDVn(DZY) + b1(kB021)*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DZY,i,j,k,i1)
        LC(L0_DZY  ) = b0(kB0  )*DUDVn(DZY) + b1(kB0  )*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DZY,i,j,k,i1)
        LC(L012_DXZ) = b0(kB012)*DUDVn(DXZ) + b1(kB012)*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DXZ,i,j,k,i1)
        LC(L1_DXZ  ) = b0(kB1  )*DUDVn(DXZ) + b1(kB1  )*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DXZ,i,j,k,i1)
        LC(L012_DYZ) = b0(kB012)*DUDVn(DYZ) + b1(kB012)*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DYZ,i,j,k,i1)
        LC(L0_DYZ  ) = b0(kB0  )*DUDVn(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYZ,i,j,k,i1)
        LC(L1_DXX  ) = b0(kB1  )*DUDVn(DXX) + b1(kB1  )*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DXX,i,j,k,i1)
        LC(L0_DYY  ) = b0(kB0  )*DUDVn(DYY) + b1(kB0  )*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYY,i,j,k,i1)
        LC(L012_DZZ) = b0(kB012)*DUDVn(DZZ) + b1(kB012)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DZZ,i,j,k,i1)
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
        real(fpp), dimension(0:5) :: b0, b1, b2, b3
        real(fpp), dimension(0:2, 0:8) :: e
        integer, dimension(0:8) :: sel

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

        ! Lijk
        e(0:2,kB012) = (/ a0, a1, a2+d2 /)
        e(0:2,kB120) = (/ a1, a2, a0+d0 /)
        e(0:2,kB021) = (/ a0, a2, a1+d1 /)
        sel = CMP_ABC
        call get_coefs_Lijk_3d(k1,k2,k0,a1,a2,a0,d1,d2,d0,b0(kB120),b1(kB120),b2(kB120),b3(kB120),sel(kB120))
        call get_coefs_Lijk_3d(k0,k2,k1,a0,a2,a1,d0,d2,d1,b0(kB021),b1(kB021),b2(kB021),b3(kB021),sel(kB021))
        call get_coefs_Lijk_3d(k0,k1,k2,a0,a1,a2,d0,d1,d2,b0(kB012),b1(kB012),b2(kB012),b3(kB012),sel(kB012))

        ! Li
        call get_coefs_Li(dom, ee, bnum, i, j, k, b0, b1)

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
                R1 = cf10*R1 + cf11*R0         + cf12*R0
                R2 = cf20*R2 + cf21*DUDVn(r) + cf22*DUDV(r)
            case (CMP_ABA)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*DUDVn(r) + cf12*DUDV(r)
                R2 = cf20*R2 + cf21*R0         + cf22*R0
            case (CMP_ABB)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*DUDVn(r) + cf12*DUDV(r)
                R2 = cf20*R2 + cf21*R1         + cf22*R1
            case (CMP_AAA)
                R0 = cf00*R0 + cf01*DUDVn(r) + cf02*DUDV(r)
                R1 = cf10*R1 + cf11*R0         + cf12*R0
                R2 = cf20*R2 + cf21*R1         + cf22*R1
            case default
                stop 1
            end select
            dom%R2_0(ee,r,i,j,k,bnum) = R0
            dom%R2_1(r,i,j,k,i1) = R1
            dom%R2_2(r,i,j,k,i2) = R2
        end do

        ! Convolve
        LC(L120_DXX) = b0(kB120)*DUDV(DXX) + b1(kB120)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DXX,i,j,k,i1) + b3(kB120)*dom%R2_2(DXX,i,j,k,i2)
        LC(L2_DYY  ) = b0(kB2  )*DUDV(DYY) + b1(kB2  )*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DYY,i,j,k,i1) + b3(kB2  )*dom%R2_2(DYY,i,j,k,i2)
        LC(L1_DZZ  ) = b0(kB1  )*DUDV(DZZ) + b1(kB1  )*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DZZ,i,j,k,i1) + b3(kB1  )*dom%R2_2(DZZ,i,j,k,i2)
        LC(L120_DYX) = b0(kB120)*DUDV(DYX) + b1(kB120)*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DYX,i,j,k,i1) + b3(kB120)*dom%R2_2(DYX,i,j,k,i2)
        LC(L2_DYX  ) = b0(kB2  )*DUDV(DYX) + b1(kB2  )*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DYX,i,j,k,i1) + b3(kB2  )*dom%R2_2(DYX,i,j,k,i2)
        LC(L120_DZX) = b0(kB120)*DUDV(DZX) + b1(kB120)*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DZX,i,j,k,i1) + b3(kB120)*dom%R2_2(DZX,i,j,k,i2)
        LC(L1_DZX  ) = b0(kB1  )*DUDV(DZX) + b1(kB1  )*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DZX,i,j,k,i1) + b3(kB1  )*dom%R2_2(DZX,i,j,k,i2)
        LC(L021_DXY) = b0(kB021)*DUDV(DXY) + b1(kB021)*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DXY,i,j,k,i1) + b3(kB021)*dom%R2_2(DXY,i,j,k,i2)
        LC(L2_DXY  ) = b0(kB2  )*DUDV(DXY) + b1(kB2  )*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DXY,i,j,k,i1) + b3(kB2  )*dom%R2_2(DXY,i,j,k,i2)
        LC(L2_DXX  ) = b0(kB2  )*DUDV(DXX) + b1(kB2  )*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DXX,i,j,k,i1) + b3(kB2  )*dom%R2_2(DXX,i,j,k,i2)
        LC(L021_DYY) = b0(kB021)*DUDV(DYY) + b1(kB021)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DYY,i,j,k,i1) + b3(kB021)*dom%R2_2(DYY,i,j,k,i2)
        LC(L0_DZZ  ) = b0(kB0  )*DUDV(DZZ) + b1(kB0  )*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DZZ,i,j,k,i1) + b3(kB0  )*dom%R2_2(DZZ,i,j,k,i2)
        LC(L021_DZY) = b0(kB021)*DUDV(DZY) + b1(kB021)*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DZY,i,j,k,i1) + b3(kB021)*dom%R2_2(DZY,i,j,k,i2)
        LC(L0_DZY  ) = b0(kB0  )*DUDV(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYZ,i,j,k,i1) + b3(kB0  )*dom%R2_2(DYZ,i,j,k,i2)
        LC(L012_DXZ) = b0(kB012)*DUDV(DXZ) + b1(kB012)*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DXZ,i,j,k,i1) + b3(kB012)*dom%R2_2(DXZ,i,j,k,i2)
        LC(L1_DXZ  ) = b0(kB1  )*DUDV(DXZ) + b1(kB1  )*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DXZ,i,j,k,i1) + b3(kB1  )*dom%R2_2(DXZ,i,j,k,i2)
        LC(L012_DYZ) = b0(kB012)*DUDV(DYZ) + b1(kB012)*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DYZ,i,j,k,i1) + b3(kB012)*dom%R2_2(DYZ,i,j,k,i2)
        LC(L0_DYZ  ) = b0(kB0  )*DUDV(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYZ,i,j,k,i1) + b3(kB0  )*dom%R2_2(DYZ,i,j,k,i2)
        LC(L1_DXX  ) = b0(kB1  )*DUDV(DXX) + b1(kB1  )*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DXX,i,j,k,i1) + b3(kB1  )*dom%R2_2(DXX,i,j,k,i2)
        LC(L0_DYY  ) = b0(kB0  )*DUDV(DYY) + b1(kB0  )*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYY,i,j,k,i1) + b3(kB0  )*dom%R2_2(DYY,i,j,k,i2)
        LC(L012_DZZ) = b0(kB012)*DUDV(DZZ) + b1(kB012)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DZZ,i,j,k,i1) + b3(kB012)*dom%R2_2(DZZ,i,j,k,i2)
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
