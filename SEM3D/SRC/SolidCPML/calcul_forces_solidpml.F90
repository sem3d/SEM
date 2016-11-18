!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"

module m_calcul_forces_solidpml
    use constants
    implicit none

    integer, parameter :: L120_DXX= 0, L2_DYY= 1, L1_DZZ= 2
    integer, parameter :: L120_DXY= 3, L2_DYX= 4
    integer, parameter :: L120_DXZ= 5, L1_DZX= 6
    integer, parameter :: L021_DYX= 7, L2_DXY= 8
    integer, parameter :: L2_DXX= 9, L021_DYY=10, L0_DZZ=11
    integer, parameter :: L021_DYZ=12, L0_DZY=13
    integer, parameter :: L012_DZX=14, L1_DXZ=15
    integer, parameter :: L012_DZY=16, L0_DYZ=17
    integer, parameter :: L1_DXX=18, L0_DYY=19, L012_DZZ=20
    integer, parameter :: dXX=0, dXY=1, dXZ=2
    integer, parameter :: dYX=3, dYY=4, dYZ=5
    integer, parameter :: dZX=6, dZY=7, dZZ=8
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
        real(fpp), dimension(0:2) :: U
        real(fpp) :: k0, d0, a0, dt, cf0,cf1
        integer :: n1, n2

        n1 = dom%I1(ee,bnum)
        n2 = dom%I2(ee,bnum)

        U = 0.5d0*Unew + 0.5d0*dom%Uold(ee,:,i,j,k,bnum)
        ! XXX valable pour ndir=1
        dt = dom%dt
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        if (n1==-1 .and. n2==-1) then
            ! Update convolution term (implicit midpoint)
            cf0 = 1d0-0.5d0*a0*dt
            cf1 = 1d0/(1d0+0.5d0*a0*dt)
            dom%R1_0(ee,0,i,j,k,bnum) = (cf0*dom%R1_0(ee,0,i,j,k,bnum) + dt*U(0))*cf1
            dom%R1_0(ee,1,i,j,k,bnum) = (cf0*dom%R1_0(ee,1,i,j,k,bnum) + dt*U(1))*cf1
            dom%R1_0(ee,2,i,j,k,bnum) = (cf0*dom%R1_0(ee,2,i,j,k,bnum) + dt*U(2))*cf1
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

        DUDV = 0.5d0*DUDVnew+0.5d0*dom%DUDVold(ee,:,i,j,k,bnum)
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

    subroutine compute_convolution_terms_1d(dom, i, j, k, bnum, ee, DUDVn, DUDV, LC)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        real(fpp), intent(in), dimension(0:8) :: DUDV, DUDVn
        real(fpp), intent(out), dimension(0:20) :: LC
        !
        integer   :: dim0, r
        real(fpp) :: k0, d0, a0, dt, cf0, cf1
        real(fpp), dimension(0:5) :: b0, b1
        real(fpp), dimension(0:8) :: cf

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
            b0(kB0) = k0
            b0(kB1) = 1.
            b0(kB2) = 1.
            b1(kB012) = k0*d0  !b1
            b1(kB021) = k0*d0  !b1
            b1(kB120) = -d0/k0 !b3
            b1(kB0) = d0*k0
            b1(kB1) = 0.
            b1(kB2) = 0.
            cf(DXX) = a0+d0
            cf(DXY) = a0+d0
            cf(DXZ) = a0+d0
        case(1)
            b0(kB012) = k0
            b0(kB021) = 1./k0
            b0(kB120) = k0
            b0(kB0) = 1.
            b0(kB1) = k0
            b0(kB2) = 1.
            b1(kB012) = k0*d0
            b1(kB021) = -d0/k0
            b1(kB120) = k0*d0
            b1(kB0) = 0.
            b1(kB1) = d0*k0
            b1(kB2) = 0.
            cf(DYX) = a0+d0
            cf(DYY) = a0+d0
            cf(DYZ) = a0+d0
        case(2)
            b0(kB012) = 1./k0
            b0(kB021) = k0
            b0(kB120) = k0
            b0(kB0) = 1.
            b0(kB1) = 1.
            b0(kB2) = k0
            b1(kB012) = -d0/k0
            b1(kB021) = k0*d0
            b1(kB120) = k0*d0
            b1(kB0) = 0.
            b1(kB1) = 0.
            b1(kB2) = d0*k0
            cf(DZX) = a0+d0
            cf(DZY) = a0+d0
            cf(DZZ) = a0+d0
        case default
            stop 1
        end select

        ! update convolution terms
        do r=0,8
            ! implicit midpoint
            cf0 = 1d0-0.5d0*cf(r)*dt
            cf1 = 1d0/(1d0+0.5d0*cf(r)*dt)
            dom%R2_0(ee,r,i,j,k,bnum) = (cf0*dom%R2_0(ee,r,i,j,k,bnum)+dt*DUDV(r))*cf1
        end do
        ! We add the terms in Dirac with the (only) convolution term
        LC(L120_DXX) = b0(kB120)*DUDVn(DXX) + b1(kB120)*dom%R2_0(ee,DXX,i,j,k,bnum)
        LC(L2_DYY  ) = b0(kB2  )*DUDVn(DYY) + b1(kB2  )*dom%R2_0(ee,DYY,i,j,k,bnum)
        LC(L1_DZZ  ) = b0(kB1  )*DUDVn(DZZ) + b1(kB1  )*dom%R2_0(ee,DZZ,i,j,k,bnum)
        LC(L120_DXY) = b0(kB120)*DUDVn(DXY) + b1(kB120)*dom%R2_0(ee,DXY,i,j,k,bnum)
        LC(L2_DYX  ) = b0(kB2  )*DUDVn(DYX) + b1(kB2  )*dom%R2_0(ee,DYX,i,j,k,bnum)
        LC(L120_DXZ) = b0(kB120)*DUDVn(DXZ) + b1(kB120)*dom%R2_0(ee,DXZ,i,j,k,bnum)
        LC(L1_DZX  ) = b0(kB1  )*DUDVn(DZX) + b1(kB1  )*dom%R2_0(ee,DZX,i,j,k,bnum)
        LC(L021_DYX) = b0(kB021)*DUDVn(DYX) + b1(kB021)*dom%R2_0(ee,DYX,i,j,k,bnum)
        LC(L2_DXY  ) = b0(kB2  )*DUDVn(DXY) + b1(kB2  )*dom%R2_0(ee,DXY,i,j,k,bnum)
        LC(L2_DXX  ) = b0(kB2  )*DUDVn(DXX) + b1(kB2  )*dom%R2_0(ee,DXX,i,j,k,bnum)
        LC(L021_DYY) = b0(kB021)*DUDVn(DYY) + b1(kB021)*dom%R2_0(ee,DYY,i,j,k,bnum)
        LC(L0_DZZ  ) = b0(kB0  )*DUDVn(DZZ) + b1(kB0  )*dom%R2_0(ee,DZZ,i,j,k,bnum)
        LC(L021_DYZ) = b0(kB021)*DUDVn(DYZ) + b1(kB021)*dom%R2_0(ee,DYZ,i,j,k,bnum)
        LC(L0_DZY  ) = b0(kB0  )*DUDVn(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum)
        LC(L012_DZX) = b0(kB012)*DUDVn(DZX) + b1(kB012)*dom%R2_0(ee,DZX,i,j,k,bnum)
        LC(L1_DXZ  ) = b0(kB1  )*DUDVn(DXZ) + b1(kB1  )*dom%R2_0(ee,DXZ,i,j,k,bnum)
        LC(L012_DZY) = b0(kB012)*DUDVn(DZY) + b1(kB012)*dom%R2_0(ee,DZY,i,j,k,bnum)
        LC(L0_DYZ  ) = b0(kB0  )*DUDVn(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum)
        LC(L1_DXX  ) = b0(kB1  )*DUDVn(DXX) + b1(kB1  )*dom%R2_0(ee,DXX,i,j,k,bnum)
        LC(L0_DYY  ) = b0(kB0  )*DUDVn(DYY) + b1(kB0  )*dom%R2_0(ee,DYY,i,j,k,bnum)
        LC(L012_DZZ) = b0(kB012)*DUDVn(DZZ) + b1(kB012)*dom%R2_0(ee,DZZ,i,j,k,bnum)

!        if (i==1.and.j==0.and.k==0.and.bnum==1.and.ee==0) then
!            write(*,*) b0(kB120), DUDVn(DXX), b1(kB120),dom%R2_0(ee,DXX,i,j,k,bnum),LC(L120_DXX)
!        end if
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
        ! First we add the terms in Dirac
        LC(L120_DXX) = DUDV(DXX)
        LC(L2_DYY  ) = DUDV(DYY)
        LC(L1_DZZ  ) = DUDV(DZZ)
        LC(L120_DXY) = DUDV(DXY)
        LC(L2_DYX  ) = DUDV(DYX)
        LC(L120_DXZ) = DUDV(DXZ)
        LC(L1_DZX  ) = DUDV(DZX)
        LC(L021_DYX) = DUDV(DYX)
        LC(L2_DXY  ) = DUDV(DXY)
        LC(L2_DXX  ) = DUDV(DXX)
        LC(L021_DYY) = DUDV(DYY)
        LC(L0_DZZ  ) = DUDV(DZZ)
        LC(L021_DYZ) = DUDV(DYZ)
        LC(L0_DZY  ) = DUDV(DYZ)
        LC(L012_DZX) = DUDV(DZX)
        LC(L1_DXZ  ) = DUDV(DXZ)
        LC(L012_DZY) = DUDV(DZY)
        LC(L0_DYZ  ) = DUDV(DYZ)
        LC(L1_DXX  ) = DUDV(DXX)
        LC(L0_DYY  ) = DUDV(DYY)
        LC(L012_DZZ) = DUDV(DZZ)
        stop 1
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
!        integer   :: dim0, r, i1, i2
!        real(fpp) :: k0, d0, a0, dt
!        real(fpp), dimension(0:5) :: b0, b1, b2, b3
!        real(fpp), dimension(0:8) :: cf
        !
        stop 1
        !
!        LC(L120_DXX) = b0(kB120)*DUDV(DXX) + b1(kB120)*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB120)*dom%R2_1(DXX,i,j,k,i1) + b3(kB120)*dom%R2_2(DXX,i,j,k,i2)
!        LC(L2_DYY  ) = b0(kB2  )*DUDV(DYY) + b1(kB2  )*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DYY,i,j,k,i1) + b3(kB2  )*dom%R2_2(DYY,i,j,k,i2)
!        LC(L1_DZZ  ) = b0(kB1  )*DUDV(DZZ) + b1(kB1  )*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DZZ,i,j,k,i1) + b3(kB1  )*dom%R2_2(DZZ,i,j,k,i2)
!        LC(L120_DXY) = b0(kB120)*DUDV(DXY) + b1(kB120)*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB120)*dom%R2_1(DXY,i,j,k,i1) + b3(kB120)*dom%R2_2(DXY,i,j,k,i2)
!        LC(L2_DYX  ) = b0(kB2  )*DUDV(DYX) + b1(kB2  )*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DYX,i,j,k,i1) + b3(kB2  )*dom%R2_2(DYX,i,j,k,i2)
!        LC(L120_DXZ) = b0(kB120)*DUDV(DXZ) + b1(kB120)*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB120)*dom%R2_1(DXZ,i,j,k,i1) + b3(kB120)*dom%R2_2(DXZ,i,j,k,i2)
!        LC(L1_DZX  ) = b0(kB1  )*DUDV(DZX) + b1(kB1  )*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DZX,i,j,k,i1) + b3(kB1  )*dom%R2_2(DZX,i,j,k,i2)
!        LC(L021_DYX) = b0(kB021)*DUDV(DYX) + b1(kB021)*dom%R2_0(ee,DYX,i,j,k,bnum) + b2(kB021)*dom%R2_1(DYX,i,j,k,i1) + b3(kB021)*dom%R2_2(DYX,i,j,k,i2)
!        LC(L2_DXY  ) = b0(kB2  )*DUDV(DXY) + b1(kB2  )*dom%R2_0(ee,DXY,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DXY,i,j,k,i1) + b3(kB2  )*dom%R2_2(DXY,i,j,k,i2)
!        LC(L2_DXX  ) = b0(kB2  )*DUDV(DXX) + b1(kB2  )*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB2  )*dom%R2_1(DXX,i,j,k,i1) + b3(kB2  )*dom%R2_2(DXX,i,j,k,i2)
!        LC(L021_DYY) = b0(kB021)*DUDV(DYY) + b1(kB021)*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB021)*dom%R2_1(DYY,i,j,k,i1) + b3(kB021)*dom%R2_2(DYY,i,j,k,i2)
!        LC(L0_DZZ  ) = b0(kB0  )*DUDV(DZZ) + b1(kB0  )*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DZZ,i,j,k,i1) + b3(kB0  )*dom%R2_2(DZZ,i,j,k,i2)
!        LC(L021_DYZ) = b0(kB021)*DUDV(DYZ) + b1(kB021)*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB021)*dom%R2_1(DYZ,i,j,k,i1) + b3(kB021)*dom%R2_2(DYZ,i,j,k,i2)
!        LC(L0_DZY  ) = b0(kB0  )*DUDV(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYZ,i,j,k,i1) + b3(kB0  )*dom%R2_2(DYZ,i,j,k,i2)
!        LC(L012_DZX) = b0(kB012)*DUDV(DZX) + b1(kB012)*dom%R2_0(ee,DZX,i,j,k,bnum) + b2(kB012)*dom%R2_1(DZX,i,j,k,i1) + b3(kB012)*dom%R2_2(DZX,i,j,k,i2)
!        LC(L1_DXZ  ) = b0(kB1  )*DUDV(DXZ) + b1(kB1  )*dom%R2_0(ee,DXZ,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DXZ,i,j,k,i1) + b3(kB1  )*dom%R2_2(DXZ,i,j,k,i2)
!        LC(L012_DZY) = b0(kB012)*DUDV(DZY) + b1(kB012)*dom%R2_0(ee,DZY,i,j,k,bnum) + b2(kB012)*dom%R2_1(DZY,i,j,k,i1) + b3(kB012)*dom%R2_2(DZY,i,j,k,i2)
!        LC(L0_DYZ  ) = b0(kB0  )*DUDV(DYZ) + b1(kB0  )*dom%R2_0(ee,DYZ,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYZ,i,j,k,i1) + b3(kB0  )*dom%R2_2(DYZ,i,j,k,i2)
!        LC(L1_DXX  ) = b0(kB1  )*DUDV(DXX) + b1(kB1  )*dom%R2_0(ee,DXX,i,j,k,bnum) + b2(kB1  )*dom%R2_1(DXX,i,j,k,i1) + b3(kB1  )*dom%R2_2(DXX,i,j,k,i2)
!        LC(L0_DYY  ) = b0(kB0  )*DUDV(DYY) + b1(kB0  )*dom%R2_0(ee,DYY,i,j,k,bnum) + b2(kB0  )*dom%R2_1(DYY,i,j,k,i1) + b3(kB0  )*dom%R2_2(DYY,i,j,k,i2)
!        LC(L012_DZZ) = b0(kB012)*DUDV(DZZ) + b1(kB012)*dom%R2_0(ee,DZZ,i,j,k,bnum) + b2(kB012)*dom%R2_1(DZZ,i,j,k,i1) + b3(kB012)*dom%R2_2(DZZ,i,j,k,i2)
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
