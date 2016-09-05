!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

!! Convolutional Perfectly Matched Layers implemented according 2 references:
!! Ref1. Improved forward wave propagation and adjoint-based sensitivity kernel calculations using a numerically stable finite-element PML
!!       Zhinan Xie, Dimitri Komatitsch, Roland Martin, Rene Matzen
!!       Geophysical Journal International, 2014, 198, 1714-1747
!! Ref2. An efficient finite element time-domain formulation for the elastic second-order wave equation: a non-split complex frequency shifted convolutional PML
!!       Rene Matzen
!!       International Journal For Numerical Methods In Engineering, 2011, 88, 951-973

!! CPML can be seen as an artificial anisotropic material according to (10) + (14) + (17c) + (19*) from Ref2.
!! Convoluting div(sigma) amounts to convoluting Cijkl : this can be seen in (17c) and in the first term at the right hand side of (24c) from Ref2.

!! This file has been initialized from the version of the Solid domain with ANISO and without ATTENUATION (which is meaningless for CPML).
!! This file differs from the version of the Solid domain by :
!! 1. Cijkl are convoluted at each time iteration.
!! 2. Cijkl has 15 terms (upper part).
!!
!!      sigij               C0ijkl                 epskl                       Conv = FT^-1(si*sj/sk)
!!
!!    | sig11 |   | l+2m l    l                   |   | eps11 |     | L231 L3   L2                                 |
!!    | sig22 |   | l    l+2m l                   |   | eps22 |     | L3   L132 L1                                 |
!!    | sig33 | = | l    l    l+2m                |   | eps33 |     | L2   L1   L123                               |
!!    | sig23 |   |                 m  m          |   | eps23 |     |                L132 L1                       |
!!    | sig32 |   |                 m  m          | * | eps32 | and |                L1   L123                     |
!!    | sig13 |   |                      m m      |   | eps13 |     |                          L231 L2             |
!!    | sig31 |   |                      m m      |   | eps31 |     |                          L2   L123           |
!!    | sig12 |   |                           m m |   | eps12 |     |                                    L231 L3   |
!!    | sig21 |   |                           m m |   | eps21 |     |                                    L3   L132 |
!!
!!    Li and Lijk (from Conv) are defined in (A.4) from Ref1.
!!
!!    We have Cijkl = C0ijkl .* Conv (.* is an element-by-element multiplication) according to (10) from Ref1.

#include "dom_solidcpml_macro.F90"

#ifndef PROCNAME
#error "Must define PROCNAME before inclusion"
#endif

#undef ngllvar
#ifndef NGLLVAL
#define ngllvar ngll
#else
#define ngllvar ngll0
#endif

#define USE_LOCAL_MEM

subroutine PROCNAME(dom,ngllvar,bnum,Fox,Foy,Foz,Depla)
    use champs_solidpml
    implicit none

    type(domain_solidpml), intent (INOUT) :: dom
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    integer, intent(in) :: bnum
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: Fox,Foz,Foy
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(in) :: Depla
    !
    integer :: i,j,k,l,e,ee
    real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real(fpp) :: xi, Li, xj, xk, Lijk, xoverl, kappa(0:2) ! solidcpml_Li, solidcpml_Lijk
    real(fpp) :: sxx,sxy,sxz,syy,syz,szz,syx,szx,szy
    real(fpp) :: EXY,EXZ,EYZ
    real(fpp) :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
    real(fpp), parameter :: zero = 0.
    real(fpp) :: t4
#ifdef USE_LOCAL_MEM
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t1,t5,t8
    ! Les indices sont reordonnes, probablement pour la localite memoire
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t2,t6,t9
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t3,t7,t10
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: dUx_dxi, dUx_deta, dUx_dzeta
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: dUy_dxi, dUy_deta, dUy_dzeta
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: dUz_dxi, dUz_deta, dUz_dzeta
    real(fpp) :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
    real(fpp) :: F1,F2,F3
#endif
    integer :: mi
    real(fpp) :: lbd, mu
    real(fpp) :: C1111,C1122,C1133,C2222,C2233,C3333,C2323,C2332,C1313,C1331,C1212,C1221,C3232,C3131,C2121

    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(e,ee)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                local_deriv_ijke(Depla,0,dUx_dxi(ee,i,j,k),dUx_deta(ee,i,j,k),dUx_dzeta(ee,i,j,k))
                END_SUBELEM_LOOP()
            end do
        end do
    end do
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(e,ee)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                local_deriv_ijke(Depla,1,dUy_dxi(ee,i,j,k),dUy_deta(ee,i,j,k),dUy_dzeta(ee,i,j,k))
                END_SUBELEM_LOOP()
            end do
        end do
    end do
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(e,ee)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                local_deriv_ijke(Depla,2,dUz_dxi(ee,i,j,k),dUz_deta(ee,i,j,k),dUz_dzeta(ee,i,j,k))
                END_SUBELEM_LOOP()
            end do
        end do
    end do

    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(e,ee)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
#define JXX dom%InvGrad_(0,0,i,j,k,bnum,ee)
#define JXY dom%InvGrad_(0,1,i,j,k,bnum,ee)
#define JXZ dom%InvGrad_(0,2,i,j,k,bnum,ee)
#define JYX dom%InvGrad_(1,0,i,j,k,bnum,ee)
#define JYY dom%InvGrad_(1,1,i,j,k,bnum,ee)
#define JYZ dom%InvGrad_(1,2,i,j,k,bnum,ee)
#define JZX dom%InvGrad_(2,0,i,j,k,bnum,ee)
#define JZY dom%InvGrad_(2,1,i,j,k,bnum,ee)
#define JZZ dom%InvGrad_(2,2,i,j,k,bnum,ee)

                dxx=dUx_dxi(ee,i,j,k)*JXX + dUx_deta(ee,i,j,k)*JXY + dUx_dzeta(ee,i,j,k)*JXZ
                dxy=dUx_dxi(ee,i,j,k)*JYX + dUx_deta(ee,i,j,k)*JYY + dUx_dzeta(ee,i,j,k)*JYZ
                dxz=dUx_dxi(ee,i,j,k)*JZX + dUx_deta(ee,i,j,k)*JZY + dUx_dzeta(ee,i,j,k)*JZZ

                dyx=dUy_dxi(ee,i,j,k)*JXX + dUy_deta(ee,i,j,k)*JXY + dUy_dzeta(ee,i,j,k)*JXZ
                dyy=dUy_dxi(ee,i,j,k)*JYX + dUy_deta(ee,i,j,k)*JYY + dUy_dzeta(ee,i,j,k)*JYZ
                dyz=dUy_dxi(ee,i,j,k)*JZX + dUy_deta(ee,i,j,k)*JZY + dUy_dzeta(ee,i,j,k)*JZZ

                dzx=dUz_dxi(ee,i,j,k)*JXX + dUz_deta(ee,i,j,k)*JXY + dUz_dzeta(ee,i,j,k)*JXZ
                dzy=dUz_dxi(ee,i,j,k)*JYX + dUz_deta(ee,i,j,k)*JYY + dUz_dzeta(ee,i,j,k)*JYZ
                dzz=dUz_dxi(ee,i,j,k)*JZX + dUz_deta(ee,i,j,k)*JZY + dUz_dzeta(ee,i,j,k)*JZZ

                mi = dom%specel(bnum*VCHUNK+ee)%mat_index
                lbd = dom%sSubDomain(mi)%DLambda
                mu = dom%sSubDomain(mi)%DMu
                solidcpml_Lijk(1,2,0,i,j,k,bnum,ee,mi); C1111 = (lbd+2*mu)*Lijk
                solidcpml_Li      (2,i,j,k,bnum,ee,mi); C1122 =  lbd      *Li
                solidcpml_Li      (1,i,j,k,bnum,ee,mi); C1133 =  lbd      *Li
                solidcpml_Lijk(0,2,1,i,j,k,bnum,ee,mi); C2222 = (lbd+2*mu)*Lijk
                solidcpml_Li      (0,i,j,k,bnum,ee,mi); C2233 =  lbd      *Li
                solidcpml_Lijk(0,1,2,i,j,k,bnum,ee,mi); C3333 = (lbd+2*mu)*Lijk
                solidcpml_Lijk(0,2,1,i,j,k,bnum,ee,mi); C2323 =        mu *Lijk
                solidcpml_Li      (0,i,j,k,bnum,ee,mi); C2332 =        mu *Li
                solidcpml_Lijk(0,1,2,i,j,k,bnum,ee,mi); C3232 =        mu *Lijk
                solidcpml_Lijk(1,2,0,i,j,k,bnum,ee,mi); C1313 =        mu *Lijk
                solidcpml_Li      (1,i,j,k,bnum,ee,mi); C1331 =        mu *Li
                solidcpml_Lijk(0,1,2,i,j,k,bnum,ee,mi); C3131 =        mu *Lijk
                solidcpml_Lijk(1,2,0,i,j,k,bnum,ee,mi); C1212 =        mu *Lijk
                solidcpml_Li      (2,i,j,k,bnum,ee,mi); C1221 =        mu *Li
                solidcpml_Lijk(0,2,1,i,j,k,bnum,ee,mi); C2121 =        mu *Lijk

                EXY=0.5*(DXY+DYX)
                EXZ=0.5*(DXZ+DZX)
                EYZ=0.5*(DYZ+DZY)

                sxx = DXX*C1111 + DYY*C1122 + DZZ*C1133
                syy = DXX*C1122 + DYY*C2222 + DZZ*C2233
                szz = DXX*C1133 + DYY*C2233 + DZZ*C3333
                sxy = EXY*C1212 + EXY*C1221
                syx = EXY*C1221 + EXY*C2121
                sxz = EXZ*C1313 + EXZ*C1331
                szx = EXZ*C1331 + EXZ*C3131
                syz = EYZ*C2323 + EYZ*C2332
                szy = EYZ*C2332 + EYZ*C3232

                !
                !=====================
                !       FX
                !=====================
                !
                xt1 = sxx*JXX+sxy*JYX+sxz*JZX
                xt2 = sxx*JXY+sxy*JYY+sxz*JZY
                xt3 = sxx*JXZ+sxy*JYZ+sxz*JZZ
                !
                !=====================
                !       FY
                !=====================
                !
                xt5 = syx*JXX+syy*JYX+syz*JZX
                xt6 = syx*JXY+syy*JYY+syz*JZY
                xt7 = syx*JXZ+syy*JYZ+syz*JZZ
                !
                !=====================
                !       FZ
                !=====================
                !
                xt8 = szx*JXX+szy*JYX+szz*JZX
                xt9 = szx*JXY+szy*JYY+szz*JZY
                xt10= szx*JXZ+szy*JYZ+szz*JZZ

                !
                !- Multiplication par le Jacobien et le poids d'integration
                !
                t4 = dom%Jacob_(i,j,k,bnum,ee) * dom%GLLw(i)
                xt1  =  xt1 * t4
                xt5  =  xt5 * t4
                xt8  =  xt8 * t4

                t4 = dom%Jacob_(i,j,k,bnum,ee) * dom%GLLw(j)
                xt2  =  xt2 * t4
                xt6  =  xt6 * t4
                xt9  =  xt9 * t4

                t4 = dom%Jacob_(i,j,k,bnum,ee) * dom%GLLw(k)
                xt3  =  xt3 * t4
                xt7  =  xt7 * t4
                xt10 = xt10 * t4
#ifdef USE_LOCAL_MEM
                t1(ee,i,j,k) = xt1
                t5(ee,i,j,k) = xt5
                t8(ee,i,j,k) = xt8

                t2(ee,j,i,k) = xt2
                t6(ee,j,i,k) = xt6
                t9(ee,j,i,k) = xt9

                t3 (ee,k,i,j) = xt3
                t7 (ee,k,i,j) = xt7
                t10(ee,k,i,j) = xt10
#else
                DO L = 0, ngll-1
                    Fox(ee,L,j,k) = Fox(ee,L,j,k) + dom%hprime(L,I)*xt1*dom%GLLw(j) * dom%GLLw(k)
                    Foy(ee,L,j,k) = Foy(ee,L,j,k) + dom%hprime(L,I)*xt5*dom%GLLw(j) * dom%GLLw(k)
                    Foz(ee,L,j,k) = Foz(ee,L,j,k) + dom%hprime(L,I)*xt8*dom%GLLw(j) * dom%GLLw(k)
                END DO
                DO L = 0, ngll-1
                    Fox(ee,i,L,k) = Fox(ee,i,L,k) + dom%hprime(L,J)*xt2*dom%GLLw(i) * dom%GLLw(k)
                    Foy(ee,i,L,k) = Foy(ee,i,L,k) + dom%hprime(L,J)*xt6*dom%GLLw(i) * dom%GLLw(k)
                    Foz(ee,i,L,k) = Foz(ee,i,L,k) + dom%hprime(L,J)*xt9*dom%GLLw(i) * dom%GLLw(k)
                END DO
                DO L = 0, ngll-1
                    Fox(ee,i,j,L) = Fox(ee,i,j,L) + dom%hprime(L,K)*xt3*dom%GLLw(i) * dom%GLLw(j)
                    Foy(ee,i,j,L) = Foy(ee,i,j,L) + dom%hprime(L,K)*xt7*dom%GLLw(i) * dom%GLLw(j)
                    Foz(ee,i,j,L) = Foz(ee,i,j,L) + dom%hprime(L,K)*xt10*dom%GLLw(i) * dom%GLLw(j)
                END DO
#endif
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo
    !
    !- Multiplication par la matrice de derivation puis par les poids
    !
#ifdef USE_LOCAL_MEM
    !=-=-=-=-=-=-=-=-=-=-
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
                !dir$ simd
                !dir$ vector aligned
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                !=-=-=-=-=-=-=-=-=-=-
                !
                t11 = dom%GLLw(j) * dom%GLLw(k)
                t12 = dom%GLLw(i) * dom%GLLw(k)
                t13 = dom%GLLw(i) * dom%GLLw(j)
                !
                t41 = zero
                t42 = zero
                t43 = zero
                t51 = zero
                t52 = zero
                t53 = zero
                t61 = zero
                t62 = zero
                t63 = zero
                !
                do l = 0,ngll-1
                    t41 = t41 + dom%htprime(l,i) * t1(ee,l,j,k)
                    t42 = t42 + dom%htprime(l,i) * t5(ee,l,j,k)
                    t43 = t43 + dom%htprime(l,i) * t8(ee,l,j,k)
                enddo

                do l = 0,ngll-1
                    t51 = t51 + dom%htprime(l,j) * t2(ee,l,i,k)
                    t52 = t52 + dom%htprime(l,j) * t6(ee,l,i,k)
                    t53 = t53 + dom%htprime(l,j) * t9(ee,l,i,k)
                enddo
                ! FX
                F1 = t41*t11 + t51*t12
                ! FY
                F2 = t42*t11 + t52*t12
                ! FZ
                F3 = t43*t11 + t53*t12
                !
                !
                do l = 0,ngll-1
                    t61 = t61 + dom%htprime(l,k) * t3 (ee,l,i,j)
                    t62 = t62 + dom%htprime(l,k) * t7 (ee,l,i,j)
                    t63 = t63 + dom%htprime(l,k) * t10(ee,l,i,j)
                enddo

                ! FX
                F1 = F1 + t61*t13
                ! FY
                F2 = F2 + t62*t13
                ! FZ
                F3 = F3 + t63*t13
                !
                Fox(ee,i,j,k) = F1
                Foy(ee,i,j,k) = F2
                Foz(ee,i,j,k) = F3

                END_SUBELEM_LOOP()
                !=-=-=-=-=-=-=-=-=-=-
            enddo
        enddo
    enddo
    !=-=-=-=-=-=-=-=-=-=-
#endif
end subroutine PROCNAME

#undef JXX
#undef JXY
#undef JXZ
#undef JYX
#undef JYY
#undef JYZ
#undef JZX
#undef JZY
#undef JZZ

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
