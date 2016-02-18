!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"

#ifndef NGLLVAL
#error "Must define NGLLVAL before inclusion"
#endif

#ifdef NGLL_GEN
#define ngllvar ngll
#else
#define ngllvar ngll0
#endif

#ifndef PROCNAME
#define CONCAT(a,b) a##b
#define PROCNAME CONCAT(calcul_forces_fluid_,NGLLVAL)
#endif

subroutine PROCNAME(dom,ngllvar,lnum,FFl,Phi)
    use sdomain
    use deriv3d
    implicit none

    type(domain_fluid), intent (INOUT) :: dom
    integer, intent(in) :: ngllvar
#ifndef NGLL_GEN
    integer, parameter :: ngll=NGLLVAL
#endif
    integer, intent(in) :: lnum
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: FFl
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Phi
    real(fpp) :: dPhi_dX,dPhi_dY,dPhi_dZ
    real(fpp) :: dPhi_dxi,dPhi_deta,dPhi_dzeta
    real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    integer :: i,j,k,l,e,ee
    real(fpp) :: sx,sy,sz,t4,F1
    real(fpp) :: t41,t11,t51,t12,t61,t13
    real(fpp) :: xt1,xt6,xt10
    real(fpp), parameter :: zero = 0.
    real(fpp) :: xdens
#ifdef ORIG
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t1,t6,t10
#endif
    !dir$ ASSUME_ALIGNED PHI: 32
    !dir$ ASSUME_ALIGNED FFL: 32
    !!        !dir$ ASSUME_ALIGNED dom%m_InvGrad: 32
    ! Assume FFl = 0d0
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#ifdef SEM_VEC
                !dir$ simd
#endif
                BEGIN_SUBELEM_LOOP(e,ee,lnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                dPhi_dxi   = 0D0
                dPhi_deta  = 0D0
                dPhi_dzeta = 0D0
                !dir$ loop count (5)
                DO L = 0, ngll-1
                    dPhi_dxi   = dPhi_dxi  +Phi(ee,L,J,K)*dom%hprime(L,I)
                    dPhi_deta  = dPhi_deta +Phi(ee,I,L,K)*dom%hprime(L,J)
                    dPhi_dzeta = dPhi_dzeta+Phi(ee,I,J,L)*dom%hprime(L,K)
                END DO

                xi1 = dom%InvGrad_(0,0,i,j,k,e)
                xi2 = dom%InvGrad_(1,0,i,j,k,e)
                xi3 = dom%InvGrad_(2,0,i,j,k,e)
                et1 = dom%InvGrad_(0,1,i,j,k,e)
                et2 = dom%InvGrad_(1,1,i,j,k,e)
                et3 = dom%InvGrad_(2,1,i,j,k,e)
                ga1 = dom%InvGrad_(0,2,i,j,k,e)
                ga2 = dom%InvGrad_(1,2,i,j,k,e)
                ga3 = dom%InvGrad_(2,2,i,j,k,e)
                !- in the physical domain
                dPhi_dx = dPhi_dxi*xi1 + dPhi_deta*et1 + dPhi_dzeta*ga1
                dPhi_dy = dPhi_dxi*xi2 + dPhi_deta*et2 + dPhi_dzeta*ga2
                dPhi_dz = dPhi_dxi*xi3 + dPhi_deta*et3 + dPhi_dzeta*ga3

                ! (fluid equivalent) stress  ( = physical velocity)
                xdens = 1d0/dom%Density_(i,j,k,e)
                sx = xdens*dPhi_dX
                sy = xdens*dPhi_dY
                sz = xdens*dPhi_dZ


                !=====================
                !       F1 
                xt1 = sx*xi1 + sy*xi2 + sz*xi3

                !=====================
                !       F2 
                xt6 = sx*et1 + sy*et2 + sz*et3

                !=====================
                !       F3 
                xt10 = sx*ga1 + sy*ga2 + sz*ga3

                !
                !- Multiply par Jacobian and weight
                !
                t4  = dom%Jacob_(i,j,k,e) * dom%GLLw(i)
                xt1 = xt1 * t4

                t4  = dom%Jacob_(i,j,k,e) * dom%GLLw(j)
                xt6 = xt6 * t4

                t4   = dom%Jacob_(i,j,k,e) * dom%GLLw(k)
                xt10 = xt10 * t4

#ifdef ORIG
                t1(ee,i,j,k) = xt1
                t6(ee,j,i,k) = xt6
                t10(ee,k,i,j) = xt10
#else
                DO L = 0, ngll-1
                    FFl(ee,L,j,k) = FFl(ee,L,j,k) + dom%hprime(L,I)*xt1*dom%GLLw(j) * dom%GLLw(k)
                    FFl(ee,i,L,k) = FFl(ee,i,L,k) + dom%hprime(L,J)*xt6*dom%GLLw(i) * dom%GLLw(k)
                    FFl(ee,i,j,L) = FFl(ee,i,j,L) + dom%hprime(L,K)*xt10*dom%GLLw(i) * dom%GLLw(j)
                END DO
#endif
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo

    !
    !- Multiplication par la matrice de derivation puis par les poids
    !
#ifdef ORIG
    !=-=-=-=-=-=-=-=-=-=-
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#ifdef SEM_VEC
                !dir$ simd
#endif
                BEGIN_SUBELEM_LOOP(e,ee,lnum)
                !=-=-=-=-=-=-=-=-=-=-
                !
                t11 = dom%GLLw(j) * dom%GLLw(k)
                t12 = dom%GLLw(i) * dom%GLLw(k)
                t13 = dom%GLLw(i) * dom%GLLw(j)
                !
                t41 = zero
                t51 = zero
                t61 = zero
                !
                !dir$ loop count (5)
                do l = 0,ngll-1
                    t41 = t41 + dom%htprime(l,i) * t1(ee,l,j,k)
                enddo

                !dir$ loop count (5)
                do l = 0,ngll-1
                    t51 = t51 + dom%htprime(l,j) * t6(ee,l,i,k)
                enddo
                ! FFl
                F1 = t41*t11 + t51*t12
                !
                !
                !dir$ loop count (5)
                do l = 0,ngll-1
                    t61 = t61 + dom%htprime(l,k) * t10(ee,l,i,j)
                enddo

                ! FX
                F1 = F1 + t61*t13
                !
                FFl(ee,i,j,k) = F1
                !=-=-=-=-=-=-=-=-=-=-
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo
#endif
    !=-=-=-=-=-=-=-=-=-=-
end subroutine PROCNAME



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
