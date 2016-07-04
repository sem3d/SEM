!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"

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

subroutine PROCNAME(ngllvar,nblocks,bnum,hprime,htprime,gllw,m_InvGrad,m_Jacob,m_IDensity,FFl,Phi)
    use sdomain
    use deriv3d
    implicit none

    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    integer, intent(in) :: bnum,nblocks
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: FFl
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Phi
    real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime,htprime
    real(fpp), dimension(0:ngll-1), intent(in) :: gllw
    real(fpp), dimension(IND_MNIJKE(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)), intent(in) :: m_InvGrad
    real(fpp), dimension(IND_IJKE(0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)), intent(in) :: m_Jacob
    real(fpp), dimension(IND_IJKE(0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)), intent(in) :: m_IDensity
    !
    real(fpp) :: dPhi_dX,dPhi_dY,dPhi_dZ
    real(fpp) :: dPhi_dxi,dPhi_deta,dPhi_dzeta
    real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    integer :: i,j,k,l,ei,ee
    real(fpp) :: sx,sy,sz,t4
    real(fpp) :: xt1,xt6,xt10
    !real(fpp), parameter :: zero = 0d0
    real(fpp) :: xdens
#ifdef USE_LOCAL_MEM
    real(fpp) :: t41,t11,t51,t12,t61,t13,F1
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t1,t6,t10
#endif
    !dir$ ASSUME_ALIGNED PHI: 32
    !dir$ ASSUME_ALIGNED FFL: 32
    !dir$ ASSUME_ALIGNED m_InvGrad: 32
    !dir$ ASSUME_ALIGNED m_Jacob: 32
    !dir$ ASSUME_ALIGNED m_IDensity: 32
    !dir$ ASSUME_ALIGNED hprime: 32
    !dir$ ASSUME_ALIGNED htprime: 32
    !dir$ ASSUME_ALIGNED gllw: 32
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(ee)
#endif
                BEGIN_SUBELEM_LOOP(ei,ee,bnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                dPhi_dxi   = 0D0
                dPhi_deta  = 0D0
                dPhi_dzeta = 0D0
                !dir$ loop count (5)
                DO L = 0, ngll-1
                    dPhi_dxi   = dPhi_dxi  +Phi(ee,L,J,K)*hprime(L,I)
                    dPhi_deta  = dPhi_deta +Phi(ee,I,L,K)*hprime(L,J)
                    dPhi_dzeta = dPhi_dzeta+Phi(ee,I,J,L)*hprime(L,K)
                END DO

                xi1 = InvGrad_(0,0,i,j,k,bnum,ee)
                xi2 = InvGrad_(1,0,i,j,k,bnum,ee)
                xi3 = InvGrad_(2,0,i,j,k,bnum,ee)
                et1 = InvGrad_(0,1,i,j,k,bnum,ee)
                et2 = InvGrad_(1,1,i,j,k,bnum,ee)
                et3 = InvGrad_(2,1,i,j,k,bnum,ee)
                ga1 = InvGrad_(0,2,i,j,k,bnum,ee)
                ga2 = InvGrad_(1,2,i,j,k,bnum,ee)
                ga3 = InvGrad_(2,2,i,j,k,bnum,ee)
                !- in the physical domain
                dPhi_dx = dPhi_dxi*xi1 + dPhi_deta*et1 + dPhi_dzeta*ga1
                dPhi_dy = dPhi_dxi*xi2 + dPhi_deta*et2 + dPhi_dzeta*ga2
                dPhi_dz = dPhi_dxi*xi3 + dPhi_deta*et3 + dPhi_dzeta*ga3

                ! (fluid equivalent) stress  ( = physical velocity)
                xdens = IDensity_(i,j,k,bnum,ee)
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
                t4  = Jacob_(i,j,k,bnum,ee) * gllw(i)
                xt1 = xt1 * t4

                t4  = Jacob_(i,j,k,bnum,ee) * gllw(j)
                xt6 = xt6 * t4

                t4   = Jacob_(i,j,k,bnum,ee) * gllw(k)
                xt10 = xt10 * t4

#ifdef USE_LOCAL_MEM
                t1(ee,i,j,k) = xt1
                t6(ee,j,i,k) = xt6
                t10(ee,k,i,j) = xt10
#else
                DO L = 0, ngll-1
                    FFl(ee,L,j,k) = FFl(ee,L,j,k) + hprime(L,I)*xt1*gllw(j) * gllw(k)
                    FFl(ee,i,L,k) = FFl(ee,i,L,k) + hprime(L,J)*xt6*gllw(i) * gllw(k)
                    FFl(ee,i,j,L) = FFl(ee,i,j,L) + hprime(L,K)*xt10*gllw(i) * gllw(j)
                ENDDO
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
!$omp simd linear(ee)
#endif
                BEGIN_SUBELEM_LOOP(ei,ee,bnum)
                !=-=-=-=-=-=-=-=-=-=-
                !
                t11 = gllw(j) * gllw(k)
                t12 = gllw(i) * gllw(k)
                t13 = gllw(i) * gllw(j)
                !
                t41 = zero
                t51 = zero
                t61 = zero
                !
                !dir$ loop count (5)
                do l = 0,ngll-1
                    t41 = t41 + htprime(l,i) * t1(ee,l,j,k)
                enddo

                !dir$ loop count (5)
                do l = 0,ngll-1
                    t51 = t51 + htprime(l,j) * t6(ee,l,i,k)
                enddo
                ! FFl
                F1 = t41*t11 + t51*t12
                !
                !
                !dir$ loop count (5)
                do l = 0,ngll-1
                    t61 = t61 + htprime(l,k) * t10(ee,l,i,j)
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
