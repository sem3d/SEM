!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
#include "index.h"
#include "loops.h"

#ifndef PROCNAME
#error "Must define PROCNAME before inclusion"
#endif

#undef ngllvar
#ifndef NGLLVAL
#define ngllvar ngll
#else
#define ngllvar ngll0
#endif

#ifdef ATTENUATION
#error "ATTENUATION defined : wrong file"
#endif
#ifdef NONLINEAR
#error "NONLINEAR defined : wrong file"
#endif
#ifdef ANISO
#error "ANISO defined : wrong file"
#endif

subroutine PROCNAME(dom,ngllvar,bnum,Fox,Foy,Foz,Depla)
    !$acc routine worker
    use champs_solid
    use nonlinear
    implicit none

    type(domain_solid), intent (INOUT) :: dom
    integer, intent(in) :: bnum
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: Fox,Foz,Foy
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(in) :: Depla
    integer :: i,j,k,l,e,ee
    real(fpp) :: wi,wj,wk,wl,ljxx,ljxy,ljxz,ljyx,ljyy,ljyz,ljzx,ljzy,ljzz
    real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real(fpp), dimension(0:VCHUNK-1) :: dUx_dxi, dUx_deta, dUx_dzeta
    real(fpp), dimension(0:VCHUNK-1) :: dUy_dxi, dUy_deta, dUy_dzeta
    real(fpp), dimension(0:VCHUNK-1) :: dUz_dxi, dUz_deta, dUz_dzeta
    real(fpp), dimension(0:VCHUNK-1) :: fxijk, fyijk, fzijk
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: sxx,sxy,sxz,syy,syz,szz

    real(fpp) :: xla2mu, xla
    real(fpp) :: xmu

    !if (ngll/=ngllvar) stop 1

    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                local_deriv_ijke(Depla,0,dUx_dxi(ee),dUx_deta(ee),dUx_dzeta(ee))
                END_SUBELEM_LOOP()
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                local_deriv_ijke(Depla,1,dUy_dxi(ee),dUy_deta(ee),dUy_dzeta(ee))
                END_SUBELEM_LOOP()
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                local_deriv_ijke(Depla,2,dUz_dxi(ee),dUz_deta(ee),dUz_dzeta(ee))
                END_SUBELEM_LOOP()

#define JXX dom%InvGrad_(0,0,i,j,k,bnum,ee)
#define JXY dom%InvGrad_(0,1,i,j,k,bnum,ee)
#define JXZ dom%InvGrad_(0,2,i,j,k,bnum,ee)
#define JYX dom%InvGrad_(1,0,i,j,k,bnum,ee)
#define JYY dom%InvGrad_(1,1,i,j,k,bnum,ee)
#define JYZ dom%InvGrad_(1,2,i,j,k,bnum,ee)
#define JZX dom%InvGrad_(2,0,i,j,k,bnum,ee)
#define JZY dom%InvGrad_(2,1,i,j,k,bnum,ee)
#define JZZ dom%InvGrad_(2,2,i,j,k,bnum,ee)

                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)

                dxx=dUx_dxi(ee)*JXX + dUx_deta(ee)*JXY + dUx_dzeta(ee)*JXZ
                dxy=dUx_dxi(ee)*JYX + dUx_deta(ee)*JYY + dUx_dzeta(ee)*JYZ
                dxz=dUx_dxi(ee)*JZX + dUx_deta(ee)*JZY + dUx_dzeta(ee)*JZZ

                dyx=dUy_dxi(ee)*JXX + dUy_deta(ee)*JXY + dUy_dzeta(ee)*JXZ
                dyy=dUy_dxi(ee)*JYX + dUy_deta(ee)*JYY + dUy_dzeta(ee)*JYZ
                dyz=dUy_dxi(ee)*JZX + dUy_deta(ee)*JZY + dUy_dzeta(ee)*JZZ

                dzx=dUz_dxi(ee)*JXX + dUz_deta(ee)*JXY + dUz_dzeta(ee)*JXZ
                dzy=dUz_dxi(ee)*JYX + dUz_deta(ee)*JYY + dUz_dzeta(ee)*JYZ
                dzz=dUz_dxi(ee)*JZX + dUz_deta(ee)*JZY + dUz_dzeta(ee)*JZZ


                ! LINEAR ELASTIC ISOTROPE
                xmu = dom%Mu_    (i,j,k,bnum,ee)
                xla = dom%Lambda_(i,j,k,bnum,ee)
                xla2mu = xla + 2. * xmu
                sxx(ee,i,j,k) = xla2mu * DXX + xla * ( DYY + DZZ )
                sxy(ee,i,j,k) = xmu * ( DXY + DYX )
                sxz(ee,i,j,k) = xmu * ( DXZ + DZX )
                syy(ee,i,j,k) = xla2mu * DYY + xla * ( DXX + DZZ )
                syz(ee,i,j,k) = xmu * ( DYZ + DZY )
                szz(ee,i,j,k) = xla2mu * DZZ + xla * ( DXX + DYY )
                END_SUBELEM_LOOP()
            end do
        end do
    end do


    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                fxijk = 0
                fyijk = 0
                fzijk = 0
                wi = dom%GLLw(i)
                wj = dom%GLLw(j)
                wk = dom%GLLw(k)
                ! Compute the terms relative to local coord xsi
                DO L = 0, ngll-1
                    wl = dom%hprime(i,L)*dom%Jacob_(L,j,k,bnum,ee)*dom%GLLw(L)*wj*wk
                    LOOP_VECTORIZE
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    LJXX = dom%InvGrad_(0,0,L,j,k,bnum,ee)
                    LJYX = dom%InvGrad_(1,0,L,j,k,bnum,ee)
                    LJZX = dom%InvGrad_(2,0,L,j,k,bnum,ee)
                    ! xt1 +xt5 + xt8
                    fxijk(ee) = fxijk(ee) + (sxx(ee,L,j,k)*ljxx+sxy(ee,L,j,k)*ljyx+sxz(ee,L,j,k)*ljzx)*wl
                    fyijk(ee) = fyijk(ee) + (sxy(ee,L,j,k)*ljxx+syy(ee,L,j,k)*ljyx+syz(ee,L,j,k)*ljzx)*wl
                    fzijk(ee) = fzijk(ee) + (sxz(ee,L,j,k)*ljxx+syz(ee,L,j,k)*ljyx+szz(ee,L,j,k)*ljzx)*wl
                    END_SUBELEM_LOOP()
                END DO


                ! Compute the terms relative to local coord eta
                DO L = 0, ngll-1
                    wl = dom%hprime(j,L)*dom%Jacob_(i,L,k,bnum,ee)*dom%GLLw(L)*wi*wk
                    LOOP_VECTORIZE
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    LJXY = dom%InvGrad_(0,1,i,L,k,bnum,ee)
                    LJYY = dom%InvGrad_(1,1,i,L,k,bnum,ee)
                    LJZY = dom%InvGrad_(2,1,i,L,k,bnum,ee)
                    ! xt2 + xt6 + xt9
                    fxijk(ee) = fxijk(ee) + (sxx(ee,i,L,k)*ljxy+sxy(ee,i,L,k)*ljyy+sxz(ee,i,L,k)*ljzy)*wl
                    fyijk(ee) = fyijk(ee) + (sxy(ee,i,L,k)*ljxy+syy(ee,i,L,k)*ljyy+syz(ee,i,L,k)*ljzy)*wl
                    fzijk(ee) = fzijk(ee) + (sxz(ee,i,L,k)*ljxy+syz(ee,i,L,k)*ljyy+szz(ee,i,L,k)*ljzy)*wl
                    END_SUBELEM_LOOP()
                END DO

                DO L = 0, ngll-1
                    wl = dom%hprime(k,L)*dom%Jacob_(i,j,L,bnum,ee)*dom%GLLw(L)*wi*wj
                    LOOP_VECTORIZE
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    LJXZ = dom%InvGrad_(0,2,i,j,L,bnum,ee)
                    LJYZ = dom%InvGrad_(1,2,i,j,L,bnum,ee)
                    LJZZ = dom%InvGrad_(2,2,i,j,L,bnum,ee)
                    fxijk(ee) = fxijk(ee) + (sxx(ee,i,j,L)*ljxz+sxy(ee,i,j,L)*ljyz+sxz(ee,i,j,L)*ljzz)*wl
                    fyijk(ee) = fyijk(ee) + (sxy(ee,i,j,L)*ljxz+syy(ee,i,j,L)*ljyz+syz(ee,i,j,L)*ljzz)*wl
                    fzijk(ee) = fzijk(ee) + (sxz(ee,i,j,L)*ljxz+syz(ee,i,j,L)*ljyz+szz(ee,i,j,L)*ljzz)*wl
                    END_SUBELEM_LOOP()
                END DO
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                Fox(ee,i,j,k) = fxijk(ee)
                Foy(ee,i,j,k) = fyijk(ee)
                Foz(ee,i,j,k) = fzijk(ee)
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo
    !
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
