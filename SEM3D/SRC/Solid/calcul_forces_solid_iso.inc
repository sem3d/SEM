!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
#include "index.h"
#include "loops.h"

#undef PROCNAME
#undef PROCNAME_ATN
#define PROCNAME PROCNAMEBASE()NGLLVAL
#define PROCNAME_ATN PROCNAMEBASE_ATN()NGLLVAL
#define PROCNAME1() PROCNAME
#define SUBPROCNAME  PROCNAME1()_subloop

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

#define x_deriv_ijke(Var,d,dS_dxi,dS_deta,dS_dzeta,a,b,c) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO LL = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Var(ee,LL,b,c,d)*dom%hprime(LL,a); \
            dS_deta  = dS_deta +Var(ee,a,LL,c,d)*dom%hprime(LL,b); \
            dS_dzeta = dS_dzeta+Var(ee,a,b,LL,d)*dom%hprime(LL,c); \
        END DO;

#define load_jacobian_matrix(a,b,c) \
        LJXX = dom%InvGrad_(0,0,a,b,c,bnum,ee); \
        LJYX = dom%InvGrad_(1,0,a,b,c,bnum,ee); \
        LJZX = dom%InvGrad_(2,0,a,b,c,bnum,ee); \
        LJXY = dom%InvGrad_(0,1,a,b,c,bnum,ee); \
        LJYY = dom%InvGrad_(1,1,a,b,c,bnum,ee); \
        LJZY = dom%InvGrad_(2,1,a,b,c,bnum,ee); \
        LJXZ = dom%InvGrad_(0,2,a,b,c,bnum,ee); \
        LJYZ = dom%InvGrad_(1,2,a,b,c,bnum,ee); \
        LJZZ = dom%InvGrad_(2,2,a,b,c,bnum,ee);

#define compute_du_deta(a,b,c) \
        x_deriv_ijke(Depla,0,dUx_dxi,dUx_deta,dUx_dzeta,a,b,c); \
        x_deriv_ijke(Depla,1,dUy_dxi,dUy_deta,dUy_dzeta,a,b,c); \
        x_deriv_ijke(Depla,2,dUz_dxi,dUz_deta,dUz_dzeta,a,b,c);

#define compute_du_dx(a,b,c) \
        dxx=dUx_dxi*LJXX + dUx_deta*LJXY + dUx_dzeta*LJXZ; \
        dxy=dUx_dxi*LJYX + dUx_deta*LJYY + dUx_dzeta*LJYZ; \
        dxz=dUx_dxi*LJZX + dUx_deta*LJZY + dUx_dzeta*LJZZ; \
        dyx=dUy_dxi*LJXX + dUy_deta*LJXY + dUy_dzeta*LJXZ; \
        dyy=dUy_dxi*LJYX + dUy_deta*LJYY + dUy_dzeta*LJYZ; \
        dyz=dUy_dxi*LJZX + dUy_deta*LJZY + dUy_dzeta*LJZZ; \
        dzx=dUz_dxi*LJXX + dUz_deta*LJXY + dUz_dzeta*LJXZ; \
        dzy=dUz_dxi*LJYX + dUz_deta*LJYY + dUz_dzeta*LJYZ; \
        dzz=dUz_dxi*LJZX + dUz_deta*LJZY + dUz_dzeta*LJZZ;

#define compute_sigma(a,b,c) \
        xmu = dom%Mu_    (a,b,c,bnum,ee); \
        xla = dom%Lambda_(a,b,c,bnum,ee); \
        xla2mu = xla + 2. * xmu; \
        sxx = xla2mu * DXX + xla * ( DYY + DZZ ); \
        sxy = xmu * ( DXY + DYX ); \
        sxz = xmu * ( DXZ + DZX ); \
        syy = xla2mu * DYY + xla * ( DXX + DZZ ); \
        syz = xmu * ( DYZ + DZY ); \
        szz = xla2mu * DZZ + xla * ( DXX + DYY );

#define get_sigma(a,b,c) \
        sxx = Sigma(ee,a,b,c,0);\
        syy = Sigma(ee,a,b,c,1);\
        szz = Sigma(ee,a,b,c,2);\
        syz = Sigma(ee,a,b,c,3);\
        sxz = Sigma(ee,a,b,c,4);\
        sxy = Sigma(ee,a,b,c,5);


subroutine PROCNAME(dom,ngllvar,bnum,Fox,Foy,Foz,Depla, Sigma)
    !$acc routine worker
    use champs_solid
    implicit none

    type(domain_solid), intent (INOUT) :: dom
    integer, intent(in) :: bnum
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: Fox,Foz,Foy
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(in) :: Depla
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:5), intent(inout) :: Sigma
    integer :: i,j,k,l,e,ee,ll
    real(fpp) :: wi,wj,wk,wl,ljxx,ljxy,ljxz,ljyx,ljyy,ljyz,ljzx,ljzy,ljzz
    real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real(fpp) :: dUx_dxi, dUx_deta, dUx_dzeta
    real(fpp) :: dUy_dxi, dUy_deta, dUy_dzeta
    real(fpp) :: dUz_dxi, dUz_deta, dUz_dzeta
    real(fpp) :: sxx,sxy,sxz,syy,syz,szz

    real(fpp) :: xla2mu, xla
    real(fpp) :: xmu
    ! TODO: CHECK
!    real(fpp), dimension(0:ngll-1,0:ngll-1) :: hprime
!    real(fpp), dimension(0:ngll-1) :: gllw
!    hprime = dom%hprime
!    gllw = dom%gllw
    !if (ngll/=ngllvar) stop 1

    !$acc loop worker collapse(3)
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                Fox(ee,i,j,k) = 0
                Foy(ee,i,j,k) = 0
                Foz(ee,i,j,k) = 0
                END_SUBELEM_LOOP()
                ! Compute sigma
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                load_jacobian_matrix(i,j,k)
                compute_du_deta(i,j,k)
                compute_du_dx(i,j,k)
                compute_sigma(i,j,k)
                ! xt1 +xt5 + xt8
                sigma(ee,i,j,k,0) = sxx
                sigma(ee,i,j,k,1) = syy
                sigma(ee,i,j,k,2) = szz
                sigma(ee,i,j,k,3) = syz
                sigma(ee,i,j,k,4) = sxz
                sigma(ee,i,j,k,5) = sxy
                END_SUBELEM_LOOP()
            end do
        end do
    end do

    !$acc loop worker collapse(3)
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                wi = dom%GLLw(i)
                wj = dom%GLLw(j)
                wk = dom%GLLw(k)
                ! Compute the terms relative to local coord xsi
                DO L = 0, ngll-1
                    LOOP_VECTORIZE
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    wl = dom%hprime(i,L)*dom%Jacob_(L,j,k,bnum,ee)*dom%GLLw(L)*wj*wk
                    load_jacobian_matrix(L,j,k)
                    get_sigma(L,j,k)
                    ! xt1 +xt5 + xt8
                    Fox(ee,i,j,k) = Fox(ee,i,j,k) + (sxx*ljxx+sxy*ljyx+sxz*ljzx)*wl
                    Foy(ee,i,j,k) = Foy(ee,i,j,k) + (sxy*ljxx+syy*ljyx+syz*ljzx)*wl
                    Foz(ee,i,j,k) = Foz(ee,i,j,k) + (sxz*ljxx+syz*ljyx+szz*ljzx)*wl
                    END_SUBELEM_LOOP()
                END DO

                ! Compute the terms relative to local coord eta
                DO L = 0, ngll-1
                    LOOP_VECTORIZE
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    wl = dom%hprime(j,L)*dom%Jacob_(i,L,k,bnum,ee)*dom%GLLw(L)*wi*wk
                    load_jacobian_matrix(i,L,k)
                    get_sigma(i,L,k)
                    ! xt2 + xt6 + xt9
                    Fox(ee,i,j,k) = Fox(ee,i,j,k) + (sxx*ljxy+sxy*ljyy+sxz*ljzy)*wl
                    Foy(ee,i,j,k) = Foy(ee,i,j,k) + (sxy*ljxy+syy*ljyy+syz*ljzy)*wl
                    Foz(ee,i,j,k) = Foz(ee,i,j,k) + (sxz*ljxy+syz*ljyy+szz*ljzy)*wl
                    END_SUBELEM_LOOP()
                END DO

                DO L = 0, ngll-1
                    LOOP_VECTORIZE
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    wl = dom%hprime(k,L)*dom%Jacob_(i,j,L,bnum,ee)*dom%GLLw(L)*wi*wj
                    load_jacobian_matrix(i,j,L)
                    get_sigma(i,j,L)
                    Fox(ee,i,j,k) = Fox(ee,i,j,k) + (sxx*ljxz+sxy*ljyz+sxz*ljzz)*wl
                    Foy(ee,i,j,k) = Foy(ee,i,j,k) + (sxy*ljxz+syy*ljyz+syz*ljzz)*wl
                    Foz(ee,i,j,k) = Foz(ee,i,j,k) + (sxz*ljxz+syz*ljyz+szz*ljzz)*wl
                    END_SUBELEM_LOOP()
                END DO
            enddo
        enddo
    enddo
    !
end subroutine PROCNAME

#undef x_deriv_ijke
#undef load_jacobian_matrix
#undef compute_du_dx
#undef compute_sigma

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
