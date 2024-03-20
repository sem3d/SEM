!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
#include "index.h"
#include "loops.h"
#include "solid_forces.h"

#undef PROCNAME
#undef PROCNAME_ATN

#undef ngllvar
!!! Generic values of ngll
#ifndef NGLLVAL
#define ngllvar ngll
#define PROCNAME PROCNAMEBASE()N
#else
!!! Compilation time specified values of ngll
#define PROCNAME PROCNAMEBASE()NGLLVAL
#define ngllvar ngll0
#endif

#define PROCNAME1() PROCNAME
#define SUBPROCNAME1  PROCNAME1()_main_acc
#define SUBPROCNAME2  PROCNAME1()_subloop_acc
#define PROCNAME_ATN  PROCNAME1()_atn_acc

#ifdef ATTENUATION
#error "ATTENUATION defined : wrong file"
#endif
#ifdef NONLINEAR
#error "NONLINEAR defined : wrong file"
#endif

subroutine PROCNAME(dom, ngllvar, var, dvdt)
    use champs_solid
    type(domain_solid), intent (INOUT) :: dom
    type(champssolid), intent(in) :: var
    type(champssolid), intent(inout) :: dvdt
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    !
    integer :: i,j,k,i_dir,ee,idx
    integer :: bnum, nblocks
    real(fpp) :: uval
    !

    nblocks = dom%nblocks

    ! ELASTIC CASE: DISPLACEMENT PREDICTION version not USE_LOCAL_MEM
    !$acc parallel loop collapse(6) async(1)
    do bnum = 0, nblocks-1
        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        do ee = 0, VCHUNK-1
                            idx = dom%Idom_(i,j,k,bnum,ee)
                            uval = var%Depla(idx,i_dir)
                            tDepla_(ee,i,j,k,i_dir,bnum) = uval
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    call SUBPROCNAME2(ngllvar,nblocks, dom%Forces, dom%Depla, dom%props, dom%hprime, dom%gllw, &
        dom%m_Jacob, dom%m_InvGrad)

    !$acc parallel loop collapse(5) async(1)
    do bnum = 0, nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        if ((i==0).or.(j==0).or.(k==0).or.(i==ngll-1).or.(j==ngll-1).or.(k==ngll-1)) then
                            idx = dom%Idom_(i,j,k,bnum,ee)
                            !$acc atomicupdate
                            dvdt%Veloc(idx,0) = dvdt%Veloc(idx,0)-tFox_(ee,i,j,k,bnum)
                            !$acc atomic update
                            dvdt%Veloc(idx,1) = dvdt%Veloc(idx,1)-tFoy_(ee,i,j,k,bnum)
                            !$acc atomic update
                            dvdt%Veloc(idx,2) = dvdt%Veloc(idx,2)-tFoz_(ee,i,j,k,bnum)
                        else
                            idx = dom%Idom_(i,j,k,bnum,ee)
                            dvdt%Veloc(idx,0) = dvdt%Veloc(idx,0)-tFox_(ee,i,j,k,bnum)
                            dvdt%Veloc(idx,1) = dvdt%Veloc(idx,1)-tFoy_(ee,i,j,k,bnum)
                            dvdt%Veloc(idx,2) = dvdt%Veloc(idx,2)-tFoz_(ee,i,j,k,bnum)
                        end if
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
end subroutine PROCNAME

#define x_deriv_ijke(Var,d,dS_dxi,dS_deta,dS_dzeta,a,b,c) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO LL = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Var(ee,LL,b,c,d,bnum)*hprime(LL,a); \
            dS_deta  = dS_deta +Var(ee,a,LL,c,d,bnum)*hprime(LL,b); \
            dS_dzeta = dS_dzeta+Var(ee,a,b,LL,d,bnum)*hprime(LL,c); \
        END DO;

#define load_jacobian_matrix(a,b,c) \
        LJXX = InvGrad(ee,0,0,a,b,c,bnum); \
        LJYX = InvGrad(ee,1,0,a,b,c,bnum); \
        LJZX = InvGrad(ee,2,0,a,b,c,bnum); \
        LJXY = InvGrad(ee,0,1,a,b,c,bnum); \
        LJYY = InvGrad(ee,1,1,a,b,c,bnum); \
        LJZY = InvGrad(ee,2,1,a,b,c,bnum); \
        LJXZ = InvGrad(ee,0,2,a,b,c,bnum); \
        LJYZ = InvGrad(ee,1,2,a,b,c,bnum); \
        LJZZ = InvGrad(ee,2,2,a,b,c,bnum);

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
        xmu = props(ee,a,b,c,cMu,bnum); \
        xla = props(ee,a,b,c,cLambda,bnum); \
        xla2mu = xla + 2. * xmu; \
        sxx = xla2mu * DXX + xla * ( DYY + DZZ ); \
        sxy = xmu * ( DXY + DYX ); \
        sxz = xmu * ( DXZ + DZX ); \
        syy = xla2mu * DYY + xla * ( DXX + DZZ ); \
        syz = xmu * ( DYZ + DZY ); \
        szz = xla2mu * DZZ + xla * ( DXX + DYY );

#define CC(d) props(ee,i,j,k,cAniso+d,bnum)
#ifdef ANISO
#define nprops 26
#else
#ifdef ATTENUATION
#define nprops 5
#else
#define nprops 2
#endif
#endif

subroutine SUBPROCNAME2(ngllvar,nblocks,Forces,Depla,props, hprime, gllw, Jacob, invgrad)
    implicit none

    integer, intent(in) :: ngllvar
    integer, intent(in) :: nblocks
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nblocks-1), intent(out) :: Forces
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nblocks-1), intent(in) :: Depla
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nprops,0:nblocks-1), intent(in) :: props
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: jacob
    real(fpp), dimension(0:VCHUNK-1,0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: invgrad
    real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
    real(fpp), dimension(0:ngll-1), intent(in) :: gllw
    integer :: i,j,k,l,ee,ll
    integer :: bnum
    real(fpp) :: wi,wj,wk,wl,ljxx,ljxy,ljxz,ljyx,ljyy,ljyz,ljzx,ljzy,ljzz
    real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real(fpp) :: dUx_dxi, dUx_deta, dUx_dzeta
    real(fpp) :: dUy_dxi, dUy_deta, dUy_dzeta
    real(fpp) :: dUz_dxi, dUz_deta, dUz_dzeta
    real(fpp) :: sxx,sxy,sxz,syy,syz,szz
    real(fpp) :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10,t4
    real(fpp) :: xla2mu, xla, jac
    real(fpp) :: xmu, xtmp, ytmp, ztmp
#ifdef ANISO
    real(fpp) :: Eyz, Exz, Exy
#endif
    ! TODO: CHECK

    !$acc parallel loop async(1) collapse(5)
    do bnum = 0, nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP1(ee)
                    Forces(ee,i,j,k,0,bnum) = 0
                    Forces(ee,i,j,k,1,bnum) = 0
                    Forces(ee,i,j,k,2,bnum) = 0
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    enddo

    !$acc parallel loop async(1) collapse(5)
    do bnum = 0, nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP1(ee)
                    compute_du_deta(i,j,k)
                    load_jacobian_matrix(i,j,k)
                    compute_du_dx(i,j,k)
#ifdef ANISO
                    Eyz = (DYZ+DZY)*M_SQRT1_2
                    Exz = (DXZ+DZX)*M_SQRT1_2
                    Exy = (DXY+DYX)*M_SQRT1_2
                    sxx = DXX*CC( 0) + DYY*CC( 1) + DZZ*CC( 2) + EYZ*CC( 3) + EXZ*CC( 4) + EXY*CC( 5)
                    syy = DXX*CC( 1) + DYY*CC( 6) + DZZ*CC( 7) + EYZ*CC( 8) + EXZ*CC( 9) + EXY*CC(10)
                    szz = DXX*CC( 2) + DYY*CC( 7) + DZZ*CC(11) + EYZ*CC(12) + EXZ*CC(13) + EXY*CC(14)
                    syz = DXX*CC( 3) + DYY*CC( 8) + DZZ*CC(12) + EYZ*CC(15) + EXZ*CC(16) + EXY*CC(17)
                    sxz = DXX*CC( 4) + DYY*CC( 9) + DZZ*CC(13) + EYZ*CC(16) + EXZ*CC(18) + EXY*CC(19)
                    sxy = DXX*CC( 5) + DYY*CC(10) + DZZ*CC(14) + EYZ*CC(17) + EXZ*CC(19) + EXY*CC(20)
                    syz = syz/M_SQRT2
                    sxz = sxz/M_SQRT2
                    sxy = sxy/M_SQRT2
#else
                    compute_sigma(i,j,k)
#endif
                    !=====================
                    !       FX
                    !=====================
                    !
                    xt1 = sxx*LJXX+sxy*LJYX+sxz*LJZX
                    xt2 = sxx*LJXY+sxy*LJYY+sxz*LJZY
                    xt3 = sxx*LJXZ+sxy*LJYZ+sxz*LJZZ
                    !
                    !=====================
                    !       FY
                    !=====================
                    !
                    xt5 = syy*LJYX+sxy*LJXX+syz*LJZX
                    xt6 = syy*LJYY+sxy*LJXY+syz*LJZY
                    xt7 = syy*LJYZ+sxy*LJXZ+syz*LJZZ
                    !
                    !=====================
                    !       FZ
                    !=====================
                    !
                    xt8 = szz*LJZX+sxz*LJXX+syz*LJYX
                    xt9 = szz*LJZY+sxz*LJXY+syz*LJYY
                    xt10= szz*LJZZ+sxz*LJXZ+syz*LJYZ
                    !
                    !- Multiplication par le Jacobien et le poids d'integration
                    !
                    wi = Gllw(i)
                    wj = Gllw(j)
                    wk = Gllw(k)
                    jac = Jacob(ee,i,j,k,bnum)
                    t4 = jac * wi * wj * wk
                    xt1  =  xt1 * t4
                    xt5  =  xt5 * t4
                    xt8  =  xt8 * t4
                    xt2  =  xt2 * t4
                    xt6  =  xt6 * t4
                    xt9  =  xt9 * t4
                    xt3  =  xt3 * t4
                    xt7  =  xt7 * t4
                    xt10 = xt10 * t4

                    DO L = 0, ngll-1
                        xtmp = hprime(L,I)*xt1
                        ytmp = hprime(L,I)*xt5
                        ztmp = hprime(L,I)*xt8
                        !$acc atomic update
                        Forces(ee,L,j,k,0,bnum) = Forces(ee,L,j,k,0,bnum) + xtmp
                        !$acc atomic update
                        Forces(ee,L,j,k,1,bnum) = Forces(ee,L,j,k,1,bnum) + ytmp
                        !$acc atomic update
                        Forces(ee,L,j,k,2,bnum) = Forces(ee,L,j,k,2,bnum) + ztmp
                    END DO
                    DO L = 0, ngll-1
                        xtmp = hprime(L,J)*xt2
                        ytmp = hprime(L,J)*xt6
                        ztmp = hprime(L,J)*xt9
                        !$acc atomic update
                        Forces(ee,i,L,k,0,bnum) = Forces(ee,i,L,k,0,bnum) + xtmp
                        !$acc atomic update
                        Forces(ee,i,L,k,1,bnum) = Forces(ee,i,L,k,1,bnum) + ytmp
                        !$acc atomic update
                        Forces(ee,i,L,k,2,bnum) = Forces(ee,i,L,k,2,bnum) + ztmp
                    END DO
                    DO L = 0, ngll-1
                        xtmp = hprime(L,K)*xt3 
                        ytmp = hprime(L,K)*xt7 
                        ztmp = hprime(L,K)*xt10
                        !$acc atomic update
                        Forces(ee,i,j,L,0,bnum) = Forces(ee,i,j,L,0,bnum) + xtmp
                        !$acc atomic update
                        Forces(ee,i,j,L,1,bnum) = Forces(ee,i,j,L,1,bnum) + ytmp
                        !$acc atomic update
                        Forces(ee,i,j,L,2,bnum) = Forces(ee,i,j,L,2,bnum) + ztmp
                    END DO
                    END_SUBELEM_LOOP()
                end do
            end do
        end do
    end do
    !
end subroutine SUBPROCNAME2

#undef x_deriv_ijke
#undef load_jacobian_matrix
#undef compute_du_dx
#undef compute_sigma


!!!! ROUTINES DE TRAITEMENT DES MIRROIRS
#define PROC_MIRROR         PROCNAME1()_mirror

subroutine PROC_MIRROR(dom, ngllvar, var, dvdt, m_dump, m_expl, m_recalc)
    use champs_solid
    type(domain_solid), intent (INOUT) :: dom
    type(champssolid), intent(in) :: var
    type(champssolid), intent(inout) :: dvdt
    integer, intent(in) :: ngllvar
    logical, intent(IN) :: m_dump, m_expl, m_recalc
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    !
    integer :: i,j,k,i_dir,ee,idx, idx_m
    integer :: bnum, nblocks, lnum
    logical :: m_load
    real(fpp),dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: Fox_m,Foy_m,Foz_m

    m_load = .not. m_dump
    nblocks = dom%nblocks

end subroutine PROC_MIRROR

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
