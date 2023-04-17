!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
#include "index.h"
#include "loops.h"


#undef PROCNUM
#undef ngllvar
#ifndef NGLLVAL
#define ngllvar ngll
#define PROCNUM() N
#else
#define ngllvar ngll0
#define PROCNUM() NGLLVAL
#endif

#undef PROCNAME
#undef PROCNAME1
#undef PROC_PRED
#undef PROC_MAIN

#define PROCNAME PROCNAMEBASE()PROCNUM()
#define PROCNAME1() PROCNAME
#define PROC_PRED  PROCNAME1()_pred
#define PROC_MAIN  PROCNAME1()_main

#undef x_deriv_ijke
#define x_deriv_ijke(Var,d,dS_dxi,dS_deta,dS_dzeta,dxx,dxy,dxz) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO L = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Var(ee,L,J,K,d)*hprime(L,I); \
            dS_deta  = dS_deta +Var(ee,I,L,K,d)*hprime(L,J); \
            dS_dzeta = dS_dzeta+Var(ee,I,J,L,d)*hprime(L,K); \
        END DO; \
        dxx = dS_dxi*InvGrad(ee,0,0,i,j,k,bnum)+dS_deta*InvGrad(ee,0,1,i,j,k,bnum)+dS_dzeta*InvGrad(ee,0,2,i,j,k,bnum); \
        dxy = dS_dxi*InvGrad(ee,1,0,i,j,k,bnum)+dS_deta*InvGrad(ee,1,1,i,j,k,bnum)+dS_dzeta*InvGrad(ee,1,2,i,j,k,bnum); \
        dxz = dS_dxi*InvGrad(ee,2,0,i,j,k,bnum)+dS_deta*InvGrad(ee,2,1,i,j,k,bnum)+dS_dzeta*InvGrad(ee,2,2,i,j,k,bnum);

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

subroutine PROCNAME(dom, ngllvar, champs1)
    use champs_solidpml
    type(domain_solidpml), intent(inout) :: dom
    type(champssolidpml), intent(inout) :: champs1
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    !
    integer :: bnum, nblocks

    nblocks = dom%nblocks
    call PROC_MAIN(ngll, nblocks, dom%nglltot, dom%dt, dom%GLLw, dom%hprime, dom%m_Idom, &
            dom%m_InvGrad, dom%m_Jacob, champs1%VelocPML, dom%m_Lambda, dom%m_Mu, &
            dom%m_Stress1, dom%m_Stress2, dom%m_Stress3, champs1%ForcesPML, &
            dom%m_PMLDumpSx, dom%m_PMLDumpSy, dom%m_PMLDumpSz)
end subroutine PROCNAME

subroutine PROC_MAIN(ngllvar,nblocks,nglltot,dt,gllw,hprime,Idom, &
    InvGrad,Jacob,VelocPML,vLambda, vMu, &
    Stress1, Stress2, Stress3, ForcesPML,PMLDumpSx, PMLDumpSy, PMLDumpSz)
    implicit none
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), intent(in), dimension(0:ngll-1,0:ngll-1) :: hprime
    real(fpp), intent(in), dimension(0:ngll-1) :: gllw
    !real(fpp), intent(out), dimension (0:VCHUNK-1,0:5,0:ngll-1,0:ngll-1,0:ngll-1) :: Stress
    real(fpp), intent(inout), dimension (0:nglltot,0:2,0:2) :: VelocPML
    integer, intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: Idom
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: vLambda, vMu
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: InvGrad
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: Jacob
    real(fpp), intent(inout), dimension(0:VCHUNK-1,0:5,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: Stress1
    real(fpp), intent(inout), dimension(0:VCHUNK-1,0:5,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: Stress2
    real(fpp), intent(inout), dimension(0:VCHUNK-1,0:5,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: Stress3
    real(fpp), intent(inout), dimension(0:nglltot,0:2,0:2) :: ForcesPML
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: PMLDumpSx
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: PMLDumpSy
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks) :: PMLDumpSz
    real(fpp), intent(in) :: dt
    integer, intent(in) :: nblocks, nglltot
    !
    real(fpp), dimension (0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2) :: Veloc
    real(fpp) :: dVx_dx, dVx_dy, dVx_dz
    real(fpp) :: dVy_dx, dVy_dy, dVy_dz
    real(fpp) :: dVz_dx, dVz_dy, dVz_dz
    real(fpp) :: dS_dxi, dS_deta, dS_dzeta
    integer :: i, j, k, l, ind, i_dir, e, ee, bnum
    real(fpp) :: lambda, mu, lam2mu, dumpsx0, dumpsx1, dumpsy0, dumpsy1, dumpsz0, dumpsz1
    real(fpp) :: wi,wj,wk,acoeff
    real(fpp) :: ljxx,ljxy,ljxz,ljyx,ljyy,ljyz,ljzx,ljzy,ljzz, jac
    real(fpp) :: sxx,sxy,sxz,syy,syz,szz

    do bnum = 0,nblocks-1
        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        BEGIN_SUBELEM_LOOP1(ee)
                        ind = Idom(ee,i,j,k,bnum)
                        Veloc(ee,i,j,k,i_dir) = VelocPML(ind,i_dir,0) + VelocPML(ind,i_dir,1) + VelocPML(ind,i_dir,2)
                        END_SUBELEM_LOOP()
                    enddo
                enddo
            enddo
        enddo
    enddo

    do bnum = 0,nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP1(ee)
                    ! partial of velocity components with respect to xi,eta,zeta
                    x_deriv_ijke(Veloc,0,dS_dxi,dS_deta,dS_dzeta,dVx_dx,dVx_dy,dVx_dz)
                    x_deriv_ijke(Veloc,1,dS_dxi,dS_deta,dS_dzeta,dVy_dx,dVy_dy,dVy_dz)
                    x_deriv_ijke(Veloc,2,dS_dxi,dS_deta,dS_dzeta,dVz_dx,dVz_dy,dVz_dz)
                    ! shortcuts
                    lambda = vLambda(ee,i,j,k,bnum)
                    mu = vMu(ee,i,j,k,bnum)
                    lam2mu = lambda+2*mu
                    dumpsx0 = PMLDumpSx(ee,0,i,j,k,bnum)
                    dumpsx1 = PMLDumpSx(ee,1,i,j,k,bnum)
                    ! Update stresses
                    ! 0:xx 1:yy 2:zz 3:xy 4:xz 5:yz
                    Stress1(ee,0,i,j,k,bnum)=dumpsx0*Stress1(ee,0,i,j,k,bnum)+dumpsx1*Dt*(lam2mu)*dVx_dx
                    Stress1(ee,1,i,j,k,bnum)=dumpsx0*Stress1(ee,1,i,j,k,bnum)+dumpsx1*Dt*(lambda)*dVx_dx
                    Stress1(ee,2,i,j,k,bnum)=dumpsx0*Stress1(ee,2,i,j,k,bnum)+dumpsx1*Dt*(lambda)*dVx_dx
                    Stress1(ee,3,i,j,k,bnum)=dumpsx0*Stress1(ee,3,i,j,k,bnum)+dumpsx1*Dt*(mu)*dVy_dx
                    Stress1(ee,4,i,j,k,bnum)=dumpsx0*Stress1(ee,4,i,j,k,bnum)+dumpsx1*Dt*(mu)*dVz_dx
                    Stress1(ee,5,i,j,k,bnum)=dumpsx0*Stress1(ee,5,i,j,k,bnum)

                    dumpsy0 = PMLDumpSy(ee,0,i,j,k,bnum)
                    dumpsy1 = PMLDumpSy(ee,1,i,j,k,bnum)
                    Stress2(ee,0,i,j,k,bnum)=dumpsy0*Stress2(ee,0,i,j,k,bnum)+dumpsy1*Dt*(lambda)*dVy_dy
                    Stress2(ee,1,i,j,k,bnum)=dumpsy0*Stress2(ee,1,i,j,k,bnum)+dumpsy1*Dt*(lam2mu)*dVy_dy
                    Stress2(ee,2,i,j,k,bnum)=dumpsy0*Stress2(ee,2,i,j,k,bnum)+dumpsy1*Dt*(lambda)*dVy_dy
                    Stress2(ee,3,i,j,k,bnum)=dumpsy0*Stress2(ee,3,i,j,k,bnum)+dumpsy1*Dt*(mu)*dVx_dy
                    Stress2(ee,4,i,j,k,bnum)=dumpsy0*Stress2(ee,4,i,j,k,bnum)
                    Stress2(ee,5,i,j,k,bnum)=dumpsy0*Stress2(ee,5,i,j,k,bnum)+dumpsy1*Dt*(mu)*dVz_dy

                    dumpsz0 = PMLDumpSz(ee,0,i,j,k,bnum)
                    dumpsz1 = PMLDumpSz(ee,1,i,j,k,bnum)
                    Stress3(ee,0,i,j,k,bnum)=dumpsz0*Stress3(ee,0,i,j,k,bnum)+dumpsz1*Dt*(lambda)*dVz_dz
                    Stress3(ee,1,i,j,k,bnum)=dumpsz0*Stress3(ee,1,i,j,k,bnum)+dumpsz1*Dt*(lambda)*dVz_dz
                    Stress3(ee,2,i,j,k,bnum)=dumpsz0*Stress3(ee,2,i,j,k,bnum)+dumpsz1*Dt*(lam2mu)*dVz_dz
                    Stress3(ee,3,i,j,k,bnum)=dumpsz0*Stress3(ee,3,i,j,k,bnum)
                    Stress3(ee,4,i,j,k,bnum)=dumpsz0*Stress3(ee,4,i,j,k,bnum)+dumpsz1*Dt*(mu)*dVx_dz
                    Stress3(ee,5,i,j,k,bnum)=dumpsz0*Stress3(ee,5,i,j,k,bnum)+dumpsz1*Dt*(mu)*dVy_dz
                    !
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    enddo
    do bnum = 0,nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP1(ee)
                    sxx = Stress1(ee,0,i,j,k,bnum)+Stress2(ee,0,i,j,k,bnum)+Stress3(ee,0,i,j,k,bnum)
                    syy = Stress1(ee,1,i,j,k,bnum)+Stress2(ee,1,i,j,k,bnum)+Stress3(ee,1,i,j,k,bnum)
                    szz = Stress1(ee,2,i,j,k,bnum)+Stress2(ee,2,i,j,k,bnum)+Stress3(ee,2,i,j,k,bnum)
                    sxy = Stress1(ee,3,i,j,k,bnum)+Stress2(ee,3,i,j,k,bnum)+Stress3(ee,3,i,j,k,bnum)
                    sxz = Stress1(ee,4,i,j,k,bnum)+Stress2(ee,4,i,j,k,bnum)+Stress3(ee,4,i,j,k,bnum)
                    syz = Stress1(ee,5,i,j,k,bnum)+Stress2(ee,5,i,j,k,bnum)+Stress3(ee,5,i,j,k,bnum)
                    load_jacobian_matrix(i,j,k)
                    wi = Gllw(i)
                    wj = Gllw(j)
                    wk = Gllw(k)
                    jac = Jacob(ee,i,j,k,bnum)
                    do l = 0,ngll-1
                        ! Forces1
                        ind = Idom(ee,l,j,k,bnum)
                        acoeff = - hprime(l,i)*wi*wj*wk*jac
                        ForcesPML(ind,0,0) = ForcesPML(ind,0,0) + acoeff*(LJXX*sxx+LJYX*sxy+LJZX*sxz)
                        ForcesPML(ind,1,0) = ForcesPML(ind,1,0) + acoeff*(LJXX*sxy+LJYX*syy+LJZX*syz)
                        ForcesPML(ind,2,0) = ForcesPML(ind,2,0) + acoeff*(LJXX*sxz+LJYX*syz+LJZX*szz)
                        !!
                        ind = Idom(ee,i,l,k,bnum)
                        acoeff = - hprime(l,j)*wi*wj*wk*jac
                        ForcesPML(ind,0,1) = ForcesPML(ind,0,1) + acoeff*(LJXY*sxx+LJYY*sxy+LJZY*sxz)
                        ForcesPML(ind,1,1) = ForcesPML(ind,1,1) + acoeff*(LJXY*sxy+LJYY*syy+LJZY*syz)
                        ForcesPML(ind,2,1) = ForcesPML(ind,2,1) + acoeff*(LJXY*sxz+LJYY*syz+LJZY*szz)
                        !!
                        ind = Idom(ee,i,j,l,bnum)
                        acoeff = - hprime(l,k)*wi*wj*wk*jac
                        ForcesPML(ind,0,2) = ForcesPML(ind,0,2) + acoeff*(LJXZ*sxx+LJYZ*sxy+LJZZ*sxz)
                        ForcesPML(ind,1,2) = ForcesPML(ind,1,2) + acoeff*(LJXZ*sxy+LJYZ*syy+LJZZ*syz)
                        ForcesPML(ind,2,2) = ForcesPML(ind,2,2) + acoeff*(LJXZ*sxz+LJYZ*syz+LJZZ*szz)
                    end do
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    enddo
end subroutine PROC_MAIN

#if 0
subroutine PROC_MAIN(ngllvar,bnum,nblocks,nglltot,ForcesPML,Idom,Stress,InvGrad,Jacob,GLLw,hprime)
    use champs_solid
    implicit none

    integer, intent(in) :: bnum, nblocks, nglltot
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), intent(inout), dimension(0:nglltot,0:2,0:2) :: ForcesPML
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:5,0:ngll-1,0:ngll-1,0:ngll-1) :: Stress
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: InvGrad
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: Jacob
    real(fpp), intent(in), dimension(0:ngll-1) :: GLLw
    real(fpp), intent(in), dimension(0:ngll-1,0:ngll-1) :: hprime
    integer, intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: Idom
    !
    integer :: i, j, k, l, ind, e, ee
    real(fpp) :: sum_vx, sum_vy, sum_vz, acoeff
    real(fpp), dimension(0:VCHUNK-1,0:2,0:ngll-1,0:ngll-1,0:ngll-1)  :: Forces1, Forces2, Forces3


    Forces1 = 0d0
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i=0,ngll-1
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                sum_vx = 0d0
                sum_vy = 0d0
                sum_vz = 0d0
                do l = 0,ngll-1
                    !ind = Idom(ee,l,j,k,bnum)
                    !ForcesPML(ind,:,0) = ForcesPML(ind,:,0) + Forces1(ee,:,l,j,k)

                    acoeff = - hprime(l,i)*GLLw(i)*GLLw(j)*GLLw(k)*Jacob(ee,l,j,k,bnum)
                    ForcesPML(ind,0,0) = ForcesPML(ind,0,0) + acoeff*InvGrad(ee,0,0,i,j,k,bnum)*sxx
                    ForcesPML(ind,0,0) = ForcesPML(ind,0,0) + acoeff*InvGrad(ee,1,0,i,j,k,bnum)*sxy
                    ForcesPML(ind,0,0) = ForcesPML(ind,0,0) + acoeff*InvGrad(ee,2,0,i,j,k,bnum)*sxz

                    ForcesPML(ind,1,0) = ForcesPML(ind,1,0) + acoeff*InvGrad(ee,0,0,i,j,k,bnum)*sxy
                    ForcesPML(ind,1,0) = ForcesPML(ind,1,0) + acoeff*InvGrad(ee,1,0,i,j,k,bnum)*syy
                    ForcesPML(ind,1,0) = ForcesPML(ind,1,0) + acoeff*InvGrad(ee,2,0,i,j,k,bnum)*syz

                    sum_vz = sum_vz + acoeff*InvGrad(ee,0,0,l,j,k,bnum)*Stress(ee,4,l,j,k) ! zx
                    sum_vz = sum_vz + acoeff*InvGrad(ee,1,0,l,j,k,bnum)*Stress(ee,5,l,j,k) ! zy
                    sum_vz = sum_vz + acoeff*InvGrad(ee,2,0,l,j,k,bnum)*Stress(ee,2,l,j,k) ! zz
                end do
                Forces1(ee,0,i,j,k) = sum_vx
                Forces1(ee,1,i,j,k) = sum_vy
                Forces1(ee,2,i,j,k) = sum_vz
                END_SUBELEM_LOOP()
            end do
        end do
    end do

    Forces2 = 0d0
    do k = 0,ngll-1
        do l = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
#if VCHUNK>1
                    !$omp simd linear(e,ee) private(sum_vx,sum_vy,sum_vz,acoeff)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    acoeff = - hprime(j,l)*GLLw(i)*GLLw(l)*GLLw(k)*Jacob(ee,i,l,k,bnum)
                    sum_vx = acoeff*(InvGrad(ee,0,1,i,l,k,bnum)*Stress(ee,0,i,l,k) + &
                                     InvGrad(ee,1,1,i,l,k,bnum)*Stress(ee,3,i,l,k) + &
                                     InvGrad(ee,2,1,i,l,k,bnum)*Stress(ee,4,i,l,k))

                    sum_vy = acoeff*(InvGrad(ee,0,1,i,l,k,bnum)*Stress(ee,3,i,l,k) + &
                                     InvGrad(ee,1,1,i,l,k,bnum)*Stress(ee,1,i,l,k) + &
                                     InvGrad(ee,2,1,i,l,k,bnum)*Stress(ee,5,i,l,k))

                    sum_vz = acoeff*(InvGrad(ee,0,1,i,l,k,bnum)*Stress(ee,4,i,l,k) + &
                                     InvGrad(ee,1,1,i,l,k,bnum)*Stress(ee,5,i,l,k) + &
                                     InvGrad(ee,2,1,i,l,k,bnum)*Stress(ee,2,i,l,k))
                    Forces2(ee,0,i,j,k) = Forces2(ee,0,i,j,k) + sum_vx
                    Forces2(ee,1,i,j,k) = Forces2(ee,1,i,j,k) + sum_vy
                    Forces2(ee,2,i,j,k) = Forces2(ee,2,i,j,k) + sum_vz
                    END_SUBELEM_LOOP()
                end do
            end do
        end do
    end do

    ! TODO reorder loops ?
    Forces3 = 0
    do l = 0,ngll-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
#if VCHUNK>1
                    !$omp simd linear(e,ee) private(sum_vx,sum_vy,sum_vz,acoeff)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    acoeff = - hprime(k,l)*GLLw(i)*GLLw(j)*GLLw(l)*Jacob(ee,i,j,l,bnum)
                    sum_vx = acoeff*(InvGrad(ee,0,2,i,j,l,bnum)*Stress(ee,0,i,j,l) + &
                                     InvGrad(ee,1,2,i,j,l,bnum)*Stress(ee,3,i,j,l) + &
                                     InvGrad(ee,2,2,i,j,l,bnum)*Stress(ee,4,i,j,l))

                    sum_vy = acoeff*(InvGrad(ee,0,2,i,j,l,bnum)*Stress(ee,3,i,j,l) + &
                                     InvGrad(ee,1,2,i,j,l,bnum)*Stress(ee,1,i,j,l) + &
                                     InvGrad(ee,2,2,i,j,l,bnum)*Stress(ee,5,i,j,l))

                    sum_vz = acoeff*(InvGrad(ee,0,2,i,j,l,bnum)*Stress(ee,4,i,j,l) + &
                                     InvGrad(ee,1,2,i,j,l,bnum)*Stress(ee,5,i,j,l) + &
                                     InvGrad(ee,2,2,i,j,l,bnum)*Stress(ee,2,i,j,l))
                    Forces3(ee,0,i,j,k) = Forces3(ee,0,i,j,k) + sum_vx
                    Forces3(ee,1,i,j,k) = Forces3(ee,1,i,j,k) + sum_vy
                    Forces3(ee,2,i,j,k) = Forces3(ee,2,i,j,k) + sum_vz
                    END_SUBELEM_LOOP()
                end do
            end do
        end do
    end do

    ! Assemblage
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ind = Idom(ee,i,j,k,bnum)
                ForcesPML(ind,:,0) = ForcesPML(ind,:,0) + Forces1(ee,:,i,j,k)
                ForcesPML(ind,:,1) = ForcesPML(ind,:,1) + Forces2(ee,:,i,j,k)
                ForcesPML(ind,:,2) = ForcesPML(ind,:,2) + Forces3(ee,:,i,j,k)
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo
end subroutine PROC_MAIN

#endif

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
