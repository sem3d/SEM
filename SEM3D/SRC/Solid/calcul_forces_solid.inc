!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

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
    real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real(fpp) :: sxx,sxy,sxz,syy,syz,szz,t4
    real(fpp) :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
    !!real(fpp), parameter :: zero = 0.
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
    real(fpp) :: xmu
#ifdef NONLINEAR
    real(fpp), dimension(0:5) :: start,trial,center,dEps_alpha,dEpsP
    real(fpp)  :: syld,ckin,kkin,biso,rinf,radius,alpha_epl
    logical    :: st_epl
#endif
#ifdef ATTENUATION
#ifdef ANISO
    real(fpp) :: Eyz, Exz, Exy
#endif
    real(fpp) :: stt, x2mu, xpression, xkappa
    integer   :: i_sls
#else
    real(fpp) :: xla2mu, xla
#endif
#ifdef NGLLVAL
    if (ngll/=ngll0) stop 1
#endif
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
#ifdef ATTENUATION
                ! Save for reuse
                dUx_dxi  (ee,i,j,k) = dxx
                dUx_deta (ee,i,j,k) = dxy
                dUx_dzeta(ee,i,j,k) = dxz
                dUy_dxi  (ee,i,j,k) = dyx
                dUy_deta (ee,i,j,k) = dyy
                dUy_dzeta(ee,i,j,k) = dyz
                dUz_dxi  (ee,i,j,k) = dzx
                dUz_deta (ee,i,j,k) = dzy
                dUz_dzeta(ee,i,j,k) = dzz
#endif

#ifdef ATTENUATION
#ifdef ANISO
                Eyz = (DYZ+DZY)*M_SQRT1_2
                Exz = (DXZ+DZX)*M_SQRT1_2
                Exy = (DXY+DYX)*M_SQRT1_2
#define CC(d)  dom%Cij_(d,i,j,k,bnum,ee)
                sxx = DXX*CC( 0) + DYY*CC( 1) + DZZ*CC( 2) + EYZ*CC( 3) + EXZ*CC( 4) + EXY*CC( 5)
                syy = DXX*CC( 1) + DYY*CC( 6) + DZZ*CC( 7) + EYZ*CC( 8) + EXZ*CC( 9) + EXY*CC(10)
                szz = DXX*CC( 2) + DYY*CC( 7) + DZZ*CC(11) + EYZ*CC(12) + EXZ*CC(13) + EXY*CC(14)
                syz = DXX*CC( 3) + DYY*CC( 8) + DZZ*CC(12) + EYZ*CC(15) + EXZ*CC(16) + EXY*CC(17)
                sxz = DXX*CC( 4) + DYY*CC( 9) + DZZ*CC(13) + EYZ*CC(16) + EXZ*CC(18) + EXY*CC(19)
                sxy = DXX*CC( 5) + DYY*CC(10) + DZZ*CC(14) + EYZ*CC(17) + EXZ*CC(19) + EXY*CC(20)
#undef CC
#else
                xmu = dom%Mu_(i,j,k,bnum,ee) * dom%onemSbeta_(i,j,k,bnum,ee)
                x2mu = 2. * xmu
                sxx = x2mu * DXX
                sxy = xmu * ( DXY + DYX )
                sxz = xmu * ( DXZ + DZX )
                syy = x2mu * DYY
                syz = xmu * ( DYZ + DZY )
                szz = x2mu * DZZ
#endif
                xkappa = dom%Kappa_(i,j,k,bnum,ee)* dom%onemPbeta_(i,j,k,bnum,ee)
                xpression = xkappa*(DXX+DYY+DZZ)
                do i_sls = 0,dom%n_sls-1
                    xpression = xpression - dom%R_vol_(i_sls,i,j,k,bnum,ee)
                    sxx = sxx - dom%R_xx_(i_sls,i,j,k,bnum,ee)
                    syy = syy - dom%R_yy_(i_sls,i,j,k,bnum,ee)
                    ! ici on utilise le fait que la trace est nulle
                    szz = szz + dom%R_xx_(i_sls,i,j,k,bnum,ee) + dom%R_yy_(i_sls,i,j,k,bnum,ee)
                    sxy = sxy - dom%R_xy_(i_sls,i,j,k,bnum,ee)
                    sxz = sxz - dom%R_xz_(i_sls,i,j,k,bnum,ee)
                    syz = syz - dom%R_yz_(i_sls,i,j,k,bnum,ee)
                enddo
                stt = (sxx + syy + szz)*M_1_3
                sxx = sxx - stt + xpression
                syy = syy - stt + xpression
                szz = szz - stt + xpression
#else
                xmu = dom%Mu_    (i,j,k,bnum,ee)
                xla = dom%Lambda_(i,j,k,bnum,ee)
                xla2mu = xla + 2. * xmu
#ifdef NONLINEAR
                
                start       =   dom%stress_(0:5,i,j,k,bnum,ee)
                center      =   dom%center_(0:5,i,j,k,bnum,ee)
                radius      =   zero !dom%radius_(i,j,k,bnum,ee)
                syld        =   dom%nl_param%LMC%syld_(i,j,k,bnum,ee)
                kkin        =   dom%nl_param%LMC%kkin_(i,j,k,bnum,ee) 
                rinf        =   zero !dom%nl_param%LMC%rinf_(i,j,k,bnum,ee) 
                biso        =   zero !dom%nl_param%LMC%biso_(i,j,k,bnum,ee) 
                ckin        =   dom%nl_param%LMC%ckin_(i,j,k,bnum,ee)
                
                dEps_alpha(0:5) = (/dxx,dyy,dzz,(dxy+dyx),(dxz+dzx),(dyz+dzy)/) 
                dEps_alpha(0:5) = dEps_alpha(0:5)*dom%dt
                dom%strain_(0:5,i,j,k,bnum,ee) = dom%strain_(0:5,i,j,k,bnum,ee) + dEps_alpha(0:5) 
                dEpsP(0:5)      = dom%plstrain_(0:5,i,j,k,bnum,ee) 
                
                sxx = xla2mu*dEps_alpha(0)+xla*sum(dEps_alpha(1:2))
                syy = xla2mu*dEps_alpha(1)+xla*(dEps_alpha(0)+dEps_alpha(2))
                szz = xla2mu*dEps_alpha(2)+xla*sum(dEps_alpha(0:1))
                sxy = xmu*dEps_alpha(3)
                sxz = xmu*dEps_alpha(4)
                syz = xmu*dEps_alpha(5)
                trial       =   (/sxx,syy,szz,sxy,sxz,syz/)
                call check_plasticity(trial, start, center,radius, &
                    syld,st_epl,alpha_epl)

                if (st_epl) then

                    dEps_alpha=(1.0d0-alpha_epl)*dEps_alpha
                    call plastic_corrector(dEps_alpha, trial, center, syld, &
                            radius, biso, rinf, ckin, kkin, xmu, xla, dEpsP)
                end if
                sxx = trial(0)
                syy = trial(1)
                szz = trial(2)
                sxy = trial(3)
                sxz = trial(4)
                syz = trial(5)
                
                dom%stress_(0:5,i,j,k,bnum,ee)   = trial(0:5)
                dom%center_(0:5,i,j,k,bnum,ee)   = center(0:5)
                dom%plstrain_(0:5,i,j,k,bnum,ee) = dEpsP(0:5)
                dom%radius_(i,j,k,bnum,ee)       = zero !radius

#else
                sxx = xla2mu * DXX + xla * ( DYY + DZZ )
                sxy = xmu * ( DXY + DYX )
                sxz = xmu * ( DXZ + DZX )
                syy = xla2mu * DYY + xla * ( DXX + DZZ )
                syz = xmu * ( DYZ + DZY )
                szz = xla2mu * DZZ + xla * ( DXX + DYY )
#endif
#endif
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
                xt5 = syy*JYX+sxy*JXX+syz*JZX
                xt6 = syy*JYY+sxy*JXY+syz*JZY
                xt7 = syy*JYZ+sxy*JXZ+syz*JZZ
                !
                !=====================
                !       FZ
                !=====================
                !
                xt8 = szz*JZX+sxz*JXX+syz*JYX
                xt9 = szz*JZY+sxz*JXY+syz*JYY
                xt10= szz*JZZ+sxz*JXZ+syz*JYZ

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
#ifdef ATTENUATION
    call PROCNAME_ATN(dom,bnum,dUx_dxi,dUx_deta,dUx_dzeta,&
        dUy_dxi,dUy_deta,dUy_dzeta, &
        dUz_dxi,dUz_deta,dUz_dzeta)
#endif
end subroutine PROCNAME

#ifdef ATTENUATION
subroutine PROCNAME_ATN(dom,bnum,mdxx,mdxy,mdxz,mdyx,mdyy,mdyz,mdzx,mdzy,mdzz)
    use champs_solid
    use deriv3d
    implicit none

    type(domain_solid), intent (INOUT) :: dom
    integer, intent(in) :: bnum
    real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: mDXX,mDXY,mDXZ
    real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: mDYX,mDYY,mDYZ
    real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: mDZX,mDZY,mDZZ
#define DXX mDXX(ee,i,j,k)
#define DXY mDXY(ee,i,j,k)
#define DXZ mDXZ(ee,i,j,k)
#define DYX mDYX(ee,i,j,k)
#define DYY mDYY(ee,i,j,k)
#define DYZ mDYZ(ee,i,j,k)
#define DZX mDZX(ee,i,j,k)
#define DZY mDZY(ee,i,j,k)
#define DZZ mDZZ(ee,i,j,k)
    !
    integer :: n_solid

    real(fpp) :: epsilondev_xx_loc,epsilondev_yy_loc
    real(fpp) :: epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc
    real(fpp) :: epsilonvol_loc
    integer   :: i_sls,i,j,k, ngll,e,ee
    real(fpp) :: factorS_loc,Sn,Snp1
    real(fpp) :: factorP_loc,Pn,Pnp1
    real(fpp) :: epsilon_trace_over_3
    real(fpp) :: factor_common_3, factor_common_P, omega_tau_s, gamma_mu, gamma_kappa
    real(fpp) :: alphaval, betaval, gammaval
    real(fpp) :: dt_tau
    n_solid = dom%n_sls

    ngll = dom%ngll
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(e,ee)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                epsilon_trace_over_3 = 0.333333333333333333333333333333d0 * (DXX + DYY + DZZ)
                epsilonvol_loc = DXX + DYY + DZZ
                epsilondev_xx_loc = DXX - epsilon_trace_over_3
                epsilondev_yy_loc = DYY - epsilon_trace_over_3
                epsilondev_xy_loc = 0.5 * (DXY + DYX)
                epsilondev_xz_loc = 0.5 * (DZX + DXZ)
                epsilondev_yz_loc = 0.5 * (DZY + DYZ)

                do i_sls = 0,n_solid-1
                    omega_tau_s = dom%omega_tau_s(i_sls)
                    gamma_kappa = dom%agamma_kappa_(i_sls,i,j,k,bnum,ee)
                    gamma_mu = dom%agamma_mu_(i_sls,i,j,k,bnum,ee)
#ifndef ANISO
                    ! get coefficients for that standard linear solid
                    factor_common_P = 2d0*gamma_kappa*omega_tau_s
                    factorP_loc   = dom%Kappa_(i,j,k,bnum,ee) * factor_common_P
                    RK4_attenu_coefs(dom%dt,omega_tau_s,alphaval, betaval, gammaval)

                    !  terme volumique
                    Pn   = factorP_loc * dom%epsilonvol_(i,j,k,bnum,ee)
                    Pnp1 = factorP_loc * epsilonvol_loc
                    dom%R_vol_(i_sls,i,j,k,bnum,ee) = alphaval * dom%R_vol_(i_sls,i,j,k,bnum,ee)  &
                        + betaval*Pn+gammaval*Pnp1
#endif
                    ! get coefficients for that standard linear solid
                    factor_common_3 = 2d0*gamma_mu*omega_tau_s

                    factorS_loc   = dom%Mu_(i,j,k,bnum,ee) * factor_common_3

                    !     termes deviatoires
                    ! term in xx
                    Sn    = factorS_loc * dom%epsilondev_xx_(i,j,k,bnum,ee)
                    Snp1  = factorS_loc * epsilondev_xx_loc
                    dom%R_xx_(i_sls,i,j,k,bnum,ee) = alphaval * dom%R_xx_(i_sls,i,j,k,bnum,ee)  &
                        + betaval*Sn+gammaval*Snp1

                    ! term in yy
                    Sn    = factorS_loc  * dom%epsilondev_yy_(i,j,k,bnum,ee)
                    Snp1  = factorS_loc  * epsilondev_yy_loc
                    dom%R_yy_(i_sls,i,j,k,bnum,ee) = alphaval * dom%R_yy_(i_sls,i,j,k,bnum,ee)  &
                        + betaval*Sn+gammaval*Snp1

                    ! term in xy
                    Sn    = factorS_loc  * dom%epsilondev_xy_(i,j,k,bnum,ee)
                    Snp1  = factorS_loc  * epsilondev_xy_loc
                    dom%R_xy_(i_sls,i,j,k,bnum,ee) = alphaval * dom%R_xy_(i_sls,i,j,k,bnum,ee) &
                        + betaval*Sn+gammaval*Snp1

                    ! term in xz
                    Sn    = factorS_loc * dom%epsilondev_xz_(i,j,k,bnum,ee)
                    Snp1  = factorS_loc * epsilondev_xz_loc
                    dom%R_xz_(i_sls,i,j,k,bnum,ee) = alphaval * dom%R_xz_(i_sls,i,j,k,bnum,ee) &
                        + betaval*Sn+gammaval*Snp1

                    ! term in yz
                    Sn    = factorS_loc * dom%epsilondev_yz_(i,j,k,bnum,ee)
                    Snp1  = factorS_loc * epsilondev_yz_loc
                    dom%R_yz_(i_sls,i,j,k,bnum,ee) = alphaval * dom%R_yz_(i_sls,i,j,k,bnum,ee)  &
                        + betaval*Sn+gammaval*Snp1

                enddo   ! end of loop on memory variables

                ! a la fin on stocke la deformation deviatorique pour le prochain
                ! pas de temps de la mise a jour de Runge-Kutta faite ci-dessus

                ! save deviatoric strain for Runge-Kutta scheme

                !     partie deviatoire
                dom%epsilondev_xx_(i,j,k,bnum,ee) = epsilondev_xx_loc
                dom%epsilondev_yy_(i,j,k,bnum,ee) = epsilondev_yy_loc
                dom%epsilondev_xy_(i,j,k,bnum,ee) = epsilondev_xy_loc
                dom%epsilondev_xz_(i,j,k,bnum,ee) = epsilondev_xz_loc
                dom%epsilondev_yz_(i,j,k,bnum,ee) = epsilondev_yz_loc
#ifndef ANISO
                !   partie isotrope
                dom%epsilonvol_   (i,j,k,bnum,ee) = epsilonvol_loc
#endif
                END_SUBELEM_LOOP()
            end do
        end do
    end do
end subroutine PROCNAME_ATN
#undef DXX
#undef DXY
#undef DXZ
#undef DYX
#undef DYY
#undef DYZ
#undef DZX
#undef DZY
#undef DZZ
#endif
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

