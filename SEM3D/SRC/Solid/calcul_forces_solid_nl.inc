!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#ifndef PROCNAME_NL
#error "Must define PROCNAME_NL before inclusion"
#endif

#undef ngllvar
#ifndef NGLLVAL
#define ngllvar ngll
#else
#define ngllvar ngll0
#endif

#define ORIG

subroutine PROCNAME_NL(dom,ngllvar,lnum,Fox,Foy,Foz,Depla)
    use champs_solid
    implicit none

    type(domain_solid), intent (INOUT) :: dom
    integer, intent(in) :: lnum
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: Fox,Foz,Foy
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(in) :: Depla
    integer :: i,j,k,l,e,ee,st_epl
    real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real(fpp) :: sxx,sxy,sxz,syy,syz,szz,t4
    real(fpp) :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
    real(fpp), parameter :: zero = 0.
    real(fpp), dimension(0:5) :: start,trial,center,dEps_alpha,dEpsP,eps0
    real(fpp)  :: syld,ckin,kkin,biso,rinf,radius,alpha_elp
#ifdef ORIG
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t1,t5,t8
    ! Les indices sont reordonnes, probablement pour la localite memoire
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t2,t6,t9
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t3,t7,t10
    real(fpp) :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
    real(fpp) :: F1,F2,F3
#endif
#ifdef ANISO
    real(fpp) :: Eyz, Exz, Exy
#endif
    real(fpp) :: xla, xmu
#ifdef ATTENUATION
    real(fpp) :: stt, x2mu, xpression, xkappa
    integer   :: i_sls
#else
    real(fpp) :: xla2mu
#endif
    real(fpp), dimension(IND_MNE(0:2,0:2,0:CHUNK-1)) :: invgrad
    real(fpp) :: dS_dxi, dS_deta, dS_dzeta
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                InvGrad(IND_MNE(0:2,0:2,0:CHUNK-1)) = dom%InvGrad_(:,:,i,j,k,lnum:lnum+CHUNK-1)
#ifdef SEM_VEC
!$omp simd linear(e,ee)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,lnum)

                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                part_deriv_ijke(Depla,0,dS_dxi,dS_deta,dS_dzeta,dxx,dxy,dxz)
                part_deriv_ijke(Depla,1,dS_dxi,dS_deta,dS_dzeta,dyx,dyy,dyz)
                part_deriv_ijke(Depla,2,dS_dxi,dS_deta,dS_dzeta,dzx,dzy,dzz)

#ifdef ATTENUATION
#ifdef ANISO
                Eyz = (DYZ+DZY)*M_SQRT1_2
                Exz = (DXZ+DZX)*M_SQRT1_2
                Exy = (DXY+DYX)*M_SQRT1_2
#define CC(d)  dom%Cij_(d,i,j,k,e)
                sxx = DXX*CC( 0) + DYY*CC( 1) + DZZ*CC( 2) + EYZ*CC( 3) + EXZ*CC( 4) + EXY*CC( 5)
                syy = DXX*CC( 1) + DYY*CC( 6) + DZZ*CC( 7) + EYZ*CC( 8) + EXZ*CC( 9) + EXY*CC(10)
                szz = DXX*CC( 2) + DYY*CC( 7) + DZZ*CC(11) + EYZ*CC(12) + EXZ*CC(13) + EXY*CC(14)
                syz = DXX*CC( 3) + DYY*CC( 8) + DZZ*CC(12) + EYZ*CC(15) + EXZ*CC(16) + EXY*CC(17)
                sxz = DXX*CC( 4) + DYY*CC( 9) + DZZ*CC(13) + EYZ*CC(16) + EXZ*CC(18) + EXY*CC(19)
                sxy = DXX*CC( 5) + DYY*CC(10) + DZZ*CC(14) + EYZ*CC(17) + EXZ*CC(19) + EXY*CC(20)
#undef CC
#else
                xmu = dom%Mu_(i,j,k,e) * dom%onemSbeta_(i,j,k,e)
                x2mu = 2. * xmu
                sxx = x2mu * DXX
                sxy = xmu * ( DXY + DYX )
                sxz = xmu * ( DXZ + DZX )
                syy = x2mu * DYY
                syz = xmu * ( DYZ + DZY )
                szz = x2mu * DZZ
#endif
                xkappa = dom%Kappa_(i,j,k,e)* dom%onemPbeta_(i,j,k,e)
                xpression = xkappa*(DXX+DYY+DZZ)
                do i_sls = 0,dom%n_sls-1
                    xpression = xpression - dom%R_vol_(i_sls,i,j,k,e)
                    sxx = sxx - dom%R_xx_(i_sls,i,j,k,e)
                    syy = syy - dom%R_yy_(i_sls,i,j,k,e)
                    ! ici on utilise le fait que la trace est nulle
                    szz = szz + dom%R_xx_(i_sls,i,j,k,e) + dom%R_yy_(i_sls,i,j,k,e)
                    sxy = sxy - dom%R_xy_(i_sls,i,j,k,e)
                    sxz = sxz - dom%R_xz_(i_sls,i,j,k,e)
                    syz = syz - dom%R_yz_(i_sls,i,j,k,e)
                enddo
                stt = (sxx + syy + szz)*M_1_3
                sxx = sxx - stt + xpression
                syy = syy - stt + xpression
                szz = szz - stt + xpression
#else
                xmu = dom%Mu_    (i,j,k,e)
                xla = dom%Lambda_(i,j,k,e)
                xla2mu = xla + 2. * xmu
 
                eps0(0:5) = dom%strain_(0:5,i,j,k,e)
                dom%strain_(0:5,i,j,k,e) = (/dxx,dyy,dzz,(dxy+dyx),(dxz+dzx),(dyz+dzy)/)
                dEps_alpha(0:5) = dom%strain_(0:5,i,j,k,e) -eps0(0:5)
                dEpsP(0:5)      = 0.d0
                
                sxx = xla2mu*dEps_alpha(0)+xla*sum(dEps_alpha(1:2))
                syy = xla2mu*dEps_alpha(1)+xla*(dEps_alpha(0)+dEps_alpha(2))
                szz = xla2mu*dEps_alpha(2)+xla*sum(dEps_alpha(0:1))
                sxy = xmu*dEps_alpha(3)
                sxz = xmu*dEps_alpha(4)
                syz = xmu*dEps_alpha(5)
#endif
                start       =   dom%stress_(0:5,i,j,k,e)
                trial       =   (/sxx,syy,szz,sxy,sxz,syz/)
                center      =   dom%center_(0:5,i,j,k,e)
                radius      =   dom%radius_(i,j,k,e)
                syld        =   dom%nl_param%LMC%syld_(i,j,k,e)
                ckin        =   dom%nl_param%LMC%Ckin_(i,j,k,e)
                kkin        =   dom%nl_param%LMC%kkin_(i,j,k,e) 
                rinf        =   dom%nl_param%LMC%rinf_(i,j,k,e) 
                biso        =   dom%nl_param%LMC%biso_(i,j,k,e) 
                call check_plasticity(trial, start, center,radius, &
                    syld,st_epl,alpha_elp)

                if (st_epl == 1) then

                    dEps_alpha=(1-alpha_elp)*dEps_alpha
                    call plastic_corrector(dEps_alpha, trial, center, syld, &
                            radius, biso, rinf, ckin, kkin, xmu, xla, dEpsP)
                end if
                sxx = trial(0)
                syy = trial(1)
                szz = trial(2)
                sxy = trial(3)
                sxz = trial(4)
                syz = trial(5)
                
                dom%stress(0:5,i,j,k,e) = trial(0:5)
                dom%center(0:5,i,j,k,e) = center(0:5)
                dom%plstrain(0:5,i,j,k,e) = dom%plstrain(0:5,i,j,k,e)+dEpsP(0:5)
                dom%radius(i,j,k,e)     = radius

                !
                !=====================
                !       FX
                !=====================
                !
                xt1 = sxx * InvGrad(IND_MNE(0,0,ee)) + sxy * InvGrad(IND_MNE(1,0,ee)) + sxz * InvGrad(IND_MNE(2,0,ee))
                xt2 = sxx * InvGrad(IND_MNE(0,1,ee)) + sxy * InvGrad(IND_MNE(1,1,ee)) + sxz * InvGrad(IND_MNE(2,1,ee))
                xt3 = sxx * InvGrad(IND_MNE(0,2,ee)) + sxy * InvGrad(IND_MNE(1,2,ee)) + sxz * InvGrad(IND_MNE(2,2,ee))
                !
                !=====================
                !       FY
                !=====================
                !
                xt5 = syy * InvGrad(IND_MNE(1,0,ee)) + sxy * InvGrad(IND_MNE(0,0,ee)) + syz * InvGrad(IND_MNE(2,0,ee))
                xt6 = syy * InvGrad(IND_MNE(1,1,ee)) + sxy * InvGrad(IND_MNE(0,1,ee)) + syz * InvGrad(IND_MNE(2,1,ee))
                xt7 = syy * InvGrad(IND_MNE(1,2,ee)) + sxy * InvGrad(IND_MNE(0,2,ee)) + syz * InvGrad(IND_MNE(2,2,ee))
                !
                !=====================
                !       FZ
                !=====================
                !
                xt8  = szz * InvGrad(IND_MNE(2,0,ee)) + sxz * InvGrad(IND_MNE(0,0,ee)) + syz * InvGrad(IND_MNE(1,0,ee))
                xt9  = szz * InvGrad(IND_MNE(2,1,ee)) + sxz * InvGrad(IND_MNE(0,1,ee)) + syz * InvGrad(IND_MNE(1,1,ee))
                xt10 = szz * InvGrad(IND_MNE(2,2,ee)) + sxz * InvGrad(IND_MNE(0,2,ee)) + syz * InvGrad(IND_MNE(1,2,ee))

                !
                !- Multiplication par le Jacobien et le poids d'integration
                !
                t4 = dom%Jacob_(i,j,k,e) * dom%GLLw(i)
                xt1  =  xt1 * t4
                xt5  =  xt5 * t4
                xt8  =  xt8 * t4

                t4 = dom%Jacob_(i,j,k,e) * dom%GLLw(j)
                xt2  =  xt2 * t4
                xt6  =  xt6 * t4
                xt9  =  xt9 * t4

                t4 = dom%Jacob_(i,j,k,e) * dom%GLLw(k)
                xt3  =  xt3 * t4
                xt7  =  xt7 * t4
                xt10 = xt10 * t4
#ifdef ORIG
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
#ifdef ORIG
    !=-=-=-=-=-=-=-=-=-=-
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#ifdef SEM_VEC
                !dir$ simd
                !dir$ vector aligned
#endif
                BEGIN_SUBELEM_LOOP(e,ee,lnum)
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
end subroutine PROCNAME_NL

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
