!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#ifndef PROCNAME
#error "Must define PROCNAME before inclusion"
#endif

#ifdef NGLL_GEN
#define ngllvar ngll
#else
#define ngllvar ngll0
#endif

#define ORIG

#define part_deriv_ijke(d,dS_dxi,dS_deta,dS_dzeta,dxx,dxy,dxz) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO L = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Depla(ee,L,J,K,d)*dom%hprime(L,I); \
            dS_deta  = dS_deta +Depla(ee,I,L,K,d)*dom%hprime(L,J); \
            dS_dzeta = dS_dzeta+Depla(ee,I,J,L,d)*dom%hprime(L,K); \
        END DO; \
        dxx = dS_dxi*InvGrad(IND_MNE(0,0,ee))+dS_deta*InvGrad(IND_MNE(0,1,ee))+dS_dzeta*InvGrad(IND_MNE(0,2,ee)); \
        dxy = dS_dxi*InvGrad(IND_MNE(1,0,ee))+dS_deta*InvGrad(IND_MNE(1,1,ee))+dS_dzeta*InvGrad(IND_MNE(1,2,ee)); \
        dxz = dS_dxi*InvGrad(IND_MNE(2,0,ee))+dS_deta*InvGrad(IND_MNE(2,1,ee))+dS_dzeta*InvGrad(IND_MNE(2,2,ee));

subroutine PROCNAME(dom,ngllvar,lnum,Fox,Foy,Foz,Depla)
    use champs_solid
    implicit none

    type(domain_solid), intent (INOUT) :: dom
    integer, intent(in) :: lnum
    integer, intent(in) :: ngllvar
#ifndef NGLL_GEN
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: Fox,Foz,Foy
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(in) :: Depla
    integer :: i,j,k,l,e,ee
    real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real(fpp) :: sxx,sxy,sxz,syy,syz,szz,t4
    real(fpp) :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
    real(fpp), parameter :: zero = 0.
#ifdef ORIG
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t1,t5,t8
    ! Les indices sont reordonnes, probablement pour la localite memoire
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t2,t6,t9
    real(fpp), dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: t3,t7,t10
    real(fpp) :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
    real(fpp) :: F1,F2,F3
#endif
    real(fpp) :: xla, xmu,xla2mu
    real(fpp), dimension(IND_MNE(0:2,0:2,0:CHUNK-1)) :: invgrad
    real(fpp) :: dS_dxi, dS_deta, dS_dzeta
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                InvGrad(IND_MNE(0:2,0:2,0:CHUNK-1)) = dom%InvGrad_(:,:,i,j,k,lnum:lnum+CHUNK-1)
#ifdef SEM_VEC
!$omp simd linear(e,ee) private(xla,xmu,xla2mu,sxx,sxy,sxz,syy,syz,szz,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10,dS_dxi, dS_deta, dS_dzeta)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,lnum)

                ! Calcul de dPhi/dx dPhi/dy dPhi/dz
!                call physical_part_deriv_ijke(ee,i,j,k,ngll,dom%hprime,InvGrad,Depla(:,:,:,:,0),dxx,dxy,dxz)
!                call physical_part_deriv_ijke(ee,i,j,k,ngll,dom%hprime,InvGrad,Depla(:,:,:,:,1),dyx,dyy,dyz)
!                call physical_part_deriv_ijke(ee,i,j,k,ngll,dom%hprime,InvGrad,Depla(:,:,:,:,2),dzx,dzy,dzz)
                part_deriv_ijke(0,dS_dxi,dS_deta,dS_dzeta,dxx,dxy,dxz)
                part_deriv_ijke(1,dS_dxi,dS_deta,dS_dzeta,dyx,dyy,dyz)
                part_deriv_ijke(2,dS_dxi,dS_deta,dS_dzeta,dzx,dzy,dzz)

                xmu = dom%Mu_    (i,j,k,e)
                xla = dom%Lambda_(i,j,k,e)
                xla2mu = xla + 2. * xmu

                sxx = xla2mu * DXX + xla * ( DYY + DZZ )
                sxy = xmu * ( DXY + DYX )
                sxz = xmu * ( DXZ + DZX )
                syy = xla2mu * DYY + xla * ( DXX + DZZ )
                syz = xmu * ( DYZ + DZY )
                szz = xla2mu * DZZ + xla * ( DXX + DYY )

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
                t4 = dom%Jacob_(i,j,k,lnum) * dom%GLLw(i)
                xt1  =  xt1 * t4
                xt5  =  xt5 * t4
                xt8  =  xt8 * t4

                t4 = dom%Jacob_(i,j,k,lnum) * dom%GLLw(j)
                xt2  =  xt2 * t4
                xt6  =  xt6 * t4
                xt9  =  xt9 * t4

                t4 = dom%Jacob_(i,j,k,lnum) * dom%GLLw(k)
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
