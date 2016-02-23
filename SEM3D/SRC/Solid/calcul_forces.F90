!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#include "index.h"

module m_calcul_forces ! wrap subroutine in module to get arg type check at build time
    use constants
    implicit none
    private :: physical_part_deriv_ijke
contains

    subroutine RK4_attenu_coefs(n_solid,dt,omega_tau_s,agamma,alphaval,betaval,gammaval)
#ifdef SEM_VEC
!$omp declare simd (RK4_attenu_coefs) uniform(dt,n_solid)
#endif
        !- routine returns the coefficients for the time integration of the
        !  relaxation function M(t)
        integer, intent(in)  :: n_solid
        real(fpp), intent(in)  :: dt
        real(fpp), intent(in) :: omega_tau_s,agamma
        real(fpp), intent(out) :: alphaval,betaval,gammaval
        !
        real(fpp) :: dt_tau

        dt_tau = -dt*omega_tau_s
        alphaval = 1d0 + dt_tau + 0.5d0 * dt_tau**2 + dt_tau**3 *(1d0/6d0) + dt_tau**4 *(1d0/24.d0)
        betaval  = dt*(0.5d0 + dt_tau * (1d0/3.d0) + dt_tau**2 *(1d0/8d0) + dt_tau**3 *(1d0/24.d0))
        gammaval = dt*(0.5d0 + dt_tau * (1d0/6.d0) + dt_tau**2 *(1d0/24d0))
    end subroutine RK4_attenu_coefs

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

    pure subroutine physical_part_deriv_ijke(e,i,j,k,ngll,hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
#ifdef SEM_VEC
!$omp declare simd (physical_part_deriv_ijke) linear(e) uniform(i,j,k,ngll,hprime,InvGrad,Scalp)
#endif
        implicit none
        integer,intent(in) :: i,j,k,e
        integer, intent(in) :: ngll
        real, dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real, dimension(IND_MNE(0:2,0:2,0:CHUNK-1)), intent(in) :: InvGrad
        real, dimension(0:CHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real, intent(out) :: dS_dx,dS_dy,dS_dz

        integer :: l
        real :: dS_dxi, dS_deta, dS_dzeta

        dS_dxi   = 0.0D+0
        dS_deta  = 0.0D+0
        dS_dzeta = 0.0D+0
        DO L = 0, ngll-1
            dS_dxi   = dS_dxi  +Scalp(e,L,J,K)*hprime(L,I)
            dS_deta  = dS_deta +Scalp(e,I,L,K)*hprime(L,J)
            dS_dzeta = dS_dzeta+Scalp(e,I,J,L)*hprime(L,K)
        END DO
        !- in the physical domain
        dS_dx = dS_dxi*InvGrad(IND_MNE(0,0,e))+dS_deta*InvGrad(IND_MNE(0,1,e))+dS_dzeta*InvGrad(IND_MNE(0,2,e))
        dS_dy = dS_dxi*InvGrad(IND_MNE(1,0,e))+dS_deta*InvGrad(IND_MNE(1,1,e))+dS_dzeta*InvGrad(IND_MNE(1,2,e))
        dS_dz = dS_dxi*InvGrad(IND_MNE(2,0,e))+dS_deta*InvGrad(IND_MNE(2,1,e))+dS_dzeta*InvGrad(IND_MNE(2,2,e))
    end subroutine physical_part_deriv_ijke

    subroutine calcul_forces_iso(dom,lnum,Fox,Foy,Foz,Depla)
        use champs_solid
        use deriv3d
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: lnum
        real(fpp), dimension(0:CHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:CHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla

        select case(dom%ngll)
        case(4)
            call calcul_forces_iso_4(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case(5)
            call calcul_forces_iso_5(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (6)
            call calcul_forces_iso_6(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (7)
            call calcul_forces_iso_7(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (8)
            call calcul_forces_iso_8(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (9)
            call calcul_forces_iso_9(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case default
            call calcul_forces_iso_n(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        end select
    end subroutine calcul_forces_iso

#undef ATTENUATION
#define NGLLVAL 4
#define PROCNAME calcul_forces_iso_4
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 5
#define PROCNAME calcul_forces_iso_5
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 6
#define PROCNAME calcul_forces_iso_6
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 7
#define PROCNAME calcul_forces_iso_7
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 8
#define PROCNAME calcul_forces_iso_8
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define NGLLVAL 9
#define PROCNAME calcul_forces_iso_9
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#define PROCNAME calcul_forces_iso_n
#include "calcul_forces_solid.inc"

    subroutine calcul_forces_iso_atn(dom,lnum,Fox,Foy,Foz,Depla)
        use champs_solid
        use deriv3d
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: lnum
        real(fpp), dimension(0:CHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real(fpp), dimension(0:CHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla

        select case(dom%ngll)
        case(4)
            call calcul_forces_iso_atn_4(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case(5)
            call calcul_forces_iso_atn_5(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (6)
            call calcul_forces_iso_atn_6(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (7)
            call calcul_forces_iso_atn_7(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (8)
            call calcul_forces_iso_atn_8(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case (9)
            call calcul_forces_iso_atn_9(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        case default
            call calcul_forces_iso_atn_n(dom,dom%ngll,lnum,Fox,Foy,Foz,Depla)
        end select
    end subroutine calcul_forces_iso_atn

#define ATTENUATION
#define NGLLVAL 4
#undef PROCNAME
#define PROCNAME calcul_forces_iso_atn_4
#define PROCNAME_ATN attenuation_update_4
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 5
#define PROCNAME calcul_forces_iso_atn_5
#define PROCNAME_ATN attenuation_update_5
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 6
#define PROCNAME calcul_forces_iso_atn_6
#define PROCNAME_ATN attenuation_update_6
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 7
#define PROCNAME calcul_forces_iso_atn_7
#define PROCNAME_ATN attenuation_update_7
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 8
#define PROCNAME calcul_forces_iso_atn_8
#define PROCNAME_ATN attenuation_update_8
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLLVAL 9
#define PROCNAME calcul_forces_iso_atn_9
#define PROCNAME_ATN attenuation_update_9
#include "calcul_forces_solid.inc"
#undef NGLLVAL
#undef PROCNAME
#undef PROCNAME_ATN
#define NGLL_GEN
#define PROCNAME calcul_forces_iso_atn_n
#define PROCNAME_ATN attenuation_update_n
#define ANISO 0
#define ATN 0
#include "calcul_forces_solid.inc"

    subroutine calcul_forces(dom,lnum,Fox,Foy,Foz,Depla,aniso,n_solid)

        use sdomain
        use attenuation_solid
        implicit none

        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: lnum
        real, dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(out) :: Fox,Foz,Foy
        real, dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2), intent(in) :: Depla
        logical :: aniso
        integer :: n_solid

        integer :: i,j,k,l
        real :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
        real :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
        real :: sxx,sxy,sxz,syy,syz,szz,t4,F1,F2,F3
        real :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
        real :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
        real, parameter :: zero = 0.
        real, dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: t1,t5,t8
        ! Les indices sont reordonnes, probablement pour la localite memoire
        real, dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: t2,t6,t9
        real, dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: t3,t7,t10
        real, dimension(:,:,:,:,:), allocatable :: C
        real, dimension(0:2,0:2) :: invgrad_ijk

        if(aniso) then
            allocate(C(0:5,0:5,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1))
            k = 0
            do i = 0,5
                do j = i,5
                    C(i,j,:,:,:) = dom%Cij_(k,:,:,:,lnum)
                    k = k + 1
                enddo
            enddo
            do i = 1,5
                do j = 0,i-1
                    C(i,j,:,:,:) = C(j,i,:,:,:)
                enddo
            enddo
        end if

        do k = 0,dom%ngll-1
            do j = 0,dom%ngll-1
                do i = 0,dom%ngll-1
                    invgrad_ijk = dom%InvGrad_(:,:,i,j,k,lnum) ! cache for performance

                    call physical_part_deriv_ijk(i,j,k,dom%ngll,dom%hprime,&
                         invgrad_ijk,Depla(:,:,:,0),dxx,dyx,dzx)
                    call physical_part_deriv_ijk(i,j,k,dom%ngll,dom%hprime,&
                         invgrad_ijk,Depla(:,:,:,1),dxy,dyy,dzy)
                    call physical_part_deriv_ijk(i,j,k,dom%ngll,dom%hprime,&
                         invgrad_ijk,Depla(:,:,:,2),dxz,dyz,dzz)

                    call calcul_sigma(dom,i,j,k,lnum,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,&
                         sxx,sxy,sxz,syy,syz,szz,aniso,C,n_solid)

                    xi1 = invgrad_ijk(0,0)
                    xi2 = invgrad_ijk(1,0)
                    xi3 = invgrad_ijk(2,0)
                    et1 = invgrad_ijk(0,1)
                    et2 = invgrad_ijk(1,1)
                    et3 = invgrad_ijk(2,1)
                    ga1 = invgrad_ijk(0,2)
                    ga2 = invgrad_ijk(1,2)
                    ga3 = invgrad_ijk(2,2)

                    !
                    !=====================
                    !       FX
                    !=====================
                    !
                    xt1 = sxx * xi1 + sxy * xi2 + sxz * xi3
                    xt2 = sxx * et1 + sxy * et2 + sxz * et3
                    xt3 = sxx * ga1 + sxy * ga2 + sxz * ga3
                    !
                    !=====================
                    !       FY
                    !=====================
                    !
                    xt5 = syy * xi2 + sxy * xi1 + syz * xi3
                    xt6 = syy * et2 + sxy * et1 + syz * et3
                    xt7 = syy * ga2 + sxy * ga1 + syz * ga3
                    !
                    !=====================
                    !       FZ
                    !=====================
                    !
                    xt8  = szz * xi3 + sxz * xi1 + syz * xi2
                    xt9  = szz * et3 + sxz * et1 + syz * et2
                    xt10 = szz * ga3 + sxz * ga1 + syz * ga2

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

                    t1(i,j,k) = xt1
                    t5(i,j,k) = xt5
                    t8(i,j,k) = xt8

                    t2(j,i,k) = xt2
                    t6(j,i,k) = xt6
                    t9(j,i,k) = xt9

                    t3(k,i,j) = xt3
                    t7(k,i,j) = xt7
                    t10(k,i,j) = xt10
                enddo
            enddo
        enddo
        if (allocated(C)) deallocate(C)

        !
        !- Multiplication par la matrice de derivation puis par les poids
        !

        !=-=-=-=-=-=-=-=-=-=-
        do k = 0,dom%ngll-1
            do j = 0,dom%ngll-1
                do i = 0,dom%ngll-1
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
                    do l = 0,dom%ngll-1
                        t41 = t41 + dom%htprime(l,i) * t1(l,j,k)
                        t42 = t42 + dom%htprime(l,i) * t5(l,j,k)
                        t43 = t43 + dom%htprime(l,i) * t8(l,j,k)
                    enddo

                    do l = 0,dom%ngll-1
                        t51 = t51 + dom%htprime(l,j) * t2(l,i,k)
                        t52 = t52 + dom%htprime(l,j) * t6(l,i,k)
                        t53 = t53 + dom%htprime(l,j) * t9(l,i,k)
                    enddo
                    ! FX
                    F1 = t41*t11 + t51*t12
                    ! FY
                    F2 = t42*t11 + t52*t12
                    ! FZ
                    F3 = t43*t11 + t53*t12
                    !
                    !
                    do l = 0,dom%ngll-1
                        t61 = t61 + dom%htprime(l,k) * t3(l,i,j)
                        t62 = t62 + dom%htprime(l,k) * t7(l,i,j)
                        t63 = t63 + dom%htprime(l,k) * t10(l,i,j)
                    enddo

                    ! FX
                    F1 = F1 + t61*t13
                    ! FY
                    F2 = F2 + t62*t13
                    ! FZ
                    F3 = F3 + t63*t13
                    !
                    Fox(i,j,k) = F1
                    Foy(i,j,k) = F2
                    Foz(i,j,k) = F3

                    !=-=-=-=-=-=-=-=-=-=-
                enddo
            enddo
        enddo
        !=-=-=-=-=-=-=-=-=-=-
        call attenuation_update(dom,lnum,dom%hprime,Depla,n_solid,aniso)
    end subroutine calcul_forces

    subroutine calcul_sigma(dom,i,j,k,lnum,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,&
                            sxx,sxy,sxz,syy,syz,szz,aniso,C,n_solid)
        use sdomain
        use attenuation_solid
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: i,j,k,lnum
        real, intent(in) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
        real, intent(out) :: sxx,sxy,sxz,syy,syz,szz
        logical aniso
        real, dimension(:,:,:,:,:), allocatable :: C
        integer :: n_solid

        real, dimension(0:5) :: eij
        real, parameter :: s2 = 1.414213562373095, s2o2 = 0.707106781186547
        real xmu, xla, xla2mu
        integer l

        if (aniso .and. allocated(C)) then
            eij(0) = DXX
            eij(1) = DYY
            eij(2) = DZZ
            eij(3) = s2o2*(DYZ+DZY)
            eij(4) = s2o2*(DXZ+DZX)
            eij(5) = s2o2*(DXY+DYX)

            sxx = 0.
            syy = 0.
            szz = 0.
            syz = 0.
            sxz = 0.
            sxy = 0.
            do l = 0,5
                sxx = sxx + C(0,l,i,j,k)*eij(l)
                syy = syy + C(1,l,i,j,k)*eij(l)
                szz = szz + C(2,l,i,j,k)*eij(l)
                syz = syz + C(3,l,i,j,k)*eij(l)
                sxz = sxz + C(4,l,i,j,k)*eij(l)
                sxy = sxy + C(5,l,i,j,k)*eij(l)
            enddo
            syz=syz/s2
            sxz=sxz/s2
            sxy=sxy/s2
        else
            if (n_solid>0) then
                call calcul_sigma_attenuation(dom,i,j,k,lnum,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,&
                     sxx,sxy,sxz,syy,syz,szz)
            else
                xmu = dom%Mu_    (i,j,k,lnum)
                xla = dom%Lambda_(i,j,k,lnum)
                xla2mu = xla + 2. * xmu

                sxx = xla2mu * DXX + xla * ( DYY + DZZ )
                sxy = xmu * ( DXY + DYX )
                sxz = xmu * ( DXZ + DZX )
                syy = xla2mu * DYY + xla * ( DXX + DZZ )
                syz = xmu * ( DYZ + DZY )
                szz = xla2mu * DZZ + xla * ( DXX + DYY )
            end if
        end if

        if (n_solid>0) then
            call sigma_attenuation(dom,i,j,k,lnum,DXX,DYY,DZZ,&
                 sxx,sxy,sxz,syy,syz,szz,n_solid)
        end if
    end subroutine calcul_sigma
end module m_calcul_forces

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
