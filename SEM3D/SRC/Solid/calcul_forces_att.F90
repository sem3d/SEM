!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file calcul_forces_att.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<
module m_calcul_forces_att ! wrap subroutine in module to get arg type check at build time
    contains
subroutine calcul_forces_att(dom,lnum,Fox,Foy,Foz, &
    htprimex,GLLwx,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
    ngllx,nglly,ngllz,n_solid)
    use sdomain

    implicit none
#include "index.h"

    type(domain_solid), intent (INOUT) :: dom
    integer, intent(in) :: lnum
    integer, intent(in) :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Fox,Foz,Foy
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: htprimex
    real, dimension(0:ngllx-1), intent(in) :: GLLwx
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    integer, intent(in) :: n_solid

    real :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    real :: sxx,sxy,sxz,syy,syz,szz,t4,F1,F2,F3
    real :: xpression , stt
    real :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
    real :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
    real, parameter :: zero = 0.
    integer :: i,j,k,l, i_sls
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: xmu,x2mu,xkappa
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t5,t8,t2,t6,t9
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t3,t7,t10

    xmu(:,:,:) = dom%Mu_(:,:,:,lnum)
    !  mu_relaxed -> mu_unrelaxed
    xmu(:,:,:) = xmu(:,:,:) * dom%onemSbeta(:,:,:,lnum)
    xkappa(:,:,:) = dom%Kappa_(:,:,:,lnum)
    !  kappa_relaxed -> kappa_unrelaxed
    xkappa(:,:,:) = xkappa(:,:,:) * dom%onemPbeta(:,:,:,lnum)
    x2mu(:,:,:) = 2. * xmu(:,:,:)

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                xpression = xkappa(i,j,k) * ( DXX(i,j,k) + &
                    DYY(i,j,k) + DZZ(i,j,k) )

                sxx = x2mu(i,j,k) * DXX(i,j,k)

                sxy = xmu(i,j,k)  *( DXY(i,j,k)  +       &
                    DYX(i,j,k)  )

                sxz = xmu(i,j,k) * ( DXZ(i,j,k)  +       &
                    DZX(i,j,k)  )

                syy = x2mu(i,j,k) * DYY(i,j,k)

                syz = xmu(i,j,k) * ( DYZ(i,j,k)  +       &
                    DZY(i,j,k)  )

                szz = x2mu(i,j,k) * DZZ(i,j,k)

                !        stt = (sxx + syy + szz)/3.
                do i_sls = 0,n_solid-1
                    xpression = xpression - dom%R_vol_(i_sls,i,j,k,lnum)
                    sxx = sxx  - dom%R_xx_(i_sls,i,j,k,lnum)
                    syy = syy  - dom%R_yy_(i_sls,i,j,k,lnum)
                    !         sxx = sxx - stt - R_xx_(i_sls,i,j,k)
                    !         syy = syy - stt - R_yy_(i_sls,i,j,k)
                    ! ici on utilise le fait que la trace est nulle
                    szz = szz  + dom%R_xx_(i_sls,i,j,k,lnum) + dom%R_yy_(i_sls,i,j,k,lnum)
                    !         szz = szz - stt + R_xx_(i_sls,i,j,k) + R_yy_(i_sls,i,j,k)
                    sxy = sxy - dom%R_xy_(i_sls,i,j,k,lnum)
                    sxz = sxz - dom%R_xz_(i_sls,i,j,k,lnum)
                    syz = syz - dom%R_yz_(i_sls,i,j,k,lnum)
                enddo
                stt = (sxx + syy + szz)/3.
                sxx = sxx - stt + xpression
                syy = syy - stt + xpression
                szz = szz - stt + xpression

                xi1 = dom%InvGrad_(0,0,i,j,k,lnum)
                xi2 = dom%InvGrad_(1,0,i,j,k,lnum)
                xi3 = dom%InvGrad_(2,0,i,j,k,lnum)
                et1 = dom%InvGrad_(0,1,i,j,k,lnum)
                et2 = dom%InvGrad_(1,1,i,j,k,lnum)
                et3 = dom%InvGrad_(2,1,i,j,k,lnum)
                ga1 = dom%InvGrad_(0,2,i,j,k,lnum)
                ga2 = dom%InvGrad_(1,2,i,j,k,lnum)
                ga3 = dom%InvGrad_(2,2,i,j,k,lnum)

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
                t4 = dom%Jacob_(i,j,k,lnum) * GLLwx(i)
                xt1  =  xt1 * t4
                xt5  =  xt5 * t4
                xt8  =  xt8 * t4

                t4 = dom%Jacob_(i,j,k,lnum) * GLLwx(j)
                xt2  =  xt2 * t4
                xt6  =  xt6 * t4
                xt9  =  xt9 * t4

                t4 = dom%Jacob_(i,j,k,lnum) * GLLwx(k)
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

    !
    !- Multiplication par la matrice de derivation puis par les poids
    !

    !=-=-=-=-=-=-=-=-=-=-
    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                !=-=-=-=-=-=-=-=-=-=-
                !
                t11 = GLLwx(j) * GLLwx(k)
                t12 = GLLwx(i) * GLLwx(k)
                t13 = GLLwx(i) * GLLwx(j)
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
                do l = 0,ngllx-1
                    t41 = t41 + htprimex(l,i) * t1(l,j,k)
                    t42 = t42 + htprimex(l,i) * t5(l,j,k)
                    t43 = t43 + htprimex(l,i) * t8(l,j,k)
                enddo

                do l = 0,nglly-1
                    t51 = t51 + htprimex(l,j) * t2(l,i,k)
                    t52 = t52 + htprimex(l,j) * t6(l,i,k)
                    t53 = t53 + htprimex(l,j) * t9(l,i,k)
                enddo
                ! FX
                F1 = t41*t11 + t51*t12
                ! FY
                F2 = t42*t11 + t52*t12
                ! FZ
                F3 = t43*t11 + t53*t12
                !
                !
                do l = 0,ngllz-1
                    t61 = t61 + htprimex(l,k) * t3(l,i,j)
                    t62 = t62 + htprimex(l,k) * t7(l,i,j)
                    t63 = t63 + htprimex(l,k) * t10(l,i,j)
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

end subroutine calcul_forces_att
end module m_calcul_forces_att

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
