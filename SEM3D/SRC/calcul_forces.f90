!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file calcul_forces.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine calcul_forces(Fox,Foy,Foz, invgrad, &
    dx,dy,dz, jac, poidsx,poidsy,poidsz, DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
    mu_,la_, ngllx,nglly,ngllz)

    use sdomain

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz
    real, dimension(0:2,0:2,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: invgrad
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: jac, mu_,la_
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Fox,Foz,Foy
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: dx
    real, dimension(0:nglly-1,0:nglly-1), intent(in) :: dy
    real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: dz
    real, dimension(0:ngllx-1), intent(in) :: poidsx
    real, dimension(0:nglly-1), intent(in) :: poidsy
    real, dimension(0:ngllz-1), intent(in) :: poidsz

    integer :: i,j,k,l

    real :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    real :: sxx,sxy,sxz,syy,syz,szz,t4,F1,F2,F3
    real :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
    real :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
    real, parameter :: zero = 0.
    real :: xmu, xla, xla2mu
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t5,t8
    ! Les indices sont reordonnes, probablement pour la localite memoire
    real, dimension(0:nglly-1,0:ngllx-1,0:ngllz-1) :: t2,t6,t9
    real, dimension(0:ngllz-1,0:ngllx-1,0:nglly-1) :: t3,t7,t10




    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1

                xmu = mu_(i,j,k)
                xla = la_(i,j,k)
                xla2mu = xla + 2. * xmu

                sxx = xla2mu * DXX(i,j,k) + xla  * ( DYY(i,j,k) + DZZ(i,j,k) )
                !---
                sxy = xmu  *( DXY(i,j,k)  + DYX(i,j,k)  )
                !---
                sxz = xmu * ( DXZ(i,j,k)  + DZX(i,j,k)  )
                !---
                syy = xla2mu * DYY(i,j,k) + xla  * ( DXX(i,j,k) + DZZ(i,j,k) )
                !---
                syz = xmu * ( DYZ(i,j,k)  + DZY(i,j,k)  )
                !---
                szz = xla2mu * DZZ(i,j,k) + xla  * ( DXX(i,j,k) + DYY(i,j,k) )
                !---

                xi1 = Invgrad(0,0,i,j,k)
                xi2 = Invgrad(1,0,i,j,k)
                xi3 = Invgrad(2,0,i,j,k)
                et1 = Invgrad(0,1,i,j,k)
                et2 = Invgrad(1,1,i,j,k)
                et3 = Invgrad(2,1,i,j,k)
                ga1 = Invgrad(0,2,i,j,k)
                ga2 = Invgrad(1,2,i,j,k)
                ga3 = Invgrad(2,2,i,j,k)

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
                t4 = jac(i,j,k) * poidsx(i)
                xt1  =  xt1 * t4
                xt5  =  xt5 * t4
                xt8  =  xt8 * t4

                t4 = jac(i,j,k) * poidsy(j)
                xt2  =  xt2 * t4
                xt6  =  xt6 * t4
                xt9  =  xt9 * t4

                t4 = jac(i,j,k) * poidsz(k)
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
                t11 = poidsy(j) * poidsz(k)
                t12 = poidsx(i) * poidsz(k)
                t13 = poidsx(i) * poidsy(j)
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
                    t41 = t41 + dx(l,i) * t1(l,j,k)
                    t42 = t42 + dx(l,i) * t5(l,j,k)
                    t43 = t43 + dx(l,i) * t8(l,j,k)
                enddo

                do l = 0,nglly-1
                    t51 = t51 + dy(l,j) * t2(l,i,k)
                    t52 = t52 + dy(l,j) * t6(l,i,k)
                    t53 = t53 + dy(l,j) * t9(l,i,k)
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
                    t61 = t61 + dz(l,k) * t3(l,i,j)
                    t62 = t62 + dz(l,k) * t7(l,i,j)
                    t63 = t63 + dz(l,k) * t10(l,i,j)
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

end subroutine calcul_forces

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
