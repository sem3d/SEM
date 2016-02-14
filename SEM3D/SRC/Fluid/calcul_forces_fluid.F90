!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module m_calcul_forces_fluid ! wrap subroutine in module to get arg type check at build time
    contains
subroutine calcul_forces_fluid(dom,lnum,FFl,hTprimex,wheix,dPhiX,dPhiY,dPhiZ,ngllx,nglly,ngllz)
    use sdomain
    implicit none
#include "index.h"

    type(domain_fluid), intent (INOUT) :: dom
    integer, intent(in) :: lnum
    integer, intent(in) :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: FFl
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
    real, dimension(0:ngllx-1), intent(in) :: wheix
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: dPhiX,dPhiY,dPhiZ

    real :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    integer :: i,j,k,l
    real :: sx,sy,sz,t4,F1
    real :: t41,t11,t51,t12,t61,t13
    real :: xt1,xt6,xt10
    real, parameter :: zero = 0.
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: xdens
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t6
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t10

    xdens(:,:,:) = 1d0/dom%Density_(:,:,:,lnum)

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1

                ! (fluid equivalent) stress  ( = physical velocity)
                sx = xdens(i,j,k)*dPhiX(i,j,k)
                sy = xdens(i,j,k)*dPhiY(i,j,k)
                sz = xdens(i,j,k)*dPhiZ(i,j,k)

                xi1 = dom%InvGrad_(0,0,i,j,k,lnum)
                xi2 = dom%InvGrad_(1,0,i,j,k,lnum)
                xi3 = dom%InvGrad_(2,0,i,j,k,lnum)
                et1 = dom%InvGrad_(0,1,i,j,k,lnum)
                et2 = dom%InvGrad_(1,1,i,j,k,lnum)
                et3 = dom%InvGrad_(2,1,i,j,k,lnum)
                ga1 = dom%InvGrad_(0,2,i,j,k,lnum)
                ga2 = dom%InvGrad_(1,2,i,j,k,lnum)
                ga3 = dom%InvGrad_(2,2,i,j,k,lnum)

                !=====================
                !       F1 
                xt1 = sx*xi1 + sy*xi2 + sz*xi3

                !=====================
                !       F2 
                xt6 = sx*et1 + sy*et2 + sz*et3

                !=====================
                !       F3 
                xt10 = sx*ga1 + sy*ga2 + sz*ga3

                !
                !- Multiply par Jacobian and weight
                !
                t4  = dom%Jacob_(i,j,k,lnum) * wheix(i)
                xt1 = xt1 * t4

                t4  = dom%Jacob_(i,j,k,lnum) * wheix(j)
                xt6 = xt6 * t4

                t4   = dom%Jacob_(i,j,k,lnum) * wheix(k)
                xt10 = xt10 * t4


                t1(i,j,k) = xt1

                t6(j,i,k) = xt6

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
                t11 = wheix(j) * wheix(k)
                t12 = wheix(i) * wheix(k)
                t13 = wheix(i) * wheix(j)
                !
                t41 = zero
                t51 = zero
                t61 = zero
                !
                do l = 0,ngllx-1
                    t41 = t41 + hTprimex(l,i) * t1(l,j,k)
                enddo

                do l = 0,nglly-1
                    t51 = t51 + hTprimex(l,j) * t6(l,i,k)
                enddo
                ! FFl
                F1 = t41*t11 + t51*t12
                !
                !
                do l = 0,ngllz-1
                    t61 = t61 + hTprimex(l,k) * t10(l,i,j)
                enddo

                ! FX
                F1 = F1 + t61*t13
                !
                FFl(i,j,k) = F1

                !=-=-=-=-=-=-=-=-=-=-
            enddo
        enddo
    enddo
    !=-=-=-=-=-=-=-=-=-=-

end subroutine calcul_forces_fluid
end module m_calcul_forces_fluid

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
