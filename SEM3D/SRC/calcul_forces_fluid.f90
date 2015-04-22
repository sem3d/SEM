subroutine calcul_forces_fluid(FFl,invgrad,      &
    hTprimex,hTprimey,hTprimez,jac,wheix,wheiy,wheiz,  &
    dPhiX,dPhiY,dPhiZ,dens_,ngllx,nglly,ngllz)

    use sdomain

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: invgrad
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: jac,dens_, dPhiX,dPhiY,dPhiZ
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: FFl
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
    real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hTprimey
    real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hTprimez
    real, dimension(0:ngllx-1), intent(in) :: wheix
    real, dimension(0:nglly-1), intent(in) :: wheiy
    real, dimension(0:ngllz-1), intent(in) :: wheiz

    real :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    integer :: i,j,k,l
    real :: sx,sy,sz,t4,F1
    real :: t41,t11,t51,t12,t61,t13
    real :: xt1,xt6,xt10
    real, parameter :: zero = 0.
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: xdens
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t6
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t10


    xdens(:,:,:) = 1d0/dens_(:,:,:)

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1

                ! (fluid equivalent) stress  ( = physical velocity)
                sx = xdens(i,j,k)*dPhiX(i,j,k)
                sy = xdens(i,j,k)*dPhiY(i,j,k)
                sz = xdens(i,j,k)*dPhiZ(i,j,k)

                xi1 = Invgrad(i,j,k,0,0)
                xi2 = Invgrad(i,j,k,1,0)
                xi3 = Invgrad(i,j,k,2,0)
                et1 = Invgrad(i,j,k,0,1)
                et2 = Invgrad(i,j,k,1,1)
                et3 = Invgrad(i,j,k,2,1)
                ga1 = Invgrad(i,j,k,0,2)
                ga2 = Invgrad(i,j,k,1,2)
                ga3 = Invgrad(i,j,k,2,2)

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
                t4 = jac(i,j,k) * wheix(i)
                xt1  =  xt1 * t4

                t4 = jac(i,j,k) * wheiy(j)
                xt6  =  xt6 * t4

                t4 = jac(i,j,k) * wheiz(k)
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
                t11 = wheiy(j) * wheiz(k)
                t12 = wheix(i) * wheiz(k)
                t13 = wheix(i) * wheiy(j)
                !
                t41 = zero
                t51 = zero
                t61 = zero
                !
                do l = 0,ngllx-1
                    t41 = t41 + hTprimex(l,i) * t1(l,j,k)
                enddo

                do l = 0,nglly-1
                    t51 = t51 + hTprimey(l,j) * t6(l,i,k)
                enddo
                ! FFl
                F1 = t41*t11 + t51*t12
                !
                !
                do l = 0,ngllz-1
                    t61 = t61 + hTprimez(l,k) * t10(l,i,j)
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
