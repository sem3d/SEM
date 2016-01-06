!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solidpml
    use constants
    use sdomain
    implicit none

contains

    subroutine forces_int_sol_pml(Elem, mat, champs1)
        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        type(champs), intent(inout) :: champs1
        !
        integer :: m1, m2, m3
        integer :: i, j, k, l, ind
        real :: sum_vx, sum_vy, sum_vz, acoeff
        real, dimension(0:2,0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1)  :: Forces1, Forces2, Forces3

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz


        do k = 0,m3-1
            do j = 0,m2-1
                do i=0,m1-1
                    sum_vx = 0d0
                    sum_vy = 0d0
                    sum_vz = 0d0
                    do l = 0,m1-1
                        acoeff = - mat%hprimex(i,l)*mat%GLLwx(l)*mat%GLLwy(j)*mat%GLLwz(k)*Elem%Jacob(l,j,k)
                        sum_vx = sum_vx + acoeff*Elem%InvGrad(0,0,l,j,k)*Elem%slpml%Diagonal_Stress(l,j,k,0)
                        sum_vx = sum_vx + acoeff*Elem%InvGrad(1,0,l,j,k)*Elem%slpml%Residual_Stress(l,j,k,0)
                        sum_vx = sum_vx + acoeff*Elem%InvGrad(2,0,l,j,k)*Elem%slpml%Residual_Stress(l,j,k,1)

                        sum_vy = sum_vy + acoeff*Elem%InvGrad(0,0,l,j,k)*Elem%slpml%Residual_Stress(l,j,k,0)
                        sum_vy = sum_vy + acoeff*Elem%InvGrad(1,0,l,j,k)*Elem%slpml%Diagonal_Stress(l,j,k,1)
                        sum_vy = sum_vy + acoeff*Elem%InvGrad(2,0,l,j,k)*Elem%slpml%Residual_Stress(l,j,k,2)

                        sum_vz = sum_vz + acoeff*Elem%InvGrad(0,0,l,j,k)*Elem%slpml%Residual_Stress(l,j,k,1)
                        sum_vz = sum_vz + acoeff*Elem%InvGrad(1,0,l,j,k)*Elem%slpml%Residual_Stress(l,j,k,2)
                        sum_vz = sum_vz + acoeff*Elem%InvGrad(2,0,l,j,k)*Elem%slpml%Diagonal_Stress(l,j,k,2)
                    end do
                    Forces1(0,i,j,k) = sum_vx
                    Forces1(1,i,j,k) = sum_vy
                    Forces1(2,i,j,k) = sum_vz
                end do
            end do
        end do

        do k = 0,m3-1
            Forces2(:,:,:,k) = 0d0
            do l = 0,m2-1
                do j = 0,m2-1
                    do i=0,m1-1
                        acoeff = - mat%hprimey(j,l)*mat%GLLwx(i)*mat%GLLwy(l)*mat%GLLwz(k)*Elem%Jacob(i,l,k)
                        sum_vx = acoeff*(Elem%InvGrad(0,1,i,l,k)*Elem%slpml%Diagonal_Stress(i,l,k,0) + &
                                         Elem%InvGrad(1,1,i,l,k)*Elem%slpml%Residual_Stress(i,l,k,0) + &
                                         Elem%InvGrad(2,1,i,l,k)*Elem%slpml%Residual_Stress(i,l,k,1))

                        sum_vy = acoeff*(Elem%InvGrad(0,1,i,l,k)*Elem%slpml%Residual_Stress(i,l,k,0) + &
                                         Elem%InvGrad(1,1,i,l,k)*Elem%slpml%Diagonal_Stress(i,l,k,1) + &
                                         Elem%InvGrad(2,1,i,l,k)*Elem%slpml%Residual_Stress(i,l,k,2))

                        sum_vz = acoeff*(Elem%InvGrad(0,1,i,l,k)*Elem%slpml%Residual_Stress(i,l,k,1) + &
                                         Elem%InvGrad(1,1,i,l,k)*Elem%slpml%Residual_Stress(i,l,k,2) + &
                                         Elem%InvGrad(2,1,i,l,k)*Elem%slpml%Diagonal_Stress(i,l,k,2))
                        Forces2(0,i,j,k) = Forces2(0,i,j,k) + sum_vx
                        Forces2(1,i,j,k) = Forces2(1,i,j,k) + sum_vy
                        Forces2(2,i,j,k) = Forces2(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        ! TODO reorder loops ?
        Forces3 = 0
        do l = 0,m3-1
            do k = 0,m3-1
                do j = 0,m2-1
                    do i=0,m1-1
                        acoeff = - mat%hprimez(k,l)*mat%GLLwx(i)*mat%GLLwy(j)*mat%GLLwz(l)*Elem%Jacob(i,j,l)
                        sum_vx = acoeff*(Elem%InvGrad(0,2,i,j,l)*Elem%slpml%Diagonal_Stress(i,j,l,0) + &
                                         Elem%InvGrad(1,2,i,j,l)*Elem%slpml%Residual_Stress(i,j,l,0) + &
                                         Elem%InvGrad(2,2,i,j,l)*Elem%slpml%Residual_Stress(i,j,l,1))

                        sum_vy = acoeff*(Elem%InvGrad(0,2,i,j,l)*Elem%slpml%Residual_Stress(i,j,l,0) + &
                                         Elem%InvGrad(1,2,i,j,l)*Elem%slpml%Diagonal_Stress(i,j,l,1) + &
                                         Elem%InvGrad(2,2,i,j,l)*Elem%slpml%Residual_Stress(i,j,l,2))

                        sum_vz = acoeff*(Elem%InvGrad(0,2,i,j,l)*Elem%slpml%Residual_Stress(i,j,l,1) + &
                                         Elem%InvGrad(1,2,i,j,l)*Elem%slpml%Residual_Stress(i,j,l,2) + &
                                         Elem%InvGrad(2,2,i,j,l)*Elem%slpml%Diagonal_Stress(i,j,l,2))
                        Forces3(0,i,j,k) = Forces3(0,i,j,k) + sum_vx
                        Forces3(1,i,j,k) = Forces3(1,i,j,k) + sum_vy
                        Forces3(2,i,j,k) = Forces3(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do


        ! Assemblage
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    ind = Elem%Idom(i,j,k)
                    champs1%ForcesPML(ind,:,0) = champs1%ForcesPML(ind,:,0) + Forces1(:,i,j,k)
                    champs1%ForcesPML(ind,:,1) = champs1%ForcesPML(ind,:,1) + Forces2(:,i,j,k)
                    champs1%ForcesPML(ind,:,2) = champs1%ForcesPML(ind,:,2) + Forces3(:,i,j,k)
                enddo
            enddo
        enddo
    end subroutine forces_int_sol_pml

    subroutine pred_sol_pml(Elem, mat, dt, champs1)
        implicit none

        type(Element), intent(inout) :: Elem
        type (subdomain), intent(IN) :: mat
        type(champs), intent(inout) :: champs1
        real, intent(in) :: dt
        !
        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVx_dx, dVx_dy, dVx_dz
        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVy_dx, dVy_dy, dVy_dz
        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVz_dx, dVz_dy, dVz_dz
        integer :: m1, m2, m3
        integer :: i, j, k, ind, i_dir
        real, dimension (:,:,:,:), allocatable :: Veloc

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        allocate(Veloc(0:m1-1,0:m2-1,0:m3-1,0:2))
        do i_dir = 0,2
            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                        ind = Elem%Idom(i,j,k)
                        Veloc(i,j,k,i_dir) = champs1%VelocPML(ind,i_dir,0) + &
                            champs1%VelocPML(ind,i_dir,1) + &
                            champs1%VelocPML(ind,i_dir,2)
                    enddo
                enddo
            enddo
        enddo

        ! partial of velocity components with respect to xi,eta,zeta
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,Elem%InvGrad, Veloc(:,:,:,0),dVx_dx,dVx_dy,dVx_dz)
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,Elem%InvGrad, Veloc(:,:,:,1),dVy_dx,dVy_dy,dVy_dz)
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,Elem%InvGrad, Veloc(:,:,:,2),dVz_dx,dVz_dy,dVz_dz)

        deallocate(Veloc)

        ! Stress_xx
        Elem%slpml%Diagonal_Stress1(:,:,:,0) = Elem%xpml%DumpSx(:,:,:,0)*Elem%slpml%Diagonal_Stress1(:,:,:,0) + Elem%xpml%DumpSx(:,:,:,1)*Dt*(Elem%Lambda+2*Elem%Mu)*dVx_dx
        Elem%slpml%Diagonal_Stress2(:,:,:,0) = Elem%xpml%DumpSy(:,:,:,0)*Elem%slpml%Diagonal_Stress2(:,:,:,0) + Elem%xpml%DumpSy(:,:,:,1)*Dt*(Elem%Lambda)*dVy_dy
        Elem%slpml%Diagonal_Stress3(:,:,:,0) = Elem%xpml%DumpSz(:,:,:,0)*Elem%slpml%Diagonal_Stress3(:,:,:,0) + Elem%xpml%DumpSz(:,:,:,1)*Dt*(Elem%Lambda)*dVz_dz

        ! Stress_yy
        Elem%slpml%Diagonal_Stress1(:,:,:,1) = Elem%xpml%DumpSx(:,:,:,0)*Elem%slpml%Diagonal_Stress1(:,:,:,1) + Elem%xpml%DumpSx(:,:,:,1)*Dt*(Elem%Lambda)*dVx_dx
        Elem%slpml%Diagonal_Stress2(:,:,:,1) = Elem%xpml%DumpSy(:,:,:,0)*Elem%slpml%Diagonal_Stress2(:,:,:,1) + Elem%xpml%DumpSy(:,:,:,1)*Dt*(Elem%Lambda+2*Elem%Mu)*dVy_dy
        Elem%slpml%Diagonal_Stress3(:,:,:,1) = Elem%xpml%DumpSz(:,:,:,0)*Elem%slpml%Diagonal_Stress3(:,:,:,1) + Elem%xpml%DumpSz(:,:,:,1)*Dt*(Elem%Lambda)*dVz_dz

        ! Stress_zz
        Elem%slpml%Diagonal_Stress1(:,:,:,2) = Elem%xpml%DumpSx(:,:,:,0)*Elem%slpml%Diagonal_Stress1(:,:,:,2) + Elem%xpml%DumpSx(:,:,:,1)*Dt*(Elem%Lambda)*dVx_dx
        Elem%slpml%Diagonal_Stress2(:,:,:,2) = Elem%xpml%DumpSy(:,:,:,0)*Elem%slpml%Diagonal_Stress2(:,:,:,2) + Elem%xpml%DumpSy(:,:,:,1)*Dt*(Elem%Lambda)*dVy_dy
        Elem%slpml%Diagonal_Stress3(:,:,:,2) = Elem%xpml%DumpSz(:,:,:,0)*Elem%slpml%Diagonal_Stress3(:,:,:,2) + Elem%xpml%DumpSz(:,:,:,1)*Dt*(Elem%Lambda+2*Elem%Mu)*dVz_dz
        Elem%slpml%Diagonal_Stress = Elem%slpml%Diagonal_Stress1 + Elem%slpml%Diagonal_Stress2 + Elem%slpml%Diagonal_Stress3

        ! Stress_xy
        Elem%slpml%Residual_Stress1 (:,:,:,0) = Elem%xpml%DumpSx(:,:,:,0)*Elem%slpml%Residual_Stress1 (:,:,:,0) + Elem%xpml%DumpSx(:,:,:,1)*Dt*(Elem%Mu)*dVy_dx
        Elem%slpml%Residual_Stress2 (:,:,:,0) = Elem%xpml%DumpSy(:,:,:,0)*Elem%slpml%Residual_Stress2 (:,:,:,0) + Elem%xpml%DumpSy(:,:,:,1)*Dt*(Elem%Mu)*dVx_dy
        Elem%slpml%Residual_Stress3 (:,:,:,0) = Elem%xpml%DumpSz(:,:,:,0)*Elem%slpml%Residual_Stress3 (:,:,:,0)

        ! Stress_xz
        Elem%slpml%Residual_Stress1 (:,:,:,1) = Elem%xpml%DumpSx(:,:,:,0)*Elem%slpml%Residual_Stress1 (:,:,:,1) + Elem%xpml%DumpSx(:,:,:,1)*Dt*(Elem%Mu)*dVz_dx
        Elem%slpml%Residual_Stress2 (:,:,:,1) = Elem%xpml%DumpSy(:,:,:,0)*Elem%slpml%Residual_Stress2 (:,:,:,1)
        Elem%slpml%Residual_Stress3 (:,:,:,1) = Elem%xpml%DumpSz(:,:,:,0)*Elem%slpml%Residual_Stress3 (:,:,:,1) + Elem%xpml%DumpSz(:,:,:,1)*Dt*(Elem%Mu)*dVx_dz

        ! Stress_yz
        Elem%slpml%Residual_Stress1 (:,:,:,2) = Elem%xpml%DumpSx(:,:,:,0)*Elem%slpml%Residual_Stress1 (:,:,:,2)
        Elem%slpml%Residual_Stress2 (:,:,:,2) = Elem%xpml%DumpSy(:,:,:,0)*Elem%slpml%Residual_Stress2 (:,:,:,2) + Elem%xpml%DumpSy(:,:,:,1)*Dt*(Elem%Mu)*dVz_dy
        Elem%slpml%Residual_Stress3 (:,:,:,2) = Elem%xpml%DumpSz(:,:,:,0)*Elem%slpml%Residual_Stress3 (:,:,:,2) + Elem%xpml%DumpSz(:,:,:,1)*Dt*(Elem%Mu)*dVy_dz

        Elem%slpml%Residual_Stress = Elem%slpml%Residual_Stress1 + Elem%slpml%Residual_Stress2 + Elem%slpml%Residual_Stress3

    end subroutine pred_sol_pml

end module dom_solidpml

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
