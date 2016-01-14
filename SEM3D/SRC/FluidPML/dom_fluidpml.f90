!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_fluidpml
    use constants
    use champs_fluidpml
    use selement
    use sdomain
    use ssubdomains
    implicit none

contains

  subroutine forces_int_flu_pml(dom, mat, champs1, Elem, lnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        type (subdomain), intent(IN) :: mat
        type(champsfluidpml), intent(inout) :: champs1
        type (Element), intent(INOUT) :: Elem
        integer :: lnum
        !
        integer :: m1, m2, m3
        integer :: i, j, k, l, ind
        real :: sum_vx, sum_vy, sum_vz, acoeff
        real, dimension(0:2,0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1)  :: ForcesFl

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        do k = 0,m3-1
            do j = 0,m2-1
                do i=0,m1-1
                    ind = Elem%Idom(i,j,k)
                    sum_vx = 0d0
                    sum_vy = 0d0
                    sum_vz = 0d0
                    do l = 0,m1-1
                        acoeff = - mat%hprimex(i,l)*mat%GLLwx(l)*mat%GLLwy(j)*mat%GLLwz(k)*Elem%Jacob(l,j,k)
                        sum_vx = sum_vx + acoeff*Elem%InvGrad(0,0,l,j,k)*dom%Veloc(l,j,k,0,lnum)
                        sum_vy = sum_vy + acoeff*Elem%InvGrad(1,0,l,j,k)*dom%Veloc(l,j,k,1,lnum)
                        sum_vz = sum_vz + acoeff*Elem%InvGrad(2,0,l,j,k)*dom%Veloc(l,j,k,2,lnum)
                    end do
                    ForcesFl(0,i,j,k) = sum_vx
                    ForcesFl(1,i,j,k) = sum_vy
                    ForcesFl(2,i,j,k) = sum_vz
                end do
            end do
        end do
        do k = 0,m3-1
            do l = 0,m2-1
                do j = 0,m2-1
                    do i=0,m1-1
                        ind = Elem%Idom(i,j,k)
                        acoeff = - mat%hprimey(j,l)*mat%GLLwx(i)*mat%GLLwy(l)*mat%GLLwz(k)*Elem%Jacob(i,l,k)
                        sum_vx = acoeff*Elem%InvGrad(0,1,i,l,k)*dom%Veloc(i,l,k,0,lnum)
                        sum_vy = acoeff*Elem%InvGrad(1,1,i,l,k)*dom%Veloc(i,l,k,1,lnum)
                        sum_vz = acoeff*Elem%InvGrad(2,1,i,l,k)*dom%Veloc(i,l,k,2,lnum)
                        ForcesFl(0,i,j,k) = ForcesFl(0,i,j,k) + sum_vx
                        ForcesFl(1,i,j,k) = ForcesFl(1,i,j,k) + sum_vy
                        ForcesFl(2,i,j,k) = ForcesFl(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        ! TODO reorder loops ?
        do l = 0,m3-1
            do k = 0,m3-1
                do j = 0,m2-1
                    do i=0,m1-1
                        ind = Elem%Idom(i,j,k)
                        acoeff = - mat%hprimez(k,l)*mat%GLLwx(i)*mat%GLLwy(j)*mat%GLLwz(l)*Elem%Jacob(i,j,l)
                        sum_vx = acoeff*Elem%InvGrad(0,2,i,j,l)*dom%Veloc(i,j,l,0,lnum)
                        sum_vy = acoeff*Elem%InvGrad(1,2,i,j,l)*dom%Veloc(i,j,l,1,lnum)
                        sum_vz = acoeff*Elem%InvGrad(2,2,i,j,l)*dom%Veloc(i,j,l,2,lnum)
                        ForcesFl(0,i,j,k) = ForcesFl(0,i,j,k) + sum_vx
                        ForcesFl(1,i,j,k) = ForcesFl(1,i,j,k) + sum_vy
                        ForcesFl(2,i,j,k) = ForcesFl(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        
        ! Assemblage
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    ind = Elem%Idom(i,j,k)
                    ! We should have atomic adds with openmp // here
                    champs1%fpml_Forces(ind,0) = champs1%fpml_Forces(ind,0) + ForcesFl(0,i,j,k)
                    champs1%fpml_Forces(ind,1) = champs1%fpml_Forces(ind,1) + ForcesFl(1,i,j,k)
                    champs1%fpml_Forces(ind,2) = champs1%fpml_Forces(ind,2) + ForcesFl(2,i,j,k)
                enddo
            enddo
        enddo

        return
    end subroutine forces_int_flu_pml

    subroutine pred_flu_pml(dom, mat, dt, champs1, Elem, lnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        type (subdomain), intent(IN) :: mat
        real, intent(in) :: dt
        type(champsfluidpml), intent(inout) :: champs1
        type(Element), intent(inout) :: Elem
        integer :: lnum
        !
        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVelPhi_dx, dVelPhi_dy, dVelPhi_dz
        integer :: m1, m2, m3
        integer :: i, j, k, ind
        real, dimension(:,:,:), allocatable :: VelPhi

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        allocate(VelPhi(0:m1-1,0:m2-1,0:m3-1))
        ! prediction in the element
        ! We do the sum V1+V2+V3 and F1+F2+F3 here
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    ind = Elem%Idom(i,j,k)
                    VelPhi(i,j,k) = champs1%fpml_VelPhi(ind,0) + &
                        champs1%fpml_VelPhi(ind,1) + &
                        champs1%fpml_VelPhi(ind,2)
                enddo
            enddo
        enddo
        ! XXX DumpS{xyz}(:,:,:,1) doit etre multiplie par 1/density
        ! d(rho*Phi)_d(xi,eta,zeta)
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,Elem%InvGrad, &
            VelPhi(:,:,:), dVelPhi_dx, dVelPhi_dy, dVelPhi_dz)

        ! prediction for (physical) velocity (which is the equivalent of a stress, here)
        ! V_x^x
        dom%Veloc(:,:,:,0,lnum) = Elem%xpml%DumpSx(:,:,:,0) * dom%Veloc(:,:,:,0,lnum) + &
            Elem%xpml%DumpSx(:,:,:,1) * Dt * dVelPhi_dx
        ! V_x^y = 0
        ! V_x^z = 0
        ! V_y^x = 0
        ! V_y^y
        dom%Veloc(:,:,:,1,lnum) = Elem%xpml%DumpSy(:,:,:,0) * dom%Veloc(:,:,:,1,lnum) + &
            Elem%xpml%DumpSy(:,:,:,1) * Dt * dVelPhi_dy
        ! V_y^z = 0
        ! V_z^x = 0
        ! V_z^y = 0
        ! V_z^z
        dom%Veloc(:,:,:,2,lnum) = Elem%xpml%DumpSz(:,:,:,0) * dom%Veloc(:,:,:,2,lnum) + &
            Elem%xpml%DumpSz(:,:,:,1) * Dt * dVelPhi_dz
        return
    end subroutine Pred_Flu_Pml

end module dom_fluidpml

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
