!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file forces_int.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<
module forces_aniso
    use deriv3d
    use sdomain
    use schamps
    implicit none

contains

    subroutine forces_int(Elem, mat, htprimex, hprimey, htprimey, hprimez, htprimez,  &
               n_solid, aniso, solid, champs1)

        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        integer, intent(IN) :: n_solid
        logical, intent(IN) :: aniso
        logical, intent(IN) :: solid   ! flag : solid or fluid element?
        type(champs), intent(inout) :: champs1

        integer :: m1,m2,m3, i,j,k,i_dir
        real :: epsilon_trace_over_3
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) ::  DXX,DXY,DXZ, &
            DYX,DYY,DYZ, &
            DZX,DZY,DZZ, &
            Fox,Foy,Foz, &
            dPhiX,dPhiY,dPhiZ,Fo_Fl

        real, dimension(:,:,:), allocatable :: epsilondev_xx_loc, epsilondev_yy_loc, &
            epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
        real, dimension(:,:,:), allocatable :: epsilonvol_loc
        real, dimension(:,:,:,:), allocatable :: Depla
        real, dimension(:,:,:), allocatable :: Phi


        m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

        if(solid)then   ! SOLID PART OF THE DOMAIN
            allocate(Depla(0:m1-1,0:m2-1,0:m3-1,0:2))
            do i_dir = 0,2
                do k = 0,m3-1
                    do j = 0,m2-1
                        do i = 0,m1-1
                            Depla(i,j,k,i_dir) = champs1%Depla(Elem%Isol(i,j,k),i_dir)
                        enddo
                    enddo
                enddo
            enddo
            
            call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Depla(:,:,:,0),dxx,dyx,dzx)
            call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Depla(:,:,:,1),dxy,dyy,dzy)
            call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Depla(:,:,:,2),dxz,dyz,dzz)
            deallocate(Depla)

            if (n_solid>0) then
                if (aniso) then
                else
                    allocate (epsilonvol_loc(0:m1-1,0:m2-1,0:m3-1))
                endif
                allocate (epsilondev_xx_loc(0:m1-1,0:m2-1,0:m3-1))
                allocate (epsilondev_yy_loc(0:m1-1,0:m2-1,0:m3-1))
                allocate (epsilondev_xy_loc(0:m1-1,0:m2-1,0:m3-1))
                allocate (epsilondev_xz_loc(0:m1-1,0:m2-1,0:m3-1))
                allocate (epsilondev_yz_loc(0:m1-1,0:m2-1,0:m3-1))
                do i = 0,m1-1
                    do j = 0,m2-1
                        do k = 0,m3-1
                            epsilon_trace_over_3 = 0.333333333333333333333333333333d0 * (DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k))
                            if (aniso) then
                            else
                                epsilonvol_loc(i,j,k) = DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k)
                            endif
                            epsilondev_xx_loc(i,j,k) = DXX(i,j,k) - epsilon_trace_over_3
                            epsilondev_yy_loc(i,j,k) = DYY(i,j,k) - epsilon_trace_over_3
                            epsilondev_xy_loc(i,j,k) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                            epsilondev_xz_loc(i,j,k) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                            epsilondev_yz_loc(i,j,k) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                        enddo
                    enddo
                enddo
            endif

            if (aniso) then
                if (n_solid>0) then
                    call calcul_forces_aniso_att(Fox,Foy,Foz, &
                        Elem%Invgrad, &
                        htprimex, htprimey, htprimez, &
                        Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                        DXX,DXY,DXZ, &
                        DYX,DYY,DYZ, &
                        DZX,DZY,DZZ, &
                        Elem%Mu, Elem%Lambda, Elem%sl%Cij, &
                        
                        m1,m2,m3, n_solid, &
                        Elem%sl%onemSbeta, Elem%sl%R_xx_, Elem%sl%R_yy_, &
                        Elem%sl%R_xy_, Elem%sl%R_xz_, Elem%sl%R_yz_)
                    call attenuation_aniso_update(Elem%sl%epsilondev_xx_,Elem%sl%epsilondev_yy_, &
                        Elem%sl%epsilondev_xy_,Elem%sl%epsilondev_xz_,Elem%sl%epsilondev_yz_, &
                        epsilondev_xx_loc,epsilondev_yy_loc, &
                        epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                        Elem%sl%R_xx_,Elem%sl%R_yy_,Elem%sl%R_xy_,Elem%sl%R_xz_,Elem%sl%R_yz_, &
                        Elem%sl%factor_common_3,Elem%sl%alphaval_3,Elem%sl%betaval_3,Elem%sl%gammaval_3, &
                        Elem%Mu, m1,m2,m3, n_solid)
                    deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
                    !      deallocate(epsilonvol_loc)
                else
                    call calcul_forces_aniso(Fox,Foy,Foz,  &
                        Elem%Invgrad, &
                        htprimex, htprimey, htprimez, &
                        Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                        DXX,DXY,DXZ, &
                        DYX,DYY,DYZ, &
                        DZX,DZY,DZZ, &
                        Elem%sl%Cij, &
                        m1,m2,m3)
                endif
            else
                if (n_solid>0) then
                    call calcul_forces_att(Fox,Foy,Foz, &
                        Elem%Invgrad, &
                        htprimex, htprimey, htprimez, &
                        Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                        DXX,DXY,DXZ, &
                        DYX,DYY,DYZ, &
                        DZX,DZY,DZZ, &
                        Elem%Mu, Elem%Kappa, &
                        m1,m2,m3, n_solid, &
                        Elem%sl%R_xx_, Elem%sl%R_yy_, &
                        Elem%sl%R_xy_, Elem%sl%R_xz_, Elem%sl%R_yz_, &
                        Elem%sl%R_vol_, &
                        Elem%sl%onemSbeta, &
                        Elem%sl%onemPbeta &
                        )
                    !                             Elem%Mu, Elem%Lambda, &
                    !                             m1,m2,m3, n_solid, &
                    !                             Elem%onemSbeta, Elem%R_xx_, Elem%R_yy_, &
                    !                             Elem%R_xy_, Elem%R_xz_, Elem%R_yz_)
                    call attenuation_update(Elem%sl%epsilondev_xx_,Elem%sl%epsilondev_yy_, &
                        Elem%sl%epsilondev_xy_,Elem%sl%epsilondev_xz_,Elem%sl%epsilondev_yz_, &
                        epsilondev_xx_loc,epsilondev_yy_loc, &
                        epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                        Elem%sl%R_xx_,Elem%sl%R_yy_,Elem%sl%R_xy_,Elem%sl%R_xz_,Elem%sl%R_yz_, &
                        Elem%sl%factor_common_3,Elem%sl%alphaval_3,Elem%sl%betaval_3,Elem%sl%gammaval_3, &
                        Elem%Mu, m1,m2,m3, n_solid, &
                        Elem%Kappa, Elem%sl%epsilonvol_,epsilonvol_loc,Elem%sl%R_vol_, &
                        Elem%sl%factor_common_P,Elem%sl%alphaval_P,Elem%sl%betaval_P,Elem%sl%gammaval_P)
                    !                              kappa_,elsilonvol_,epsilon_vol_loc,R_vol_, &
                    !                              factor_common_P_,alphaval_P_,betaval_P_,gammaval_P_)
                    deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
                    deallocate(epsilonvol_loc)
                else
                    call calcul_forces(Fox,Foy,Foz,  &
                        Elem%Invgrad, &
                        htprimex, htprimey, htprimez, &
                        Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                        DXX,DXY,DXZ, &
                        DYX,DYY,DYZ, &
                        DZX,DZY,DZZ, &
                        Elem%Mu, Elem%Lambda, &
                        m1,m2,m3)
                endif
            endif

            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                        champs1%Forces(Elem%Isol(i,j,k),0) = champs1%Forces(Elem%Isol(i,j,k),0)-Fox(i,j,k)
                        champs1%Forces(Elem%Isol(i,j,k),1) = champs1%Forces(Elem%Isol(i,j,k),1)-Foy(i,j,k)
                        champs1%Forces(Elem%Isol(i,j,k),2) = champs1%Forces(Elem%Isol(i,j,k),2)-Foz(i,j,k)
                    enddo
                enddo
            enddo
   !---------------------------------
        else      ! FLUID PART OF THE DOMAIN
            ! d(rho*Phi)_dX
            ! d(rho*Phi)_dY
            ! d(rho*Phi)_dZ
            allocate(Phi(0:m1-1,0:m2-1,0:m3-1))
            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                        Phi(i,j,k) = champs1%Phi(Elem%Iflu(i,j,k))
                    enddo
                enddo
            enddo
            call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad, Phi, &
                dPhiX, dPhiY, dPhiZ)
            deallocate(Phi)

            ! internal forces
            call calcul_forces_fluid(Fo_Fl,                &
                         Elem%Invgrad, &
                         htprimex,htprimey,htprimez, &
                         Elem%Jacob,mat%GLLwx,mat%GLLwy,mat%GLLwz, &
                         dPhiX,dPhiY,dPhiZ,       &
                         Elem%Density,            &
                         m1,m2,m3)

            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                       champs1%ForcesFl(Elem%Iflu(i,j,k)) = champs1%ForcesFl(Elem%Iflu(i,j,k))-Fo_Fl(i,j,k)
                    enddo
                enddo
            enddo

        end if

        return
    end subroutine forces_int

    subroutine forces_int_flu_pml(Elem, mat, champs1)
        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        type(champs), intent(inout) :: champs1
        !
        integer :: m1, m2, m3
        integer :: i, j, k, l, ind
        real :: sum_vx, sum_vy, sum_vz, acoeff
        real, dimension(0:2,0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1)  :: ForcesFl

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        do k = 0,m3-1
            do j = 0,m2-1
                do i=0,m1-1
                    ind = Elem%flpml%IFluPml(i,j,k)
                    sum_vx = 0d0
                    sum_vy = 0d0
                    sum_vz = 0d0
                    do l = 0,m1-1
                        acoeff = - mat%hprimex(i,l)*mat%GLLwx(l)*mat%GLLwy(j)*mat%GLLwz(k)*Elem%Jacob(l,j,k)
                        sum_vx = sum_vx + acoeff*Elem%InvGrad(l,j,k,0,0)*Elem%flpml%Veloc(l,j,k,0)
                        sum_vy = sum_vy + acoeff*Elem%InvGrad(l,j,k,1,0)*Elem%flpml%Veloc(l,j,k,1)
                        sum_vz = sum_vz + acoeff*Elem%InvGrad(l,j,k,2,0)*Elem%flpml%Veloc(l,j,k,2)
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
                        ind = Elem%flpml%IFluPml(i,j,k)
                        acoeff = - mat%hprimey(j,l)*mat%GLLwx(i)*mat%GLLwy(l)*mat%GLLwz(k)*Elem%Jacob(i,l,k)
                        sum_vx = acoeff*Elem%InvGrad(i,l,k,0,1)*Elem%flpml%Veloc(i,l,k,0)
                        sum_vy = acoeff*Elem%InvGrad(i,l,k,1,1)*Elem%flpml%Veloc(i,l,k,1)
                        sum_vz = acoeff*Elem%InvGrad(i,l,k,2,1)*Elem%flpml%Veloc(i,l,k,2)
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
                        ind = Elem%flpml%IFluPml(i,j,k)
                        acoeff = - mat%hprimez(k,l)*mat%GLLwx(i)*mat%GLLwy(j)*mat%GLLwz(l)*Elem%Jacob(i,j,l)
                        sum_vx = acoeff*Elem%InvGrad(i,j,l,0,2)*Elem%flpml%Veloc(i,j,l,0)
                        sum_vy = acoeff*Elem%InvGrad(i,j,l,1,2)*Elem%flpml%Veloc(i,j,l,1)
                        sum_vz = acoeff*Elem%InvGrad(i,j,l,2,2)*Elem%flpml%Veloc(i,j,l,2)
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
                    ind = Elem%flpml%IFluPml(i,j,k)
                    ! We should have atomic adds with openmp // here
                    champs1%fpml_Forces(ind+0) = champs1%fpml_Forces(ind+0) + ForcesFl(0,i,j,k)
                    champs1%fpml_Forces(ind+1) = champs1%fpml_Forces(ind+1) + ForcesFl(1,i,j,k)
                    champs1%fpml_Forces(ind+2) = champs1%fpml_Forces(ind+2) + ForcesFl(2,i,j,k)
                enddo
            enddo
        enddo

        return
    end subroutine forces_int_flu_pml


    subroutine pred_flu_pml(Elem, bega, dt, mat, champs0, champs1)
        ! same as previously, but for fluid part
        implicit none

        type(Element), intent(inout) :: Elem
        type (subdomain), intent(IN) :: mat
        type(champs), intent(in) :: champs0
        type(champs), intent(inout) :: champs1
        real, intent(in) :: bega, dt
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
                    ind = Elem%flpml%IFluPml(i,j,k)
                    VelPhi(i,j,k) = champs0%fpml_VelPhi(ind)+ &
                        champs0%fpml_VelPhi(ind+1)+champs0%fpml_VelPhi(ind+2) &
                        + dt*(0.5-bega)*(champs1%fpml_Forces(ind)+&
                        champs1%fpml_Forces(ind+1)+champs1%fpml_Forces(ind+2))
                enddo
            enddo
        enddo
        ! XXX DumpS{xyz}(:,:,:,1) doit etre multiplie par 1/density

        ! d(rho*Phi)_d(xi,eta,zeta)
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,Elem%InvGrad, &
            VelPhi(:,:,:), dVelPhi_dx, dVelPhi_dy, dVelPhi_dz)

        ! prediction for (physical) velocity (which is the equivalent of a stress, here)
        ! V_x^x
        Elem%flpml%Veloc(:,:,:,0) = Elem%xpml%DumpSx(:,:,:,0) * Elem%flpml%Veloc(:,:,:,0) + &
            Elem%xpml%DumpSx(:,:,:,1) * Dt * dVelPhi_dx
        ! V_x^y = 0
        ! V_x^z = 0
        ! V_y^x = 0
        ! V_y^y
        Elem%flpml%Veloc(:,:,:,1) = Elem%xpml%DumpSy(:,:,:,0) * Elem%flpml%Veloc(:,:,:,1) + &
            Elem%xpml%DumpSy(:,:,:,1) * Dt * dVelPhi_dy
        ! V_y^z = 0
        ! V_z^x = 0
        ! V_z^y = 0
        ! V_z^z
        Elem%flpml%Veloc(:,:,:,2) = Elem%xpml%DumpSz(:,:,:,0) * Elem%flpml%Veloc(:,:,:,2) + &
            Elem%xpml%DumpSz(:,:,:,1) * Dt * dVelPhi_dz

        return
    end subroutine Pred_Flu_Pml


end module forces_aniso

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
