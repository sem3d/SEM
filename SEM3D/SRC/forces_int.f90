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

contains

    subroutine forces_int(Elem, mat, htprimex, hprimey, htprimey, hprimez, htprimez,  &
               n_solid, aniso, solid)


        use sdomain

        implicit none

        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        integer, intent(IN) :: n_solid
        logical, intent(IN) :: aniso
        logical, intent(IN) :: solid   ! flag : solid or fluid element?

        integer :: m1,m2,m3, i,j,k
        real :: epsilon_trace_over_3
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) ::  DXX,DXY,DXZ, &
            DYX,DYY,DYZ, &
            DZX,DZY,DZZ, &
            Fox,Foy,Foz, &
            dPhiX,dPhiY,dPhiZ,Fo_Fl

        real, dimension(:,:,:), allocatable :: epsilondev_xx_loc, epsilondev_yy_loc, &
            epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
        real, dimension(:,:,:), allocatable :: epsilonvol_loc


        m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

        if(solid)then   ! SOLID PART OF THE DOMAIN

        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Elem%sl%Forces(:,:,:,0),dxx,dyx,dzx)
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Elem%sl%Forces(:,:,:,1),dxy,dyy,dzy)
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Elem%sl%Forces(:,:,:,2),dxz,dyz,dzz)

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
                    Elem%Invgrad(:,:,:,0,0), &
                    Elem%Invgrad(:,:,:,1,0), &
                    Elem%Invgrad(:,:,:,2,0), &
                    Elem%Invgrad(:,:,:,0,1), &
                    Elem%Invgrad(:,:,:,1,1), &
                    Elem%Invgrad(:,:,:,2,1), &
                    Elem%Invgrad(:,:,:,0,2), &
                    Elem%Invgrad(:,:,:,1,2), &
                    Elem%Invgrad(:,:,:,2,2), &
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
                    Elem%Invgrad(:,:,:,0,0), &
                    Elem%Invgrad(:,:,:,1,0), &
                    Elem%Invgrad(:,:,:,2,0), &
                    Elem%Invgrad(:,:,:,0,1), &
                    Elem%Invgrad(:,:,:,1,1), &
                    Elem%Invgrad(:,:,:,2,1), &
                    Elem%Invgrad(:,:,:,0,2), &
                    Elem%Invgrad(:,:,:,1,2), &
                    Elem%Invgrad(:,:,:,2,2), &
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
                    Elem%Invgrad(:,:,:,0,0), &
                    Elem%Invgrad(:,:,:,1,0), &
                    Elem%Invgrad(:,:,:,2,0), &
                    Elem%Invgrad(:,:,:,0,1), &
                    Elem%Invgrad(:,:,:,1,1), &
                    Elem%Invgrad(:,:,:,2,1), &
                    Elem%Invgrad(:,:,:,0,2), &
                    Elem%Invgrad(:,:,:,1,2), &
                    Elem%Invgrad(:,:,:,2,2), &
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
                    Elem%Invgrad(:,:,:,0,0), &
                    Elem%Invgrad(:,:,:,1,0), &
                    Elem%Invgrad(:,:,:,2,0), &
                    Elem%Invgrad(:,:,:,0,1), &
                    Elem%Invgrad(:,:,:,1,1), &
                    Elem%Invgrad(:,:,:,2,1), &
                    Elem%Invgrad(:,:,:,0,2), &
                    Elem%Invgrad(:,:,:,1,2), &
                    Elem%Invgrad(:,:,:,2,2), &
                    htprimex, htprimey, htprimez, &
                    Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Mu, Elem%Lambda, &
                    m1,m2,m3)
            endif
        endif

        Elem%sl%Forces(:,:,:,0) = -Fox
        Elem%sl%Forces(:,:,:,1) = -Foy
        Elem%sl%Forces(:,:,:,2) = -Foz

   !---------------------------------
        else      ! FLUID PART OF THE DOMAIN
            ! d(rho*Phi)_dX
            ! d(rho*Phi)_dY
            ! d(rho*Phi)_dZ
            call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Elem%fl%ForcesFl, &
                dPhiX, dPhiY, dPhiZ)

            ! internal forces
            call calcul_forces_fluid(Fo_Fl,                &
                         Elem%Invgrad(:,:,:,0,0), &
                         Elem%Invgrad(:,:,:,1,0), &
                         Elem%Invgrad(:,:,:,2,0), &
                         Elem%Invgrad(:,:,:,0,1), &
                         Elem%Invgrad(:,:,:,1,1), &
                         Elem%Invgrad(:,:,:,2,1), &
                         Elem%Invgrad(:,:,:,0,2), &
                         Elem%Invgrad(:,:,:,1,2), &
                         Elem%Invgrad(:,:,:,2,2), &
                         htprimex,htprimey,htprimez, &
                         Elem%Jacob,mat%GLLwx,mat%GLLwy,mat%GLLwz, &
                         dPhiX,dPhiY,dPhiZ,       &
                         Elem%Density,            &
                         m1,m2,m3)

            Elem%fl%ForcesFl(:,:,:) = -Fo_Fl(:,:,:)

        end if

        return
    end subroutine forces_int

end module forces_aniso
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
