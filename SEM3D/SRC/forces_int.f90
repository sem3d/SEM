!>
!! \file forces_int.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<
module forces_aniso

    interface
       subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
         CHARACTER*1        TRANSA, TRANSB
         INTEGER            M, N, K, LDA, LDB, LDC
         DOUBLE PRECISION   ALPHA, BETA
         DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
       end subroutine DGEMM
    end interface

contains

    subroutine forces_int(Elem, htprimex, hprimey, htprimey, hprimez, htprimez,  &
               n_solid, aniso, solid)


        use sdomain

        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        integer, intent(IN) :: n_solid
        logical, intent(IN) :: aniso
        logical, intent(IN) :: solid   ! flag : solid or fluid element?

        integer :: n_z, m1,m2,m3, i,j,k
        real :: epsilon_trace_over_3
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dUx_dxi, dUx_deta, dUx_dzeta, &
            dUy_dxi, dUy_deta, dUy_dzeta, &
            dUz_dxi, dUz_deta, dUz_dzeta, &
            DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
            Fox,Foy,Foz,                         &
            dPhi_dxi, dPhi_deta, dPhi_dzeta, dPhiX,dPhiY,dPhiZ,Fo_Fl

        real, dimension(:,:,:), allocatable :: epsilondev_xx_loc, epsilondev_yy_loc, &
            epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
        real, dimension(:,:,:), allocatable :: epsilonvol_loc


        m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

        if(solid)then   ! SOLID PART OF THE DOMAIN

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1, Elem%Forces(:,:,:,0), m1, 0., dUx_dxi, m1)
        do n_z = 0,Elem%ngllz-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,0), m1, hprimey, m2, 0., dUx_deta(:,:,n_z), m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,0), m1*m2, hprimez, m3, 0., dUx_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1, Elem%Forces(:,:,:,1), m1, 0., dUy_dxi, m1)
        do n_z = 0,Elem%ngllz-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,1), m1, hprimey, m2, 0., dUy_deta(:,:,n_z), m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,1), m1*m2, hprimez, m3, 0., dUy_dzeta, m1*m2)

        call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1, Elem%Forces(:,:,:,2), m1, 0., dUz_dxi, m1)
        do n_z = 0,Elem%ngllz-1
            call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%Forces(:,:,n_z,2), m1, hprimey, m2, 0., dUz_deta(:,:,n_z), m1)
        enddo
        call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%Forces(:,:,:,2), m1*m2, hprimez, m3, 0., dUz_dzeta, m1*m2)

        do i = 0,m1-1
            do j = 0,m2-1
                do k = 0,m3-1

                    dxx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
                    dyy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)
                    dzz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

                    dyx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)
                    dzx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

                    dxy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
                    dzy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

                    dxz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
                    dyz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)

                enddo
            enddo
        enddo

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
                    Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Mu, Elem%Lambda, Elem%Cij, &
                    
                    m1,m2,m3, n_solid, &
                    Elem%onemSbeta, Elem%R_xx_, Elem%R_yy_, &
                    Elem%R_xy_, Elem%R_xz_, Elem%R_yz_)
                call attenuation_aniso_update(Elem%epsilondev_xx_,Elem%epsilondev_yy_, &
                    Elem%epsilondev_xy_,Elem%epsilondev_xz_,Elem%epsilondev_yz_, &
                    epsilondev_xx_loc,epsilondev_yy_loc, &
                    epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                    Elem%R_xx_,Elem%R_yy_,Elem%R_xy_,Elem%R_xz_,Elem%R_yz_, &
                    Elem%factor_common_3,Elem%alphaval_3,Elem%betaval_3,Elem%gammaval_3, &
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
                    Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Cij, &
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
                    Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Mu, Elem%Kappa, &
                    m1,m2,m3, n_solid, &
                                !     Elem%onemSbeta, &
                    Elem%R_xx_, Elem%R_yy_, &
                    Elem%R_xy_, Elem%R_xz_, Elem%R_yz_, &
                                !      Elem%onemPbeta, &
                    Elem%R_vol_)
                !                             Elem%Mu, Elem%Lambda, &
                !                             m1,m2,m3, n_solid, &
                !                             Elem%onemSbeta, Elem%R_xx_, Elem%R_yy_, &
                !                             Elem%R_xy_, Elem%R_xz_, Elem%R_yz_)
                call attenuation_update(Elem%epsilondev_xx_,Elem%epsilondev_yy_, &
                    Elem%epsilondev_xy_,Elem%epsilondev_xz_,Elem%epsilondev_yz_, &
                    epsilondev_xx_loc,epsilondev_yy_loc, &
                    epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                    Elem%R_xx_,Elem%R_yy_,Elem%R_xy_,Elem%R_xz_,Elem%R_yz_, &
                    Elem%factor_common_3,Elem%alphaval_3,Elem%betaval_3,Elem%gammaval_3, &
                    Elem%Mu, m1,m2,m3, n_solid, &
                    Elem%Kappa, Elem%epsilonvol_,epsilonvol_loc,Elem%R_vol_, &
                    Elem%factor_common_P,Elem%alphaval_P,Elem%betaval_P,Elem%gammaval_P)
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
                    Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Mu, Elem%Lambda, &
                    m1,m2,m3)
            endif
        endif

        Elem%Forces(:,:,:,0) = -Fox
        Elem%Forces(:,:,:,1) = -Foy
        Elem%Forces(:,:,:,2) = -Foz

   !---------------------------------
        else      ! FLUID PART OF THE DOMAIN

        !- gradients at GLLs points
        ! d(rho*Phi)_dxi
            call DGEMM ('N', 'N', m1, m2*m3, m1, 1., htprimex, m1, Elem%ForcesFl(:,:,:), m1, 0., dPhi_dxi, m1)
        ! d(rho*Phi)_deta
            do n_z = 0,Elem%ngllz-1
                call DGEMM ('N', 'N', m1, m2, m2, 1., Elem%ForcesFl(:,:,n_z), m1, hprimey, m2, 0., dPhi_deta(:,:,n_z), m1)
            enddo
        ! d(rho*Phi)_dzeta
            call DGEMM ('N', 'N', m1*m2, m3, m3, 1., Elem%ForcesFl(:,:,:), m1*m2, hprimez, m3, 0., dPhi_dzeta, m1*m2)

        ! d(rho*Phi)_dX 
            dPhiX(:,:,:) = dPhi_dxi(:,:,:)*Elem%InvGrad(:,:,:,0,0) + dPhi_deta(:,:,:)*Elem%InvGrad(:,:,:,0,1) +   &
                           dPhi_dzeta(:,:,:)*Elem%InvGrad(:,:,:,0,2)
        ! d(rho*Phi)_dY 
            dPhiY(:,:,:) = dPhi_dxi(:,:,:)*Elem%InvGrad(:,:,:,1,0) + dPhi_deta(:,:,:)*Elem%InvGrad(:,:,:,1,1) +   &
                           dPhi_dzeta(:,:,:)*Elem%InvGrad(:,:,:,1,2)
        ! d(rho*Phi)_dZ 
            dPhiZ(:,:,:) = dPhi_dxi(:,:,:)*Elem%InvGrad(:,:,:,2,0) + dPhi_deta(:,:,:)*Elem%InvGrad(:,:,:,2,1) +   &
                           dPhi_dzeta(:,:,:)*Elem%InvGrad(:,:,:,2,2)

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
                         Elem%Jacob,Elem%wgtx,Elem%wgty,Elem%wgtz, &
                         dPhiX,dPhiY,dPhiZ,       &
                         Elem%Density,            &
                         m1,m2,m3)

            Elem%ForcesFl(:,:,:) = -Fo_Fl(:,:,:)

           
        end if

        return
    end subroutine forces_int

end module forces_aniso
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
