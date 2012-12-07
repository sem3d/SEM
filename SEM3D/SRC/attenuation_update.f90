!>
!! \file attenuation_update.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine attenuation_update(epsilondev_xx_,epsilondev_yy_, &
    epsilondev_xy_,epsilondev_xz_,epsilondev_yz_, &
    epsilondev_xx_loc,epsilondev_yy_loc, &
    epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
    R_xx_,R_yy_,R_xy_,R_xz_,R_yz_, &
    factor_common_3_,alphaval_3_,betaval_3_,gammaval_3_, &
    mu_, ngllx,nglly,ngllz, n_solid, &
    kappa_,epsilonvol_,epsilonvol_loc,R_vol_, &
    factor_common_P_,alphaval_P_,betaval_P_,gammaval_P_)

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz, n_solid
    !   partie deviatoire
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: mu_
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: epsilondev_xx_,epsilondev_yy_, &
        epsilondev_xy_,epsilondev_xz_,epsilondev_yz_, &
        epsilondev_xx_loc,epsilondev_yy_loc, &
        epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: R_xx_,R_yy_,R_xy_,R_xz_,R_yz_
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: factor_common_3_, &
        alphaval_3_,betaval_3_,gammaval_3_

    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: factorS_loc,alphavalS_loc,betavalS_loc,gammavalS_loc,Sn,Snp1
    integer :: i_sls
    !    partie isotrope
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: kappa_
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: epsilonvol_,epsilonvol_loc
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: R_vol_
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: factor_common_P_, &
        alphaval_P_,betaval_P_,gammaval_P_

    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: factorP_loc,alphavalP_loc,betavalP_loc,gammavalP_loc,Pn,Pnp1


    do i_sls = 0,n_solid-1

        ! get coefficients for that standard linear solid
        factorS_loc  (:,:,:) = mu_(:,:,:) *factor_common_3_(i_sls,:,:,:)
        alphavalS_loc(:,:,:) = alphaval_3_(i_sls,:,:,:)
        betavalS_loc (:,:,:) = betaval_3_(i_sls,:,:,:)
        gammavalS_loc(:,:,:) = gammaval_3_(i_sls,:,:,:)

        ! get coefficients for that standard linear solid
        factorP_loc  (:,:,:) = kappa_(:,:,:) *factor_common_P_(i_sls,:,:,:)
        alphavalP_loc(:,:,:) = alphaval_P_(i_sls,:,:,:)
        betavalP_loc (:,:,:) = betaval_P_(i_sls,:,:,:)
        gammavalP_loc(:,:,:) = gammaval_P_(i_sls,:,:,:)

        !  terme volumique
        Pn  (:,:,:)  = factorP_loc(:,:,:) * epsilonvol_(:,:,:)
        Pnp1(:,:,:)  = factorP_loc(:,:,:) * epsilonvol_loc(:,:,:)
        R_vol_(i_sls,:,:,:) = alphavalP_loc(:,:,:) * R_vol_(i_sls,:,:,:)  &
            + betavalP_loc(:,:,:)*Pn(:,:,:)+gammavalP_loc(:,:,:)*Pnp1(:,:,:)


        !     termes deviatoires
        ! term in xx
        Sn  (:,:,:)  = factorS_loc(:,:,:) * epsilondev_xx_(:,:,:)
        Snp1(:,:,:)  = factorS_loc(:,:,:) * epsilondev_xx_loc(:,:,:)
        R_xx_(i_sls,:,:,:) = alphavalS_loc(:,:,:) * R_xx_(i_sls,:,:,:)  &
            + betavalS_loc(:,:,:)*Sn(:,:,:)+gammavalS_loc(:,:,:)*Snp1(:,:,:)

        ! term in yy
        Sn  (:,:,:)  = factorS_loc(:,:,:)  * epsilondev_yy_(:,:,:)
        Snp1(:,:,:)  = factorS_loc(:,:,:)  * epsilondev_yy_loc(:,:,:)
        R_yy_(i_sls,:,:,:) = alphavalS_loc(:,:,:) * R_yy_(i_sls,:,:,:)  &
            + betavalS_loc(:,:,:)*Sn(:,:,:)+gammavalS_loc(:,:,:)*Snp1(:,:,:)

        ! term in xy
        Sn  (:,:,:)  = factorS_loc(:,:,:)  * epsilondev_xy_(:,:,:)
        Snp1(:,:,:)  = factorS_loc(:,:,:)  * epsilondev_xy_loc(:,:,:)
        R_xy_(i_sls,:,:,:) = alphavalS_loc(:,:,:) * R_xy_(i_sls,:,:,:) &
            + betavalS_loc(:,:,:)*Sn(:,:,:)+gammavalS_loc(:,:,:)*Snp1(:,:,:)

        ! term in xz
        Sn  (:,:,:)  = factorS_loc(:,:,:) * epsilondev_xz_(:,:,:)
        Snp1(:,:,:)  = factorS_loc(:,:,:) * epsilondev_xz_loc(:,:,:)
        R_xz_(i_sls,:,:,:) = alphavalS_loc(:,:,:) * R_xz_(i_sls,:,:,:) &
            + betavalS_loc(:,:,:)*Sn(:,:,:)+gammavalS_loc(:,:,:)*Snp1(:,:,:)

        ! term in yz
        Sn  (:,:,:)  = factorS_loc(:,:,:) * epsilondev_yz_(:,:,:)
        Snp1(:,:,:)  = factorS_loc(:,:,:) * epsilondev_yz_loc(:,:,:)
        R_yz_(i_sls,:,:,:) = alphavalS_loc(:,:,:) * R_yz_(i_sls,:,:,:)  &
            + betavalS_loc(:,:,:)*Sn(:,:,:)+gammavalS_loc(:,:,:)*Snp1(:,:,:)

    enddo   ! end of loop on memory variables

    ! a la fin on stocke la deformation deviatorique pour le prochain
    ! pas de temps de la mise a jour de Runge-Kutta faite ci-dessus

    ! save deviatoric strain for Runge-Kutta scheme

    !     partie deviatoire
    epsilondev_xx_(:,:,:) = epsilondev_xx_loc(:,:,:)
    epsilondev_yy_(:,:,:) = epsilondev_yy_loc(:,:,:)
    epsilondev_xy_(:,:,:) = epsilondev_xy_loc(:,:,:)
    epsilondev_xz_(:,:,:) = epsilondev_xz_loc(:,:,:)
    epsilondev_yz_(:,:,:) = epsilondev_yz_loc(:,:,:)
    !   partie isotrope
    epsilonvol_(:,:,:) = epsilonvol_loc(:,:,:)

end subroutine attenuation_update
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
