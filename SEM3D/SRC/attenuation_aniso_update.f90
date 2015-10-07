!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file attenuation_aniso_update.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine attenuation_aniso_update(epsilondev_xx_,epsilondev_yy_, &
    epsilondev_xy_,epsilondev_xz_,epsilondev_yz_, &
    epsilondev_xx_loc,epsilondev_yy_loc, &
    epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
    R_xx_,R_yy_,R_xy_,R_xz_,R_yz_, &
    factor_common_3_,alphaval_3_,betaval_3_,gammaval_3_, &
    mu_, ngllx,nglly,ngllz, n_solid)

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz, n_solid
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: mu_
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: epsilondev_xx_,epsilondev_yy_, &
        epsilondev_xy_,epsilondev_xz_,epsilondev_yz_, &
        epsilondev_xx_loc,epsilondev_yy_loc, &
        epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: R_xx_,R_yy_,R_xy_,R_xz_,R_yz_
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: factor_common_3_, &
        alphaval_3_,betaval_3_,gammaval_3_

    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: factor_loc,alphaval_loc,betaval_loc,gammaval_loc,Sn,Snp1
    integer :: i_sls


    do i_sls = 0,n_solid-1

        ! get coefficients for that standard linear solid
        factor_loc  (:,:,:) = mu_(:,:,:) *factor_common_3_(i_sls,:,:,:)
        alphaval_loc(:,:,:) = alphaval_3_(i_sls,:,:,:)
        betaval_loc (:,:,:) = betaval_3_(i_sls,:,:,:)
        gammaval_loc(:,:,:) = gammaval_3_(i_sls,:,:,:)

        ! term in xx
        Sn  (:,:,:)  = factor_loc(:,:,:) * epsilondev_xx_(:,:,:)
        Snp1(:,:,:)  = factor_loc(:,:,:) * epsilondev_xx_loc(:,:,:)
        R_xx_(i_sls,:,:,:) = alphaval_loc(:,:,:) * R_xx_(i_sls,:,:,:)  &
            + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)

        ! term in yy
        Sn  (:,:,:)  = factor_loc(:,:,:)  * epsilondev_yy_(:,:,:)
        Snp1(:,:,:)  = factor_loc(:,:,:)  * epsilondev_yy_loc(:,:,:)
        R_yy_(i_sls,:,:,:) = alphaval_loc(:,:,:) * R_yy_(i_sls,:,:,:)  &
            + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)

        ! term in xy
        Sn  (:,:,:)  = factor_loc(:,:,:)  * epsilondev_xy_(:,:,:)
        Snp1(:,:,:)  = factor_loc(:,:,:)  * epsilondev_xy_loc(:,:,:)
        R_xy_(i_sls,:,:,:) = alphaval_loc(:,:,:) * R_xy_(i_sls,:,:,:) &
            + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)

        ! term in xz
        Sn  (:,:,:)  = factor_loc(:,:,:) * epsilondev_xz_(:,:,:)
        Snp1(:,:,:)  = factor_loc(:,:,:) * epsilondev_xz_loc(:,:,:)
        R_xz_(i_sls,:,:,:) = alphaval_loc(:,:,:) * R_xz_(i_sls,:,:,:) &
            + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)

        ! term in yz
        Sn  (:,:,:)  = factor_loc(:,:,:) * epsilondev_yz_(:,:,:)
        Snp1(:,:,:)  = factor_loc(:,:,:) * epsilondev_yz_loc(:,:,:)
        R_yz_(i_sls,:,:,:) = alphaval_loc(:,:,:) * R_yz_(i_sls,:,:,:)  &
            + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)

    enddo   ! end of loop on memory variables

    ! a la fin on stocke la deformation deviatorique pour le prochain
    ! pas de temps de la mise a jour de Runge-Kutta faite ci-dessus

    ! save deviatoric strain for Runge-Kutta scheme

    epsilondev_xx_(:,:,:) = epsilondev_xx_loc(:,:,:)
    epsilondev_yy_(:,:,:) = epsilondev_yy_loc(:,:,:)
    epsilondev_xy_(:,:,:) = epsilondev_xy_loc(:,:,:)
    epsilondev_xz_(:,:,:) = epsilondev_xz_loc(:,:,:)
    epsilondev_yz_(:,:,:) = epsilondev_yz_loc(:,:,:)

end subroutine attenuation_aniso_update

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
