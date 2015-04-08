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

    !    partie isotrope
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: kappa_
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: epsilonvol_,epsilonvol_loc
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: R_vol_
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: factor_common_P_, &
        alphaval_P_,betaval_P_,gammaval_P_

    integer :: i_sls,i,j,k
    real :: factorS_loc,alphavalS_loc,betavalS_loc,gammavalS_loc,Sn,Snp1
    real :: factorP_loc,alphavalP_loc,betavalP_loc,gammavalP_loc,Pn,Pnp1

    do k = 0, ngllz - 1
        do j = 0, nglly - 1
            do i = 0, ngllx-1
                do i_sls = 0,n_solid-1

                    ! get coefficients for that standard linear solid
                    factorS_loc   = mu_(i,j,k) *factor_common_3_(i_sls,i,j,k)
                    alphavalS_loc = alphaval_3_(i_sls,i,j,k)
                    betavalS_loc  = betaval_3_(i_sls,i,j,k)
                    gammavalS_loc = gammaval_3_(i_sls,i,j,k)

                    ! get coefficients for that standard linear solid
                    factorP_loc   = kappa_(i,j,k) *factor_common_P_(i_sls,i,j,k)
                    alphavalP_loc = alphaval_P_(i_sls,i,j,k)
                    betavalP_loc  = betaval_P_(i_sls,i,j,k)
                    gammavalP_loc = gammaval_P_(i_sls,i,j,k)

                    !  terme volumique
                    Pn   = factorP_loc * epsilonvol_(i,j,k)
                    Pnp1 = factorP_loc * epsilonvol_loc(i,j,k)
                    R_vol_(i_sls,i,j,k) = alphavalP_loc * R_vol_(i_sls,i,j,k)  &
                        + betavalP_loc*Pn+gammavalP_loc*Pnp1


                    !     termes deviatoires
                    ! term in xx
                    Sn    = factorS_loc * epsilondev_xx_(i,j,k)
                    Snp1  = factorS_loc * epsilondev_xx_loc(i,j,k)
                    R_xx_(i_sls,i,j,k) = alphavalS_loc * R_xx_(i_sls,i,j,k)  &
                        + betavalS_loc*Sn+gammavalS_loc*Snp1

                    ! term in yy
                    Sn    = factorS_loc  * epsilondev_yy_(i,j,k)
                    Snp1  = factorS_loc  * epsilondev_yy_loc(i,j,k)
                    R_yy_(i_sls,i,j,k) = alphavalS_loc * R_yy_(i_sls,i,j,k)  &
                        + betavalS_loc*Sn+gammavalS_loc*Snp1

                    ! term in xy
                    Sn    = factorS_loc  * epsilondev_xy_(i,j,k)
                    Snp1  = factorS_loc  * epsilondev_xy_loc(i,j,k)
                    R_xy_(i_sls,i,j,k) = alphavalS_loc * R_xy_(i_sls,i,j,k) &
                        + betavalS_loc*Sn+gammavalS_loc*Snp1

                    ! term in xz
                    Sn    = factorS_loc * epsilondev_xz_(i,j,k)
                    Snp1  = factorS_loc * epsilondev_xz_loc(i,j,k)
                    R_xz_(i_sls,i,j,k) = alphavalS_loc * R_xz_(i_sls,i,j,k) &
                        + betavalS_loc*Sn+gammavalS_loc*Snp1

                    ! term in yz
                    Sn    = factorS_loc * epsilondev_yz_(i,j,k)
                    Snp1  = factorS_loc * epsilondev_yz_loc(i,j,k)
                    R_yz_(i_sls,i,j,k) = alphavalS_loc * R_yz_(i_sls,i,j,k)  &
                        + betavalS_loc*Sn+gammavalS_loc*Snp1

                enddo   ! end of loop on memory variables

                ! a la fin on stocke la deformation deviatorique pour le prochain
                ! pas de temps de la mise a jour de Runge-Kutta faite ci-dessus

                ! save deviatoric strain for Runge-Kutta scheme

                !     partie deviatoire
                epsilondev_xx_(i,j,k) = epsilondev_xx_loc(i,j,k)
                epsilondev_yy_(i,j,k) = epsilondev_yy_loc(i,j,k)
                epsilondev_xy_(i,j,k) = epsilondev_xy_loc(i,j,k)
                epsilondev_xz_(i,j,k) = epsilondev_xz_loc(i,j,k)
                epsilondev_yz_(i,j,k) = epsilondev_yz_loc(i,j,k)
                !   partie isotrope
                epsilonvol_(i,j,k) = epsilonvol_loc(i,j,k)

            end do
        end do
    end do

end subroutine attenuation_update
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
