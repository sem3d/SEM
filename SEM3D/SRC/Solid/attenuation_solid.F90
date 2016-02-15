!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file attenuation_update.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

module attenuation_solid

    implicit none

contains

  subroutine attenuation_update(dom,lnum,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,ngll,n_solid,aniso)

        use sdomain
        implicit none
#include "index.h"

        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: ngll,n_solid
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
        logical aniso

        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: epsilondev_xx_loc,epsilondev_yy_loc
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: epsilonvol_loc

        integer :: i_sls,i,j,k
        real :: factorS_loc,alphavalS_loc,betavalS_loc,gammavalS_loc,Sn,Snp1
        real :: factorP_loc,alphavalP_loc,betavalP_loc,gammavalP_loc,Pn,Pnp1
        real epsilon_trace_over_3

        if (n_solid>0) then
            do i = 0,ngll-1
                do j = 0,ngll-1
                    do k = 0,ngll-1
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

        do k = 0, ngll - 1
            do j = 0, ngll - 1
                do i = 0, ngll-1
                    do i_sls = 0,n_solid-1
    
                        ! get coefficients for that standard linear solid
                        factorS_loc   = dom%Mu_(i,j,k,lnum) * dom%factor_common_3(i_sls,i,j,k,lnum)
                        alphavalS_loc = dom%alphaval_3(i_sls,i,j,k,lnum)
                        betavalS_loc  = dom%betaval_3 (i_sls,i,j,k,lnum)
                        gammavalS_loc = dom%gammaval_3(i_sls,i,j,k,lnum)
    
                        ! get coefficients for that standard linear solid
                        factorP_loc   = dom%Kappa_(i,j,k,lnum) * dom%factor_common_P(i_sls,i,j,k,lnum)
                        alphavalP_loc = dom%alphaval_P(i_sls,i,j,k,lnum)
                        betavalP_loc  = dom%betaval_P (i_sls,i,j,k,lnum)
                        gammavalP_loc = dom%gammaval_P(i_sls,i,j,k,lnum)
    
                        !  terme volumique
                        Pn   = factorP_loc * dom%epsilonvol_(i,j,k,lnum)
                        Pnp1 = factorP_loc * epsilonvol_loc(i,j,k)
                        dom%R_vol_(i_sls,i,j,k,lnum) = alphavalP_loc * dom%R_vol_(i_sls,i,j,k,lnum)  &
                            + betavalP_loc*Pn+gammavalP_loc*Pnp1
    
                        !     termes deviatoires
                        ! term in xx
                        Sn    = factorS_loc * dom%epsilondev_xx_(i,j,k,lnum)
                        Snp1  = factorS_loc * epsilondev_xx_loc(i,j,k)
                        dom%R_xx_(i_sls,i,j,k,lnum) = alphavalS_loc * dom%R_xx_(i_sls,i,j,k,lnum)  &
                            + betavalS_loc*Sn+gammavalS_loc*Snp1
    
                        ! term in yy
                        Sn    = factorS_loc  * dom%epsilondev_yy_(i,j,k,lnum)
                        Snp1  = factorS_loc  * epsilondev_yy_loc(i,j,k)
                        dom%R_yy_(i_sls,i,j,k,lnum) = alphavalS_loc * dom%R_yy_(i_sls,i,j,k,lnum)  &
                            + betavalS_loc*Sn+gammavalS_loc*Snp1
    
                        ! term in xy
                        Sn    = factorS_loc  * dom%epsilondev_xy_(i,j,k,lnum)
                        Snp1  = factorS_loc  * epsilondev_xy_loc(i,j,k)
                        dom%R_xy_(i_sls,i,j,k,lnum) = alphavalS_loc * dom%R_xy_(i_sls,i,j,k,lnum) &
                            + betavalS_loc*Sn+gammavalS_loc*Snp1
    
                        ! term in xz
                        Sn    = factorS_loc * dom%epsilondev_xz_(i,j,k,lnum)
                        Snp1  = factorS_loc * epsilondev_xz_loc(i,j,k)
                        dom%R_xz_(i_sls,i,j,k,lnum) = alphavalS_loc * dom%R_xz_(i_sls,i,j,k,lnum) &
                            + betavalS_loc*Sn+gammavalS_loc*Snp1
    
                        ! term in yz
                        Sn    = factorS_loc * dom%epsilondev_yz_(i,j,k,lnum)
                        Snp1  = factorS_loc * epsilondev_yz_loc(i,j,k)
                        dom%R_yz_(i_sls,i,j,k,lnum) = alphavalS_loc * dom%R_yz_(i_sls,i,j,k,lnum)  &
                            + betavalS_loc*Sn+gammavalS_loc*Snp1
    
                    enddo   ! end of loop on memory variables
    
                    ! a la fin on stocke la deformation deviatorique pour le prochain
                    ! pas de temps de la mise a jour de Runge-Kutta faite ci-dessus
    
                    ! save deviatoric strain for Runge-Kutta scheme
    
                    !     partie deviatoire
                    dom%epsilondev_xx_(i,j,k,lnum) = epsilondev_xx_loc(i,j,k)
                    dom%epsilondev_yy_(i,j,k,lnum) = epsilondev_yy_loc(i,j,k)
                    dom%epsilondev_xy_(i,j,k,lnum) = epsilondev_xy_loc(i,j,k)
                    dom%epsilondev_xz_(i,j,k,lnum) = epsilondev_xz_loc(i,j,k)
                    dom%epsilondev_yz_(i,j,k,lnum) = epsilondev_yz_loc(i,j,k)
                    !   partie isotrope
                    dom%epsilonvol_   (i,j,k,lnum) = epsilonvol_loc(i,j,k)
    
                end do
            end do
        end do
    
    end subroutine attenuation_update

    subroutine attenuation_aniso_update(dom,lnum,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,&
        ngll,n_solid)

        use sdomain
        implicit none
#include "index.h"

        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: ngll,n_solid
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ

        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: epsilondev_xx_loc,epsilondev_yy_loc
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: factor_loc,alphaval_loc,betaval_loc,gammaval_loc,Sn,Snp1
        integer :: i_sls, i, j, k
        real epsilon_trace_over_3

        if (n_solid>0) then
            do i = 0,ngll-1
                do j = 0,ngll-1
                    do k = 0,ngll-1
                        epsilon_trace_over_3 = 0.333333333333333333333333333333d0 * (DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k))
                        epsilondev_xx_loc(i,j,k) = DXX(i,j,k) - epsilon_trace_over_3
                        epsilondev_yy_loc(i,j,k) = DYY(i,j,k) - epsilon_trace_over_3
                        epsilondev_xy_loc(i,j,k) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                        epsilondev_xz_loc(i,j,k) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                        epsilondev_yz_loc(i,j,k) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                    enddo
                enddo
            enddo
        endif

        do i_sls = 0,n_solid-1

            ! get coefficients for that standard linear solid
            factor_loc  (:,:,:) = dom%Mu_(:,:,:,lnum) * dom%factor_common_3(i_sls,:,:,:,lnum)
            alphaval_loc(:,:,:) = dom%alphaval_3(i_sls,:,:,:,lnum)
            betaval_loc (:,:,:) = dom%betaval_3 (i_sls,:,:,:,lnum)
            gammaval_loc(:,:,:) = dom%gammaval_3(i_sls,:,:,:,lnum)
    
            ! term in xx
            Sn  (:,:,:)  = factor_loc(:,:,:) * dom%epsilondev_xx_(:,:,:,lnum)
            Snp1(:,:,:)  = factor_loc(:,:,:) * epsilondev_xx_loc(:,:,:)
            dom%R_xx_(i_sls,:,:,:,lnum) = alphaval_loc(:,:,:) * dom%R_xx_(i_sls,:,:,:,lnum)  &
                + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)
    
            ! term in yy
            Sn  (:,:,:)  = factor_loc(:,:,:)  * dom%epsilondev_yy_(:,:,:,lnum)
            Snp1(:,:,:)  = factor_loc(:,:,:)  * epsilondev_yy_loc(:,:,:)
            dom%R_yy_(i_sls,:,:,:,lnum) = alphaval_loc(:,:,:) * dom%R_yy_(i_sls,:,:,:,lnum)  &
                + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)
    
            ! term in xy
            Sn  (:,:,:)  = factor_loc(:,:,:)  * dom%epsilondev_xy_(:,:,:,lnum)
            Snp1(:,:,:)  = factor_loc(:,:,:)  * epsilondev_xy_loc(:,:,:)
            dom%R_xy_(i_sls,:,:,:,lnum) = alphaval_loc(:,:,:) * dom%R_xy_(i_sls,:,:,:,lnum) &
                + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)
    
            ! term in xz
            Sn  (:,:,:)  = factor_loc(:,:,:) * dom%epsilondev_xz_(:,:,:,lnum)
            Snp1(:,:,:)  = factor_loc(:,:,:) * epsilondev_xz_loc(:,:,:)
            dom%R_xz_(i_sls,:,:,:,lnum) = alphaval_loc(:,:,:) * dom%R_xz_(i_sls,:,:,:,lnum) &
                + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)
    
            ! term in yz
            Sn  (:,:,:)  = factor_loc(:,:,:) * dom%epsilondev_yz_(:,:,:,lnum)
            Snp1(:,:,:)  = factor_loc(:,:,:) * epsilondev_yz_loc(:,:,:)
            dom%R_yz_(i_sls,:,:,:,lnum) = alphaval_loc(:,:,:) * dom%R_yz_(i_sls,:,:,:,lnum)  &
                + betaval_loc(:,:,:)*Sn(:,:,:)+gammaval_loc(:,:,:)*Snp1(:,:,:)
    
        enddo   ! end of loop on memory variables
    
        ! a la fin on stocke la deformation deviatorique pour le prochain
        ! pas de temps de la mise a jour de Runge-Kutta faite ci-dessus
    
        ! save deviatoric strain for Runge-Kutta scheme
    
        dom%epsilondev_xx_(:,:,:,lnum) = epsilondev_xx_loc(:,:,:)
        dom%epsilondev_yy_(:,:,:,lnum) = epsilondev_yy_loc(:,:,:)
        dom%epsilondev_xy_(:,:,:,lnum) = epsilondev_xy_loc(:,:,:)
        dom%epsilondev_xz_(:,:,:,lnum) = epsilondev_xz_loc(:,:,:)
        dom%epsilondev_yz_(:,:,:,lnum) = epsilondev_yz_loc(:,:,:)
    
    end subroutine attenuation_aniso_update

end module attenuation_solid

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
