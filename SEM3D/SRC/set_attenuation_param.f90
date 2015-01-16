!>
!! \file set_attenuation_param.f90
!! \brief
!!
!<

module attenuation
    use sdomain
    use constants, only : M_PI
    implicit none
contains

    subroutine set_attenuation_param ( Tdomain )
        implicit none
        type(domain),target, intent (INOUT) :: Tdomain

        call set_attenuation_param_solid_nopml (Tdomain)

    end subroutine set_attenuation_param

    subroutine set_attenuation_param_solid_nopml (Tdomain)
        implicit none
        type(domain),target, intent (INOUT) :: Tdomain

        integer :: n, i,j,k, n_solid, i_count, ngllx,nglly,ngllz, mat
        real :: factor_scale_mu0, w_c_source, a_val,b_val, factorS_scale, f0_source, dt, &
            factor_scale_mu, one_minus_sum_betaS, Q_mu, T1_attenuation,T2_attenuation, f0_modele
        !   pour la partie isotrope
        real :: one_minus_sum_betaP, Q_P, factor_scale_P0, aP_val, bP_val,factor_scale_P
        real :: factorP_scale
        ! Dernieres valeurs de Q_mu, Q_P utilisees pour eviter de rappeler compute_iso_Q avec toujours
        ! les meme valeurs
        real :: Q_mu_old, Q_P_old
        real, dimension(:), allocatable :: tau_mu,tau_sigma, betaS, tauinv
        real, dimension(:), allocatable :: P_mu,P_sigma, betaP, Pinv

        ! NL
        interface
           subroutine compute_iso_Q(T1,T2,n,Q,mu,sigma)
               implicit none
               integer :: n
               real :: T1, T2, Q
               real, dimension(n) :: mu, sigma
           end subroutine compute_iso_Q
        end interface
        ! NL
        ! Q_mu ne peut pas etre negatif
        Q_mu_old = -1.d8
        Q_P_old = -1.d8

        f0_modele = 1.d0 / Tdomain%T0_modele

        T1_attenuation = Tdomain%T1_att;   T2_attenuation = Tdomain%T2_att
        f0_source = exp((log(1.d0/T1_attenuation)+log(1.d0/T2_attenuation))/2.d0)
        w_c_source = 2.d0*M_PI*f0_source

        n_solid = Tdomain%n_sls
        allocate (tau_mu(0:n_solid-1))
        allocate (tau_sigma(0:n_solid-1))
        allocate (betaS(0:n_solid-1))
        allocate (P_mu(0:n_solid-1))
        allocate (P_sigma(0:n_solid-1))
        allocate (betaP(0:n_solid-1))

        do n = 0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%PML) cycle
            if (.NOT.(Tdomain%specel(n)%solid)) cycle
            !! if (Tdomain%specel(n)%PML==.false.) then
            !   on met la meme loi partout sinon reflexion a la jonction avec PML

            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz

            do i = 0,ngllx-1
                do j = 0,nglly-1
                    do k = 0,ngllz-1

                        Q_mu = Tdomain%specel(n)%sl%Qs(i,j,k)
                        Q_P = Tdomain%specel(n)%sl%Qp(i,j,k)
                        if (Q_mu .ne. Q_mu_old) then
                            call compute_iso_Q(T1_attenuation,T2_attenuation,n_solid,Q_mu,tau_mu,tau_sigma)
                            Q_mu_old = Q_mu
                        end if
                        if (Q_P .ne. Q_P_old) then
                            call compute_iso_Q(T1_attenuation,T2_attenuation,n_solid,Q_P,P_mu,P_sigma)
                            Q_P_old = Q_P
                        end if

                        betaS(:) = 1.d0 - tau_mu(:) / tau_sigma(:)
                        betaP(:) = 1.d0 - P_mu(:) / P_sigma(:)

                        factor_scale_mu0 = 1.d0
                        factor_scale_P0 = 1.d0
                        !          factor_scale_mu0 = 1.d0 + 2.d0 * log(f0_source/f0_modele) / (PI*Q_mu)
                        !          factor_scale_P0 = 1.d0 + 2.d0 * log(f0_source/f0_modele) / (PI*Q_P)
                        a_val = 0.d0
                        !          a_val = 1.d0
                        b_val = 0.d0
                        aP_val = 0.d0
                        !          aP_val = 1.d0
                        bP_val = 0.d0
                        one_minus_sum_betaS = 1.d0
                        one_minus_sum_betaP = 1.d0

                        do i_count = 0,n_solid-1
                            !   partie deviatoire
                            a_val = a_val + w_c_source * w_c_source * tau_mu(i_count) * &
                                (tau_mu(i_count) - tau_sigma(i_count)) / &
                                (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
                            b_val = b_val + w_c_source * (tau_mu(i_count) - tau_sigma(i_count)) / &
                                (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
                            one_minus_sum_betaS = one_minus_sum_betaS - betaS(i_count)
                            !   partie isotrope
                            aP_val = aP_val + w_c_source * w_c_source * P_mu(i_count) * &
                                (P_mu(i_count) - P_sigma(i_count)) / &
                                (1.d0 + w_c_source*w_c_source*P_mu(i_count)*P_mu(i_count))
                            bP_val = b_val + w_c_source * (P_mu(i_count) - P_sigma(i_count)) / &
                                (1.d0 + w_c_source*w_c_source*P_mu(i_count)*P_mu(i_count))
                            one_minus_sum_betaP = one_minus_sum_betaP - betaP(i_count)
                        enddo

                        Tdomain%specel(n)%sl%onemSbeta(i,j,k) = one_minus_sum_betaS
                        Tdomain%specel(n)%sl%onemPbeta(i,j,k) = one_minus_sum_betaP

                        a_val = 1.d0 - a_val
                        aP_val = 1.d0 - aP_val

                        !--- quantity by which to scale mu to get mu_relaxed
                        factor_scale_mu = a_val*( (sqrt(1.d0 + b_val*b_val/(a_val*a_val))+1.d0))/2.
                        factor_scale_P = aP_val*( (sqrt(1.d0 + bP_val*bP_val/(aP_val*aP_val))+1.d0))/2.

                        !--- total factor by which to scale mu0
                        factorS_scale = factor_scale_mu * factor_scale_mu0
                        factorP_scale = factor_scale_P * factor_scale_P0

                        !          print*,' one_minus_sum_betaS ',one_minus_sum_betaS
                        !          print*,' one_minus_sum_betaP ',one_minus_sum_betaP
                        !          print*,' b_val  ',b_val
                        !          print*,' a_val  ',a_val
                        !          print*,' w_c_source  ',w_c_source
                        !          print*,' factorS_scale ',factorS_scale
                        !          stop
                        !--- check that the correction factor is close to one
                        if (factorS_scale < 0.1d0 .or. factorS_scale > 1.1d0) then
                            print*,' one_minus_sum_betaS ',one_minus_sum_betaS
                            print*,' one_minus_sum_betaP ',one_minus_sum_betaP
                            print*,' b_val  ',b_val
                            print*,' a_val  ',a_val
                            print*,' w_c_source  ',w_c_source
                            print*,' factorS_scale ',factorS_scale
                            stop 'incorrect correction factor in attenuation model deviatoric part'
                        endif
                        if (factorP_scale < 0.1d0 .or. factorP_scale > 1.1d0) then
                            stop 'incorrect correction factor in attenuation model isotropic part'
                        endif

                        ! on corrige mu, mais avant on remplace lambda par kappa
                        ! (ATTENTION a ne pas se planter dans le suite) :
                        !          Tdomain%specel(n)%Lambda(i,j,k) = Tdomain%specel(n)%Lambda(i,j,k) + &
                        !                                            2.d0/3.d0 * Tdomain%specel(n)%Mu(i,j,k)
                        Tdomain%specel(n)%Mu(i,j,k) = factorS_scale * Tdomain%specel(n)%Mu(i,j,k)
                        Tdomain%specel(n)%Kappa(i,j,k) = factorP_scale * Tdomain%specel(n)%Kappa(i,j,k)

                        !on initialise les coefs Runge-Kutta:
                        mat = Tdomain%specel(n)%mat_index
                        dt = Tdomain%sSubdomain(mat)%dt
                        allocate (tauinv(0:n_solid-1))
                        allocate (Pinv(0:n_solid-1))
                        tauinv(:) = - 1.d0 / tau_sigma(:)
                        Pinv(:) = - 1.d0 / P_sigma(:)
                        !!         partie deviatorique
                        Tdomain%specel(n)%sl%factor_common_3(:,i,j,k) = 2.d0 * betaS(:) * tauinv(:)
                        Tdomain%specel(n)%sl%alphaval_3(:,i,j,k) = 1.d0 + dt*tauinv(:) + dt**2*tauinv(:)**2 / 2.d0 + &
                            dt**3*tauinv(:)**3 / 6.d0 + dt**4*tauinv(:)**4 / 24.d0
                        Tdomain%specel(n)%sl%betaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 3.d0 + &
                            dt**3*tauinv(:)**2 / 8.d0 + dt**4*tauinv(:)**3 / 24.d0
                        Tdomain%specel(n)%sl%gammaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 6.d0 + &
                            dt**3*tauinv(:)**2 / 24.d0
                        deallocate (tauinv)
                        !!         partie istotrope
                        Tdomain%specel(n)%sl%factor_common_P(:,i,j,k) = 2.d0 * betaP(:) * Pinv(:)
                        Tdomain%specel(n)%sl%alphaval_P(:,i,j,k) = 1.d0 + dt*Pinv(:) + dt**2*Pinv(:)**2 / 2.d0 + &
                            dt**3*Pinv(:)**3 / 6.d0 + dt**4*Pinv(:)**4 / 24.d0
                        Tdomain%specel(n)%sl%betaval_P(:,i,j,k) = dt / 2.d0 + dt**2*Pinv(:) / 3.d0 + &
                            dt**3*Pinv(:)**2 / 8.d0 + dt**4*Pinv(:)**3 / 24.d0
                        Tdomain%specel(n)%sl%gammaval_P(:,i,j,k) = dt / 2.d0 + dt**2*Pinv(:) / 6.d0 + &
                            dt**3*Pinv(:)**2 / 24.d0
                        deallocate (Pinv)
                        !         partie deviatorique
                        !          Tdomain%specel(n)%factor_common_3(:,i,j,k) = 0.
                        !          Tdomain%specel(n)%alphaval_3(:,i,j,k) = 0.
                        !          Tdomain%specel(n)%betaval_3(:,i,j,k) = 0.
                        !          Tdomain%specel(n)%gammaval_3(:,i,j,k) = 0.
                        !          deallocate (tauinv)
                        !         partie istotrope
                        !          Tdomain%specel(n)%factor_common_P(:,i,j,k) = 0.
                        !          Tdomain%specel(n)%alphaval_P(:,i,j,k) = 0.
                        !          Tdomain%specel(n)%betaval_P(:,i,j,k) =0.
                        !          Tdomain%specel(n)%gammaval_P(:,i,j,k) = 0.
                        !          deallocate (Pinv)
                    enddo
                enddo
            enddo
        enddo

        deallocate (tau_mu, tau_sigma, betaS)
        deallocate (P_mu, P_sigma, betaP)
    end subroutine set_attenuation_param_solid_nopml

    subroutine set_attenuation_aniso_param (Tdomain)

        ! Written by Paul Cupillard 29/06/2006


        use sdomain
        use constants, only : M_PI

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain

        integer :: n, i,j,k, n_solid, i_count, ngllx,nglly,ngllz, mat
        real :: factor_scale_mu0, w_c_source, a_val,b_val, big_omega, factor_scale, f0_source, dt, &
            factor_scale_mu, one_minus_sum_beta, Q_mu, T1_attenuation,T2_attenuation, f0_modele
        real, dimension(:), allocatable :: tau_mu,tau_sigma, beta, tauinv

        f0_modele = 1.d0 / Tdomain%T0_modele

        T1_attenuation = Tdomain%T1_att;   T2_attenuation = Tdomain%T2_att
        f0_source = exp((log(1.d0/T1_attenuation)+log(1.d0/T2_attenuation))/2.d0)
        w_c_source = 2.d0*M_PI*f0_source

        n_solid = Tdomain%n_sls
        allocate (tau_mu(0:n_solid-1))
        allocate (tau_sigma(0:n_solid-1))
        allocate (beta(0:n_solid-1))

        do n = 0,Tdomain%n_elem-1
            !! if (Tdomain%specel(n)%PML==.false.) then
            if (.NOT.(Tdomain%specel(n)%PML)) then

                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz

                do i = 0,ngllx-1
                    do j = 0,nglly-1
                        do k = 0,ngllz-1

                            Q_mu = Tdomain%specel(n)%sl%Q(i,j,k)
                            call compute_constant_Q(T1_attenuation,T2_attenuation,n_solid,Q_mu,tau_mu,tau_sigma)

                            beta(:) = 1.d0 - tau_mu(:) / tau_sigma(:)

                            factor_scale_mu0 = 1.d0 + 2.d0 * log(f0_source/f0_modele) / (M_PI*Q_mu)
                            a_val = 1.d0
                            b_val = 0.d0
                            one_minus_sum_beta = 1.d0

                            do i_count = 0,n_solid-1
                                a_val = a_val - w_c_source * w_c_source * tau_mu(i_count) * &
                                    (tau_mu(i_count) - tau_sigma(i_count)) / &
                                    (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
                                b_val = b_val + w_c_source * (tau_mu(i_count) - tau_sigma(i_count)) / &
                                    (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
                                one_minus_sum_beta = one_minus_sum_beta - beta(i_count)
                            enddo

                            Tdomain%specel(n)%sl%onemSbeta(i,j,k) = one_minus_sum_beta

                            big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0);

                            !--- quantity by which to scale mu to get mu_relaxed
                            factor_scale_mu = b_val * b_val / (2.d0 * big_omega)

                            !--- total factor by which to scale mu0
                            factor_scale = factor_scale_mu * factor_scale_mu0

                            !--- check that the correction factor is close to one
                            if (factor_scale < 0.9d0 .or. factor_scale > 1.1d0) &
                                stop 'incorrect correction factor in attenuation model'

                            ! on corrige mu, mais avant on remplace lambda par kappa
                            ! (ATTENTION a ne pas se planter dans le suite) :
                            Tdomain%specel(n)%Lambda(i,j,k) = Tdomain%specel(n)%Lambda(i,j,k) + &
                                2.d0/3.d0 * Tdomain%specel(n)%Mu(i,j,k)
                            Tdomain%specel(n)%Mu(i,j,k) = factor_scale * Tdomain%specel(n)%Mu(i,j,k)

                            !on initialise les coefs Runge-Kutta:
                            mat = Tdomain%specel(n)%mat_index
                            dt = Tdomain%sSubdomain(mat)%dt
                            allocate (tauinv(0:n_solid-1))
                            tauinv(:) = - 1.d0 / tau_sigma(:)
                            Tdomain%specel(n)%sl%factor_common_3(:,i,j,k) = 2.d0 * beta(:) * tauinv(:)
                            Tdomain%specel(n)%sl%alphaval_3(:,i,j,k) = 1.d0 + dt*tauinv(:) + dt**2*tauinv(:)**2 / 2.d0 + &
                                dt**3*tauinv(:)**3 / 6.d0 + dt**4*tauinv(:)**4 / 24.d0
                            Tdomain%specel(n)%sl%betaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 3.d0 + &
                                dt**3*tauinv(:)**2 / 8.d0 + dt**4*tauinv(:)**3 / 24.d0
                            Tdomain%specel(n)%sl%gammaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 6.d0 + &
                                dt**3*tauinv(:)**2 / 24.d0
                            deallocate (tauinv)
                        enddo
                    enddo
                enddo

            endif
        enddo

        deallocate (tau_mu, tau_sigma, beta)

    end subroutine set_attenuation_aniso_param

end module attenuation
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
