!>
!! \file set_attenuation_param.f90
!! \brief
!!
!<

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!--  linear anelasticity : association of standard linear solids (Zener) in
!    parallel. See Dahlen & Tromp (98) or Moczo & Krystek (2005) for
!    clear and correct explanations about tinker-toy models.
!--  these authors gave correct formulae, as Lombard & Piraux (2011) did too. What
!    follows is strongly influenced by that last paper.
!-------------------------------------------------------------------
!-------------------------------------------------------------------
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


    !---------------------------------------------------------
    !---------------------------------------------------------
    subroutine set_attenuation_param_solid_nopml (Tdomain)
        type(domain),target, intent (INOUT) :: Tdomain
        integer    :: n_solid,n_freq_q
        real       :: f_c_source,f_ref,dt
        real, dimension(:), allocatable   :: omega_tau_s,omega_Q
        real, dimension(:), allocatable   :: agamma_mu,agamma_kappa
        integer :: n,i,j,k,ngllx,nglly,ngllz,mat
        real :: Q_mu, Q_kappa, Q_mu_old, Q_kappa_old

        !- last values of Q factors (allows to save time if same Qs)
        Q_mu_old = -1.d8
        Q_kappa_old = -1.d8

        !- initialization of general parameters, valid for each GLL point, and mechanisms
        !   related to shear or compression
        call init_attenu(Tdomain,n_solid,n_freq_q,f_ref,f_c_source,omega_tau_s,omega_Q)

        !-----------------------------------------------------------------------------------
        !- here we will look (at each GLL where a given value of Q is given) for the values
        !    of gamma, a combination of the stress and strain relaxation
        !    times associated to each standard linear (Zener) solid of our system:
        !         gamma = 1/N(tau_eps/tau_s-1),
        !    where N is the number of SLS.
        !- this gamma is the same as the Lombard and Piraux's kappa, or equivalent to the gamma of
        !    Peyrusse et al., 2014, modulo a factor (relaxed modulus/unrelaxed modulus)
        !- in fact these gammas are the modulus defect (normalized by the relaxed modulus) for each SLS;
        !-  these gamma parameters are determinated following the assumption of a constant Q on
        !   the seismic band of interest. To be generalized.

        !- once we get these gamma values, everything's necessary for taking account of anelasticity,
        !    can be calculated: the relaxed modulus, the modulus defect.
        !-----------------------------------------------------------------------------------

        allocate(agamma_mu(0:n_solid-1))
        allocate(agamma_kappa(0:n_solid-1))

        do n = 0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%PML) cycle
            if (.NOT.(Tdomain%specel(n)%solid)) cycle

            ! a faire dans PMLs...
                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz

                do i = 0,ngllx-1
                    do j = 0,nglly-1
                        do k = 0,ngllz-1

                            Q_mu = Tdomain%specel(n)%sl%Qs(i,j,k)
                            Q_kappa = Tdomain%specel(n)%sl%Qp(i,j,k)

                            !- from Qs to gammas
                            if (Q_mu .ne. Q_mu_old) then
                                call inv_gamma_Q_const(n_solid,n_freq_q,Q_mu,omega_Q,omega_tau_s,agamma_mu)
                                Q_mu_old = Q_mu
                            end if
                            if (Q_kappa .ne. Q_kappa_old) then
                                call inv_gamma_Q_const(n_solid,n_freq_q,Q_kappa,omega_Q,omega_tau_s,agamma_kappa)
                                Q_kappa_old = Q_kappa
                            end if


                            !- getting the values of the relaxed moduli
                            Tdomain%specel(n)%Mu(i,j,k) = Tdomain%specel(n)%Mu(i,j,k)*   &
                                get_relaxed_modulus(n_solid,Q_mu,f_ref,   &
                                f_c_source,omega_tau_s,agamma_mu)

                            Tdomain%specel(n)%Kappa(i,j,k) = Tdomain%specel(n)%Kappa(i,j,k)*  &
                                get_relaxed_modulus(n_solid,Q_kappa,f_ref,     &
                                f_c_source,omega_tau_s,agamma_kappa)


                            !- factor to get from relaxed to unrelaxed modulus: M_U = M_R*(1+delta_M/M_R)
                            Tdomain%specel(n)%sl%onemSbeta(i,j,k) = 1d0+get_modulus_defect(n_solid, agamma_mu)
                            Tdomain%specel(n)%sl%onemPbeta(i,j,k) = 1d0+get_modulus_defect(n_solid, agamma_kappa)



                            !- Runge-kutta parameters for the time integration of the terms related to the relaxation function
                            mat = Tdomain%specel(n)%mat_index
                            dt = Tdomain%sSubdomain(mat)%dt
                            call RK4_attenu_coefficients(n_solid,dt,omega_tau_s,agamma_mu,   &
                                Tdomain%specel(n)%sl%factor_common_3(:,i,j,k),             &
                                Tdomain%specel(n)%sl%alphaval_3(:,i,j,k),                  &
                                Tdomain%specel(n)%sl%betaval_3(:,i,j,k),                   &
                                Tdomain%specel(n)%sl%gammaval_3(:,i,j,k))
                            call RK4_attenu_coefficients(n_solid,dt,omega_tau_s,agamma_kappa,   &
                                Tdomain%specel(n)%sl%factor_common_P(:,i,j,k),             &
                                Tdomain%specel(n)%sl%alphaval_P(:,i,j,k),                  &
                                Tdomain%specel(n)%sl%betaval_P(:,i,j,k),                   &
                                Tdomain%specel(n)%sl%gammaval_P(:,i,j,k))


                        enddo
                    enddo
                enddo

        enddo


        deallocate(agamma_mu,agamma_kappa,omega_tau_s,omega_Q)



    end subroutine set_attenuation_param_solid_nopml

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    subroutine init_attenu(Tdomain,n_solid,n_freq_q,f_ref,f_source,omega_tau_s,omega_Q)
        implicit none
        type(domain), intent(in)  :: Tdomain
        integer, intent(out)  :: n_solid, n_freq_q
        real, intent(out)     :: f_ref,f_source
        real, dimension(:), allocatable, intent(out) :: omega_tau_s,omega_Q
        real  :: f_min,f_max, f_tmp

        !- number of standard linear solids
        n_solid = Tdomain%n_sls
        !- seismic periods of interest - bounds of the range => frequency range
        f_max = 1d0/Tdomain%T2_att
        f_min = 1d0/Tdomain%T1_att
        if (f_max<f_min) then
            f_tmp = f_max
            f_max = f_min
            f_min = f_tmp
        endif

        f_source = exp((log(f_min)+log(f_max))/2.d0)

        f_ref = 1d0/Tdomain%T0_modele

        !- f_ref must be in the attenuation band
        if((f_ref < f_min) .or. (f_ref > f_max))  &
            stop "  --> In set_attenuation_param.f90 : Bad value of T0_modele, which is not in the attenuation band"

        !- frequency sampling for 1/stress relaxation times (1/tau_s)
        allocate(omega_tau_s(0:n_solid-1))
        call freq_sampl(n_solid,f_min,f_max,omega_tau_s)

        !- frequency sampling for Q
        !- number of frequencies for which Q will be given
        !- here it is the value chosen by Lombard & Piraux (in fact taking
        !   more values does not change the result).
        n_freq_q = 2*n_solid-1
        allocate(omega_Q(0:n_freq_q-1))
        call freq_sampl(n_freq_q,f_min,f_max,omega_Q)

    end subroutine init_attenu
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    subroutine freq_sampl(nval,fmin,fmax,omega)
        ! returns the values of frequencies on range [fmin,fmax], logarithmic sampling
        implicit none
        integer, intent(in) :: nval
        real, intent(in) :: fmin,fmax
        real, dimension(0:nval-1), intent(out) :: omega
        integer   :: i

        if(nval == 1)then
            omega(0) = dsqrt(fmin*fmax)
        else
            do i = 0,nval-1
                omega(i) = 2d0*M_PI*fmin*(fmax/fmin)**(i/(nval-1d0))
            end do
        end if

    end subroutine freq_sampl
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    subroutine inv_gamma_Q_const(n_solid,n_freq_q,Qref,omega_Q,omega_tau_s,agamma)
        !-- from Qref (data) constant and tau_s (a priori fixed), gives tau_e for
        !   each standard linear solid
        !-- for the time being: only a classical least-squares fit, that means linear
        !   operators in the inversion process. Reliability to check thoroughly.
        use mleastsq, only : cg_inv
        implicit none
        integer, intent(in)  :: n_solid,n_freq_q
        real, intent(in)     :: Qref
        real, dimension(0:n_freq_q-1), intent(in) :: omega_q
        real, dimension(0:n_solid-1), intent(in) :: omega_tau_s
        real, dimension(0:n_solid-1), intent(out) :: agamma
        real, dimension(0:n_freq_q-1) :: data_vec
        real, dimension(0:n_freq_q-1,0:n_solid-1) :: Amat
        integer  :: nf,ns

        !- coefficients of the forward problem matrix
        do nf = 0,n_freq_q-1
            call G_gamma(n_solid,Qref,omega_q(nf),omega_tau_s(0:),Amat(nf,0:))
        end do

        !- singular value decomposition, and inverse problem solution for gamma coefficients
        ! data vector: Q-1 in each slot
        data_vec(0:) = 1d0/Qref
        call cg_inv(n_freq_q,n_solid,Amat(0:,0:),data_vec(0:),agamma(0:))

        ! call cg_inv(n_solid,n_solid,MATMUL(TRANSPOSE(Amat(0:,0:)),Amat(0:,0:)),MATMUL(TRANSPOSE(Amat(0:,0:)),data_vec(0:)),agamma(0:))

        !- values of gamma physically correct? (they must be positive)
        do ns = 0,n_solid-1
            if(agamma(ns) <= 0d0)then
                print*, "  --> One of the gamma values is negative: this is not physically plausible."
                stop "    --> Try another run with less linear solids."
            end if
        end do

        !- now we could get values of tau_sigma and tau_epsilon: the 1st is just the inverse
        !   of omega_tau_s, the second comes from: kappa = 1/N(tau_eps/tau_sigma-1)
        !  tau_s(0:) = 1d0/omega_tau_s(0:)
        !  tau_e(0:) = tau_s(0:)*(1d0+n_solid*kappa(0:))
        !- but in fact it is unnecessary: everything's needed in the SEM's attenuation process
        !  is caught with gamma and omega_tau_s

        !- graphical checking
        !  call print_Q_model(n_solid,omega_tau_s,agamma,omega_tau_s(0)/2d0/M_PI,omega_tau_s(n_solid-1)/2d0/M_PI)

    end subroutine inv_gamma_Q_const
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    subroutine G_gamma(nsol,Q,omega_q,omega_tau,Amat)
        ! coefficients of the forward problem (omega_tau -> Q)
        implicit none
        integer, intent(in)  :: nsol
        real, intent(in)     :: Q,omega_q
        real, dimension(0:nsol-1),intent(in) :: omega_tau
        real, dimension(0:nsol-1),intent(out) :: Amat
        real    :: xx,yy
        integer :: i

        do i = 0,nsol-1
            xx = omega_q*(omega_tau(i)-omega_q/Q)
            yy = omega_q*omega_q+omega_tau(i)*omega_tau(i)
            Amat(i) = xx/yy
        end do

    end subroutine G_gamma
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    subroutine Q_direct(nsol,omega_tau,agamma,omega_sampl,qinv)
        ! getting the Q-value, given the sampling and the Gamma values (= values
        !   of tau_s and tau_e for each SLS)
        implicit none
        integer, intent(in)  :: nsol
        real, intent(in)     :: omega_sampl
        real, dimension(0:nsol-1), intent(in) :: omega_tau,agamma
        real, intent(out) :: qinv
        real    :: xx,yy
        integer :: i

        xx = 0d0 ; yy = 1d0

        do i = 0,nsol-1
            xx = xx+(omega_sampl*omega_tau(i)*agamma(i))/(omega_tau(i)**2+omega_sampl**2)
            yy = yy+(omega_sampl**2*agamma(i))/(omega_tau(i)**2+omega_sampl**2)
        end do
        qinv = xx/yy


    end subroutine Q_direct
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    real function get_relaxed_modulus(n_sol,Qref,f_ref,f_source,omega_s,agamma)
        implicit none
        integer, intent(in)  :: n_sol
        real, intent(in)  :: f_ref,f_source,Qref
        real, dimension(0:n_sol-1), intent(in) :: omega_s,agamma
        integer  :: ns
        real :: theta_1,theta_2,R,ws,r_m_part,omega_r

        ws = 2d0*M_PI*f_source

        !-  factor to get from the modulus M at f_ref (at the outset in SEM, from the material file)
        !    to the relaxed modulus (zero frequency)

        !- from M(f_ref) to M(f_source) (at the center of the frequency band of interest)
        r_m_part = 1d0+2d0*log(f_source/f_ref)/(M_PI*Qref)

        !- from M(f_source) to M_relaxed: see for instance expressions of M(omega) in Moczo and Kryztek,
        !   and procedure to get M_relaxed in Liu, Anderson and Kanamori (1976) or M_unrelaxed in Peyrusse et al.
        !- note that in Liu et al. this is not a GZB that is dealt

        !- (1/v(omega)) = real part (density/M(omega))^0.5)
        theta_1 = 1d0 ; theta_2 = 0d0 ; R = 0d0

        do ns = 0,n_sol-1
            omega_r = ws/omega_s(ns)
            theta_1 = theta_1 + (omega_r*omega_r*agamma(ns))/(1d0+omega_r*omega_r)
            theta_2 = theta_2 + (omega_r*agamma(ns))/(1d0+omega_r*omega_r)
        end do

        R = dsqrt(theta_1*theta_1 + theta_2*theta_2)

        get_relaxed_modulus = r_m_part*(theta_1+R)/(2d0*R*R)
        if(get_relaxed_modulus > 1d0)   &
            stop "  --> Unphysical value for the relaxed modulus."

    end function get_relaxed_modulus
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    real function get_modulus_defect(n_solid, agamma)
        implicit none
        !- this function returns the value of the modulus defect, relative to the
        !    relaxed modulus:
        !    delta_M/M_relaxed = (M_unrelaxed - M_relaxed)/M_relaxed = sum_i gamma_i
        integer, intent(in) :: n_solid
        real, dimension(0:n_solid-1), intent(in) :: agamma

        get_modulus_defect = SUM(agamma)

    end function get_modulus_defect
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    subroutine RK4_attenu_coefficients(n_solid,dt,omega_tau_s,agamma,   &
        factor_common,alphaval,betaval,gammaval)
        !- routine returns the coefficients for the time integration of the
        !  relaxation function M(t)
        integer, intent(in)  :: n_solid
        real, intent(in)  :: dt
        real, dimension(0:n_solid-1), intent(in) :: omega_tau_s,agamma
        real, dimension(0:n_solid-1), intent(out) :: factor_common,alphaval,betaval,gammaval
        real, dimension(0:n_solid-1) :: invtau

        invtau(:) = -1d0*omega_tau_s(:)
        factor_common = 2d0 * agamma(:) * omega_tau_s(:)
        alphaval = 1d0 + dt*invtau(:) + dt**2*invtau(:)**2 / 2.d0 +   &
            dt**3*invtau(:)**3 / 6.d0 + dt**4*invtau(:)**4 / 24.d0
        betaval = dt / 2.d0 + dt**2*invtau(:) / 3.d0 + &
            dt**3*invtau(:)**2 / 8.d0 + dt**4*invtau(:)**3 / 24.d0
        gammaval = dt / 2.d0 + dt**2*invtau(:) / 6.d0 + &
            dt**3*invtau(:)**2 / 24.d0


    end subroutine RK4_attenu_coefficients
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    subroutine print_Q_model(n_solid,omega_tau_s,agamma,fmin,fmax)

        implicit none
        integer, intent(in)  :: n_solid
        real, dimension(0:n_solid-1), intent(in) :: omega_tau_s
        real, dimension(0:n_solid-1), intent(in) :: agamma
        real, intent(in)  :: fmin,fmax
        real  :: qinv
        integer  :: i
        integer, parameter  :: n_sampl = 5000
        real, dimension(:), allocatable :: omega_sampl

        allocate(omega_sampl(0:n_sampl-1))
        open(10,file="Q_verif",action="write",status="replace")
        call freq_sampl(n_sampl,fmin/10,fmax*10,omega_sampl)
        do i = 0,size(omega_sampl)-1
            call Q_direct(n_solid,omega_tau_s,agamma,omega_sampl(i),qinv)
            write(10,*) omega_sampl(i)/2d0/M_PI,qinv
        end do
        close(10)
        deallocate(omega_sampl)

    end subroutine print_Q_model



    subroutine set_attenuation_aniso_param (Tdomain)
        !! TODO...
        use sdomain
        use constants, only : M_PI

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain

        call set_attenuation_param(Tdomain)
!        integer :: n, i,j,k, n_solid, i_count, ngllx,nglly,ngllz, mat
!        real :: factor_scale_mu0, w_c_source, a_val,b_val, big_omega, factor_scale, f0_source, dt, &
!            factor_scale_mu, one_minus_sum_beta, Q_mu, T1_attenuation,T2_attenuation, f0_modele
!        real, dimension(:), allocatable :: tau_mu,tau_sigma, beta, tauinv
!
!        f0_modele = 1.d0 / Tdomain%T0_modele
!
!        T1_attenuation = Tdomain%T1_att;   T2_attenuation = Tdomain%T2_att
!        f0_source = exp((log(1.d0/T1_attenuation)+log(1.d0/T2_attenuation))/2.d0)
!        w_c_source = 2.d0*M_PI*f0_source
!
!        n_solid = Tdomain%n_sls
!        allocate (tau_mu(0:n_solid-1))
!        allocate (tau_sigma(0:n_solid-1))
!        allocate (beta(0:n_solid-1))
!
!        do n = 0,Tdomain%n_elem-1
!            !! if (Tdomain%specel(n)%PML==.false.) then
!            if (.NOT.(Tdomain%specel(n)%PML)) then
!
!                ngllx = Tdomain%specel(n)%ngllx
!                nglly = Tdomain%specel(n)%nglly
!                ngllz = Tdomain%specel(n)%ngllz
!
!                do i = 0,ngllx-1
!                    do j = 0,nglly-1
!                        do k = 0,ngllz-1
!
!                            Q_mu = Tdomain%specel(n)%sl%Q(i,j,k)
!                            call compute_constant_Q(T1_attenuation,T2_attenuation,n_solid,Q_mu,tau_mu,tau_sigma)
!
!                            beta(:) = 1.d0 - tau_mu(:) / tau_sigma(:)
!
!                            factor_scale_mu0 = 1.d0 + 2.d0 * log(f0_source/f0_modele) / (M_PI*Q_mu)
!                            a_val = 1.d0
!                            b_val = 0.d0
!                            one_minus_sum_beta = 1.d0
!
!                            do i_count = 0,n_solid-1
!                                a_val = a_val - w_c_source * w_c_source * tau_mu(i_count) * &
!                                    (tau_mu(i_count) - tau_sigma(i_count)) / &
!                                    (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
!                                b_val = b_val + w_c_source * (tau_mu(i_count) - tau_sigma(i_count)) / &
!                                    (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
!                                one_minus_sum_beta = one_minus_sum_beta - beta(i_count)
!                            enddo
!
!                            Tdomain%specel(n)%sl%onemSbeta(i,j,k) = one_minus_sum_beta
!
!                            big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0);
!
!                            !--- quantity by which to scale mu to get mu_relaxed
!                            factor_scale_mu = b_val * b_val / (2.d0 * big_omega)
!
!                            !--- total factor by which to scale mu0
!                            factor_scale = factor_scale_mu * factor_scale_mu0
!
!                            !--- check that the correction factor is close to one
!                            if (factor_scale < 0.9d0 .or. factor_scale > 1.1d0) &
!                                stop 'incorrect correction factor in attenuation model'
!
!                            ! on corrige mu, mais avant on remplace lambda par kappa
!                            ! (ATTENTION a ne pas se planter dans le suite) :
!                            Tdomain%specel(n)%Lambda(i,j,k) = Tdomain%specel(n)%Lambda(i,j,k) + &
!                                2.d0/3.d0 * Tdomain%specel(n)%Mu(i,j,k)
!                            Tdomain%specel(n)%Mu(i,j,k) = factor_scale * Tdomain%specel(n)%Mu(i,j,k)
!
!                            !on initialise les coefs Runge-Kutta:
!                            mat = Tdomain%specel(n)%mat_index
!                            dt = Tdomain%sSubdomain(mat)%dt
!                            allocate (tauinv(0:n_solid-1))
!                            tauinv(:) = - 1.d0 / tau_sigma(:)
!                            Tdomain%specel(n)%sl%factor_common_3(:,i,j,k) = 2.d0 * beta(:) * tauinv(:)
!                            Tdomain%specel(n)%sl%alphaval_3(:,i,j,k) = 1.d0 + dt*tauinv(:) + dt**2*tauinv(:)**2 / 2.d0 + &
!                                dt**3*tauinv(:)**3 / 6.d0 + dt**4*tauinv(:)**4 / 24.d0
!                            Tdomain%specel(n)%sl%betaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 3.d0 + &
!                                dt**3*tauinv(:)**2 / 8.d0 + dt**4*tauinv(:)**3 / 24.d0
!                            Tdomain%specel(n)%sl%gammaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 6.d0 + &
!                                dt**3*tauinv(:)**2 / 24.d0
!                            deallocate (tauinv)
!                        enddo
!                    enddo
!                enddo
!
!            endif
!        enddo
!
!        deallocate (tau_mu, tau_sigma, beta)
!
    end subroutine set_attenuation_aniso_param

end module attenuation
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
