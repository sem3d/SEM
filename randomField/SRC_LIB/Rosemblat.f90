module rosemblat

    use mpi
    implicit none

    private

    ! Public transformation subroutines
    public :: transform_Gaussian
    public :: transform_Lognormal
    public :: transform_Uniform
    public :: transform_Weibull
    public :: transform_Gamma

    ! Public auxiliary statistical functions
    public :: std_normal_cdf
    public :: std_normal_ccdf
    public :: digamma
    public :: solve_weibull_params
    public :: solve_gamma_params

    ! Numerical constants
    double precision, parameter :: PI = 3.141592653589793238462643383279502884d0
    double precision, parameter :: SQRT_2 = 1.414213562373095048801688724209698079d0
    double precision, parameter :: SQRT_3 = 1.732050807568877293527446341505872367d0
    double precision, parameter :: EPS_TAIL = 1.0d-300
    double precision, parameter :: ONE_MINUS_EPS = 0.9999999999999999d0

    ! External dcdflib routine for Gamma distribution
    interface
        subroutine cdfgam(which, p, q, x, shape, scale, status, bound)
            integer, intent(in) :: which
            double precision, intent(inout) :: p, q, x, shape, scale, bound
            integer, intent(out) :: status
        end subroutine cdfgam
    end interface

contains

    !-----------------------------------------------------------------------------------------------
    ! Cumulative Distribution Function of Standard Normal Distribution: Phi(z)
    !-----------------------------------------------------------------------------------------------
    pure function std_normal_cdf(z) result(p)
        double precision, intent(in) :: z
        double precision :: p
        p = 0.5d0 * erfc(-z / SQRT_2)
    end function std_normal_cdf

    !-----------------------------------------------------------------------------------------------
    ! Complementary CDF of Standard Normal Distribution: 1 - Phi(z) = Phi(-z)
    ! Numerically stable for positive extreme tails without catastrophic cancellation.
    !-----------------------------------------------------------------------------------------------
    pure function std_normal_ccdf(z) result(q)
        double precision, intent(in) :: z
        double precision :: q
        q = 0.5d0 * erfc(z / SQRT_2)
    end function std_normal_ccdf

    !-----------------------------------------------------------------------------------------------
    ! Digamma function: psi(x) = d/dx ln(Gamma(x))
    ! Computes digamma to machine precision using recurrence and asymptotic expansion.
    !-----------------------------------------------------------------------------------------------
    pure function digamma(x) result(psi)
        double precision, intent(in) :: x
        double precision :: psi, y, r

        y = x
        psi = 0.0d0
        ! Recurrence relation: psi(y+1) = psi(y) + 1/y
        do while (y < 6.0d0)
            psi = psi - 1.0d0 / y
            y = y + 1.0d0
        end do

        ! Asymptotic expansion for y >= 6
        r = 1.0d0 / (y * y)
        psi = psi + log(y) - 0.5d0 / y - r * (1.0d0/12.0d0 - r * (1.0d0/120.0d0 &
              - r * (1.0d0/252.0d0 - r * (1.0d0/240.0d0))))
    end function digamma

    !-----------------------------------------------------------------------------------------------
    ! Abort calculation immediately via MPI upon error
    !-----------------------------------------------------------------------------------------------
    subroutine rosemblat_abort(message)
        character(len=*), intent(in) :: message
        integer :: mpi_err

        write(*,*) ""
        write(*,*) "================================================================================"
        write(*,*) "FATAL ERROR in module rosemblat:"
        write(*,*) trim(message)
        write(*,*) "Aborting MPI calculation immediately..."
        write(*,*) "================================================================================"
        write(*,*) ""
        call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_err)
    end subroutine rosemblat_abort

    !-----------------------------------------------------------------------------------------------
    ! Transform to Gaussian Distribution
    ! X = fieldAvg + sqrt(fieldVar) * Z
    !-----------------------------------------------------------------------------------------------
    subroutine transform_Gaussian(fieldAvg, fieldVar, randField)
        double precision, intent(in) :: fieldAvg, fieldVar
        double precision, dimension(:,:), intent(inout) :: randField

        if (fieldVar < 0.0d0) then
            call rosemblat_abort("transform_Gaussian: fieldVar must be non-negative.")
        end if

        randField(:,:) = randField(:,:) * sqrt(fieldVar) + fieldAvg
    end subroutine transform_Gaussian

    !-----------------------------------------------------------------------------------------------
    ! Transform to Lognormal Distribution
    ! X = exp(mu_Y + sigma_Y * Z)
    !-----------------------------------------------------------------------------------------------
    subroutine transform_Lognormal(fieldAvg, fieldVar, randField)
        double precision, intent(in) :: fieldAvg, fieldVar
        double precision, dimension(:,:), intent(inout) :: randField
        double precision :: normalVar, normalAvg

        if (fieldAvg <= 0.0d0) then
            call rosemblat_abort("transform_Lognormal: fieldAvg must be strictly positive (> 0).")
        end if
        if (fieldVar <= 0.0d0) then
            call rosemblat_abort("transform_Lognormal: fieldVar must be strictly positive (> 0).")
        end if

        normalVar = log(1.0d0 + fieldVar / (fieldAvg**2))
        normalAvg = log(fieldAvg) - normalVar / 2.0d0

        randField(:,:) = exp(randField(:,:) * sqrt(normalVar) + normalAvg)
    end subroutine transform_Lognormal

    !-----------------------------------------------------------------------------------------------
    ! Transform to Uniform Distribution: X ~ U(a, b)
    ! a = fieldAvg - sqrt(3 * fieldVar)
    ! b = fieldAvg + sqrt(3 * fieldVar)
    ! X = a + (b - a) * Phi(Z) = fieldAvg + sqrt(3 * fieldVar) * erf(Z / sqrt(2))
    !-----------------------------------------------------------------------------------------------
    subroutine transform_Uniform(fieldAvg, fieldVar, randField)
        double precision, intent(in) :: fieldAvg, fieldVar
        double precision, dimension(:,:), intent(inout) :: randField
        double precision :: half_width

        if (fieldVar <= 0.0d0) then
            call rosemblat_abort("transform_Uniform: fieldVar must be strictly positive (> 0).")
        end if

        half_width = sqrt(3.0d0 * fieldVar)

        ! Using standard Fortran erf(z / sqrt(2)): 2*Phi(Z) - 1 = erf(Z / sqrt(2))
        ! X = fieldAvg + half_width * erf(Z / sqrt(2))
        randField(:,:) = fieldAvg + half_width * erf(randField(:,:) / SQRT_2)
    end subroutine transform_Uniform

    !-----------------------------------------------------------------------------------------------
    ! Solve Weibull parameters (k = shape, lambda = scale) from mean and variance
    ! Mean = lambda * Gamma(1 + 1/k)
    ! Var  = lambda^2 * [ Gamma(1 + 2/k) - Gamma(1 + 1/k)^2 ]
    ! CV^2 = Var / Mean^2 = Gamma(1 + 2/k) / Gamma(1 + 1/k)^2 - 1
    !-----------------------------------------------------------------------------------------------
    subroutine solve_weibull_params(mean, var, k_shape, lambda_scale, status)
        double precision, intent(in) :: mean, var
        double precision, intent(out) :: k_shape, lambda_scale
        integer, intent(out) :: status
        double precision :: cv, cv2, target, k, g, dg, delta
        integer :: iter

        status = 0
        if (mean <= 0.0d0 .or. var <= 0.0d0) then
            status = -1
            return
        end if

        cv = sqrt(var) / mean
        cv2 = cv * cv
        target = log(1.0d0 + cv2)

        ! Menon's initial approximation
        k = cv**(-1.086d0)
        if (k < 0.01d0) k = 0.01d0
        if (k > 100.0d0) k = 100.0d0

        ! Newton-Raphson solver for g(k) = ln Gamma(1 + 2/k) - 2 ln Gamma(1 + 1/k) - target = 0
        do iter = 1, 100
            g = log_gamma(1.0d0 + 2.0d0/k) - 2.0d0 * log_gamma(1.0d0 + 1.0d0/k) - target
            dg = (2.0d0 / (k * k)) * (digamma(1.0d0 + 1.0d0/k) - digamma(1.0d0 + 2.0d0/k))

            if (abs(dg) < 1.0d-30) then
                status = -2
                return
            end if

            delta = g / dg

            ! Safeguard step to maintain k in positive domain
            if (k - delta <= 0.0001d0) then
                delta = 0.5d0 * k
            end if

            k = k - delta

            if (abs(delta) < 1.0d-13 * k .or. abs(g) < 1.0d-13) exit
        end do

        if (iter > 100) then
            status = -3
            return
        end if

        k_shape = k
        lambda_scale = mean / exp(log_gamma(1.0d0 + 1.0d0/k))
    end subroutine solve_weibull_params

    !-----------------------------------------------------------------------------------------------
    ! Transform to Weibull Distribution: X ~ Weibull(k, lambda)
    ! F_X(x) = 1 - exp(-(x/lambda)^k)
    ! F_X^-1(u) = lambda * [ -ln(1 - u) ]^(1/k)
    ! Using u = Phi(Z) => 1 - u = Phi(-Z) = 0.5 * erfc(Z / sqrt(2))
    ! X = lambda * [ -ln( 0.5 * erfc(Z / sqrt(2)) ) ]^(1/k)
    !-----------------------------------------------------------------------------------------------
    subroutine transform_Weibull(fieldAvg, fieldVar, randField)
        double precision, intent(in) :: fieldAvg, fieldVar
        double precision, dimension(:,:), intent(inout) :: randField
        double precision :: k_shape, lambda_scale, q_tail, inv_k
        integer :: status, i, j
        character(len=256) :: err_msg

        if (fieldAvg <= 0.0d0) then
            call rosemblat_abort("transform_Weibull: fieldAvg must be strictly positive (> 0).")
        end if
        if (fieldVar <= 0.0d0) then
            call rosemblat_abort("transform_Weibull: fieldVar must be strictly positive (> 0).")
        end if

        call solve_weibull_params(fieldAvg, fieldVar, k_shape, lambda_scale, status)
        if (status /= 0) then
            write(err_msg, '(A, I0, A, ES12.4, A, ES12.4)') &
                "transform_Weibull: failed to converge on parameters (status=", status, &
                ", fieldAvg=", fieldAvg, ", fieldVar=", fieldVar, ")"
            call rosemblat_abort(trim(err_msg))
        end if

        inv_k = 1.0d0 / k_shape

        do j = 1, size(randField, 2)
            do i = 1, size(randField, 1)
                ! q_tail = 1 - Phi(z) = Phi(-z) = 0.5 * erfc(z / sqrt(2))
                q_tail = 0.5d0 * erfc(randField(i, j) / SQRT_2)
                ! Bound tail probability to prevent log(0) or zero division
                if (q_tail < EPS_TAIL) q_tail = EPS_TAIL
                if (q_tail > ONE_MINUS_EPS) q_tail = ONE_MINUS_EPS

                randField(i, j) = lambda_scale * ((-log(q_tail))**inv_k)
            end do
        end do
    end subroutine transform_Weibull

    !-----------------------------------------------------------------------------------------------
    ! Solve Gamma parameters: shape alpha and rate lambda_rate
    ! Mean = alpha / lambda_rate => alpha = Mean^2 / Var
    ! Var  = alpha / lambda_rate^2 => lambda_rate = Mean / Var, scale beta = Var / Mean
    !-----------------------------------------------------------------------------------------------
    subroutine solve_gamma_params(mean, var, shape_alpha, scale_rate)
        double precision, intent(in) :: mean, var
        double precision, intent(out) :: shape_alpha, scale_rate

        if (mean <= 0.0d0 .or. var <= 0.0d0) then
            call rosemblat_abort("solve_gamma_params: mean and var must be strictly positive.")
        end if

        shape_alpha = (mean * mean) / var
        scale_rate  = mean / var
    end subroutine solve_gamma_params

    !-----------------------------------------------------------------------------------------------
    ! Transform to Gamma Distribution: X ~ Gamma(alpha, beta)
    ! Uses cdfgam from dcdflib to compute quantiles with machine accuracy.
    !-----------------------------------------------------------------------------------------------
    subroutine transform_Gamma(fieldAvg, fieldVar, randField)
        double precision, intent(in) :: fieldAvg, fieldVar
        double precision, dimension(:,:), intent(inout) :: randField
        double precision :: shape_alpha, scale_rate
        double precision :: p, q, x, bound, z
        integer :: which, status, i, j
        character(len=256) :: err_msg

        if (fieldAvg <= 0.0d0) then
            call rosemblat_abort("transform_Gamma: fieldAvg must be strictly positive (> 0).")
        end if
        if (fieldVar <= 0.0d0) then
            call rosemblat_abort("transform_Gamma: fieldVar must be strictly positive (> 0).")
        end if

        call solve_gamma_params(fieldAvg, fieldVar, shape_alpha, scale_rate)

        which = 2 ! Compute X from P, Q, SHAPE, SCALE in cdfgam

        do j = 1, size(randField, 2)
            do i = 1, size(randField, 1)
                z = randField(i, j)

                ! Evaluate standard normal CDF and CCDF
                p = 0.5d0 * erfc(-z / SQRT_2)
                q = 0.5d0 * erfc(z / SQRT_2)

                ! Ensure p and q are strictly in (0, 1) and sum to 1.0
                if (p < EPS_TAIL) then
                    p = EPS_TAIL
                    q = 1.0d0 - EPS_TAIL
                else if (q < EPS_TAIL) then
                    q = EPS_TAIL
                    p = 1.0d0 - EPS_TAIL
                end if

                call cdfgam(which, p, q, x, shape_alpha, scale_rate, status, bound)

                if (status /= 0) then
                    write(err_msg, '(A, I0, A, ES12.4, A, ES12.4)') &
                        "transform_Gamma: cdfgam returned error (status=", status, &
                        ", p=", p, ", z=", z, ")"
                    call rosemblat_abort(trim(err_msg))
                end if

                randField(i, j) = x
            end do
        end do
    end subroutine transform_Gamma

end module rosemblat
