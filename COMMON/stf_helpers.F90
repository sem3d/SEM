!! Shared module for Source Time Functions
!! Consolidates all math waveforms for both 2D and 3D solvers.
MODULE stf_helpers
    USE constants
    IMPLICIT none

CONTAINS

    real(fpp) function Gaussian_2D (time, tau, f0)
        implicit none
        real(fpp), intent(in) :: time, tau, f0
        real(fpp) :: sigma
        sigma = (M_PI * f0 * (time - tau))**2
        Gaussian_2D = (time - tau) * exp(-sigma)
    end function Gaussian_2D

    real(fpp) function Gaussian_3D (time, ts, tau)
        implicit none
        real(fpp), intent(in) :: time, ts, tau
        if ((time - ts) < 8.0_fpp * tau) then
            Gaussian_3D = -2.0_fpp * (time - ts) * exp(-(time - ts)**2 / tau**2)
        else
            Gaussian_3D = 0.0_fpp
        endif
    end function Gaussian_3D

    real(fpp) function Ricker_2D (time, tau, f0)
        implicit none
        real(fpp), intent(in) :: time, tau, f0
        real(fpp) :: sigma
        sigma = (M_PI * f0 * (time - tau))**2
        Ricker_2D = (1.0_fpp - 2.0_fpp * sigma) * exp(-sigma)
        if (sigma > 50.0_fpp) then
            Ricker_2D = 0.0_fpp
        endif
    end function Ricker_2D

    real(fpp) function Ricker_3D (time, tau, f0)
        implicit none
        real(fpp), intent(in) :: time, tau, f0
        real(fpp) :: sigma, alpha
        alpha = -1.0_fpp * M_PI**2 * f0**2
        if (time < 2.5_fpp * tau) then
            sigma = alpha * (time - tau)**2
            Ricker_3D = 2.0_fpp * alpha * (1.0_fpp + 2.0_fpp * sigma) * exp(sigma)
        else
            Ricker_3D = 0.0_fpp
        endif
    end function Ricker_3D

    real(fpp) function Gabor (time, tau, fp, gamma, ts)
        implicit none
        real(fpp), intent(in) :: time, tau, fp, gamma, ts
        real(fpp) :: sigma, xomega, xval1, xval2
        xomega = M_PI * 0.5_fpp
        if (time < 32.0_fpp) then
            sigma = 2.0_fpp * M_PI * fp * (time - ts)
            xval1 = cos(sigma + xomega)
            sigma = sigma / gamma
            sigma = sigma**2
            if (sigma < 100.0_fpp) then
                xval2 = exp(-sigma)
            else
                xval2 = 0.0_fpp
            endif
            Gabor = xval2 * xval1 * tau
        else
            Gabor = 0.0_fpp
        endif
    end function Gabor

    real(fpp) function Source_sinewave_math (time, ts, f0)
        implicit none
        real(fpp), intent(in) :: time, ts, f0
        Source_sinewave_math = sin(2.0_fpp * M_PI * f0 * (time - ts))
    end function Source_sinewave_math

    real(fpp) function Source_square_math (t, ts, dt, k)
        implicit none
        real(fpp), intent(in) :: t, ts, dt, k
        real(fpp) :: w0, winf
        w0 = (log(cosh(k * (-ts))) - log(cosh(k * (ts + dt)))) / k
        winf = dt
        Source_square_math = (tanh(k * (t - ts)) + tanh(k * (ts + dt - t))) / (winf - w0)
    end function Source_square_math

    real(fpp) function Source_tanh_math (time, ts, k)
        implicit none
        real(fpp), intent(in) :: time, ts, k
        Source_tanh_math = 0.5_fpp * (tanh(k * (time - ts)) + 1.0_fpp)
    end function Source_tanh_math

    real(fpp) function Source_Spice_Bench_math (time, ts, cutoff_freq, gamma)
        implicit none
        real(fpp), intent(in) :: time, ts, cutoff_freq, gamma
        real(fpp) :: T, k, s
        if (time < ts) then
            Source_Spice_Bench_math = 0.0_fpp
            return
        endif
        T = 1.0_fpp / cutoff_freq
        k = gamma
        if (k < 1.0_fpp) k = 1.0_fpp
        s = ((time - ts) / T)**k
        Source_Spice_Bench_math = (1.0_fpp - (1.0_fpp + s) * exp(-s))
    end function Source_Spice_Bench_math

    real(fpp) function Ricker_Fl (time, tau, f0)
        implicit none
        real(fpp), intent(in) :: time, tau, f0
        real(fpp) :: alpha, sigma, Ricker_val
        alpha = -1.0_fpp * M_PI**2 * f0**2
        if (time < 2.5_fpp * tau) then
            sigma = alpha * (time - tau)**2
            Ricker_val = 2.0_fpp * alpha * (1.0_fpp + 2.0_fpp * sigma) * exp(sigma)
            Ricker_Fl = 2.0_fpp * alpha * (time - tau) * Ricker_val + &
                        8.0_fpp * alpha * alpha * (time - tau) * exp(sigma)
        else
            Ricker_Fl = 0.0_fpp
        endif
    end function Ricker_Fl

    real(fpp) function Triangle (time, tau)
        implicit none
        real(fpp), intent(in) :: time, tau
        if (time < 0.005_fpp) then
            Triangle = -time * tau * 1.0e10_fpp
        elseif (time < 0.01_fpp) then
            Triangle = -(-time + 0.01_fpp) * tau * 1.0e10_fpp
        else
            Triangle = 0.0_fpp
        endif
    end function Triangle

    real(fpp) function HSF (time, tau)
        implicit none
        real(fpp), intent(in) :: time, tau
        if (time < 0.00000001_fpp) then
            HSF = 0.0_fpp
        endif
        if (time == 0.0_fpp) then
            HSF = -0.5_fpp * tau
        endif
        if (time > 0.00000001_fpp) then
            HSF = -1.0_fpp * tau
        endif
    end function HSF

    real(fpp) function DM (time, tau, Q, X, Y, L, v, d, a)
        implicit none
        real(fpp), intent(in) :: time, tau, Q, X, Y, L, v, d, a
        DM = Q * Y / 2.0_fpp * (X**(((v * (time - tau) - a)**2 / d**2)) + &
                                X**(((v * (time - tau) - a - L)**2 / d**2)))
    end function DM

    real(fpp) function Ormsby (time, tau, band)
        implicit none
        real(fpp), intent(in) :: time, tau
        real(fpp), dimension(0:3), intent(in) :: band
        real(fpp) :: t, f1, f2, f3, f4, den1, den2
        
        f1 = band(0)
        f2 = band(1)
        f3 = band(2)
        f4 = band(3)
        t = time - tau
        
        if (abs(t) .lt. 1.0e-12_fpp) then
            Ormsby = M_PI * (f4 + f3) - M_PI * (f2 + f1)
        else
            den1 = M_PI**2 * (f4 - f3) * t**2
            den2 = M_PI**2 * (f2 - f1) * t**2
            Ormsby = (sin(M_PI * f4 * t)**2 - sin(M_PI * f3 * t)**2) / den1 - &
                     (sin(M_PI * f2 * t)**2 - sin(M_PI * f1 * t)**2) / den2
        endif
    end function Ormsby

END MODULE stf_helpers
