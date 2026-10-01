module spectra_RF
    use displayCarvalhol
    use math_RF
    use write_Log_File
    use constants_RF
    use type_RF
    use special_functions

    implicit none
contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kMaxND(corrMod, kMax, hurst_nu_in, sphere_D_in);
        implicit none
        !INPUT
        integer, intent(in) :: corrMod;
        double precision, optional, intent(in) :: hurst_nu_in, sphere_D_in

        !OUTPUT
        double precision, dimension(:),   intent(out) :: kMax;

        !LOCAL VARIABLES
        integer          :: nDim
        double precision :: nu, sphere_D, a, val

        nDim = size(kMax)
        nu = 0.5d0
        if(present(hurst_nu_in)) nu = hurst_nu_in
        if(nu <= 0.0d0) nu = 0.5d0

        sphere_D = 1.0d0
        if(present(sphere_D_in)) sphere_D = sphere_D_in
        if(sphere_D <= 0.0d0) sphere_D = 1.0d0

        select case(corrMod)
            case(cm_GAUSSIAN)
                select case(nDim)
                    case(1)
                        kMax(:) = 6.457D0 !Value to cover 99% Spectra area
                    case(2)
                        kMax(:) = 7.035D0
                    case(3)
                        kMax(:) = 7.355D0
                end select

            case(cm_COS)
                kMax(:) = 12.0D0

            case(cm_KARMAN)
                select case(nDim)
                    case(1)
                        a = 2.0d0 * sqrt(PI) * exp(log_gamma(nu + 0.5d0) - log_gamma(nu))
                        val = a * sqrt(max(1.0d0, (0.01d0)**(-1.0d0/max(0.05d0, nu)) - 1.0d0))
                        kMax(:) = max(8.0d0, min(val, 80.0d0))
                    case(2)
                        a = 2.0d0 * sqrt(max(0.01d0, nu))
                        val = a * sqrt(max(1.0d0, (0.01d0)**(-1.0d0/max(0.05d0, nu)) - 1.0d0))
                        kMax(:) = max(8.0d0, min(val, 50.0d0))
                    case(3)
                        a = (6.0d0 * sqrt(PI) * exp(log_gamma(nu + 1.5d0) - log_gamma(nu)))**(1.0d0/3.0d0)
                        val = a * sqrt(max(1.0d0, (0.01d0)**(-1.0d0/max(0.05d0, nu + 0.5d0)) - 1.0d0))
                        kMax(:) = max(8.0d0, min(val, 40.0d0))
                end select

            case(cm_EXPONENTIAL)
                select case(nDim)
                    case(1)
                        kMax(:) = 50.0D0
                    case(2)
                        kMax(:) = 14.14D0
                    case(3)
                        kMax(:) = 14.0D0
                end select

            case(cm_POWERLAW)
                select case(nDim)
                    case(1)
                        kMax(:) = 10.43D0
                    case(2)
                        kMax(:) = 7.5D0
                    case(3)
                        kMax(:) = 6.13D0
                end select

            case(cm_LOWPASS)
                select case(nDim)
                    case(1)
                        kMax(:) = PI
                    case(2)
                        kMax(:) = 2.0D0
                    case(3)
                        kMax(:) = 1.5D0 * PI
                end select

            case(cm_TRIANGULAR)
                select case(nDim)
                    case(1)
                        kMax(:) = PI
                    case(2)
                        kMax(:) = 1.5D0 * PI
                    case(3)
                        kMax(:) = 2.0D0 * PI
                end select

            case(cm_MONODISPERSE)
                kMax(:) = 20.0D0 * PI / sphere_D

            case default
                kMax(:) = 10.0D0
        end select

    end subroutine set_kMaxND

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kPoints(RDF, xStep);
        implicit none

        !INPUT OUTPUT
        type(RF) :: RDF
        !INPUT
        double precision, dimension(:), intent(in) :: xStep

        !LOCAL
        integer :: i, j
        integer(kind=8) :: i_long, kNLocal, kN_unsigned
        double precision :: kAdjust    = 1.0D0 !"kNStep minimum" multiplier
        double precision :: periodMult = 1.1D0 !"range" multiplier
        double precision :: rAdjust    = 5.0D0 !"kNStep minimum" multiplier (isotropic)
        integer, dimension(RDF%nDim) :: kInitial, kFinal, kNStepLocal
        double precision, dimension(:,:) , allocatable :: kSign, kPoints_base;
        integer(kind=8) :: countInit, countEnd

        if(allocated(RDF%kPoints)) deallocate(RDF%kPoints)

        call set_kMaxND(RDF%corrMod, RDF%kMax, RDF%hurst_nu, RDF%sphere_D) !Defining kMax according to corrMod
        call wLog(" RDF%kMax = ")
        call wLog(RDF%kMax)

        select case (RDF%method)
            case(ISOTROPIC)
                RDF%kDelta(1) = 2.0D0*PI/(periodMult*sqrt(sum(RDF%xRange**2))) !Diagonal
                RDF%kNStep(1) = 1 + int(rAdjust*(ceiling(maxval(RDF%kMax)/RDF%kDelta(1)))); !Number of points in k
                RDF%kDelta(1) = maxval(RDF%kMax)/(RDF%kNStep(1)-1); !Redefining kDelta after ceiling and adjust
                RDF%kNTotal   = RDF%kNStep(1);

                allocate(RDF%kPoints(1, RDF%kNTotal))

                do i_long = 1, RDF%kNTotal
                    call get_Permutation(i, [RDF%kMax(1)], [RDF%kNStep(1)], RDF%kPoints(:, i), snapExtremes = .true.);
                end do

            case(SHINOZUKA, RANDOMIZATION)
                RDF%kDelta(:) = 2.0D0*PI/(periodMult*RDF%xRange)
                !write(*,*)" RDF%kDelta = ", RDF%kDelta
                RDF%kNStep(:)   = 1 + int(kAdjust*(ceiling(RDF%kMax/RDF%kDelta(:)))); !Number of points in k
                !write(*,*)" RDF%kNStep = ", RDF%kNStep
                RDF%kDelta(:) = (RDF%kMax)/(RDF%kNStep-1); !Redefining kDelta after ceiling and adjust
                write(*,*)" RDF%kDelta = ", RDF%kDelta

                allocate(kSign (2**(RDF%nDim-1), RDF%nDim));
                call set_kSign(kSign) !Set the sign permutations for kVec

                !write(*,*)" kSign = ", kSign

                kN_unsigned = product(int(RDF%kNStep,8));
                RDF%kNTotal = kN_unsigned*int(size(kSign,1),8)

                call wLog(" RDF%kNStep = ")
                call wLog(RDF%kNStep)
                call wLog(" RDF%kDelta = ")
                call wLog(RDF%kDelta)
                call wLog(" RDF%kNTotal = ")
                call wLog(RDF%kNTotal)

                write(*,*)" kN_unsigned = ", kN_unsigned
                write(*,*)" RDF%kNTotal = ", RDF%kNTotal


                allocate(RDF%kPoints (RDF%nDim, RDF%kNTotal));
                allocate(kPoints_base(RDF%nDim, kN_unsigned));

                if(RDF%method == SHINOZUKA) then
                    call setGrid(kPoints_base(:,:), dble(kInitial)*0d0, RDF%kDelta, RDF%kNStep)
                else if (RDF%method == RANDOMIZATION) then
                    call random_number(kPoints_base(:,:))
                    do i = 1, RDF%nDim
                        kPoints_base(i,:) = kPoints_base(i,:) * RDF%kMax(i)
                    end do
                end if

                !Copy kPoints changing the sign
                do i = 1, size(kSign,1)
                    countInit = (i-1)*kN_unsigned + 1
                    countEnd  = countInit + kN_unsigned - 1
                    do j = 1, RDF%nDim
                        RDF%kPoints(j,countInit:countEnd) = kPoints_base(j,:)*kSign(i,j)
                    end do
                end do

            case(FFT)
                RDF%kNStep(:) = find_xNStep(xMaxExt = RDF%xRange, xStep=xStep)
                RDF%kDelta(:) = 2.0D0*PI/(periodMult*RDF%xRange)
                RDF%kMax(:)   = ((RDF%kNStep(:) - 1) * RDF%kDelta(:))/1.0D0
                RDF%kNTotal   = product(int(RDF%kNStep,8));


                    kInitial = 0
                    kInitial(RDF%nDim)=RDF%kNInit-1 !-1 so k starts at 0
                    kFinal = 0
                    kFinal(RDF%nDim)= RDF%kNInit + RDF%kNEnd -1
                    call wLog("RDF%kNInit = ")
                    call wLog(RDF%kNInit)
                    call wLog("RDF%kNEnd = ")
                    call wLog(RDF%kNEnd)
                    kNStepLocal = RDF%kNStep
                    kNStepLocal(RDF%nDim)=RDF%kNEnd - RDF%kNInit + 1
                    kNLocal = product(kNStepLocal)
                    call wLog("kNLocal = ")
                    call wLog(kNLocal)
                    call wLog("HERE !!!!!!!!!!!!")
                    allocate(RDF%kPoints(RDF%nDim, kNLocal))
                    call wLog("shape(RDF%kPoints) = ")
                    call wLog(shape(RDF%kPoints))
                    call wLog("Making kPoints= ")
                    call setGrid(RDF%kPoints, dble(kInitial)*RDF%kDelta, RDF%kDelta, kNStepLocal)


        end select

        if(allocated(kSign)) deallocate(kSign)
        if(allocated(kPoints_base)) deallocate(kPoints_base)

    end subroutine set_kPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_SkVec(RDF, corrL_in);
        implicit none

        !INPUT OUTPUT
        type(RF) :: RDF
        double precision, dimension(:), optional, intent(in) :: corrL_in
        !LOCAL
        integer :: i, j, freqK = 6
        double precision, dimension(RDF%nDim) :: corrL
        double precision, dimension(:), allocatable :: kk2, z_k
        double precision :: nu, a, prefactor, prod_corrL
        double precision :: eta_p, D_p
        logical :: use_vw

        corrL(:) = 1.0D0
        if(present(corrL_in)) corrL = corrL_in
        write(*,*) 'correlation length',corrL
        write(*,*) "size(RDF%kPoints,2) = ", size(RDF%kPoints,2)

        if(allocated(RDF%SkVec)) deallocate(RDF%SkVec)
        allocate(RDF%SkVec(size(RDF%kPoints,2)))

        prod_corrL = product(corrL(1:RDF%nDim))

        allocate(kk2(size(RDF%SkVec,1)))
        allocate(z_k(size(RDF%SkVec,1)))
        kk2(:) = 0.0d0
        do i = 1, size(RDF%kPoints, 1)
            kk2(:) = kk2(:) + (RDF%kPoints(i,:) * corrL(min(i, RDF%nDim)))**2.0D0
        end do
        z_k(:) = sqrt(kk2(:))

        call wLog("RDF%corrMod")
        call wLog(RDF%corrMod)

        select case(RDF%corrMod)

            case(cm_GAUSSIAN)
                if(RDF%rang == 0) write(*,*) "Gaussian Correlation Model"
                call wLog("cm_GAUSSIAN")
                RDF%SkVec(:) = prod_corrL * exp(-(z_k(:)**2.0d0) / (4.0d0 * PI))

            case(cm_COS)
                call wLog("cm_COS")
                RDF%SkVec = 0.0d0
                if(size(RDF%SkVec) > freqK) then
                    RDF%SkVec(freqK) = 1.0d0
                    call wLog("kPoint used in cosinus CM: ")
                    call wLog(RDF%kPoints(:,freqK))
                    if(RDF%rang == 0) write(*,*) "kPoint used in cosinus CM: ", RDF%kPoints(:,freqK)
                else
                    stop ("When using the cosinus correlation Model you should have more than 'freqK' kPoints")
                end if

            case(cm_KARMAN)
                call wLog("cm_KARMAN")
                nu = RDF%hurst_nu
                if(nu <= 0.0d0) nu = 0.5d0

                select case(RDF%nDim)
                    case(1)
                        a = 2.0d0 * sqrt(PI) * exp(log_gamma(nu + 0.5d0) - log_gamma(nu))
                        prefactor = corrL(1) / (2.0d0 * PI)
                        RDF%SkVec(:) = prefactor / ((1.0d0 + (z_k(:)/a)**2.0D0)**(nu + 0.5d0))
                    case(2)
                        a = 2.0d0 * sqrt(nu)
                        prefactor = prod_corrL / (4.0d0 * PI)
                        RDF%SkVec(:) = prefactor / ((1.0d0 + (z_k(:)/a)**2.0D0)**(nu + 1.0d0))
                    case(3)
                        a = (6.0d0 * sqrt(PI) * exp(log_gamma(nu + 1.5d0) - log_gamma(nu)))**(1.0d0/3.0d0)
                        prefactor = prod_corrL / (6.0d0 * PI**2)
                        RDF%SkVec(:) = prefactor / ((1.0d0 + (z_k(:)/a)**2.0D0)**(nu + 1.5d0))
                end select

            case(cm_EXPONENTIAL)
                call wLog("cm_EXPONENTIAL")
                select case(RDF%nDim)
                    case(1)
                        a = 2.0d0
                        prefactor = corrL(1) / (2.0d0 * PI)
                        RDF%SkVec(:) = prefactor / (1.0d0 + (z_k(:)/a)**2.0D0)
                    case(2)
                        a = sqrt(2.0d0)
                        prefactor = prod_corrL / (4.0d0 * PI)
                        RDF%SkVec(:) = prefactor / ((1.0d0 + (z_k(:)/a)**2.0D0)**1.5d0)
                    case(3)
                        a = 6.0d0**(1.0d0/3.0d0)
                        prefactor = prod_corrL / (6.0d0 * PI**2)
                        RDF%SkVec(:) = prefactor / ((1.0d0 + (z_k(:)/a)**2.0D0)**2.0d0)
                end select

            case(cm_POWERLAW)
                call wLog("cm_POWERLAW")
                select case(RDF%nDim)
                    case(1)
                        a = PI / 2.0d0
                        prefactor = corrL(1) / (2.0d0 * PI)
                        RDF%SkVec(:) = prefactor * (1.0d0 + z_k(:)/a) * exp(-z_k(:)/a)
                    case(2)
                        prefactor = prod_corrL / (4.0d0 * PI)
                        do j = 1, size(RDF%SkVec, 1)
                            if (z_k(j) < 1.0d-8) then
                                RDF%SkVec(j) = prefactor
                            else
                                RDF%SkVec(j) = prefactor * z_k(j) * BESSK1(z_k(j))
                            end if
                        end do
                    case(3)
                        a = (3.0d0 * PI / 4.0d0)**(1.0d0/3.0d0)
                        prefactor = prod_corrL / (6.0d0 * PI**2)
                        RDF%SkVec(:) = prefactor * exp(-z_k(:)/a)
                end select

            case(cm_LOWPASS)
                call wLog("cm_LOWPASS")
                RDF%SkVec(:) = 0.0d0
                select case(RDF%nDim)
                    case(1)
                        a = PI
                        prefactor = corrL(1) / (2.0d0 * PI)
                        where (z_k(:) <= a) RDF%SkVec(:) = prefactor
                    case(2)
                        a = 2.0d0
                        prefactor = prod_corrL / (2.0d0 * PI**2)
                        where (z_k(:) <= a) RDF%SkVec(:) = prefactor
                    case(3)
                        a = 1.5d0 * PI
                        prefactor = (2.0d0 * prod_corrL) / (9.0d0 * PI**4)
                        where (z_k(:) <= a) RDF%SkVec(:) = prefactor
                end select

            case(cm_TRIANGULAR)
                call wLog("cm_TRIANGULAR")
                RDF%SkVec(:) = 0.0d0
                select case(RDF%nDim)
                    case(1)
                        a = PI
                        prefactor = corrL(1) / (2.0d0 * PI)
                        where (z_k(:) <= a) RDF%SkVec(:) = prefactor * (1.0d0 - z_k(:)/a)
                    case(2)
                        a = 1.5d0 * PI
                        prefactor = prod_corrL / (2.0d0 * PI**2)
                        where (z_k(:) <= a) RDF%SkVec(:) = prefactor * (1.0d0 - z_k(:)/a)
                    case(3)
                        a = 2.0d0 * PI
                        prefactor = (3.0d0 * prod_corrL) / (8.0d0 * PI**4)
                        where (z_k(:) <= a) RDF%SkVec(:) = prefactor * (1.0d0 - z_k(:)/a)
                end select

            case(cm_MONODISPERSE)
                call wLog("cm_MONODISPERSE")
                eta_p = RDF%sphere_eta
                if (eta_p <= 0.0d0 .or. eta_p >= 0.74d0) eta_p = 0.2d0
                D_p = RDF%sphere_D
                if (D_p <= 0.0d0) D_p = 1.0d0
                use_vw = RDF%use_verlet_weis

                do j = 1, size(RDF%SkVec, 1)
                    RDF%SkVec(j) = eval_monodisperse_psd(z_k(j), RDF%nDim, eta_p, D_p, use_vw)
                end do

            case default
                write(*,*) "ERROR - Unknown correlation model ID in set_SkVec: ", RDF%corrMod
                stop ("Unknown correlation model ID in set_SkVec")

        end select

        deallocate(kk2)
        deallocate(z_k)

    end subroutine set_SkVec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    ! Helper: Evaluates Power Spectral Density for Monodisperse Spheres (Percus-Yevick + Verlet-Weis)
    !-----------------------------------------------------------------------------------------------
    function eval_monodisperse_psd(k_in, nDim, eta_in, D_in, use_vw_in) result(phi_out)
        implicit none
        double precision, intent(in) :: k_in, eta_in, D_in
        integer,          intent(in) :: nDim
        logical,          intent(in) :: use_vw_in
        double precision             :: phi_out

        ! Local variables
        double precision :: k_val, eta, D, R, u, S_py, S_tot, m_k, rho
        double precision :: eta_w, D_w, R_w, rho_w, l1, l2, c_val
        double precision :: g_cs, g_py, A, mu, denom, I_cos, I_sin, dh
        double precision :: j1_2u, c2d, j1_val

        k_val = abs(k_in)
        eta = eta_in
        D = D_in
        R = D / 2.0d0

        select case(nDim)
            case(1)
                rho = eta / D
                u = k_val * D / 2.0d0
                if (u < 1.0d-6) then
                    m_k = D
                else
                    m_k = D * sin(u) / u
                end if
                S_tot = (1.0d0 - eta**2) / max(1.0d-8, (1.0d0 - 2.0d0*eta*cos(k_val*D) + eta**2))
                phi_out = rho * (m_k**2) * S_tot / (2.0d0 * PI * eta * (1.0d0 - eta))

            case(2)
                rho = eta / (PI * R**2)
                u = k_val * R
                if (u < 1.0d-6) then
                    m_k = PI * R**2
                    c2d = -PI * (R**2) * (1.0d0 + 2.0d0*eta) / ((1.0d0 - eta)**3)
                else
                    j1_val = bessel_j1(u)
                    m_k = 2.0d0 * PI * (R**2) * j1_val / u
                    j1_2u = bessel_j1(2.0d0 * u)
                    c2d = -2.0d0 * PI * (R**2) * (j1_2u / u) * &
                          ((1.0d0 + eta)/((1.0d0 - eta)**3) - (3.0d0*eta/((1.0d0 - eta)**3)) * (j1_2u / (2.0d0*u)))
                end if
                S_tot = max(1.0d-8, 1.0d0 / (1.0d0 - rho * c2d))
                phi_out = rho * (m_k**2) * S_tot / ((2.0d0 * PI)**2 * eta * (1.0d0 - eta))

            case(3)
                rho = 3.0d0 * eta / (4.0d0 * PI * R**3)
                if (use_vw_in .and. eta > 0.05d0) then
                    eta_w = eta * (1.0d0 - eta / 16.0d0)
                    D_w   = D * (eta_w / eta)**(1.0d0 / 3.0d0)
                    R_w   = D_w / 2.0d0
                    rho_w = 3.0d0 * eta_w / (4.0d0 * PI * R_w**3)
                else
                    eta_w = eta
                    D_w   = D
                    R_w   = R
                    rho_w = rho
                end if

                l1 = (1.0d0 + 2.0d0*eta_w)**2 / ((1.0d0 - eta_w)**4)
                l2 = -(1.0d0 + eta_w/2.0d0)**2 / ((1.0d0 - eta_w)**4)

                u = k_val * R_w
                if (u < 1.0d-6) then
                    c_val = -8.0d0 * PI * (R_w**3) / 3.0d0 * ((4.0d0 + eta_w)*l1 + 18.0d0*eta_w*l2)
                else
                    c_val = -4.0d0 * PI * (R_w**3) / (u**3) * ( &
                        l1 * (sin(2.0d0*u) - 2.0d0*u*cos(2.0d0*u)) + &
                        3.0d0*eta_w*l2/u * (4.0d0*u*sin(2.0d0*u) + (2.0d0 - 4.0d0*u**2)*cos(2.0d0*u) - 2.0d0) + &
                        eta_w*l1/(2.0d0*u**3) * ((6.0d0*u**2 - 3.0d0 - 2.0d0*u**4)*cos(2.0d0*u) + &
                        (4.0d0*u**3 - 6.0d0*u)*sin(2.0d0*u) + 3.0d0) )
                end if
                S_py = 1.0d0 / max(1.0d-8, 1.0d0 - rho_w * c_val)

                if (use_vw_in .and. eta > 0.05d0) then
                    g_cs = (1.0d0 - eta/2.0d0) / ((1.0d0 - eta)**3)
                    g_py = (1.0d0 + eta_w/2.0d0) / ((1.0d0 - eta_w)**2)
                    A    = D * (g_cs - g_py)
                    mu   = (24.0d0 * (A / D) / (eta_w * g_py)) / D

                    denom = 4.0d0 * mu**4 + k_val**4
                    I_cos = mu * (2.0d0 * mu**2 + k_val**2) / denom
                    I_sin = k_val * (2.0d0 * mu**2 - k_val**2) / denom

                    if (k_val < 1.0d-8) then
                        dh = 4.0d0 * PI * A * (D * I_cos + I_sin)
                    else
                        dh = (4.0d0 * PI * A / k_val) * (sin(k_val * D) * I_cos + cos(k_val * D) * I_sin)
                    end if
                    S_tot = max(1.0d-8, S_py + rho * dh)
                else
                    S_tot = S_py
                end if

                ! Form factor of physical sphere
                u = k_val * R
                if (u < 1.0d-6) then
                    m_k = 4.0d0 * PI * (R**3) / 3.0d0
                else
                    m_k = 4.0d0 * PI * (R**3) * (sin(u) - u * cos(u)) / (u**3)
                end if

                phi_out = rho * (m_k**2) * S_tot / ((2.0d0 * PI)**3 * eta * (1.0d0 - eta))

            case default
                phi_out = 0.0d0
        end select

    end function eval_monodisperse_psd

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_rMax(corrMod, rMax, corrL);
        implicit none
        !INPUT
        integer,                intent(in) :: corrMod;
        double precision,   dimension(:),  intent(in), optional :: corrL;

        !OUTPUT
        double precision, dimension(:),   intent(out) :: rMax;

        !LOCAL VARIABLES
        double precision :: pi = 3.1415926535898
        double precision, dimension(:), allocatable:: corrL_effec

        allocate(corrL_effec (size(rMax)))

        if (present(corrL)) corrL_effec = corrL
        if (.not. present(corrL)) corrL_effec = 1
        !        do i = 1, 100
        !            kMax = i/10.0 * corrL(:)
        !            write(*,*) "kMax = ", kMax
        !            write(*,*) "Spectrum = ", get_SpectrumND(kMax, corrMod, corrL)
        !            call DispCarvalhol (kMax, "kMax")
        !            call DispCarvalhol (get_SpectrumND(kMax, corrMod, corrL), "Spectrum")
        !        end do

        select case(corrMod)
        case(cm_GAUSSIAN)
            rMax(:) = 2*pi*corrL_effec(:); !CRITERIA STILL TO BE TESTED
        case default
            rMax(:) = 4.0d0 * pi * corrL_effec(:)
        end select

        deallocate(corrL_effec)

    end subroutine set_rMax

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kSign(kSign);
        implicit none
        !INPUT and OUTPUT
        double precision,   dimension(:, :), intent(inout) :: kSign;

        !LOCAL VARIABLES
        integer:: i, j, pos, seedStep, nDim;

        nDim = size(kSign, 2);

        kSign(:,:) = 1d0

        do i = 2, nDim
            do j = 1, size(kSign,1)
                seedStep = 2**(nDim-i);
                pos      = cyclicMod(int((j-0.9)/seedStep), 2)
                if (mod(pos,2) == 1) kSign(j,i) = -1d0
            end do
        end do
    end subroutine set_kSign

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kArray_rand(corrMod, kArray_rand, nDim);
        implicit none
        !INPUT
        integer, intent(in) :: corrMod;
        integer :: nDim

        !OUTPUT
        double precision, dimension(:),  intent(out) :: kArray_rand;

        !LOCAL VARIABLES
        integer :: i
        double precision ::  bound,mean,cdf_x,q,sd,uniRand
        integer :: st, which
        !real :: av = 0.0
        !real :: gennor
        !real :: std = 1.0
        !real :: normRand

!        mean = 0.0D0
!        sd = 1.0D0
!        test2 = gennor (real(mean), real(sd))
!        write(*,*) "gennor (0.0, 1.0) = ", test2

!        call random_number(kArray_rand)
!        kArray_rand = kMaxScal*kArray_rand

!TODO RANDOMIZATION
        which = 2 !Find x, for a given cdf_x

        !write(*,*) "corrMod = ", corrMod

        select case(corrMod)
            case(cm_GAUSSIAN)
                !write(*,*) "Gaussian kArray"
                mean = 0.0D0
                sd = SQRT_2PI
                !sd = 1.0D0/SQRT_2PI
                do i = 1, size(kArray_rand)
                    call random_number(uniRand)
                    cdf_x = (0.5D0 + uniRand/2.0D0)
                    q = 1-cdf_x
                    call cdfnor(which,cdf_x,q,kArray_rand(i),mean,sd,st,bound) !Normal CDF
                    kArray_rand = kArray_rand**(dble(nDim))
                end do

!                sd = SQRT_2PI
!                !sd = 1.0D0/SQRT_2PI
!                do i = 1, size(kArray_rand)
!                    cdf_x = 1.0D0
!                    do j = 1, nDim
!                        call random_number(uniRand)
!                        cdf_x = cdf_x*(0.5D0 + uniRand/2.0D0)
!                    end do
!                    q = 1-cdf_x
!                    call cdfnor(which,cdf_x,q,kArray_rand(i),mean,sd,st,bound) !Normal CDF
!                end do


                !write(*,*) "gennor (mean, sd) = ",gennor (0.5, 1.0)
                ! PROBLEM, how to put the seed
!                mean = 0.0D0
!                sd = 1.0D0/SQRT_2PI
!
!                kArray_rand(:) = 1
!                do i = 1, size(kArray_rand)
!                    do j = 1, nDim
!                        normRand = gennor (real(mean), real(sd))
!                        kArray_rand(i) = kArray_rand(i)*abs(dble(normRand))
!                    end do
!                end do
        end select



        !call reorder_vector(kArray_rand)

        !kArray_rand = kArray_rand**(dble(nDim))

        !call dispCarvalhol(kArray, "kArray")

    end subroutine set_kArray_rand

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kDelta_rand(kArray_rand, kDelta_rand);
        implicit none
        !INPUT
        double precision, dimension(:),  intent(in) :: kArray_rand;

        !OUTPUT
        double precision, dimension(:),  intent(out) :: kDelta_rand;

        !LOCAL VARIABLES
        integer :: i
        double precision :: inc

        kDelta_rand = 0.0D0

        !Correction for the first term
        inc = (kArray_rand(1) - 0.0D0)
        kDelta_rand(1) = kDelta_rand(1) + inc

        !Calculating distance between numbers (we supose they're ordered)
        do i = 1, size(kArray_rand)-1
            inc = (kArray_rand(i+1) - kArray_rand(i))/2.0D0
            kDelta_rand(i)   = kDelta_rand(i)   + inc
            kDelta_rand(i+1) = kDelta_rand(i+1) + inc
        end do

    end subroutine set_kDelta_rand
!***********************************************************
!* This program tests the subroutine BESSK to calculate the*
!* modified Bessel function of the third kind of order N   *
!* for any real positive argument X.                       *
!* ------------------------------------------------------- *
    function BESSK(N,X)
        implicit none
        integer :: N,J
        double precision :: X,BESSK,TOX,BK,BKM,BKP
! ------------------------------------------------------------------------
!     CE SOUS-PROGRAMME CALCULE LA FONCTION BESSEL MODIFIFIEE 3E ESPECE
!     D'ORDRE N ENTIER POUR TOUT X REEL POSITIF > 0.  ON UTILISE ICI LA
!     FORMULE DE RECURRENCE CLASSIQUE EN PARTANT DE BESSK0 ET BESSK1.
!
!     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD
!     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
!     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
! ------------------------------------------------------------------------ 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
! ------------------------------------------------------------------------
        IF (N.EQ.0) THEN
            write(*,*) 'DEBUG OK'
            BESSK = BESSK0(X)
            RETURN
        ENDIF
        IF (N.EQ.1) THEN
            BESSK = BESSK1(X)
            RETURN
        ENDIF
        IF (X.EQ.0.D0) THEN
            BESSK = 1.D30
            RETURN
        ENDIF
        TOX = 2.D0/X
        BK  = BESSK1(X)
        BKM = BESSK0(X)
        DO 11 J=1,N-1
        BKP = BKM+DFLOAT(J)*TOX*BK
        BKM = BK
        BK  = BKP
     11 CONTINUE
        BESSK = BK
        RETURN
    END
! ----------------------------------------------------------------------
    FUNCTION BESSK0(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
!     POUR TOUT X REEL NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
        implicit none
        double precision :: X,BESSK0,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
        DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0, &
        0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
        DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1, & 
        -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/

        IF(X.EQ.0.D0) THEN
            BESSK0=1.D30
            RETURN
        ENDIF
        IF(X.LE.2.D0) THEN
            Y=X*X/4.D0
            AX=-LOG(X/2.D0)*BESSI0(X)
            BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
        ELSE
            Y=(2.D0/X)
            AX=EXP(-X)/SQRT(X)
            BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
        ENDIF
        RETURN
    END
! ----------------------------------------------------------------------
    FUNCTION BESSK1(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
!     POUR TOUT X REEL POSITF NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
        implicit none
        double precision :: X,BESSK1,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
        DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,  &
        -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
        DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1, &
        0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/

        IF(X.EQ.0.D0) THEN
            BESSK1=1.D32
            RETURN
        ENDIF
        IF(X.LE.2.D0) THEN
            Y=X*X/4.D0
            AX=LOG(X/2.D0)*BESSI1(X)
            BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
        ELSE
            Y=(2.D0/X)
            AX=EXP(-X)/DSQRT(X)
            BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
        ENDIF
        RETURN
    END
!
!     Bessel Function of the 1st kind of order zero.
!
    FUNCTION BESSI0(X)
        IMPLICIT NONE
        double precision :: X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
        DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
        0.2659732D0,0.360768D-1,0.45813D-2/
        DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,  &
        0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
        0.2635537D-1,-0.1647633D-1,0.392377D-2/

        IF(ABS(X).LT.3.75D0) THEN
            Y=(X/3.75D0)**2
            BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
        ELSE
            AX=ABS(X)
            Y=3.75D0/AX
            BX=EXP(AX)/DSQRT(AX)
            AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
            BESSI0=AX*BX
        ENDIF
        RETURN
    END
!
!     Bessel Function of the 1st kind of order one.
!
    FUNCTION BESSI1(X)
        IMPLICIT NONE
        double precision :: X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
        DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
        0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
        DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
        -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,        &
        -0.2895312D-1,0.1787654D-1,-0.420059D-2/

        IF(ABS(X).LT.3.75D0) THEN
            Y=(X/3.75D0)**2
            BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
        ELSE
            AX=ABS(X)
            Y=3.75D0/AX
            BX=EXP(AX)/DSQRT(AX)
            AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
            BESSI1=AX*BX
        ENDIF
        RETURN
      END

end module spectra_RF
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
