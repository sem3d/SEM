%% =========================================================================
% TEST_CONVERGENCE_RATES.M
% Convergence Study of the Empirical Autocorrelation via the Law of Large Numbers
%
% Checks whether the RMSE of the autocorrelation function decays proportionally to:
%                   RMSE ~ 1 / sqrt(N)
% as predicted by the Central Limit Theorem and the Ergodic Theorem.
%
% Tests in 1D, 2D and 3D by varying the sample size / number of realizations
% by factors of 4 (e.g., N/4, N, 4*N or M = 1, 4, 16, 64 realizations).
% =========================================================================

clear; close all; clc;

fprintf('=========================================================================\n');
fprintf('  CONVERGENCE STUDY OF THE AUTOCORRELATION ESTIMATOR (1D, 2D, 3D)       \n');
fprintf('  Validation of the asymptotic rate: RMSE ~ 1 / sqrt(N) (CLT)           \n');
fprintf('=========================================================================\n\n');

%% 1. 1D CONVERGENCE: Varying the Domain Length L (or Number of Points N)
fprintf('--- 1. 1D CONVERGENCE (Exponential Model, \\ell = 2.0) ---\n');
model_1d = 'exponential';
corrL_1d = 2.0;
extra_1d.nu = 0.5;

% Domain size factors: L = 20, 80, 320 (factors of 4!)
% N = 1024, 4096, 16384 (same spacing dx = 20/1024)
dx1 = 20.0 / 1024;
factors_1d = [0.25, 1.0, 4.0, 16.0];
N_base_1d = 2048;
N_list_1d = round(N_base_1d * factors_1d);
rmse_1d = zeros(size(N_list_1d));
r2_1d   = zeros(size(N_list_1d));

% Fix the seed for reproducibility
rng(42);

for idx = 1:length(N_list_1d)
    N = N_list_1d(idx);
    L = N * dx1;
    
    % Average of independent runs to smooth the variance of the error estimator
    n_runs = 5;
    rmse_runs = zeros(n_runs, 1);
    
    for r = 1:n_runs
        u = generate_rf_1d(N, L, model_1d, corrL_1d, extra_1d);
        max_lag = round(N / 8);
        [r_lags, R_emp] = measure_corr_1d(u, dx1, max_lag);
        R_theo = get_theo_corr(r_lags, model_1d, corrL_1d, 1, extra_1d);
        rmse_runs(r) = sqrt(mean((R_emp(:) - R_theo(:)).^2));
    end
    rmse_1d(idx) = mean(rmse_runs);
    fprintf('  N = %5d (factor %5.2f): RMSE = %.4e\n', N, factors_1d(idx), rmse_1d(idx));
end

% Log-log decay rate:
p_1d = polyfit(log(factors_1d), log(rmse_1d), 1);
slope_1d = p_1d(1);
fprintf('  >> Observed 1D convergence rate: N^(%.3f) [Theory: N^(-0.500)]\n\n', slope_1d);


%% 2. 2D CONVERGENCE: Varying the Sampled Area (L^2 with factors of 4)
fprintf('--- 2. 2D CONVERGENCE (von Karman Model, \\ell = 2.5) ---\n');
model_2d = 'karman';
corrL_2d = 2.5;
extra_2d.nu = 0.6;

% N = 64, 128, 256 (Areas N^2 = 4096, 16384, 65536 -> factors of 4 in area!)
N_list_2d = [64, 128, 256];
area_factors_2d = (N_list_2d / 128).^2; % [0.25, 1.0, 4.0]
dx2 = 20.0 / 128;
rmse_2d = zeros(size(N_list_2d));

for idx = 1:length(N_list_2d)
    N = N_list_2d(idx);
    L = N * dx2;
    
    n_runs = 4;
    rmse_runs = zeros(n_runs, 1);
    for r = 1:n_runs
        u = generate_rf_2d(N, L, model_2d, corrL_2d, extra_2d);
        max_lag = round(N / 4);
        [r_lags, R_radial] = measure_corr_2d(u, dx2, dx2, max_lag);
        R_theo = get_theo_corr(r_lags, model_2d, corrL_2d, 2, extra_2d);
        rmse_runs(r) = sqrt(mean((R_radial(:) - R_theo(:)).^2));
    end
    rmse_2d(idx) = mean(rmse_runs);
    fprintf('  N = %3dx%3d (Area x %4.2f): RMSE = %.4e\n', N, N, area_factors_2d(idx), rmse_2d(idx));
end

p_2d = polyfit(log(area_factors_2d), log(rmse_2d), 1);
slope_2d = p_2d(1);
fprintf('  >> 2D convergence rate (vs Area): Area^(%.3f) [Theory: Area^(-0.500)]\n\n', slope_2d);


%% 3. 3D CONVERGENCE: Ensemble Average / Volume (M = 1, 4, 16 realizations)
fprintf('--- 3. 3D CONVERGENCE (Monodisperse Model, diameter D = 1.5, \\eta = 0.30) ---\n');
model_3d = 'monodisperse';
corrL_3d = 2.0;
extra_3d.D = 1.5;
extra_3d.eta = 0.30;
extra_3d.use_vw = true;

N3 = 48;
L3 = 16.0;
dx3 = L3 / N3;
max_lag_3d = round(N3 / 4);

% By the Law of Large Numbers, the mean of M independent realizations
% should have a standard deviation that decays exactly as 1 / sqrt(M)
M_list = [1, 4, 16];
rmse_3d = zeros(size(M_list));

% Generate a pool of 16 3D realizations
n_total_realiz = 16;
fprintf('  Generating a pool of %d 3D realizations of size %dx%dx%d...\n', n_total_realiz, N3, N3, N3);
R_cube_pool = cell(n_total_realiz, 1);

for m = 1:n_total_realiz
    u = generate_rf_3d(N3, L3, model_3d, corrL_3d, extra_3d);
    [r3d, R_sph] = measure_corr_3d(u, dx3, dx3, dx3, max_lag_3d);
    R_cube_pool{m} = R_sph(:);
end

R3d_theo = get_theo_corr(r3d, model_3d, corrL_3d, 3, extra_3d);
R3d_theo = R3d_theo(:);

for idx = 1:length(M_list)
    M = M_list(idx);
    % Combine the first M realizations
    R_mean = zeros(size(R3d_theo));
    for m = 1:M
        R_mean = R_mean + R_cube_pool{m};
    end
    R_mean = R_mean / M;
    
    rmse_3d(idx) = sqrt(mean((R_mean - R3d_theo).^2));
    fprintf('  Ensemble M = %2d realizations (factor x%2d): RMSE = %.4e\n', M, M, rmse_3d(idx));
end

p_3d = polyfit(log(M_list), log(rmse_3d), 1);
slope_3d = p_3d(1);
fprintf('  >> 3D convergence rate (vs M): M^(%.3f) [Theory: M^(-0.500)]\n\n', slope_3d);


%% 4. CONVERGENCE PLOT AND COMPARISON WITH 1/sqrt(N)
fig = figure('Name', 'Validation of the Law of Large Numbers (1/sqrt(N))', ...
             'Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.45]);

% --- Subplot 1: 1D ---
subplot(1, 3, 1);
loglog(factors_1d, rmse_1d, 'o-', 'Color', [0.0 0.45 0.74], 'LineWidth', 2.0, 'MarkerSize', 8, 'DisplayName', 'Measured RMSE');
hold on;
% 1/sqrt(N) reference line
ref_1d = rmse_1d(2) * (factors_1d / factors_1d(2)).^(-0.5);
loglog(factors_1d, ref_1d, '--k', 'LineWidth', 1.5, 'DisplayName', 'Theory ~ N^{-0.5}');
grid on; box on;
xlabel('Sample Size Factor N / N_0'); ylabel('Autocorrelation RMSE');
title(sprintf('1D: Rate = N^{%.2f}\n(Theory = N^{-0.50})', slope_1d), 'FontWeight', 'bold');
legend('Location', 'northeast');

% --- Subplot 2: 2D ---
subplot(1, 3, 2);
loglog(area_factors_2d, rmse_2d, 's-', 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'MarkerSize', 8, 'DisplayName', 'Measured RMSE');
hold on;
ref_2d = rmse_2d(2) * (area_factors_2d / area_factors_2d(2)).^(-0.5);
loglog(area_factors_2d, ref_2d, '--k', 'LineWidth', 1.5, 'DisplayName', 'Theory ~ Area^{-0.5}');
grid on; box on;
xlabel('Sampled Area Factor A / A_0'); ylabel('Autocorrelation RMSE');
title(sprintf('2D: Rate = Area^{%.2f}\n(Theory = Area^{-0.50})', slope_2d), 'FontWeight', 'bold');
legend('Location', 'northeast');

% --- Subplot 3: 3D ---
subplot(1, 3, 3);
loglog(M_list, rmse_3d, '^-', 'Color', [0.47 0.67 0.19], 'LineWidth', 2.0, 'MarkerSize', 8, 'DisplayName', 'Measured RMSE');
hold on;
ref_3d = rmse_3d(1) * (M_list / M_list(1)).^(-0.5);
loglog(M_list, ref_3d, '--k', 'LineWidth', 1.5, 'DisplayName', 'Theory ~ M^{-0.5}');
grid on; box on;
xlabel('Number of Realizations M (Ensemble)'); ylabel('Autocorrelation RMSE');
title(sprintf('3D: Rate = M^{%.2f}\n(Theory = M^{-0.50})', slope_3d), 'FontWeight', 'bold');
legend('Location', 'northeast');

saveas(fig, 'convergence_clt_verification.png');
fprintf('Figure saved successfully: convergence_clt_verification.png\n');


%% =========================================================================
%% CONCISE AUXILIARY FUNCTIONS
%% =========================================================================

function [r, R_emp] = measure_corr_1d(u, dx, max_lag)
    u = u(:) - mean(u);
    N = length(u);
    var_u = mean(u.^2);
    U = fft(u, 2*N);
    R_full = real(ifft(abs(U).^2));
    lags = (0:max_lag)';
    weights = N - lags;
    R_emp = (R_full(1:max_lag+1) ./ weights) / var_u;
    r = lags * dx;
end

function [r, R_radial] = measure_corr_2d(u, dx, dy, max_lag)
    [Nx, Ny] = size(u);
    u = u - mean(u(:));
    var_u = mean(u(:).^2);
    U = fft2(u, 2*Nx, 2*Ny);
    R_full = real(fftshift(ifft2(abs(U).^2)));
    cx = Nx + 1; cy = Ny + 1;
    lag_x = -max_lag:max_lag; lag_y = -max_lag:max_lag;
    [LX, LY] = meshgrid(lag_x, lag_y);
    W = (Nx - abs(LX)) .* (Ny - abs(LY)); W(W <= 0) = 1;
    R_map = R_full(cy + lag_y, cx + lag_x) ./ (W * var_u);
    [X_d, Y_d] = meshgrid(lag_x * dx, lag_y * dy);
    Dist = sqrt(X_d.^2 + Y_d.^2);
    dr = min(dx, dy);
    r_edges = 0:dr:(max_lag * dr);
    r = [0, (r_edges(1:end-1) + r_edges(2:end)) / 2];
    R_radial = zeros(size(r)); R_radial(1) = 1.0;
    for k = 2:length(r)
        mask = (Dist >= r_edges(k-1)) & (Dist < r_edges(k));
        if any(mask(:)); R_radial(k) = mean(R_map(mask)); else; R_radial(k) = NaN; end
    end
    nan_idx = isnan(R_radial);
    if any(nan_idx); R_radial(nan_idx) = interp1(r(~nan_idx), R_radial(~nan_idx), r(nan_idx), 'linear', 'extrap'); end
end

function [r, R_sph] = measure_corr_3d(u, dx, dy, dz, max_lag)
    [Nx, Ny, Nz] = size(u);
    u = u - mean(u(:));
    var_u = mean(u(:).^2);
    nfft = [2*Nx, 2*Ny, 2*Nz];
    U = fftn(u, nfft);
    R_full = real(fftshift(fftn(abs(U).^2))) / prod(nfft);
    cx = Nx + 1; cy = Ny + 1; cz = Nz + 1;
    lag_x = -max_lag:max_lag;
    [LX, LY, LZ] = meshgrid(lag_x, lag_x, lag_x);
    W = (Nx - abs(LX)) .* (Ny - abs(LY)) .* (Nz - abs(LZ)); W(W <= 0) = 1;
    sub_R = R_full(cy + lag_x, cx + lag_x, cz + lag_x) ./ (W * var_u);
    Dist3D = sqrt((LX * dx).^2 + (LY * dy).^2 + (LZ * dz).^2);
    dr = min([dx, dy, dz]);
    r_edges = 0:dr:(max_lag * dr);
    r = [0, (r_edges(1:end-1) + r_edges(2:end)) / 2];
    R_sph = zeros(size(r)); R_sph(1) = 1.0;
    for k = 2:length(r)
        mask = (Dist3D >= r_edges(k-1)) & (Dist3D < r_edges(k));
        if any(mask(:)); R_sph(k) = mean(sub_R(mask)); else; R_sph(k) = NaN; end
    end
    nan_idx = isnan(R_sph);
    if any(nan_idx); R_sph(nan_idx) = interp1(r(~nan_idx), R_sph(~nan_idx), r(nan_idx), 'linear', 'extrap'); end
end

function u = generate_rf_1d(N, L, model, corrL, extra)
    dk = 2.0 * pi / L;
    kx = (0:N-1)'; kx(kx > N/2) = kx(kx > N/2) - N;
    k_mag = abs(kx) * dk;
    Sk = eval_sk(k_mag, 1, model, corrL, extra);
    u = real(ifft(fft(randn(N, 1)) .* sqrt(Sk)));
    u = (u - mean(u)) / std(u);
end

function u = generate_rf_2d(N, L, model, corrL, extra)
    dk = 2.0 * pi / L;
    kx = (0:N-1); kx(kx > N/2) = kx(kx > N/2) - N;
    [KX, KY] = meshgrid(kx, kx);
    K_mag = sqrt(KX.^2 + KY.^2) * dk;
    Sk = eval_sk(K_mag, 2, model, corrL, extra);
    u = real(ifft2(fft2(randn(N, N)) .* sqrt(Sk)));
    u = (u - mean(u(:))) / std(u(:));
end

function u = generate_rf_3d(N, L, model, corrL, extra)
    dk = 2.0 * pi / L;
    kx = (0:N-1); kx(kx > N/2) = kx(kx > N/2) - N;
    [KX, KY, KZ] = meshgrid(kx, kx, kx);
    K_mag = sqrt(KX.^2 + KY.^2 + KZ.^2) * dk;
    Sk = eval_sk(K_mag, 3, model, corrL, extra);
    u = real(ifftn(fftn(randn(N, N, N)) .* sqrt(Sk)));
    u = (u - mean(u(:))) / std(u(:));
end

function Sk = eval_sk(k_mag, nDim, model, corrL, extra)
    switch lower(model)
        case 'exponential'
            switch nDim
                case 1; a = 2.0; Sk = (corrL / (2*pi)) ./ (1.0 + (k_mag*corrL/a).^2);
                case 2; a = sqrt(2); Sk = (corrL^2 / (4*pi)) ./ ((1.0 + (k_mag*corrL/a).^2).^1.5);
                case 3; a = 6^(1/3); Sk = (corrL^3 / (6*pi^2)) ./ ((1.0 + (k_mag*corrL/a).^2).^2.0);
            end
        case 'karman'
            nu = 0.5; if isfield(extra, 'nu'); nu = extra.nu; end
            switch nDim
                case 1
                    a = 2.0 * sqrt(pi) * exp(gammaln(nu + 0.5) - gammaln(nu));
                    Sk = (corrL / (2*pi)) ./ ((1.0 + (k_mag*corrL/a).^2).^(nu + 0.5));
                case 2
                    a = 2.0 * sqrt(nu);
                    Sk = (corrL^2 / (4*pi)) ./ ((1.0 + (k_mag*corrL/a).^2).^(nu + 1.0));
                case 3
                    a = (6.0 * sqrt(pi) * exp(gammaln(nu + 1.5) - gammaln(nu)))^(1/3);
                    Sk = (corrL^3 / (6*pi^2)) ./ ((1.0 + (k_mag*corrL/a).^2).^(nu + 1.5));
            end
        case 'monodisperse'
            D = extra.D; eta = extra.eta;
            % 3D Percus-Yevick
            R_sph = D / 2.0;
            eta_w = eta * (1.0 - eta / 16.0);
            D_w = D * (eta_w / eta)^(1/3);
            rho_w = 6.0 * eta_w / (pi * D_w^3);
            l1 = -(1.0 + 2.0*eta_w)^2 / (1.0 - eta_w)^4;
            l2 = 6.0 * eta_w * (1.0 + 0.5*eta_w)^2 / (1.0 - eta_w)^4;
            c_val = -0.5 * eta_w * (1.0 + 2.0*eta_w)^2 / (1.0 - eta_w)^4;
            Sk = zeros(size(k_mag));
            for i = 1:numel(k_mag)
                k = k_mag(i);
                u = k * D_w;
                if u < 1e-5
                    c_k = (4.0*pi*D_w^3) * (l1/3.0 + l2/4.0 + c_val/6.0);
                    m_k = (4.0/3.0) * pi * R_sph^3;
                else
                    c_k = (4.0*pi / (k^3)) * ( ...
                        l1 * (sin(u) - u*cos(u)) + ...
                        l2 * (2.0*u*sin(u) - (u^2 - 2.0)*cos(u) - 2.0)/u + ...
                        c_val * ((4.0*u^3 - 24.0*u)*sin(u) - (u^4 - 12.0*u^2 + 24.0)*cos(u) + 24.0)/(u^3) );
                    m_k = (4.0 * pi / (k^3)) * (sin(k*R_sph) - (k*R_sph)*cos(k*R_sph));
                end
                S_py = 1.0 / max(1e-8, 1.0 - rho_w * c_k);
                g_cs = (1.0 - 0.5*eta) / (1.0 - eta)^3;
                g_py = (1.0 + 0.5*eta_w) / (1.0 - eta_w)^2;
                A_vw = (g_cs - g_py) * D_w;
                mu = 24.0 * A_vw / D_w;
                denom = max(1e-8, mu^2 + k^2);
                I_cos = (mu * cos(k * D_w) - k * sin(k * D_w)) / denom;
                I_sin = (mu * sin(k * D_w) + k * cos(k * D_w)) / denom;
                dh = (4.0 * pi * A_vw / k) * (D_w * I_sin + (2.0*mu*I_cos + (mu^2 - k^2)*I_sin/k)/denom);
                S_tot = max(1e-8, S_py + rho_w * dh);
                rho_p = 6.0 * eta / (pi * D^3);
                Sk(i) = rho_p * (m_k^2) * S_tot / ((2*pi)^3 * eta * (1.0 - eta));
            end
        otherwise
            Sk = (corrL^nDim) * exp(-(k_mag .* corrL).^2 / 8.0);
    end
end

function R_theo = get_theo_corr(r, model, corrL, nDim, extra)
    r = r(:);
    switch lower(model)
        case 'exponential'
            switch nDim
                case 1; a = 2.0;
                case 2; a = sqrt(2.0);
                case 3; a = 6.0^(1/3);
            end
            R_theo = exp(- a * r / corrL);
        case 'karman'
            nu = extra.nu;
            switch nDim
                case 1; a = 2.0 * sqrt(pi) * exp(gammaln(nu + 0.5) - gammaln(nu));
                case 2; a = 2.0 * sqrt(nu);
                case 3; a = (6.0 * sqrt(pi) * exp(gammaln(nu + 1.5) - gammaln(nu)))^(1/3);
            end
            z = (a * r) / corrL;
            R_theo = zeros(size(r));
            R_theo(z < 1e-6) = 1.0;
            pos = z >= 1e-6;
            pref = (2.0^(1.0 - nu)) / gamma(nu);
            R_theo(pos) = pref * (z(pos).^nu) .* besselk(nu, z(pos));
        case 'monodisperse'
            D = extra.D; eta = extra.eta;
            k_pts = linspace(1e-4, 25.0 * pi / D, 1000);
            dk = k_pts(2) - k_pts(1);
            Sk_pts = eval_sk(k_pts, 3, 'monodisperse', corrL, extra);
            R_theo = zeros(size(r));
            for ir = 1:length(r)
                rr = r(ir);
                if rr < 1e-6; sinc_t = ones(size(k_pts)); else; sinc_t = sin(k_pts * rr) ./ (k_pts * rr); end
                R_theo(ir) = sum(Sk_pts .* sinc_t .* (k_pts.^2) * (4.0*pi)) * dk;
            end
            if R_theo(1) > 0; R_theo = R_theo / R_theo(1); end
    end
end
