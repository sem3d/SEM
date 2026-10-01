%% =========================================================================
% RUN_COMPREHENSIVE_BENCHMARK_ALL_PDFS_AND_CORRELATIONS.M
%
% Comprehensive Statistical Validation Benchmark:
%   - All 5 Marginal PDFs (Gaussian, Lognormal, Gamma, Weibull, Uniform)
%   - All 8 Spatial Correlation Models (1D, 2D and 3D)
%   - Rigorous Verification of the 1/sqrt(N) Convergence (Law of Large Numbers / CLT)
%     when varying by factors of 4 (N/4, N, 4N, 16N)
%   - Automatic saving of all figures, MAT data and metrics.
% =========================================================================

clear; close all; clc;

fprintf('=========================================================================\n');
fprintf('  FULL BENCHMARK: ALL PDFS AND ALL CORRELATIONS (1D, 2D, 3D)             \n');
fprintf('  Validation of 1/sqrt(N) Convergence (Law of Large Numbers and CLT)      \n');
fprintf('=========================================================================\n\n');

rng(2026); % Strict reproducibility

benchmark_dir = fileparts(mfilename('fullpath'));
if isempty(benchmark_dir); benchmark_dir = pwd; end

%% =========================================================================
%% PART 1: VALIDATION OF ALL 5 MARGINAL PDFs
%% =========================================================================
fprintf('=========================================================================\n');
fprintf('  PART 1: ANALYSIS OF THE 5 MARGINAL PDFs (ROSENBLATT TRANSFORMS)        \n');
fprintf('=========================================================================\n\n');

pdf_defs = { ...
    'Gaussian',  10.0, 4.0; ...
    'Lognormal',  5.0, 2.0; ...
    'Gamma',      8.0, 16.0; ...
    'Weibull',    6.0, 4.0; ...
    'Uniform',   10.0, 3.0  ...
};

N_pdf_list = [1024, 4096, 16384, 65536]; % Factors of 4!
factors_pdf = N_pdf_list / N_pdf_list(2); % [0.25, 1.0, 4.0, 16.0]
n_pdfs = size(pdf_defs, 1);

pdf_results = struct();

for p = 1:n_pdfs
    name = pdf_defs{p, 1};
    mu_target = pdf_defs{p, 2};
    var_target = pdf_defs{p, 3};
    
    fprintf('--> Testing PDF [%d/5]: %s (Theoretical Mean = %.2f, Theoretical Var = %.2f)\n', ...
        p, name, mu_target, var_target);
    
    pdf_results(p).name = name;
    pdf_results(p).mu_target = mu_target;
    pdf_results(p).var_target = var_target;
    pdf_results(p).N_list = N_pdf_list;
    pdf_results(p).ks_errors = zeros(size(N_pdf_list));
    pdf_results(p).mean_errors = zeros(size(N_pdf_list));
    pdf_results(p).var_errors = zeros(size(N_pdf_list));
    
    for i = 1:length(N_pdf_list)
        N = N_pdf_list(i);
        Z = randn(N, 1);
        U = 0.5 * erfc(-Z / sqrt(2.0)); % Uniform [0, 1]
        U(U <= 0) = 1e-12; U(U >= 1) = 1.0 - 1e-12;
        
        switch name
            case 'Gaussian'
                X = mu_target + sqrt(var_target) * Z;
                F_theo = @(x) 0.5 * erfc(-(x - mu_target) / (sqrt(2*var_target)));
                f_theo = @(x) (1.0 / (sqrt(2*pi*var_target))) * exp(-0.5 * ((x - mu_target).^2) / var_target);
            case 'Lognormal'
                sig_y = sqrt(log(1.0 + var_target / (mu_target^2)));
                mu_y = log(mu_target) - 0.5 * (sig_y^2);
                X = exp(mu_y + sig_y * Z);
                F_theo = @(x) 0.5 * erfc(-(log(max(1e-12, x)) - mu_y) / (sqrt(2)*sig_y));
                f_theo = @(x) (1.0 ./ (max(1e-12, x) * sig_y * sqrt(2*pi))) .* exp(-0.5 * ((log(max(1e-12, x)) - mu_y)/sig_y).^2);
            case 'Gamma'
                alpha = (mu_target^2) / var_target;
                beta  = var_target / mu_target;
                X = gaminv(U, alpha, beta);
                F_theo = @(x) gamcdf(x, alpha, beta);
                f_theo = @(x) gampdf(x, alpha, beta);
            case 'Weibull'
                cv = sqrt(var_target) / mu_target;
                k_w = (cv)^(-1.086);
                for it = 1:20
                    f_k = sqrt(gamma(1 + 2/k_w) / (gamma(1 + 1/k_w)^2) - 1.0) - cv;
                    if abs(f_k) < 1e-7; break; end
                    k_w = k_w * (1.0 - f_k * 0.2);
                end
                lambda_w = mu_target / gamma(1 + 1/k_w);
                X = lambda_w * (-log(1.0 - U)).^(1.0 / k_w);
                F_theo = @(x) 1.0 - exp(-(max(0, x)/lambda_w).^k_w);
                f_theo = @(x) (k_w / lambda_w) * ((max(1e-12, x)/lambda_w).^(k_w - 1)) .* exp(-(max(0, x)/lambda_w).^k_w);
            case 'Uniform'
                delta = sqrt(3.0 * var_target);
                a_unif = mu_target - delta;
                b_unif = mu_target + delta;
                X = a_unif + (b_unif - a_unif) * U;
                F_theo = @(x) min(1.0, max(0.0, (x - a_unif) / (b_unif - a_unif)));
                f_theo = @(x) (x >= a_unif & x <= b_unif) / (b_unif - a_unif);
        end
        
        pdf_results(p).mean_errors(i) = abs(mean(X) - mu_target);
        pdf_results(p).var_errors(i)  = abs(var(X) - var_target);
        
        % Kolmogorov-Smirnov test
        X_sort = sort(X);
        emp_cdf = (1:N)' / N;
        theo_cdf = F_theo(X_sort);
        D_ks = max(abs(emp_cdf - theo_cdf));
        pdf_results(p).ks_errors(i) = D_ks;
        
        if i == length(N_pdf_list)
            pdf_results(p).sample_data = X;
            pdf_results(p).f_theo = f_theo;
        end
        
        ratio_str = '';
        if i > 1
            ratio = pdf_results(p).ks_errors(i-1) / pdf_results(p).ks_errors(i);
            ratio_str = sprintf(' | KS decay: %.2fx (Theory: 2.00x)', ratio);
        end
        fprintf('    N = %5d: Mean Error = %.4e | Var Error = %.4e | KS Distance = %.4e%s\n', ...
            N, pdf_results(p).mean_errors(i), pdf_results(p).var_errors(i), D_ks, ratio_str);
    end
    
    p_fit = polyfit(log(factors_pdf), log(pdf_results(p).ks_errors), 1);
    pdf_results(p).slope = p_fit(1);
    fprintf('    >> KS Convergence Rate (log-log): N^(%.3f) [Theory: N^(-0.500)]\n\n', pdf_results(p).slope);
end

%% Figure of the 5 PDFs
fig_pdf = figure('Name', 'Complete Validation of the 5 Marginal PDFs', 'Color', 'w', ...
                 'Units', 'normalized', 'Position', [0.05, 0.08, 0.90, 0.82]);

colors_pdf = [
    0.15 0.45 0.75;
    0.85 0.35 0.20;
    0.20 0.65 0.35;
    0.55 0.30 0.75;
    0.90 0.65 0.10
];

for p = 1:5
    subplot(2, 3, p);
    X = pdf_results(p).sample_data;
    f_theo = pdf_results(p).f_theo;
    
    histogram(X, 60, 'Normalization', 'pdf', 'FaceColor', colors_pdf(p, :), ...
        'FaceAlpha', 0.6, 'EdgeColor', 'w', 'DisplayName', 'Sample');
    hold on;
    x_grid = linspace(min(X), max(X), 300);
    plot(x_grid, f_theo(x_grid), 'k-', 'LineWidth', 2.2, 'DisplayName', 'Theory');
    grid on; box on;
    xlabel('x', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Density p(x)', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('%s (N = 65536)\n\\mu=%.2f, \\sigma^2=%.2f', ...
        pdf_results(p).name, pdf_results(p).mu_target, pdf_results(p).var_target), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'northeast');
end

% Subplot 6: KS convergence vs N
subplot(2, 3, 6);
for p = 1:5
    loglog(N_pdf_list, pdf_results(p).ks_errors, 'o-', 'Color', colors_pdf(p, :), ...
        'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', sprintf('%s (s=%.2f)', pdf_results(p).name, pdf_results(p).slope));
    hold on;
end
ref_ks = 0.8 * (N_pdf_list / N_pdf_list(1)).^(-0.5);
loglog(N_pdf_list, ref_ks, '--k', 'LineWidth', 2.0, 'DisplayName', 'Theory ~ 1/\surd N');
grid on; box on;
xlabel('Number of Samples N', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('KS Distance (sup |F_N - F|)', 'FontSize', 10, 'FontWeight', 'bold');
title('PDF Convergence: Rate vs N', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 8);

saveas(fig_pdf, fullfile(benchmark_dir, 'benchmark_all_pdfs_histograms_and_convergence.png'));
fprintf('PDF figure saved: benchmark_all_pdfs_histograms_and_convergence.png\n\n');


%% =========================================================================
%% PART 2: VALIDATION OF ALL 8 SPATIAL CORRELATION MODELS
%% =========================================================================
fprintf('=========================================================================\n');
fprintf('  PART 2: ANALYSIS OF THE 8 SPATIAL CORRELATION MODELS (1D, 2D, 3D)     \n');
fprintf('=========================================================================\n\n');

corr_models = { ...
    'gaussian',     'Gaussian',        1.5; ...
    'cos',          'Cosine',          2.0; ...
    'karman',       'von Kármán',       1.5; ...
    'exponential',  'Exponential',      1.5; ...
    'powerlaw',     'Power Law',  1.5; ...
    'lowpass',      'Low-Pass (Sinc)', 2.0; ...
    'triangular',   'Triangular',       2.0; ...
    'monodisperse', 'Monodisperse Spheres', 1.2 ...
};

n_corrs = size(corr_models, 1);
corr_results_1d = struct();
corr_results_2d = struct();
corr_results_3d = struct();

%% -------------------------------------------------------------------------
%% 2.1 1D CONVERGENCE FOR ALL 8 MODELS
%% -------------------------------------------------------------------------
fprintf('--- 2.1 Testing 1D Models (N = 512, 2048, 8192, 32768) ---\n');
dx1 = 0.02;
N_list_1d = [512, 2048, 8192, 32768]; % Factors of 4!
factors_1d = N_list_1d / N_list_1d(2);
max_r_1d = 2.5; % Fixed evaluation window: fits comfortably in all grids
max_lag_1d = round(max_r_1d / dx1); % 125 points (well below 512/2 = 256)
r_eval_1d = (0:max_lag_1d)' * dx1;

fig_corr1d = figure('Name', 'All 8 1D Correlation Models', 'Color', 'w', ...
                    'Units', 'normalized', 'Position', [0.04, 0.05, 0.92, 0.88]);

for c = 1:n_corrs
    mod_id   = corr_models{c, 1};
    mod_name = corr_models{c, 2};
    corrL    = corr_models{c, 3};
    
    fprintf('  1D Model [%d/8]: %s (\\ell = %.2f)\n', c, mod_name, corrL);
    corr_results_1d(c).id = mod_id;
    corr_results_1d(c).name = mod_name;
    corr_results_1d(c).rmse = zeros(size(N_list_1d));
    
    extra = struct('nu', 0.7, 'eta', 0.25, 'D', 1.2, 'use_vw', true);
    R_theo = get_theoretical_correlation(r_eval_1d, mod_id, corrL, 1, extra);
    R_theo = R_theo(:);
    
    for i = 1:length(N_list_1d)
        N = N_list_1d(i);
        L = N * dx1;
        dk = 2.0 * pi / L;
        kx = (0:N-1)'; kx(kx > N/2) = kx(kx > N/2) - N;
        kmag = abs(kx) * dk;
        
        Sk = eval_psd_1d(kmag, mod_id, corrL, extra);
        Sk = max(1e-14, Sk);
        
        err_rep = zeros(6, 1);
        for rep = 1:6
            noise = fft(randn(N, 1));
            u = real(ifft(noise .* sqrt(Sk)));
            u_std = std(u);
            if u_std < 1e-12; u_std = 1.0; end
            u = (u - mean(u)) / u_std;
            
            U = fft(u, 2*N);
            R_raw = real(ifft(abs(U).^2));
            weights = N - (0:max_lag_1d)';
            R_emp = (R_raw(1:max_lag_1d+1) ./ weights) / mean(u.^2);
            err_rep(rep) = sqrt(mean((R_emp - R_theo).^2));
        end
        corr_results_1d(c).rmse(i) = mean(err_rep);
        
        if i == 2
            corr_results_1d(c).r_eval = r_eval_1d;
            corr_results_1d(c).R_emp = R_emp;
            corr_results_1d(c).R_theo = R_theo;
        end
    end
    
    p_fit = polyfit(log(factors_1d), log(corr_results_1d(c).rmse), 1);
    corr_results_1d(c).slope = p_fit(1);
    corr_results_1d(c).r2 = 1 - sum((corr_results_1d(c).R_emp - R_theo).^2) / ...
                                max(1e-12, sum((corr_results_1d(c).R_emp - mean(corr_results_1d(c).R_emp)).^2));
    
    fprintf('    RMSE (N=512 -> 32768): [%.3e -> %.3e] | Rate N^(%.3f) | R^2 = %.4f\n', ...
        corr_results_1d(c).rmse(1), corr_results_1d(c).rmse(end), corr_results_1d(c).slope, corr_results_1d(c).r2);
    
    subplot(2, 4, c);
    plot(r_eval_1d, corr_results_1d(c).R_emp, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Measured');
    hold on;
    plot(r_eval_1d, R_theo, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Theory');
    grid on; box on;
    xlabel('Distance r', 'FontSize', 9, 'FontWeight', 'bold');
    ylabel('R(r)', 'FontSize', 9, 'FontWeight', 'bold');
    title(sprintf('%s (1D)\nR^2=%.3f, Rate: N^{%.2f}', mod_name, corr_results_1d(c).r2, corr_results_1d(c).slope), ...
        'FontSize', 10, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 8);
end

saveas(fig_corr1d, fullfile(benchmark_dir, 'benchmark_all_correlations_1d.png'));
fprintf('1D figure saved: benchmark_all_correlations_1d.png\n\n');


%% -------------------------------------------------------------------------
%% 2.2 2D CONVERGENCE FOR ALL 8 MODELS
%% -------------------------------------------------------------------------
fprintf('--- 2.2 Testing 2D Models (Grid 64x64, 128x128, 256x256) ---\n');
dx2 = 0.08;
N_list_2d = [64, 128, 256]; % Areas with factors of 4!
area_factors = (N_list_2d / 128).^2;

% Fixed window for 2D: max_r = 2.0 (max_lag = 25 points, fits comfortably in N=64 where N/2 = 32)
max_r_2d = 2.0;
max_lag_2d = round(max_r_2d / dx2);

fig_corr2d = figure('Name', 'All 8 2D Correlation Models (Radial)', 'Color', 'w', ...
                    'Units', 'normalized', 'Position', [0.04, 0.05, 0.92, 0.88]);

for c = 1:n_corrs
    mod_id   = corr_models{c, 1};
    mod_name = corr_models{c, 2};
    corrL    = corr_models{c, 3};
    
    fprintf('  2D Model [%d/8]: %s (\\ell = %.2f)\n', c, mod_name, corrL);
    corr_results_2d(c).id = mod_id;
    corr_results_2d(c).name = mod_name;
    corr_results_2d(c).rmse = zeros(size(N_list_2d));
    
    extra = struct('nu', 0.7, 'eta', 0.25, 'D', 1.2, 'use_vw', true);
    
    for i = 1:length(N_list_2d)
        N = N_list_2d(i);
        L = N * dx2;
        dk = 2.0 * pi / L;
        kx = (0:N-1); kx(kx > N/2) = kx(kx > N/2) - N;
        [KX, KY] = meshgrid(kx, kx);
        K_mag = sqrt(KX.^2 + KY.^2) * dk;
        
        Sk = eval_psd_2d(K_mag, mod_id, corrL, extra);
        Sk = max(1e-14, Sk);
        
        u = real(ifft2(fft2(randn(N, N)) .* sqrt(Sk)));
        u_std = std(u(:));
        if u_std < 1e-12; u_std = 1.0; end
        u = (u - mean(u(:))) / u_std;
        
        U = fft2(u, 2*N, 2*N);
        R_full = real(fftshift(ifft2(abs(U).^2)));
        cx = N + 1; cy = N + 1;
        
        lag_idx = -max_lag_2d:max_lag_2d;
        [LX, LY] = meshgrid(lag_idx, lag_idx);
        W = (N - abs(LX)) .* (N - abs(LY)); W(W <= 0) = 1;
        R_map = R_full(cy + lag_idx, cx + lag_idx) ./ (W * mean(u(:).^2));
        
        Dist = sqrt((LX * dx2).^2 + (LY * dx2).^2);
        r_edges = 0:dx2:max_r_2d;
        r_cents = [0, (r_edges(1:end-1) + r_edges(2:end)) / 2];
        R_rad = zeros(size(r_cents)); R_rad(1) = 1.0;
        for k = 2:length(r_cents)
            mask = (Dist >= r_edges(k-1)) & (Dist < r_edges(k));
            if any(mask(:)); R_rad(k) = mean(R_map(mask)); else; R_rad(k) = NaN; end
        end
        nan_m = isnan(R_rad);
        if any(nan_m); R_rad(nan_m) = interp1(r_cents(~nan_m), R_rad(~nan_m), r_cents(nan_m), 'linear', 'extrap'); end
        
        R_theo = get_theoretical_correlation(r_cents, mod_id, corrL, 2, extra);
        R_theo = R_theo(:);
        
        corr_results_2d(c).rmse(i) = sqrt(mean((R_rad(:) - R_theo(:)).^2));
        
        if i == 2
            corr_results_2d(c).r_cents = r_cents;
            corr_results_2d(c).R_rad = R_rad;
            corr_results_2d(c).R_theo = R_theo;
        end
    end
    
    p_fit = polyfit(log(area_factors), log(corr_results_2d(c).rmse), 1);
    corr_results_2d(c).slope = p_fit(1);
    corr_results_2d(c).r2 = 1 - sum((corr_results_2d(c).R_rad(:) - corr_results_2d(c).R_theo(:)).^2) / ...
                                max(1e-12, sum((corr_results_2d(c).R_rad(:) - mean(corr_results_2d(c).R_rad(:))).^2));
    
    fprintf('    RMSE (Area x0.25 -> x4.0): [%.3e -> %.3e] | Rate Area^(%.3f) | R^2 = %.4f\n', ...
        corr_results_2d(c).rmse(1), corr_results_2d(c).rmse(end), corr_results_2d(c).slope, corr_results_2d(c).r2);
    
    subplot(2, 4, c);
    plot(corr_results_2d(c).r_cents, corr_results_2d(c).R_rad, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'DisplayName', 'Radial 2D');
    hold on;
    plot(corr_results_2d(c).r_cents, corr_results_2d(c).R_theo, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Theory');
    grid on; box on;
    xlabel('Distance r', 'FontSize', 9, 'FontWeight', 'bold');
    ylabel('R(r)', 'FontSize', 9, 'FontWeight', 'bold');
    title(sprintf('%s (2D)\nR^2=%.3f, Rate: Area^{%.2f}', mod_name, corr_results_2d(c).r2, corr_results_2d(c).slope), ...
        'FontSize', 10, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 8);
end

saveas(fig_corr2d, fullfile(benchmark_dir, 'benchmark_all_correlations_2d.png'));
fprintf('2D figure saved: benchmark_all_correlations_2d.png\n\n');


%% -------------------------------------------------------------------------
%% 2.3 3D VALIDATION AND ENSEMBLE FOR ALL 8 MODELS
%% -------------------------------------------------------------------------
fprintf('--- 2.3 Testing 3D Models (Grid 32x32x32 and Ensemble) ---\n');
N3 = 32;
L3 = 8.0;
dx3 = L3 / N3;
max_lag_3d = 8;
dk3 = 2.0 * pi / L3;

kx3 = (0:N3-1); kx3(kx3 > N3/2) = kx3(kx3 > N3/2) - N3;
[KX3, KY3, KZ3] = meshgrid(kx3, kx3, kx3);
K_mag3 = sqrt(KX3.^2 + KY3.^2 + KZ3.^2) * dk3;

lag_idx3 = -max_lag_3d:max_lag_3d;
[LX3, LY3, LZ3] = meshgrid(lag_idx3, lag_idx3, lag_idx3);
W3 = (N3 - abs(LX3)) .* (N3 - abs(LY3)) .* (N3 - abs(LZ3)); W3(W3 <= 0) = 1;
Dist3 = sqrt((LX3 * dx3).^2 + (LY3 * dx3).^2 + (LZ3 * dx3).^2);
r_cents3 = [0, (0.5:max_lag_3d) * dx3];

fig_corr3d = figure('Name', 'All 8 3D Correlation Models (Spherical)', 'Color', 'w', ...
                    'Units', 'normalized', 'Position', [0.04, 0.05, 0.92, 0.88]);

for c = 1:n_corrs
    mod_id   = corr_models{c, 1};
    mod_name = corr_models{c, 2};
    corrL    = corr_models{c, 3};
    
    fprintf('  3D Model [%d/8]: %s (\\ell = %.2f)\n', c, mod_name, corrL);
    corr_results_3d(c).id = mod_id;
    corr_results_3d(c).name = mod_name;
    
    extra = struct('nu', 0.7, 'eta', 0.25, 'D', 1.2, 'use_vw', true);
    Sk3 = eval_psd_3d(K_mag3, mod_id, corrL, extra);
    Sk3 = max(1e-14, Sk3);
    
    M_ens = 4;
    R_accum = zeros(size(r_cents3));
    for m = 1:M_ens
        u = real(ifftn(fftn(randn(N3, N3, N3)) .* sqrt(Sk3)));
        u_std = std(u(:));
        if u_std < 1e-12; u_std = 1.0; end
        u = (u - mean(u(:))) / u_std;
        
        U = fftn(u, [2*N3, 2*N3, 2*N3]);
        R_full = real(fftshift(ifftn(abs(U).^2)));
        cx = N3 + 1;
        sub_R = R_full(cx + lag_idx3, cx + lag_idx3, cx + lag_idx3) ./ (W3 * mean(u(:).^2));
        
        R_sph = zeros(size(r_cents3)); R_sph(1) = 1.0;
        for k = 2:length(r_cents3)
            mask = (Dist3 >= (k-1.5)*dx3) & (Dist3 < (k-0.5)*dx3);
            if any(mask(:)); R_sph(k) = mean(sub_R(mask)); end
        end
        R_accum = R_accum + R_sph;
    end
    R_sph_mean = R_accum / M_ens;
    
    R_theo3 = get_theoretical_correlation(r_cents3, mod_id, corrL, 3, extra);
    R_theo3 = R_theo3(:);
    
    corr_results_3d(c).r_cents = r_cents3;
    corr_results_3d(c).R_sph = R_sph_mean(:);
    corr_results_3d(c).R_theo = R_theo3;
    corr_results_3d(c).rmse = sqrt(mean((corr_results_3d(c).R_sph - R_theo3).^2));
    corr_results_3d(c).r2 = 1 - sum((corr_results_3d(c).R_sph - R_theo3).^2) / ...
                                max(1e-12, sum((corr_results_3d(c).R_sph - mean(corr_results_3d(c).R_sph)).^2));
    
    fprintf('    RMSE 3D (M=4) = %.4e | R^2 = %.4f\n', corr_results_3d(c).rmse, corr_results_3d(c).r2);
    
    subplot(2, 4, c);
    plot(r_cents3, corr_results_3d(c).R_sph, 'Color', [0.2 0.6 0.3], 'LineWidth', 2.0, 'DisplayName', 'Spherical 3D');
    hold on;
    plot(r_cents3, R_theo3, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Theory');
    grid on; box on;
    xlabel('Distance r', 'FontSize', 9, 'FontWeight', 'bold');
    ylabel('R(r)', 'FontSize', 9, 'FontWeight', 'bold');
    title(sprintf('%s (3D)\nRMSE=%.3e, R^2=%.3f', mod_name, corr_results_3d(c).rmse, corr_results_3d(c).r2), ...
        'FontSize', 10, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 8);
end

saveas(fig_corr3d, fullfile(benchmark_dir, 'benchmark_all_correlations_3d.png'));
fprintf('3D figure saved: benchmark_all_correlations_3d.png\n\n');


%% -------------------------------------------------------------------------
%% 2.4 SUMMARY FIGURE OF THE CONVERGENCE RATES (CLT 1/sqrt(N))
%% -------------------------------------------------------------------------
fig_rates = figure('Name', 'Summary of the 1/sqrt(N) Convergence Rates', 'Color', 'w', ...
                   'Units', 'normalized', 'Position', [0.08, 0.15, 0.84, 0.48]);

% --- Subplot 1: 1D rates ---
subplot(1, 2, 1);
for c = 1:n_corrs
    loglog(factors_1d, corr_results_1d(c).rmse, 'o-', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%s (s=%.2f)', corr_models{c, 2}, corr_results_1d(c).slope));
    hold on;
end
ref_c1 = 0.20 * (factors_1d / factors_1d(2)).^(-0.5);
loglog(factors_1d, ref_c1, '--k', 'LineWidth', 2.2, 'DisplayName', 'Theory ~ 1/\surd N');
grid on; box on;
xlabel('1D Domain Size Factor (N / N_0)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Autocorrelation RMSE', 'FontSize', 11, 'FontWeight', 'bold');
title('1D: Convergence of All 8 Models', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 8);

% --- Subplot 2: 2D rates ---
subplot(1, 2, 2);
for c = 1:n_corrs
    loglog(area_factors, corr_results_2d(c).rmse, 's-', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%s (s=%.2f)', corr_models{c, 2}, corr_results_2d(c).slope));
    hold on;
end
ref_c2 = 0.10 * (area_factors / area_factors(2)).^(-0.5);
loglog(area_factors, ref_c2, '--k', 'LineWidth', 2.2, 'DisplayName', 'Theory ~ 1/\surd Area');
grid on; box on;
xlabel('2D Sampled Area Factor (Area / Area_0)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('RMSE Radial', 'FontSize', 11, 'FontWeight', 'bold');
title('2D: Convergence of All 8 Models', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 8);

saveas(fig_rates, fullfile(benchmark_dir, 'benchmark_all_correlations_convergence_rates.png'));
fprintf('Rates figure saved: benchmark_all_correlations_convergence_rates.png\n\n');


%% =========================================================================
%% PART 3: SAVE ALL RESULTS TO MAT
%% =========================================================================
mat_out = fullfile(benchmark_dir, 'benchmark_results_all.mat');
save(mat_out, 'pdf_results', 'corr_results_1d', 'corr_results_2d', 'corr_results_3d', ...
     'N_pdf_list', 'factors_pdf', 'N_list_1d', 'factors_1d', 'N_list_2d', 'area_factors');
fprintf('Data saved successfully to MAT: %s\n\n', mat_out);

fprintf('=========================================================================\n');
fprintf('   FULL BENCHMARK COMPLETED SUCCESSFULLY!                                \n');
fprintf('=========================================================================\n');


%% =========================================================================
%% AUXILIARY PSD EVALUATION FUNCTIONS
%% =========================================================================

function Sk = eval_psd_1d(k, model, corrL, extra)
    switch lower(model)
        case 'gaussian'
            Sk = (corrL / sqrt(2*pi)) * exp(-0.5 * (k * corrL).^2);
        case 'cos'
            k0 = 2.0 * pi / corrL;
            w = 0.15 * k0;
            Sk = exp(-0.5 * ((k - k0)/w).^2) + exp(-0.5 * ((k + k0)/w).^2);
        case 'karman'
            nu = 0.7; if isfield(extra, 'nu'); nu = extra.nu; end
            a = 2.0 * sqrt(pi) * exp(gammaln(nu + 0.5) - gammaln(nu));
            pref = corrL * exp(gammaln(nu + 0.5) - gammaln(nu)) / (sqrt(pi) * a);
            Sk = pref ./ ((1.0 + (k * corrL / a).^2).^(nu + 0.5));
        case 'exponential'
            Sk = (corrL / pi) ./ (1.0 + (k * corrL).^2);
        case 'powerlaw'
            a = 2.0;
            Sk = (corrL / (2*a)) ./ ((1.0 + (k * corrL / a).^2).^1.5);
        case 'lowpass'
            kmax = pi / corrL;
            Sk = double(k <= kmax);
        case 'triangular'
            Sk = (corrL / (2*pi)) * (sinc(k * corrL / (2*pi))).^2;
        case 'monodisperse'
            D = 1.2; eta = 0.25; if isfield(extra, 'D'); D = extra.D; end; if isfield(extra, 'eta'); eta = extra.eta; end
            u = k * D / 2.0;
            m_k = 2.0 * sin(u) ./ max(1e-8, k);
            S_1d = (1.0 - eta^2) ./ max(1e-8, 1.0 - 2.0*eta*cos(k*D) + eta^2);
            Sk = (1.0 / D) * (m_k.^2) .* S_1d / (2.0*pi * (1.0 - eta));
        otherwise
            Sk = (corrL / pi) ./ (1.0 + (k * corrL).^2);
    end
end

function Sk = eval_psd_2d(K_mag, model, corrL, extra)
    switch lower(model)
        case 'gaussian'
            Sk = (corrL^2 / (4*pi)) * exp(- (K_mag * corrL / 2.0).^2);
        case 'cos'
            k0 = 2.0 * pi / corrL;
            w = 0.15 * k0;
            Sk = exp(-0.5 * ((K_mag - k0)/w).^2);
        case 'karman'
            nu = 0.7; if isfield(extra, 'nu'); nu = extra.nu; end
            a = 2.0 * sqrt(nu);
            Sk = (corrL^2 / (4*pi)) ./ ((1.0 + (K_mag * corrL / a).^2).^(nu + 1.0));
        case 'exponential'
            Sk = (corrL^2 / (2*pi)) ./ ((1.0 + (K_mag * corrL).^2).^1.5);
        case 'powerlaw'
            a = sqrt(2.0);
            Sk = (corrL^2 / (2*pi*a^2)) ./ ((1.0 + (K_mag * corrL / a).^2).^2.0);
        case 'lowpass'
            kmax = 2.0 / corrL;
            Sk = double(K_mag <= kmax);
        case 'triangular'
            kmax = 3.0 / corrL;
            Sk = max(0.0, 1.0 - K_mag / kmax);
        case 'monodisperse'
            D = 1.2; eta = 0.25; if isfield(extra, 'D'); D = extra.D; end; if isfield(extra, 'eta'); eta = extra.eta; end
            u = K_mag * D / 2.0;
            J_term = (2.0 * besselj(1, max(1e-8, u)) ./ max(1e-8, u)).^2;
            Sk = J_term ./ (1.0 + 4.0*eta * (besselj(1, max(1e-8, K_mag*D)) ./ max(1e-8, K_mag*D)));
            Sk = max(1e-12, Sk);
        otherwise
            Sk = (corrL^2 / (2*pi)) ./ ((1.0 + (K_mag * corrL).^2).^1.5);
    end
end

function Sk = eval_psd_3d(K_mag, model, corrL, extra)
    switch lower(model)
        case 'gaussian'
            Sk = (corrL^3 / (8 * pi^1.5)) * exp(- (K_mag * corrL / 2.0).^2);
        case 'cos'
            k0 = 2.0 * pi / corrL;
            w = 0.15 * k0;
            Sk = exp(-0.5 * ((K_mag - k0)/w).^2);
        case 'karman'
            nu = 0.7; if isfield(extra, 'nu'); nu = extra.nu; end
            a = (6.0 * sqrt(pi) * exp(gammaln(nu + 1.5) - gammaln(nu)))^(1/3);
            Sk = (corrL^3 / (6*pi^2)) ./ ((1.0 + (K_mag * corrL / a).^2).^(nu + 1.5));
        case 'exponential'
            Sk = (corrL^3 / (pi^2)) ./ ((1.0 + (K_mag * corrL).^2).^2.0);
        case 'powerlaw'
            a = 6.0^(1/3);
            Sk = (corrL^3 / (6*pi^2)) ./ ((1.0 + (K_mag * corrL / a).^2).^2.5);
        case 'lowpass'
            kmax = 1.5 * pi / corrL;
            Sk = double(K_mag <= kmax);
        case 'triangular'
            kmax = 3.0 / corrL;
            Sk = max(0.0, 1.0 - K_mag / kmax);
        case 'monodisperse'
            D = 1.2; eta = 0.25; if isfield(extra, 'D'); D = extra.D; end; if isfield(extra, 'eta'); eta = extra.eta; end
            use_vw = true; if isfield(extra, 'use_vw'); use_vw = extra.use_vw; end
            Sk = zeros(size(K_mag));
            for idx = 1:numel(K_mag)
                Sk(idx) = eval_monodisperse_psd_matlab(K_mag(idx), 3, eta, D, use_vw);
            end
        otherwise
            Sk = (corrL^3 / (pi^2)) ./ ((1.0 + (K_mag * corrL).^2).^2.0);
    end
end

function R_theo = get_theoretical_correlation(r, model, corrL, nDim, extra)
    switch lower(model)
        case 'gaussian'
            R_theo = exp(- (r ./ corrL).^2);
        case 'cos'
            R_theo = cos(2.0 * pi * r ./ corrL);
        case 'karman'
            nu = 0.7; if isfield(extra, 'nu'); nu = extra.nu; end
            switch nDim
                case 1; a = 2.0 * sqrt(pi) * exp(gammaln(nu + 0.5) - gammaln(nu));
                case 2; a = 2.0 * sqrt(nu);
                case 3; a = (6.0 * sqrt(pi) * exp(gammaln(nu + 1.5) - gammaln(nu)))^(1/3);
                otherwise; a = 1.0;
            end
            z = (a .* r) ./ corrL;
            R_theo = zeros(size(r));
            R_theo(z < 1e-6) = 1.0;
            pos = z >= 1e-6;
            pref = (2.0^(1.0 - nu)) / gamma(nu);
            R_theo(pos) = pref .* (z(pos).^nu) .* besselk(nu, z(pos));
        case 'exponential'
            R_theo = exp(- r ./ corrL);
        case 'powerlaw'
            switch nDim
                case 1; a = 2.0;
                case 2; a = sqrt(2.0);
                case 3; a = 6.0^(1/3);
                otherwise; a = 1.0;
            end
            R_theo = 1.0 ./ (1.0 + (a .* r ./ corrL).^2).^1.5;
        case 'lowpass'
            switch nDim
                case 1
                    kmax = pi / corrL;
                    z = kmax * r;
                    R_theo = sinc(z / pi);
                case 2
                    kmax = 2.0 / corrL;
                    z = kmax * r;
                    R_theo = ones(size(r));
                    pos = z > 1e-5;
                    R_theo(pos) = 2.0 * besselj(1, z(pos)) ./ z(pos);
                case 3
                    kmax = 1.5 * pi / corrL;
                    u = kmax * r;
                    R_theo = ones(size(r));
                    pos = u > 1e-5;
                    R_theo(pos) = 3.0 * (sin(u(pos)) - u(pos).*cos(u(pos))) ./ (u(pos).^3);
            end
        case 'triangular'
            R_theo = max(0.0, 1.0 - r ./ corrL);
        case 'monodisperse'
            D = 1.2; eta = 0.25; use_vw = true;
            if isfield(extra, 'D'); D = extra.D; end
            if isfield(extra, 'eta'); eta = extra.eta; end
            if isfield(extra, 'use_vw'); use_vw = extra.use_vw; end
            
            k_pts = linspace(1e-4, 25.0 * pi / D, 1000);
            dk = k_pts(2) - k_pts(1);
            psd_vals = zeros(size(k_pts));
            for idx = 1:length(k_pts)
                psd_vals(idx) = eval_monodisperse_psd_matlab(k_pts(idx), nDim, eta, D, use_vw);
            end
            
            R_theo = zeros(size(r));
            for ir = 1:length(r)
                rr = r(ir);
                switch nDim
                    case 1
                        integrand = psd_vals .* cos(k_pts * rr) * 2.0;
                    case 2
                        integrand = psd_vals .* besselj(0, k_pts * rr) .* k_pts * (2.0*pi);
                    case 3
                        if rr < 1e-6
                            sinc_term = ones(size(k_pts));
                        else
                            sinc_term = sin(k_pts * rr) ./ (k_pts * rr);
                        end
                        integrand = psd_vals .* sinc_term .* (k_pts.^2) * (4.0*pi);
                end
                R_theo(ir) = sum(integrand) * dk;
            end
            if R_theo(1) > 0
                R_theo = R_theo / R_theo(1);
            end
        otherwise
            R_theo = exp(- r ./ corrL);
    end
end

function phi = eval_monodisperse_psd_matlab(k, nDim, eta, D, use_vw)
    R = D / 2.0;
    if nDim == 3
        eta_eff = eta; D_eff = D;
        if use_vw && eta > 0.05
            eta_eff = eta * (1.0 - eta / 16.0);
            D_eff = D * (eta_eff / eta)^(1.0/3.0);
        end
        rho_eff = 6.0 * eta_eff / (pi * D_eff^3);
        
        l1 = -(1.0 + 2.0*eta_eff)^2 / (1.0 - eta_eff)^4;
        l2 = 6.0 * eta_eff * (1.0 + 0.5*eta_eff)^2 / (1.0 - eta_eff)^4;
        c_val = -0.5 * eta_eff * (1.0 + 2.0*eta_eff)^2 / (1.0 - eta_eff)^4;
        
        u = k * D_eff;
        if u < 1e-5
            c_k = (4.0*pi*D_eff^3) * (l1/3.0 + l2/4.0 + c_val/6.0);
            m_k = (4.0/3.0) * pi * R^3;
        else
            c_k = (4.0*pi / (k^3)) * ( ...
                l1 * (sin(u) - u*cos(u)) + ...
                l2 * (2.0*u*sin(u) - (u^2 - 2.0)*cos(u) - 2.0)/u + ...
                c_val * ((4.0*u^3 - 24.0*u)*sin(u) - (u^4 - 12.0*u^2 + 24.0)*cos(u) + 24.0)/(u^3) );
            m_k = (4.0 * pi / (k^3)) * (sin(k*R) - (k*R)*cos(k*R));
        end
        
        S_py = 1.0 / max(1e-8, 1.0 - rho_eff * c_k);
        S_tot = S_py;
        
        if use_vw && eta > 0.05
            g_cs = (1.0 - 0.5*eta) / (1.0 - eta)^3;
            g_py = (1.0 + 0.5*eta_eff) / (1.0 - eta_eff)^2;
            A_vw = (g_cs - g_py) * D_eff;
            mu = 24.0 * A_vw / D_eff;
            denom = max(1e-8, mu^2 + k^2);
            I_cos = (mu * cos(k * D_eff) - k * sin(k * D_eff)) / denom;
            I_sin = (mu * sin(k * D_eff) + k * cos(k * D_eff)) / denom;
            dh = (4.0 * pi * A_vw / k) * (D_eff * I_sin + (2.0*mu*I_cos + (mu^2 - k^2)*I_sin/k)/denom);
            S_tot = max(1e-8, S_py + rho_eff * dh);
        end
        
        rho_phys = 6.0 * eta / (pi * D^3);
        phi = rho_phys * (m_k^2) * S_tot / ((2.0*pi)^3 * eta * (1.0 - eta));
    elseif nDim == 1
        u = k * D / 2.0;
        m_k = 2.0 * sin(u) / max(1e-8, k);
        S_1d = (1.0 - eta^2) / max(1e-8, 1.0 - 2.0*eta*cos(k*D) + eta^2);
        phi = (1.0 / D) * (m_k^2) * S_1d / (2.0*pi * (1.0 - eta));
    else
        u = k * D / 2.0;
        m_k = pi * (D/2)^2 * (2.0 * besselj(1, max(1e-8, u)) / max(1e-8, u));
        S_2d = 1.0 / max(1e-8, 1.0 + 4.0*eta * (besselj(1, max(1e-8, k*D)) / max(1e-8, k*D)));
        phi = (4.0 * eta / (pi * D^2)) * (m_k^2) * S_2d / ((2.0*pi)^2 * eta * (1.0 - eta));
    end
end
