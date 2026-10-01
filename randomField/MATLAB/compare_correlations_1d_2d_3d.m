%% =========================================================================
% COMPARE_CORRELATIONS_1D_2D_3D.M
% Measurement, Analysis and Validation of Spatial Correlation Functions in 1D, 2D and 3D
%
% This script computes the empirical spatial autocorrelation function from
% random fields (1D, 2D and 3D) and compares it quantitatively and
% visually with the analytical formulations of the models:
%   1. Gaussian (cm_GAUSSIAN)
%   2. Cosine (cm_COS)
%   3. von Karman (cm_KARMAN) - Fractal with Hurst exponent nu
%   4. Exponential (cm_EXPONENTIAL)
%   5. Power Law (cm_POWERLAW)
%   6. Low-Pass / Sinc (cm_LOWPASS)
%   7. Triangular (cm_TRIANGULAR)
%   8. Monodisperse Spheres (cm_MONODISPERSE) - Percus-Yevick / Verlet-Weis
%
% Measurement Methodology:
%   - Unbiased estimator via FFT (Wiener-Khinchin) with double zero-padding
%     to avoid aliasing / artificial circular periodicity.
%   - Isotropic measurement via radial average in 2D and spherical shells in 3D.
%   - Directional decomposition (X, Y, Z axes and diagonals) to assess
%     the degree of spatial isotropy.
%   - Metrics: RMSE, coefficient of determination R^2, maximum error and
%     empirical correlation length (r at R(r)=1/e).
%
% =========================================================================

clear; close all; clc;

%% 0. Global Settings and Model Selection
fprintf('=========================================================================\n');
fprintf('   AUTOCORRELATION MEASUREMENT AND VALIDATION SYSTEM (1D, 2D, 3D)       \n');
fprintf('=========================================================================\n\n');

% Test options:
% Available models: 'gaussian', 'exponential', 'karman', 'cos', 'powerlaw', 'lowpass', 'triangular', 'monodisperse'
model_1d = 'exponential';
model_2d = 'karman';
model_3d = 'monodisperse';

% Physical parameters of the fields:
corrL_1d = 2.0;    % Correlation length in 1D
corrL_2d = 3.0;    % Correlation length in 2D
corrL_3d = 2.5;    % Correlation length in 3D

% Additional model-specific parameters:
extra_1d.nu = 0.5;                % Hurst exponent for von Karman
extra_2d.nu = 0.7;                % Hurst exponent for von Karman
extra_3d.eta = 0.30;              % Packing fraction (Monodisperse)
extra_3d.D   = 1.5;               % Sphere diameter (Monodisperse)
extra_3d.use_vw = true;           % Verlet-Weis correction

% Size of the test grids:
N1 = 2048;  L1 = 40.0;            % 1D: 2048 points, domain [0, 40]
N2 = 256;   L2 = 30.0;            % 2D: 256 x 256 points, domain [0, 30]^2
N3 = 64;    L3 = 20.0;            % 3D: 64 x 64 x 64 points, domain [0, 20]^3

dx1 = L1 / N1;
dx2 = L2 / N2; dy2 = L2 / N2;
dx3 = L3 / N3; dy3 = L3 / N3; dz3 = L3 / N3;

%% 1. DATA GENERATION OR LOADING

% --- 1.1 1D Field ---
fprintf('-> Generating 1D field (%s, corrL = %.2f, N = %d)...\n', model_1d, corrL_1d, N1);
u1d = generate_synthetic_rf(1, N1, L1, model_1d, corrL_1d, extra_1d);

% --- 1.2 2D Field ---
% Check whether real data saved on disk exists; otherwise generate
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir); script_dir = pwd; end
mat2d_file = fullfile(script_dir, 'random_fields_data.mat');

if exist(mat2d_file, 'file') && false % Change to true to load the saved file
    fprintf('-> Loading 2D field from file %s...\n', mat2d_file);
    loaded_data = load(mat2d_file);
    u2d = loaded_data.field_gaussian;
    N2 = size(u2d, 1);
    dx2 = double(loaded_data.Lx) / N2; dy2 = double(loaded_data.Ly) / size(u2d, 2);
    model_2d = 'gaussian'; corrL_2d = 1.0;
else
    fprintf('-> Generating 2D field (%s, corrL = %.2f, N = %dx%d)...\n', model_2d, corrL_2d, N2, N2);
    u2d = generate_synthetic_rf(2, N2, L2, model_2d, corrL_2d, extra_2d);
end

% --- 1.3 3D Field ---
fprintf('-> Generating 3D field (%s, corrL = %.2f, N = %dx%dx%d)...\n', model_3d, corrL_3d, N3, N3, N3);
u3d = generate_synthetic_rf(3, N3, L3, model_3d, corrL_3d, extra_3d);


%% =========================================================================
%% 2. AUTOCORRELATION MEASUREMENT IN 1D
%% =========================================================================
fprintf('\n-------------------------------------------------------------------------\n');
fprintf(' 1. 1D AUTOCORRELATION MEASUREMENT AND VALIDATION\n');
fprintf('-------------------------------------------------------------------------\n');

max_lag_1d = round(N1 / 4); % Evaluate up to 1/4 of the domain for high sampling reliability
[r1d, R1d_emp] = measure_correlation_1d(u1d, dx1, max_lag_1d);

% Theoretical curve
R1d_theo = get_theoretical_correlation(r1d, model_1d, corrL_1d, 1, extra_1d);
R1d_theo = R1d_theo(:);
R1d_emp  = R1d_emp(:);

% Error metrics
rmse_1d = sqrt(mean((R1d_emp - R1d_theo).^2));
r2_1d   = 1 - sum((R1d_emp - R1d_theo).^2) / sum((R1d_emp - mean(R1d_emp)).^2);
max_err_1d = max(abs(R1d_emp - R1d_theo));
l_emp_1d  = estimate_correlation_length(r1d, R1d_emp);
l_theo_1d = estimate_correlation_length(r1d, R1d_theo);

fprintf('  Model 1D: %s\n', upper(model_1d));
fprintf('  RMSE:                 %.4e\n', rmse_1d);
fprintf('  R^2 (Determination):    %.4f\n', r2_1d);
fprintf('  Max Absolute Error:   %.4e\n', max_err_1d);
fprintf('  Theoretical Length (1/e): %.3f | Measured (1/e): %.3f\n', l_theo_1d, l_emp_1d);


%% =========================================================================
%% 3. AUTOCORRELATION MEASUREMENT IN 2D (Radial Average and Directional Cuts)
%% =========================================================================
fprintf('\n-------------------------------------------------------------------------\n');
fprintf(' 2. 2D AUTOCORRELATION MEASUREMENT AND VALIDATION\n');
fprintf('-------------------------------------------------------------------------\n');

max_lag_2d = round(N2 / 3);
[r2d, R2d_radial, R2d_map, lags_x2, lags_y2, cuts_2d] = measure_correlation_2d(u2d, dx2, dy2, max_lag_2d);

R2d_theo   = get_theoretical_correlation(r2d, model_2d, corrL_2d, 2, extra_2d);
R2d_theo   = R2d_theo(:);
R2d_radial = R2d_radial(:);

rmse_2d = sqrt(mean((R2d_radial - R2d_theo).^2));
r2_2d   = 1 - sum((R2d_radial - R2d_theo).^2) / sum((R2d_radial - mean(R2d_radial)).^2);
max_err_2d = max(abs(R2d_radial - R2d_theo));
l_emp_2d  = estimate_correlation_length(r2d, R2d_radial);
l_theo_2d = estimate_correlation_length(r2d, R2d_theo);

fprintf('  Model 2D: %s\n', upper(model_2d));
fprintf('  RMSE (Radial Average):  %.4e\n', rmse_2d);
fprintf('  R^2 (Determination):    %.4f\n', r2_2d);
fprintf('  Max Absolute Error:   %.4e\n', max_err_2d);
fprintf('  Theoretical Length (1/e): %.3f | Measured (1/e): %.3f\n', l_theo_2d, l_emp_2d);


%% =========================================================================
%% 4. AUTOCORRELATION MEASUREMENT IN 3D (Spherical Average and Cuts)
%% =========================================================================
fprintf('\n-------------------------------------------------------------------------\n');
fprintf(' 3. 3D AUTOCORRELATION MEASUREMENT AND VALIDATION\n');
fprintf('-------------------------------------------------------------------------\n');

max_lag_3d = round(N3 / 3);
[r3d, R3d_spherical, cuts_3d, slice_center_3d, lags_x3, lags_y3] = measure_correlation_3d(u3d, dx3, dy3, dz3, max_lag_3d);

R3d_theo      = get_theoretical_correlation(r3d, model_3d, corrL_3d, 3, extra_3d);
R3d_theo      = R3d_theo(:);
R3d_spherical = R3d_spherical(:);

rmse_3d = sqrt(mean((R3d_spherical - R3d_theo).^2));
r2_3d   = 1 - sum((R3d_spherical - R3d_theo).^2) / sum((R3d_spherical - mean(R3d_spherical)).^2);
max_err_3d = max(abs(R3d_spherical - R3d_theo));
l_emp_3d  = estimate_correlation_length(r3d, R3d_spherical);
l_theo_3d = estimate_correlation_length(r3d, R3d_theo);

fprintf('  Model 3D: %s\n', upper(model_3d));
fprintf('  RMSE (Spherical Average): %.4e\n', rmse_3d);
fprintf('  R^2 (Determination):    %.4f\n', r2_3d);
fprintf('  Max Absolute Error:   %.4e\n', max_err_3d);
fprintf('  Theoretical Length (1/e): %.3f | Measured (1/e): %.3f\n', l_theo_3d, l_emp_3d);


%% =========================================================================
%% 5. COMPLETE GRAPHICAL VISUALIZATION
%% =========================================================================

% Modern color palette
c_emp = [0.0, 0.45, 0.74];   % Dark blue
c_theo = [0.85, 0.33, 0.1];  % Vibrant orange
c_res = [0.47, 0.67, 0.19];  % Moss green
c_cut1 = [0.49, 0.18, 0.56]; % Purple
c_cut2 = [0.93, 0.69, 0.13]; % Golden yellow

% -------------------------------------------------------------------------
% FIGURE 1: 1D VALIDATION
% -------------------------------------------------------------------------
fig1 = figure('Name', '1D Correlation Validation', 'Color', 'w', ...
       'Units', 'normalized', 'Position', [0.05, 0.52, 0.42, 0.38]);

subplot(2, 2, [1, 2]);
x_coords = (0:N1-1)*dx1;
plot(x_coords, u1d, 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.0);
grid on; box on;
title(sprintf('1D Random Field (%s, L = %.1f)', upper(model_1d), corrL_1d), 'FontWeight', 'bold');
xlabel('Position x'); ylabel('u(x)');
xlim([0, min(L1, 20)]);

subplot(2, 2, 3);
plot(r1d, R1d_emp, 'o-', 'Color', c_emp, 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'Empirical (FFT)');
hold on;
plot(r1d, R1d_theo, '--', 'Color', c_theo, 'LineWidth', 2.0, 'DisplayName', 'Exact Theory');
yline(exp(-1), ':', 'Color', [0.5 0.5 0.5], 'DisplayName', 'R = 1/e');
grid on; box on;
xlabel('Lag Distance r'); ylabel('Autocorrelation R(r)');
title(sprintf('1D Autocorrelation (R^2 = %.3f, RMSE = %.2e)', r2_1d, rmse_1d), 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, max(r1d)]); ylim([-0.2, 1.05]);

subplot(2, 2, 4);
plot(r1d, R1d_emp - R1d_theo, '-', 'Color', c_res, 'LineWidth', 1.5);
grid on; box on;
xlabel('Lag Distance r'); ylabel('Residual [R_{emp} - R_{theo}]');
title('Correlation Residual Error', 'FontWeight', 'bold');
yline(0, '--k');
xlim([0, max(r1d)]);


% -------------------------------------------------------------------------
% FIGURE 2: 2D VALIDATION
% -------------------------------------------------------------------------
fig2 = figure('Name', '2D Correlation Validation', 'Color', 'w', ...
       'Units', 'normalized', 'Position', [0.50, 0.52, 0.45, 0.38]);

subplot(2, 2, 1);
imagesc((0:N2-1)*dx2, (0:N2-1)*dy2, u2d);
axis image; colormap(subplot(2,2,1), viridis_fallback()); colorbar;
title(sprintf('2D Field (%s, L = %.1f)', upper(model_2d), corrL_2d), 'FontWeight', 'bold');
xlabel('x'); ylabel('y');

subplot(2, 2, 2);
imagesc(lags_x2, lags_y2, R2d_map);
axis image; colormap(subplot(2,2,2), jet_fallback()); colorbar;
title('2D Map R(\Delta x, \Delta y)', 'FontWeight', 'bold');
xlabel('\Delta x'); ylabel('\Delta y');

subplot(2, 2, 3);
plot(r2d, R2d_radial, 'o-', 'Color', c_emp, 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'Radial Average');
hold on;
plot(r2d, R2d_theo, '--', 'Color', c_theo, 'LineWidth', 2.0, 'DisplayName', 'Exact Theory');
yline(exp(-1), ':', 'Color', [0.5 0.5 0.5], 'DisplayName', 'R = 1/e');
grid on; box on;
xlabel('Lag Radius r'); ylabel('Autocorrelation R(r)');
title(sprintf('Radial Average (R^2 = %.3f)', r2_2d), 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, max(r2d)]); ylim([-0.2, 1.05]);

subplot(2, 2, 4);
plot(r2d, R2d_theo, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Theory');
hold on;
plot(cuts_2d.r_x, cuts_2d.R_x, '-', 'Color', c_emp, 'LineWidth', 1.2, 'DisplayName', 'X Cut (\Delta y=0)');
plot(cuts_2d.r_y, cuts_2d.R_y, '-.', 'Color', c_cut1, 'LineWidth', 1.2, 'DisplayName', 'Y Cut (\Delta x=0)');
plot(cuts_2d.r_diag, cuts_2d.R_diag, ':', 'Color', c_cut2, 'LineWidth', 1.5, 'DisplayName', 'Diagonal (45^\circ)');
grid on; box on;
xlabel('Lag Radius r'); ylabel('R(r)');
title('Isotropy Test (Cuts)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, max(r2d)]); ylim([-0.2, 1.05]);


% -------------------------------------------------------------------------
% FIGURE 3: 3D VALIDATION
% -------------------------------------------------------------------------
fig3 = figure('Name', '3D Correlation Validation', 'Color', 'w', ...
       'Units', 'normalized', 'Position', [0.20, 0.06, 0.55, 0.40]);

subplot(1, 3, 1);
% Central slices of the 3D volume
slice_z = round(N3/2);
imagesc((0:N3-1)*dx3, (0:N3-1)*dy3, u3d(:,:,slice_z));
axis image; colormap(subplot(1,3,1), viridis_fallback()); colorbar;
title(sprintf('Central Z Slice of 3D Field\n(%s)', upper(model_3d)), 'FontWeight', 'bold');
xlabel('x'); ylabel('y');

subplot(1, 3, 2);
imagesc(lags_x3, lags_y3, slice_center_3d);
axis image; colormap(subplot(1,3,2), jet_fallback()); colorbar;
title('Central Slice R_{3D}(\Delta x, \Delta y, 0)', 'FontWeight', 'bold');
xlabel('\Delta x'); ylabel('\Delta y');

subplot(1, 3, 3);
plot(r3d, R3d_spherical, 'o-', 'Color', c_emp, 'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', 'Spherical Average');
hold on;
plot(r3d, R3d_theo, '--', 'Color', c_theo, 'LineWidth', 2.0, 'DisplayName', 'Exact Theory');
plot(cuts_3d.r_x, cuts_3d.R_x, ':', 'Color', c_cut1, 'LineWidth', 1.1, 'DisplayName', 'X Axis');
plot(cuts_3d.r_y, cuts_3d.R_y, '-.', 'Color', c_cut2, 'LineWidth', 1.1, 'DisplayName', 'Y Axis');
plot(cuts_3d.r_z, cuts_3d.R_z, '--', 'Color', [0.3 0.6 0.3], 'LineWidth', 1.1, 'DisplayName', 'Z Axis');
grid on; box on;
xlabel('Lag Radius r'); ylabel('Autocorrelation R(r)');
title(sprintf('3D Spherical Average\n(R^2 = %.3f, RMSE = %.2e)', r2_3d, rmse_3d), 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, max(r3d)]); ylim([-0.25, 1.05]);

fprintf('\n>>> Graphical validation completed successfully! (3 figures open)\n\n');


%% =========================================================================
%% AUXILIARY AUTOCORRELATION MEASUREMENT FUNCTIONS
%% =========================================================================

function [r, R_emp] = measure_correlation_1d(u, dx, max_lag)
    % Measures the 1D autocorrelation function using an unbiased FFT
    u = u(:);
    N = length(u);
    if nargin < 3 || isempty(max_lag)
        max_lag = floor(N / 2);
    end
    
    u_cent = u - mean(u);
    var_u = mean(u_cent.^2);
    if var_u < 1e-12; var_u = 1.0; end
    
    % Double zero-padding to avoid circular periodicity
    nfft = 2 * N;
    U = fft(u_cent, nfft);
    R_full = real(ifft(abs(U).^2));
    
    % Unbiased estimator: normalizes by the effective number of pairs at each lag
    lags = (0:max_lag)';
    weights = N - lags;
    
    R_emp = (R_full(1:max_lag+1) ./ weights) / var_u;
    r = lags * dx;
end

function [r, R_radial, R_map, lags_x, lags_y, cuts] = measure_correlation_2d(u, dx, dy, max_lag)
    % Measures the 2D autocorrelation matrix and its isotropic radial average
    [Nx, Ny] = size(u);
    if nargin < 4 || isempty(max_lag)
        max_lag = min(floor(Nx/2), floor(Ny/2));
    end
    
    u_cent = u - mean(u(:));
    var_u = mean(u_cent(:).^2);
    if var_u < 1e-12; var_u = 1.0; end
    
    nfft_x = 2 * Nx;
    nfft_y = 2 * Ny;
    
    U = fft2(u_cent, nfft_x, nfft_y);
    R_full = real(fftshift(ifft2(abs(U).^2)));
    
    % Coordinates and triangular overlap weights
    cx = Nx + 1;
    cy = Ny + 1;
    
    lag_x_idx = -max_lag:max_lag;
    lag_y_idx = -max_lag:max_lag;
    
    [LX, LY] = meshgrid(lag_x_idx, lag_y_idx);
    W = (Nx - abs(LX)) .* (Ny - abs(LY));
    W(W <= 0) = 1;
    
    R_map = R_full(cy + lag_y_idx, cx + lag_x_idx) ./ (W * var_u);
    
    lags_x = lag_x_idx * dx;
    lags_y = lag_y_idx * dy;
    
    % Radial average over rings
    [X_dist, Y_dist] = meshgrid(lags_x, lags_y);
    Dist = sqrt(X_dist.^2 + Y_dist.^2);
    
    dr = min(dx, dy);
    max_r = max_lag * dr;
    r_edges = 0:dr:max_r;
    r = (r_edges(1:end-1) + r_edges(2:end)) / 2;
    r = [0, r];
    
    R_radial = zeros(size(r));
    R_radial(1) = 1.0; % By definition at lag 0
    
    for k = 2:length(r)
        mask = (Dist >= r_edges(k-1)) & (Dist < r_edges(k));
        if any(mask(:))
            R_radial(k) = mean(R_map(mask));
        else
            R_radial(k) = NaN;
        end
    end
    
    % Interpolate any empty NaNs
    nan_idx = isnan(R_radial);
    if any(nan_idx)
        R_radial(nan_idx) = interp1(r(~nan_idx), R_radial(~nan_idx), r(nan_idx), 'linear', 'extrap');
    end
    
    % Directional cuts
    idx_zero_y = find(lag_y_idx == 0);
    idx_zero_x = find(lag_x_idx == 0);
    
    cuts.r_x = lags_x(lags_x >= 0);
    cuts.R_x = R_map(idx_zero_y, lag_x_idx >= 0);
    
    cuts.r_y = lags_y(lags_y >= 0);
    cuts.R_y = R_map(lag_y_idx >= 0, idx_zero_x)';
    
    diag_lags = 0:max_lag;
    cuts.r_diag = diag_lags * sqrt(dx^2 + dy^2);
    cuts.R_diag = zeros(size(diag_lags));
    for d = 1:length(diag_lags)
        cuts.R_diag(d) = R_map(idx_zero_y + diag_lags(d), idx_zero_x + diag_lags(d));
    end
end

function [r, R_spherical, cuts, slice_center, lags_x, lags_y] = measure_correlation_3d(u, dx, dy, dz, max_lag)
    % Measures the 3D autocorrelation and computes the average over spherical shells
    [Nx, Ny, Nz] = size(u);
    if nargin < 5 || isempty(max_lag)
        max_lag = min([floor(Nx/2), floor(Ny/2), floor(Nz/2)]);
    end
    
    u_cent = u - mean(u(:));
    var_u = mean(u_cent(:).^2);
    if var_u < 1e-12; var_u = 1.0; end
    
    nfft = [2*Nx, 2*Ny, 2*Nz];
    U = fftn(u_cent, nfft);
    R_full = real(fftshift(fftn(abs(U).^2))) / prod(nfft); % Inverse via normalized fftn
    
    cx = Nx + 1; cy = Ny + 1; cz = Nz + 1;
    
    lag_x = -max_lag:max_lag;
    lag_y = -max_lag:max_lag;
    lag_z = -max_lag:max_lag;
    
    [LX, LY, LZ] = meshgrid(lag_x, lag_y, lag_z);
    W = (Nx - abs(LX)) .* (Ny - abs(LY)) .* (Nz - abs(LZ));
    W(W <= 0) = 1;
    
    sub_R = R_full(cy + lag_y, cx + lag_x, cz + lag_z);
    sub_R = sub_R ./ (W * var_u);
    
    % Central slice in the Z = 0 plane
    slice_center = sub_R(:, :, max_lag + 1);
    lags_x = lag_x * dx;
    lags_y = lag_y * dy;
    
    % 3D distances for the spherical average
    Dist3D = sqrt((LX * dx).^2 + (LY * dy).^2 + (LZ * dz).^2);
    
    dr = min([dx, dy, dz]);
    max_r = max_lag * dr;
    r_edges = 0:dr:max_r;
    r = [0, (r_edges(1:end-1) + r_edges(2:end)) / 2];
    
    R_spherical = zeros(size(r));
    R_spherical(1) = 1.0;
    
    for k = 2:length(r)
        mask = (Dist3D >= r_edges(k-1)) & (Dist3D < r_edges(k));
        if any(mask(:))
            R_spherical(k) = mean(sub_R(mask));
        else
            R_spherical(k) = NaN;
        end
    end
    
    nan_idx = isnan(R_spherical);
    if any(nan_idx)
        R_spherical(nan_idx) = interp1(r(~nan_idx), R_spherical(~nan_idx), r(nan_idx), 'linear', 'extrap');
    end
    
    % 3D cuts
    pos_idx = max_lag + 1;
    cuts.r_x = (0:max_lag) * dx;
    cuts.R_x = squeeze(sub_R(pos_idx, pos_idx:end, pos_idx));
    
    cuts.r_y = (0:max_lag) * dy;
    cuts.R_y = squeeze(sub_R(pos_idx:end, pos_idx, pos_idx))';
    
    cuts.r_z = (0:max_lag) * dz;
    cuts.R_z = squeeze(sub_R(pos_idx, pos_idx, pos_idx:end));
end


%% =========================================================================
%% ANALYTICAL THEORETICAL CORRELATION FUNCTIONS
%% =========================================================================

function R_theo = get_theoretical_correlation(r, model_name, corrL, nDim, extra)
    % Returns the exact theoretical autocorrelation function for the model and dimension
    r = r(:)';
    switch lower(model_name)
        case 'gaussian'
            % In the Fortran code: SkVec ~ exp(- (k*corrL)^2 / 8)
            % The normalized inverse transform gives R(r) = exp(- 2 * (r / corrL)^2)
            R_theo = exp(- 2.0 * (r ./ corrL).^2);
            
        case 'exponential'
            % C(r) = exp(- a * r / corrL)
            switch nDim
                case 1; a = 2.0;
                case 2; a = sqrt(2.0);
                case 3; a = 6.0^(1/3);
                otherwise; a = 1.0;
            end
            R_theo = exp(- a * r ./ corrL);
            
        case 'karman'
            % von Karman: R(r) = (2^(1-nu) / Gamma(nu)) * (a*r/L)^nu * K_nu(a*r/L)
            nu = 0.5;
            if isfield(extra, 'nu'); nu = extra.nu; end
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
            
            % Evaluation via besselk
            prefactor = (2.0^(1.0 - nu)) / gamma(nu);
            R_theo(pos) = prefactor .* (z(pos).^nu) .* besselk(nu, z(pos));
            R_theo(isnan(R_theo) | isinf(R_theo)) = 0.0;
            
        case 'cos'
            % Quasi-periodic cosine
            R_theo = cos(2.0 * pi * r ./ corrL);
            
        case 'powerlaw'
            % Power law
            switch nDim
                case 1; a = 2.0;
                case 2; a = sqrt(2.0);
                case 3; a = 6.0^(1/3);
                otherwise; a = 1.0;
            end
            R_theo = 1.0 ./ (1.0 + (a .* r ./ corrL).^2).^1.5;
            
        case 'lowpass'
            % Analytical sinc per dimension
            switch nDim
                case 1
                    % Sinc 1D: sin(kmax*r)/(kmax*r)
                    kmax = pi / corrL;
                    z = kmax * r;
                    R_theo = sinc(z / pi);
                case 2
                    % 2D: 2 * J1(kmax*r) / (kmax*r)
                    kmax = 2.0 / corrL;
                    z = kmax * r;
                    R_theo = ones(size(r));
                    pos = z > 1e-5;
                    R_theo(pos) = 2.0 * besselj(1, z(pos)) ./ z(pos);
                case 3
                    % 3D: 3 * (sin(u) - u*cos(u)) / u^3
                    kmax = 1.5 * pi / corrL;
                    u = kmax * r;
                    R_theo = ones(size(r));
                    pos = u > 1e-5;
                    R_theo(pos) = 3.0 * (sin(u(pos)) - u(pos).*cos(u(pos))) ./ (u(pos).^3);
            end
            
        case 'triangular'
            % Bartlett / Triangular in 1D
            R_theo = max(0, 1.0 - r ./ corrL);
            
        case 'monodisperse'
            % Percus-Yevick with Verlet-Weis: inverse transform of the analytical PSD
            D = 1.0; eta = 0.25; use_vw = true;
            if isfield(extra, 'D'); D = extra.D; end
            if isfield(extra, 'eta'); eta = extra.eta; end
            if isfield(extra, 'use_vw'); use_vw = extra.use_vw; end
            
            % Fourier-Bessel spectral integration for exact precision
            k_pts = linspace(1e-4, 25.0 * pi / D, 1500);
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
                R_theo = R_theo / R_theo(1); % Normalization R(0) = 1
            end
            
        otherwise
            R_theo = exp(-r ./ corrL);
    end
end

function phi = eval_monodisperse_psd_matlab(k, nDim, eta, D, use_vw)
    % Evaluates the analytical Percus-Yevick + Verlet-Weis spectral density
    R = D / 2.0;
    if nDim == 3
        eta_eff = eta;
        D_eff = D;
        if use_vw && eta > 0.05
            eta_eff = eta * (1.0 - eta / 16.0);
            D_eff = D * (eta_eff / eta)^(1.0/3.0);
        end
        R_eff = D_eff / 2.0;
        rho_eff = 6.0 * eta_eff / (pi * D_eff^3);
        
        % Wertheim-Thiele coefficients
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
        rho = eta / D;
        u = k * D / 2.0;
        if u < 1e-6
            m_k = D;
        else
            m_k = D * sin(u) / u;
        end
        S_tot = (1.0 - eta^2) / max(1e-8, (1.0 - 2.0*eta*cos(k*D) + eta^2));
        phi = rho * (m_k^2) * S_tot / (2.0 * pi * eta * (1.0 - eta));
    else
        % 2D Disks
        rho = eta / (pi * R^2);
        u = k * R;
        if u < 1e-6
            m_k = pi * R^2;
            j1_term = R / 2.0;
        else
            m_k = 2.0 * pi * R * besselj(1, u) / k;
            j1_term = besselj(1, 2.0*u) / (2.0*k);
        end
        c2d = - (1.0 + 2.0*eta) / (1.0 - eta)^3 * (2.0 * pi * R * j1_term);
        S_tot = 1.0 / max(1e-8, 1.0 - rho * c2d);
        phi = rho * (m_k^2) * S_tot / ((2.0*pi)^2 * eta * (1.0 - eta));
    end
end


%% =========================================================================
%% FAST SPECTRAL RANDOM FIELD GENERATOR FOR TESTING
%% =========================================================================

function u = generate_synthetic_rf(nDim, N, L, model_name, corrL, extra)
    % Generates an isotropic random field via the spectral filtering method
    dx = L / N;
    dk = 2.0 * pi / L;
    
    switch nDim
        case 1
            kx = (0:N-1)';
            kx(kx > N/2) = kx(kx > N/2) - N;
            k_mag = abs(kx) * dk;
            
            % Evaluate the theoretical spectrum S(k)
            Sk = eval_spectrum_nd(k_mag, nDim, model_name, corrL, extra);
            
            % Linear filtering of real Gaussian white noise
            white_noise = randn(N, 1);
            U_k = fft(white_noise) .* sqrt(Sk(:));
            u = real(ifft(U_k));
            
        case 2
            kx = (0:N-1); kx(kx > N/2) = kx(kx > N/2) - N;
            ky = (0:N-1); ky(ky > N/2) = ky(ky > N/2) - N;
            [KX, KY] = meshgrid(kx, ky);
            K_mag = sqrt(KX.^2 + KY.^2) * dk;
            
            Sk = eval_spectrum_nd(K_mag, nDim, model_name, corrL, extra);
            
            white_noise = randn(N, N);
            U_k = fft2(white_noise) .* sqrt(Sk);
            u = real(ifft2(U_k));
            
        case 3
            kx = (0:N-1); kx(kx > N/2) = kx(kx > N/2) - N;
            [KX, KY, KZ] = meshgrid(kx, kx, kx);
            K_mag = sqrt(KX.^2 + KY.^2 + KZ.^2) * dk;
            
            Sk = eval_spectrum_nd(K_mag, nDim, model_name, corrL, extra);
            
            white_noise = randn(N, N, N);
            U_k = fftn(white_noise) .* sqrt(Sk);
            u = real(ifftn(U_k));
    end
    
    % Normalization to zero mean and unit variance
    u = u - mean(u(:));
    std_u = std(u(:));
    if std_u > 1e-12; u = u / std_u; end
end

function Sk = eval_spectrum_nd(k_mag, nDim, model_name, corrL, extra)
    % Evaluates the analytical power spectrum S(k)
    switch lower(model_name)
        case 'gaussian'
            Sk = (corrL^nDim) * exp(-(k_mag .* corrL).^2 / 8.0);
        case 'exponential'
            switch nDim
                case 1; a = 2.0; Sk = (corrL / (2.0*pi)) ./ (1.0 + (k_mag*corrL/a).^2);
                case 2; a = sqrt(2.0); Sk = (corrL^2 / (4.0*pi)) ./ ((1.0 + (k_mag*corrL/a).^2).^1.5);
                case 3; a = 6.0^(1/3); Sk = (corrL^3 / (6.0*pi^2)) ./ ((1.0 + (k_mag*corrL/a).^2).^2.0);
            end
        case 'karman'
            nu = 0.5; if isfield(extra, 'nu'); nu = extra.nu; end
            switch nDim
                case 1
                    a = 2.0 * sqrt(pi) * exp(gammaln(nu + 0.5) - gammaln(nu));
                    Sk = (corrL / (2.0*pi)) ./ ((1.0 + (k_mag*corrL/a).^2).^(nu + 0.5));
                case 2
                    a = 2.0 * sqrt(nu);
                    Sk = (corrL^2 / (4.0*pi)) ./ ((1.0 + (k_mag*corrL/a).^2).^(nu + 1.0));
                case 3
                    a = (6.0 * sqrt(pi) * exp(gammaln(nu + 1.5) - gammaln(nu)))^(1/3);
                    Sk = (corrL^3 / (6.0*pi^2)) ./ ((1.0 + (k_mag*corrL/a).^2).^(nu + 1.5));
            end
        case 'cos'
            k0 = 2.0 * pi / corrL;
            sigma_k = 0.1 * k0;
            Sk = exp(-((k_mag - k0).^2) / (2.0 * sigma_k^2));
        case 'monodisperse'
            D = 1.0; eta = 0.25; use_vw = true;
            if isfield(extra, 'D'); D = extra.D; end
            if isfield(extra, 'eta'); eta = extra.eta; end
            if isfield(extra, 'use_vw'); use_vw = extra.use_vw; end
            Sk = zeros(size(k_mag));
            for i = 1:numel(k_mag)
                Sk(i) = eval_monodisperse_psd_matlab(k_mag(i), nDim, eta, D, use_vw);
            end
        otherwise
            Sk = (corrL^nDim) * exp(-(k_mag .* corrL).^2 / 8.0);
    end
end

function l = estimate_correlation_length(r, R)
    % Estimates the correlation length where R(r) reaches 1/e (~0.3679)
    target = exp(-1);
    idx = find(R <= target, 1, 'first');
    if isempty(idx) || idx == 1
        l = NaN;
    else
        % Fine linear interpolation
        r1 = r(idx-1); r2 = r(idx);
        R1 = R(idx-1); R2 = R(idx);
        l = r1 + (target - R1) * (r2 - r1) / (R2 - R1);
    end
end

function cmap = viridis_fallback()
    % Modern blue-green-yellow colormap
    try
        cmap = viridis(256);
    catch
        cmap = parula(256);
    end
end

function cmap = jet_fallback()
    cmap = jet(256);
end
