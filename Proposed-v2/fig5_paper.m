%% Fig.5 (paper) reproduction: EH power vs transmit power
% This script generates 4 curves:
%   - Proposed 6DMHS (rotation + translation)
%   - 6DMHS with rotation only
%   - 6DMHS with translation only
%   - FPA (fixed position & rotation)
%
% Notes:
% - Requires CVX for fun_run_R00.m (and optimize_system.m for translation schemes).
% - This script can be time-consuming because it solves many convex programs.

clc; clear; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

rng(1, 'twister');

%% Sweep setup (match the paper figure style)
Pt_range = 30:1:45;              % dBm
Pt_ref_for_stage2 = 40;          % dBm (used to design translation/holographic beamforming once)

%% System parameters (paper)
K_rice = 10;                     % linear K-factor
kappa0 = 1;                      % sensing spacing control (2^kappa0)
N_fft = 256;                     % FFT size in sensing
sigma0 = 1e-12;                  % noise power (W)
R00 = 1;                         % throughput constraint (bit/s/Hz), used in fun_run_R00

fc = 30e9;
c = 3e8;
lambda = c / fc;
dx = lambda / 4;
dy = lambda / 4;

Mx = 32; My = 32; M = Mx * My;
K = 4;  B = 4;                   % paper uses K=B

Ptx_ref = 10^((Pt_ref_for_stage2 - 30) / 10);
dmin = 0.1;

%% Geometry (RHS elements & sensing ring)
[geom] = build_geometry(Mx, My, dx, dy, lambda, kappa0);

%% Stage I: initial deployment (FPA-like) + channel + sensing
[Coor_Ele0, normal0, Rot0, ~, Coor_Ele_r0, Coor_Ele_c0, ...
    Coor_Ele_x_d0, Coor_Ele_x_u0, Coor_Ele_y_l0, Coor_Ele_y_r0] = ...
    Orientation_Initial(B, geom.Coor_Ele_init, geom.Coor_Ele_init_r, geom.Coor_Ele_init_c, 0, ...
    geom.Coor_Ele_init_x_d, geom.Coor_Ele_init_x_u, geom.Coor_Ele_init_y_l, geom.Coor_Ele_init_y_r);

[~, e0, e2, ~, ~, ~, ~, Iota, eta, Omega, h_r, h_c, h_sen, ...
    h_x_d, h_x_u, h_y_l, h_y_r, Omega0] = ...
    Channel_Generation_init(K, B, M, Mx, My, geom.M_x_r, geom.M_y_c, normal0, K_rice, ...
    Coor_Ele0, Coor_Ele_r0, Coor_Ele_c0, lambda, Rot0, ...
    Coor_Ele_x_d0, Coor_Ele_x_u0, Coor_Ele_y_l0, Coor_Ele_y_r0);

% Sensing_Algorithm returns [est_G, est_local, est_error]; use est_G (global) as in existing scripts.
[est, ~, ~] = Sensing_Algorithm(K, B, Iota, Mx, My, geom.M_x_r, geom.M_y_c, ...
    h_r, h_c, Rot0, e0, K_rice, lambda, e2, sigma0, Omega, eta, N_fft, geom.d, kappa0, ...
    h_sen, h_x_d, h_x_u, h_y_l, h_y_r, ...
    geom.Coor_Ele_init_x_d, geom.Coor_Ele_init_x_u, geom.Coor_Ele_init_y_l, geom.Coor_Ele_init_y_r, ...
    Omega0);

%% Build 4 configurations once (Stage II): {s, mm, Rot, normal_rot, Coor_Ele_noTrans} -> H_eff

% (1) Fixed (FPA) orientation
[Rot_fixed, normal_fixed, Coor_Ele_fixed] = build_rot_and_elements_fixed(Rot0, geom.Coor_Ele_init);

% (2) Rotation alignment (for proposed & rotation-only)
[Rot_align, normal_align, Coor_Ele_align] = build_rot_and_elements_align(est, geom.Coor_Ele_init);

% Fixed position selection for "FPA" and "rotation only"
s_fixed = eye(B, B); % B1=B=4, each RHS picks one fixed point on the sphere

mm_fpa = build_mm_for_fixed_positions(s_fixed, Rot_fixed, normal_fixed, est, geom, lambda);
H_eff_fpa = generate_effective_channel(s_fixed, K, B, M, Iota, Coor_Ele_fixed, ...
    est, normal_fixed, Rot_fixed, geom.Coor_Ele_init, e0, e2, eta, K_rice, Omega0, ...
    lambda, geom.delta_x, geom.delta_y, dx, dy, mm_fpa, geom.V_F);

mm_rot_only = build_mm_for_fixed_positions(s_fixed, Rot_align, normal_align, est, geom, lambda);
H_eff_rot_only = generate_effective_channel(s_fixed, K, B, M, Iota, Coor_Ele_align, ...
    est, normal_align, Rot_align, geom.Coor_Ele_init, e0, e2, eta, K_rice, Omega0, ...
    lambda, geom.delta_x, geom.delta_y, dx, dy, mm_rot_only, geom.V_F);

% Translation-only: optimize positions with fixed rotation
[s_trans_only, ~, mm_trans_only, ~] = optimize_system(B, K, M, Mx, My, est, ...
    normal_fixed, Rot_fixed, geom.Coor_Ele_init, lambda, dx, dy, ...
    geom.delta_x, geom.delta_y, dmin, geom.V_F0, geom.V_F, Ptx_ref);

H_eff_trans_only = generate_effective_channel(s_trans_only, K, B, M, Iota, Coor_Ele_fixed, ...
    est, normal_fixed, Rot_fixed, geom.Coor_Ele_init, e0, e2, eta, K_rice, Omega0, ...
    lambda, geom.delta_x, geom.delta_y, dx, dy, mm_trans_only, geom.V_F);

% Proposed: optimize positions with aligned rotation
[s_prop, ~, mm_prop, ~] = optimize_system(B, K, M, Mx, My, est, ...
    normal_align, Rot_align, geom.Coor_Ele_init, lambda, dx, dy, ...
    geom.delta_x, geom.delta_y, dmin, geom.V_F0, geom.V_F, Ptx_ref);

H_eff_prop = generate_effective_channel(s_prop, K, B, M, Iota, Coor_Ele_align, ...
    est, normal_align, Rot_align, geom.Coor_Ele_init, e0, e2, eta, K_rice, Omega0, ...
    lambda, geom.delta_x, geom.delta_y, dx, dy, mm_prop, geom.V_F);

% Normalize once (kept consistent with existing scripts)
H_eff_fpa = H_eff_fpa / sqrt(sigma0);
H_eff_rot_only = H_eff_rot_only / sqrt(sigma0);
H_eff_trans_only = H_eff_trans_only / sqrt(sigma0);
H_eff_prop = H_eff_prop / sqrt(sigma0);

%% Stage III: SWIPT optimization (vary Pt)
P_EH_fpa = zeros(size(Pt_range));
P_EH_rot_only = zeros(size(Pt_range));
P_EH_trans_only = zeros(size(Pt_range));
P_EH_prop = zeros(size(Pt_range));

cacheFile = fullfile(thisDir, 'fig5_paper_results.mat');

for idx = 1:numel(Pt_range)
    Pt = Pt_range(idx);
    Ptx = 10^((Pt - 30) / 10);

    % Deterministic initialization per Pt for reproducibility
    rng(1000 + idx, 'twister');
    psi1 = 0.1 * rand(K, K);
    psi2 = 0.01 * rand(K, 1) + 0.01j * rand(K, 1);
    rho = 0.1 * rand(K, 1);

    fprintf('[%02d/%02d] Pt=%g dBm ...\n', idx, numel(Pt_range), Pt);

    P_EH_prop(idx) = fun_run_R00(R00, H_eff_prop, rho, psi1, psi2, sigma0, Ptx, B, K);
    P_EH_rot_only(idx) = fun_run_R00(R00, H_eff_rot_only, rho, psi1, psi2, sigma0, Ptx, B, K);
    P_EH_trans_only(idx) = fun_run_R00(R00, H_eff_trans_only, rho, psi1, psi2, sigma0, Ptx, B, K);
    P_EH_fpa(idx) = fun_run_R00(R00, H_eff_fpa, rho, psi1, psi2, sigma0, Ptx, B, K);

    save(cacheFile, 'Pt_range', 'P_EH_prop', 'P_EH_rot_only', 'P_EH_trans_only', 'P_EH_fpa');
end

%% Plot
figure('Position', [100 100 860 520]);
hold on;
plot(Pt_range, P_EH_prop, '-o', 'LineWidth', 2, 'MarkerSize', 7, ...
    'Color', [0 0.4470 0.7410]);
plot(Pt_range, P_EH_rot_only, '-v', 'LineWidth', 2, 'MarkerSize', 7, ...
    'Color', [0.8500 0.3250 0.0980]);
plot(Pt_range, P_EH_trans_only, '-d', 'LineWidth', 2, 'MarkerSize', 7, ...
    'Color', [0.4940 0.1840 0.5560]);
plot(Pt_range, P_EH_fpa, '-s', 'LineWidth', 2, 'MarkerSize', 7, ...
    'Color', [0.4660 0.6740 0.1880]);
grid on;
box on;
xlabel('Transmit power (dBm)');
ylabel('EH power (W)');
legend('Proposed 6DMHS', '6DMHS with rotation only', '6DMHS with translation only', 'FPA', ...
    'Location', 'northwest');

ax = gca;
ax.YAxis.Exponent = -3;

saveas(gcf, fullfile(thisDir, 'fig5_paper.fig'));
saveas(gcf, fullfile(thisDir, 'fig5_paper.png'));

%% Local helpers
function geom = build_geometry(Mx, My, dx, dy, lambda, kappa0)
    M = Mx * My;

    delta_x = (2 * (0:(Mx-1)).' - Mx + 1) / 2;
    delta_y = (2 * (0:(My-1)).' - My + 1) / 2;
    Coor_Ele_init = [kron(delta_x * dx, ones(My, 1)), kron(ones(Mx, 1), delta_y * dy), zeros(M, 1)];

    Coor_Ele_init_r = [delta_x * dx, (-Mx + 1) / 2 * (lambda / 2) * ones(Mx, 1), zeros(Mx, 1)];
    Coor_Ele_init_c = [(My - 1) / 2 * dx * ones(My, 1), delta_y * dy, zeros(My, 1)];

    delta_ele = 2^kappa0;
    d = dx * delta_ele;
    Coor_Ele_init_r = Coor_Ele_init_r(1:delta_ele:end, :);
    Coor_Ele_init_c = Coor_Ele_init_c(1:delta_ele:end, :);
    M_x_r = size(Coor_Ele_init_r, 1);
    M_y_c = size(Coor_Ele_init_c, 1);

    delta_x_r = [(-M_x_r + 1) / 2 : 1 : (M_x_r - 1) / 2].';
    delta_y_c = [(-M_y_c + 1) / 2 : 1 : (M_y_c - 1) / 2].';

    Coor_Ele_init_x_d = [delta_x_r * (lambda / 2), (-M_y_c + 1) / 2 * (lambda / 2) * ones(M_x_r, 1), zeros(M_x_r, 1)];
    Coor_Ele_init_x_u = [delta_x_r * (lambda / 2), (M_y_c - 1) / 2 * (lambda / 2) * ones(M_x_r, 1), zeros(M_x_r, 1)];
    Coor_Ele_init_y_l = [(-M_x_r + 1) / 2 * (lambda / 2) * ones(M_y_c, 1), delta_y_c * (lambda / 2), zeros(M_y_c, 1)];
    Coor_Ele_init_y_r = [(M_x_r - 1) / 2 * (lambda / 2) * ones(M_y_c, 1), delta_y_c * (lambda / 2), zeros(M_y_c, 1)];

    % Feed response model
    Coor_Feed = [0.002984305671304, 0.001586131606604, 0];
    Dis_Feed2Ele = sqrt((Coor_Ele_init - Coor_Feed).^2 * ones(3, 1));
    V_F0 = exp(-1j * 2 * pi * sqrt(3) / lambda * Dis_Feed2Ele);
    eta0 = 8 / 3 / M;
    V_F = sqrt(eta0) * V_F0;

    geom = struct();
    geom.delta_x = delta_x;
    geom.delta_y = delta_y;
    geom.Coor_Ele_init = Coor_Ele_init;
    geom.Coor_Ele_init_r = Coor_Ele_init_r;
    geom.Coor_Ele_init_c = Coor_Ele_init_c;
    geom.Coor_Ele_init_x_d = Coor_Ele_init_x_d;
    geom.Coor_Ele_init_x_u = Coor_Ele_init_x_u;
    geom.Coor_Ele_init_y_l = Coor_Ele_init_y_l;
    geom.Coor_Ele_init_y_r = Coor_Ele_init_y_r;
    geom.M_x_r = M_x_r;
    geom.M_y_c = M_y_c;
    geom.d = d;
    geom.V_F0 = V_F0;
    geom.V_F = V_F;
end

function [Rot_fixed, normal_fixed, Coor_Ele_fixed] = build_rot_and_elements_fixed(Rot0, Coor_Ele_init)
    B = numel(Rot0);
    Rot_fixed = Rot0;
    normal_fixed = cell(B, 1);
    Coor_Ele_fixed = cell(B, 1);
    for b = 1:B
        normal_fixed{b} = Rot_fixed{b} * [0; 0; 1];
        Coor_Ele_fixed{b} = (Rot_fixed{b} * Coor_Ele_init.').';
    end
end

function [Rot_align, normal_align, Coor_Ele_align] = build_rot_and_elements_align(est, Coor_Ele_init)
    K = numel(est);
    Rot_align = cell(K, 1);
    normal_align = cell(K, 1);
    Coor_Ele_align = cell(K, 1);

    n1 = [0.02; 0.01; sqrt(1 - 0.02^2 - 0.01^2)];
    n1 = n1 / norm(n1);

    for k = 1:K
        n2 = est{k} / norm(est{k});
        u = cross(n1, n2);
        if norm(u) < 1e-12
            Rot_align{k} = eye(3);
        else
            u = u / norm(u);
            alph = acos(max(-1, min(1, n1.' * n2)));
            U = [0, -u(3), u(2); ...
                 u(3), 0, -u(1); ...
                 -u(2), u(1), 0];
            Rot_align{k} = cos(alph) * eye(3) + (1 - cos(alph)) * (u * u.') + sin(alph) * U;
        end
        normal_align{k} = Rot_align{k} * [0; 0; 1];
        Coor_Ele_align{k} = (Rot_align{k} * Coor_Ele_init.').';
    end
end

function mm = build_mm_for_fixed_positions(s, Rot_Matrix, normal_vector_rot, est, geom, lambda)
    % Build holographic beamforming weights mm for fixed positions (s) and given rotation matrices.
    B1 = size(s, 1);
    B = size(s, 2);

    [X, Y, Z] = Orientation_uniformSpherePoints(B1);
    q = [X.', Y.', Z.'];

    M = size(geom.Coor_Ele_init, 1);
    mm = zeros(M, B);

    for b = 1:B
        pos = s(:, b).' * q; % 1x3
        a1 = pos + (Rot_Matrix{b} * geom.Coor_Ele_init.').';

        if normal_vector_rot{b}.' * est{b} > 0
            a_st = exp(1j * 2 * pi / lambda * (a1 * est{b}));
        else
            a_st = zeros(M, 1);
        end

        mm(:, b) = real(conj(geom.V_F0) .* conj(a_st) + 1) / 2;
    end
end
