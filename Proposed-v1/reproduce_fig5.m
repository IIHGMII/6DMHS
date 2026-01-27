%% reproduce_fig5.m
% Reproduce paper Fig.5: EH power vs Transmit power (4 curves).
% Note: This script is written for the codebase in `Proposed-v1/`.
%
% 中文说明：
% - 运行位置：在 `Proposed-v1` 目录下运行本脚本。
% - 依赖：CVX（建议配合 MOSEK），以及本目录已有函数：
%   `Orientation_uniformSpherePoints.m`、`Channel_Generation_init.m`、`sig.m`、`fun_run_R00.m`
% - 为了让脚本可在合理时间内跑完，这里对“平移优化”采用了简化的离散穷举：
%   仅根据单用户全息增益挑选平移点，再把等效CSI输入到 `fun_run_R00` 做下行数字波束与功分优化。
% - 若你需要严格复现论文(P1)的交替优化，可把 `select_translation_positions(...)` 替换为你完整的(P1)实现。

clc; clear; close all;
rng(1, 'twister');

%% ======================= 可调参数 =======================
Pt_dBm_list = 30:45;      % 发射功率(dBm)

% 论文仿真参数：K_R = 10（线性值，对应 10 dB）
K_rice_dB   = 10;

R00         = 1;          % 吞吐阈值(bit/s/Hz)

% 论文仿真参数：sigma0^2=-100 dBm, sigma_cov^2=-50 dBm
sigma0_dBm     = -100;    % 天线噪声
sigma_cov_dBm  = -50;     % 转换噪声
sigma0      = 10^((sigma0_dBm - 30) / 10);
sigma_cov   = 10^((sigma_cov_dBm - 30) / 10);

K = 4;                    % 用户数（论文中K=B）
B = 4;                    % RHS数

B1   = 50;                % 平移候选点数（论文：M1=50）
dmin = 0.1;               % 平移点最小距离约束（单位与q一致，当前q在单位球面上）

% RHS最大增益方向（“未旋转”的局部坐标系中）
% 论文在“Directional Beamforming Gain”示例中给出的最大增益方向为 [0.175, -0.275, 0.945]^T
u_max = [0.175; -0.275; 0.945];
u_max = u_max / norm(u_max, 2);

% FP迭代次数（fun_run_R00 外层迭代）；越大越慢
fp_opts = struct('max_iter', 30, 'quiet', true, 'do_plot', false, 'sigma_cov', sigma_cov);

% 平移选点(P1)优化的参考发射功率（按论文参数 Pt=40 dBm）
Pt_ref_dBm = 40;
Ptx_ref = 10^((Pt_ref_dBm - 30)/10);

% 平移选点优化配置：越大越接近论文的交替优化，但也越慢
pos_opts = struct('outer_iter', 3, 'greedy_iter', 2, 'fp_iter', 10, 'quiet', true);

%% ======================= 固定系统参数 =======================
fc = 30e9;
c  = 3e8;
lambda = c / fc;
dx = lambda / 4;
dy = lambda / 4;

Mx = 32; My = 32; M = Mx * My;

delta_x = (2 * (0:(Mx-1)).' - Mx + 1) / 2;
delta_y = (2 * (0:(My-1)).' - My + 1) / 2;
Coor_Ele_rel = [kron(delta_x * dx, ones(My, 1)), kron(ones(Mx, 1), delta_y * dy), zeros(M, 1)];

% 馈源模型（与原始脚本保持一致）
Coor_Feed = [0.002984305671304, 0.001586131606604, 0];
Dis_Feed2Ele = sqrt((Coor_Ele_rel - Coor_Feed).^2 * ones(3, 1));
V_F0 = exp(-1j * 2*pi*sqrt(3)/lambda * Dis_Feed2Ele);
eta0 = 8/3/M;
V_F  = sqrt(eta0) * V_F0;

K_rice = 10^(K_rice_dB/10);

%% ======================= 生成“环境参数”（方向/散射/损耗）=======================
% 这里借用 Channel_Generation_init 返回的 e0/e2/eta/Omega0（文件内已硬编码一组场景）
% 注意：本脚本后续会基于不同的(q,R)重新计算下行信道相位，但方向/散射系数沿用该场景。
[q_init, R_init] = init_positions_and_rotations(B);
[Coor_Ele_init_cell, normal_init] = build_global_geometry(q_init, R_init, Coor_Ele_rel);
[Coor_Ele_r, Coor_Ele_c, Coor_Ele_x_d, Coor_Ele_x_u, Coor_Ele_y_l, Coor_Ele_y_r, M_x_r, M_y_c] = ...
    build_sensing_geometry(q_init, R_init, Mx, My, dx, dy, lambda, 1);

[~, e0, e2, ~, ~, ~, ~, Iota, eta, ~, ~, ~, ~, ~, ~, ~, ~, Omega0] = ...
    Channel_Generation_init(K, B, M, Mx, My, M_x_r, M_y_c, normal_init, K_rice, ...
        Coor_Ele_init_cell, Coor_Ele_r, Coor_Ele_c, lambda, R_init, ...
        Coor_Ele_x_d, Coor_Ele_x_u, Coor_Ele_y_l, Coor_Ele_y_r);

% Fig.5 一般用于展示“不同发射机配置”的收益，默认按“无感知误差”处理
est = e0;

%% ======================= 预先构造4种配置(q,R) =======================
% FPA：不旋转、不平移
cfg.fpa.q = q_init;
cfg.fpa.R = R_init;

% Rotation only：旋转对齐最大增益方向，但平移固定
cfg.rot_only.q = q_init;
cfg.rot_only.R = cell(B, 1);
for b = 1:B
    cfg.rot_only.R{b} = rot_from_to(u_max, est{b});
end

% Translation only：平移优化，但不旋转（保持初始姿态）
[q_cand, ~] = init_positions_and_rotations(B1);
[cfg.trans_only.q, cfg.trans_only.R] = optimize_translation_positions_trans_only_p1(q_cand, est, Coor_Ele_rel, V_F0, V_F, lambda, dmin, Ptx_ref, pos_opts);

% Proposed 6DMHS：旋转 + 平移优化
cfg.proposed.R = cfg.rot_only.R;
cfg.proposed.q = optimize_translation_positions_p1(q_cand, cfg.proposed.R, est, Coor_Ele_rel, V_F0, V_F, lambda, dmin, Ptx_ref, pos_opts);

cfg_order = {'proposed', 'rot_only', 'trans_only', 'fpa'};
cfg_label = {'Proposed 6DMHS', '6DMHS with rotation only', '6DMHS with translation only', 'FPA'};
cfg_style = {'-o', '-v', '-d', '-s'};

%% ======================= 扫 Pt，计算 EH =======================
eh_min = zeros(numel(Pt_dBm_list), numel(cfg_order));

for ip = 1:numel(Pt_dBm_list)
    Pt_dBm = Pt_dBm_list(ip);
    Ptx = 10^((Pt_dBm - 30)/10);

    fprintf('[%02d/%02d] Pt=%g dBm (Ptx=%g W)\n', ip, numel(Pt_dBm_list), Pt_dBm, Ptx);

    for ic = 1:numel(cfg_order)
        key = cfg_order{ic};
        q_use = cfg.(key).q;
        R_use = cfg.(key).R;

        [Coor_Ele_cell, normal_vec] = build_global_geometry(q_use, R_use, Coor_Ele_rel);
        mm = build_holographic_beamformer(Coor_Ele_cell, est, V_F0, lambda);
        H_eff = build_effective_channel(Coor_Ele_cell, normal_vec, mm, V_F, e0, e2, eta, Omega0, Iota, K_rice, lambda);

        best_val = -inf;
        max_restart = 3;
        for rr = 1:max_restart
            rho0  = 0.1 * rand(K, 1);
            psi10 = 0.05 * (randn(K, K) + 1j * randn(K, K));

            % 先通过“吞吐率子问题”迭代得到一个更稳健的psi2初值（避免FP下界为负导致CVX不可行）
            psi2_fixed = init_psi2_fixed(H_eff, rho0, sigma0, sigma_cov, Ptx, B, K);

            val = fun_run_R00(R00, H_eff, rho0, psi10, psi2_fixed, sigma0, Ptx, B, K, fp_opts);
            if isfinite(val) && val > best_val
                best_val = val;
            end
        end
        eh_min(ip, ic) = best_val;
    end
end

%% ======================= 画图 =======================
figure('Color', 'w');
hold on;
for ic = 1:numel(cfg_order)
    plot(Pt_dBm_list, eh_min(:, ic), cfg_style{ic}, 'LineWidth', 2);
end
grid on;
xlabel('Transmit power (dBm)');
ylabel('EH power (W)');
legend(cfg_label, 'Location', 'northwest');
title('Fig.5：EH power vs Transmit power');

%% ======================= 保存结果 =======================
out_dir = fullfile(pwd, 'Data_record');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save(fullfile(out_dir, 'fig5_results.mat'), 'Pt_dBm_list', 'eh_min', 'cfg_order', 'cfg_label', 'K_rice_dB', 'R00', 'sigma0');
saveas(gcf, fullfile(out_dir, 'fig5_reproduce.png'));

fprintf('已保存：%s\n', fullfile(out_dir, 'fig5_results.mat'));
fprintf('已保存：%s\n', fullfile(out_dir, 'fig5_reproduce.png'));

%% ======================= 本脚本使用的局部函数 =======================

function [q, R] = init_positions_and_rotations(N)
    [X, Y, Z] = Orientation_uniformSpherePoints(N);
    q = [X(:), Y(:), Z(:)];
    R = cell(N, 1);
    n1 = [0; 0; 1];
    for i = 1:N
        n2 = q(i, :).';
        n2 = n2 / norm(n2, 2);
        R{i} = rot_from_to(n1, n2);
    end
end

function R = rot_from_to(a, b)
    a = a / max(norm(a, 2), eps);
    b = b / max(norm(b, 2), eps);
    c = dot(a, b);
    v = cross(a, b);
    s = norm(v, 2);
    if s < 1e-12
        if c >= 0
            R = eye(3);
            return;
        end
        % 180°：选一个与a不共线的轴
        tmp = [1; 0; 0];
        if abs(dot(a, tmp)) > 0.9
            tmp = [0; 1; 0];
        end
        axis = cross(a, tmp);
        axis = axis / max(norm(axis, 2), eps);
        A = skew(axis);
        R = eye(3) + 2 * (A * A);
        return;
    end
    axis = v / s;
    theta = atan2(s, c);
    A = skew(axis);
    R = eye(3) * cos(theta) + (1 - cos(theta)) * (axis * axis.') + sin(theta) * A;
end

function A = skew(u)
    A = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
end

function [Coor_Ele_cell, normal_vec] = build_global_geometry(q, R, Coor_Ele_rel)
    B = size(q, 1);
    Coor_Ele_cell = cell(B, 1);
    normal_vec = cell(B, 1);
    for b = 1:B
        Coor_Ele_cell{b} = q(b, :) + Coor_Ele_rel * R{b}.';
        normal_vec{b} = R{b} * [0; 0; 1];
    end
end

function mm = build_holographic_beamformer(Coor_Ele_cell, est, V_F0, lambda)
    % Q=1 情况下：Psi = (Re{conj(V_F0).*conj(a)} + 1)/2
    B = numel(Coor_Ele_cell);
    M = size(V_F0, 1);
    mm = zeros(M, B);
    for b = 1:B
        a = exp(1j * 2*pi/lambda * (Coor_Ele_cell{b} * est{b}));
        mm(:, b) = (real(conj(V_F0) .* conj(a)) + 1) / 2;
    end
end

function H_eff = build_effective_channel(Coor_Ele_cell, normal_vec, mm, V_F, e0, e2, eta, Omega0, Iota, K_rice, lambda)
    % H_eff(b,k) = h_{k,b}^T diag(mm_b) V_F
    B = numel(Coor_Ele_cell);
    K = numel(e0);
    M = size(mm, 1);
    H_eff = zeros(B, K);
    for b = 1:B
        v_b = mm(:, b) .* V_F; % diag(mm_b)*V_F
        for k = 1:K
            h_kb = zeros(M, 1);

            % LoS
            h_kb = h_kb + sig(normal_vec{b}.' * e0{k}) ...
                * sqrt(K_rice/(1 + K_rice)) * sqrt(Omega0(k, b)) ...
                * exp(1j * 2*pi/lambda * (Coor_Ele_cell{b} * e0{k}));

            % NLoS
            for iota = 1:Iota
                dir = e2{k, iota};
                h_kb = h_kb + sig(normal_vec{b}.' * dir) ...
                    * sqrt(1/(1 + K_rice)) * sqrt(Omega0(k, b)) * eta(k, iota) ...
                    * exp(1j * 2*pi/lambda * (Coor_Ele_cell{b} * dir));
            end

            H_eff(b, k) = h_kb.' * v_b;
        end
    end
end

function q_best = select_translation_positions(q_cand, R_use, est, Coor_Ele_rel, V_F0, V_F, lambda, dmin)
    % 简化的平移选择：对每个RHS，先计算在每个候选点上的“单用户全息增益”，再做穷举最大化最小值
    B = numel(R_use);
    if B ~= 4
        error('当前实现为便于复现，仅支持 B=4（你可自行扩展）。');
    end

    B1 = size(q_cand, 1);
    g_tbl = zeros(B, B1);
    n_tbl = zeros(3, B);
    for b = 1:B
        n_tbl(:, b) = R_use{b} * [0; 0; 1];
    end

    for b = 1:B
        for i = 1:B1
            Coor_tmp = q_cand(i, :) + Coor_Ele_rel * R_use{b}.';
            a = exp(1j * 2*pi/lambda * (Coor_tmp * est{b}));
            mm = (real(conj(V_F0) .* conj(a)) + 1) / 2;
            g_tbl(b, i) = abs(a.' * (mm .* V_F))^2;
        end
    end

    % 候选集合：与用户方向同半球
    C_set = cell(B, 1);
    for b = 1:B
        % 论文约束：n_b^T q_b >= 0（避免遮挡）；并要求在用户方向半球
        C_set{b} = find((q_cand * est{b}) > 0 & (q_cand * n_tbl(:, b)) >= 0);
        if isempty(C_set{b})
            error('平移候选集合为空：请增大B1或检查est方向。');
        end
    end

    % 预计算距离矩阵
    D = zeros(B1, B1);
    for i = 1:B1
        for j = 1:B1
            D(i, j) = norm(q_cand(i, :) - q_cand(j, :), 2);
        end
    end

    best_val = -inf;
    best_idx = [C_set{1}(1), C_set{2}(1), C_set{3}(1), C_set{4}(1)];

    for i1 = C_set{1}.'
        for i2 = C_set{2}.'
            if i2 == i1, continue; end
            if D(i1, i2) < dmin, continue; end
            % 论文反射约束(近似)：n_1^T q_2 < n_1^T q_1，n_2^T q_1 < n_2^T q_2
            if (n_tbl(:, 1).'*q_cand(i2, :).') >= (n_tbl(:, 1).'*q_cand(i1, :).'), continue; end
            if (n_tbl(:, 2).'*q_cand(i1, :).') >= (n_tbl(:, 2).'*q_cand(i2, :).'), continue; end
            for i3 = C_set{3}.'
                if any(i3 == [i1, i2]), continue; end
                if D(i1, i3) < dmin || D(i2, i3) < dmin, continue; end
                if (n_tbl(:, 1).'*q_cand(i3, :).') >= (n_tbl(:, 1).'*q_cand(i1, :).'), continue; end
                if (n_tbl(:, 3).'*q_cand(i1, :).') >= (n_tbl(:, 3).'*q_cand(i3, :).'), continue; end
                if (n_tbl(:, 2).'*q_cand(i3, :).') >= (n_tbl(:, 2).'*q_cand(i2, :).'), continue; end
                if (n_tbl(:, 3).'*q_cand(i2, :).') >= (n_tbl(:, 3).'*q_cand(i3, :).'), continue; end
                for i4 = C_set{4}.'
                    if any(i4 == [i1, i2, i3]), continue; end
                    if D(i1, i4) < dmin || D(i2, i4) < dmin || D(i3, i4) < dmin, continue; end
                    if (n_tbl(:, 1).'*q_cand(i4, :).') >= (n_tbl(:, 1).'*q_cand(i1, :).'), continue; end
                    if (n_tbl(:, 4).'*q_cand(i1, :).') >= (n_tbl(:, 4).'*q_cand(i4, :).'), continue; end
                    if (n_tbl(:, 2).'*q_cand(i4, :).') >= (n_tbl(:, 2).'*q_cand(i2, :).'), continue; end
                    if (n_tbl(:, 4).'*q_cand(i2, :).') >= (n_tbl(:, 4).'*q_cand(i4, :).'), continue; end
                    if (n_tbl(:, 3).'*q_cand(i4, :).') >= (n_tbl(:, 3).'*q_cand(i3, :).'), continue; end
                    if (n_tbl(:, 4).'*q_cand(i3, :).') >= (n_tbl(:, 4).'*q_cand(i4, :).'), continue; end

                    val = min([g_tbl(1, i1), g_tbl(2, i2), g_tbl(3, i3), g_tbl(4, i4)]);
                    if val > best_val
                        best_val = val;
                        best_idx = [i1, i2, i3, i4];
                    end
                end
            end
        end
    end

    q_best = q_cand(best_idx, :);
end

function [q_best, R_best] = select_translation_positions_trans_only(q_cand, est, Coor_Ele_rel, V_F0, V_F, lambda, dmin)
    % Translation-only benchmark:
    % - RHS positions are optimized on a radius-1 sphere
    % - For each RHS, its orientation is fixed to the normal direction from BS to its center position

    B = numel(est);
    if B ~= 4
        error('当前实现为便于复现，仅支持 B=4（你可自行扩展）。');
    end

    B1 = size(q_cand, 1);
    R_cand = cell(B1, 1);
    for i = 1:B1
        R_cand{i} = rot_from_to([0; 0; 1], q_cand(i, :).');
    end

    g_tbl = zeros(B, B1);
    for b = 1:B
        for i = 1:B1
            Coor_tmp = q_cand(i, :) + Coor_Ele_rel * R_cand{i}.';
            a = exp(1j * 2*pi/lambda * (Coor_tmp * est{b}));
            mm = (real(conj(V_F0) .* conj(a)) + 1) / 2;
            g_tbl(b, i) = abs(a.' * (mm .* V_F))^2;
        end
    end

    % 候选集合：与用户方向同半球（论文中由可行性约束给出）
    C_set = cell(B, 1);
    for b = 1:B
        C_set{b} = find((q_cand * est{b}) > 0);
        if isempty(C_set{b})
            error('平移候选集合为空：请增大B1或检查est方向。');
        end
    end

    % 预计算距离矩阵
    D = zeros(B1, B1);
    for i = 1:B1
        for j = 1:B1
            D(i, j) = norm(q_cand(i, :) - q_cand(j, :), 2);
        end
    end

    best_val = -inf;
    best_idx = [C_set{1}(1), C_set{2}(1), C_set{3}(1), C_set{4}(1)];

    for i1 = C_set{1}.'
        for i2 = C_set{2}.'
            if i2 == i1, continue; end
            if D(i1, i2) < dmin, continue; end
            for i3 = C_set{3}.'
                if any(i3 == [i1, i2]), continue; end
                if D(i1, i3) < dmin || D(i2, i3) < dmin, continue; end
                for i4 = C_set{4}.'
                    if any(i4 == [i1, i2, i3]), continue; end
                    if D(i1, i4) < dmin || D(i2, i4) < dmin || D(i3, i4) < dmin, continue; end
                    val = min([g_tbl(1, i1), g_tbl(2, i2), g_tbl(3, i3), g_tbl(4, i4)]);
                    if val > best_val
                        best_val = val;
                        best_idx = [i1, i2, i3, i4];
                    end
                end
            end
        end
    end

    q_best = q_cand(best_idx, :);
    R_best = cell(B, 1);
    for b = 1:B
        R_best{b} = R_cand{best_idx(b)};
    end
end

function q_best = optimize_translation_positions_p1(q_cand, R_use, est, Coor_Ele_rel, V_F0, V_F, lambda, dmin, Ptx_ref, opts)
    % Proposed 平移选点：固定旋转R_use，交替优化位置与数字波束（对应论文(P1)的思想）

    if nargin < 9 || isempty(Ptx_ref)
        Ptx_ref = 10;
    end
    if nargin < 10 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'outer_iter'),  opts.outer_iter  = 3; end
    if ~isfield(opts, 'greedy_iter'), opts.greedy_iter = 2; end
    if ~isfield(opts, 'fp_iter'),     opts.fp_iter     = 10; end
    if ~isfield(opts, 'quiet'),       opts.quiet       = true; end

    B = numel(R_use);
    K = numel(est);
    if B ~= 4 || K ~= 4
        error('当前实现为便于复现，仅支持 B=4, K=4（你可自行扩展）。');
    end

    B1 = size(q_cand, 1);

    % 固定法向量
    n_tbl = zeros(3, B);
    for b = 1:B
        n_tbl(:, b) = R_use{b} * [0; 0; 1];
    end

    % 候选集合：与对应用户同半球 + 避免阻塞（n_b^T q_b >= 0）
    C_set = cell(B, 1);
    for b = 1:B
        C_set{b} = find((q_cand * est{b}) > 0 & (q_cand * n_tbl(:, b)) >= 0);
        if isempty(C_set{b})
            error('平移候选集合为空：请增大B1或检查est方向。');
        end
    end

    % 距离矩阵
    D = zeros(B1, B1);
    for i = 1:B1
        for j = 1:B1
            D(i, j) = norm(q_cand(i, :) - q_cand(j, :), 2);
        end
    end

    % 预计算 g_{k,b}（每个b、每个候选点、每个用户方向）
    H_gain = zeros(B, B1, K);
    for b = 1:B
        n_b = n_tbl(:, b);
        for i = 1:B1
            Coor_tmp = q_cand(i, :) + Coor_Ele_rel * R_use{b}.';

            % Psi_b：按用户b方向设计（与论文一致）
            if (n_b.' * est{b}) <= 0
                mm_b = zeros(size(V_F0));
            else
                a_bb = exp(1j * 2*pi/lambda * (Coor_tmp * est{b}));
                mm_b = (real(conj(V_F0) .* conj(a_bb)) + 1) / 2;
            end
            v_b = mm_b .* V_F;

            for k = 1:K
                if (n_b.' * est{k}) <= 0
                    H_gain(b, i, k) = 0;
                else
                    a_bk = exp(1j * 2*pi/lambda * (Coor_tmp * est{k}));
                    H_gain(b, i, k) = a_bk.' * v_b;
                end
            end
        end
    end

    % 初始化：找一个可行组合
    idx_sel = find_initial_positions_p1(C_set, q_cand, D, n_tbl, dmin);

    % 初始化数字波束
    Xtr = randn(B, K) + 1j * randn(B, K);
    Xtr = Xtr / max(norm(Xtr, 'fro'), eps) * sqrt(Ptx_ref);

    % 交替优化：位置(离散)+数字波束(连续)
    for outer = 1:opts.outer_iter
        for gg = 1:opts.greedy_iter
            for b = 1:B
                best_idx = idx_sel(b);
                best_val = eval_min_gain_from_tbl(H_gain, idx_sel, Xtr);

                for cand = C_set{b}.'
                    if cand == idx_sel(b)
                        continue;
                    end
                    idx_try = idx_sel;
                    idx_try(b) = cand;
                    if ~is_feasible_positions_p1(idx_try, q_cand, D, n_tbl, dmin)
                        continue;
                    end
                    val = eval_min_gain_from_tbl(H_gain, idx_try, Xtr);
                    if val > best_val
                        best_val = val;
                        best_idx = cand;
                    end
                end
                idx_sel(b) = best_idx;
            end
        end

        H_sel = gather_gain_matrix(H_gain, idx_sel);
        [Xtr_new, ~] = solve_gain_fp(H_sel, Ptx_ref, opts.fp_iter, opts.quiet);
        if ~isempty(Xtr_new)
            Xtr = Xtr_new;
        end
    end

    q_best = q_cand(idx_sel, :);
end

function [q_best, R_best] = optimize_translation_positions_trans_only_p1(q_cand, est, Coor_Ele_rel, V_F0, V_F, lambda, dmin, Ptx_ref, opts)
    % Translation-only 平移选点：R与q绑定（法向对齐BS->中心），交替优化位置与数字波束

    if nargin < 8 || isempty(Ptx_ref)
        Ptx_ref = 10;
    end
    if nargin < 9 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'outer_iter'),  opts.outer_iter  = 3; end
    if ~isfield(opts, 'greedy_iter'), opts.greedy_iter = 2; end
    if ~isfield(opts, 'fp_iter'),     opts.fp_iter     = 10; end
    if ~isfield(opts, 'quiet'),       opts.quiet       = true; end

    B = numel(est);
    K = numel(est);
    if B ~= 4 || K ~= 4
        error('当前实现为便于复现，仅支持 B=4, K=4（你可自行扩展）。');
    end

    B1 = size(q_cand, 1);

    % 每个候选点对应一个旋转（法向对齐BS->中心）
    R_cand = cell(B1, 1);
    n_cand = zeros(B1, 3);
    for i = 1:B1
        R_cand{i} = rot_from_to([0; 0; 1], q_cand(i, :).');
        n_cand(i, :) = (R_cand{i} * [0; 0; 1]).';
    end

    % 候选集合：与对应用户同半球
    C_set = cell(B, 1);
    for b = 1:B
        C_set{b} = find((q_cand * est{b}) > 0);
        if isempty(C_set{b})
            error('平移候选集合为空：请增大B1或检查est方向。');
        end
    end

    % 距离矩阵
    D = zeros(B1, B1);
    for i = 1:B1
        for j = 1:B1
            D(i, j) = norm(q_cand(i, :) - q_cand(j, :), 2);
        end
    end

    % 预计算 g_{k,b}
    H_gain = zeros(B, B1, K);
    for b = 1:B
        for i = 1:B1
            R_i = R_cand{i};
            n_i = n_cand(i, :).';
            Coor_tmp = q_cand(i, :) + Coor_Ele_rel * R_i.';

            % Psi_b：按用户b方向设计
            if (n_i.' * est{b}) <= 0
                mm_b = zeros(size(V_F0));
            else
                a_bb = exp(1j * 2*pi/lambda * (Coor_tmp * est{b}));
                mm_b = (real(conj(V_F0) .* conj(a_bb)) + 1) / 2;
            end
            v_b = mm_b .* V_F;

            for k = 1:K
                if (n_i.' * est{k}) <= 0
                    H_gain(b, i, k) = 0;
                else
                    a_bk = exp(1j * 2*pi/lambda * (Coor_tmp * est{k}));
                    H_gain(b, i, k) = a_bk.' * v_b;
                end
            end
        end
    end

    % 初始化：找一个可行组合（translation-only主要是dmin）
    idx_sel = find_initial_positions_trans_only(C_set, D, dmin);

    % 初始化数字波束
    Xtr = randn(B, K) + 1j * randn(B, K);
    Xtr = Xtr / max(norm(Xtr, 'fro'), eps) * sqrt(Ptx_ref);

    for outer = 1:opts.outer_iter
        for gg = 1:opts.greedy_iter
            for b = 1:B
                best_idx = idx_sel(b);
                best_val = eval_min_gain_from_tbl(H_gain, idx_sel, Xtr);

                for cand = C_set{b}.'
                    if cand == idx_sel(b)
                        continue;
                    end
                    idx_try = idx_sel;
                    idx_try(b) = cand;
                    if ~is_feasible_positions_trans_only(idx_try, D, dmin)
                        continue;
                    end
                    val = eval_min_gain_from_tbl(H_gain, idx_try, Xtr);
                    if val > best_val
                        best_val = val;
                        best_idx = cand;
                    end
                end
                idx_sel(b) = best_idx;
            end
        end

        H_sel = gather_gain_matrix(H_gain, idx_sel);
        [Xtr_new, ~] = solve_gain_fp(H_sel, Ptx_ref, opts.fp_iter, opts.quiet);
        if ~isempty(Xtr_new)
            Xtr = Xtr_new;
        end
    end

    q_best = q_cand(idx_sel, :);
    R_best = cell(B, 1);
    for b = 1:B
        R_best{b} = R_cand{idx_sel(b)};
    end
end

function idx_sel = find_initial_positions_p1(C_set, q_cand, D, n_tbl, dmin)
    B = numel(C_set);
    if B ~= 4
        error('当前实现仅支持B=4。');
    end

    for i1 = C_set{1}.'
        for i2 = C_set{2}.'
            if any(i2 == i1), continue; end
            if D(i1, i2) < dmin, continue; end
            if (n_tbl(:, 1).'* (q_cand(i2, :).'-q_cand(i1, :).')) >= -1e-12, continue; end
            if (n_tbl(:, 2).'* (q_cand(i1, :).'-q_cand(i2, :).')) >= -1e-12, continue; end

            for i3 = C_set{3}.'
                if any(i3 == [i1, i2]), continue; end
                if D(i1, i3) < dmin || D(i2, i3) < dmin, continue; end
                for i4 = C_set{4}.'
                    idx_try = [i1, i2, i3, i4];
                    if is_feasible_positions_p1(idx_try, q_cand, D, n_tbl, dmin)
                        idx_sel = idx_try;
                        return;
                    end
                end
            end
        end
    end

    error('未找到满足约束的初始可行平移点组合，请增大B1或放宽dmin。');
end

function idx_sel = find_initial_positions_trans_only(C_set, D, dmin)
    B = numel(C_set);
    if B ~= 4
        error('当前实现仅支持B=4。');
    end

    for i1 = C_set{1}.'
        for i2 = C_set{2}.'
            if any(i2 == i1), continue; end
            if D(i1, i2) < dmin, continue; end
            for i3 = C_set{3}.'
                if any(i3 == [i1, i2]), continue; end
                if D(i1, i3) < dmin || D(i2, i3) < dmin, continue; end
                for i4 = C_set{4}.'
                    idx_try = [i1, i2, i3, i4];
                    if is_feasible_positions_trans_only(idx_try, D, dmin)
                        idx_sel = idx_try;
                        return;
                    end
                end
            end
        end
    end

    error('未找到满足约束的初始可行平移点组合，请增大B1或放宽dmin。');
end

function ok = is_feasible_positions_p1(idx_sel, q_cand, D, n_tbl, dmin)
    B = numel(idx_sel);
    if numel(unique(idx_sel)) < B
        ok = false;
        return;
    end

    % 距离约束
    for b1 = 1:B
        for b2 = b1+1:B
            if D(idx_sel(b1), idx_sel(b2)) < dmin
                ok = false;
                return;
            end
        end
    end

    % 阻塞约束：n_b^T q_b >= 0
    for b = 1:B
        if (n_tbl(:, b).'*q_cand(idx_sel(b), :).') < 0
            ok = false;
            return;
        end
    end

    % 反射约束（按原代码库符号）：n_b^T(q_{b2}-q_b) < 0, ∀b2≠b
    for b1 = 1:B
        q1 = q_cand(idx_sel(b1), :).';
        for b2 = 1:B
            if b2 == b1, continue; end
            q2 = q_cand(idx_sel(b2), :).';
            if (n_tbl(:, b1).'* (q2 - q1)) >= -1e-12
                ok = false;
                return;
            end
        end
    end

    ok = true;
end

function ok = is_feasible_positions_trans_only(idx_sel, D, dmin)
    B = numel(idx_sel);
    if numel(unique(idx_sel)) < B
        ok = false;
        return;
    end
    for b1 = 1:B
        for b2 = b1+1:B
            if D(idx_sel(b1), idx_sel(b2)) < dmin
                ok = false;
                return;
            end
        end
    end
    ok = true;
end

function H_sel = gather_gain_matrix(H_gain, idx_sel)
    [B, ~, K] = size(H_gain);
    H_sel = zeros(B, K);
    for b = 1:B
        H_sel(b, :) = reshape(H_gain(b, idx_sel(b), :), 1, K);
    end
end

function val = eval_min_gain_from_tbl(H_gain, idx_sel, Xtr)
    H_sel = gather_gain_matrix(H_gain, idx_sel);
    K = size(H_sel, 2);
    val_k = zeros(K, 1);
    for k = 1:K
        y = H_sel(:, k).' * Xtr; % 1xK
        val_k(k) = sum(abs(y).^2);
    end
    val = min(val_k);
end

function [Xtr_best, best_val] = solve_gain_fp(H_sel, Ptx, max_iter, quiet)
    % 给定方向增益通道H_sel，优化数字波束Xtr使 min_k ||H(:,k)^T Xtr||^2 最大
    B = size(H_sel, 1);
    K = size(H_sel, 2);

    psi = 0.05 * (randn(K, K) + 1j * randn(K, K));
    Xtr_best = [];
    best_val = -inf;

    for it = 1:max_iter
        if quiet
            cvx_begin quiet
        else
            cvx_begin
        end
        cvx_solver mosek
        variable Xtr(B, K) complex
        variable P0(1)
        expressions P_temp(K, K) G_lb(K, 1)

        for k = 1:K
            for k1 = 1:K
                P_temp(k1, k) = H_sel(:, k).' * Xtr(:, k1);
            end
            G_lb(k) = 2 * real(psi(:, k)' * P_temp(:, k)) - abs(psi(:, k)' * psi(:, k));
        end

        maximize P0
        subject to
        G_lb >= P0;
        sum(sum(pow_abs(Xtr, 2))) <= Ptx;
        cvx_end

        if ~(strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved'))
            break;
        end

        for k = 1:K
            psi(:, k) = P_temp(:, k);
        end

        if cvx_optval > best_val
            best_val = cvx_optval;
            Xtr_best = Xtr;
        end
    end
end

function [Coor_Ele_r, Coor_Ele_c, Coor_Ele_x_d, Coor_Ele_x_u, Coor_Ele_y_l, Coor_Ele_y_r, M_x_r, M_y_c] = ...
    build_sensing_geometry(q, R, Mx, My, dx, dy, lambda, kappa0)
    % 仅用于喂给 Channel_Generation_init，以获取 e0/e2/eta/Omega0 等场景参数
    delta_x = (2 * (0:(Mx-1)).' - Mx + 1) / 2;
    delta_y = (2 * (0:(My-1)).' - My + 1) / 2;

    Coor_Ele_init_r = [delta_x * dx, (-Mx + 1)/2 * (lambda/2) * ones(Mx, 1), zeros(Mx, 1)];
    Coor_Ele_init_c = [(Mx - 1)/2 * dx * ones(My, 1), delta_y * dy, zeros(My, 1)];

    delta_ele = 2^kappa0;
    Coor_Ele_init_r = Coor_Ele_init_r(1:delta_ele:end, :);
    Coor_Ele_init_c = Coor_Ele_init_c(1:delta_ele:end, :);

    M_x_r = size(Coor_Ele_init_r, 1);
    M_y_c = size(Coor_Ele_init_c, 1);

    delta_x_r = [(-M_x_r + 1)/2 : 1 : (M_x_r - 1)/2].';
    delta_y_c = [(-M_y_c + 1)/2 : 1 : (M_y_c - 1)/2].';

    Coor_Ele_init_x_d = [delta_x_r * (lambda/2), (-M_y_c + 1)/2 * (lambda/2) * ones(M_x_r, 1), zeros(M_x_r, 1)];
    Coor_Ele_init_x_u = [delta_x_r * (lambda/2), ( M_y_c - 1)/2 * (lambda/2) * ones(M_x_r, 1), zeros(M_x_r, 1)];
    Coor_Ele_init_y_l = [(-M_x_r + 1)/2 * (lambda/2) * ones(M_y_c, 1), delta_y_c * (lambda/2), zeros(M_y_c, 1)];
    Coor_Ele_init_y_r = [( M_x_r - 1)/2 * (lambda/2) * ones(M_y_c, 1), delta_y_c * (lambda/2), zeros(M_y_c, 1)];

    B = size(q, 1);
    Coor_Ele_r  = cell(B, 1);
    Coor_Ele_c  = cell(B, 1);
    Coor_Ele_x_d = cell(B, 1);
    Coor_Ele_x_u = cell(B, 1);
    Coor_Ele_y_l = cell(B, 1);
    Coor_Ele_y_r = cell(B, 1);
    for b = 1:B
        Coor_Ele_r{b}  = q(b, :) + Coor_Ele_init_r * R{b}.';
        Coor_Ele_c{b}  = q(b, :) + Coor_Ele_init_c * R{b}.';
        Coor_Ele_x_d{b} = q(b, :) + Coor_Ele_init_x_d * R{b}.';
        Coor_Ele_x_u{b} = q(b, :) + Coor_Ele_init_x_u * R{b}.';
        Coor_Ele_y_l{b} = q(b, :) + Coor_Ele_init_y_l * R{b}.';
        Coor_Ele_y_r{b} = q(b, :) + Coor_Ele_init_y_r * R{b}.';
    end
end

function psi2_fixed = init_psi2_fixed(H_eff, rho, sigma0, sigma_cov, Ptx, B, K)
    % 用与 RUN_OPT_Protocol 类似的“先最大化吞吐率”步骤，得到一个可行的psi2初值
    % 这样可以显著降低 fun_run_R00 第一次CVX求解就不可行的概率。

    psi2_fixed = 0.01 * (randn(K, 1) + 1j * randn(K, 1));

    max_restart = 3;
    max_iter = 15;

    for rr = 1:max_restart
        psi2 = psi2_fixed;
        ok = true;

        for it = 1:max_iter
            if any(isnan(psi2))
                ok = false;
                break;
            end

            cvx_begin quiet
            cvx_solver mosek
            variable Xtr(B, K) complex
            variable R0(1)
            expressions R_temp(K, K) R_thro(K, 1)

            for k  = 1:K
                for k1 = 1:K
                    if k ~= k1
                        R_temp(k1, k) = quad_form(H_eff(:, k).' * Xtr(:, k1), 1);
                    else
                        R_temp(k1, k) = 0;
                    end
                end
            end

            for k = 1:K
                % log2(1+·) = log(1+·)/log(2)
                noise_term = (1 - rho(k)) * sigma0 + sigma_cov;
                R_thro(k) = log( 1 + 2 * sqrt(1 - rho(k)) * real( psi2(k)' * Xtr(:, k)' * conj(H_eff(:, k)) ) ...
                    - abs(psi2(k))^2 * ( (1 - rho(k)) * sum(R_temp(:, k)) + noise_term ) ) / log(2);
            end

            maximize R0
            subject to
            R_thro >= R0;
            sum(sum(pow_abs(Xtr, 2))) <= Ptx;
            cvx_end

            if ~(strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved'))
                ok = false;
                break;
            end

            for k = 1:K
                denom = (1 - rho(k)) * sum(R_temp(:, k)) + (1 - rho(k)) * sigma0 + sigma_cov;
                psi2(k) = sqrt(1 - rho(k)) * Xtr(:, k)' * conj(H_eff(:, k)) / denom;
            end
        end

        if ok && ~any(isnan(psi2))
            psi2_fixed = psi2;
            return;
        end

        % 重新随机初始化，重试
        psi2_fixed = 0.01 * (randn(K, 1) + 1j * randn(K, 1));
    end
end
