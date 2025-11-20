%% 绘制 EH Power Performance versus Transmit Power
clc; clear; close all;

%% 参数设置
Pt_range = 30:2:50;  % 发射功率范围 (dBm)
num_trials = length(Pt_range);
K_rice_dB = -10;
K_rice = 10^(K_rice_dB/10);

P_EH = zeros(num_trials, 1);

%% 系统参数初始化
kappa0 = 1;
N = 256;
sigma0 = 1e-12;
R00 = 1;

fc = 30e9; c = 3e8; lambda = c/fc;
dx = lambda/4; dy = lambda/4;
Mx = 6; My = 6; M = Mx * My;
K = 4; B = 4;

delta_x = (2*[0:(Mx-1)].' - Mx + 1)/2;
delta_y = (2*[0:(My-1)].' - My + 1)/2;
Coor_Ele_init = [kron(delta_x*dx, ones(My,1)), kron(ones(Mx,1), delta_y*dy), zeros(M,1)];

Coor_Ele_init_r = [delta_x*dx, (-Mx+1)/2*(lambda/2)*ones(Mx,1), zeros(Mx,1)];
Coor_Ele_init_c = [(My-1)/2*dx*ones(My,1), delta_y*dy, zeros(My,1)];

delta_ele = 2^kappa0; d = dx * delta_ele;
Coor_Ele_init_r = Coor_Ele_init_r(1:delta_ele:end,:);
Coor_Ele_init_c = Coor_Ele_init_c(1:delta_ele:end,:);
M_x_r = size(Coor_Ele_init_r, 1);
M_y_c = size(Coor_Ele_init_c, 1);

delta_x_r = [(-M_x_r+1)/2:1:(M_x_r-1)/2].';
delta_y_c = [(-M_y_c+1)/2:1:(M_y_c-1)/2].';
Coor_Ele_init_x_d = [delta_x_r*(lambda/2), (-M_y_c+1)/2*(lambda/2)*ones(M_y_c,1), zeros(M_x_r,1)];
Coor_Ele_init_x_u = [delta_x_r*(lambda/2), (M_y_c-1)/2*(lambda/2)*ones(M_y_c,1), zeros(M_x_r,1)];
Coor_Ele_init_y_l = [(-M_x_r+1)/2*(lambda/2)*ones(M_x_r,1), delta_y_c*(lambda/2), zeros(M_x_r,1)];
Coor_Ele_init_y_r = [(M_x_r-1)/2*(lambda/2)*ones(M_x_r,1), delta_y_c*(lambda/2), zeros(M_x_r,1)];

Coor_Feed = [0.002984305671304, 0.001586131606604, 0];
Dis_Feed2Ele = sqrt((Coor_Ele_init - Coor_Feed).^2 * ones(3,1));
V_F0 = exp(-1j*2*pi*sqrt(3)/lambda * Dis_Feed2Ele);
eta0 = 8/3/M;
V_F = sqrt(eta0) * V_F0;

dmin = 0.1;

%% 主循环
fprintf('开始仿真...\n========================================\n');

for idx = 1:num_trials
    Pt = Pt_range(idx);
    Ptx = 10^((Pt - 30)/10);
    fprintf('进度: %d/%d | Pt = %d dBm\n', idx, num_trials, Pt);
    
    try
        % Step 1-3: 初始化、信道生成、感知
        [Coor_Ele, normal_vector, Rot_Matrix, ~, Coor_Ele_r, Coor_Ele_c, ...
         Coor_Ele_x_d, Coor_Ele_x_u, Coor_Ele_y_l, Coor_Ele_y_r] = ...
            Orientation_Initial(B, Coor_Ele_init, Coor_Ele_init_r, ...
                               Coor_Ele_init_c, 0, Coor_Ele_init_x_d, ...
                               Coor_Ele_init_x_u, Coor_Ele_init_y_l, ...
                               Coor_Ele_init_y_r);
        
        [h, e0, e2, ~, ~, ~, ~, Iota, eta, Omega, h_r, h_c, h_sen, ...
         h_x_d, h_x_u, h_y_l, h_y_r, Omega0] = ...
            Channel_Generation_init(K, B, M, Mx, My, M_x_r, M_y_c, ...
                                   normal_vector, K_rice, Coor_Ele, ...
                                   Coor_Ele_r, Coor_Ele_c, lambda, Rot_Matrix, ...
                                   Coor_Ele_x_d, Coor_Ele_x_u, ...
                                   Coor_Ele_y_l, Coor_Ele_y_r);
        
        [est, ~, ~] = Sensing_Algorithm(K, B, Iota, Mx, My, M_x_r, M_y_c, ...
                                       h_r, h_c, Rot_Matrix, e0, K_rice, ...
                                       lambda, e2, sigma0, Omega, eta, N, d, ...
                                       kappa0, h_sen, h_x_d, h_x_u, h_y_l, ...
                                       h_y_r, Coor_Ele_init_x_d, ...
                                       Coor_Ele_init_x_u, Coor_Ele_init_y_l, ...
                                       Coor_Ele_init_y_r, Omega0);
        
        % Step 4: 旋转优化
        for k = 1:K
            n1 = [0.02; 0.01; sqrt(1 - 0.02^2 - 0.01^2)];
            n2 = est{k};
            u = cross(n1, n2); u = u / norm(u);
            alph = acos(n1.' * n2);
            U = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
            Rot_Matrix{k} = cos(alph)*eye(3) + (1-cos(alph))*(u*u.') + sin(alph)*U;
            Coor_Ele{k} = Coor_Ele_init * Rot_Matrix{k}.';
            normal_vector_rot{k} = Rot_Matrix{k} * [0; 0; 1];
        end
        
        % Step 5-7: 优化与SWIPT
        [s, a_st, mm, Xtr] = optimize_system(B, K, M, Mx, My, est, ...
                                             normal_vector_rot, Rot_Matrix, ...
                                             Coor_Ele_init, lambda, dx, dy, ...
                                             delta_x, delta_y, dmin, V_F0, V_F, Ptx);
        
        H_eff = generate_effective_channel(s, K, B, M, Iota, Coor_Ele, est, ...
                                          normal_vector_rot, Rot_Matrix, ...
                                          Coor_Ele_init, e0, e2, eta, K_rice, ...
                                          Omega0, lambda, delta_x, delta_y, ...
                                          dx, dy, mm, V_F);
        H_eff = H_eff / sqrt(sigma0);
        
        psi1 = 0.1 * rand(K, K);
        psi2 = 0.01 * rand(K, 1) + 0.01j * rand(K, 1);
        rho = 0.1 * rand(K, 1);
        
        P_EH(idx) = fun_run_R00(R00, H_eff, rho, psi1, psi2, sigma0, Ptx, B, K);
        fprintf('  -> EH Power = %.6f W\n', P_EH(idx));
        
    catch ME
        fprintf('  -> 错误: %s\n', ME.message);
        P_EH(idx) = NaN;
    end
    fprintf('----------------------------------------\n');
end

%% 绘图
figure('Position', [100, 100, 800, 600]);
plot(Pt_range, P_EH*1000, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('Transmit Power (dBm)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('EH Power (mW)', 'FontSize', 14, 'FontWeight', 'bold');
title('Energy Harvesting Performance vs Transmit Power', 'FontSize', 16);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

save('EH_vs_Pt_results.mat', 'Pt_range', 'P_EH');
saveas(gcf, 'EH_vs_Pt.png');
saveas(gcf, 'EH_vs_Pt.fig');

fprintf('\n仿真完成！\n');