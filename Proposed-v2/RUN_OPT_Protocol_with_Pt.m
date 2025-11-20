function [P_out] = RUN_OPT_Protocol_with_Pt(Pt, K_rice)
% 修改版的RUN_OPT_Protocol函数，接受Pt作为参数

% 转换发射功率
Ptx = 10^((Pt - 30) / 10);

% 设置其他参数
kappa0 = 1; % 0 1 2 3 4
N = 1024;
sigma0 = 1e-12;
R00 = 1;

format long

% 初始化参数
for nummmmmm = 1
    fc = 30e9; %载波频率
    c = 3e8; %光速
    lambda = c / fc; %波长
    dx = lambda / 4; %x轴方向天线间隔
    dy = lambda / 4; %y轴方向天线间隔
    Mx = 32;
    My = 32;
    M = Mx * My;
    Q = 1; %每个全息超表面馈源个数
    Dx = Mx * dx;
    Dy = My * dy;
    K = 4; %用户个数
    B = 4; %全息超表面个数，每个全息超表面上配备有M个元素
    NUM_Position = 100; %6DMA空间位置
    delta_x = (2 * [0:(Mx - 1)].' - Mx + 1) / 2;
    delta_y = (2 * [0 : (My - 1)].' - My + 1) / 2;
    Coor_Ele_init = [kron(delta_x*dx, ones(My, 1)), kron(ones(Mx, 1), delta_y*dy), zeros(M, 1)];

    Coor_Ele_init_r = [delta_x * dx, (-Mx + 1) / 2 * (lambda / 2) * ones(Mx, 1), zeros(Mx, 1)];
    Coor_Ele_init_c = [(My - 1) / 2 * dx * ones(My, 1), delta_y * dy, zeros(My, 1)];

    delta_ele = 2^kappa0;
    d = dx * delta_ele;
    Coor_Ele_init_r = Coor_Ele_init_r(1:delta_ele:end, :);
    Coor_Ele_init_c = Coor_Ele_init_c(1:delta_ele:end, :);
    M_x_r = size(Coor_Ele_init_r, 1);
    M_y_c = size(Coor_Ele_init_c, 1);
    dmin = 0.1;

    % 引入全息超表面传输模型
    Coor_Feed = [0.002984305671304, 0.001586131606604, 0]; %馈源坐标
    Dis_Feed2Ele = sqrt((Coor_Ele_init - Coor_Feed).^2*ones(3, 1)); %馈源到元素的距离
    V_F0 = exp(-1j*2*pi*sqrt(3)/lambda*Dis_Feed2Ele); %馈源响应向量
    eta0 = 8 / 3 / M;
    V_F = sqrt(eta0) * V_F0;
end

for num = 1:1
    % 6DMA初始空间姿态
    [Coor_Ele, normal_vector, Rot_Matrix, q, Coor_Ele_x, Coor_Ele_y, Coor_Ele_x_d, Coor_Ele_x_u, Coor_Ele_y_l, Coor_Ele_y_r] = Orientation_Initial(B, Coor_Ele_init, Coor_Ele_init_r, Coor_Ele_init_c, 0, [], [], [], []);

    % 针对初始6DMA位置，生成无线信道
    [h, e0, e2, theta0, phi0, theta_scatterer0, phi_scatterer0, Iota, eta, Omega, h_r, h_c] = Channel_Generation_init(K, B, M, Mx, My, M_x_r, M_y_c, normal_vector, K_rice, Coor_Ele, Coor_Ele_r, Coor_Ele_c, lambda, Rot_Matrix);

    % 生成全息感知
    [est, est_L, est_error] = Sensing_Algorithm(K, B, Iota, Mx, My, M_x_r, M_y_c, h_r, h_c, Rot_Matrix, e0, K_rice, lambda, e2, sigma0, Omega, eta, N, d, kappa0);
end

% 调整旋转角度
for k = 1:K
    n1 = [0.02; 0.01; sqrt(1 - 0.02^2 - 0.01^2)]; %最大增益方向
    n2 = est{k}; %用户方向
    u = cross(n1, n2);
    u = u / norm(u);
    % 计算旋转角度（弧度单位）
    alph = acos(n1.'*n2);
    % 计算矩阵U
    U = [0, -u(3), u(2); ...
        u(3), 0, -u(1); ...
        -u(2), u(1), 0];
    % 计算罗德格里斯旋转矩阵
    Rot_Matrix{k} = cos(alph) * eye(3) + (1 - cos(alph)) * (u * u.') + sin(alph) * U;
    Coor_Ele{k} = Coor_Ele_init * Rot_Matrix{k}.'; %旋转到第k个用户对应方向后的用户坐标
    normal_vector_rot{k} = Rot_Matrix{k} * [0; 0; 1]; % 旋转后的法向量
end

% 生成B个离散空间点的波束增益向量
B1 = 32; %一共有64个空间离散点位
for k = 1:B, e(:, k) = est{k}; end

[X, Y, Z] = Orientation_uniformSpherePoints(B1);
q = [X.', Y.', Z.'];

% 第k个用户支持的离散位置集合
C_set = cell(K, 1);
for k = 1:K
    C_set{k} = [];
    for b0 = 1:B1
        if q(b0, :) * e(:, k) > 0
            C_set{k} = [C_set{k}, b0];
        end
    end
end

for k = 1:K
    %注意，dx = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2
    Px{k} = dx * Rot_Matrix{k} * [1; 0; 0];
    Py{k} = dy * Rot_Matrix{k} * [0; 1; 0];
end

% RHS对准第k0个用户，位于第b0个位置，在第k1个用户方向的导向矢量
for k0 = 1:K
    for b0 = 1:B1
        for k1 = 1:K
            a1 = q(b0, :) + (Rot_Matrix{k0} * Coor_Ele_init.').';
            a_st_set{k0, b0, k1} = exp(1j*2*pi/lambda*(a1 * e(:, k1)));
        end
    end
end

% 生成方向矩阵和距离矩阵
U_direc = zeros(B1, K);
D_dis = zeros(B1, B1);

for i = 1:B1
    for j = 1:K
        U_direc(i , j) = (Rot_Matrix{j} * [0; 0; 1]).' * q(i, :).'; %此处根据导向矢量求
    end
end
for i = 1:B1
    for j = 1:B1
        D_dis(i, j) = norm(q(i, :)-q(j, :), 2);
    end
end

% 初始空间姿态
for b = 1:B
    s(:, b) = [zeros(b-1, 1); 1; zeros(B1-(b - 1)-1, 1)];
end
g = eye(B, B);

% 制备b个平移向量的可行区间,搜索初始可行解
s_dic = zeros(B1, B1, B, B);
for i1 = 1:B
    for i2 = 1:B
        if i1 == i2, s_dic(:, 1, i1, i2) = ones(B1, 1);
        else, s_dic(:, :, i1, i2) = ones(B1, B1);
        end
    end
end

% 优化部分
psi = 0.1 * rand(K, K);
psi1 = 0.05 * randn(K, K) + 0.05j * randn(K, K);
psi2 = 0.01 * rand(K, 1) + 0.01j * rand(K, 1);
rho = 0.01 * rand(K, 1); %功率分割因子,rho用于传能，(1-rho)用于通信
P_out_ite = [];

% 计算有效信道
H_eff = zeros(B, K);
for k = 1:K
    for b = 1:B
        H_eff(b, k) = V_F' * h{b, k};
    end
end

% 使用fun_run_R00函数计算EH功率
P_out = fun_run_R00(R00, H_eff, rho, psi1, psi2, sigma0, Ptx, B, K);

end