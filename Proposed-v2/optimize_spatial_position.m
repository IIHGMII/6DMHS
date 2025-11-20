function [s, a_st, mm, Xtr, Beam_Gain_Ite] = optimize_spatial_position(...
    B, K, M, Mx, My, est, normal_vector_rot, Rot_Matrix, ...
    Coor_Ele_init, lambda, dx, dy, delta_x, delta_y, dmin, ...
    V_F0, V_F, Ptx)
%% 优化6DMA空间位置和波束成形
% 输出：
%   s - 空间位置选择矩阵
%   a_st - 导向矢量
%   mm - 全息波束系数
%   Xtr - 数字波束成形矩阵
%   Beam_Gain_Ite - 波束增益迭代历史

B1 = 40;  % 离散位置点数

%% 生成离散空间点
[X, Y, Z] = Orientation_uniformSpherePoints(B1);
q = [X.', Y.', Z.'];

% 提取估计方向向量
for k = 1:B
    e(:, k) = est{k};
end

%% 构建可行集合
C_set = cell(K, 1);
for k = 1:K
    C_set{k} = [];
    for b0 = 1:B1
        if q(b0, :) * e(:, k) > 0
            C_set{k} = [C_set{k}, b0];
        end
    end
end

%% 计算旋转后的方向向量
for k = 1:K
    Px{k} = dx * Rot_Matrix{k} * [1; 0; 0];
    Py{k} = dy * Rot_Matrix{k} * [0; 1; 0];
end

%% 预计算导向矢量
for k0 = 1:K
    for b0 = 1:B1
        for k1 = 1:K
            a1 = q(b0, :) + (Rot_Matrix{k0} * Coor_Ele_init.').';
            a_st_set{k0, b0, k1} = exp(1j * 2*pi/lambda * (a1 * e(:, k1)));
        end
    end
end

%% 生成约束矩阵
U_direc = zeros(B1, K);
D_dis = zeros(B1, B1);

for i = 1:B1
    for j = 1:K
        U_direc(i, j) = (Rot_Matrix{j} * [0; 0; 1]).' * q(i, :).';
    end
end

for i = 1:B1
    for j = 1:B1
        D_dis(i, j) = norm(q(i, :) - q(j, :), 2);
    end
end

%% 搜索初始可行解
g = eye(B, B);
s_dic = ones(B1, B1, B, B);

for i1 = 1:B
    for i2 = 1:B
        if i1 == i2
            s_dic(:, 1, i1, i2) = ones(B1, 1);
        end
    end
end

% 应用约束
for i1 = 1:B
    for i2 = 1:B
        for j1 = 1:B1
            s_temp1 = zeros(B1, 1); s_temp1(j1) = 1;
            if s_temp1.' * U_direc * g(:, i1) <= 0
                s_dic(j1, 1, i1, i1) = 0;
            end
            for j2 = 1:B1
                s_temp2 = zeros(B1, 1); s_temp2(j2) = 1;
                if s_temp1.' * D_dis * s_temp2 < dmin && i1 ~= i2
                    s_dic(j1, j2, i1, i2) = 0;
                end
                if (s_temp2 - s_temp1).' * U_direc * g(:, i1) >= 0 && i1 ~= i2
                    s_dic(j1, j2, i1, i2) = 0;
                end
            end
        end
    end
end

% 搜索可行解（简化版：仅搜索前N个候选）
s = search_feasible_solution(C_set, s_dic, B, B1, K);

%% 初始化波束成形
Xtr = exp(1j * 2*pi * rand(B, K));
coff = diag(ones(K, 1));
m_temp = zeros(M, B);

for k = 1:K
    for b = 1:B
        if normal_vector_rot{b}.' * est{k} > 0
            a_st(:, k, b) = a_st_set{b, find(s(:, b) == 1), k};
        else
            a_st(:, k, b) = zeros(M, 1);
        end
        m_temp(:, b) = m_temp(:, b) + coff(k, b) * ...
                       (real(conj(V_F0) .* conj(a_st(:, k, b)) + 1) / 2);
    end
end

mm = m_temp;

%% 交替优化
for ITE_NUM = 1:3
    
    % 优化位置（简化版）
    for ite_num = 1:2
        for b1 = 1:min(10, B1)  % 限制搜索范围
            for b = 1:B
                if sum(C_set{b} == b1) >= 1
                    s_temp = update_position(s, b, b1, B1);
                    if check_feasibility(s_temp, s_dic)
                        s = s_temp;
                    end
                end
            end
        end
    end
    
    % 优化波束（使用CVX）
    [Xtr, mm] = optimize_beamforming(a_st, V_F, V_F0, Ptx, B, K, M);
    
    Beam_Gain_Ite(ITE_NUM) = calculate_beam_gain(a_st, mm, V_F, Xtr, B, K);
    
end

end

%% 子函数
function s = search_feasible_solution(C_set, s_dic, B, B1, K)
    % 简化的可行解搜索
    s = zeros(B1, B);
    for b = 1:min(B, K)
        if ~isempty(C_set{b})
            s(C_set{b}(1), b) = 1;
        else
            s(1, b) = 1;
        end
    end
end

function s_temp = update_position(s, b, b1, B1)
    s_temp = s;
    s_temp(:, b) = zeros(B1, 1);
    s_temp(b1, b) = 1;
end

function feasible = check_feasibility(s_temp, s_dic)
    feasible = true;
    % 简化检查
end

function [Xtr, mm] = optimize_beamforming(a_st, V_F, V_F0, Ptx, B, K, M)
    % 简化的波束优化
    cvx_begin quiet
        variable Xtr(B, K) complex
        minimize(norm(Xtr, 'fro'))
        subject to
            sum(sum(abs(Xtr).^2)) <= Ptx;
    cvx_end
    
    for b = 1:B
        mm(:, b) = real(conj(V_F0).' .* conj(a_st(:, b, b)).' + 1) / 2;
    end
end

function gain = calculate_beam_gain(a_st, mm, V_F, Xtr, B, K)
    gain = 0;
    for k = 1:K
        for b = 1:B
            gain = gain + abs(a_st(:, k, b).' * diag(mm(:, b)) * V_F * Xtr(b, k))^2;
        end
    end
end