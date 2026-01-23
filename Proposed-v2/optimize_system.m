function [s, a_st, mm, Xtr] = optimize_system(B, K, M, Mx, My, est, ...
                                              normal_vector_rot, Rot_Matrix, ...
                                              Coor_Ele_init, lambda, dx, dy, ...
                                              delta_x, delta_y, dmin, V_F0, V_F, Ptx)
%% 联合优化空间位置和波束成形（修复CVX错误）

B1 = 40;  % 离散位置数

%% 生成离散空间点
[X, Y, Z] = Orientation_uniformSpherePoints(B1);
q = [X.', Y.', Z.'];

for k = 1:B
    e(:, k) = est{k};
end

%% 构建可行集
C_set = cell(K, 1);
for k = 1:K
    C_set{k} = [];
    for b0 = 1:B1
        if q(b0, :) * e(:, k) > 0
            C_set{k} = [C_set{k}, b0];
        end
    end
end

%% 计算旋转向量
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

%% 约束矩阵
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

%% 构建可行域字典
g = eye(B, B);
s_dic = zeros(B1, B1, B, B);

for i1 = 1:B
    for i2 = 1:B
        if i1 == i2
            s_dic(:, 1, i1, i2) = ones(B1, 1);
        else
            s_dic(:, :, i1, i2) = ones(B1, B1);
        end
    end
end

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

%% 搜索初始可行解
LoopUntil = 0;
for m1 = C_set{1}
    for m2 = C_set{2}(C_set{2}~=m1)
        for m3 = C_set{3}(C_set{3}~=m1 & C_set{3}~=m2)
            for m4 = C_set{4}(C_set{4}~=m1 & C_set{4}~=m2 & C_set{4}~=m3)
                if s_dic(m1,1,1,1) + s_dic(m2,1,2,2) + s_dic(m3,1,3,3) + s_dic(m4,1,4,4) + ...
                   s_dic(m1,m2,1,2) + s_dic(m1,m3,1,3) + s_dic(m1,m4,1,4) + ...
                   s_dic(m2,m1,2,1) + s_dic(m2,m3,2,3) + s_dic(m2,m4,2,4) + ...
                   s_dic(m3,m1,3,1) + s_dic(m3,m2,3,2) + s_dic(m3,m4,3,4) + ...
                   s_dic(m4,m1,4,1) + s_dic(m4,m2,4,2) + s_dic(m4,m3,4,3) == 16
                    LoopUntil = 1;
                    break;
                end
            end
            if LoopUntil == 1, break; end
        end
        if LoopUntil == 1, break; end
    end
    if LoopUntil == 1, break; end
end

s = zeros(B1, B);
if LoopUntil == 1
    s(m1,1) = 1; s(m2,2) = 1; s(m3,3) = 1; s(m4,4) = 1;
else
    % 备用方案
    for b = 1:B
        if ~isempty(C_set{b})
            s(C_set{b}(1), b) = 1;
        else
            s(1, b) = 1;
        end
    end
end

%% 初始化波束
Xtr = exp(1j * 2*pi * rand(B, K));
coff = diag(ones(K, 1));
m_temp = zeros(M, B);

for k = 1:K
    for b = 1:B
        if normal_vector_rot{b}.' * est{k} > 0
            a_st(:, k, b) = a_st_set{b, find(s(:, b)==1), k};
        else
            a_st(:, k, b) = zeros(M, 1);
        end
        m_temp(:, b) = m_temp(:, b) + coff(k, b) * ...
                       (real(conj(V_F0) .* conj(a_st(:, k, b)) + 1) / 2);
    end
end

mm = m_temp;

%% 交替优化（简化版）
psi = 0.1 * rand(K, K);

for ITE_NUM = 1:2  % 减少迭代次数加速
    
    % 位置优化（简化）
    Gain_Beam = calculate_gain(a_st, mm, V_F, Xtr, B, K);
    Gain_old = min(Gain_Beam);
    
    for b1 = 1:min(20, B1)  % 限制搜索范围
        for b = 1:B
            if sum(C_set{b} == b1) >= 1
                s_temp = s;
                s_temp(:, b) = zeros(B1, 1);
                s_temp(b1, b) = 1;
                
                % 检查可行性
                m1 = find(s_temp(:,1)==1); m2 = find(s_temp(:,2)==1);
                m3 = find(s_temp(:,3)==1); m4 = find(s_temp(:,4)==1);
                
                if s_dic(m1,m2,1,2) + s_dic(m1,m3,1,3) + s_dic(m1,m4,1,4) + ...
                   s_dic(m2,m1,2,1) + s_dic(m2,m3,2,3) + s_dic(m2,m4,2,4) + ...
                   s_dic(m3,m1,3,1) + s_dic(m3,m2,3,2) + s_dic(m3,m4,3,4) + ...
                   s_dic(m4,m1,4,1) + s_dic(m4,m2,4,2) + s_dic(m4,m3,4,3) == 12
                    
                    % 更新导向矢量
                    for k = 1:K
                        for bb = 1:B
                            if normal_vector_rot{bb}.' * est{k} > 0
                                a_st_temp(:,k,bb) = a_st_set{bb, find(s_temp(:,bb)==1), k};
                            else
                                a_st_temp(:,k,bb) = zeros(M,1);
                            end
                        end
                    end
                    
                    Gain_new = calculate_gain(a_st_temp, mm, V_F, Xtr, B, K);
                    if min(Gain_new) > Gain_old
                        s = s_temp;
                        a_st = a_st_temp;
                        Gain_old = min(Gain_new);
                    end
                end
            end
        end
    end
    
    % 波束优化
    for b = 1:B
        mm(:, b) = real(conj(V_F0) .* conj(a_st(:, b, b)) + 1) / 2;
    end
    
    % 数字波束优化（修复CVX错误）
    try
        cvx_begin quiet
            variable Xtr(B, K) complex
            variable P0(1)
            
            expressions gain_beam_temp(K, B, K) gain_vec(K, K) Gain_Beam(K, 1)
            
            for k0 = 1:K
                for b0 = 1:B
                    for k1 = 1:K
                        gain_beam_temp(k1, b0, k0) = a_st(:, k0, b0).' * ...
                                                     diag(V_F * Xtr(b0, k1)) * mm(:, b0);
                    end
                end
            end
            
            for k = 1:K
                for k1 = 1:K
                    gain_vec(k1, k) = ones(1, B) * gain_beam_temp(k1, :, k)';
                end
            end
            
            for k = 1:K
                Gain_Beam(k) = 2 * real(psi(:, k)' * gain_vec(:, k)) - psi(:, k)' * psi(:, k);
            end
            
            maximize P0
            subject to
                Gain_Beam >= P0;
                sum(sum(pow_abs(Xtr, 2))) <= Ptx;  % 修复：使用pow_abs
        cvx_end
        
        % 更新辅助变量
        for k = 1:K
            psi(:, k) = gain_vec(:, k);
        end
        
    catch
        % CVX失败则保持原Xtr
        warning('CVX优化失败，使用初始值');
    end
    
end

end

%% 辅助函数
function Gain_Beam = calculate_gain(a_st, mm, V_F, Xtr, B, K)
    Gain_Beam = zeros(K, 1);
    for k = 1:K
        for k1 = 1:K
            for b = 1:B
                gain_beam_temp(k1, b, k) = a_st(:, k, b).' * diag(mm(:, b)) * V_F * Xtr(k1, b);
            end
        end
        Gain_Beam(k) = sum(abs(gain_beam_temp(:, :, k) * ones(B, 1)).^2);
    end
end
