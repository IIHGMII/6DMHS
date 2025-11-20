function H_eff = generate_effective_channel(s, K, B, M, Iota, Coor_Ele, ...
                                            est, normal_vector_rot, ...
                                            Rot_Matrix, Coor_Ele_init, ...
                                            e0, e2, eta, K_rice, Omega0, ...
                                            lambda, delta_x, delta_y, dx, dy, ...
                                            mm, V_F)
%% 生成等效信道
% 输出：H_eff - B×K等效信道矩阵

% 离散空间点
B1 = size(s, 1);
[X, Y, Z] = Orientation_uniformSpherePoints(B1);
q = [X.', Y.', Z.'];

% 计算旋转向量
for k = 1:K
    Px{k} = dx * Rot_Matrix{k} * [1; 0; 0];
    Py{k} = dy * Rot_Matrix{k} * [0; 1; 0];
end

% 生成完整信道
H = zeros(B*M, K);

for k = 1:K
    for iota = 0:Iota
        for b = 1:B
            a1 = s(:, b).' * q + Coor_Ele{b};
            
            if iota == 0
                if normal_vector_rot{b}.' * e0{k} > 0
                    H((b-1)*M+1:b*M, k) = H((b-1)*M+1:b*M, k) + ...
                        sqrt(K_rice/(1+K_rice)) * sqrt(Omega0(k,b)) * ...
                        exp(1j * 2*pi/lambda * (a1 * e0{k}));
                end
            else
                if normal_vector_rot{b}.' * e2{k, iota} > 0
                    H((b-1)*M+1:b*M, k) = H((b-1)*M+1:b*M, k) + ...
                        sqrt(1/(1+K_rice)) * sqrt(Omega0(k,b)) * ...
                        eta(k, iota) * exp(1j * 2*pi/lambda * (a1 * e2{k, iota}));
                end
            end
        end
    end
end

% 计算等效信道
H_eff = zeros(B, K);
for b = 1:B
    for k = 1:K
        H_eff(b, k) = H((b-1)*M+1:b*M, k).' * diag(mm(:, b)) * V_F;
    end
end

end