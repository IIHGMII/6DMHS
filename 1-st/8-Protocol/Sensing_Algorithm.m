function [est_G,est,est_error] = Sensing_Algorithm(K,B,Iota,Mx,My,M_x_r,M_y_c,h_r,h_c,Rot_Matrix,e0,K_rice,lambda,e2,sigma0,Omega,eta,N,d,kappa0)

% N = 32;    %FFT采样点数
P_s = sqrt( 1e-8 ); %全息感知器功率









delta_n = [ (-N+1) / 2 : 1 : (N-1) /2 ].';
F = exp( - 1j * 2 * pi * ( delta_n * delta_n.' ) / N   );
% % % % % % % % % % % % % % % % % % % % F_x =  exp( -1j * 2*pi * delta_x_r * delta_x_r.' / M_x_r ); %用于"行"感知单元感知信息提取的FFT矩阵
% % % % % % % % % % % % % % % % % % % % F_y =  exp( -1j * 2*pi * delta_y_c * delta_y_c.' / M_y_c ); %用于"列"感知单元感知信息提取的FFT矩阵



%% 扩展到和FFT矩阵相同维度
for k = 1:K
    for b = 1:B
        hx{k,b} = [ zeros((N-M_x_r)/2,1) ; h_r{k,b} ; zeros((N-M_x_r)/2,1) ];
        hy{k,b} = [ zeros((N-M_y_c)/2,1) ; h_c{k,b} ; zeros((N-M_y_c)/2,1) ];
    end
end

for k = 1:K
    for b = 1:B
        P_h(k,b) = sum(abs(hx{k,b}) + abs(hy{k,b}));
    end
    index_k(k) = find( P_h(k,:) == max(P_h(k,:)) ) ; %第k个用户最大增益信道对应的第b个全息超表面
end

% 求第k个用户在对应最大增益全息超表面的局部坐标系中的坐标
for k = 1:K
    e0_L{k} = Rot_Matrix{ index_k(k) }.' * e0{k};
    h_L_x{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda * (lambda/2)  * ( [ (-M_x_r+1)/2 : 1 : (M_x_r-1)/2 ].' * e0_L{k}(1) + (-M_y_c+1)/2 * e0_L{k}(2) )  );
    h_L_y{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda * (lambda/2)  * ( (M_x_r-1)/2 * e0_L{k}(1) + [ (-M_y_c+1)/2 : 1 : (M_y_c-1)/2 ].' * e0_L{k}(2)  ) );
    for iota = 1:Iota
        e2_L{k,iota} = Rot_Matrix{ index_k(k) }.' * e2{k,iota};
        if e2_L{k,iota}(3) >= 0
            h_L_x{k} = h_L_x{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda * (lambda/2) * ( [ (-M_x_r+1)/2 : 1 : (M_x_r-1)/2 ].' * e2_L{k,iota}(1) + (-M_y_c+1)/2 * e2_L{k,iota}(2) )  );
            h_L_y{k} = h_L_y{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda * (lambda/2) * ( (M_x_r-1)/2 * e2_L{k,iota}(1) + [ (-M_y_c+1)/2 : 1 : (M_y_c-1)/2 ].' * e2_L{k,iota}(2)  ) );
        end
    end
end

clear hx hy

% 当元素间隔超过半波长的时候，就会产生孪生波束，导致感知失败
% for p = 1:2:length(h_L_x{k})
%     h_L_x{k}(p) = 0;
%     h_L_y{k}(p) = 0;
% end

%% 添加全息感知

for k = 1:K
    for b = 1:B
        S_Holo_x{k,b} = [ zeros((N-M_x_r)/2,1) ; sqrt(P_s) * exp(1j * 2 * pi * rand(M_x_r,1)); zeros((N-M_x_r)/2,1) ];
        S_Holo_y{k,b} = [ zeros((N-M_x_r)/2,1) ; sqrt(P_s) * exp(1j * 2 * pi * rand(M_y_c,1)); zeros((N-M_x_r)/2,1) ];
    end
end


for k = 1:K
        hx{k} = [ zeros((N-M_x_r)/2,1) ; h_L_x{k}; zeros((N-M_x_r)/2,1) ];
        hy{k} = [ zeros((N-M_y_c)/2,1) ; h_L_y{k}; zeros((N-M_y_c)/2,1) ];

        yx{k} = hx{k} + sqrt(sigma0) * [ zeros((N-M_x_r)/2,1) ; 1/sqrt(2) * ( randn(M_x_r,1) + 1j * randn(M_x_r,1) ); zeros((N-M_x_r)/2,1) ];
        yy{k} = hy{k} + sqrt(sigma0) * [ zeros((N-M_y_c)/2,1) ; 1/sqrt(2) * ( randn(M_y_c,1) + 1j * randn(M_y_c,1) ); zeros((N-M_y_c)/2,1) ];

        Holo_Graph_x{k} = 2 * real( conj(yx{k}) .* S_Holo_x{k,index_k(k)} );
        Holo_Graph_y{k} = 2 * real( conj(yy{k}) .* S_Holo_y{k,index_k(k)} );

        Holo_Graph_Recov_x{k} = Holo_Graph_x{k} .* S_Holo_x{k,index_k(k)};
        Holo_Graph_Recov_y{k} = Holo_Graph_y{k} .* S_Holo_y{k,index_k(k)};
end



for k = 1:K
%     Hx{k} = F * hx{k};  Hy{k} = F * hy{k};
    Hx{k} = F * Holo_Graph_Recov_x{k};  Hy{k} = F * Holo_Graph_Recov_y{k};
    index_x_Sens(k) = find(abs(Hx{k}) == max(abs(Hx{k})));
    index_y_Sens(k) = find(abs(Hy{k}) == max(abs(Hy{k})));
    est_temp{k}(1) = delta_n(index_x_Sens(k)) / N * 2;
    est_temp{k}(2) = delta_n(index_y_Sens(k)) / N * 2;
    est{k} = [est_temp{k}(1) ; est_temp{k}(2) ; sqrt( max( [ 1 - est_temp{k}(1)^2 - est_temp{k}(2)^2 , 0 ] ) )   ];
end

for k = 1:K
est_error(k,1) = norm( Rot_Matrix{ index_k(k) } * est{k} - e0{k} , 2);
est_G{k} = Rot_Matrix{ index_k(k) } * est{k};
end




end