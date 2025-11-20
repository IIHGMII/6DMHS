function [est_G,est,est_error] = Sensing_Algorithm(K,B,Iota,Mx,My,M_x_r,M_y_c,h_r,h_c,Rot_Matrix,e0,K_rice,lambda,e2,sigma0,Omega,eta,N,d,kappa0,h_sen,h_x_d , h_x_u , h_y_l , h_y_r , Coor_Ele_init_x_d , Coor_Ele_init_x_u , Coor_Ele_init_y_l , Coor_Ele_init_y_r)

% N = 32;    %FFT采样点数
P_s = sqrt( 1e-4 ); %全息感知器功率








% % N=128;
delta_n = [ (-N+1) / 2 : 1 : (N-1) /2 ].';
F = exp( - 1j * 2 * pi * ( delta_n * delta_n.' ) / N   );
% % % % % % % % % % % % % % % % % % % % F_x =  exp( -1j * 2*pi * delta_x_r * delta_x_r.' / M_x_r ); %用于"行"感知单元感知信息提取的FFT矩阵
% % % % % % % % % % % % % % % % % % % % F_y =  exp( -1j * 2*pi * delta_y_c * delta_y_c.' / M_y_c ); %用于"列"感知单元感知信息提取的FFT矩阵

for k = 1:K
    for b = 1:B
        P_h(k,b) = sum( abs(h_sen{k,b}) , 'all' );
    end
    index_k(k) = ( find( P_h(k,:) == max(P_h(k,:)) ) ); %第k个用户最大增益信道对应的第b个全息超表面
end

%% 变换到局部坐标系
%% x,d
for k = 1:K
    e0_L{k} = Rot_Matrix{ index_k(k) }.' * e0{k};
    h_x_d_L{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda   * ( Coor_Ele_init_x_d * e0_L{k} )  );
    h_x_u_L{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda   * ( Coor_Ele_init_x_u * e0_L{k}  ) );
    h_y_l_L{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda   * ( Coor_Ele_init_y_l * e0_L{k} )  );
    h_y_r_L{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda   * ( Coor_Ele_init_y_r * e0_L{k} )  );

    for iota = 1:Iota
        e2_L{k,iota} = Rot_Matrix{ index_k(k) }.' * e2{k,iota};
        if e2_L{k,iota}(3) >= 0
            h_x_d_L{k} = h_x_d_L{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda  * ( Coor_Ele_init_x_d * e2_L{k,iota}  ) );
            h_x_u_L{k} = h_x_u_L{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda  * ( Coor_Ele_init_x_u * e2_L{k,iota}  ) );
            h_y_l_L{k} = h_y_l_L{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda  * ( Coor_Ele_init_y_l * e2_L{k,iota}  ) );
            h_y_r_L{k} = h_y_r_L{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda  * ( Coor_Ele_init_y_r * e2_L{k,iota}  ) );
        end
    end
end

for k = 1:K
%     h_x_d_L_exp{k} = [ zeros((N-M_x_r)/2,1) ; h_x_d_L{k} ; zeros((N-M_x_r)/2,1) ];
%     h_x_u_L_exp{k} = [ zeros((N-M_x_r)/2,1) ; h_x_u_L{k} ; zeros((N-M_x_r)/2,1) ];
%     h_y_l_L_exp{k} = [ zeros((N-M_x_r)/2,1) ; h_y_l_L{k} ; zeros((N-M_x_r)/2,1) ];
%     h_y_r_L_exp{k} = [ zeros((N-M_x_r)/2,1) ; h_y_r_L{k} ; zeros((N-M_x_r)/2,1) ];
    H_Sen{k} = sparse(M_x_r,M_x_r);
    H_Sen{k}(1,:)   = h_x_d_L{k}.';
    H_Sen{k}(end,:) = h_x_u_L{k}.';
    H_Sen{k}(:,1)   = h_y_l_L{k};
    H_Sen{k}(:,end) = h_y_r_L{k};

    H_SEN{k} = zeros(N);
    kk = (N - 2)/2;    % k=2
    H_SEN{k}(kk+1:kk+M_x_r, kk+1:kk+M_x_r) = H_Sen{k};


    A{k} = abs( F * H_SEN{k} *F.' );
    index_x_Sen(k) = find( A{k} == max(max( A{k} )) );
end


for k = 1:K
    est{k}(1,1) = delta_n( floor( index_x_Sen(k) / N ) ) /N * 2;
    est{k}(2,1)    = delta_n( index_x_Sen(k) - floor( index_x_Sen(k) / N )*N) /N *2;
    est{k}(3,1) = sqrt( max( [ 1 - est{k}(1)^2 - est{k}(2)^2 , 0 ] ) );
    est_error(k,1) = norm( Rot_Matrix{ index_k(k) } * est{k} - e0{k} , 2);
    est_G{k} = Rot_Matrix{ index_k(k) } * est{k};
end



end