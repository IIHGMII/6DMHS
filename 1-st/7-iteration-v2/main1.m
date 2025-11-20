% 6DMA,使用全息感知方法
clc;clear;

tic

kappa0 = 3; % 2 3 5 9 17



Pt = 40;
Ptx = 10^( (Pt - 30)/10 ); sigma0 = 1e-12; R00 = 1;
% kappa0 = 3; % 2 3 5 9 17

format long

% %% ****************************************初始化参数************************************************************************************************************************************************************

for nummmmmm = 1

    % Ptx = 10;
    fc = 30e9;   %载波频率
    c = 3e8;      %光速
    lambda = c/fc;   %波长
    dx = lambda/4;   %x轴方向天线间隔
    dy = lambda/4;   %y轴方向天线间隔
    Mx = 32; My = 32; M = Mx * My;
    Q  = 1; %每个全息超表面馈源个数
    Dx = Mx * dx;   Dy = My * dy;
    K_rice = 10;
    K = 8; %用户个数
    B = 8; %全息超表面个数，每个全息超表面上配备有M个元素
    NUM_Position = 100; %6DMA空间位置
    delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [ 0 : ( My - 1 ) ].' - My + 1 )/2;
    % delta_x = [0:(Mx-1)].';   delta_y = [0:(My-1)].';
    % Coor_Ele_init = [ kron( ones(My,1) , delta_x*dx ) , kron( delta_y*dy , ones(Mx,1) ) , zeros(M,1)  ];  %辐射元素坐标
    Coor_Ele_init = [ kron( delta_x*dx , ones(My,1) ) , kron( delta_x*dx , ones(My,1) ) , zeros(M,1)  ];  %辐射元素坐标
    % 标记横轴为x轴，纵轴为y轴
    M_x_r = Mx/2;      delta_x_r = [ (-M_x_r + 1)/2 : 1 : ( M_x_r - 1 )/2 ].' ;
    M_y_c = My/2;      delta_y_c = [ (-M_y_c + 1)/2 : 1 : ( M_y_c - 1 )/2 ].' ;
    Coor_Ele_init_r = [ delta_x_r * (lambda)/2 , (-M_x_r + 1)/2 * (lambda/2) * ones(M_x_r,1) , zeros(M_x_r,1) ];    %第0行辐射元素的坐标
    Coor_Ele_init_c = [ ( M_y_c - 1 )/2 * (lambda/2) * ones(M_y_c,1) , delta_y_c * (lambda/2) , zeros(M_y_c,1) ];    %第0列辐射元素的坐标
%     Coor_Ele_init_r = Coor_Ele_init_r_temp(1:2:end,:);
%     Coor_Ele_init_c = Coor_Ele_init_c_temp(1:2:end,:);
    dmin = 0.1;
% figure(1); plot( Coor_Ele_init(:,1) , Coor_Ele_init(:,2) , 'o' )
% % 
%% 引入全息超表面传输模型
Coor_Feed    = [ 0.002984305671304 , 0.001586131606604 , 0 ];   %馈源坐标
% Coor_Feed    = [ 0 , 3.838239570983593e-04 , 0 ];
Dis_Feed2Ele = sqrt( ( Coor_Ele_init - Coor_Feed ).^2 * ones(3,1) );    %馈源到元素的距离
V_F0          = exp( - 1j * 2*pi*sqrt(3)/lambda * Dis_Feed2Ele );    %馈源响应向量
eta0 = 8/3/M;
V_F = sqrt(eta0) * V_F0;

end



for  num = 1 : 1
%% 6DMA初始空间姿态
% Coor_Ele表示6DMA天线坐标，normal_vector是UPA的法向量
[ Coor_Ele , normal_vector , Rot_Matrix , Coor_q_init , Coor_Ele_r,Coor_Ele_c] = Orientation_Initial(B,Coor_Ele_init,Coor_Ele_init_r,Coor_Ele_init_c,0);
%% 针对初始6DMA位置，生成无线信道
[ h , e0 , e2 , theta0, phi0 , theta_scatterer0 , phi_scatterer0 , Iota , eta , Omega , h_r , h_c ] = Channel_Generation_init( K,B,M,Mx,My,normal_vector,K_rice,Coor_Ele,Coor_Ele_r,Coor_Ele_c,lambda , Rot_Matrix );
%% 生成全息感知
% [ error_theta , error_phi , error_varphi , est_direc_vector_rec , Rot_Matrix , normal_vector_rot ] = Sensing_Holographic_HalfWave1(B,M,Mx,My,lambda,delta_x,delta_y,dx,dy,K,h,Rot_Matrix,e0,kappa0);
% est_direc_vector_rec：指向用户的方向向量
% normal_vector_rot：全息超表面的法向量
% Rot_Matrix：旋转矩阵
end

% error_theta
% error_phi
% error_varphi
% E = [error_theta;error_phi;error_varphi];
% MSE = (sum( E.^2 , 'all' )) / 1 / 3 / K 

N = 128;    %FFT采样点数
delta_n = [ (-N+1) / 2 : 1 : (N-1) /2 ].';
F = exp( - 1j * 2 * pi * delta_n * delta_n.' / N   );
% F_x =  exp( -1j * 2*pi * delta_x_r * delta_x_r.' / M_x_r ); %用于"行"感知单元感知信息提取的FFT矩阵
% F_y =  exp( -1j * 2*pi * delta_y_c * delta_y_c.' / M_y_c ); %用于"列"感知单元感知信息提取的FFT矩阵



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
    index_k(k) = find( P_h(k,:) == max(P_h(k,:)) ); %第k个用户最大增益信道对应的第b个全息超表面
end

% 求第k个用户在对应最大增益全息超表面的局部坐标系中的坐标
for k = 1:K
    e0_L{k} = Rot_Matrix{ index_k(k) }.' * e0{k};
    h_L_x{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda * (lambda/2) * ( [ (-M_x_r+1)/2 : 1 : (M_x_r-1)/2 ].' * e0_L{k}(1) + (-M_y_c+1)/2 * e0_L{k}(2) )  );
    h_L_y{k} = sqrt( K_rice / (K_rice + 1) ) * sqrt(Omega(1)) * exp( 1j * 2 * pi / lambda * (lambda/2) * ( (-M_x_r+1)/2 * e0_L{k}(1) + [ (-M_y_c+1)/2 : 1 : (M_y_c-1)/2 ].' * e0_L{k}(2)  ) );
    for iota = 1:Iota
        e2_L{k,iota} = Rot_Matrix{ index_k(k) }.' * e2{k,iota};
        if e2_L{k,iota}(3) >= 0
            h_L_x{k} = h_L_x{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda * (lambda/2) * ( [ (-M_x_r+1)/2 : 1 : (M_x_r-1)/2 ].' * e2_L{k,iota}(1) + (-M_y_c+1)/2 * e2_L{k,iota}(2) )  );
            h_L_y{k} = h_L_y{k} + sqrt( 1 / (K_rice + 1) ) * eta(k,iota) * exp( 1j * 2 * pi / lambda * (lambda/2) * ( (-M_x_r+1)/2 * e2_L{k,iota}(1) + [ (-M_y_c+1)/2 : 1 : (M_y_c-1)/2 ].' * e2_L{k,iota}(2)  ) );
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
P_s = sqrt( 0.0001 ); %全息感知器功率

for k = 1:K
    for b = 1:B
        S_Holo_x{k,b} = [ zeros((N-M_x_r)/2,1) ; sqrt(P_s) * exp(1j * 2 * pi * rand(M_x_r,1)); zeros((N-M_x_r)/2,1) ];
        S_Holo_y{k,b} = [ zeros((N-M_x_r)/2,1) ; sqrt(P_s) * exp(1j * 2 * pi * rand(M_y_c,1)); zeros((N-M_x_r)/2,1) ];
    end
end

sigma0 = 1e-12;

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
    est{k} = [est_temp{k}(1) ; est_temp{k}(2) ; sqrt( 1 - est_temp{k}(1)^2 - est_temp{k}(2)^2 )];
end

for k = 1:K
%     disp('下一个')

norm( Rot_Matrix{ index_k(k) } * est{k} - e0{k} , 2)

end




toc