function [P_out] = untitled( K_rice )
K_rice = 3 ;
K_rice = 10^( K_rice/10 );

tic

% for i = 1:1
%     [P_out(i)] = RUN_OPT_Protocol( K_rice )
% end
kappa0 = 1; % 0 1 2 3 4
N = 256;
 

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
%     K_rice = 10;
    K = 4; %用户个数
    B = 4; %全息超表面个数，每个全息超表面上配备有M个元素
    NUM_Position = 100; %6DMA空间位置
    delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [ 0 : ( My - 1 ) ].' - My + 1 )/2;
    Coor_Ele_init = [ kron( delta_x*dx , ones(My,1)  ) , kron( ones(Mx,1) , delta_y*dy ) , zeros(M,1)  ];
    % 标记横轴为x轴，纵轴为y轴
    M_x_r = Mx/2;      delta_x_r = [ (-M_x_r + 1)/2 : 1 : ( M_x_r - 1 )/2 ].' ;
    M_y_c = My/2;      delta_y_c = [ (-M_y_c + 1)/2 : 1 : ( M_y_c - 1 )/2 ].' ;
%     Coor_Ele_init_r = [ delta_x_r * (lambda)/2 , (-M_x_r + 1)/2 * (lambda/2) * ones(M_x_r,1) , zeros(M_x_r,1) ];    %第0行辐射元素的坐标
%     Coor_Ele_init_c = [ ( M_y_c - 1 )/2 * (lambda/2) * ones(M_y_c,1) , delta_y_c * (lambda/2) , zeros(M_y_c,1) ];    %第0列辐射元素的坐标

    Coor_Ele_init_r = [ delta_x * dx , (-Mx + 1)/2 * (lambda/2) * ones(Mx,1) , zeros(Mx,1) ];    %第0行辐射元素的坐标
    Coor_Ele_init_c = [ ( Mx - 1 )/2 * dx * ones(My,1) , delta_y * dy , zeros(My,1) ];    %第0列辐射元素的坐标

    delta_ele = 2^kappa0;   d = dx * delta_ele;
    Coor_Ele_init_r = Coor_Ele_init_r(1:delta_ele:end,:);
    Coor_Ele_init_c = Coor_Ele_init_c(1:delta_ele:end,:);
    M_x_r = size(Coor_Ele_init_r,1);
    M_y_c = size(Coor_Ele_init_c,1);

    Coor_Ele_init_x_d = [ delta_x_r * (lambda/2) , (-M_y_c + 1)/2 * (lambda/2) * ones(M_y_c,1) , zeros(M_x_r,1) ];
    Coor_Ele_init_x_u = [ delta_x_r * (lambda/2) , (M_y_c  - 1)/2 * (lambda/2) * ones(M_y_c,1) , zeros(M_x_r,1) ];
    Coor_Ele_init_y_l = [ (-M_x_r+1) / 2 * (lambda/2) * ones(M_x_r,1) , delta_y_c * (lambda/2) , zeros(M_x_r,1) ];
    Coor_Ele_init_y_r = [ ( M_x_r-1) / 2 * (lambda/2) * ones(M_x_r,1) , delta_y_c * (lambda/2) , zeros(M_x_r,1) ];

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
[ Coor_Ele , normal_vector , Rot_Matrix , Coor_q_init , Coor_Ele_r,Coor_Ele_c,Coor_Ele_x_d,Coor_Ele_x_u,Coor_Ele_y_l,Coor_Ele_y_r] = Orientation_Initial(B,Coor_Ele_init,Coor_Ele_init_r,Coor_Ele_init_c,0,Coor_Ele_init_x_d , Coor_Ele_init_x_u , Coor_Ele_init_y_l , Coor_Ele_init_y_r);
%% 针对初始6DMA位置，生成无线信道
[ h , e0 , e2 , theta0, phi0 , theta_scatterer0 , phi_scatterer0 , Iota , eta , Omega , h_r , h_c , h_sen , h_x_d , h_x_u , h_y_l , h_y_r , Omega0 ] = Channel_Generation_init( K,B,M,Mx,My,M_x_r , M_y_c , normal_vector,K_rice,Coor_Ele,Coor_Ele_r,Coor_Ele_c,lambda , Rot_Matrix , Coor_Ele_x_d,Coor_Ele_x_u,Coor_Ele_y_l,Coor_Ele_y_r);
%% 生成全息感知

% 所提感知算法
[est,est_L,est_error] = Sensing_Algorithm(K,B,Iota,Mx,My,M_x_r,M_y_c,h_r,h_c,Rot_Matrix,e0,K_rice,lambda,e2,sigma0,Omega,eta,N,d,kappa0,h_sen,h_x_d , h_x_u , h_y_l , h_y_r , Coor_Ele_init_x_d , Coor_Ele_init_x_u , Coor_Ele_init_y_l , Coor_Ele_init_y_r , Omega0);
% 
% est_error


% est{1} = [   -0.253906250000000
%    0.746093750000000
%   -0.615528823388373];
%   
% est{2} = [      0.863281250000000
%   -0.519531250000000
%                    0];
% est{3} = [    0.480468750000000
%    0.746093750000000
%   -0.460970602624913  ];
% est{4} = [     0.402343750000000
%    0.753906250000000
%    0.519369688224943];

% Benchmark 1
% [ek, error_ben1(:,ite_num)] = Benchmark_Sensing(B,M,K,V_F,h,lambda,eta0,delta_x,Mx,Coor_Ele,e0,N);

% est_direc_vector_rec：指向用户的方向向量
% normal_vector_rot：全息超表面的法向量
% Rot_Matrix：旋转矩阵
end


est{1}
[e0{1} e2{1,1} e2{1,2} e2{1,3} ]

est{2}
[e0{2} e2{2,1} e2{2,2} e2{2,3} ]

est{3}
[e0{3} e2{3,1} e2{3,2} e2{3,3} ]

est{4}
[e0{4} e2{4,1} e2{4,2} e2{4,3} ]






for k = 1:K
n1 = [ 0.02 ; 0.01 ; sqrt( 1 - 0.02^2 - 0.01^2 ) ]; %最大增益方向
n2 = est{k}; %用户方向
u = cross(n1,n2);   u = u / norm(u);
% 计算旋转角度（弧度单位）
alph = acos( n1.' * n2 );
% 计算矩阵U
U = [ 0    , -u(3) , u(2) ; ...
      u(3) , 0     , -u(1); ...
     -u(2), u(1)  , 0          ];
% 计算罗德格里斯旋转矩阵
Rot_Matrix{k} = cos(alph) * eye(3) + ( 1 - cos(alph) ) * ( u * u.' ) + sin(alph) * U;
Coor_Ele{k} = Coor_Ele_init * Rot_Matrix{k}.'; %旋转到第k个用户对应方向后的用户坐标
normal_vector_rot{k} = Rot_Matrix{k} * [0;0;1] ; % 旋转后的法向量

end

%% 生成B个离散空间点的波束增益向量
B1 = 40; %一共有64个空间离散点位
% e = [ error_theta , error_phi , error_varphi ].'; % 估计方向向量：e \in e*K
for k = 1:B, e(:,k) = est{k}; end



[ X , Y , Z ] = Orientation_uniformSpherePoints( B1 );
q = [ X.' , Y.' , Z.' ];
n1 = [ 0 ; 0 ; 1 ];

% 第k个用户支持的离散位置集合
C_set = cell(K,1);
for k = 1:K
    C_set{k} = [];
    for b0 = 1:B1
        if q(b0,:) * e(:,k) > 0
            C_set{k} = [ C_set{k} , b0 ];
        end
    end
end


for k = 1:K
    %注意，dx = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2
    Px{k}    = dx * Rot_Matrix{k} * [1;0;0];
    Py{k}    = dy * Rot_Matrix{k} * [0;1;0];
end

%%%%%%%%
%RHS对准第k0个用户，位于第b0个位置，在第k1个用户方向的导向矢量
for k0 = 1:K
    for b0 = 1:B1
        for k1 = 1:K
            a1 =  q(b0,:) +   ( Rot_Matrix{k0} * Coor_Ele_init.' ).';
            a_st_set{k0,b0,k1} = exp( 1j * 2*pi/lambda * ( a1 * e(:,k1) ) ) ;
%             a_st_set{k0,b0,k1} = exp( 1j * 2*pi/lambda * ( q(b0,:) * e(:,k1) +  kron( delta_x * ( Px{k0}.' * e(:,k1) ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{k0}.' * e(:,k1) ) )  ) ) ;
        end
    end
end
% 
% a1 =  q(b0,:) +   ( Rot_Matrix{k0} * Coor_Ele_init.' ).'
% a2 =  q(b0,:) + kron( delta_x * ( Px{k0}.' ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{k0}.'  ) )  
% a3 = [];
% for mx = -15.5:15.5
%     for my = -15.5:15.5
%         a3 = [ a3 ; q(b0,:) + [ mx *dx * Rot_Matrix{k0} * [1;0;0]  + my *dy *  Rot_Matrix{k0} * [0;1;0] ].'  ];
%     end
% end

























% 生成方向矩阵和距离矩阵
U_direc = zeros( B1 , K );   D_dis = zeros( B1 , B1 );

for i = 1:B1
    for j = 1:K
        U_direc( i , j ) = ( Rot_Matrix{j} * [0;0;1] ).' * q(i,:).'   ; %此处根据导向矢量求
    end
end
for i = 1:B1
    for j = 1:B1
        D_dis(i,j) = norm( q(i,:) - q(j,:) , 2 ) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 初始空间姿态
for b = 1:B
    s(:,b) = [ zeros(b-1,1) ; 1 ; zeros( B1 - (b-1) - 1 , 1 ) ];
end
g = eye(B,B);


%% Alternating Optimizing Method
% 制备b个平移向量的可行区间,搜索初始可行解

s_dic = zeros(B1,B1,B,B);  
for i1 = 1:B
    for i2 = 1:B
        if i1 == i2, s_dic(:,1,i1,i2) = ones(B1,1); else, s_dic(:,:,i1,i2) = ones(B1,B1); end
    end
end

for i1 = 1:B
    for i2 = 1:B
        for j1 = 1:B1
                s_temp1 = zeros(B1,1); s_temp1(j1) = 1;
                if s_temp1.' * U_direc * g(:,i1) <= 0, s_dic(j1,1,i1,i1) = 0 ; end
            for j2 = 1:B1
                s_temp2 = zeros(B1,1); s_temp2(j2) = 1;
                if s_temp1.' * D_dis * s_temp2 < dmin & i1 ~= i2, s_dic(j1,j2,i1,i2) = 0; end
                if (s_temp2 - s_temp1).' * U_direc * g(:,i1) >= 0 & i1 ~= i2, s_dic(j1,j2,i1,i2) = 0; end
            end
        end
    end
end

LoopUntil = 0;

Bset = 1:B1;

for m1 = C_set{1}
    for m2 = C_set{2}(C_set{2}~=m1)
        for m3 = C_set{3}(C_set{3}~=m1&C_set{3}~=m2)
            for m4 = C_set{4}(C_set{4}~=m1 & C_set{4}~=m2 & C_set{4} ~=m3)
                                if s_dic(m1,1,1,1)  + s_dic(m2,1,2,2) + s_dic(m3,1,3,3) + s_dic(m4,1,4,4) + ...
                                   s_dic(m1,m2,1,2) + s_dic(m1,m3,1,3) + s_dic(m1,m4,1,4) + ...
                                   s_dic(m2,m1,2,1) + s_dic(m2,m3,2,3) + s_dic(m2,m4,2,4) + ...
                                   s_dic(m3,m1,3,1) + s_dic(m3,m2,3,2) + s_dic(m3,m4,3,4) + ...
                                   s_dic(m4,m1,4,1) + s_dic(m4,m2,4,2) + s_dic(m4,m3,4,3)  ...
                                   == 16
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
s = zeros(B1,B);

if LoopUntil == 1
    s(m1,1) = 1;s(m2,2) = 1;s(m3,3) = 1;s(m4,4) = 1;
else 
    error('ERROR'); %必须存在可行点，才能进一步优化，否则不能进行优化
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 改成新的模型
Xtr = exp(1j*2*pi*randn(B,K));
coff   = diag( rand(K,1) ); %全息超表面波束系数
% 初始化波束增益，模拟波束
m_temp = zeros(M,B); %模拟波束增益


for k = 1:K
    for b = 1:B
        if normal_vector_rot{b}.' * est{k} > 0
            % a_st(:,k,b) =  exp( 1j * 2*pi/lambda * ( s(:,b).' * q * e(:,k) +  kron( delta_x * ( Px{b}.' * e(:,k) ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b}.' * e(:,k) ) )  ) ) ;
            a_st(:,k,b) = a_st_set{b,find(s(:,b)==1),k};
        else
            a_st(:,k,b)  = zeros(M,1);
        end
        m_temp(:,b)    = m_temp(:,b)  +  coff(k,b) * (  real( conj(V_F0) .* conj(a_st(:,k,b)) + 1 ) / 2 )  ;
    end
end

mm = m_temp;

for k = 1:K
%     Gain_Beam(k,1) = 0 ;
    for k1 = 1:K
    for b = 1:B
%         Gain_Beam(k,1) = Gain_Beam(k,1) + abs( a_st{k,b}.' * diag( m_temp(:,b) ) *  V_F * Xtr(k) )^2 ;
        gain_beam_temp(k1,b,k) = a_st(:,k,b).' * diag( m_temp(:,b) ) *  V_F * Xtr(k1,b);
    end
    end
    Gain_Beam(k,1) = sum( abs( gain_beam_temp(:,:,k) * ones(B,1) ).^2 ); 
end
% %     for b = 1:B
% %         Gain_Beam(b,1) = abs( a_st{b,b}.' * diag( m_temp(:,b) ) *  V_F * Xtr(b) )^2 ; 
% %     end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = 0.1 * rand(K,K);

for ITE_NUM = 1:5


Gain_Beam_Fair = min( Gain_Beam ) ; %公平波束增益
Gain_Beam_Fair_ITE(1) = Gain_Beam_Fair;     Gain_Beam_Fair_old =Gain_Beam_Fair;
time1 = 0;    time2 = 0;                        
%没有在每次迭代完更新导向矢量

for ite_num = 1:3
    
    Gain_Beam_Fair_ITE(ite_num) = Gain_Beam_Fair_old;


for b1 = 1 : B1
    
    Gain_Beam_Fair_temp1 = zeros(B,1);

    for b = 1:B
        if sum( C_set{b} == b1 ) >= 1
        % m_temp = mm; %mm_temp(:,b) = 0;
        m_temp = zeros(M,B);
        s_temp = s;
        s_temp(:,b) = zeros(B1,1); s_temp(b1,b) = 1;


        for k0 = 1:K
            for b0 = 1:B
                if normal_vector_rot{b0}.' * est{k0} > 0
                    % a_st_temp(:,k0,b0) =  exp( 1j * 2*pi/lambda * ( (s_temp(:,b0).') * q * e(:,k0) +  kron( delta_x * ( Px{b0}.' * e(:,k0) ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b0}.' * e(:,k0) ) )  ) ) ;            
                    a_st_temp(:,k0,b0) = a_st_set{ b0 , find(s_temp(:,b0)==1) , k0 };
                else
                    a_st_temp(:,k0,b0)  = zeros(M,1);
                end
                m_temp(:,b0)    = m_temp(:,b0)  +  coff(k0,b0) * (  real( conj(V_F0) .* conj(a_st_temp(:,k0,b0)) + 1 ) / 2 )  ;
            end
        end

        for k0 = 1:K


            mm1 = sparse( diag( m_temp(:,1) ) ); mm2 = sparse( diag( m_temp(:,2) ) ); mm3 = sparse( diag( m_temp(:,3) ) ); mm4 = sparse( diag( m_temp(:,4) ) ); 
            gain_beam_temp(:,:,k0) = [  a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(1,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(2,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(3,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(4,1) 
                                        a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(1,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(2,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(3,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(4,2) 
                                        a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(1,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(2,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(3,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(4,3) 
                                        a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(1,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(2,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(3,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(4,4) 
                                        ];

            Gain_Beam(k0,1) = sum( abs( gain_beam_temp(:,:,k0) * ones(B,1) ).^2 ); 

        end



     

        % m1 = find(s(:,1) == 1); m2 = find(s(:,2) == 1); m3 = find(s(:,3) == 1); m4 = find(s(:,4) == 1); m5 = find(s(:,5) == 1); m6 = find(s(:,6) == 1);
        m1 = find(s_temp(:,1) == 1); m2 = find(s_temp(:,2) == 1); m3 = find(s_temp(:,3) == 1); m4 = find(s_temp(:,4) == 1); 

        if                 s_dic(m1,m2,1,2)  + s_dic(m1,m3,1,3)  + s_dic(m1,m4,1,4)  + ...
                           s_dic(m2,m1,2,1)  + s_dic(m2,m3,2,3)  + s_dic(m2,m4,2,4)  + ...
                           s_dic(m3,m1,3,1)  + s_dic(m3,m2,3,2)  + s_dic(m3,m4,3,4)  + ...
                           s_dic(m4,m1,4,1)  + s_dic(m4,m2,4,2)  + s_dic(m4,m3,4,3)    ...
                             == 12

                            Gain_Beam_Fair_temp1(b,1) = min( Gain_Beam );
                           if Gain_Beam_Fair_temp1(b,1) > Gain_Beam_Fair_old
                              s = s_temp;
                              mm = m_temp;
                              Gain_Beam_Fair_old = Gain_Beam_Fair_temp1(b,1);
                              a_st = a_st_temp;
                           end

        
        end






    
        end
   
    
    end

end

end

% figure; plot( 1:length(Gain_Beam_Fair_ITE) , Gain_Beam_Fair_ITE , '-o' )

% 更新的结果似乎并没有用于空间姿态的设计上，还需要重新检查



%% 优化全息波束和数字波束
for b = 1:B
    mm(:,b) = real( conj(V_F0).' .* conj(a_st(:,b,b)).' + 1 )/2;
end



% 
% for b = 1:B
%     for k = 1:K
%         A_ST(k,:,b) = real( conj(V_F0).' .* conj(a_st(:,k,b)).' + 1 )/2 ;
%     end
% end

for ite_num = 1:3
% % 优化全息波束
% cvx_begin
% variable coff(K,B)
% variable P0(1)
% expressions mm(M,B) Gain_Beam(K,1) gain_beam_temp(K,B,K) gain_vec(K,K)
% 
% for b = 1:B, mm(:,b) = coff(:,b).' * A_ST(:,:,b); end
% 
% for k0 = 1:K
%     for b0 = 1:B
%         for k1 = 1:K
%             gain_beam_temp(k1,b0,k0) = a_st(:,k0,b0).' * diag(V_F * Xtr(b0,k1) ) * mm(:,b0);
%         end
%     end
% end
% 
% for k = 1:K
%     for k1 = 1:K
%         gain_vec(k1,k) = ones(1,B) * gain_beam_temp(k1,:,k)';
%     end
% end
% 
% for k = 1:K
%     Gain_Beam(k) = 2 * real( psi(:,k)' * gain_vec(:,k) ) - psi(:,k)' * psi(:,k);
% end
% 
% maximize P0
% subject to
% 
% Gain_Beam >= P0;
% 
% for b0 = 1:B
%     sum( coff(:,b0) ) <= 1;
% end
% 0<= coff <= 1;
%  
% cvx_end
% 
% for k = 1:K
%     psi(:,k) = gain_vec(:,k);
% end



%优化数字波束
% for b = 1:B, mm(:,b) = coff(:,b).' * A_ST(:,:,b); end

% 优化全息波束
cvx_begin
variable Xtr(B,K) complex
variable P0(1)
expressions Gain_Beam(K,1) gain_beam_temp(K,B,K) gain_vec(K,K)


for k0 = 1:K
    for b0 = 1:B
        for k1 = 1:K
            gain_beam_temp(k1,b0,k0) = a_st(:,k0,b0).' * diag(V_F * Xtr(b0,k1) ) * mm(:,b0);
        end
    end
end

for k = 1:K
    for k1 = 1:K
        gain_vec(k1,k) = ones(1,B) * gain_beam_temp(k1,:,k)';
    end
end

for k = 1:K
    Gain_Beam(k) = 2 * real( psi(:,k)' * gain_vec(:,k) ) - psi(:,k)' * psi(:,k);
end

maximize P0
subject to

Gain_Beam >= P0;

for b0 = 1:B
   sum( sum( pow_abs( Xtr , 2 ) ) ) <= Ptx;
end
 
cvx_end

for k = 1:K
    psi(:,k) = gain_vec(:,k);
end
% Beam_Gain_Ite(ite_num) = P0;

% fprintf('ite_num=%d, size(gain_beam_temp)=[%d %d %d], size(gain_vec)=[%d %d]\n', ...
%         ite_num, size(gain_beam_temp,1), size(gain_beam_temp,2), size(gain_beam_temp,3), ...
%         size(gain_vec,1), size(gain_vec,2));


end

Beam_Gain_Ite(ITE_NUM) = P0;






end

Pt = 40 ;
Ptx = 10^( (Pt - 30)/10 ); sigma0 = 1e-8;

% R00 = 18 ;





%% 生成信道
% Omega0 = Omega(1,1);

% 生成无线信道
H = zeros( B*M , K );
for k = 1:K
    for iota = 0:Iota
        for b = 1:B
            a1 = s(:,b).' * q + Coor_Ele{b};    %第b个全息超表面对准第b个用户
            if iota == 0
                    H( (b-1)*M+1 : b*M , k ) = H( (b-1)*M+1 : b*M , k ) + sig( normal_vector_rot{b}.' * e0{k} ) * sqrt( K_rice/(1+K_rice) ) * sqrt(Omega0(k,b)) * exp( 1j * 2*pi/lambda * ( a1 * e0{k}  ) ) ;
            else
                    H( (b-1)*M+1 : b*M , k ) = H( (b-1)*M+1 : b*M , k ) + sig( normal_vector_rot{b}.' * e2{k,iota} ) * sqrt(Omega0(k,b))* sqrt( 1/(1+K_rice) ) *  eta(k,iota) * exp( 1j * 2*pi/lambda * ( a1 * e2{k,iota} ) ) ;
            end
        end
    end
end



%% 不优化模拟波束，只优化数字波束
% 生成等效信道
for b = 1:B
    for k = 1:K
        H_eff(b,k) = H( (b-1)*M+1 : b*M , k ).' * diag( mm(:,b) ) * V_F; %等效于有B个射频链路，K个用户的MISO信道
    end
end

H_eff = H_eff /sqrt(sigma0);

% 初始化优化参数
norm_coff = 1;
psi1 = 0.05 * randn(K,K) + 0.05j * randn(K,K);   psi2 = 0.01 * rand(K,1) + 0.01j * rand(K,1) ;
rho  = 0.01 * rand(K,1); %功率分割因子,rho用于传能，(1-rho)用于通信
P_out_ite = [];
for ite_num = 1:20

%% 优化数能性能
cvx_begin
cvx_solver mosek
variable Xtr(B,K) complex
variable R0(1)
expressions R_temp(K,K) R_thro(K,1)

for k  = 1:K
    for k1 = 1:K
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
        end
    end
end

for k = 1:K
    R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + 1 ) ) / log(2);
end

maximize R0
subject to
R_thro >= R0;
sum( sum( pow_abs( Xtr , 2 ) ) ) <= Ptx;

cvx_end

psi_old = psi2;

for k = 1:K
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + 1 ) ;
end

P_out_ite(ite_num) = R0 * norm_coff;

norm_coff = norm(psi2,2);

if sum( isnan(psi2) ) >= 1
    psi2 = psi_old;
end


end

psi2_fixed = psi2;


num = 1;        rho = 0.1000000 * rand(K,1);

for  iiiii = 2

    R00 = (iiiii-1) * 1;

    [P_ave(iiiii)] = fun_run_R00( R00 , H_eff , rho , psi1 , psi2_fixed , sigma0 , Ptx , B , K );


end

P_out = P_ave(end)*sigma0


toc





















end