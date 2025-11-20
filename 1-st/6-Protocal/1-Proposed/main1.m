% 6DMA,使用全息感知方法
clc;clear;

Pt = 40;
Ptx = 10^( (Pt - 30)/10 ); sigma0 = 1e-7; R00 = 0;
kappa0 = 2; % 2 3 5 9 17

format long

% %% ****************************************初始化参数************************************************************************************************************************************************************
%Mx = 2; Dx = ( Mx - 1 ) * ( lambda / 4 ); Mx1 = Dx / dx + 1;
for nummmmmm = 1
    % Ptx = 10;
    fc = 30e9;   %载波频率
    c = 3e8;      %光速
    lambda = c/fc;   %波长
    dx = lambda/4;   %x轴方向天线间隔
    dy = lambda/4;  %y轴方向天线间隔
    Mx = 32; My = 32; M = Mx * My;
    Q  = 1; %每个全息超表面馈源个数
    Dx = Mx * dx;   Dy = My * dy;
    K_rice = 10^( -6 / 10 );
    K = 6; %用户个数
    B = 6; %全息超表面个数，每个全息超表面上配备有M个元素
    NUM_Position = 100; %6DMA空间位置
    delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [ 0 : ( My - 1 ) ].' - My + 1 )/2;
    % delta_x = [0:(Mx-1)].';   delta_y = [0:(My-1)].';
    % Coor_Ele_init = [ kron( ones(My,1) , delta_x*dx ) , kron( delta_y*dy , ones(Mx,1) ) , zeros(M,1)  ];  %辐射元素坐标
    Coor_Ele_init = [ kron( delta_x*dx , ones(My,1) ) , kron( ones(Mx,1) , delta_y*dy ) , zeros(M,1)  ];  %辐射元素坐标
    dmin = 0.1;
% figure(1); plot( Coor_Ele_init(:,1) , Coor_Ele_init(:,2) , 'o' )
% % 
%% 引入全息超表面传输模型
% Coor_Feed    = [ 0.002984305671304 , 0.001586131606604 , 0 ];   %馈源坐标
% Coor_Feed    = [ 0 , 3.838239570983593e-04 , 0 ];
% Coor_Feed = [ -0.0784 -0.0776 0 ];
% Coor_Feed = [ -0.023434 , 0.020677 0 ];
% Coor_Feed = [9.662364468102914e-04  0.002891155032136 0]; %最优馈源位置
% Coor_Feed = [0.010771   -0.020238 0 ; 0 0 0]; %随机馈源位置
%1个馈源
Coor_Feed = [ -0.0021   -0.0011 0 ];
% Coor_Feed = [  -0.000174821899166  -0.002230978260917  0 ];
% Coor_Feed = [0.002199321278631   0.000608108596301 0];
% Coor_Feed = [0.001862061370379  -0.002782682785307 0];
% Coor_Feed = [ -0.000479262690508   0.003326898762819  0 ];
%2个馈源
% Coor_Feed = [0.039979341818787   0.039967547200585 0 ; -0.000303984029534   0.018993235831849 0 ];
% Coor_Feed = [  -0.003059467310324  -0.000605847834945 0 ;  -0.003059467310324  -0.000605847834945 0  ];
% Coor_Feed = [ 0.006254478919187   0.031651406732093 0 ;    0.022812579049173  -0.022827586593938 0  ];
% Coor_Feed = [ -0.010223106017578   0.006912049032291 0 ;  -0.011969875142873   0.004567084538087 0 ];
%4个馈源
% Coor_Feed = [-0.001359734542402  -0.001447323745485 0 ; 0.002099163641678   0.000481926394073 0 ; -0.000147263827002  -0.001675380973294 0 ; -0.000081629756206   0.001767954310867 0 ];
% Coor_Feed = [-0.005373640461791  -0.011251513360779 0 ; -0.004434219428950   0.011376867810633 0 ; -0.005732142715214  -0.010916931655205 0 ; -0.002361468879039  -0.011842218780711 0 ];
% Coor_Feed = [ 0.067295693781457   0.067376813191020 0 ; 0.067329666144528   0.067333687477053 0 ;  0.067291583727011   0.067301061710797 0 ; 0.067294251284837   0.067388194672966 0 ];

for qq = 1:Q
    Dis_Feed2Ele(:,qq) = sqrt( ( Coor_Ele_init - Coor_Feed(qq,:) ).^2 * ones(3,1) );    %馈源到元素的距离
    V_F0(:,qq)         = exp( - 1j * 2*pi*sqrt(3)/lambda * Dis_Feed2Ele(:,qq) );    %馈源响应向量
end

eta0 = 8/3/M;
% eta0 = 1/M;

V_F = sqrt(eta0) * V_F0;

end



for  num = 1 : 1
%% 6DMA初始空间姿态
% Coor_Ele表示6DMA天线坐标，normal_vector是UPA的法向量
[ Coor_Ele , normal_vector , Rot_Matrix ] = Orientation_Initial(B,Coor_Ele_init,0);
%% 针对初始6DMA位置，生成无线信道
[ h , e0 , theta0, phi0 , theta_scatterer0 , phi_scatterer0 , Iota , eta , Omega ] = Channel_Generation_init( K,B,M,Mx,My,normal_vector,K_rice,Coor_Ele,lambda , Rot_Matrix );
%% 生成全息感知
% [ error_theta , error_phi , error_varphi , est_direc_vector_rec , Rot_Matrix , normal_vector_rot ] = Sensing_Holographic_HalfWave1(B,M,Mx,My,lambda,delta_x,delta_y,dx,dy,K,h,Rot_Matrix,e0,kappa0,normal_vector,Coor_Ele);
% [ error_theta , error_phi , error_varphi , est_direc_vector_rec , Rot_Matrix , normal_vector_rot ] = Sensing_Holographic_HalfWave1(B,M,Mx,My,lambda,delta_x,delta_y,dx,dy,K,h,Rot_Matrix,e0);
[ error_theta , error_phi , error_varphi , est_direc_vector_rec , Rot_Matrix , normal_vector_rot ] = Sensing_Holographic_HalfWave1(B,M,Mx,My,lambda,delta_x,delta_y,dx,dy,K,h,Rot_Matrix,e0,kappa0,normal_vector,Coor_Ele);
% est_direc_vector_rec：指向用户的方向向量
% normal_vector_rot：全息超表面的法向量
% Rot_Matrix：旋转矩阵
end

error_theta
error_phi
error_varphi


%% 生成B个离散空间点的波束增益向量
B1 = 50; %一共有64个空间离散点位
% e = [ error_theta , error_phi , error_varphi ].'; % 估计方向向量：e \in e*K
for k = 1:B, e(:,k) = est_direc_vector_rec{k}; end



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
            a_st_set{k0,b0,k1} = exp( 1j * 2*pi/lambda * ( q(b0,:) * e(:,k1) +  kron( delta_x * ( Px{k0}.' * e(:,k1) ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{k0}.' * e(:,k1) ) )  ) ) ;
        end
    end
end



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
tic
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
                for m5 = C_set{5}(C_set{5}~=m1 & C_set{5}~=m2 & C_set{5} ~=m3 & C_set{5}~=m4 )
                    for m6 = C_set{6}(C_set{6}~=m1 & C_set{6}~=m2 & C_set{6} ~=m3 & C_set{6}~=m4 & C_set{6}~=m5 )
                                if s_dic(m1,1,1,1)  + s_dic(m2,1,2,2) + s_dic(m3,1,3,3) + s_dic(m4,1,4,4) + s_dic(m5,1,5,5) + s_dic(m6,1,6,6) + ...
                                   s_dic(m1,m2,1,2) + s_dic(m1,m3,1,3) + s_dic(m1,m4,1,4) + s_dic(m1,m5,1,5) + s_dic(m1,m6,1,6)  + ...
                                   s_dic(m2,m1,2,1) + s_dic(m2,m3,2,3) + s_dic(m2,m4,2,4) + s_dic(m2,m5,2,5) + s_dic(m2,m6,2,6)  + ...
                                   s_dic(m3,m1,3,1) + s_dic(m3,m2,3,2) + s_dic(m3,m4,3,4) + s_dic(m3,m5,3,5) + s_dic(m3,m6,3,6)  + ...
                                   s_dic(m4,m1,4,1) + s_dic(m4,m2,4,2) + s_dic(m4,m3,4,3) + s_dic(m4,m5,4,5) + s_dic(m4,m6,4,6)  + ...
                                   s_dic(m5,m1,5,1) + s_dic(m5,m2,5,2) + s_dic(m5,m3,5,3) + s_dic(m5,m4,5,4) + s_dic(m5,m6,5,6)  + ...
                                   s_dic(m6,m1,6,1) + s_dic(m6,m2,6,2) + s_dic(m6,m3,6,3) + s_dic(m6,m4,6,4) + s_dic(m6,m5,6,5)   ...
                                   == 36
                                   LoopUntil = 1;
                                   break;
                        end
                        if LoopUntil == 1, break; end
                    end
                    if LoopUntil == 1, break; end
                end
                if LoopUntil == 1, break; end
            end
            if LoopUntil == 1, break; end
        end
        if LoopUntil == 1, break; end
    end
    if LoopUntil == 1, break; end
end
s = zeros(B1,B);

if LoopUntil == 1
    s(m1,1) = 1;s(m2,2) = 1;s(m3,3) = 1;s(m4,4) = 1;s(m5,5) = 1;s(m6,6) = 1;%s(m7,7) = 1;
else 
    error('ERROR'); %必须存在可行点，才能进一步优化，否则不能进行优化
end

toc





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 改成新的模型
Xtr = ones(Q,B,K);
coff   = ones(Q,K,B); %全息超表面波束系数
% 初始化波束增益，模拟波束
m_temp = zeros(M,B,Q); %模拟波束增益

tic
for k = 1:K
    for b = 1:B
        if normal_vector_rot{b}.' * est_direc_vector_rec{k} > 0
            % a_st(:,k,b) =  exp( 1j * 2*pi/lambda * ( s(:,b).' * q * e(:,k) +  kron( delta_x * ( Px{b}.' * e(:,k) ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b}.' * e(:,k) ) )  ) ) ;
            a_st(:,k,b) = a_st_set{b,find(s(:,b)==1),k};
        else
            a_st(:,k,b)  = zeros(M,1);
        end
        for qq = 1:Q
            m_temp(:,b)    = m_temp(:,b)  +  coff(qq,k,b) * (  real( conj(V_F0(:,qq)) .* conj(a_st(:,k,b)) + 1 ) / 2 )  ;
        end
    end
end
toc

mm = m_temp;

for k = 1:K
%     Gain_Beam(k,1) = 0 ;
    for k1 = 1:K
    for b = 1:B
%         Gain_Beam(k,1) = Gain_Beam(k,1) + abs( a_st{k,b}.' * diag( m_temp(:,b) ) *  V_F * Xtr(k) )^2 ;
        %% 波束增益
        gain_beam_temp(k1,b,k) = a_st(:,k,b).' * diag( m_temp(:,b) ) *  V_F * Xtr(:,b,k1); %在b个超表面上，第k1个用户的符号在第k个用户方向处的波束增益
    end
    end
    Gain_Beam(k,1) = sum( abs( gain_beam_temp(:,:,k) * ones(B,1) ).^2 ); 
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = 0.1 * rand(K,K);
for ITE_NUM = 1:2


Gain_Beam_Fair = min( Gain_Beam ) ; %公平波束增益
Gain_Beam_Fair_ITE(1) = Gain_Beam_Fair;     Gain_Beam_Fair_old =Gain_Beam_Fair;
time1 = 0;    time2 = 0;                        
% 没有在每次迭代完更新导向矢量

for ite_num = 1:2
    
    Gain_Beam_Fair_ITE(ite_num) = Gain_Beam_Fair_old;


for b1 = 1 : B1
    
    Gain_Beam_Fair_temp1 = zeros(B,1);

    for b = 1:B
        if sum( C_set{b} == b1 ) >= 1
 tic
        % m_temp = mm; %mm_temp(:,b) = 0;
        m_temp = zeros(M,B);
        s_temp = s;
        s_temp(:,b) = zeros(B1,1); s_temp(b1,b) = 1;
toc
tic
        for k0 = 1:K
            for b0 = 1:B
                if normal_vector_rot{b0}.' * est_direc_vector_rec{k0} > 0
                    a_st_temp(:,k0,b0) = a_st_set{ b0 , find(s_temp(:,b0)==1) , k0 };
                else
                    a_st_temp(:,k0,b0)  = zeros(M,1);
                end
                for qq = 1:Q
                    m_temp(:,b0)    = m_temp(:,b0)  +  coff(qq,k0,b0) * (  real( conj(V_F0(:,qq)) .* conj(a_st_temp(:,k0,b0)) + 1 ) / 2 )  ;
                end
            end
        end
toc
        for k0 = 1:K


            %% 波束增益
            mm1 = sparse( diag( m_temp(:,1) ) ); mm2 = sparse( diag( m_temp(:,2) ) ); mm3 = sparse( diag( m_temp(:,3) ) ); mm4 = sparse( diag( m_temp(:,4) ) ); mm5 = sparse( diag( m_temp(:,5) ) ); mm6 = sparse( diag( m_temp(:,6) ) ); %mm7 = sparse( diag( m_temp(:,7) ) ); 
            gain_beam_temp(:,:,k0) = [  a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(:,1,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(:,2,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(:,3,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(:,4,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(:,5,1) , a_st_temp(:,k0,1).' * mm1 *  V_F * Xtr(:,6,1)  
                                        a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(:,1,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(:,2,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(:,3,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(:,4,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(:,5,2) , a_st_temp(:,k0,2).' * mm2 *  V_F * Xtr(:,6,2) 
                                        a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(:,1,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(:,2,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(:,3,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(:,4,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(:,5,3) , a_st_temp(:,k0,3).' * mm3 *  V_F * Xtr(:,6,3) 
                                        a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(:,1,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(:,2,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(:,3,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(:,4,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(:,5,4) , a_st_temp(:,k0,4).' * mm4 *  V_F * Xtr(:,6,4) 
                                        a_st_temp(:,k0,5).' * mm5 *  V_F * Xtr(:,1,5) , a_st_temp(:,k0,5).' * mm5 *  V_F * Xtr(:,2,5) , a_st_temp(:,k0,5).' * mm5 *  V_F * Xtr(:,3,5) , a_st_temp(:,k0,5).' * mm5 *  V_F * Xtr(:,4,5) , a_st_temp(:,k0,5).' * mm5 *  V_F * Xtr(:,5,5) , a_st_temp(:,k0,5).' * mm5 *  V_F * Xtr(:,6,5) 
                                        a_st_temp(:,k0,6).' * mm6 *  V_F * Xtr(:,1,6) , a_st_temp(:,k0,6).' * mm6 *  V_F * Xtr(:,2,6) , a_st_temp(:,k0,6).' * mm6 *  V_F * Xtr(:,3,6) , a_st_temp(:,k0,6).' * mm6 *  V_F * Xtr(:,4,6) , a_st_temp(:,k0,6).' * mm6 *  V_F * Xtr(:,5,6) , a_st_temp(:,k0,6).' * mm6 *  V_F * Xtr(:,6,6) 
                                        ];

            Gain_Beam(k0,1) = sum( abs( gain_beam_temp(:,:,k0) * ones(B,1) ).^2 ); 

        end


disp(toc)
     
tic
        % m1 = find(s(:,1) == 1); m2 = find(s(:,2) == 1); m3 = find(s(:,3) == 1); m4 = find(s(:,4) == 1); m5 = find(s(:,5) == 1); m6 = find(s(:,6) == 1);
        m1 = find(s_temp(:,1) == 1); m2 = find(s_temp(:,2) == 1); m3 = find(s_temp(:,3) == 1); m4 = find(s_temp(:,4) == 1); m5 = find(s_temp(:,5) == 1); m6 = find(s_temp(:,6) == 1); %m7 = find(s_temp(:,7) == 1);

        if                 s_dic(m1,m2,1,2)  + s_dic(m1,m3,1,3)  + s_dic(m1,m4,1,4)  + s_dic(m1,m5,1,5)  + s_dic(m1,m6,1,6)  + ...
                           s_dic(m2,m1,2,1)  + s_dic(m2,m3,2,3)  + s_dic(m2,m4,2,4)  + s_dic(m2,m5,2,5)  + s_dic(m2,m6,2,6)  + ...
                           s_dic(m3,m1,3,1)  + s_dic(m3,m2,3,2)  + s_dic(m3,m4,3,4)  + s_dic(m3,m5,3,5)  + s_dic(m3,m6,3,6)  + ...
                           s_dic(m4,m1,4,1)  + s_dic(m4,m2,4,2)  + s_dic(m4,m3,4,3)  + s_dic(m4,m5,4,5)  + s_dic(m4,m6,4,6)  + ...
                           s_dic(m5,m1,5,1)  + s_dic(m5,m2,5,2)  + s_dic(m5,m3,5,3)  + s_dic(m5,m4,5,4)  + s_dic(m5,m6,5,6)  + ...
                           s_dic(m6,m1,6,1)  + s_dic(m6,m2,6,2)  + s_dic(m6,m3,6,3)  + s_dic(m6,m4,6,4)  + s_dic(m6,m5,6,5)    == 30

                            Gain_Beam_Fair_temp1(b,1) = min( Gain_Beam );
                           if Gain_Beam_Fair_temp1(b,1) > Gain_Beam_Fair_old
                              s = s_temp;
                              mm = m_temp;
                              Gain_Beam_Fair_old = Gain_Beam_Fair_temp1(b,1);
                              a_st = a_st_temp;
                           end

        
        end



toc
time2 = time2 + toc;


    
        end
   
    
    end

end

end

% figure; plot( 1:length(Gain_Beam_Fair_ITE) , Gain_Beam_Fair_ITE , '-o' )

% 更新的结果似乎并没有用于空间姿态的设计上，还需要重新检查



%% 优化全息波束和数字波束
A_st = zeros( M , Q , K , B );

for b = 1:B
    for k = 1:K
        for qq = 1:Q
            A_ST(:,qq,k,b) = real( V_F0(:,qq).' .* a_st(:,k,b).' + 1 ) / 2;
        end
    end
end


% for b = 1:B
%     for k = 1:K
%         for qq = 1:Q
%             A_ST(:,k,b,qq) = real( V_F0(:,qq).' .* a_st(:,k,b).' + 1 ) / 2 ;
%         end
%     end
% end

for ite_num = 1:6
% 优化全息波束
cvx_begin
variable coff(Q,K,B)
variable P0(1)
expressions mm_temp1(M,K,B) mm(M,B) Gain_Beam(K,1) gain_beam_temp(K,B,K) gain_vec(K,K)

for k = 1:K
for b = 1:B
%     for qq = 1:Q
        mm_temp1(:,k,b) =  A_ST(:,:,k,b) * coff(:,k,b);
%         mm_temp1(:,k,b,qq) = coff(k,b,qq) * A_ST(:,k,b,qq); 
%     end
end
end

for b = 1:B
    mm(:,b) = mm_temp1(:,:,b) * ones(K,1);
end

%% 波束增益
for k0 = 1:K
    for b0 = 1:B
        for k1 = 1:K
            gain_beam_temp(k1,b0,k0) = a_st(:,k0,b0).' * diag(V_F * Xtr(:,b0,k1) ) * mm(:,b0);
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
    sum( sum( coff(:,:,b0) ) ) == 1;
end

0<= coff <= 1;
 
cvx_end

for k = 1:K
    psi(:,k) = gain_vec(:,k);
end
% Beam_Gain_Ite(ite_num) = P0;


%优化数字波束

% 优化全息波束
cvx_begin
variable Xtr(Q,B,K) complex
variable P0(1)
expressions Gain_Beam(K,1) gain_beam_temp(K,B,K) gain_vec(K,K)


for k0 = 1:K
    for b0 = 1:B
        for k1 = 1:K
            gain_beam_temp(k1,b0,k0) = a_st(:,k0,b0).' * diag(V_F * Xtr(:,b0,k1) ) * mm(:,b0);
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
   sum( sum( sum( pow_abs( Xtr , 2 ) ) ) ) <= Ptx;
end
 
cvx_end

for k = 1:K
    psi(:,k) = gain_vec(:,k);
end
% Beam_Gain_Ite(ite_num) = P0;
end

Beam_Gain_Ite(ITE_NUM) = P0;

end

% figure; plot( log(Beam_Gain_Ite) , '-o' )





fprintf('\n');