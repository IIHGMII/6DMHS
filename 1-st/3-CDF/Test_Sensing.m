clc;clear;



Pt = 40;
Ptx = 10^( (Pt - 30)/10 ); sigma0 = 1e-8; R00 = 1;


format long

% %% ****************************************初始化参数************************************************************************************************************************************************************
for nummmmmm = 1
    % Ptx = 10;
    fc = 30e9;   %载波频率
    c = 3e8;      %光速
    lambda = c/fc;   %波长；
    dx = lambda/4;   %x轴方向天线间隔
    dy = lambda/4;  %y轴方向天线间隔
    Mx = 32; My = 32; M = Mx * My;
    Q  = 1; %每个全息超表面馈源个数
    Dx = Mx * dx;   Dy = My * dy;
    K_rice = 5;
    K = 7; %用户个数
    B = 7; %全息超表面个数，每个全息超表面上配备有M个元素
    NUM_Position = 100; %6DMA空间位置
    delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [ 0 : ( My - 1 ) ].' - My + 1 )/2;
    % delta_x = [0:(Mx-1)].';   delta_y = [0:(My-1)].';
    % Coor_Ele_init = [ kron( ones(My,1) , delta_x*dx ) , kron( delta_y*dy , ones(Mx,1) ) , zeros(M,1)  ];  %辐射元素坐标
    Coor_Ele_init = [ kron( delta_x*dx , ones(My,1) ) , kron( ones(Mx,1) , delta_y*dy ) , zeros(M,1)  ];  %辐射元素坐标
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



tic

parfor ite_num = 1:120

for  num = 1 : 1
%% 6DMA初始空间姿态
% Coor_Ele表示6DMA天线坐标，normal_vector是UPA的法向量
[ Coor_Ele , normal_vector , Rot_Matrix ] = Orientation_Initial(B,Coor_Ele_init,0);
%% 针对初始6DMA位置，生成无线信道
[ h , e0 , theta0, phi0 , theta_scatterer0 , phi_scatterer0 , Iota , eta , Omega ] = Channel_Generation_init( K,B,M,Mx,My,normal_vector,K_rice,Coor_Ele,lambda , Rot_Matrix );
% h{k,b}

end

[ek, error_ben1(:,ite_num),error_theta,error_phi,error_varphi] = Ben_Sensing1(B,M,K,V_F,h,lambda,eta0,delta_x,Mx,Coor_Ele,e0);

E(:,:,num) = [error_theta , error_phi , error_varphi];

end

% sum(sum( error_ben1.^2 )) / numel( error_ben1 )

toc


MSE = (sum( E.^2 , 'all' )) / NUM / 3 / K ;
disp( MSE );
fprintf('\n');
disp('Running Over!');
fprintf('\n');
% []

E_vec = reshape(E,3*B*NUM,1);
f0=0;
f=[];
for num = 1:500
    f(1,num) = length(  find( abs(E_vec) <= ( num * 0.001 )  ) );
%     f0 = f(1,num);
end
f = [0,f]/(3*B*NUM);






