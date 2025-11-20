function [theta_error,phi_error] = Run_Sensing(aaaaaaa)
ttttttttttttt=aaaaaaa;
% 在本函数中，所有全息超表面元素都配备一个功率计
%% ****************************************初始化参数************************************************************************************************************************************************************
fc = 30e9;   %载波频率
c = 3e8;      %光速
lambda = c/fc;   %波长；
dx = lambda/4;   %x轴方向天线间隔
dy = lambda/4;  %y轴方向天线间隔
Mx = 64; My = 64; M = Mx * My;
Dx = Mx * dx;   Dy = My * dy;
K_rice = 2;
K = 1; %用户个数
L = 6; %散射体个数
delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [0:( My-1)].' - My + 1 )/2;
Coor_Ele = [ kron( ones(My,1) , delta_x*dx ) , kron( delta_y*dy , ones(Mx,1) ) , zeros(M,1)  ];  %辐射元素坐标
% % % figure; plot( Coor_Ele(:,1) , Coor_Ele(:,2) , 'o' )
% 



% %% ****************************************生成无线信道************************************************************************************************************************************************************
% %接收机坐标
Dis_min_x = 0.5 * sqrt( Dx^3/lambda );
% 求用户的余弦归一化的水平角和俯仰角，二者取值范围均为[-1,1]
for k = 1:K
    while true
        Coor_Rec(k,:) = [ (rand(1)-0.5)*10 , (rand(1)-0.5)*10 , (rand(1)-0.5)*4 ];
        if norm(Coor_Rec(k,:),2) > Dis_min_x, break; end
    end
    theta_Rec(k) = Coor_Rec(k,1) / norm(Coor_Rec(k,:),2);
    phi_Rec(k) = Coor_Rec(k,3) / norm(Coor_Rec(k,:),2);
end
% 
% %散射体坐标
% 求散射体的余弦归一化的水平角和俯仰角，二者取值范围均为[-1,1]
for k = 1:K
for l = 1:L
    % while true
        Coor_Sca(l,:,k) = [ (rand(1)-0.5)*10 , (rand(1)-0.5)*10 , (rand(1)-0.5)*4 ];
        % if norm(Coor_Sca(l,:),2) > 3, break; end
    % end
    theta_Sca(l,k) = Coor_Sca(1) / norm(Coor_Sca(l,:),2);
    phi_Sca(l,k) = Coor_Sca(3) / norm(Coor_Sca(l,:),2);
end
end
% 散射体增益
alpha_Sca = 1/sqrt(2*L) * randn(L,K) + 1j/sqrt(2*L) * randn(L,K);
% 定义导向矢量：
b = @(theta,phi) 1/sqrt(M) * kron( exp( 1j * 2*pi/lambda * delta_x * dx * theta ) , exp( 1j * 2*pi/lambda * delta_y * dy * phi  ) );
%第k个用户的无线信道增益可以表示为
for k = 1:K
    h{k} = 0 ;
    for l = 0:L
        if l ==0
            h{k} = h{k} + sqrt(M) * sqrt(K_rice/(1+K_rice)) * b( theta_Rec(k) , phi_Rec(k) );
        else
            h{k} = h{k} + sqrt(M) * alpha_Sca(l,k) * sqrt(1/(1+K_rice)) * b( theta_Sca(l,k) , phi_Sca(l,k) );
        end
    end
    h{k} = 1e-6 * h{k};
end

%% 生成全息图
x_ref = 1e-3 * exp( 1j * rand(M,1) ); %全息干扰值
for k = 1:K, H_holographic{k} = 2 * real( conj(h{k}) .* x_ref ); end %全息图，存储在计算机中
theta_grid = ( 2*[0:Mx-1] - Mx +1 ) / Mx; phi_grid = ( 2*[0:My-1] - My +1 ) / My;

for k = 1:K
for num1 = 1:Mx
    for num2 = 1:My
        H_recover(num1,num2) = sqrt(M) * b( theta_grid(num1) , phi_grid(num2) )' .* x_ref.' * H_holographic{k} ;
    end
end
    [maxColumnValues, rowIndices] = max(abs(H_recover));  % 每列的最大值及其行索引
    [maxValue, colIndex] = max(maxColumnValues);  % 找出最大值及其列索引
    rowIndex = rowIndices(colIndex);  % 最大值的行索引
    theta_est(k) = theta_grid(rowIndex);
    phi_est(k)   = phi_grid(colIndex)  ;
end
% 误差记录
theta_error = theta_Rec - theta_est;
phi_error   = phi_Rec   - phi_est  ;
end