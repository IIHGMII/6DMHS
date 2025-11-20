function [Coor_Ele,normal_vector,Rot_Matrix,q,Coor_Ele_x,Coor_Ele_y] = Orientation_Initial(B,Coor_Ele_init,Coor_Ele_init_r,Coor_Ele_init_c,flag)

%此处，在空间中均匀生成B个点数，让全息超表面均匀的在自由空间中排列
[X , Y , Z] = Orientation_uniformSpherePoints(B);

% 初始法向量
n1 = [0;0;1];   am1 = [ 0.02 ; 0.01 ; sqrt( 1 - 0.02^2 -0.01^2 ) ];

for b = 1:B
    % 平移变量
    q(b,:) = [X(b) , Y(b) , Z(b) ];
    % 目标单位的法向量,并归一化
    n2 = [X(b) ; Y(b) ; Z(b) ]; n2 = n2 / norm(n2);
    if norm(n2 - n1) <= 1e-3
        Coor_Ele{b} = q(b,:) + Coor_Ele_init;   %此时UPA的法向量为[0;0;1]
        normal_vector{b} = q(b,:).';
        Rot_Matrix{b} = eye(3);
        am2{b} = Rot_Matrix{b} * am1;
    elseif norm(n2 + n1) <= 1e-3
        Coor_Ele{b} = q(b,:) + Coor_Ele_init;   %此时UPA的法向量为[0;0;-1]
        normal_vector{b} = q(b,:).'; 
        Rot_Matrix{b} = [ 1,0,0;0,1,0;0,0,-1 ];
        am2{b} = Rot_Matrix{b} * am1;
    else
    % 计算旋转轴(叉乘),并归一化
    u = cross(n1,n2);   u = u / norm(u);
    % 计算旋转角度（弧度单位）
    alph = acos( n1.' * n2 );
    % 计算矩阵U
    U = [ 0    , -u(3) , u(2) ; ...
          u(3) , 0     , -u(1); ...
          -u(2), u(1)  , 0          ];
    % 计算罗德格里斯旋转矩阵
    Rot_Matrix{b} = cos(alph) * eye(3) + ( 1 - cos(alph) ) * ( u * u.' ) + sin(alph) * U;
    Coor_Ele{b} = q(b,:) + Coor_Ele_init * Rot_Matrix{b}.';
    normal_vector{b} = q(b,:).'; 
    am2{b} = Rot_Matrix{b} * am1;
    end
end


for b = 1:B
    % 目标单位的法向量,并归一化
    n2 = [X(b) ; Y(b) ; Z(b) ]; n2 = n2 / norm(n2);
    
    if norm(n2-n1) <= 1e-3
        Coor_Ele_x{b} = q(b,:) + Coor_Ele_init_r;
        normal_vector{b} = q(b,:).';
        Rot_Matrix{b} = eye(3);
        am2{b} = Rot_Matrix{b} * am1;
    elseif norm(n2+n1) <= 1e-3
        Coor_Ele_x{b} = q(b,:) + Coor_Ele_init_r;
        normal_vector{b} = q(b,:).';
        Rot_Matrix{b} = [1,0,0;0,1,0;0,0,-1];
        am2{b} = Rot_Matrix{b} * am1;
    else
    
    % 计算旋转轴(叉乘),并归一化
    u = cross(n1,n2);   u = u / norm(u);
    % 计算旋转角度（弧度单位）
    alph = acos( n1.' * n2 );
    % 计算矩阵U
    U = [ 0    , -u(3) , u(2) ; ...
          u(3) , 0     , -u(1); ...
          -u(2), u(1)  , 0          ];
    % 计算罗德格里斯旋转矩阵
    Rot_Matrix{b} = cos(alph) * eye(3) + ( 1 - cos(alph) ) * ( u * u.' ) + sin(alph) * U;
    Coor_Ele_x{b} = q(b,:) + Coor_Ele_init_r * Rot_Matrix{b}.';
    normal_vector{b} = q(b,:).'; 
    am2{b} = Rot_Matrix{b} * am1;
    
    end
end

% 验证旋转后的x-轴方向的感知元素坐标表达式是否正确
% for b = 1:B
% Coor_Ele_x_Check{b} =  q(b,:).' + Rot_Matrix{b} * Coor_Ele_init_r.';
% Coor_Ele_x_Check{b} - Coor_Ele_x{b}.'
% end
% 经检验，旋转后的坐标正确



for b = 1:B
    % 目标单位的法向量,并归一化
    n2 = [X(b) ; Y(b) ; Z(b) ]; n2 = n2 / norm(n2);
    
    if norm(n2-n1) <= 1e-3
        Coor_Ele_y{b} = q(b,:) + Coor_Ele_init_c;
%         normal_vector{b} = q(b,:).';
%         Rot_Matrix{b} = eye(3);
%         am2{b} = Rot_Matrix{b} * am1;
    elseif norm(n2+n1) <= 1e-3
        Coor_Ele_y{b} = q(b,:) + Coor_Ele_init_c;
%         normal_vector{b} = q(b,:).';
%         Rot_Matrix{b} = [1,0,0;0,1,0;0,0,-1];
%         am2{b} = Rot_Matrix{b} * am1;
    else
    
    % 计算旋转轴(叉乘),并归一化
    u = cross(n1,n2);   u = u / norm(u);
    % 计算旋转角度（弧度单位）
    alph = acos( n1.' * n2 );
    % 计算矩阵U
    U = [ 0    , -u(3) , u(2) ; ...
          u(3) , 0     , -u(1); ...
          -u(2), u(1)  , 0          ];
    % 计算罗德格里斯旋转矩阵
%     Rot_Matrix{b} = cos(alph) * eye(3) + ( 1 - cos(alph) ) * ( u * u.' ) + sin(alph) * U;
    Coor_Ele_y{b} = q(b,:) + Coor_Ele_init_c * Rot_Matrix{b}.';
%     normal_vector{b} = q(b,:).'; 
%     am2{b} = Rot_Matrix{b} * am1;
    end
end

% 验证旋转后的x-轴方向的感知元素坐标表达式是否正确
% for b = 1:B
% Coor_Ele_y_Check{b} =  q(b,:).' + Rot_Matrix{b} * Coor_Ele_init_c.';
% Coor_Ele_y_Check{b} - Coor_Ele_y{b}.'
% end
% 经检验，结果正确

if flag == 1

% % 展示6DMA空间图

figure(2);

for b =1:B
    scatter3(Coor_Ele{b}(:,1), Coor_Ele{b}(:,2), Coor_Ele{b}(:,3), 'filled')
    hold on;
end

end





end