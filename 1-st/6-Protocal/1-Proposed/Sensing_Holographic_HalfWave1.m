function [ error_theta , error_phi , error_varphi , est_direc_vector_rec , Rot_Matrix , normal_vector_rot ] = Sensing_Holographic_HalfWave1(B,M,Mx,My,lambda,delta_x,delta_y,dx,dy,K,h0,Rot_Matrix,e0,kappa0,normal_vector,Coor_Ele)
%以半波长为间隔，生成全息图
% kappa0 = 33;
Pt_UL = 5e-4;
I1 = sparse(0,0);
for num = 1:(Mx*My/4)
    I1 = blkdiag( I1 , [1,0] );
end

I2 = sparse(0,0);
for num = 1:Mx/2
    I2 = blkdiag( I2 , [ eye(My) , zeros(My,My) ] );
end
I = I1 * I2;

% I0 = kron( repmat([1;0],Mx/2,1) , repmat([1;0],My/2,1).' );
% I0 = ones(Mx,Mx);

Mx0 = 1:(kappa0-1):Mx; My0 = 1:(kappa0-1):My;

for k = 1:K
    for b = 1:B
        h1 = reshape( h0{k,b} , Mx , My );
%         h_temp = I * h0{k,b};
%         h{k,b} = sqrt(Pt_UL) * reshape( h_temp , Mx/2 , My/2 ).' ;
        h{k,b} = sqrt(Pt_UL) * h1(Mx0,My0);
    end
end

%% 生成全息感知
%生成参考信号
for k = 1:K
% for b = 1:B, x_ref{k,b} = sqrt(1e-2) * I0 .* exp( 1j * rand(Mx,My)*2*pi ); end
for b = 1:B, x_ref{k,b} = sqrt(1e-2) * exp( 1j * rand( length(Mx0) , length(My0) )*2*pi ); end
end

%全息图，存储在计算机中
for k = 1:K
    for b = 1:B
        n{k,b} = sqrt(1e-12) * ( randn(length(Mx0),length(My0)) + 1j * randn(length(Mx0),length(My0)) );
        H_holographic{k,b} =  2 * real( conj(h{k,b}) .* x_ref{k,b} + n{k,b} )  ; 
    end
end 

for b = 1:B
    b_steering{b} = @(e)  exp( -1j * 2*pi/lambda * Coor_Ele{b} * e ) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%使用穷举搜索，找到感知角度。
angle_set_0 = -1:0.01:1;    angle_set_1 = -1:0.02:1;
for k = 1 : K
    for b = 1 : B
        T{k,b}   = x_ref{k,b} .* H_holographic{k,b};
        num1 = 1;
        for d_x0 = angle_set_0
            num2 = 1;
            for d_y0 = angle_set_1
                if d_x0^2 + d_y0^2 <=1

                    e_erg0 = [ d_x0 ; d_y0 ; sqrt( max( 1 - d_x0^2 - d_y0^2 , 0 ) ) ];
                    if normal_vector{b}.' * e_erg0 >=0
                        % 第三项为正数
                        H_rec0{b,k}(num1,num2) = abs( b_steering{b}(e_erg0).' * vec( T{k,b} ) );
                    end

                    e_erg1 = [ d_x0 ; d_y0 ; -sqrt( max( 1 - d_x0^2 - d_y0^2 , 0 ) ) ];
                    if normal_vector{b}.' * e_erg1 >=0
                        % 第三项为负数
                        H_rec1{b,k}(num1,num2) = abs( b_steering{b}(e_erg1).' * vec( T{k,b} ) );
                    end

                end

                num2 = num2 + 1;

            end
            num1 = num1 + 1;
        end

    end
end

for k = 1 : K
    for b = 1 : B
        H_max{k}(b,:) = [ max(max( H_rec0{b,k} )) , max(max( H_rec1{b,k} )) ];
    end
end

% 搜索增益最大得超表面
for k = 1:K
    [maxGain{k}, idx{k}] = max( H_max{k}(:) );
    [maxRow{k}, maxCol{k}] = ind2sub(size(H_max{k}), idx{k});
    if maxCol{k} == 1
        H_angle{k} = H_rec0{maxRow{k},k};
        flag(k,1) = 1;
    elseif maxCol{k} == 2
        H_angle{k} = H_rec1{maxRow{k},k};
        flag(k,1) = 2;
    end
end

% 搜索每个用户的最大角度
for k = 1:K
    [maxGain0{k}, idx0{k}] = max( H_angle{k}(:) );
    [maxRow0{k}, maxCol0{k}] = ind2sub(size(H_angle{k}), idx0{k});
    theta_est(k) = angle_set_0(maxRow0{k}); phi_est(k) = angle_set_1(maxCol0{k});
    if flag(k,1) == 1
        varphi_est(k) =  sqrt( max( 1 - theta_est(k)^2 - phi_est(k)^2 , 0 ) );
    elseif flag(k,1) == 2
        varphi_est(k) = -sqrt( max( 1 - theta_est(k)^2 - phi_est(k)^2 , 0 ) );
    end

    est_direc_vector_rec{k} = [ theta_est(k) ; phi_est(k) ; varphi_est(k) ];
    % % est_direc_vector_rec{k} = e0{k};
end






for k = 1:K
    error_theta(k,1) = est_direc_vector_rec{k}(1,1) - e0{k}(1);
    error_phi(k,1)   = est_direc_vector_rec{k}(2,1) - e0{k}(2);
    error_varphi(k,1)= est_direc_vector_rec{k}(3,1) - e0{k}(3);
end


















% 旋转最大增益角
for k = 1:K
% n1 = [ 0.02 ; 0.01 ; sqrt( 1 - 0.02^2 - 0.01^2 ) ]; %最大增益方向
% n1 = [ -0.278445348362427 ; 0.236305977207617 ; 0.930928393116936 ]; %优化馈源位置，最佳方向
% n1 = [[ -0.285894092386310 ; 0.236512162150458 ; 0.928615402141017 ]]; %随机馈源位置，最佳方向
% n1 = [0;0;1];
% n1 = [  -0.477630009460290
%    0.843521030655574
%    0.245645771192427];
n1 = [   -0.3799
    0.8718
    0.3092];

n2 = est_direc_vector_rec{k}; %用户方向
u = cross(n1,n2);   u = u / norm(u);
% 计算旋转角度（弧度单位）
alph = acos( n1.' * n2 );
% 计算矩阵U
U = [ 0    , -u(3) , u(2) ; ...
      u(3) , 0     , -u(1); ...
     -u(2), u(1)  , 0          ];
% 计算罗德格里斯旋转矩阵
Rot_Matrix{k} = cos(alph) * eye(3) + ( 1 - cos(alph) ) * ( u * u.' ) + sin(alph) * U;

normal_vector_rot{k} = Rot_Matrix{k} * [0;0;1] ; % 旋转后的法向量

end





end


