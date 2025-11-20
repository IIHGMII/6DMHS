function [ error_theta , error_phi , theta_est_global , phi_est_global , est_normal_vector_global_out , Rot_Matrix ] = Sensing_Holographic_HalfWave(B,M,Mx,My,lambda,delta_x,delta_y,dx,dy,K,h,Rot_Matrix,e0)
%以半波长为间隔，生成全息图
Pt_UL = 1e-3;
I1 = sparse(0,0);
for num = 1:(Mx*My/4)
    I1 = blkdiag( I1 , [1,0] );
end

I2 = sparse(0,0);
for num = 1:Mx/2
    I2 = blkdiag( I2 , [ eye(My) , zeros(My,My) ] );
end
I = I1 * I2;

for k = 1:K
    for b = 1:B
        h{k,b} = Pt_UL * sparse(I * h{k,b});
    end
end

%% 生成全息感知
%生成参考信号
for k = 1:K
for b = 1:B, x_ref{k,b} = sqrt(1e-3) * exp( 1j * rand(M/4,1)*2*pi ); end
end

%全息图，存储在计算机中
for k = 1:K
    for b = 1:B
        H_holographic{k,b} =  2 * real( conj(h{k,b}) .* x_ref{k,b} )  ; 
    end
end 

theta_grid = ( 2*[0:Mx-1]' - Mx +1 ) / (Mx); phi_grid = ( 2*[0:My-1]' - My +1 ) / (My);
% 导向矢量
b_steering = @(theta,phi)  1/sqrt(M) * I * kron( exp( 1j * 2*pi/lambda * delta_x * dx * theta ) , exp( 1j * 2*pi/lambda * delta_y * dy * phi  ) ) ;


for k = 1:K
    H_recover = zeros(Mx,My); maxColumnValues=[]; rowIndices=[]; maxValue=[]; colIndex=[]; rowIndex=[];
    for b = 1:B
        for num1 = 1:Mx
            for num2 = 1:My
                H_recover(num1,num2) = sqrt(M) * b_steering( theta_grid(num1) , phi_grid(num2) )' .* x_ref{k,b}.' * H_holographic{k,b} ;
            end
        end
        [maxColumnValues, rowIndices] = max(abs(H_recover));  % 每列的最大z值及其行索引
        [maxValue, colIndex] = max(maxColumnValues);  % 找出最大值及其列索引
        rowIndex = rowIndices(colIndex);  % 最大值的行索引
        theta_est_local(k,b) = theta_grid(rowIndex);
        phi_est_local(k,b)   = phi_grid(colIndex)  ;
        est_normal_vector_local{k,b} = [ theta_est_local(k,b) ; phi_est_local(k,b) ; sqrt( max( 1 - ( theta_est_local(k,b)^2 + phi_est_local(k,b)^2 ) , 0 ) ) ];
        est_normal_vector_global{k,b}= Rot_Matrix{b} * est_normal_vector_local{k,b};
        % theta_error(k,b) = sin(phi(k)) * cos(theta(k)) - theta_est(k,b);
        % phi_error(k,b)   = sin(phi(k)) * sin(theta(k)) - phi_est(k,b);
    end

end


for k = 1:K
    theta_est_global_temp = [];  phi_est_global_temp = [];
    for b = 1:B
        A(b,1) = abs(sum(h{k,b}));
    end
    for b=1:B
        if A(b,1) == max(A)
            % est_normal_vector_global{k,b};
            theta_est_global_temp = [ theta_est_global_temp ; est_normal_vector_global{k,b}(1) ] ;
            phi_est_global_temp   = [ phi_est_global_temp ; est_normal_vector_global{k,b}(2) ] ;
            est_normal_vector_global_out(:,k) = est_normal_vector_global{k,b};
        end
    end
    theta_est_global(k,1) = mean(theta_est_global_temp);
    phi_est_global(k,1)   = mean(phi_est_global_temp);
    error_theta(k,1) = theta_est_global(k,1) - e0{k}(1);
    error_phi(k,1)   = phi_est_global(k,1)   - e0{k}(2);
end

for k = 1:K
n1 = [0;0;1]; n2 = est_normal_vector_global_out(:,k);
u = cross(n1,n2);   u = u / norm(u);
% 计算旋转角度（弧度单位）
alph = acos( n1.' * n2 );
% 计算矩阵U
U = [ 0    , -u(3) , u(2) ; ...
      u(3) , 0     , -u(1); ...
     -u(2), u(1)  , 0          ];
% 计算罗德格里斯旋转矩阵
Rot_Matrix{k} = cos(alph) * eye(3) + ( 1 - cos(alph) ) * ( u * u.' ) + sin(alph) * U;
end




end