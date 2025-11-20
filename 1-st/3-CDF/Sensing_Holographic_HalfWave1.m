function [ error_theta , error_phi , error_varphi , est_direc_vector_rec , Rot_Matrix , normal_vector_rot ] = Sensing_Holographic_HalfWave1(B,M,Mx,My,lambda,delta_x,delta_y,dx,dy,K,h0,Rot_Matrix,e0)
%以半波长为间隔，生成全息图
kappa0 = 17;
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

Mx0 = 1:(kappa0-1):Mx;My0 = 1:(kappa0-1):My;

for k = 1:K
    for b = 1:B
        h1 = reshape( h0{k,b} , Mx , My );
%         h_temp = I * h0{k,b};
%         h{k,b} = sqrt(Pt_UL) * reshape( h_temp , Mx/2 , My/2 ).' ;
        h{k,b} = sqrt(Pt_UL) * h1(Mx0,My0).';
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
        H_holographic{k,b} =  2 * real( conj(h{k,b}) .* x_ref{k,b} + n{k,b}  )  ; 
    end
end 

theta_grid = ( 2*[0:Mx-1]' - Mx +1 ) / (Mx); phi_grid = ( 2*[0:My-1]' - My +1 ) / (My);
% 导向矢量
b_steering = @(theta,phi)  1/sqrt(M) * kron( exp( 1j * 2*pi/lambda * delta_x * dx * theta ) , exp( 1j * 2*pi/lambda * delta_y * dy * phi  ) ) ;

% 构造2D-FFT矩阵
% delta_x = delta_x(1:2:end); delta_y = delta_y(1:2:end); 
delta_x = ( 2 * [0:((length(Mx0))-1)].' - (length(My0)) + 1 )/2; delta_y = ( 2 * [0:( (length(My0))-1)].' - (length(My0)) + 1 )/2;
F =  exp( -1j * 2*pi * delta_x * delta_y.' / (length(Mx0)) );

for k = 1:K
    for b = 1:B
        Px{b}    = dx * Rot_Matrix{b} * [1;0;0];
        Py{b}    = dy * Rot_Matrix{b} * [0;1;0];
        T{k,b}   = x_ref{k,b} .* H_holographic{k,b};
        % T{b,k}     = exp( 1j * 2*pi/lambda * ( kron( delta_x * ( Px{b}.' * e0{k} ) , ones(1,My) ) + kron( ones(Mx,1) , delta_y.' * ( Py{b}.' * e0{k} ) ) ) );
        Sen_T{k,b} = abs( F * T{k,b} * F.' );
        [maxValues, rowIndices] = max(Sen_T{k,b});
        % 从列最大值中找到最大值及其列索引
        [maxValue, colIndex{k,b}] = max(maxValues);
        rowIndex{k,b} = rowIndices(colIndex{k,b});
        A(k,b)  = sum( sum(abs(h{k,b})));
    end
end

for k = 1:K
    A0(k,:) = sort(A(k,:),'descend');
    Qx=[];Qy=[];cx=[];cy=[];
    for b = 1:B
        if A(k,b) == A0(k,1) || A(k,b) == A0(k,2) || A(k,b) == A0(k,3)
            Qx = [ Qx ; Px{b}.' ];
            Qy = [ Qy ; Py{b}.' ];
            cx = [ cx ; delta_x(rowIndex{k,b}) / Mx * lambda ];
            cy = [ cy ; delta_y(colIndex{k,b}) / My * lambda ];
        end
    end
Q = [Qx;Qy];    c0 = [cx;cy];
t = Q\c0;
%有误差
est_direc_vector_rec{k} = t / norm(t);
%无误差
% est_direc_vector_rec{k} = e0{k};



error_theta(k,1) = est_direc_vector_rec{k}(1,1) - e0{k}(1);
error_phi(k,1)   = est_direc_vector_rec{k}(2,1) - e0{k}(2);
error_varphi(k,1)= est_direc_vector_rec{k}(3,1) - e0{k}(3);

end



% for k = 1:K
%     H_recover = zeros(Mx,My); maxColumnValues=[]; rowIndices=[]; maxValue=[]; colIndex=[]; rowIndex=[];
%     for b = 1:B
%         for num1 = 1:Mx
%             for num2 = 1:My
%                 H_recover(num1,num2) = sqrt(M) * b_steering( theta_grid(num1) , phi_grid(num2) )' .* x_ref{k,b}.' * H_holographic{k,b} ;
%             end
%         end
%         [maxColumnValues, rowIndices] = max(abs(H_recover));  % 每列的最大z值及其行索引
%         [maxValue, colIndex] = max(maxColumnValues);  % 找出最大值及其列索引
%         rowIndex = rowIndices(colIndex);  % 最大值的行索引
%         theta_est_local(k,b) = theta_grid(rowIndex);
%         phi_est_local(k,b)   = phi_grid(colIndex)  ;
%         est_normal_vector_local{k,b} = [ theta_est_local(k,b) ; phi_est_local(k,b) ; sqrt( max( 1 - ( theta_est_local(k,b)^2 + phi_est_local(k,b)^2 ) , 0 ) ) ];
%         est_normal_vector_global{k,b}= Rot_Matrix{b} * est_normal_vector_local{k,b};
%         % theta_error(k,b) = sin(phi(k)) * cos(theta(k)) - theta_est(k,b);
%         % phi_error(k,b)   = sin(phi(k)) * sin(theta(k)) - phi_est(k,b);
%     end
% 
% end
% 
% 
% for k = 1:K
%     theta_est_global_temp = [];  phi_est_global_temp = [];
%     for b = 1:B
%         A(b,1) = abs(sum(h{k,b}));
%     end
%     for b=1:B
%         if A(b,1) == max(A)
%             % est_normal_vector_global{k,b};
%             theta_est_global_temp = [ theta_est_global_temp ; est_normal_vector_global{k,b}(1) ] ;
%             phi_est_global_temp   = [ phi_est_global_temp ; est_normal_vector_global{k,b}(2) ] ;
%             est_normal_vector_global_out(:,k) = est_normal_vector_global{k,b};
%         end
%     end
%     theta_est_global(k,1) = mean(theta_est_global_temp);
%     phi_est_global(k,1)   = mean(phi_est_global_temp);
%     error_theta(k,1) = theta_est_global(k,1) - e0{k}(1);
%     error_phi(k,1)   = phi_est_global(k,1)   - e0{k}(2);
% end

% for k = 1:K
% n1 = [0;0;1]; n2 = est_normal_vector_global_out(:,k);
% u = cross(n1,n2);   u = u / norm(u);
% % 计算旋转角度（弧度单位）
% alph = acos( n1.' * n2 );
% % 计算矩阵U
% U = [ 0    , -u(3) , u(2) ; ...
%       u(3) , 0     , -u(1); ...
%      -u(2), u(1)  , 0          ];
% % 计算罗德格里斯旋转矩阵
% Rot_Matrix{k} = cos(alph) * eye(3) + ( 1 - cos(alph) ) * ( u * u.' ) + sin(alph) * U;
% end

% 旋转最大增益角
for k = 1:K
n1 = [ 0.02 ; 0.01 ; sqrt( 1 - 0.02^2 - 0.01^2 ) ]; %最大增益方向
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


