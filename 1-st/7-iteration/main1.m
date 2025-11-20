% 6DMA,使用全息感知方法
clc;clear;

tic

kappa0 = 3; % 2 3 5 9 17



Pt = 40;
Ptx = 10^( (Pt - 30)/10 ); sigma0 = 1e-7; R00 = 1;
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
    Mx = 64; My = 64; M = Mx * My;
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

N = 1024;    %FFT采样点数
delta_n = [ (-N+1) / 2 : 1 : (N-1) /2 ].';
F = exp( - 1j * 2 * pi * delta_n * delta_n.' / N   );
% F_x =  exp( -1j * 2*pi * delta_x_r * delta_x_r.' / M_x_r ); %用于"行"感知单元感知信息提取的FFT矩阵
% F_y =  exp( -1j * 2*pi * delta_y_c * delta_y_c.' / M_y_c ); %用于"列"感知单元感知信息提取的FFT矩阵

for k = 1:K
    for b = 1:B
        hx{k,b} = [ zeros((N-M_x_r)/2,1) ; h_r{k,b} ; zeros((N-M_x_r)/2,1) ];
        hy{k,b} = [ zeros((N-M_y_c)/2,1) ; h_c{k,b} ; zeros((N-M_y_c)/2,1) ];
    end
end

for k = 1:K
    for b = 1:B
        P_x_sum(k,b) = sum(abs(hx{k,b}));
        if P_x_sum(k,b) > 1e-2
            flag_x{k,b} = 1;
            e_est_index{k,b}(1) = find( abs( F * hx{k,b} ) == max(abs( F * hx{k,b} )) );
            e_est{k,b}(1) = delta_n(e_est_index{k,b}(1)) / N * 2;
        else
            flag_x{k,b} = 0;
            e_est{k,b}(1) = nan;
        end

        P_y_sum(k,b) = sum( abs(hy{k,b}) );
        if P_y_sum(k,b) > 1e-2
            flag_y{k,b} = 1;
            e_est_index{k,b}(2) = find( abs( F * hy{k,b} ) == max(abs( F * hy{k,b} )) );
            e_est{k,b}(2) = delta_n(e_est_index{k,b}(2)) / N * 2;
        else
            flag_y{k,b} = 0;
            e_est{k,b}(2) = nan;
        end        

    end
end




for k = 1:K
    index_h_x(k) = find(P_x_sum(k,:) == max(P_x_sum(k,:)));
    index_h_y(k) = find(P_y_sum(k,:) == max(P_y_sum(k,:)));
    index_h_x_min(k) = find(P_x_sum(k,:) == min(P_x_sum(k,:)));
    index_h_y_min(k) = find(P_y_sum(k,:) == min(P_y_sum(k,:)));
end

for k = 1:K
    e0{k}' * ( Rot_Matrix{index_h_x(k)} * [1;0;0] )   - e_est{k,index_h_x(k)}(1);
    e0{k}' * ( Rot_Matrix{index_h_y(k)} * [0;1;0] )   - e_est{k,index_h_y(k)}(2);
end

%% 估计方向向量的第三个元素

e_est0{k}  = [ e_est{k,index_h_x(k)}(1) ; e_est{k,index_h_x(k)}(2) ; sqrt( 1 - e_est{k,index_h_x(k)}(1)^2 - e_est{k,index_h_x(k)}(2)^2 )  ];
% Coor_q_init( index_h_x(k) , : ) * e_est0{k} ;
% Coor_q_init( index_h_x_min(k) , : ) * e_est0{k};
e1{k} = Rot_Matrix{index_h_x(k)} \ e_est0{k};
Coor_q_init( index_h_x(k) , : ) * e1{k};
Coor_q_init( index_h_x_min(k) , : ) * e1{k};

e_est0{k} = [ e_est{k,index_h_x(k)}(1) ; e_est{k,index_h_x(k)}(2) ; -sqrt( 1 - e_est{k,index_h_x(k)}(1)^2 - e_est{k,index_h_x(k)}(2)^2 )  ];
% Coor_q_init( index_h_x(k) , : ) * e_est0{k} ;
% Coor_q_init( index_h_x_min(k) , : ) * e_est0{k} ;
e2{k} = Rot_Matrix{index_h_x(k)} \ e_est0{k};
Coor_q_init( index_h_x(k) , : ) * e2{k};
Coor_q_init( index_h_x_min(k) , : ) * e2{k};

A = [Rot_Matrix{index_h_x(k)}(:,1).' ; Rot_Matrix{index_h_x(k)}(:,2).' ];
y = [e_est{k,index_h_x(k)}(1);e_est{k,index_h_y(k)}(2)];
% x = sphere_only_solution(A, y)
cvx_begin
variable x(3,1)

minimize ( norm( y - [x.' * Rot_Matrix{index_h_x(k)} * [ 1 ; 0 ; 0] ; x.' * Rot_Matrix{index_h_x(k)} * [ 0 ; 1 ; 0]  ] ) )

% minimize norm( y - A *x , 2 )
subject to
norm(x,2) <= 1;
Coor_q_init( index_h_x(k) , : ) * x > 0;
cvx_end
norm(x,2)
x / norm(x,2)
% e_est0{k} 
e0{k}
%最小二乘法
% A{k} = [];
% A{k} = [(Rot_Matrix{index_h_x(k)} * [1;0;0])' ; (Rot_Matrix{index_h_y(k)} * [0;1;0])'];
% C{k} = [ e_est{k,index_h_x(k)}(1) ;e_est{k,index_h_y(k)}(2) ];
% A{k}\C{k}
% e0{k}
% e0{1}' * ( Rot_Matrix{1} * [1;0;0] )   - e_est{1,1}
% e0{1}' * ( Rot_Matrix{2} * [1;0;0] )   - e_est{1,2}
% e0{1}' * ( Rot_Matrix{3} * [1;0;0] )   - e_est{1,3}
% e0{1}' * ( Rot_Matrix{4} * [1;0;0] )   - e_est{1,4}
% e0{1}' * ( Rot_Matrix{5} * [1;0;0] )   - e_est{1,5}(1)
% e0{1}' * ( Rot_Matrix{6} * [1;0;0] )   - e_est{1,6}







toc
function x = sphere_only_solution(A, y)
% Solve: min ||y - A x||  s.t. ||x||=1
% Robust secular-equation solver on lambda in (-sigma_min^2, +inf)

    y = y(:);
    [U,S,V] = svd(A,'econ');
    sigma = diag(S);             % r x 1
    r = numel(sigma);
    b = A.'*y;
    c = V.'*b;                   % 与右奇异向量对齐的系数
    s2 = sigma.^2;

    % 只挑选 |c_i| 有效分量，对应的最小 sigma^2 作为左端点
    tol_c = 1e-12 * norm(c);
    idx = find(abs(c(1:r)) > tol_c);       % 有能量的奇异分量
    if isempty(idx)
        % A^T y 全零：目标退化为 min_{||x||=1} x^T A^T A x
        % 取 A^T A 的最小特征向量（即最小奇异向量）
        x = V(:, end);                      % 对应最小奇异值
        x = x / norm(x);
        return;
    end
    s2_active_min = min(s2(idx));
    % 设定搜索区间：左端稍大于 -s2_active_min，右端很大正数
    lamL = -s2_active_min + max(1e-15, 1e-12 * max(1, s2_active_min)); 
    lamR = max(1, 1e6 * max(s2) + 1);

    % 世俗方程（只对有效分量求和也可；对所有分量更稳健）
    f = @(lam) sum( (c(1:r).^2) ./ (s2 + lam).^2 ) - 1;

    % 确保端点异号：f(lamL) 应该 >0（趋于 +inf），f(lamR) < 0（-> -1）
    fL = f(lamL);
    % 若数值上没变正，再把 lamL 向左逼近（但不能越过 -s2_active_min）
    tries = 0;
    while ~(isfinite(fL) && fL > 0) && tries < 5
        lamL = (lamL + (-s2_active_min)) / 2;   % 向左靠近
        fL = f(lamL); tries = tries + 1;
    end

    % 右端点直到为负
    fR = f(lamR);
    while ~(isfinite(fR) && fR < 0)
        lamR = lamR * 2;
        fR = f(lamR);
        if lamR > 1e20, break; end
    end

    % 若仍未异号，用二分在 (lamL, lamR) 上保守求根
    if ~(fL > 0 && fR < 0)
        % 退化情况下用安全二分（单调递减保证可用）
        lam = bisect_monotone(f, lamL, lamR, 200, 1e-12);
    else
        lam = fzero(f, [lamL, lamR]);      % 成功异号则用 fzero
    end

    % 回代得到 x，并归一化
    z = zeros(size(V,2),1);
    z(1:r) = c(1:r) ./ (s2 + lam);
    x = V * z;
    x = x / norm(x);
end

function lam = bisect_monotone(f, a, b, maxit, tol)
% 对单调递减函数 f 在 [a,b] 求根（无需严格异号）
    fa = f(a); fb = f(b);
    for k = 1:maxit
        m = 0.5*(a+b); fm = f(m);
        if ~isfinite(fm), b = m; continue; end
        % 目标是 f(m)=0；若 fm>0，根在右侧（因 f 递减）
        if fm > 0
            a = m; fa = fm;
        else
            b = m; fb = fm;
        end
        if abs(b-a) <= tol*(1+abs(m)), break; end
    end
    lam = 0.5*(a+b);
end
