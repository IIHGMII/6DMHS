function [ h , e0 , e2 , theta, phi , theta_scatterer , phi_scatterer , IOTA , eta , Omega , h_r , h_c , h_sen , h_x_d , h_x_u , h_y_l , h_y_r ] = Channel_Generation_init( K,B,M,Mx,My,M_x_r , M_y_c ,normal_vector,K_rice,Coor_Ele,Coor_Ele_r,Coor_Ele_c,lambda , Rot_Matrix , Coor_Ele_x_d , Coor_Ele_x_u , Coor_Ele_y_l , Coor_Ele_y_r )

    dx = lambda/4;   %x轴方向天线间隔
    dy = lambda/4;  %y轴方向天线间隔
    delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [0:( My-1)].' - My + 1 )/2;
%% 针对初始6DMA位置，生成无线信道
PL    = 32.4 + 17.3 * log10( 5 ) + 20 * log10( 30 ) - 25 ; % Path-loss,包含20dB天线增益
Omega = 10^( -PL/10 ) * ones(K,B);  % path-loss
% Omega = ones(K,B);
IOTA  = 3;         % 散射体数量
theta =  rand(K,1) * 2 * pi; phi = rand(K,1) * 2 * pi; % K个用户的水平角和俯仰角
% theta =  (rand(K,1)-0.5) * pi; phi = rand(K,1) * 2 * pi; % K个用户的水平角和俯仰角
% theta = 5.263863993147139 ; phi =  1.057064847365129;
theta_scatterer = rand(K,IOTA) * 2 * pi; phi_scatterer = rand(K,IOTA) * 2 * pi; % K个用户的水平角和俯仰角
% NLoS链路信道
for k = 1:K
    for iota = 1:IOTA
        eta(k,iota) = 1/sqrt(IOTA) * sqrt( Omega(k) ) * 1/sqrt(2) * ( randn(1) + 1j*randn(1) );
    end
end

theta = [1.45524585485558
1.41504973364835
2.35835054090831
5.43852179143626];
phi = [5.83563651866559
2.02237603845018
0.0488728727564714
4.43117058629477];
theta_scatterer = [1.56165821938995	0.0891024192032021	4.97921394752858
2.21108342743673	0.787690933376525	3.74825942684378
4.40679921072753	0.638612512771587	3.13896018436890
3.63580537416366	6.14856924227454	5.07064531671004];
phi_scatterer = [2.69668281857476	1.21137536400753	4.64457857472923
6.11843287903727	2.30257164625589	2.92921398237037
4.01675984308597	2.54356244257540	1.59404072941367
5.60346946308695	1.59961515297204	1.83698635921796];

eta = [0.000162840258005646 - 0.00410263506222009i	0.00141051775190082 - 0.00356861079602921i	0.00101196371517660 + 6.34354465897823e-05i
1.32018770279289e-05 - 0.00116493108313416i	0.000188490299469544 + 0.00190905595007172i	-0.00300202684589278 + 0.000658481193069772i
0.00104360111647310 + 5.04415928924031e-05i	0.000150090015563367 - 0.000851331624431704i	0.00205548167310716 - 0.000936050547026122i
0.000794415601009792 - 0.000778357233264543i	0.00120470376911885 - 0.00126613881231049i	-0.00380369315367250 - 0.00117017556283066i];


%% 生成无线信道
% LoS链路信道
for k = 1:K
    for b = 1:B
        h{k,b} = sparse(M,1);%%%%%%%%%%%%%%%%%%%%%注意
            e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h{k,b} = sparse( h{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(M) * sqrt(Omega(k,b)) * Lambda_LoS(k,b) * 1/sqrt(M) * exp( 1j * 2*pi/lambda * Coor_Ele{b} * e0{k} ) );
    end
end

%% 测试
for b = 1:B
    Px{b}    = dx * Rot_Matrix{b} * [1;0;0];
    Py{b}    = dy * Rot_Matrix{b} * [0;1;0];
    T{b}     = exp( 1j * 2*pi/lambda * ( kron( delta_x * ( Px{b}.' * e0{1} ) , ones(1,My) ) + kron( ones(Mx,1) , delta_y.' * ( Py{b}.' * e0{1} ) ) ) );

end








for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            % eta(k,iota,b) = 1/sqrt(IOTA) * Omega(k,b) * 1/sqrt(2) * ( randn(1) + 1j*randn(1) );
            e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            % e = Rot_Matrix{b}.' * e;
            if dot( normal_vector{b} , e1{k,b} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h{k,b} = sparse( h{k,b} + sqrt(1/(K_rice+1)) * sqrt(M) * Lambda(k,iota,b) * eta(k,iota) * 1/sqrt(M) * exp( 1j * 2*pi/lambda * Coor_Ele{b} * e1{k,b} ) );
        end
    end
end
% for k = 1:K
%     H{k} = [];
%     for b = 1:B
%         H{k} = [H{k} ; h{k,b} ];
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第b个RHS的"行"感知元素元素的信道
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_r{k,b} = sparse(M_x_r,1);%%%%%%%%%%%%%%%%%%%%%注意
            e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_r{k,b} = sparse( h_r{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(Mx/2) * sqrt(Omega(k,b)) * Lambda_LoS(k,b) * 1/sqrt(Mx/2) * exp( 1j * 2*pi/lambda * Coor_Ele_r{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            e2{k,iota} = e1{k,b};
            if dot( normal_vector{b} , e1{k,b} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_r{k,b} = sparse( h_r{k,b} + sqrt(1/(K_rice+1)) * sqrt(Mx/2) * Lambda(k,iota,b) * eta(k,iota) * 1/sqrt(Mx/2) * exp( 1j * 2*pi/lambda * Coor_Ele_r{b} * e1{k,b} ) );
        end
    end
end

%% 第b个RHS的"列"感知元素元素的信道
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_c{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
            e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_c{k,b} = sparse( h_c{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(My/2) * sqrt(Omega(k,b)) * Lambda_LoS(k,b) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_c{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,b} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_c{k,b} = sparse( h_c{k,b} + sqrt(1/(K_rice+1)) * sqrt(My/2) * Lambda(k,iota,b) * eta(k,iota) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_c{b} * e1{k,b} ) );
        end
    end
end

%% x,d
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_x_d{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
            e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_x_d{k,b} = sparse( h_x_d{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(My/2) * sqrt(Omega(k,b)) * Lambda_LoS(k,b) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_x_d{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,b} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_x_d{k,b} = sparse( h_x_d{k,b} + sqrt(1/(K_rice+1)) * sqrt(My/2) * Lambda(k,iota,b) * eta(k,iota) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_x_d{b} * e1{k,b} ) );
        end
    end
end


%% x,u
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_x_u{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
            e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_x_u{k,b} = sparse( h_x_u{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(My/2) * sqrt(Omega(k,b)) * Lambda_LoS(k,b) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_x_u{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,b} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_x_u{k,b} = sparse( h_x_u{k,b} + sqrt(1/(K_rice+1)) * sqrt(My/2) * Lambda(k,iota,b) * eta(k,iota) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_x_u{b} * e1{k,b} ) );
        end
    end
end

%% y,l
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_y_l{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
            e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_y_l{k,b} = sparse( h_y_l{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(My/2) * sqrt(Omega(k,b)) * Lambda_LoS(k,b) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_y_l{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,b} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_y_l{k,b} = sparse( h_y_l{k,b} + sqrt(1/(K_rice+1)) * sqrt(My/2) * Lambda(k,iota,b) * eta(k,iota) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_y_l{b} * e1{k,b} ) );
        end
    end
end


%% y,r
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_y_r{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
            e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_y_r{k,b} = sparse( h_y_r{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(My/2) * sqrt(Omega(k,b)) * Lambda_LoS(k,b) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_y_r{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,b} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_y_r{k,b} = sparse( h_y_r{k,b} + sqrt(1/(K_rice+1)) * sqrt(My/2) * Lambda(k,iota,b) * eta(k,iota) * 1/sqrt(My/2) * exp( 1j * 2*pi/lambda * Coor_Ele_y_r{b} * e1{k,b} ) );
        end
    end
end


% 把四个合并起来
for k = 1:K
    for b = 1:B
        h_sen{k,b} = zeros( M_x_r , M_y_c );
        h_sen{k,b}(1,:) = h_x_d{k,b}.';
        h_sen{k,b}(end,:) = h_x_u{k,b}.';
        h_sen{k,b}(:,1) = h_y_l{k,b};
        h_sen{k,b}(:,end) = h_y_r{k,b};
    end
end

end