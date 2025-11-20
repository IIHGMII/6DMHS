function [ h , e0 , e2 , theta, phi , theta_scatterer , phi_scatterer , IOTA , eta , Omega , h_r , h_c ] = Channel_Generation_init( K,B,M,Mx,My,M_x_r , M_y_c ,normal_vector,K_rice,Coor_Ele,Coor_Ele_r,Coor_Ele_c,lambda , Rot_Matrix )

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

theta = [5.65324751373651
0.634257254182679
3.77645183619374
4.38871458749764];
phi = [6.24705677088034
1.18057722500655
3.39885463460762
1.96987960831687];
theta_scatterer = [0.607163484461980	2.10633549812748	5.99823033411429
2.48252740624671	1.69203352022247	2.47406702277745
4.04986358450860	1.90096131481041	1.09482003332386
2.56429383379772	5.10666530962600	5.46840304235066];
phi_scatterer = [3.70417630103166	2.77936846405697	2.65092733203256
2.28396320963008	3.69426329237485	3.88172729270198
2.04250098610203	1.55919163922539	2.87291258043158
2.23391571266596	4.05454413013264	4.53450222339513];

eta = [-3.16112961080035e-05 + 0.00227272321894752i	-0.00195163133032960 + 0.00104791282451845i	-0.00166488740695727 - 0.000782693970061389i
0.000637654914307428 + 0.00164528932966888i	0.000540567805446713 + 0.00304097907835682i	-0.000348157744648930 - 0.000951338442146517i
0.000984095034622721 + 0.00139973823257161i	0.00257268710051583 + 0.00170606263502061i	-0.000253635283657737 + 0.00152587973355016i
-0.000877278955484061 - 0.00261314663279099i	-0.000709397522920676 - 0.00145843620681477i	-0.000257029462448349 - 0.00100029852083689i];


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









end