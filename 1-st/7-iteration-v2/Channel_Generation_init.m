function [ h , e0 , e2 , theta, phi , theta_scatterer , phi_scatterer , IOTA , eta , Omega , h_r , h_c ] = Channel_Generation_init( K,B,M,Mx,My,normal_vector,K_rice,Coor_Ele,Coor_Ele_r,Coor_Ele_c,lambda , Rot_Matrix )

    dx = lambda/4;   %x轴方向天线间隔
    dy = lambda/4;  %y轴方向天线间隔
    delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [0:( My-1)].' - My + 1 )/2;
%% 针对初始6DMA位置，生成无线信道
PL    = 32.4 + 17.3 * log10( 6.5 ) + 20 * log10( 30 ) - 15 ; % Path-loss,包含20dB天线增益
Omega = 10^( -PL/10 ) * ones(K,B);  % path-loss
% Omega = ones(K,B);
IOTA  = 20;         % 散射体数量
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
        h_r{k,b} = sparse(Mx/2,1);%%%%%%%%%%%%%%%%%%%%%注意
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
        h_c{k,b} = sparse(My/2,1);%%%%%%%%%%%%%%%%%%%%%注意
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