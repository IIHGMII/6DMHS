for ITE_NUMMMM = 1:5
% 完美信道估计
Pt     = 40 ;
Ptx    = 10^( (Pt - 30)/10 ); sigma0 = 1e-8; R00 = 0;
L      = 10;
% K_rice = 10^( -6.5 / 10 );

% R00 = 18 ;

% K_rice = 10^( 10 / 10 )




%% 生成信道
Omega0 = Omega(1,1);

% 第k个用户的散射角度，包含有Iota个散射链路
es = zeros(3,K,Iota); theta_scatterer = zeros(K,Iota); phi_scatterer = zeros(K,Iota);
theta = theta0; phi = phi0;
theta_scatterer = theta_scatterer0;     phi_scatterer = phi_scatterer0;
for k = 1:K
    for iota = 1 : Iota
%         theta_scatterer(k,iota) = 2 * ( rand(1) - 0.5 ) * 2 * pi; %水平角
%         phi_scatterer(k,iota)   = 2 * ( rand(1) - 0.5 ) * 2 * pi; %俯仰角
        es(:,k,iota)              = [ sin( phi_scatterer(k,iota) ) * cos( theta_scatterer(k,iota) ) ; sin( phi_scatterer(k,iota) ) * sin( theta_scatterer(k,iota) ) ; cos( phi_scatterer(k,iota) ) ];
    end
end

% 生成无线信道
H = zeros( B*M , K );
for k = 1:K
    for iota = 0:Iota
        for b = 1:B
            if iota == 0
                if normal_vector_rot{b}.' * e0{k} > 0
                    H( (b-1)*M+1 : b*M , k ) = H( (b-1)*M+1 : b*M , k ) + sqrt( K_rice/(1+K_rice) ) * sqrt(Omega0) * exp( 1j * 2*pi/lambda * ( s(:,b).' * q * e0{k} +  kron( delta_x * ( Px{b}.' * e0{k} ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b}.' * e0{k} ) )  ) ) ;
                end
            else
                if normal_vector_rot{b}.' * es(:,k,iota) > 0
                    H( (b-1)*M+1 : b*M , k ) = H( (b-1)*M+1 : b*M , k ) + sqrt( 1/(1+K_rice) ) * sqrt(Omega0) * eta(k,iota) * exp( 1j * 2*pi/lambda * ( s(:,b).' * q * es(:,k,iota) +  kron( delta_x * ( Px{b}.' * es(:,k,iota) ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b}.' * es(:,k,iota) ) )  ) ) ;
                end
            end
        end
    end
end
% for k = 1:K
% for l = 1:L
%     addr_max_gain(k) = [  ]
% end
% end


%% 不优化模拟波束，只优化数字波束
% H_eff = [];
% % 生成等效信道
% for b = 1:B
%     for k = 1:K
%         H_eff( (b-1)*Q+1 : b*Q , k ) = ( H( (b-1)*M+1 : b*M , k ).' * diag( mm(:,b) ) * V_F ).' ;
% %         H_eff{b,k} = H( (b-1)*M+1 : b*M , k ).' * diag( mm(:,b) ) * V_F; %等效于有B个射频链路，K个用户的MISO信道
%     end
% end


% 初始化优化参数
norm_coff = 1;
psi1 = 0.1 * randn(K,K) + 0.1j * randn(K,K);   psi2 = 0.01 * rand(K,1) + 0.01j * rand(K,1) ;
rho  = 0.0 * rand(K,1); %功率分割因子,rho用于传能，(1-rho)用于通信
P_out_ite = [];





for k = 1:K

    e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
    
    for iota = 1:L
        e1{k,iota} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
    end

end


for k = 1:K
    for b = 1:B
        for l = 1:(L+1)
            if l == 1
                a_st0(:,l,k,b) = exp( 1j * 2*pi/lambda * ( q(find(s(:,b)==1),:) * e0{k} +  kron( delta_x * ( Px{b}.' * e0{k} ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b}.' * e0{k} ) )  ) );
            elseif l>=2
                a_st0(:,l,k,b) = exp( 1j * 2*pi/lambda * ( q(find(s(:,b)==1),:) * e1{k,l-1} +  kron( delta_x * ( Px{b}.' * e1{k,l-1} ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b}.' * e1{k,l-1} ) )  ) );
            end
        end
    end
end

A_ST = zeros( M , L+1 , K , B );

for b = 1:B
    for k = 1:K
        for l = 1:(L+1)
            A_ST(:,l,k,b) = real( V_F0(:,1).' .* a_st0(:,l,k,b).' + 1 ) / 2;
        end
    end
end

for b = 1:B
    for k = 1:K
        H_eff( (b-1)*Q+1 : b*Q , k ) = ( H( (b-1)*M+1 : b*M , k ).' * diag( mm(:,b) ) * V_F ).' ;
    end
end
 Xtr = randn(Q*B,K) + 1j * randn(Q*B,K);

for ite_num = 1:2

%% 优化模拟波束
cvx_begin
cvx_solver mosek
variable coff(L+1,K,B) 
variable R0(1)
expressions mm_temp1(M,K,B) mm(M,B) R_temp(K,K) R_thro(K,1) R_temp0(K,1) H_eff(B*Q,K)

for k = 1:K
    for b = 1:B
        mm_temp1(:,k,b) =  A_ST(:,:,k,b) * coff(:,k,b);
    end
end

for b = 1:B
    mm(:,b) = mm_temp1(:,:,b) * ones(K,1);
end

for b = 1:B
    for k = 1:K
        H_eff( (b-1)*Q+1 : b*Q , k ) = ( H( (b-1)*M+1 : b*M , k ).' * diag( mm(:,b) ) * V_F ).' ;
    end
end

for k  = 1:K
    for k1 = 1:K
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
            R_temp0(k) = H_eff(:,k).' * Xtr(:,k);
        end
    end
end

for k = 1:K
    % R_thro(k)      = ( 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/ norm_coff;
    R_thro(k,1)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / log(2);
    % R_thro(k,1)      = ( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / norm_coff ;
    % R_thro(k,1)      = (  2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / norm_coff ;
    %     R_thro(k,1)      = ( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/ norm_coff;
end

maximize R0
subject to

R_thro >= R0;

for b = 1:B
    sum( sum( coff(:,:,b) ) ) == 1;
end

for b = 1:B
        sum( coff(:,b,b) ) == 1;
end

0<= coff <= 1;


cvx_end

for k = 1:K
    psi2(k) = sqrt(1-rho(k)) * R_temp0(k) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end

%% 优化数能性能
cvx_begin
cvx_solver mosek
variable Xtr(Q*B,K) complex
variable R0(1)
expressions R_temp(K,K) R_thro(K,1) R_temp0(K,1)

for k  = 1:K
    for k1 = 1:K
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
            R_temp0(k) = H_eff(:,k).' * Xtr(:,k); 
        end
    end
end

for k = 1:K
    % R_thro(k)      = ( 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/ norm_coff;
    R_thro(k,1)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / log(2);
%     R_thro(k,1)      = ( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/ norm_coff;
end

maximize R0
subject to
R_thro >= R0;
sum( sum( pow_abs( Xtr , 2 ) ) ) <= Ptx;

cvx_end

psi_old = psi2;

for k = 1:K
%     psi1(:,k) = P_temp(:,k);
%     psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
    psi2(k) = sqrt(1-rho(k)) * R_temp0(k) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end

P_out_ite(ite_num) = R0 * norm_coff;

norm_coff = norm(psi2,2);

if sum( isnan(psi2) ) >= 1
    psi2 = psi_old;
end


end

psi2_fixed = psi2;


num = 1;        rho = 0.0 * rand(K,1);

for  iiiii = 11

    Ptx = 10^( (30+iiiii-1-30)/10 );

%     psi2 = psi2_fixed;

    [P_ave(iiiii)] = fun_run_R00( R00 , H_eff , rho , psi1 , psi2_fixed , sigma0 , Ptx , B , K , Q , M , H , A_ST , V_F , Xtr , L );

    % P_ave(iiiii) = fun_out;

end

% figure(1); plot( P_ave , '-o' )
P_ITE(ITE_NUMMMM) = max(P_ave);
end


