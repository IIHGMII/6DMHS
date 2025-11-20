%% 生成信道
Omega0 = Omega(1,1);

% 第k个用户的散射角度，包含有Iota个散射链路
es = zeros(3,K,Iota); theta_scatterer = zeros(K,Iota); phi_scatterer = zeros(K,Iota);
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
                    H( (b-1)*M+1 : b*M , k ) = H( (b-1)*M+1 : b*M , k ) + sqrt( 1/(1+K_rice) ) *  eta(k,iota) * exp( 1j * 2*pi/lambda * ( s(:,b).' * q * es(:,k,iota) +  kron( delta_x * ( Px{b}.' * es(:,k,iota) ) , ones(My,1) ) + kron( ones(Mx,1) , delta_y * ( Py{b}.' * es(:,k,iota) ) )  ) ) ;
                end
            end
        end
    end
end

%% 不优化模拟波束，只优化数字波束
% 生成等效信道
for b = 1:B
    for k = 1:K
        H_eff(b,k) = H( (b-1)*M+1 : b*M , k ).' * diag( mm(:,b) ) * V_F; %等效于有B个射频链路，K个用户的MISO信道
    end
end

% 初始化优化参数
norm_coff = 1;
psi1 = 0.005 * randn(K,K) + 0.005j * randn(K,K);   psi2 = 0.001 * rand(K,1) + 0.001j * rand(K,1) ;
rho  = 0.0000001 * rand(K,1); %功率分割因子,rho用于传能，(1-rho)用于通信
P_out_ite = [];
for ite_num = 1:20

%% 优化数能性能
cvx_begin
cvx_solver mosek
variable Xtr(B,K) complex
variable R0(1)
expressions R_temp(K,K) R_thro(K,1)

for k  = 1:K
    for k1 = 1:K
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
        end
    end
end

for k = 1:K
    % R_thro(k)      = ( 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/ norm_coff;
    R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / log(2);

end

maximize R0
subject to
R_thro >= R0;
sum( sum( pow_abs( Xtr , 2 ) ) ) <= Ptx;

cvx_end

psi_old = psi2;

for k = 1:K
%     psi1(:,k) = P_temp(:,k);
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end

P_out_ite(ite_num) = R0 * norm_coff;

norm_coff = norm(psi2,2);

if sum( isnan(psi2) ) >= 1
    psi2 = psi_old;
end


end

psi2_fixed = psi2;


num = 1;        rho = 0.000000000000 * rand(K,1);

for  iiiii = 1:1

    R00 = (iiiii-1) * 1;

%     psi2 = psi2_fixed;

    [P_ave(iiiii)] = fun_run_R00( R00 , H_eff , rho , psi1 , psi2_fixed , sigma0 , Ptx , B , K );

    % P_ave(iiiii) = fun_out;

end

P_out = P_ave;

figure;plot(P_out)

%无误差
P1 = [0.00230555049009775	0.00225010292567540	0.00212225689264841	0.00200780886558690	0.00192265678721690	0.00186173621485415	0.00181738349260611	0.00178216367354857	0.00176503146597107	0.00174695108703070	0.00173180279032037	0.00171197141757212	0.00169283026857875	0.00165081989590456	0.00155271685739499	8.99798073181712e-05];


%MSE = 3.57e-4,RMSE = 0.018894443627691
P2 = [0.00195464499778940	0.00191938688411638	0.00182922546819305	0.00175076359727757	0.00169025502404063	0.00164584174499974	0.00161408811238635	0.00159070340510185	0.00157390561290566	0.00156018769832748	0.00154486950958741	0.00152451279238886	0.00151211337305289	0.00146793971249838	0.00138412782674177	0];

% MSE = 0.0011; RMSE = 0.031;
P3 = [0.00189678889412516	0.00187013322780967	0.00178053281938157	0.00169274497314632	0.00162351379558286	0.00157306070416075	0.00153528494562230	0.00150859205072502	0.00149060571655986	0.00147609066055030	0.00146186699166452	0.00144546685935212	0.00141873563988119	0.00138132109523695	0.00128222792144648	0];



% MSE = 0.007; RMSE = 0.083;
P4 = [0.000114307417387126	0.000107355437764615	7.09641374919613e-05	5.46686036104717e-05	4.69783411051652e-05	4.39806991293168e-05	4.13526543941319e-05	4.04341011649568e-05	3.87742373292339e-05	3.57851563629057e-05	3.07566306943398e-05	1.99577287507151e-05	0	0	0	0];

figure;plot(1:16,10*log10(P1/1e-3),1:16,10*log10(P2/1e-3),1:16,10*log10(P3/1e-3),1:16,10*log10(P4/1e-3))
figure;plot(1:16,10*log10(P1/1e-3),1:16,10*log10(P2/1e-3),1:16,10*log10(P3/1e-3))
axis([1,15,1, 4])
