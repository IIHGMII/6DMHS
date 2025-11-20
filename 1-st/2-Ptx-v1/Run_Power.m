%% 初始化优化参数
Ptx = 10; sigma0 = 1e-7; R00 = 1; norm_coff = 1;
psi1 = 0.01 * rand(K,K);   psi2 = 0.1 * rand(K,1);
rho  = 0.1 * rand(K,1); %功率分割因子,rho用于传能，(1-rho)用于通信
P_out_ite = [];
for ite_num = 1:10
psi_old = psi2;
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
    R_thro(k)      = ( 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/norm_coff;
end

maximize R0
subject to
R_thro >= R0;
sum( sum( pow_abs( Xtr , 2 ) ) ) <= Ptx;

cvx_end

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









for ite_num = 1:20

%% 优化数能性能
cvx_begin
cvx_solver mosek
variable Xtr(B,K) complex
variable P0(1)
expressions P_temp(K,K) P_EH(K,1) R_temp(K,K) R_thro(K,1)

for k  = 1:K
    for k1 = 1:K
        P_temp(k1,k) =  Xtr(:,k1)' * conj( H_eff(:,k) ) ;
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
        end
    end
end

for k = 1:K
    P_EH(k,1) =  2 * sqrt(rho(k)) * real(  psi1(:,k)' * P_temp(:,k)  ) - abs( psi1(:,k)' * psi1(:,k) ) ;
    
    R_thro(k)      = 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;

end



maximize P0
subject to
P_EH >= P0;
R_thro    >= R00 ;
sum( sum( pow_abs( Xtr , 2 ) ) ) <= Ptx;
cvx_end

for k = 1:K
    psi1(:,k) =  sqrt(rho(k)) * P_temp(:,k);
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end


%% 
cvx_begin
cvx_solver mosek
variable rho(K,1) 
variable P0(1)
expressions P_temp(K,K) P_EH(K,1) R_temp(K,K) R_thro(K,1)

for k  = 1:K
    for k1 = 1:K
        P_temp(k1,k) =  Xtr(:,k1)' * conj( H_eff(:,k) ) ;
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
        end
    end
end

for k = 1:K
    P_EH(k,1) =  2 * sqrt(rho(k)) * real(  psi1(:,k)' * P_temp(:,k)  ) - abs( psi1(:,k)' * psi1(:,k) ) ;
    R_thro(k)      = 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end

maximize P0
subject to
P_EH    >= P0;
R_thro  >= R00 ;
0 <= rho <= 1;
cvx_end

for k = 1:K
    psi1(:,k) =  sqrt(rho(k)) * P_temp(:,k);
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end



P_out_ite(ite_num) = cvx_optval;

end



figure(2); plot( P_out_ite , '-o' )
