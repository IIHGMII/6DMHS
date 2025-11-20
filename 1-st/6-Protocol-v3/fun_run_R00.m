function [fun_out] = fun_run_R00( R00 , H_eff , rho , psi1 , psi2 , sigma0 , Ptx , B , K  )

% H_eff0 = 1/sqrt(sigma0) * H_eff;
H_eff0 =  H_eff;

ite_num = 1;
while true

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
            R_temp(k1,k) = quad_form( H_eff0(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
        end
    end
end

for k = 1:K
    P_EH(k,1) =  2 * sqrt(rho(k)) * real(  psi1(:,k)' * P_temp(:,k)  ) - abs( psi1(:,k)' * psi1(:,k) ) ;
    % if ite_num == 1
    %    R_thro(k)      =   2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 );
    % else
    %     R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/log(2);
    % end
    R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff0(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + 1 ) ) / log(2);

end



maximize P0
subject to
P_EH >= P0;
% R_thro  >= (2^R00-1)/norm_coff ;
% R_thro  >= (2^R00-1) ;
R_thro >= R00;
sum( sum( pow_abs( Xtr , 2 ) ) ) <= Ptx;
cvx_end


for k = 1:K
    psi1(:,k) =  sqrt(rho(k)) * P_temp(:,k);
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff0(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + 1 ) ;
end
norm_coff = norm(psi2,2);
if sum( isnan( psi2 ) ) >= 1 , P_out_ite(ite_num) = 0 ; break ; end


%% 
cvx_begin
cvx_solver mosek
variable rho(K,1) nonnegative
variable P0(1)
expressions P_temp(K,K) P_EH(K,1) R_temp(K,K) R_thro(K,1)

for k  = 1:K
    for k1 = 1:K
        P_temp(k1,k) =  Xtr(:,k1)' * conj( H_eff(:,k) ) ;
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff0(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
        end
    end
end

for k = 1:K
    P_EH(k,1) =  2 * sqrt(rho(k)) * real(  psi1(:,k)' * P_temp(:,k)  ) - abs( psi1(:,k)' * psi1(:,k) ) ;
    % if ite_num == 1
    %    R_thro(k)      =   2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 );
    % else
    %     R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/log(2);
    % end
    R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff0(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + 1 ) ) / log(2);

end

maximize P0
subject to
P_EH    >= P0;
% R_thro  >= (2^R00-1) ;
R_thro >= R00;
0 <= rho <= 1;
cvx_end


for k = 1:K
    psi1(:,k) =  sqrt(rho(k)) * P_temp(:,k);
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff0(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + 1 ) ;
end
norm_coff = norm(psi2,2);

if sum( isnan( psi2 ) ) >= 1 , P_out_ite(ite_num) = 0 ; break ; end


P_out_ite(ite_num) = cvx_optval;

if ite_num >=25, break; end

ite_num = ite_num + 1;

end

fun_out = max(P_out_ite) * sigma0;

end