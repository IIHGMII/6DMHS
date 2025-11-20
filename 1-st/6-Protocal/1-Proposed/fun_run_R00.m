function [fun_out] = fun_run_R00( R00 , H_eff , rho , psi1 , psi2 , sigma0 , Ptx , B , K ,Q , M , H , A_ST , V_F , Xtr , L )

% L=5;

fac = 0; mm_old = zeros(M,B);

ite_num = 1;

while true

    for iii = 1:1

        %% 优化模拟波束
cvx_begin
cvx_solver mosek
variable coff(L+1,K,B) 
variable P0(1)
expressions mm_temp1(M,K,B) mm(M,B) R_temp(K,K) R_thro(K,1) R_temp0(K,1) H_eff(B*Q,K) P_temp(K,K) P_EH(K,1)

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
        P_temp(k1,k) = H_eff(:,k).' * Xtr(:,k1) ;
        if k ~= k1
            R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
        else
            R_temp(k1,k) = 0;
            R_temp0(k) = H_eff(:,k).' * Xtr(:,k);
        end
    end
end

for k = 1:K
    P_EH(k,1) =  2 * sqrt(rho(k)) * real(  psi1(:,k)' * P_temp(:,k) ) - abs( psi1(:,k)' * psi1(:,k) ) ;
    % R_thro(k)      = ( 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/ norm_coff;
    R_thro(k,1)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / log(2);
%     R_thro(k,1)      = ( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/ norm_coff;
end

maximize P0
subject to
P_EH >= P0;
R_thro >= R00;
for b = 1:B
    sum( sum( coff(:,:,b) ) ) == 1;
end
% 
for b = 1:B
    sum( coff(:,b,b) ) == 1;
end

% for k = 1:K
    % sum(sum( sum( coff(:,k,:) ) ) ) == 1;
% end

0 <= coff <= 1;

cvx_end

for k = 1:K
%     % % % % % % % psi1(:,k) =  sqrt(rho(k)) * P_temp(:,k);
    psi1(:,k) = sqrt(rho(k)) * P_temp(:,k);
    psi2(k) = sqrt(1-rho(k)) * R_temp0(k) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end
% % norm_coff = norm(psi2,2);
% % if sum( isnan( psi2 ) ) >= 1 , P_out_ite(ite_num) = 0 ; break ; end
% % % % 重球加速
mm = mm + fac * ( mm - mm_old );
mm = mm/max(max(mm));
mm_old = mm;
fac = max( (ite_num-2)/(ite_num+1) , 0 );

    end

%% 优化数能性能
cvx_begin
cvx_solver mosek
variable Xtr(Q*B,K) complex
variable P0(1)
expressions P_temp(K,K) P_EH(K,1) R_temp(K,K) R_thro(K,1) R_temp0(K,1)

for k  = 1:K
    for k1 = 1:K
            P_temp(k1,k) = H_eff(:,k).' * Xtr(:,k1) ;
            if k ~= k1
                R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
            else
                R_temp(k1,k) = 0;
                R_temp0(k) = H_eff(:,k).' * Xtr(:,k); 
            end
    end
end

for k = 1:K
    P_EH(k,1) =  2 * sqrt(rho(k)) * real(  psi1(:,k)' * P_temp(:,k) ) - abs( psi1(:,k)' * psi1(:,k) ) ;
    % if ite_num == 1
    %    R_thro(k)      =   2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 );
    % else
    %     R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/log(2);
    % end
    R_thro(k,1)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / log(2);

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
%     psi1(:,k) =  sqrt(rho(k)) * P_temp(:,k);
    psi1(:,k) = sqrt(rho(k)) * P_temp(:,k);
    psi2(k) = sqrt(1-rho(k)) * R_temp0(k) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end
norm_coff = norm(psi2,2);
if sum( isnan( psi2 ) ) >= 1 , P_out_ite(ite_num) = 0 ; break ; end


%% 
cvx_begin
cvx_solver mosek
variable rho(K,1) nonnegative
variable P0(1)
expressions P_temp(K,K) P_EH(K,1) R_temp(K,K) R_thro(K,1) R_temp0(K,1)

for k  = 1:K
    for k1 = 1:K
            P_temp(k1,k) = H_eff(:,k).' * Xtr(:,k1) ;
            if k ~= k1
                R_temp(k1,k) = quad_form( H_eff(:,k).' * Xtr(:,k1) , 1 );
            else
                R_temp(k1,k) = 0;
                R_temp0(k) = H_eff(:,k).' * Xtr(:,k); 
            end
    end
end

for k = 1:K
    P_EH(k,1) =  2 * sqrt(rho(k)) * real(  psi1(:,k)' * P_temp(:,k) ) - abs( psi1(:,k)' * psi1(:,k) ) ;
    % if ite_num == 1
    %    R_thro(k)      =   2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 );
    % else
    %     R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/log(2);
    % end
    R_thro(k,1)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * R_temp0(k) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ) / log(2);

end

maximize P0
subject to
P_EH    >= P0;
% R_thro  >= (2^R00-1) ;
R_thro >= R00;
0 <= rho <= 1;
cvx_end

for k = 1:K
%     psi1(:,k) =  sqrt(rho(k)) * P_temp(:,k);
    psi1(:,k) = sqrt(rho(k)) * P_temp(:,k);
    psi2(k) = sqrt(1-rho(k)) * R_temp0(k) / ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) ;
end

norm_coff = norm(psi2,2);

if sum( isnan( psi2 ) ) >= 1 , P_out_ite(ite_num) = 0 ; break ; end


P_out_ite(ite_num) = cvx_optval;

if ite_num >= 50, break; end

ite_num = ite_num + 1;

end

fun_out = max(P_out_ite);

end