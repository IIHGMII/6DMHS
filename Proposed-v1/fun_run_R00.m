function [fun_out, history] = fun_run_R00( R00 , H_eff , rho , psi1 , psi2 , sigma0 , Ptx , B , K , varargin )
%FUN_RUN_R00 采用FP方法在给定吞吐阈值R00下最大化最小EH功率
% 兼容旧接口：仅给前9个参数时，行为与旧版一致（但默认不再弹出迭代曲线图）。
%
% 可选参数（以结构体opts形式传入）：
%   opts.max_iter (默认 50)      FP外层最大迭代次数
%   opts.quiet    (默认 true)    使用cvx_begin quiet
%   opts.do_plot  (默认 false)   是否绘制P_out_ite收敛曲线
%
% 输出：
%   fun_out  : 本次运行的最优(迭代中最大)目标值
%   history  : 结构体，包含P_out_ite等中间量

opts = struct('max_iter', 50, 'quiet', true, 'do_plot', false, 'sigma_cov', 0);
if ~isempty(varargin) && isstruct(varargin{1})
    user_opts = varargin{1};
    fns = fieldnames(user_opts);
    for i = 1:numel(fns)
        opts.(fns{i}) = user_opts.(fns{i});
    end
end

% 预分配，避免在早期异常(例如CVX不可行导致psi2为NaN)时引用未定义变量
P_out_ite = -inf(1, opts.max_iter);

ite_num = 1;
while true

%% 优化数能性能（给定rho）
if opts.quiet
    cvx_begin quiet
else
    cvx_begin
end
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
    % if ite_num == 1
    %    R_thro(k)      =   2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 );
    % else
    %     R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/log(2);
    % end
    % log2(1+·) = log(1+·)/log(2)
    noise_term = (1 - rho(k)) * sigma0 + opts.sigma_cov;
    R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) ...
        - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + noise_term ) ) / log(2);

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
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + (1-rho(k)) * sigma0 + opts.sigma_cov ) ;
end
norm_coff = norm(psi2,2);
if sum( isnan( psi2 ) ) >= 1
    if ite_num > 1
        P_out_ite(ite_num) = max(P_out_ite(1:ite_num-1));
    else
        P_out_ite(ite_num) = -inf;
    end
    break;
end


%% 
if opts.quiet
    cvx_begin quiet
else
    cvx_begin
end
cvx_solver mosek
variable rho(K,1) nonnegative
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
    % if ite_num == 1
    %    R_thro(k)      =   2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 );
    % else
    %     R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + sigma0 ) )/log(2);
    % end
    noise_term = (1 - rho(k)) * sigma0 + opts.sigma_cov;
    R_thro(k)      = log( 1 + 2 * sqrt(1-rho(k)) * real( psi2(k)' * Xtr(:,k)' * conj( H_eff(:,k) ) ) ...
        - abs(psi2(k))^2 * ( (1-rho(k)) * sum( R_temp(:,k) ) + noise_term ) ) / log(2);

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
    psi2(k)   = sqrt(1-rho(k)) * Xtr(:,k)' * conj( H_eff(:,k) ) / ( (1-rho(k)) * sum( R_temp(:,k) ) + (1-rho(k)) * sigma0 + opts.sigma_cov ) ;
end
norm_coff = norm(psi2,2);

if sum( isnan( psi2 ) ) >= 1
    if ite_num > 1
        P_out_ite(ite_num) = max(P_out_ite(1:ite_num-1));
    else
        P_out_ite(ite_num) = -inf;
    end
    break;
end


P_out_ite(ite_num) = cvx_optval;

if ite_num >= opts.max_iter, break; end

ite_num = ite_num + 1;

end
history = struct();
history.P_out_ite = P_out_ite;

if opts.do_plot
    figure; plot(P_out_ite);
    xlabel('迭代次数');
    ylabel('目标值(最小EH功率)');
    grid on;
end

fun_out = max(P_out_ite);

end
