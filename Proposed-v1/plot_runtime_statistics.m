%% plot_runtime_statistics.m
% 读取 Monte Carlo 统计结果，绘制：
% - Fig A: Stage I/II/III/Total 运行时间箱线图（含 T_c 红色虚线）
% - Fig B: Stage I/II/III 运行时间 CDF（含 T_c 垂直虚线）

clc; clear;

in_file = fullfile('Data_record', 'runtime_statistics_mc100.mat');
S = load(in_file);

time_stage1_all = S.time_stage1_all(:);
time_stage2_all = S.time_stage2_all(:);
time_stage3_all = S.time_stage3_all(:);
time_total_all  = S.time_total_all(:);

% 相干时间 T_c（ms）
if isfield(S, 'coherence') && isfield(S.coherence, 'Tc_ms')
    Tc_ms = S.coherence.Tc_ms;
    Tc_note = sprintf('T_c=%.3f ms', Tc_ms);
else
    % 兜底：按论文给定 T_c=120 symbols，并假设用户速度/载频计算
    fc = 30e9; v_kmh = 3; Tc_symbols = 120;
    c = 3e8; lambda = c/fc; v = v_kmh/3.6; fD = v/lambda;
    Tc_ms = (0.423/fD) * 1e3;
    symbol_rate = Tc_symbols / (Tc_ms/1e3); %#ok<NASGU>
    Tc_note = sprintf('T_c=%.3f ms (fallback)', Tc_ms);
end

%% Fig A: Box Plot
figure('Color', 'w');
data = [time_stage1_all; time_stage2_all; time_stage3_all; time_total_all];
group = [repmat({'Stage I'}, numel(time_stage1_all), 1); ...
         repmat({'Stage II'}, numel(time_stage2_all), 1); ...
         repmat({'Stage III'}, numel(time_stage3_all), 1); ...
         repmat({'Total'}, numel(time_total_all), 1)];

boxplot(data, group);
grid on;
ylabel('Runtime (ms)');
title('Fig A: Runtime Box Plot (MC=100)');
hold on;
yline(Tc_ms, 'r--', Tc_note, 'LineWidth', 1.5);

%% Fig B: CDF
figure('Color', 'w');
hold on; grid on;

plot_empirical_cdf(time_stage1_all, 'Stage I');
plot_empirical_cdf(time_stage2_all, 'Stage II');
plot_empirical_cdf(time_stage3_all, 'Stage III');

xline(Tc_ms, 'r--', Tc_note, 'LineWidth', 1.5);
xlabel('Runtime (ms)');
ylabel('CDF');
title('Fig B: Runtime CDF (MC=100)');
legend('Location', 'southeast');

%% ====== 本地函数 ======
function plot_empirical_cdf(x, label_name)
    x = x(:);
    x = x(isfinite(x));
    if isempty(x)
        return;
    end
    x = sort(x);
    n = numel(x);
    y = (1:n).' / n;
    plot(x, y, 'LineWidth', 2, 'DisplayName', label_name);
end
