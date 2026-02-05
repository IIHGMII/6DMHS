%% run_monte_carlo_runtime_statistics.m
% 目的：运行 100 次 Monte Carlo，统计 Stage I/II/III 运行时间与迭代次数，并保存 .mat 结果。
% 注意：计时已在 RUN_OPT_Protocol.m 内部完成（仅覆盖算法核心计算），本脚本不把数据生成/保存/绘图计入计时。

clc; clear;

%% ========================= 可调参数 =========================
num_mc = 100;
K_rice_dB = 10;                 % 论文常用：10 dB
K_rice = 10^(K_rice_dB/10);

out_dir = 'Data_record';
out_file = fullfile(out_dir, 'runtime_statistics_mc100.mat');

% 相干时间（论文：T_c = 120 symbols）
coherence.fc = 30e9;
coherence.v_kmh = 3;            % 假设用户速度（可改：如 30 表示 30 km/h）
coherence.Tc_symbols = 120;
coherence.Tc_formula = 'Tc_sec = 0.423/fD, fD=v/lambda';

% solver 信息（不改变算法，仅用于记录）
solver_info = struct();
solver_info.name = 'CVX/MOSEK';
solver_info.tolerance = 1e-4;
solver_info.max_iters = 50;     % fun_run_R00 默认 opts.max_iter=50（Stage III FP）

use_parfor = false;             % 如有 Parallel Computing Toolbox 可改为 true
%% ==========================================================

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% 运行时间 (单位: ms)
time_stage1_all = NaN(1, num_mc);
time_stage2_all = NaN(1, num_mc);
time_stage3_all = NaN(1, num_mc);

% 迭代次数
iter_stage2_AO_all = NaN(1, num_mc);
iter_stage3_AO_all = NaN(1, num_mc);
iter_stage3_FP_all = NaN(1, num_mc);

% CVX 调用次数（如可获取）
cvx_calls_stage2_all = NaN(1, num_mc);
cvx_calls_stage3_all = NaN(1, num_mc);
cvx_calls_total_all  = NaN(1, num_mc);

% Monte Carlo
if use_parfor
    parfor mc = 1:num_mc
        [time_stage1_all(mc), time_stage2_all(mc), time_stage3_all(mc), ...
            iter_stage2_AO_all(mc), iter_stage3_AO_all(mc), iter_stage3_FP_all(mc), ...
            cvx_calls_stage2_all(mc), cvx_calls_stage3_all(mc), cvx_calls_total_all(mc)] = ...
            one_run(mc, K_rice);
    end
else
    for mc = 1:num_mc
        [time_stage1_all(mc), time_stage2_all(mc), time_stage3_all(mc), ...
            iter_stage2_AO_all(mc), iter_stage3_AO_all(mc), iter_stage3_FP_all(mc), ...
            cvx_calls_stage2_all(mc), cvx_calls_stage3_all(mc), cvx_calls_total_all(mc)] = ...
            one_run(mc, K_rice);
    end
end

time_total_all = time_stage1_all + time_stage2_all + time_stage3_all;

%% 统计量：平均/最大/最小/标准差
stats = struct();
stats.time.stage1 = basic_stats(time_stage1_all);
stats.time.stage2 = basic_stats(time_stage2_all);
stats.time.stage3 = basic_stats(time_stage3_all);
stats.time.total  = basic_stats(time_total_all);

stats.iter.stage2_AO = basic_stats(iter_stage2_AO_all);
stats.iter.stage3_AO = basic_stats(iter_stage3_AO_all);
stats.iter.stage3_FP = basic_stats(iter_stage3_FP_all);

stats.cvx.stage2 = basic_stats(cvx_calls_stage2_all);
stats.cvx.stage3 = basic_stats(cvx_calls_stage3_all);
stats.cvx.total  = basic_stats(cvx_calls_total_all);

%% 硬件信息
hardware_info = get_hardware_info();

%% 相干时间换算（ms）
coherence = finalize_coherence(coherence);

%% Table II（用于论文）
table2_stats = build_table2(time_stage1_all, time_stage2_all, time_stage3_all, time_total_all, ...
    iter_stage2_AO_all, iter_stage3_AO_all);

disp('Table II（Runtime Statistics）：');
disp(table2_stats);
fprintf('相干时间限制 T_c = %.3f ms（v=%.1f km/h, fc=%.0f GHz, T_c=%d symbols）\n', ...
    coherence.Tc_ms, coherence.v_kmh, coherence.fc/1e9, coherence.Tc_symbols);

%% 保存
save(out_file, ...
    'time_stage1_all','time_stage2_all','time_stage3_all', ...
    'iter_stage2_AO_all','iter_stage3_AO_all','iter_stage3_FP_all', ...
    'hardware_info','solver_info', ...
    'cvx_calls_stage2_all','cvx_calls_stage3_all','cvx_calls_total_all', ...
    'time_total_all','stats','table2_stats','coherence');

fprintf('已保存：%s\n', out_file);

%% ========================= 本地函数 =========================
function [t1,t2,t3,it2,it3,fp3,cvx2,cvx3,cvxT] = one_run(mc_seed, K_rice)
    % 每次 MC 的随机种子（不影响算法逻辑，仅确保可复现）
    rng(mc_seed, 'twister');
    try
        [~,~,run_stats] = RUN_OPT_Protocol(K_rice);
        t1 = run_stats.time_stage1_ms;
        t2 = run_stats.time_stage2_ms;
        t3 = run_stats.time_stage3_ms;
        it2 = run_stats.iter_stage2_AO;
        it3 = run_stats.iter_stage3_AO;
        fp3 = run_stats.iter_stage3_FP;
        cvx2 = run_stats.cvx_calls_stage2;
        cvx3 = run_stats.cvx_calls_stage3;
        cvxT = run_stats.cvx_calls_total;
    catch
        t1 = NaN; t2 = NaN; t3 = NaN;
        it2 = NaN; it3 = NaN; fp3 = NaN;
        cvx2 = NaN; cvx3 = NaN; cvxT = NaN;
    end
end

function s = basic_stats(x)
    x = x(:);
    s.avg = mean(x, 'omitnan');
    s.max = max(x, [], 'omitnan');
    s.min = min(x, [], 'omitnan');
    s.std = std(x, 0, 'omitnan');
end

function hw = get_hardware_info()
    hw = struct();

    % MATLAB 版本
    hw.MATLAB_version = version;

    % RAM（尽量获取物理内存）
    try
        m = memory;
        hw.RAM = sprintf('%.1f GB', m.PhysicalMemory.Total/1024^3);
    catch
        hw.RAM = 'unknown';
    end

    % CPU 型号（Windows 优先 wmic，失败则 unknown）
    hw.CPU = 'unknown';
    try
        [status,out] = system('wmic cpu get Name /value');
        if status == 0
            out = strtrim(out);
            % 解析：Name=...
            k = strfind(out, 'Name=');
            if ~isempty(k)
                val = strtrim(out(k(1)+5:end));
                val = regexprep(val, '[\r\n]+', '');
                if ~isempty(val), hw.CPU = val; end
            end
        end
    catch
    end
end

function coh = finalize_coherence(coh)
    c = 3e8;
    coh.lambda = c / coh.fc;
    coh.v_mps = coh.v_kmh / 3.6;
    coh.fD_Hz = coh.v_mps / coh.lambda;
    coh.Tc_sec = 0.423 / coh.fD_Hz;                 % 常用近似
    coh.symbol_rate = coh.Tc_symbols / coh.Tc_sec;  % Rs = Tc_symbols / Tc_sec
    coh.symbol_period = 1 / coh.symbol_rate;
    coh.Tc_ms = coh.Tc_sec * 1e3;
end

function tbl = build_table2(t1,t2,t3,tt, it2, it3)
    stage = {'I';'II';'III';'Total'};
    avg_t = [mean(t1,'omitnan'); mean(t2,'omitnan'); mean(t3,'omitnan'); mean(tt,'omitnan')];
    max_t = [max(t1,[],'omitnan'); max(t2,[],'omitnan'); max(t3,[],'omitnan'); max(tt,[],'omitnan')];
    std_t = [std(t1,0,'omitnan'); std(t2,0,'omitnan'); std(t3,0,'omitnan'); std(tt,0,'omitnan')];

    avg_iter_str = strings(4,1);
    avg_iter_str(1) = "N/A";
    avg_iter_str(2) = string(mean(it2,'omitnan'));
    avg_iter_str(3) = string(mean(it3,'omitnan'));
    avg_iter_str(4) = "-";

    tbl = table(stage, avg_t, max_t, std_t, avg_iter_str, ...
        'VariableNames', {'Stage','Avg_Time_ms','Max_Time_ms','Std_ms','Avg_Iterations'});
end
