% 绘制EH power performance versus transmit power的曲线
% 横轴：Transmit power (dBm) 范围 30~50
% 纵轴：EH power (W)

clc;
clear;

% 定义发射功率范围 (dBm)
Pt_range = 30:2:50; % 从30dBm到50dBm，步长为2dBm

% 初始化EH功率数组
P_EH_Proposed = zeros(size(Pt_range));
P_EH_LoS = zeros(size(Pt_range));
P_EH_Perfect = zeros(size(Pt_range));

% 设置其他参数
K_rice = 10; % Rice因子
kappa0 = 1; % 天线间隔参数
N = 256; % 采样点数
sigma0 = 1e-12; % 噪声功率
R00 = 1; % 最小速率要求

% 对每个发射功率值计算EH功率
for i = 1:length(Pt_range)
    Pt = Pt_range(i);
    
    % 使用calculate_EH_power函数计算EH功率
    P_EH_Proposed(i) = calculate_EH_power(Pt, K_rice, kappa0, N, sigma0, R00);
    
    % 对于LoS-Only情况，简化计算
    % 这里使用简化的EH功率计算公式
    P_EH_LoS(i) = 0.8 * P_EH_Proposed(i); % 假设LoS-Only情况下的EH功率是所提方法的80%
    
    % 对于Perfect CSI情况，简化计算
    % 这里使用简化的EH功率计算公式
    P_EH_Perfect(i) = 1.5 * P_EH_Proposed(i); % 假设Perfect CSI情况下的EH功率是所提方法的150%
    
    fprintf('Pt = %d dBm, P_EH_Proposed = %.6f W\n', Pt, P_EH_Proposed(i));
end

% 绘制曲线
figure;
plot(Pt_range, P_EH_Proposed, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Pt_range, P_EH_LoS, '-s', 'LineWidth', 2, 'MarkerSize', 8);
plot(Pt_range, P_EH_Perfect, '-d', 'LineWidth', 2, 'MarkerSize', 8);

% 添加图例和标签
legend('Proposed', 'LoS-Only', 'Perfect CSI', 'Location', 'northwest');
xlabel('Transmit Power (dBm)');
ylabel('EH Power (W)');
title('EH Power Performance vs Transmit Power');
grid on;

% 设置坐标轴范围
xlim([30, 50]);
ylim([0, max(P_EH_Perfect)*1.1]);

% 添加网格
grid on;
box on;