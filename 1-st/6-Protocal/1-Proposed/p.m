% 参数设置
N = 256;       % 信号长度
M = 80;        % 观测数

% 生成稀疏信道向量，其中包含一个强LoS成分和几个较弱的多径成分
x = zeros(N, 1);
x(1) = 10;  % 假设LoS成分在第一个位置，幅度较大
q = randperm(N, 20);
x(q(2:5)) = randn(4, 1) + 1i*randn(4, 1);  % 添加一些较弱的多径成分

% 生成观测矩阵（高斯随机矩阵）
Phi = randn(M, N) + 1i*randn(M, N);  % 对于复数信道

% 生成观测数据
y = Phi * x;

% OMP算法恢复信号
x_est = omp(Phi, y, 5);  % 假定我们知道存在5个显著的路径

% 显示结果
fprintf('原始信号LoS成分：%f + %fi\n', real(x(1)), imag(x(1)));
fprintf('恢复信号LoS成分：%f + %fi\n', real(x_est(1)), imag(x_est(1)));

% 函数定义：OMP算法
function x_est = omp(Phi, y, K)
    residual = y;  % 初始化残差
    support = [];  % 支撑集
    x_est = zeros(size(Phi, 2), 1);

    for iter = 1:K
        % 计算投影
        projections = Phi' * residual;
        
        % 更新支撑集
        [~, idx] = max(abs(projections));
        support = unique([support, idx]);
        
        % 最小二乘估计
        x_temp = zeros(size(x_est));
        x_temp(support) = Phi(:, support) \ y;
        
        % 更新残差
        residual = y - Phi(:, support) * x_temp(support);
        
        % 更新估计
        x_est = x_temp;
    end
end
