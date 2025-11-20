function [xmin, fmin, iter] = B_goldenSectionSearch(func, a, b, tol)
    % func: 目标函数
    % a, b: 初始搜索区间
    % tol: 精度阈值

    % 定义黄金分割比率
    ratio = (sqrt(5) - 1) / 2;

    % 初始化
    c = b - ratio * (b - a);
    d = a + ratio * (b - a);
    fc = func(c);
    fd = func(d);
    iter = 0;

    % 主循环
    while abs(b - a) > tol
        iter = iter + 1;
        if fc < fd
            b = d;
            d = c;
            fd = fc;
            c = b - ratio * (b - a);
            fc = func(c);
        else
            a = c;
            c = d;
            fc = fd;
            d = a + ratio * (b - a);
            fd = func(d);
        end
    end

    % 确定最小值位置和函数值
    xmin = (a + b) / 2;
    fmin = func(xmin);

    % fprintf('Minimum occurs at x = %.4f\n', xmin);
    % fprintf('f(x) at minimum is: %.4f\n', fmin);
    % fprintf('Total iterations: %d\n', iter);
end
