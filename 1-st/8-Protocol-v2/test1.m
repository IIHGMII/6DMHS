% 原始法向量 a、新法向量 b
a = [0; 0; 1];  % 原来法向是 z+
b = [1; 1; 1];  % 新的法向
a = a / norm(a);
b = b / norm(b);

% 计算旋转轴和角度
v = cross(a, b);
s = norm(v);
c = dot(a, b);

if s < 1e-12
    % 特殊情况：平行或反平行
    if c > 0
        R = eye(3);  % 平行，不旋转
    else
        % 反平行：取任意垂直向量做旋转轴
        perp = null(a.'); % 任意一个与 a 垂直的单位向量
        k = perp(:,1);
        theta = pi;
        K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
        R = eye(3)*cos(theta) + (1-cos(theta))*(k*k.') + K*sin(theta);
    end
else
    k = v / s;
    theta = atan2(s, c);
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3)*cos(theta) + (1-cos(theta))*(k*k.') + K*sin(theta);
end

% 测试：旋转后的 a 是否等于 b
a_rot = R * a;
disp([a_rot b])

% 旋转平面上的所有点
% 例如：原来平面是 z=0，正方形顶点：
square_pts = [-1 -1 0; 1 -1 0; 1 1 0; -1 1 0]';
square_rot = R * square_pts;

disp('旋转后的顶点坐标：')
disp(square_rot')
