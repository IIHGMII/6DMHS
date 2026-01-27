function [X, Y, Z] = Orientation_uniformSpherePoints(N)
%% 输入参数：
% N: 空间球面采样点数
%% 输出参数
% 空间采样点坐标(X,Y,Z)
%% 方法
% 利用黄金分割法来生成采样点
    indices = 0:N-1;  % Array of indices
    phi = pi * (3 - sqrt(5));  % Golden angle in radians
    
    % Spherical coordinates
    z = 1 - (2 * indices / (N - 1));  % Heights are evenly spaced
    radius = sqrt(1 - z.^2);  % Radius at each height
    theta = phi * indices;  % Golden angle increment
    
    % Conversion to Cartesian coordinates
    X = radius .* cos(theta);
    Y = radius .* sin(theta);
    Z = z;
    
    % Plotting the points for visualization
    % figure;
    % scatter3(X, Y, Z, 'filled');
    % axis equal;
    % title('Uniform Distribution of Points on Sphere');
    % xlabel('X');
    % ylabel('Y');
    % zlabel('Z');
end