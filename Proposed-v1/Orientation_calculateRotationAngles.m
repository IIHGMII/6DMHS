function [alpha, beta, gamma] = Orientation_calculateRotationAngles(x0, y0, z0)
%% 输入参数：
% 全局坐标系中的方向向量(x0,y0,z0)
%% 输出参数：
% 将向量(0,0,1)旋转到(x0,y0,z0),所需的旋转向量(alpha,beta,gamma)
%% 归一化方向向量(x0,y0,z0)
    norm_v = sqrt(x0^2 + y0^2 + z0^2);
    x0 = x0 / norm_v;
    y0 = y0 / norm_v;
    z0 = z0 / norm_v;
%% 计算旋转向量(x0,y0,z0)
    % Calculate angles
    gamma = atan2(y0, x0);  % Around z-axis
    beta = acos(z0);        % Around y-axis after z-rotation
    alpha = 0;              % Around x-axis, typically 0 since no further rotation is needed to align z-axis

    % Display the angles in radians
    % disp(['Alpha (around x-axis): ', num2str(alpha), ' radians']);
    % disp(['Beta (around y-axis): ', num2str(beta), ' radians']);
    % disp(['Gamma (around z-axis): ', num2str(gamma), ' radians']);
end
