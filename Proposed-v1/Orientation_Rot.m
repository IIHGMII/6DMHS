function R = Orientation_Rot(alpha,beta,gamma)
%% 输入参数：
% alpha: 沿x轴旋转角度，正视图方向逆时针为正
% beta : 沿y轴旋转角度，左视图方向顺时针为正
% gamma: 沿z轴旋转角度，俯视图方向逆时针为正
%% 输出参数
% 旋转矩阵R
%%
R = [ cos(gamma) * cos(beta), cos(gamma)*cos(beta) - sin(gamma)*cos(alpha), cos(gamma)*sin(beta)*cos(alpha)+sin(gamma)*sin(alpha) ; ...
      sin(gamma)*cos(beta), sin(gamma)*sin(beta)*sin(alpha) + cos(gamma)*cos(alpha), sin(gamma)*sin(beta)*cos(alpha)-cos(gamma)*sin(alpha) ; ...
      -sin(beta), cos(beta)*sin(alpha), cos(beta)*cos(alpha) ];

end