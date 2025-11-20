%_______________________________________________________________________________________%
%  Aquila Optimizer (AO) source codes (version 1.0)                                     %
%                                                                                       %
%  Developed in MATLAB R2015a (7.13)                                                    %
%  Author and programmer: Laith Abualigah                                               %
%  Abualigah, L, Yousri, D, Abd Elaziz, M, Ewees, A, Al-qaness, M, Gandomi, A.          %
%         e-Mail: Aligah.2020@gmail.com                                                 %
%       Homepage:                                                                       %
%         1- https://scholar.google.com/citations?user=39g8fyoAAAAJ&hl=en               %
%         2- https://www.researchgate.net/profile/Laith_Abualigah                       %
%                                                                                       %
%   Main paper:                                                                         %
%_____________Aquila Optimizer: A novel meta-heuristic optimization algorithm___________%
%__Main paper: please, cite it as follws:_______________________________________________%

%Abualigah, L., Yousri, D., Elaziz, M.A., Ewees, A.A., A. Al-qaness, M.A., Gandomi, A.H.%
%Aquila Optimizer: A novel meta-heuristic optimization Algorithm,          
%Computers & Industrial Engineering (2021), doi:https://doi.org/10.1016/j.cie.2021.107250 
% Computers & Industrial Engineering.
%_______________________________________________________________________________________%

clear all ;
clc ;
NUM = 200000;
Solution_no = 4;  
F_name = 'F24';    
M_Iter = NUM;    

[LB,UB,Dim,F_obj]=Get_F(F_name); 
[Best_FF,Best_P,conv]=AO(Solution_no,M_Iter,LB,UB,Dim,F_obj);  

display(['The best-obtained solution by AO is : ', num2str(Best_P)]);
display(['The best optimal values of the objective funciton found by AO is : ', num2str(Best_FF)]);



figure;plot( 1 : NUM , conv );

%% 计算波束增益
fc = 30e9 ;
c  = 3e8 ;
lambda = c/fc ; dx = lambda / 4 ; dy = lambda / 4 ;
Q = 1 ;
Mx = 64 ; My = 8 ; M = Mx * My ;
D  = Mx * dx*1.2;
gamma = sqrt(3) ;
% b = @(theta,phi) 1/sqrt(M) * exp( 1j*2*pi/lambda * ( kron([0:Mx-1],ones(1,My)).'*theta*dx + kron(ones(1,Mx),[0:My-1]).'*phi*dy ) );

% Coor_Feed = [0.02 ; 0.08 ; 0];
Coor_Ele  = [ kron([0:Mx-1],ones(1,My)) * dx ; kron(ones(1,Mx),[0:My-1])*dy ; zeros(1,M) ];
% for i = 1:M
%     D(i,1) = norm( Coor_Ele(:,i) - Coor_Feed(Q) , 2 );
% end

% Best_P(3) = 0.16608; Best_P(4) = 0.002617;
% Best_P(3) = 0; Best_P(4) = 0;
%% 1个馈源
% o = @(theta,phi) - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( sin(theta*pi/2) * cos(phi*pi)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(theta*pi/2) * sin(phi*pi) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - ( Best_P(3)  * D / 2 )^2 - ( Best_P(4) * D / 2 )^2 ) ).' )  ) ).^2  ) ) ;
o = @(theta,phi) - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( cos(theta*pi) * cos(phi*pi/2)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(theta*pi) * cos(phi*pi/2) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - ( Best_P(3)  * D / 2 )^2 - ( Best_P(4) * D / 2 )^2 ) ).' )  ) ).^2  ) ) ;

%% 2个馈源
% o = @(theta,phi) - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( sin(0.95 * theta*pi/2) * cos(0.95 * phi*pi)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(0.95 * theta*pi/2) * sin(0.95 * phi*pi) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - (  Best_P(3) * D / 2 )^2 - ( Best_P(4) * D / 2 )^2 ) ).' )  ) ).^2  ) ) ...
%                  - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( sin(0.95 * theta*pi/2) * cos(0.95 * phi*pi)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(0.95 * theta*pi/2) * sin(0.95 * phi*pi) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - (  Best_P(5) * D / 2 )^2 - ( Best_P(6) * D / 2 )^2 ) ).' )  ) ).^2  ) );
%% 4个馈源
% o = @(theta,phi) - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( sin(0.95 * theta*pi/2) * cos(0.95 * phi*pi)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(0.95 * theta*pi/2) * sin(0.95 * phi*pi) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - ( Best_P(3) * D / 2 )^2 - ( Best_P(4) * D / 2 )^2 ) ).' )  ) ).^2  ) ) ...
%                  - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( sin(0.95 * theta*pi/2) * cos(0.95 * phi*pi)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(0.95 * theta*pi/2) * sin(0.95 * phi*pi) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - ( Best_P(5) * D / 2 )^2 - ( Best_P(6) * D / 2 )^2 ) ).' )  ) ).^2  ) ) ...
%                  - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( sin(0.95 * theta*pi/2) * cos(0.95 * phi*pi)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(0.95 * theta*pi/2) * sin(0.95 * phi*pi) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - ( Best_P(7) * D / 2 )^2 - ( Best_P(8) * D / 2 )^2 ) ).' )  ) ).^2  ) ) ...
%                  - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( sin(0.95 * theta*pi/2) * cos(0.95 * phi*pi)  )*dx + kron(ones(1,Mx),[0:My-1]).'*( sin(0.95 * theta*pi/2) * sin(0.95 * phi*pi) )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - ( Best_P(9) * D / 2 )^2 - ( Best_P(10) * D / 2 )^2 ) ).' )  ) ).^2  ) );

%% 1个馈源，馈源位置
% [ Best_P(3)  * D / 2 ,  Best_P(4)  * D / 2]

%% 2个馈源，馈源位置
% [( Best_P(3) )  * D / 2 ,  Best_P(4)    * D / 2]
% [( Best_P(5) )  * D / 2 ,  Best_P(6)    * D / 2]
%% 4个馈源，馈源位置
% [ Best_P(3)  * D / 2 ,  Best_P(4)  * D / 2]
% [ Best_P(5)  * D / 2 ,  Best_P(6)  * D / 2]
% [ Best_P(7)  * D / 2 ,  Best_P(8)  * D / 2]
% [ Best_P(9)  * D / 2 ,  Best_P(10)  * D / 2]


o1 = @(theta,phi) - 1/4/sqrt(M) * abs( sum( ( 1 + exp( -1j * 2*pi/lambda * (( kron([0:Mx-1],ones(1,My)).'*( theta  )*dx + kron(ones(1,Mx),[0:My-1]).'*( phi )*dy ) - gamma * ( sqrt( Coor_Ele(1,:).^2 + Coor_Ele(2,:).^2 - ( Best_P(3)  * D / 2 )^2 - ( Best_P(4) * D / 2 )^2 ) ).' )  ) ).^2  ) ) ;

num1=1;
for theta0 = -1:0.005:1
    num2 = 1;
    for phi0 = -1:0.005:1
        if theta0^2 + phi0^2 <= 1
            Beam_Gain(num1,num2) = -o1(theta0,phi0);
        else
            Beam_Gain(num1,num2) = nan;
            
        end
        num2 = num2 + 1;
    end
    num1 = num1 + 1;
end


x = linspace(-1, 1, size(Beam_Gain, 2));
y = linspace(-1, 1, size(Beam_Gain, 1));

[X, Y] = meshgrid( x , y );
figure(1);
mesh(X,Y,Beam_Gain );  % 绘制矩阵 A 的表面图，不显示网格线
view(2);            % 将视图设置为二维
colorbar;           % 显示颜色条，颜色表示值的大小
% title('2D Heatmap Style Plot of Matrix A using surf');  % 添加标题
% xlabel('Column Index');  % X轴标签
% ylabel('Row Index');     % Y轴标签
colormap turbo;            % 设置颜色映射，'jet' 是一种流行的颜色映射
grid off;
[maxVal, ind] = max(Beam_Gain(:));
[row, col] = ind2sub(size(Beam_Gain), ind);
hold on;  % 保持当前图像，以便添加标记
plot3(x(col), y(row), maxVal, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'k');  % 在最大值位置放置一个黑色的标记
hold off;
axis([-1 1 -1 1]);
box on;
xlabel('$f_1(\theta,\phi)$' , 'FontSize' , 12 , 'FontName' , 'timesnewroman' , 'Interpreter' , 'latex' );
ylabel('$f_2(\theta,\phi)$' , 'FontSize' , 12 ,'FontName' , 'timesnewroman' , 'Interpreter' , 'latex' );
hText = text( 0 , 0 , 0, ['(0.16,-0.273)'],'HorizontalAlignment','center', 'VerticalAlignment','top' );
uistack(hText, 'top');






% num1 = 1; num2 = 1;
% for theta0 = -1:0.002:1
%     num2 = 1;
%     for phi0 = -1:0.004:1
% %         if theta0^2 + phi0^2 <= 1
%             Beam_Gain(num1,num2) = ( -o(theta0,phi0) );
% %         else
% %             Beam_Gain(num1,num2) = 0;
% %         end
%         num2 = num2 + 1;
%     end
%     num1 = num1 + 1;
% end
% 
% % 创建X和Y坐标网格
% x = linspace(-1*pi/2, 1*pi/2, size(Beam_Gain, 2));
% y = linspace(-1*pi, 1*pi, size(Beam_Gain, 1));
% 
% [X, Y] = meshgrid( x , y );
% 
% % 绘制矩阵 A 的三维表面图，作为二维热图使用
% figure;             % 创建新的图形窗口
% % surf(X,Y,Beam_Gain, 'EdgeColor', 'none');  % 绘制矩阵 A 的表面图，不显示网格线
% mesh(X,Y,Beam_Gain );  % 绘制矩阵 A 的表面图，不显示网格线
% view(2);            % 将视图设置为二维
% colorbar;           % 显示颜色条，颜色表示值的大小
% % title('2D Heatmap Style Plot of Matrix A using surf');  % 添加标题
% % xlabel('Column Index');  % X轴标签
% % ylabel('Row Index');     % Y轴标签
% colormap jet;            % 设置颜色映射，'jet' 是一种流行的颜色映射
% 
% % 找到并标记最大值的位置
% [maxVal, ind] = max(Beam_Gain(:));
% [row, col] = ind2sub(size(Beam_Gain), ind);
% hold on;  % 保持当前图像，以便添加标记
% plot3(x(col), y(row), maxVal, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'k');  % 在最大值位置放置一个黑色的标记
% hold off;
% axis([-pi/2 pi/2 -pi pi]);
% box on;
% xlabel('The elevation angle $\theta$' , 'FontSize' , 12 , 'FontName' , 'timesnewroman' , 'Interpreter' , 'latex' );
% ylabel('The azimuth angle $\phi$' , 'FontSize' , 12 ,'FontName' , 'timesnewroman' , 'Interpreter' , 'latex' );
% text( 1 , 1 , 0, '(-0.34,-0.73)' );
% % uistack(hText, 'top');
% hText = text(0, 0, 0, 'This is the top text', ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
%     'FontSize', 12, 'Color', 'red');
% uistack(hText, 'top');



























% theta = y(row); phi = x(col);
% n1 = [ sin(0.95 * theta * pi/2) * cos(0.95 * phi * pi) ; sin(0.95 * theta * pi/2) * sin(0.95 * phi*pi) ; cos(0.95 * theta*pi/2) ]


