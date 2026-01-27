function [ h , e0 , e2 , theta, phi , theta_scatterer , phi_scatterer , IOTA , eta , Omega , h_r , h_c , h_sen , h_x_d , h_x_u , h_y_l , h_y_r , Omega0 ] = Channel_Generation_init( K,B,M,Mx,My,M_x_r , M_y_c ,normal_vector,K_rice,Coor_Ele,Coor_Ele_r,Coor_Ele_c,lambda , Rot_Matrix , Coor_Ele_x_d , Coor_Ele_x_u , Coor_Ele_y_l , Coor_Ele_y_r )

    dx = lambda/4;   %x轴方向天线间隔
    dy = lambda/4;  %y轴方向天线间隔
    delta_x = ( 2 * [0:(Mx-1)].' - Mx + 1 )/2; delta_y = ( 2 * [0:( My-1)].' - My + 1 )/2;
%% 针对初始6DMA位置，生成无线信道
PL    = 32.4 + 17.3 * log10( 5 ) + 20 * log10( 30 ) - 25 ; % Path-loss,包含20dB天线增益
Omega = 10^( -PL/10 ) * ones(K,B);  % path-loss
% Omega = ones(K,B);
IOTA  = 3;         % 散射体数量
% theta =  rand(K,1) * 2 * pi; phi = rand(K,1) * 2 * pi; % K个用户的水平角和俯仰角
% theta =  (rand(K,1)-0.5) * pi; phi = rand(K,1) * 2 * pi; % K个用户的水平角和俯仰角
% theta = 5.263863993147139 ; phi =  1.057064847365129;
% theta_scatterer = rand(K,IOTA) * 2 * pi; phi_scatterer = rand(K,IOTA) * 2 * pi; % K个用户的水平角和俯仰角
% NLoS链路信道
% for k = 1:K
%     for iota = 1:IOTA
%         eta(k,iota) = 1/sqrt(IOTA) * sqrt( Omega(k) ) * 1/sqrt(2) * ( randn(1) + 1j*randn(1) );
%     end
% end

% theta = [1.45524585485558
% 1.41504973364835
% 2.35835054090831
% 5.43852179143626];
% phi = [5.83563651866559
% 2.02237603845018
% 0.0488728727564714
% 4.43117058629477];
% theta_scatterer = [1.56165821938995	0.0891024192032021	4.97921394752858
% 2.21108342743673	0.787690933376525	3.74825942684378
% 4.40679921072753	0.638612512771587	3.13896018436890
% 3.63580537416366	6.14856924227454	5.07064531671004];
% phi_scatterer = [2.69668281857476	1.21137536400753	4.64457857472923
% 6.11843287903727	2.30257164625589	2.92921398237037
% 4.01675984308597	2.54356244257540	1.59404072941367
% 5.60346946308695	1.59961515297204	1.83698635921796];
% 
% eta = [0.000162840258005646 - 0.00410263506222009i	0.00141051775190082 - 0.00356861079602921i	0.00101196371517660 + 6.34354465897823e-05i
% 1.32018770279289e-05 - 0.00116493108313416i	0.000188490299469544 + 0.00190905595007172i	-0.00300202684589278 + 0.000658481193069772i
% 0.00104360111647310 + 5.04415928924031e-05i	0.000150090015563367 - 0.000851331624431704i	0.00205548167310716 - 0.000936050547026122i
% 0.000794415601009792 - 0.000778357233264543i	0.00120470376911885 - 0.00126613881231049i	-0.00380369315367250 - 0.00117017556283066i];


% Loc_Rec = 2 * [     -0.455360759292362   0.936069924833439  -2.257524203591861
%   -0.474371387063496  -2.254673933889913   0.915358234303231
%    1.203201306715649  -2.450389020600264  -1.935736297706728];
%  
% Loc_Sca(:,:,1) = [-3.18000793268770	-3.10480777675081	-1.75813656859536
% -2.22762077650595	-1.58711193668651	-2.18100312792807
% 1.41121921033561	-1.38433345774810	0.968648439403465];
% 
% Loc_Sca(:,:,2) = [1.53388933190742	-1.51758114957103	3.82361104918924
% -2.14729645599737	3.71437679759213	2.35791856724873
% 3.55425228612239	-0.726824778327774	3.11436511974883];
% 
% Loc_Sca(:,:,3) = [-3.42822583920652	-0.589331358795705	-2.32028205337480
% -3.50954175334392	1.15399752609044	-2.70361577060424
% -3.08483937805276	-3.98908857542391	-2.66434659940792];
% 
% eta = [ -0.042298152048777 + 0.435197206553587i -0.321655014542619 + 0.356156309166077i  0.461267624548972 + 0.955378222632790i
%   0.794253652750950 - 0.140875686998121i -0.099973739973389 + 0.406673088435238i -0.538712370725737 + 0.614602815141328i
%   0.186727093751285 + 0.473243221941504i -0.587982556918943 + 0.603481476938166i -0.227080390951307 + 0.247357707423056i];

Loc_Rec = [
   -0.809152890890314   2.400140296822853  -1.970621624264194
   2.965467557427142  -1.827964708554828  -0.197932110191203
   1.061295612592333   1.623630544975729  -1.040407270208264
   1.094988407791883   1.991820259910715   1.298719792580405];
Loc_Sca(:,:,1) = [   3.307006849112155   1.058873969803276  -3.219676760004724
   3.660054683474381   3.719108281594212  -2.739095346579614
  -0.116994810217270   2.402243751110401  -2.864909290982277]; 
Loc_Sca(:,:,2) = [  -3.714306571406484   2.793034446950217   3.471945982060404
   1.945059744999329  -0.862183843726655   1.243823121420453
  -3.745337228980635  -1.784616120312880  -3.630628874950768];
Loc_Sca(:,:,3) = [  -3.724431355976730  -0.490045122748814  -0.947532343255933
  -2.505019163564971  -0.081884833694152  -0.435310394312804
   2.037493455858887  -1.791799384011373   1.437621414829398];
Loc_Sca(:,:,4) = [  -1.276914186670934   0.682142007838219  -2.209504484070904
   0.047656413321139   1.592613781253488   3.127226020286388
  -2.891004457370567  -2.805647955527540  -1.939933967010108];
eta = [-0.177013201865064 + 0.139875852735497i	1.23900231681864 + 0.296145034729511i	-0.0506817179300599 + 0.608166501686522i
0.199590045948221 + 0.422411652453875i	-0.321406858746061 + 0.362685997920621i	-1.20199897577424 + 0.587216295578675i
0.130315605568860 + 0.127723987234153i	0.256260427045311 + 0.446323840408201i	0.245236511313696 - 0.495661207479030i
0.456158733510309 - 0.444608636660307i	0.630421866453449 + 0.0350812382164900i	0.959570143507271 - 0.251318415715197i];
Omega0 = 1e-4 * [   0.092390160651137   0.078648952175001   0.108924279234190   0.139836456218170
   0.089926964433591   0.073542831246887   0.089292024835166   0.122363802402645
   0.075988029097874   0.060805743593390   0.078900123019043   0.089111896788634
   0.085198940615822   0.063896774310169   0.074571957476659   0.086591282450251
];
%% 生成无线信道
% LoS链路信道
for k = 1:K
    theta(k) = atan2(Loc_Rec(k,2) , Loc_Rec(k,1));
    phi(k) = atan2( Loc_Rec(k,3) , sqrt( Loc_Rec(k,1)^2 + Loc_Rec(k,2)^2 ) );
    for b = 1:B
        h{k,b} = sparse(M,1);%%%%%%%%%%%%%%%%%%%%%注意
%             e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            e0{k} = Loc_Rec(k,:).' / norm( Loc_Rec(k,:) , 2 );
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h{k,b} = sparse( h{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda_LoS(k,b) * exp( 1j * 2*pi/lambda * Coor_Ele{b} * e0{k} ) );
    end
end

%% 测试
for b = 1:B
    Px{b}    = dx * Rot_Matrix{b} * [1;0;0];
    Py{b}    = dy * Rot_Matrix{b} * [0;1;0];
    T{b}     = exp( 1j * 2*pi/lambda * ( kron( delta_x * ( Px{b}.' * e0{1} ) , ones(1,My) ) + kron( ones(Mx,1) , delta_y.' * ( Py{b}.' * e0{1} ) ) ) );

end








for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
            theta_scatterer(k,iota) = atan2( Loc_Sca(iota,2,k) , Loc_Sca(iota,1,k));
            phi_scatterer(k,iota)   = atan2( Loc_Sca(iota,3,k) , sqrt( Loc_Sca(iota,1,k)^2 + Loc_Sca(iota,2,k)^2 ) );


            % eta(k,iota,b) = 1/sqrt(IOTA) * Omega(k,b) * 1/sqrt(2) * ( randn(1) + 1j*randn(1) );
%             e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            e1{k,iota} = Loc_Sca(iota,:,k).' / norm( Loc_Sca(iota,:,k) , 2 );
%           e = Rot_Matrix{b}.' * e;
            if dot( normal_vector{b} , e1{k,iota} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h{k,b} = sparse( h{k,b} + sqrt(1/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda(k,iota,b) * eta(k,iota) * exp( 1j * 2*pi/lambda * Coor_Ele{b} * e1{k,iota} ) );
        end
    end
end
% for k = 1:K
%     H{k} = [];
%     for b = 1:B
%         H{k} = [H{k} ; h{k,b} ];
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第b个RHS的"行"感知元素元素的信道
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_r{k,b} = sparse(M_x_r,1);%%%%%%%%%%%%%%%%%%%%%注意
%             e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_r{k,b} = sparse( h_r{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda_LoS(k,b) * exp( 1j * 2*pi/lambda * Coor_Ele_r{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
%             e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            e2{k,iota} = e1{k,iota};
            if dot( normal_vector{b} , e1{k,iota} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_r{k,b} = sparse( h_r{k,b} + sqrt(1/(K_rice+1)) * Lambda(k,iota,b) * sqrt(Omega0(k,b))  * eta(k,iota) * exp( 1j * 2*pi/lambda * Coor_Ele_r{b} * e1{k,iota} ) );
        end
    end
end

%% 第b个RHS的"列"感知元素元素的信道
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_c{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
%             e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_c{k,b} = sparse( h_c{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda_LoS(k,b) * exp( 1j * 2*pi/lambda * Coor_Ele_c{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
%             e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,iota} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_c{k,b} = sparse( h_c{k,b} + sqrt(1/(K_rice+1)) * Lambda(k,iota,b) * sqrt(Omega0(k,b))  * eta(k,iota) * exp( 1j * 2*pi/lambda * Coor_Ele_c{b} * e1{k,iota} ) );
        end
    end
end

%% x,d
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_x_d{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
%             e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_x_d{k,b} = sparse( h_x_d{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda_LoS(k,b) * exp( 1j * 2*pi/lambda * Coor_Ele_x_d{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
%             e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,iota} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_x_d{k,b} = sparse( h_x_d{k,b} + sqrt(1/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda(k,iota,b) * eta(k,iota) * exp( 1j * 2*pi/lambda * Coor_Ele_x_d{b} * e1{k,iota} ) );
        end
    end
end


%% x,u
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_x_u{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
%             e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_x_u{k,b} = sparse( h_x_u{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda_LoS(k,b) * exp( 1j * 2*pi/lambda * Coor_Ele_x_u{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
%             e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,iota} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_x_u{k,b} = sparse( h_x_u{k,b} + sqrt(1/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda(k,iota,b) * eta(k,iota) * exp( 1j * 2*pi/lambda * Coor_Ele_x_u{b} * e1{k,iota} ) );
        end
    end
end

%% y,l
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_y_l{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
%             e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_y_l{k,b} = sparse( h_y_l{k,b} + sqrt(K_rice/(K_rice+1)) * sqrt(Omega0(k,b)) * Lambda_LoS(k,b) * exp( 1j * 2*pi/lambda * Coor_Ele_y_l{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
%             e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,iota} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_y_l{k,b} = sparse( h_y_l{k,b} + sqrt(1/(K_rice+1)) *  sqrt(Omega0(k,b)) * Lambda(k,iota,b) * eta(k,iota) * exp( 1j * 2*pi/lambda * Coor_Ele_y_l{b} * e1{k,iota} ) );
        end
    end
end


%% y,r
% LoS链路信道
for k = 1:K
    for b = 1:B
        h_y_r{k,b} = sparse(M_y_c,1);%%%%%%%%%%%%%%%%%%%%%注意
%             e0{k} = [ sin(phi(k)) * cos(theta(k)) ; sin(phi(k)) * sin(theta(k)) ; cos(phi(k)) ];
            % e = Rot_Matrix{b}.' * e0;
            if dot( normal_vector{b} , e0{k} ) > 0
                Lambda_LoS(k,b) = 1;
            else
                Lambda_LoS(k,b) = 0;
            end
            h_y_r{k,b} = sparse( h_y_r{k,b} + sqrt(K_rice/(K_rice+1))  * sqrt(Omega0(k,b)) * Lambda_LoS(k,b) * exp( 1j * 2*pi/lambda * Coor_Ele_y_r{b} * e0{k} ) );
    end
end
% NLoS链路信道
for k = 1:K
    for b = 1:B
        % h{k,b} = 0;%%%%%%%%%%%%%%%%%%%%%注意
        for iota = 1:IOTA
%             e1{k,b} = [ sin(phi_scatterer(k,iota)) * cos(theta_scatterer(k,iota)) ; sin(phi_scatterer(k,iota)) * sin(theta_scatterer(k,iota)) ; cos(phi_scatterer(k,iota)) ];
            if dot( normal_vector{b} , e1{k,iota} ) > 0
                Lambda(k,iota,b) = 1;
            else
                Lambda(k,iota,b) = 0;
            end
            h_y_r{k,b} = sparse( h_y_r{k,b} + sqrt(1/(K_rice+1)) *  sqrt(Omega0(k,b)) * Lambda(k,iota,b) * eta(k,iota) *  exp( 1j * 2*pi/lambda * Coor_Ele_y_r{b} * e1{k,iota} ) );
        end
    end
end


% 把四个合并起来
for k = 1:K
    for b = 1:B
        h_sen{k,b} = zeros( M_x_r , M_y_c );
        h_sen{k,b}(1,:) = h_x_d{k,b}.';
        h_sen{k,b}(end,:) = h_x_u{k,b}.';
        h_sen{k,b}(:,1) = h_y_l{k,b};
        h_sen{k,b}(:,end) = h_y_r{k,b};
    end
end

end