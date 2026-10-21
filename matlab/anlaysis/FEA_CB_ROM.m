%%
%系统模型信息
stlFile = 'old_wafer_matlab.STL';
gm  = fegeometry("old_wafer_matlab.stl");

act_loadedVertex = [0.0,0.0,0.0;
                    0.15*sqrt(2)/2,-0.15*sqrt(2)/2,0;
                    -0.15*sqrt(2)/2,-0.15*sqrt(2)/2,0;
                    0,0.15,0];%输入N*3维数的矩阵，N个坐标
[gm,act_loadedVertex] = ...
    addVertex(gm,Coordinates=act_loadedVertex);

model = femodel(AnalysisType="structuralModal", ...
                Geometry = gm);

model.MaterialProperties = ...
    materialProperties(YoungsModulus=70e9, ...
                       PoissonsRatio=0.33, ...
                       MassDensity=2700);

model = generateMesh(model,Hmin=0.02);
% figure;
% pdemesh(model);
% title("Mesh with Quadratic Tetrahedral Elements");

% 
% result_modal = solve(model,FrequencyRange=[-0.1 1e5]);
% 
% natuarlfre_rad = result_modal.NaturalFrequencies;
% natuarlfre_Hz = natuarlfre_rad/2/pi;
% 
% figure;
% for i = 1:8
%     subplot(4,2,i);
%     pdeplot3D(result_modal.Mesh,ColorMapData=result_modal.ModeShapes.uz(:,i+6));
%     axis equal
%     title(sprintf(['Mode=%d, z-displacement\n', ...
%     'Frequency(Hz): FEM=%g'], ...
%     i,result_modal.NaturalFrequencies(i+6)/2/pi));
% end

FEM = assembleFEMatrices(model);

M_full = FEM.M;    % 3N x 3N
K_full = FEM.K;    
M_full = full(M_full);
K_full = full(K_full);
maxModesSolve = 20;   % 用 eigs 求解时最多求的特征对数（按问题规模调整）

[phai,W2] = eigs(K_full,M_full,maxModesSolve, 'SM');
f = sqrt(abs(diag(W2))) / (2*pi);
nflex = 14;
nrigid = 6;

Tn = phai(:,1:nrigid+nflex);

node_coodinfo = model.Geometry.Mesh.Nodes;
numNodes = size(node_coodinfo,2);
findnodetol = 1e-6;
bottomfacenode_idx = find(abs(node_coodinfo(3,:)) < findnodetol);
bottomnode_coodinfo = node_coodinfo(:,bottomfacenode_idx);


% bottomnode_coodinfo: 3 × N (you already have this)
xy = bottomnode_coodinfo(1:2,:).';   % N × 2, 转置便于排序

% 先 x 后 y 排序
[~, idx_sort] = sortrows(xy, [1 2]);

% 得到重新编号后的底面节点
bottomnode_coodinfo_sorted = bottomnode_coodinfo(:, idx_sort);

% 同时对节点索引和自由度索引一致重排
bottomfacenode_idx_sorted = bottomfacenode_idx(idx_sort);

% 工具函数保持不变
getDoF = @(x,numNodes) [x(:); x(:) + numNodes; x(:) + 2*numNodes];

bottomfaceDoFs_sorted = getDoF(bottomfacenode_idx_sorted, numNodes);
bottomfacezdof_sorted = bottomfacenode_idx_sorted(:) + 2*numNodes;
Tn_bottomface_sorted = Tn(bottomfacezdof_sorted,:);


% 
%===== 绘制底面节点 xy 散点并标号 =====
% figure;
% plot(bottomnode_coodinfo_sorted(1,:), bottomnode_coodinfo_sorted(2,:), 'o', 'MarkerSize', 4);
% hold on; axis equal; grid on;
% title('Bottom Face Nodes (sorted by x then y)');
% xlabel('X'); ylabel('Y');
% 
% % 标号
% for k = 1:length(bottomfacenode_idx_sorted)
%     text(bottomnode_coodinfo_sorted(1,k), ...
%          bottomnode_coodinfo_sorted(2,k), ...
%          sprintf('%d',k), ...      % 标号k
%          'FontSize', 6, ...
%          'VerticalAlignment','bottom', ...
%          'HorizontalAlignment','right');
% end
% 
% hold off;

%Craig–Bampton降阶方法

retained = bottomfacezdof_sorted(:);
allDOF = (1:size(K_full,1))';
interior = setdiff(allDOF, retained);

Kii = K_full(interior, interior);
Kib = K_full(interior, retained);
Kbi = K_full(retained, interior);
Kbb = K_full(retained, retained);

Mii = M_full(interior, interior);
Mib = M_full(interior, retained);
Mbb = M_full(retained, retained);






[phai_i_CB, W2_i_CB] = eigs(Kii, Mii, maxModesSolve, 'SM');  % 固定接口内模
f_i_CB = sqrt(abs(diag(W2_i_CB)))/2/pi;

% 约束模
Psi = - (Kii \ Kib);   % interior x interface  (每列为 interior 上的约束模式)

% 构造 T （size: (Ni+Nb) x (nmodes + Nb) ）
T_top = [phai_i_CB, Psi];          % interior rows
T_bottom = [zeros(length(retained), maxModesSolve), eye(length(retained))]; % interface rows
T = [T_top; T_bottom];

% 投影（注意 K_full 行列顺序必须是 [interior; interface]）
perm = [interior; retained];
K_perm = K_full(perm, perm);
M_perm = M_full(perm, perm);
[phai_sort,W2_sort] = eigs(K_perm,M_perm,maxModesSolve, 'SM');
f_sort = sqrt(abs(diag(W2_sort))) / (2*pi);
%将phai_full按照[interior; retained]顺序重排观察相关度
phai_full = phai(perm,:);
for i = 7:maxModesSolve
    MAC_full1(i) = (phai_full(:,i)'*phai_sort(:,i))^2/...
        ((phai_full(:,i)'*phai_full(:,i))*(phai_sort(:,i)'*phai_sort(:,i)));
end
aa1 = phai_sort'*M_perm*phai_sort;
aa2 = phai_full'*M_perm*phai_full;

K_red = T' * K_perm * T;
M_red = T' * M_perm * T;

%非常重要的强制对称
K_red = (K_red + K_red')/2;
M_red = (M_red + M_red')/2;

[phai_red_CB, W2_r_CB] = eigs(K_red, M_red, maxModesSolve, 'SM');   % nred x nev_compare
freq_red = sqrt(abs(diag(W2_r_CB))) / (2*pi);

aa3 = phai_red_CB'*M_red*phai_red_CB;

phi_red2full = T * phai_red_CB;

for i = 7:maxModesSolve
    MAC_full2(i) = (phai_full(:,i)'*phi_red2full(:,i))^2/...
        ((phai_full(:,i)'*phai_full(:,i))*(phi_red2full(:,i)'*phi_red2full(:,i)));
end

phai_b_CB = phai_red_CB(maxModesSolve+1:end,1:maxModesSolve);

phai_i = phai(interior,:);
phai_b = phai(retained,:);
for i = 7:maxModesSolve
    MAC_phaib(i) = (phai_b(:,i)'*phai_b_CB(:,i))^2/...
        ((phai_b(:,i)'*phai_b(:,i))*(phai_b_CB(:,i)'*phai_b_CB(:,i)));
end


%对于重根的振型子空间问题
bb1 = sum(phai_full(:,12).^2+phai_full(:,13).^2);
bb2 = sum(phai_sort(:,12).^2+phai_sort(:,13).^2);

bb3 = sum(phai_b(:,12).^2+phai_b(:,13).^2);
bb4 = sum(phai_b_CB(:,12).^2+phai_b_CB(:,13).^2);

bb5 = sum(phai_b(:,1).^2+phai_b(:,2).^2 + phai_b(:,3).^2+phai_b(:,4).^2+phai_b(:,5).^2+phai_b(:,6).^2);
bb6 = sum(phai_b_CB(:,1).^2+phai_b_CB(:,2).^2 + phai_b_CB(:,3).^2+phai_b_CB(:,4).^2+phai_b_CB(:,5).^2+phai_b_CB(:,6).^2);

%Z-Rx-Ry刚体模态分离,自己构造的刚体向量,为了统一弧度与米的量纲，除以半径
%但是不能简单的只除以半径，仅仅用"几何半径"RR（或特征长度Lc）去归一化，
% 确实会导致Z轴模态和Rx,Rx,Ry模态在能量上差一倍（针对圆形均质薄板）。

%正六边形，Rx，Ry的转动惯量为5/24*m *R^2
%若想要让z方向与Rx，Ry方向有相同量纲，需要将其统一，即phai_Ry，phai_Rx除以R*sqrt(5/24)
phai_Ry = -bottomnode_coodinfo_sorted(1,:)/(0.2*sqrt(5/24));
phai_Rx = bottomnode_coodinfo_sorted(2,:)/(0.2*sqrt(5/24));
%z部分使用单位向量
phai_Z = ones(size(bottomnode_coodinfo_sorted(3,:)));

phai_rigid_self = [phai_Z',phai_Rx',phai_Ry'];

rank_exam1 = [phai_rigid_self,phai_b(:,1:6)];
rank1 = rank(rank_exam1,1e-3);

rank_exam2 = [phai_rigid_self,Tn_bottomface_sorted(:,1:6)];
rank2 = rank(rank_exam2,1e-3);

[U1,sigma_rank_exam,V1] = svd(phai_rigid_self);



% phai_rigid_self = U1(:,1:3);

AA = phai_rigid_self'*phai_rigid_self;

