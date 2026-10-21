flex_target_mode_num = 2;
%柔性溢出模态数量
flex_spillover_mode_num = 3;
%组装的全振型矩阵，其具有3个刚体模态与目标加溢出模态总和的柔性模态个数
phai_rigid_bottom_z = phai_rigid_self;
phai_construct = [phai_rigid_bottom_z,Tn_bottomface_sorted(:,7:6+flex_target_mode_num+flex_spillover_mode_num)];


act_indice = [act_idx1,act_idx2,act_idx3,act_idx4];
sen_indice = [sen_idx1,sen_idx2,sen_idx3,sen_idx4];
act_location = bottomnode_coodinfo_sorted(:,act_indice);
sen_location = bottomnode_coodinfo_sorted(:,sen_indice);

plot(bottomnode_coodinfo_sorted(1,:), bottomnode_coodinfo_sorted(2,:), 'o', 'MarkerSize', 4);
hold on; axis equal; grid on;
plot(act_location(1,:), act_location(2,:), 'x', 'MarkerSize', 10, 'LineWidth', 3);
hold on;
plot(sen_location(1,:), sen_location(2,:), 'o', 'MarkerSize', 10, 'LineWidth', 3);
title('Bottom Face Nodes (sorted by x then y)');
xlabel('X'); ylabel('Y');

phai_act = phai_construct(act_indice,:);
%转化为振宇论文中的排列
Q_act = phai_act'*phai_act;
act_cond_rigid = cond(Q_act(1:3,1:3));
Br = phai_act(:,1:3)';
Bf = phai_act(:,4:end)';
Pbnull = null(Br);
B_equal = Bf*Pbnull;
Q_BT_Pnull = B_equal*B_equal';

phai_sen = phai_construct(sen_indice,:);
%转化为振宇论文中的排列
Q_sen = phai_sen'*phai_sen;
sen_cond_rigid = cond(Q_sen(1:3,1:3));
Cr = phai_sen(:,1:3);
Cf = phai_sen(:,4:end);
Pcnull = null(Cr')';
C_equal = Pcnull*Cf;
Q_Pcnull_C = C_equal'*C_equal;


% %逻辑轴输入输出端口的传递函数
% Tu = pinv(Br);
% Ty = pinv(Cr')';
% 
% %冗余输入输出空间端口的传递函数

