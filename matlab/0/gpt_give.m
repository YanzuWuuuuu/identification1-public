%% 演示 G1,G2,G3 的有理逼近与可辨识性
clear; clc; close all;

% 频率点 (模拟频域数据)
w = linspace(0.1, 5, 200);   % rad/s
s = 1j*w;

%% 1. 定义真实系统 (频域传递函数形式)
G1_true = (1 + 0.5*s) ./ (1 + 0.2*s);
G2_true = (0.5 + 0.2*s) ./ (1 - 0.1*s);
G3_true = (0.2 + 0.1*s) ./ (1 + 0.05*s);

%% 输入设计
% ---- 方案 A: 独立输入 ----
U1_A = ones(size(s));
U2_A = sin(w);
U3_A = cos(w);

% ---- 方案 B: 相关输入 (退化) ----
U1_B = ones(size(s));
U2_B = sin(w);
U3_B = U1_B + U2_B;   % 与前两者相关

% ---- 方案 C: 优化输入 (正交/独立) ----
U1_C = ones(size(s));
U2_C = exp(-0.2*w);
U3_C = w;  % 线性函数

% 输出信号
Y_A = G1_true.*U1_A + G2_true.*U2_A + G3_true.*U3_A;
Y_B = G1_true.*U1_B + G2_true.*U2_B + G3_true.*U3_B;
Y_C = G1_true.*U1_C + G2_true.*U2_C + G3_true.*U3_C;

%% 2. 点对点恢复 G1,G2,G3
% 定义一个函数：点对点解最小二乘
recover = @(U1,U2,U3,Y) arrayfun(@(k) ...
    (pinv([U1(k),U2(k),U3(k)])*Y(k)), 1:length(Y), 'UniformOutput', false);

% A 案例
G_A = cell2mat(recover(U1_A,U2_A,U3_A,Y_A));
G1_A = G_A(1,:); G2_A = G_A(2,:); G3_A = G_A(3,:);

% B 案例
G_B = cell2mat(recover(U1_B,U2_B,U3_B,Y_B));
G1_B = G_B(1,:); G2_B = G_B(2,:); G3_B = G_B(3,:);

% C 案例
G_C = cell2mat(recover(U1_C,U2_C,U3_C,Y_C));
G1_C = G_C(1,:); G2_C = G_C(2,:); G3_C = G_C(3,:);

%% 3. 有理函数拟合 ([1/1] 和 [2/2])
% 使用 MATLAB 内置函数 invfreqs 进行拟合
n_order = [1 2]; % [1/1] 与 [2/2]
w_fit = w;       % 拟合频率点

fit_rational = @(Gdata,order) ...
    invfreqs(Gdata, s, order, order); % 返回传递函数分子/分母系数

[b1,a1] = fit_rational(G1_A,1);
[b2,a2] = fit_rational(G1_A,2);

G1_fit_11 = freqs(b1,a1,w);
G1_fit_22 = freqs(b2,a2,w);

%% 4. 输入矩阵条件数比较
Umat_A = [U1_A(:),U2_A(:),U3_A(:)];
Umat_B = [U1_B(:),U2_B(:),U3_B(:)];
Umat_C = [U1_C(:),U2_C(:),U3_C(:)];

cond_A = cond(Umat_A);
cond_B = cond(Umat_B);
cond_C = cond(Umat_C);

fprintf('条件数对比:\n');
fprintf('独立输入 A: %.2e\n', cond_A);
fprintf('相关输入 B: %.2e\n', cond_B);
fprintf('优化输入 C: %.2e\n', cond_C);

%% 5. 画图比较
figure;
subplot(3,2,1);
plot(w,real(G1_true),'-',w,real(G1_A),'--');
title('G1 真值 vs A 恢复 (实部)'); legend('True','Recovered');

subplot(3,2,2);
plot(w,real(G1_true),'-',w,real(G1_B),'--');
title('G1 真值 vs B 恢复 (实部)'); legend('True','Recovered');

subplot(3,2,3);
plot(w,real(G1_true),'-',w,real(G1_C),'--');
title('G1 真值 vs C 恢复 (实部)'); legend('True','Recovered');

subplot(3,2,4);
plot(w,abs(G1_true),'-',w,abs(G1_fit_11),'--',w,abs(G1_fit_22),':');
title('G1 拟合 (A 案例)'); legend('True','[1/1]','[2/2]');

subplot(3,2,5);
semilogy(w,abs(G1_true-G1_fit_11),'-',w,abs(G1_true-G1_fit_22),'--');
title('G1 拟合误差'); legend('[1/1]','[2/2]');

subplot(3,2,6);
bar([cond_A cond_B cond_C]);
set(gca,'XTickLabel',{'独立 A','相关 B','优化 C'});
title('输入矩阵条件数对比');
