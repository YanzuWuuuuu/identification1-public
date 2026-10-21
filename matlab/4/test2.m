% Levenberg-Marquardt 算法实现（指数函数拟合）
clear; clc;

% 步骤1: 生成带噪声的模拟数据
x = linspace(0, 2, 50)';
y_true = 2.5*exp(1.3*x) + 0.5;    % 真实参数 [a=2.5, b=1.3, c=0.5]
rng(1);
y_noise = y_true + 0.5*randn(size(x)); % 添加高斯噪声

% 步骤2: 算法参数设置
params = [1; 0.1; 0.1];   % 初始猜测参数 [a; b; c]
max_iter = 100;            % 最大迭代次数
tol = 1e-6;                % 收敛容差
lambda = 0.01;             % 初始阻尼因子
update_factor = 10;        % 阻尼因子调整系数

% 步骤3: 主循环
for iter = 1:max_iter
    % 计算当前残差和雅可比矩阵
    [r, J] = residual_jacobian(x, params, y_noise);
    
    % 检查收敛条件
    if norm(r) < tol
        break;
    end
    
    % 构建线性系统 (J'*J + lambda*diag(J'*J)) * delta = -J'*r
    H = J'*J;
    g = J'*r;
    H_lm = H + lambda*diag(diag(H));  % 添加阻尼项
    
    % 求解增量
    delta = -H_lm \ g;
    
    % 计算新参数和残差
    params_new = params + delta;
    r_new = residual_jacobian(x, params_new, y_noise);
    
    % 判断是否接受更新
    if norm(r_new) < norm(r)
        lambda = lambda / update_factor;
        params = params_new;
    else
        lambda = lambda * update_factor;
    end
end

% 步骤4: 显示结果
fprintf('算法在 %d 次迭代后收敛\n', iter);
fprintf('优化参数:\n a = %.4f\n b = %.4f\n c = %.4f\n', params);

% 步骤5: 绘制结果
y_fit = params(1)*exp(params(2)*x) + params(3);
figure;
plot(x, y_noise, 'bo', 'DisplayName', '带噪声数据');
hold on;
plot(x, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', '拟合曲线');
plot(x, y_true, 'k--', 'DisplayName', '真实曲线');
xlabel('x'); ylabel('y');
legend('Location', 'northwest');
title('Levenberg-Marquardt算法拟合结果');

% 残差和雅可比计算函数
function [r, J] = residual_jacobian(x, params, y)
    a = params(1);
    b = params(2);
    c = params(3);
    
    % 计算模型预测值
    y_pred = a*exp(b*x) + c;
    
    % 残差向量
    r = y_pred - y;
    
    % 雅可比矩阵
    J = zeros(length(x), 3);
    J(:,1) = exp(b*x);         % dr/da
    J(:,2) = a*x.*exp(b*x);    % dr/db
    J(:,3) = 1;                % dr/dc
end