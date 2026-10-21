clc;
clear;
close all;
%初值估计模块


%调试所选频响范围的df即可 5:1:1000;中间那个数为待调试数
no = 4;
ni = 4;
% 选择所选频响范围,待修改变量
f_select = 40:20:300;
f_pri = [285.25];
isubsystem = length(f_pri);
% frf导入
model;
load_data;
load('belta_init.mat');
f = f';
[frf_select,w_select] = select_frf(vecG,f_select',f);
mean_frf = frf_select;
N_select = length(w_select);
%获得协方差逆的权矩阵
% Weight = weight_matrix_cal(N_select,no,ni,count,frf_select);
%系统频域参数nior估计
%构建schur权矩阵
v = [];
W = [];
for ifrf = 1:N_select
    for j = 1:no*ni
        v_temp(j) = 1/(abs(mean_frf(j,ifrf))+1e-10);
    end
    w = diag(v_temp);
    w = cat_diag(w,w);
    W = cat_diag(W,w);
end
Weight = W;



belta_est_real = zeros((3*no*ni+2),1);
belta_est_real(1,1) = belta_real(no*ni+1,2);
belta_est_real(2,1) = belta_real(no*ni+2,2);
B0_real = belta_real(1:no*ni)';
B1_real = belta_est_real(1,1)*belta_real(1:no*ni)';
B2_real = belta_est_real(2,1)*belta_real(1:no*ni)'+belta_real(no*ni+3:2*no*ni+2)';
belta_est_real(3:no*ni+2,1) = B0_real;
belta_est_real(no*ni+2+1:2*no*ni+2,1) = B1_real;
belta_est_real(2*no*ni+2+1:3*no*ni+2,1) = B2_real;

belta0 = zeros((no*ni)*(isubsystem+1),1);
belta0(1:no*ni,1) = ones(no*ni,1)';
den_para0 = [];
for i = 1:length(f_pri)
    B0{i} = ones(1,16);
    theta0{i} = [1 0.02, B0{i}]'/((2*pi*f_pri(i))^2);
    den_para0{i} =  [1 2*0.02*2*pi*f_pri(i)]/((2*pi*f_pri(i))^2);
end
% 
for i = 2:isubsystem+1
    belta0((i-1)*(no*ni)+1:i*(no*ni),1) = B0{i-1}'/((2*pi*f_pri(i-1))^2);
end

e = 1e-10;
M = 100;
[num_init,belta_init] = init_estimate(f_select,mean_frf,den_para0,Weight,belta0,no*ni,e,M,isubsystem,0);
belta_est0 = zeros((no*ni+2)*(isubsystem-1)+(3*no*ni+2),isubsystem);
belta_est0(1,1) = den_para0{1}(2);
belta_est0(2,1) = den_para0{1}(1);
B0_temp = belta_init(1:no*ni);
B1_temp = den_para0{1}(2)*belta_init(1:no*ni);
B2_temp = den_para0{1}(1)*belta_init(1:no*ni)+belta_init(no*ni+1:2*no*ni);

belta_est0(3:no*ni+2,1) = B0_temp;
belta_est0(no*ni+2+1:2*no*ni+2,1) = B1_temp;
belta_est0(2*no*ni+2+1:3*no*ni+2,1) = B2_temp;
for i = 2:isubsystem
    belta_est0(3*no*ni+2+(i-2)*(no*ni+2)+1,i) = den_para0{i}(2);
    belta_est0(3*no*ni+2+(i-2)*(no*ni+2)+2,i) = den_para0{i}(1);
    for j = 1:ni
        for k = 1:no
    belta_est0((3*no*ni+2)+(i-2)*(no*ni+2)+2+(j-1)*no+k,i) = belta_init((i)*(no*ni)+(j-1)*no+k);
        end
    end
end

[num_est,belta_est] = SRIVC1(f_select,mean_frf,Weight,belta_est0,no*ni,e,M,isubsystem);


M_est_init = belta_init(1:no*ni);
for i = 1:isubsystem
    B_est{i} = belta_init(no*ni*i+1:no*ni*(i+1),1);
end

for i = 1:N_select
    s = 1i*w_select(i);
    frf_est_init(:,i) = M_est_init/(s^2);
    for j = 1:isubsystem
        A = 1 + den_para0{j}(2)*s + den_para0{j}(1)*s^2;
        frf_est_init(:,i) = frf_est_init(:,i)+ B_est{j}/A;
    end
end

B0_est = belta_est(3:no*ni+2,1);
B1_est = belta_est(no*ni+2+1:2*no*ni+2,1);
B2_est = belta_est(2*no*ni+2+1:3*no*ni+2,1);

for i = 1:N_select
    s = 1i*w_select(i);
    A = s^2*(1+belta_est(1,1)*s+belta_est(2,1)*s^2);
    frf_est(:,i) = (B0_est + B1_est*s + B2_est*s^2)/A;
    for j = 2:isubsystem
         A = 1 + belta_est0(3*no*ni+2+(j-2)*(no*ni+2)+1,j)*s + belta_est0(3*no*ni+2+(j-2)*(no*ni+2)+2,j)*s^2;
         B = belta_est0(3*no*ni+2+(j-2)*(no*ni+2)+3:3*no*ni+2+(j-1)*(no*ni+2))';
         frf_est(:,i) = frf_est(:,i) + B/A;
    end
end

%将刚体模态与第一阶柔性模态分开
Brb_est = B0_est;
Bfb1_est = B2_est-belta_est(2,1)*B0_est;
belta_divide_est = zeros(size(belta_real));
belta_divide_est(1:no*ni,1) = Brb_est;
belta_divide_est(no*ni+1,2) = belta_est(1,1);
belta_divide_est(no*ni+2,2) = belta_est(2,1);
belta_divide_est(no*ni+3:2*no*ni+2,2) = Bfb1_est;
%画图

err_frf = abs(mean_frf)-abs(frf_est);
mag_err_frf = 20*log10(abs(err_frf));

mag_frf_est_init = 20*log10(abs(frf_est_init));
mag_mean_frf = 20*log10(abs(mean_frf));
mag_frf_est = 20*log10(abs(frf_est));


figure(1);
for i = 1:ni
    for j = 1:no
    subplot(4,4,(j-1)*ni+i);
%     semilogx(f_select,mag_frf_est_init((i-1)*no+j,:),'k','LineWidth',2);
%     hold on;
    semilogx(f_select,mag_mean_frf((i-1)*no+j,:),'g','LineWidth',2);
    hold on;
    semilogx(f_select,mag_frf_est((i-1)*no+j,:),'b','LineWidth',2);
    hold on;
    semilogx(f_select,mag_err_frf((i-1)*no+j,:),'r','LineWidth',2);
    end
end
norm_err = norm(belta_divide_est-belta_real)/norm(belta_real)

N_f = length(f);
w = 2*pi*f;
for ifrf = 1:N_f
    s = 1i*w(ifrf);
    A = 1+belta_est(1,1)*s+belta_est(2,1)*s^2;
    frf_minus(:,ifrf) = Brb_est/s^2 + Bfb1_est/A;
end

% for ifrf = 1:N_f
%     s = 1i*w(ifrf);
%     A = 1+belta_real(no*ni+1,2)*s+belta_real(no*ni+2,2)*s^2;
%     frf_minus(:,ifrf) = belta_real(1:no*ni,1)/s^2 + belta_real(no*ni+3:2*no*ni+2,2)/A;
% end


errand_frf = LRMfrf - frf_minus;

% mag_errand_frf = 20*log10(abs(errand_frf));
% figure(2);
% for i = 1:ni
%     for j = 1:no
%     subplot(4,4,(j-1)*ni+i);
%     semilogx(f,mag_errand_frf((i-1)*no+j,:),'LineWidth',2);
%     end
% end

% save('second_step_data',"errand_frf",'f','belta_divide_est',...
%     'belta_real','count');