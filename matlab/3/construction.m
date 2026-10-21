load("vecM_est.mat");
load_data;
f_select = 20:2:1000;
N = length(f_select);
[frf_select,w_select] = select_frf(vecG,f_select',f');
est_frf1 = zeros(no*ni,N);
est_frf2 = zeros(no*ni,N);
for ifrf = 1:N
    s = 1i*w_select(ifrf);
    B_frf1 = vecM_est(1:no*ni,1);
    est_frf1(:,ifrf) =  B_frf1/s^2;
end
for ifrf = 1:N
    s = 1i*w_select(ifrf);
    est_frf2(:,ifrf) = zeros(no*ni,1);
    for j = 1:isubsystem
    A = 1 + belta_est((no*ni+2)*(j-1)+1,j)*s + belta_est((no*ni+2)*(j-1)+2,j)*s^2;
    B = belta_est((no*ni+2)*(j-1)+3:(no*ni+2)*j,j);
    est_frf2(:,ifrf) = est_frf2(:,ifrf) + B/A;
    end
end

est_total_frf = est_frf1 + est_frf2;
err_frf = frf_select - est_total_frf;

mag_est_frf_body = 20*log10(abs(est_frf1));
mag_est_frf = 20*log10(abs(est_total_frf));
mag_frf = 20*log10(abs(frf_select));
mag_err_frf = 20*log10(abs(err_frf));

figure;
for i = 1:ni
    for j = 1:no
    subplot(4,4,(j-1)*ni+i);
    semilogx(f_select,mag_frf((i-1)*no+j,:),'g','LineWidth',2);
    hold on;
    semilogx(f_select,mag_err_frf((i-1)*no+j,:),'r','LineWidth',2);
    hold on;
    semilogx(f_select,mag_est_frf((i-1)*no+j,:),'b','LineWidth',2);
    hold on;
    end
end

belta_real_total_temp = cat_diag(belta_real,belta_real0);
belta_est_total_temp = cat_diag(belta_divide_est,belta_est);

belta_real_total = zeros(length(belta_real_total_temp(:,1)),1);
belta_est_total = zeros(length(belta_est_total_temp(:,1)),1);

belta_real_total(1:no*ni,1) = belta_real_total_temp(1:no*ni,1);
belta_est_total(1:no*ni,1) = belta_est_total_temp(1:no*ni,1);
for i = 1:isubsystem+1
    belta_real_total(no*ni+(i-1)*(no*ni+2)+1:no*ni+i*(no*ni+2),1) = belta_real_total_temp(no*ni+(i-1)*(no*ni+2)+1:no*ni+i*(no*ni+2),i+1);
    belta_est_total(no*ni+(i-1)*(no*ni+2)+1:no*ni+i*(no*ni+2),1) = belta_est_total_temp(no*ni+(i-1)*(no*ni+2)+1:no*ni+i*(no*ni+2),i+1);
end

norm_err = norm(belta_real_total-belta_est_total)/norm(belta_real_total)

isubsystem = isubsystem+1;
% save('additive_model_para','belta_est_total','belta_real_total','isubsystem');