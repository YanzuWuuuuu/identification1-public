function [Weight_total] = weight_matrix_cal(N_select,no,ni,count,frf_select)
%
%   此处显示详细说明
Weight_total = zeros(2*N_select*no*ni,2*N_select*no*ni);
C_vecG = zeros(no*ni,no*ni,N_select);
for i = 1:count
vecG_select(:,:,i) = frf_select{i};
end
for i  = 1:N_select
    temp2 = squeeze(vecG_select(:,i,:));
    temp1 = temp2';
    for j = 1:count
        temp1(:,j) = temp1(:,j) - mean(temp1(:,j));
    end
    C_vecG = cov(temp1);
    %权矩阵
    Weight = inv(C_vecG);
    [temp3] = diag_real_imag_divide(Weight);
    Weight_total((i-1)*2*no*ni+1:i*2*no*ni,(i-1)*2*no*ni+1:i*2*no*ni) = temp3;
end
end

