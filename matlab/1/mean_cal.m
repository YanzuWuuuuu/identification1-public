function [mean_frf] = mean_cal(vecG,count)
%MEAN_CAL 平均值计算
%   输入：vecG三维frf矩阵，第一维为图个数，第二维为frf，第三维为第几次frf, count为输入次数
temp = 0;
for i = 1:count
    temp =temp + vecG{i};
end
mean_frf = temp/count;
end

