function [frf_select,w_select] = select_frf(frf,f_select,f,count)
%SELECT_FRF 在非参数频响中获得所选范围频响
% 输入参数：  frf：全频响 f_select:所选频段 f：全频率
% 输出参数：  frf_select:所选频段频响 w_select:所选角频率

% 计算欧几里得距离矩阵（高效向量化方法）
distSq = sum(f_select.^2, 2) + sum(f.^2, 2)' - 2 * (f_select * f');

% 找到每个data2点的最近邻索引
[~, indices] = min(distSq, [], 2);
fidex = f(indices);
w_select = 2*pi*fidex;
    for i = 1:count
    frf_select{i} = frf{i}(:,indices);
    end
end

