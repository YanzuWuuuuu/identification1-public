function [A,B,G_piao] = obtain_Gpiao(G,ifrf,s,belta,noni,K)
%OBTAIN_GPIAO  获得 G_piao(i) = G_frf - sum(B/A)(除了第i个的其他估计量的总和)
% 输入K:子系统个数 G：扫频frf  ifrf:频点序号 s：当前频点jw值 belta:估计参数
% 输出 G_piao

frf = G;
for i = 1:K
    A(i) = 1+belta((i-1)*(noni+2)+1,i)*s+belta((i-1)*(noni+2)+2,i)*s^2;
    B{i} = belta((i-1)*(noni+2)+3:(i)*(noni+2),i);
end
for i = 1:K
    temp = G(:,ifrf);
    for j = 1:K
        if j~=i
          temp = temp-B{j}/A(j);
        end
    end
        G_piao{i} = temp;
end
end