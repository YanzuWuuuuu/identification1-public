function [A,B,G_piao] = obtain_Gpiao1(G,ifrf,s,belta,noni,K)
%OBTAIN_GPIAO  获得 G_piao(i) = G_frf - sum(B/A)(除了第i个的其他估计量的总和)
% 输入K:子系统个数 G：扫频frf  ifrf:频点序号 s：当前频点jw值 belta:估计参数
% 输出 G_piao

A(1) = 1+belta(1,1)*s+belta(2,1)*s^2;
B{1} = (belta(3:noni+2)+ belta(3+noni:2*noni+2)*s+belta(3+2*noni:3*noni+2)*s^2)';
for i = 2:K
    A(i) = 1+belta((3*noni+2)+(i-2)*(noni+2)+1,i)*s+belta((3*noni+2)+(i-2)*(noni+2)+2,i)*s^2;
    B{i} = belta((3*noni+2)+(i-2)*(noni+2)+3:(3*noni+2)+(i-1)*(noni+2),i);
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