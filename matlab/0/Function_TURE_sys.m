function [P,ls4] = Function_TURE_sys(NFFT,m)
%%
% NFFT:离散数据长度
% m:质量
% ls4：输出堆积

fs=5000;
%%
%对象----------
s=tf('s');
K1 = 0.013;  w1 = 2*pi*150; xi1 = 0.018;   %第一个谐振点  xi1 = 0.040;
K2 = 0.0023; w2 = 2*pi*292; xi2 = 0.018;   %第二个谐振点
K3 = 0.0038; w3 = 2*pi*510; xi3 = 0.027;   %第三个谐振点

P0 = 1/(m*s^2);
P1 = K1/(1*s^2+2*xi1*w1*s+w1^2);
P2 = K2/(1*s^2+2*xi2*w2*s+w2^2);
P3 = K3/(1*s^2+2*xi3*w3*s+w3^2); 
 P =  P0+ P1 + P2 + P3; 
% P =  P0 + P1;
 
sys=P;
    j2=0;
   % NFFT=(length(x));
 for w=fs/2*linspace(0,1,NFFT/2+1)*(2*pi)
     res = freqresp(sys, w);
     amp_true = abs(res);
     amp_db = 20*log10(abs(res));
     angle_true = angle(res)*180/pi;
     j2=j2+1;
     ls4(j2,:)=[j2,res,amp_true,amp_db,angle_true];
end