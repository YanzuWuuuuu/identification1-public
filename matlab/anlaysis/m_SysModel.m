
%********************************************************
%Author:	Fazhi SONG, Harbin Institute of Technology, Harbin 150001, China
%Date:      2022-07-09
%********************************************************

%-----------------------------对象--------------------------------------
B1 = [0.0972,0.0990,-0.1405,-0.0024;
      0.1004,-0.1002,9.8523e-4,0.0011;
      0.9717,0.9893,0.9564,1.0781];
B2 = [-0.5405,0.5613,0.0015,-0.0292;
      0.0725,0.0316,0.6871,-0.1560;
      0.2364,0.2200,0.0601,-0.8658;
      -0.2834,-0.3269,-0.5538,-0.0619;
      -0.2543,0.1529,-0.0934,-0.0085];
B = [B1;B2];

C1 = [0.1714,0.0861,1.0014;
     0.1697,-0.0756,0.9986;
     -0.1708,-0.1108,1.0012;
     -0.1673,0.1034,0.9989];
C2 = [-0.5297,0.5215,0.6477,0.4121,-0.9896;
     0.5585,0.5490,0.5350,0.0562,1.0118;
     -0.6198,0.4402,0.6328,-0.0175,1.0857;
     0.6258,0.4527,0.6219,0.3133,-0.9832];
C = [C1,C2];
td = 2e-4;
Tf = B1'*inv((B1*B1'));
Tcb = inv((C1'*C1))*C1';
s = tf('s');
%%Jrx,Jry,mz
Jrx = 0.0439;
Jry = 0.0384;
mz = 4.2809;
phai_body = [1/(Jrx*s^2),1/(Jry*s^2),1/(mz*s^2)];
phai1 = tf(3,3);
for i = 1:length(phai_body)
    phai1(i,i) = phai_body(i);
end

%%internal mode
zeta = 0.001;
f = [286.43,344.58,488.27,705.67,738.83];
phai2 = tf(length(kesai),length(kesai));
for i = 1:length(kesai)
    phai2(i,i) = 1/(s^2+2*kesai(i)*(2*pi*f(i))*s+(2*pi*f(i))^2);
end

% phai2(1,1) = 0;

phai2(4,4) = 0;
phai2(5,5) = 0;

phai = blkdiag(phai1,phai2);
% D = exp(-td*s);
% phai = D.* phai;

% [num_P1 den_P1] = tfdata(phai(1,1),'v');
% [num_P2 den_P2] = tfdata(phai(2,2),'v');
% [num_P3 den_P3] = tfdata(phai(3,3),'v');
% [num_P4 den_P4] = tfdata(phai(4,4),'v');
% [num_P5 den_P5] = tfdata(phai(5,5),'v');
% [num_P6 den_P6] = tfdata(phai(6,6),'v');
% [num_P7 den_P7] = tfdata(phai(7,7),'v');
% [num_P8 den_P8] = tfdata(phai(8,8),'v');
%-----------------------------------------------

nrigidcontrol = 3;
nflexcontrol = 5;

f = [zeros(1,3),f];
wrad = 2*pi*f;
w2_ss = zeros(length(wrad),length(wrad));
zeta_w_2_ss = zeros(length(wrad),length(wrad));

for i = 1:length(wrad)
    w2_ss(i,i) = wrad(i)^2;
    zeta_w_2_ss(i,i) = 2*zeta*wrad(i);
end

Pbnull = null(B1);
Pcnull = null(C1')';

Tu = pinv(B1);
Ty = pinv(C1')';

[GG,GG11,GG12,GG21,GG22] = construct_rigidflex_decple_controlmodel(nrigidcontrol,nflexcontrol,w2_ss,zeta_w_2_ss,...
    B1,C1,B2,C2,Tu,Ty,Pcnull,Pbnull);

figure;
P = C*phai*B;
bode(P);

%逻辑轴输入输出端口的传递函数
figure;
bode(GG11);
xlim([10,1e4]);
%冗余输入输出空间端口的传递函数

figure;
bode(GG22);

Qb = B1*B1';
Qc = C1'*C1;

B1_piao = B1;
C1_piao = C1;

Rgx = sqrt(Jrx/mz);
Rgy = sqrt(Jry/mz);
B1_piao(1,:) = B1(1,:)/Rgx;
B1_piao(2,:) = B1(2,:)/Rgy;
B1_piao(3,:) = B1(3,:);

C1_piao(:,1) = C1(:,1)/Rgx;
C1_piao(:,2) = C1(:,2)/Rgy;
C1_piao(:,3) = C1(:,3);

Qb_piao = B1_piao*B1_piao';
Qc_piao = C1_piao'*C1_piao;

cond_Qb = cond(Qb_piao);
cond_Qc = cond(Qc_piao);