% close all;
% clc;
% clear;


ts = 0.0002;  
fs = 1/ts;

% load('Rz.mat');
% load('X.mat');
% load('Y.mat');
% load('ux.mat');
% load('uy.mat');
% load('urz.mat');
% load('ux1.mat');
% load('uy1.mat');
% load('urz1.mat');


% [f, AX]  = Function_SFZ_FFT((X),fs);
% [f, AY] = Function_SFZ_FFT((Y),fs);
% [f, ARz]  = Function_SFZ_FFT((Rz),fs);
% [f, AUx] = Function_SFZ_FFT((ux),fs);
% [f, AUy] = Function_SFZ_FFT((uy),fs);
% [f, AUrz] = Function_SFZ_FFT((urz),fs);
[f, AU1] = Function_SFZ_FFT((u1),fs);
[f, AU2] = Function_SFZ_FFT((u2),fs);
[f, AU3] = Function_SFZ_FFT((u3),fs);
[f, AU4] = Function_SFZ_FFT((u4),fs);
[f, AY1] = Function_SFZ_FFT((y1),fs);
[f, AY2] = Function_SFZ_FFT((y2),fs);
[f, AY3] = Function_SFZ_FFT((y3),fs);
[f, AY4] = Function_SFZ_FFT((y4),fs);

[f, AURX] = Function_SFZ_FFT((ux),fs);
[f, AURY] = Function_SFZ_FFT((uy),fs);
[f, AUZ] = Function_SFZ_FFT((urz),fs);

fs=5000;
NFFT=(length(ux));
ny=fft(X,NFFT); 
nu1=fft(ux,NFFT);
nu2=fft(uy,NFFT);
nu3=fft(urz,NFFT);
f=fs/2*linspace(0,1,NFFT/2+1); 
%B; %数据初始拟合位置  
abs_y=(ny(1:NFFT/2+1)); abs_x1=(nu1(1:NFFT/2+1)); abs_x2=(nu2(1:NFFT/2+1)); abs_x3=(nu3(1:NFFT/2+1));

AYY1  = mag2db(abs(AY1));
AYY2  = mag2db(abs(AY2));
AYY3  = mag2db(abs(AY3));
AYY4  = mag2db(abs(AY4));

AUU1  = mag2db(abs(AU1));
AUU2  = mag2db(abs(AU2));
AUU3  = mag2db(abs(AU3));
AUU4  = mag2db(abs(AU4));


U = [AU1.';
     AU2.';
     AU3.';
     AU4.'];

U11 = [AURX.';
     AURY.';
     AUZ.'];

[tf_magitude_tr1_LRM,Tr1,theta1_store,P1,Pdb_1_LRM,Pphase_1_LRM] = LRM_MISO(AY1,U,9,1);
[tf_magitude_tr2_LRM,Tr2,theta2_store,P2,Pdb_2_LRM,Pphase_2_LRM] = LRM_MISO(AY2,U,9,1);
[tf_magitude_tr3_LRM,Tr3,theta3_store,P3,Pdb_3_LRM,Pphase_3_LRM] = LRM_MISO(AY3,U,9,1);
[tf_magitude_tr4_LRM,Tr4,theta4_store,P4,Pdb_4_LRM,Pphase_4_LRM] = LRM_MISO(AY4,U,9,1);

% P1 = conj(P1);
% P2 = conj(P2);
% P3 = conj(P3);
% P4 = conj(P4);

T_U11 = Tr1(1,:)./AU1';
T_U12 = Tr1(2,:)./AU2';
T_U13 = Tr1(3,:)./AU3';
T_U14 = Tr1(4,:)./AU4';


T_U21 = Tr2(1,:)./AU1';
T_U22 = Tr2(2,:)./AU2';
T_U23 = Tr2(3,:)./AU3';
T_U24 = Tr2(4,:)./AU4';

T_U31 = Tr3(1,:)./AU1';
T_U32 = Tr3(2,:)./AU2';
T_U33 = Tr3(3,:)./AU3';
T_U34 = Tr3(4,:)./AU4';

T_U41 = Tr4(1,:)./AU1';
T_U42 = Tr4(2,:)./AU2';
T_U43 = Tr4(3,:)./AU3';
T_U44 = Tr4(4,:)./AU4';

magTU11 = mag2db(abs(T_U11));
magTU12 = mag2db(abs(T_U12));
magTU13 = mag2db(abs(T_U13));
magTU14 = mag2db(abs(T_U14));

magTU21 = mag2db(abs(T_U21));
magTU22 = mag2db(abs(T_U22));
magTU23 = mag2db(abs(T_U23));
magTU24 = mag2db(abs(T_U24));

magTU31 = mag2db(abs(T_U31));
magTU32 = mag2db(abs(T_U32));
magTU33 = mag2db(abs(T_U33));
magTU34 = mag2db(abs(T_U34));

magTU41 = mag2db(abs(T_U41));
magTU42 = mag2db(abs(T_U42));
magTU43 = mag2db(abs(T_U43));
magTU44 = mag2db(abs(T_U44));


[magreal11,phareal11,wout] = bode(P44(1,1),2*pi*f);
[magreal12,phareal12,wout] = bode(P44(1,2),2*pi*f);
[magreal13,phareal13,wout] = bode(P44(1,3),2*pi*f);
[magreal14,phareal14,wout] = bode(P44(1,4),2*pi*f);

[magreal21,phareal21,wout] = bode(P44(2,1),2*pi*f);
[magreal22,phareal22,wout] = bode(P44(2,2),2*pi*f);
[magreal23,phareal23,wout] = bode(P44(2,3),2*pi*f);
[magreal24,phareal24,wout] = bode(P44(2,4),2*pi*f);

[magreal31,phareal31,wout] = bode(P44(3,1),2*pi*f);
[magreal32,phareal32,wout] = bode(P44(3,2),2*pi*f);
[magreal33,phareal33,wout] = bode(P44(3,3),2*pi*f);
[magreal34,phareal34,wout] = bode(P44(3,4),2*pi*f);

[magreal41,phareal41,wout] = bode(P44(4,1),2*pi*f);
[magreal42,phareal42,wout] = bode(P44(4,2),2*pi*f);
[magreal43,phareal43,wout] = bode(P44(4,3),2*pi*f);
[magreal44,phareal44,wout] = bode(P44(4,4),2*pi*f);

frfreal11 = squeeze(magreal11.*exp(1i*deg2rad(phareal11)));
frfreal12 = squeeze(magreal12.*exp(1i*deg2rad(phareal12)));
frfreal13 = squeeze(magreal13.*exp(1i*deg2rad(phareal13)));
frfreal14 = squeeze(magreal14.*exp(1i*deg2rad(phareal14)));

frfreal21 = squeeze(magreal21.*exp(1i*deg2rad(phareal21)));
frfreal22 = squeeze(magreal22.*exp(1i*deg2rad(phareal22)));
frfreal23 = squeeze(magreal23.*exp(1i*deg2rad(phareal23)));
frfreal24 = squeeze(magreal24.*exp(1i*deg2rad(phareal24)));

frfreal31 = squeeze(magreal31.*exp(1i*deg2rad(phareal31)));
frfreal32 = squeeze(magreal32.*exp(1i*deg2rad(phareal32)));
frfreal33 = squeeze(magreal33.*exp(1i*deg2rad(phareal33)));
frfreal34 = squeeze(magreal34.*exp(1i*deg2rad(phareal34)));

frfreal41 = squeeze(magreal41.*exp(1i*deg2rad(phareal41)));
frfreal42 = squeeze(magreal42.*exp(1i*deg2rad(phareal42)));
frfreal43 = squeeze(magreal43.*exp(1i*deg2rad(phareal43)));
frfreal44 = squeeze(magreal44.*exp(1i*deg2rad(phareal44)));

err11 = (P1(1,:) - frfreal11.').';
err12 = (P1(2,:) - frfreal12.')';
err13 = (P1(3,:) - frfreal13.')';
err14 = (P1(4,:) - frfreal14.')';

err21 = (P2(1,:) - frfreal21.').';
err22 = (P2(2,:) - frfreal22.').';
err23 = (P2(3,:) - frfreal23.').';
err24 = (P2(4,:) - frfreal24.').';

err31 = (P3(1,:) - frfreal31.').';
err32 = (P3(2,:) - frfreal32.').';
err33 = (P3(3,:) - frfreal33.').';
err34 = (P3(4,:) - frfreal34.').';

err41 = (P4(1,:) - frfreal41.').';
err42 = (P4(2,:) - frfreal42.').';
err43 = (P4(3,:) - frfreal43.').';
err44 = (P4(4,:) - frfreal44.').';

magerr11 = mag2db(abs(err11));
magerr12 = mag2db(abs(err12));
magerr13 = mag2db(abs(err13));
magerr14 = mag2db(abs(err14));

magerr21 = mag2db(abs(err21));
magerr22 = mag2db(abs(err22));
magerr23 = mag2db(abs(err23));
magerr24 = mag2db(abs(err24));

magerr31 = mag2db(abs(err31));
magerr32 = mag2db(abs(err32));
magerr33 = mag2db(abs(err33));
magerr34 = mag2db(abs(err34));

magerr41 = mag2db(abs(err41));
magerr42 = mag2db(abs(err42));
magerr43 = mag2db(abs(err43));
magerr44 = mag2db(abs(err44));


figure;
subplot(4,4,1);
semilogx(f,Pdb_1_LRM(1,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr11,'r');
hold on;
bode(P44(1,1),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');


subplot(4,4,2);
semilogx(f,Pdb_1_LRM(2,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr12,'r');
hold on;
bode(P44(1,2),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,3);
semilogx(f,Pdb_1_LRM(3,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr13,'r');
hold on;
bode(P44(1,3),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,4);
semilogx(f,Pdb_1_LRM(4,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr14,'r');
hold on;
bode(P44(1,4),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,5);
semilogx(f,Pdb_2_LRM(1,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr21,'r');
hold on;
bode(P44(2,1),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,6);
semilogx(f,Pdb_2_LRM(2,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr22,'r');
hold on;
bode(P44(2,2),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,7);
semilogx(f,Pdb_2_LRM(3,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr23,'r');
hold on;
bode(P44(2,3),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,8);
semilogx(f,Pdb_2_LRM(4,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr24,'r');
hold on;
bode(P44(2,4),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,9);
semilogx(f,Pdb_3_LRM(1,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr31,'r');
hold on;
bode(P44(3,1),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,10);
semilogx(f,Pdb_3_LRM(2,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr32,'r');
hold on;
bode(P44(3,2),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,11);
semilogx(f,Pdb_3_LRM(3,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr33,'r');
hold on;
bode(P44(3,3),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,12);
semilogx(f,Pdb_3_LRM(4,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr34,'r');
hold on;
bode(P44(3,4),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,13);
semilogx(f,Pdb_4_LRM(1,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr41,'r');
hold on;
bode(P44(4,1),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,14);
semilogx(f,Pdb_4_LRM(2,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr42,'r');
hold on;
bode(P44(4,2),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,15);
semilogx(f,Pdb_4_LRM(3,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr43,'r');
hold on;
bode(P44(4,3),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

subplot(4,4,16);
semilogx(f,Pdb_4_LRM(4,:)','LineWidth',2,'Color','g');
hold on;
semilogx(f,magerr44,'r');
hold on;
bode(P44(4,4),bdoption);
xlabel('频率(Hz)');
ylabel('幅值(dB)');

