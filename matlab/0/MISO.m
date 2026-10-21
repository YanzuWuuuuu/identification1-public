clear; 
clc; 
close all;
fs=5000;
%%
%输入数据
% ux=load('ux.mat');uy=load('uy.mat');urz=load('urz.mat'); 
% ux=ux.ux; uy=uy.uy; urz=urz.urz;
% % ux1=load('ux1.mat');uy1=load('uy1.mat');urz1=load('urz1.mat'); 
% % ux=ux1.ux1; uy=uy1.uy1; urz=urz1.urz1;
% %输出数据
% X=load('X.mat'); Y=load('Y.mat'); Rz=load('Rz.mat');
% X=X.X;Y=Y.Y;Rz=Rz.Rz;

%%
B=8;  %所得数据拟合的起始位置
% [ls2,f,tf_magitude_y11,tf_phase_y11,tf_magitude_y12,tf_phase_y12,tf_magitude_y13,tf_phase_y13,tf_magitude_y_ploy11,tf_phase_y_ploy11,tf_magitude_y_ploy12,tf_phase_y_ploy12,tf_magitude_y_ploy13,tf_phase_y_ploy13]=FunctionMOSO_FFT_LRM_X(Ux,Uy,Urz,x,B);
% %clear ls2
% [ls2,f,tf_magitude_y21,tf_phase_y21,tf_magitude_y22,tf_phase_y22,tf_magitude_y23,tf_phase_y23,tf_magitude_y_ploy21,tf_phase_y_ploy21,tf_magitude_y_ploy22,tf_phase_y_ploy22,tf_magitude_y_ploy23,tf_phase_y_ploy23]=FunctionMOSO_FFT_LRM_X(Ux,Uy,Urz,y,B);
% %clear ls2
% [ls2,f,tf_magitude_y31,tf_phase_y31,tf_magitude_y32,tf_phase_y32,tf_magitude_y33,tf_phase_y33,tf_magitude_y_ploy31,tf_phase_y_ploy31,tf_magitude_y_ploy32,tf_phase_y_ploy32,tf_magitude_y_ploy33,tf_phase_y_ploy33]=FunctionMOSO_FFT_LRM_X(Ux,Uy,Urz,rz,B);
% %clear ls2

[ls2,f,tf_magitude_y_ploy11,tf_phase_y_ploy11,tf_magitude_y_ploy12,tf_phase_y_ploy12,tf_magitude_y_ploy13,tf_phase_y_ploy13]=FunctionMOSO_FFT_LRM_X(ux,uy,urz,X,B);
%clear ls2
[ls2,f,tf_magitude_y_ploy21,tf_phase_y_ploy21,tf_magitude_y_ploy22,tf_phase_y_ploy22,tf_magitude_y_ploy23,tf_phase_y_ploy23]=FunctionMOSO_FFT_LRM_X(ux,uy,urz,Y,B);
%clear ls2
[ls2,f,tf_magitude_y_ploy31,tf_phase_y_ploy31,tf_magitude_y_ploy32,tf_phase_y_ploy32,tf_magitude_y_ploy33,tf_phase_y_ploy33]=FunctionMOSO_FFT_LRM_X(ux,uy,urz,Rz,B);
%clear ls2

%%
%获取真实控制对象的频响数据
NFFT=(length(X));
Mx=60;     %质量
My=146;    %质量
Jz=30;     %Rz转动惯量
[PX,ls4] = Function_TURE_sys(NFFT,Mx);
[PY,ls5] = Function_TURE_sys(NFFT,My);
[PRZ,ls6] = Function_TURE_sys(NFFT,Jz);

%%
%rel_BJ=abs(ls2(7:7168,2))-abs(ls3(7:7168,2));
%relWU=20*log10(rel_BJ);
%figure(4)
%semilogx(f(1:length(relWU)), relWU,'k');hold on;

figure(5)
subplot(3,3,1);
% semilogx(f,tf_magitude_y11,'b');hold on;
% semilogx(f, (ls4(:,4)),'k');hold on;
semilogx(f(1:length(tf_magitude_y_ploy11)), tf_magitude_y_ploy11,'r');hold on;
bode(P(1,1),bdoption);
subplot(3,3,2);
% semilogx(f,tf_magitude_y12,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy12)), tf_magitude_y_ploy12,'r');hold on;
bode(P(1,2),bdoption);
subplot(3,3,3);
% semilogx(f,tf_magitude_y13,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy13)), tf_magitude_y_ploy13,'r');hold on;
bode(P(1,3),bdoption);

subplot(3,3,4);
% semilogx(f,tf_magitude_y21,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy21)), tf_magitude_y_ploy21,'r');hold on;
bode(P(2,1),bdoption);
subplot(3,3,5);
% semilogx(f,tf_magitude_y22,'b');hold on;
% semilogx(f, (ls5(:,4)),'k');hold on;
semilogx(f(1:length(tf_magitude_y_ploy22)), tf_magitude_y_ploy22,'r');hold on;
bode(P(2,2),bdoption);
subplot(3,3,6);
% semilogx(f,tf_magitude_y23,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy23)), tf_magitude_y_ploy23,'r');hold on;
bode(P(2,3),bdoption);

subplot(3,3,7);
% semilogx(f,tf_magitude_y11,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy31)), tf_magitude_y_ploy31,'r');hold on;
bode(P(3,1),bdoption);
subplot(3,3,8);
% semilogx(f,tf_magitude_y12,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy32)), tf_magitude_y_ploy32,'r');hold on;
bode(P(3,2),bdoption);
subplot(3,3,9);
% semilogx(f,tf_magitude_y13,'b');hold on;
% semilogx(f, (ls6(:,4)),'k');hold on;
semilogx(f(1:length(tf_magitude_y_ploy33)), tf_magitude_y_ploy33,'r');hold on;
bode(P(3,3),bdoption);
% %%
figure(6)
subplot(3,3,1);
% semilogx(f,tf_phase_y11,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy11)), tf_phase_y_ploy11,'r');hold on;
subplot(3,3,2);
% semilogx(f,tf_phase_y12,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy12)), tf_phase_y_ploy12,'r');hold on;
subplot(3,3,3);
% semilogx(f,tf_phase_y13,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy13)), tf_phase_y_ploy13,'r');hold on;

subplot(3,3,4);
% semilogx(f,tf_phase_y21,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy21)), tf_phase_y_ploy21,'r');hold on;
subplot(3,3,5);
% semilogx(f,tf_phase_y22,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy22)), tf_phase_y_ploy22,'r');hold on;
subplot(3,3,6);
% semilogx(f,tf_phase_y23,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy23)), tf_phase_y_ploy23,'r');hold on;

subplot(3,3,7);
% semilogx(f,tf_phase_y11,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy31)), tf_phase_y_ploy31,'r');hold on;
subplot(3,3,8);
% semilogx(f,tf_phase_y12,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy32)), tf_phase_y_ploy32,'r');hold on;
subplot(3,3,9);
% semilogx(f,tf_phase_y13,'b');hold on;
semilogx(f(1:length(tf_magitude_y_ploy33)), tf_phase_y_ploy33,'r');hold on;
