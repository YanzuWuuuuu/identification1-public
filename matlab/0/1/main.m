close all;
clc;
clear;


ts = 0.0002;  
fs = 1/ts;

load('u1.mat');
load('u2.mat');
load('u3.mat');
load('u4.mat');
load('y1.mat');
load('y2.mat');
load('y3.mat');
load('y4.mat');


[f, AU1] = Function_SFZ_FFT((u1),fs);
[f, AU2] = Function_SFZ_FFT((u2),fs);
[f, AU3] = Function_SFZ_FFT((u3),fs);
[f, AU4] = Function_SFZ_FFT((u4),fs);
[f, AY1] = Function_SFZ_FFT((y1),fs);
[f, AY2] = Function_SFZ_FFT((y2),fs);
[f, AY3] = Function_SFZ_FFT((y3),fs);
[f, AY4] = Function_SFZ_FFT((y4),fs);

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

[tf_magitude_tr1_LRM,Tr1,theta1_store,P1,Pdb_1_LRM,Pphase_1_LRM] = LRM_MISO(AY1,U,9,1);
[tf_magitude_tr2_LRM,Tr2,theta2_store,P2,Pdb_2_LRM,Pphase_2_LRM] = LRM_MISO(AY2,U,9,1);
[tf_magitude_tr3_LRM,Tr3,theta3_store,P3,Pdb_3_LRM,Pphase_3_LRM] = LRM_MISO(AY3,U,9,1);
[tf_magitude_tr4_LRM,Tr4,theta4_store,P4,Pdb_4_LRM,Pphase_4_LRM] = LRM_MISO(AY4,U,9,1);

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

figure;
subplot(4,4,1);
semilogx(f,Pdb_1_LRM(1,:)','g');
hold on;
semilogx(f,magTU11,'r');

subplot(4,4,2);
semilogx(f,Pdb_1_LRM(2,:)','g');
hold on;
semilogx(f,magTU12,'r');

subplot(4,4,3);
semilogx(f,Pdb_1_LRM(3,:)','g');
hold on;
semilogx(f,magTU13,'r');

subplot(4,4,4);
semilogx(f,Pdb_1_LRM(4,:)','g');
hold on;
semilogx(f,magTU14,'r');

subplot(4,4,5);
semilogx(f,Pdb_2_LRM(1,:)','g');
hold on;
semilogx(f,magTU21,'r');

subplot(4,4,6);
semilogx(f,Pdb_2_LRM(2,:)','g');
hold on;
semilogx(f,magTU22,'r');

subplot(4,4,7);
semilogx(f,Pdb_2_LRM(3,:)','g');
hold on;
semilogx(f,magTU23,'r');

subplot(4,4,8);
semilogx(f,Pdb_2_LRM(4,:)','g');
hold on;
semilogx(f,magTU24,'r');

subplot(4,4,9);
semilogx(f,Pdb_3_LRM(1,:)','g');
hold on;
semilogx(f,magTU31,'r');

subplot(4,4,10);
semilogx(f,Pdb_3_LRM(2,:)','g');
hold on;
semilogx(f,magTU32,'r');

subplot(4,4,11);
semilogx(f,Pdb_3_LRM(3,:)','g');
hold on;
semilogx(f,magTU33,'r');

subplot(4,4,12);
semilogx(f,Pdb_3_LRM(4,:)','g');
hold on;
semilogx(f,magTU34,'r');

subplot(4,4,13);
semilogx(f,Pdb_4_LRM(1,:)','g');
hold on;
semilogx(f,magTU41,'r');

subplot(4,4,14);
semilogx(f,Pdb_4_LRM(2,:)','g');
hold on;
semilogx(f,magTU42,'r');

subplot(4,4,15);
semilogx(f,Pdb_4_LRM(3,:)','g');
hold on;
semilogx(f,magTU43,'r');

subplot(4,4,16);
semilogx(f,Pdb_4_LRM(4,:)','g');
hold on;
semilogx(f,magTU44,'r');


