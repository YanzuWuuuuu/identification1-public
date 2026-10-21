Mx=60;       %质量
My=146;      %质量
Jz=30;       %Rz转动惯量

s = tf('s');
Jrx = 0.0439;
Jry = 0.0384;
mz = 4.2809;

phai_body = [1/(Jrx*s^2),0,0;
             0,1/(Jry*s^2),0;
             0,0,1/(mz*s^2)];

kesai = [0.01,0.01,0.01,0.01,0.01];
ff = [286.43,344.58,488.27,705.67,738.83];
% ff = [286.43,344.58,488.27,705.67];
phai2 = tf(length(kesai),length(kesai));
for i = 1:length(kesai)
    phai2(i,i) = 1/(s^2+2*kesai(i)*(2*pi*ff(i))*s+(2*pi*ff(i))^2);
end
phai_flx = [phai2(1,1),0,0,0,0;
            0,phai2(2,2),0,0,0;
            0,0,phai2(3,3),0,0;
            0,0,0,phai2(4,4),0;
            0,0,0,0,phai2(5,5)];

phai = [phai_body,zeros(3,5);
        zeros(5,3),phai_flx];

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

P44 = C*phai*B;

vecM = C1*[1/Jrx,0,0;
           0,1/Jry,0;
           0,0,1/mz]*B1;
vecM = vecM(:);

vecB10 = C2(:,1)*B2(1,:);
vecB10 = vecB10(:)/(2*pi*ff(1))^2;

vecB20 = C2(:,2)*B2(2,:);
vecB20 = vecB20(:)/(2*pi*ff(2))^2;

vecB30 = C2(:,3)*B2(3,:);
vecB30 = vecB30(:)/(2*pi*ff(3))^2;

vecB40 = C2(:,4)*B2(4,:);
vecB40 = vecB40(:)/(2*pi*ff(4))^2;

vecB50 = C2(:,5)*B2(5,:);
vecB50 = vecB50(:)/(2*pi*ff(5))^2;

vecD = zeros(no*ni,1);


% belta_real = [vecM;vecB10;vecB20;vecB30;vecB40;vecB50];

vecB{1} = vecB10;
vecB{2} = vecB20;
vecB{3} = vecB30;
vecB{4} = vecB40;
vecB{5} = vecB40;
zeta_real = 0.01;

for i = 1:length(ff)
    B_real{i} = vecB{i};
    theta_real{i} = [1 zeta_real, B_real{i}']'/((2*pi*ff(i))^2);
    den_para_real{i} =  [1 2*zeta_real*(2*pi*ff(i))]/((2*pi*ff(i))^2);
end


theta_real1 = [flip(den_para_real{1}');vecB10];
theta_real2 = [flip(den_para_real{2}');vecB20];
theta_real3 = [flip(den_para_real{3}');vecB30];
theta_real4 = [flip(den_para_real{4}');vecB40];
theta_real5 = [flip(den_para_real{5}');vecB50];

belta_real = blkdiag(theta_real1,theta_real2,theta_real3,theta_real4,theta_real5);