
function [P_LRM,tf_magitude_y_LRM,tf_phase_y_LRM] = LRM_singleDOF(Y,U,n,B)

% B = 30;%%该参数为初始位置
nw = B-1;%%窗宽度
E = n-nw;%%该参数为结束位置
R = 4;%%拟合多项式阶数

P_LRM = zeros(n,1);
K_temp = zeros(2*nw+1,2*(R+1)+R);
for k = B:E
    Y_temp = Y(k-nw:k+nw);
    for r = -nw:nw
% R = 4
    K_temp(r+nw+1,:) = [U(k+r),r*U(k+r),r^2*U(k+r),r^3*U(k+r),r^4*U(k+r),1,r,r^2,r^3,r^4,-r*Y(k+r),-r^2*Y(k+r),-r^3*Y(k+r),-r^4*Y(k+r)];
% % R = 2
%     K_temp(r+nw+1,:) = [U(k+r),r*U(k+r),r^2*U(k+r),1,r,r^2,-r*Y(k+r),-r^2*Y(k+r)];
    end
    thetak =inv(K_temp'*K_temp)*K_temp'*Y_temp;
    P_LRM(B+k-1) = thetak(1);
end
tf_magitude_y_LRM=20*log10(abs(P_LRM)); 
tf_phase_y_LRM=rad2deg(angle(P_LRM));

for ii=1:length(tf_phase_y_LRM)
    if tf_phase_y_LRM(ii)>0
        tf_phase_y_LRM(ii)=tf_phase_y_LRM(ii)-360;
    end
end
% subplot(2,1,1);
% semilogx(f,tf_magitude_y_LRM,'g');
% hold on;
% grid on;
% xlabel('Frequency [Hz]','FontSize',10);
% ylabel('Magnitude (dB)','FontSize',10); set(gca,'FontSize',10);
% subplot(2,1,2); 
% semilogx(f, tf_phase_y_LRM,'g');hold on;
% grid on;
% xlabel('Frequency [Hz]','FontSize',10);
% ylabel('Phase[deg]','FontSize',10); set(gca,'FontSize',10);
end