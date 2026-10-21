function [ls2,f,tf_magitude_y_ploy1,tf_phase_y_ploy1,tf_magitude_y_ploy2,tf_phase_y_ploy2,tf_magitude_y_ploy3,tf_phase_y_ploy3] = FunctionMOSO_FFT_LRM_X(u1,u2,u3,y,B)

fs=5000;
%%
NFFT=(length(u1));
ny=fft(y,NFFT); 
nu1=fft(u1,NFFT);
nu2=fft(u2,NFFT);
nu3=fft(u3,NFFT);
f=fs/2*linspace(0,1,NFFT/2+1);  
% %%
% Abs_P1=(ny(1:NFFT/2+1))./(nu1(1:NFFT/2+1)); tf_magitude_y1=20*log10(abs(Abs_P1)); tf_phase_y1=angle(Abs_P1);
% for i=1:length(tf_phase_y1)
%     if tf_phase_y1(i)>pi/2
%         tf_phase_y1(i)=tf_phase_y1(i)-2*pi;
%     end
% end
% tf_phase_y1=tf_phase_y1.*180./pi;
% %
% Abs_P2=(ny(1:NFFT/2+1))./(nu2(1:NFFT/2+1)); tf_magitude_y2=20*log10(abs(Abs_P2)); tf_phase_y2=angle(Abs_P2);
% for i=1:length(tf_phase_y2)
%     if tf_phase_y2(i)>pi/2
%         tf_phase_y2(i)=tf_phase_y2(i)-2*pi;
%     end
% end
% tf_phase_y2=tf_phase_y2.*180./pi;
% %
% Abs_P3=(ny(1:NFFT/2+1))./(nu1(1:NFFT/2+1)); tf_magitude_y3=20*log10(abs(Abs_P3)); tf_phase_y3=angle(Abs_P3);
% for i=1:length(tf_phase_y3)
%     if tf_phase_y3(i)>pi/2
%         tf_phase_y3(i)=tf_phase_y3(i)-2*pi;
%     end
% end
% tf_phase_y3=tf_phase_y3.*180./pi;
%%
%B; %数据初始拟合位置  
A=NFFT/2-B; %数据最后拟合位置  
abs_y=(ny(1:NFFT/2+1)); abs_x1=(nu1(1:NFFT/2+1)); abs_x2=(nu2(1:NFFT/2+1)); abs_x3=(nu3(1:NFFT/2+1));
%nw
nw=B-1; %nw为左右选取数据宽度，不包括中间处。 nw=B-2; nw=B-3;都行，最小二乘所堆积的数据为2*nw+1。
 j=B-1; %j为存取数据表格初始行
 
 for k=B:A
     j=j+1;
     Y=abs_y(k-nw:k+nw);       
     j1=0;
     for r=-nw:1:nw
%    R=4        
%        varphi=[abs_x1(k+r), r*abs_x1(k+r), r^2*abs_x1(k+r),r^3*abs_x1(k+r),r^4*abs_x1(k+r),1,r,r^2,r^3,r^4,...  
%            abs_x2(k+r), r*abs_x2(k+r), r^2*abs_x2(k+r),r^3*abs_x2(k+r),r^4*abs_x2(k+r),1,r,r^2,r^3,r^4,...
%            abs_x3(k+r), r*abs_x3(k+r),
%            r^2*abs_x3(k+r),r^3*abs_x3(k+r),r^4*abs_x3(k+r),1,r,r^2,r^3,r^4,...
%            -r*abs_y(k+r), -r^2*abs_y(k+r),-r^3*abs_y(k+r),-r^4*abs_y(k+r)];
%   R=2
        varphi=[abs_x1(k+r), r*abs_x1(k+r), r^2*abs_x1(k+r),1,r,r^2,...  
          abs_x2(k+r), r*abs_x2(k+r), r^2*abs_x2(k+r),1,r,r^2,...
          abs_x3(k+r), r*abs_x3(k+r), r^2*abs_x3(k+r),1,r,r^2,...
          -r*abs_y(k+r), -r^2*abs_y(k+r)];
      j1=j1+1;
      Phi(j1,:)=varphi;
    end
    par=Phi\Y; % par=(inv(Phi'*Phi))*Phi'*Y;
    ls3(j,:)=[j,par'];
%ls2;
ls2=ls3;
end
  
%%
%此部分为上面提供的数据，用梯度在寻优一次，结果差别不

% jj=B-1;
% initial_para1=ls2(:,2:7);  initial_para2=ls2(:,8:13); initial_para3=ls2(:,14:19); 
% for k=B:A
% Y=abs_y(k-nw:k+nw);  
% jj=jj+1;
% par4=initial_para1(jj,1:6)'; par5=initial_para2(jj,1:6)';  par6=initial_para3(jj,1:6)';
%  par1=par4;  par2=par5; par3=par6;
%   for kk=1:50
%        j1=0;j2=0;j3=0;
%       
%       for r=-nw:1:nw
%          varphi1=[abs_x1(k+r), r*abs_x1(k+r), r^2*abs_x1(k+r),1,r,r^2];
%          j1=j1+1;
%         Phi1(j1,:)=varphi1;
%       end
%       for r=-nw:1:nw
%          varphi2=[abs_x2(k+r), r*abs_x2(k+r), r^2*abs_x2(k+r),1,r,r^2];
%         j2=j2+1;
%         Phi2(j2,:)=varphi2;
%       end
%       for r=-nw:1:nw
%          varphi3=[abs_x3(k+r), r*abs_x3(k+r), r^2*abs_x3(k+r),1,r,r^2];
%          j3=j3+1;
%          Phi3(j3,:)=varphi3;
%       end
%     par1=par1+Phi1'*(Y-Phi1*par4-Phi2*par5-Phi3*par6)*1.5/max(eig(Phi1'*Phi1));
%     par2=par2+Phi2'*(Y-Phi1*par4-Phi2*par5-Phi3*par6)*1.5/max(eig(Phi2'*Phi2));
%     par3=par3+Phi3'*(Y-Phi1*par4-Phi2*par5-Phi3*par6)*1.5/max(eig(Phi3'*Phi3));
% %      par1=Phi1\(Y-Phi1*par4-Phi2*par5-Phi3*par6);
% %      par2=Phi2\(Y-Phi1*par4-Phi2*par5-Phi3*par6);
% %      par3=Phi3\(Y-Phi1*par4-Phi2*par5-Phi3*par6);
%      par4=par1;  par5=par2; par6=par3; 
%   end
%       ls3(jj,:)=[jj,par4',par5',par6'];
% end

%%
tf_magitude_y_ploy1=20*log10(abs(ls3(:,2)));
tf_phase_y_ploy1=angle(ls3(:,2));
for ii=1:length(tf_phase_y_ploy1)
    if tf_phase_y_ploy1(ii)>pi/2
        tf_phase_y_ploy1(ii)=tf_phase_y_ploy1(ii)-2*pi;
    end
end
tf_phase_y_ploy1=tf_phase_y_ploy1.*180./pi;

%  此部分在不同功能可选择注释
for ii=1:length(tf_phase_y_ploy1)
    if tf_phase_y_ploy1(ii)>0
        tf_phase_y_ploy1(ii)=tf_phase_y_ploy1(ii)-2.*tf_phase_y_ploy1(ii);
    end
    if tf_phase_y_ploy1(ii)<-180
        tf_phase_y_ploy1(ii)=tf_phase_y_ploy1(ii)+2.*(-180-tf_phase_y_ploy1(ii));
    end
end
%%
tf_magitude_y_ploy2=20*log10(abs(ls3(:,8)));
tf_phase_y_ploy2=angle(ls3(:,8));
for ii=1:length(tf_phase_y_ploy2)
    if tf_phase_y_ploy2(ii)>pi/2
        tf_phase_y_ploy2(ii)=tf_phase_y_ploy2(ii)-2*pi;
    end
end
tf_phase_y_ploy2=tf_phase_y_ploy2.*180./pi;

%  此部分在不同功能可选择注释
for ii=1:length(tf_phase_y_ploy2)
    if tf_phase_y_ploy2(ii)>0
        tf_phase_y_ploy2(ii)=tf_phase_y_ploy2(ii)-2.*tf_phase_y_ploy2(ii);
    end
    if tf_phase_y_ploy2(ii)<-180
        tf_phase_y_ploy2(ii)=tf_phase_y_ploy2(ii)+2.*(-180-tf_phase_y_ploy2(ii));
    end
end
%%
tf_magitude_y_ploy3=20*log10(abs(ls3(:,14)));
tf_phase_y_ploy3=angle(ls3(:,14));
for ii=1:length(tf_phase_y_ploy3)
    if tf_phase_y_ploy3(ii)>pi/2
        tf_phase_y_ploy3(ii)=tf_phase_y_ploy3(ii)-2*pi;
    end
end
tf_phase_y_ploy3=tf_phase_y_ploy3.*180./pi;

%  此部分在不同功能可选择注释
for ii=1:length(tf_phase_y_ploy3)
    if tf_phase_y_ploy3(ii)>0
        tf_phase_y_ploy3(ii)=tf_phase_y_ploy3(ii)-2.*tf_phase_y_ploy3(ii);
    end
    if tf_phase_y_ploy3(ii)<-180
        tf_phase_y_ploy3(ii)=tf_phase_y_ploy3(ii)+2.*(-180-tf_phase_y_ploy3(ii));
    end
end

end