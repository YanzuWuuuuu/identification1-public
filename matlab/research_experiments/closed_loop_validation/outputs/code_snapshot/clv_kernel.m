function [H,dH]=clv_kernel(w,zeta,a,b,delay,z,Ts,modal,wantJ)
% Exact ZOH frequency response of a mechanical second-order block followed
% by 1/(1+a*s+b*s^2), including an arbitrary nonnegative fractional delay.
% Matrix exponentials remain regular at repeated filter roots.
% dH coordinates: log(w), log(zeta), log(a), log(b), delay in seconds.
wr=2*pi*1000; A=zeros(4); A(1:2,1:2)=[0 wr;-1/(b*wr) -a/b];
A(3:4,3:4)=[0 wr;-w*w/wr -2*zeta*w]; A(4,1)=1/wr;
B=[0;1/(b*wr);0;0]; E=[A B;zeros(1,5)];
if modal
    root=sqrt(1-zeta^2); C=[0 0 2*zeta*w 2*wr;0 0 -2*w*root 0];
else, C=[0 0 1 0];end
m=floor(max(delay,0)/Ts); h=Ts-(delay-m*Ts);
FT=expm(E*Ts); if abs(h-Ts)<1e-15,Fh=FT;else,Fh=expm(E*h);end
Ad=FT(1:4,1:4); Bd=FT(1:4,5); B0=Fh(1:4,5);
n=numel(z); zp=reshape(z,1,1,n); phase=zp.^(-m);
R=pagemldivide(zp.*eye(4)-Ad,repmat(eye(4),1,1,n));
X=pagemtimes(R,(B0+reshape(Bd-B0,4,1,1)./zp).*phase);
H=reshape(pagemtimes(C,X),size(C,1),n);
if ~wantJ, dH=[];return;end
dH=complex(zeros(size(C,1),n,5));
for k=1:4
    if ~modal && k<=2,continue;end
    dE=zeros(5);dC=zeros(size(C));
    switch k
        case 1
            dE(4,3)=-2*w*w/wr;dE(4,4)=-2*zeta*w;
            dC(1,3)=2*zeta*w;dC(2,3)=-2*w*root;
        case 2
            dE(4,4)=-2*zeta*w;dC(1,3)=2*zeta*w;dC(2,3)=2*w*zeta^2/root;
        case 3,dE(2,2)=-a/b;
        case 4,dE(2,1)=1/(b*wr);dE(2,2)=a/b;dE(2,5)=-1/(b*wr);
    end
    DT=frechet(E,dE,Ts);if abs(h-Ts)<1e-15,Dh=DT;else,Dh=frechet(E,dE,h);end
    dAd=DT(1:4,1:4);dBd=DT(1:4,5);dB0=Dh(1:4,5);
    rhs=pagemtimes(dAd,X)+(dB0+reshape(dBd-dB0,4,1,1)./zp).*phase;
    dX=pagemtimes(R,rhs);
    dH(:,:,k)=reshape(pagemtimes(C,dX)+pagemtimes(dC,X),size(C,1),n);
end
dB0=-Fh(1:4,1:4)*B;
dX=pagemtimes(R,(dB0-reshape(dB0,4,1,1)./zp).*phase);
dH(:,:,5)=reshape(pagemtimes(C,dX),size(C,1),n);
end
function D=frechet(E,dE,t)
F=expm([E dE;zeros(5) E]*t);D=F(1:5,6:10);
end
