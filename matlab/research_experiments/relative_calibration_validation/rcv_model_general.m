function [G,J,p] = rcv_model_general(x,c,s)
%RCV_MODEL_GENERAL Four-channel model with any number of modal pole pairs.
% K = numel(c.rs), with 70+16*K real parameters. The first 31 parameters
% encode rank-three M and real D. Pole frequencies and damping use K each,
% rank-one complex residues use 14*K, delays use 7, and denominators use 32.
% This preserves the pilot's four-mode chart and analytic Jacobian ordering.
x=x(:); s=s(:).'; N=numel(s); wantJ=nargout>1;
nModes=numel(c.rs); nParameters=70+16*nModes;
assert(numel(x)==nParameters,'Expected %d parameters for %d modes.',nParameters,nModes);
poleOffset=31; dampingOffset=31+nModes;
residueOffset=31+2*nModes;
delayOffset=31+16*nModes;
aOffset=delayOffset+7; bOffset=delayOffset+23;

A=reshape(x(1:9),3,3); b=x(10:12); row=x(13:15).';
v=A\b; w=row/A;
rigidCore=[A,b;row,row*v]; p.M=c.ms*c.U*rigidCore*c.V';
p.D=c.ds*reshape(x(16:31),4,4);
omega=2*pi*exp(x(poleOffset+(1:nModes)));
zeta=exp(x(dampingOffset+(1:nModes)));
p.lambda=-zeta.*omega+1i*omega.*sqrt(1-zeta.^2);
p.R=complex(zeros(4,4,nModes));
P=p.M(:)*(1./s.^2)+p.D(:)*ones(1,N);
if wantJ, JP=complex(zeros(16,N,nParameters)); end
if wantJ
    for k=1:15
        dK=zeros(4);
        if k<=9
            E=zeros(3); E(k)=1; dK(1:3,1:3)=E; dK(4,4)=-w*E*v;
        elseif k<=12
            e=zeros(3,1); e(k-9)=1; dK(1:3,4)=e; dK(4,4)=w*e;
        else
            e=zeros(1,3); e(k-12)=1; dK(4,1:3)=e; dK(4,4)=e*v;
        end
        dM=c.ms*c.U*dK*c.V'; JP(:,:,k)=dM(:)*(1./s.^2);
    end
    for k=1:16, JP(k,:,15+k)=c.ds; end
end

for i=1:nModes
    ix=residueOffset+(i-1)*14+(1:14); order=c.order(:,i);
    left=ones(4,1); left(order(2:4))=x(ix(1:3))+1i*x(ix(4:6));
    right=c.rs(i)*(x(ix(7:10))+1i*x(ix(11:14))).';
    R=left*right; p.R(:,:,i)=R; lam=p.lambda(i);
    f=1./(s-lam); fc=1./(s-conj(lam));
    P=P+R(:)*f+conj(R(:))*fc;
    if wantJ
        dl=[lam,-omega(i)*zeta(i)-1i*omega(i)*zeta(i)^2/sqrt(1-zeta(i)^2)];
        JP(:,:,poleOffset+i)=R(:)*(dl(1)*f.^2)+conj(R(:))*(conj(dl(1))*fc.^2);
        JP(:,:,dampingOffset+i)=R(:)*(dl(2)*f.^2)+conj(R(:))*(conj(dl(2))*fc.^2);
        for k=1:14
            if k<=6
                e=zeros(4,1); e(order(2+mod(k-1,3)))=1;
                if k>3, e=1i*e; end
                dR=e*right;
            else
                e=zeros(1,4); e(1+mod(k-7,4))=c.rs(i);
                if k>10, e=1i*e; end
                dR=left*e;
            end
            JP(:,:,ix(k))=dR(:)*f+conj(dR(:))*fc;
        end
    end
end

p.ty=[0;x(delayOffset+(1:3))]*c.ts;
p.tu=x(delayOffset+(4:7))*c.ts;
p.a=reshape(exp(x(aOffset+(1:16)))*c.as,4,4);
p.b=reshape(exp(x(bOffset+(1:16)))*c.bs,4,4);
ell=p.ty+p.tu.'; den=1+p.a(:)*s+p.b(:)*s.^2;
H=exp(-ell(:)*s)./den; G=P.*H;
if wantJ
    JP=JP.*H;
    for j=2:4
        mask=zeros(4); mask(j,:)=1;
        JP(:,:,delayOffset+j-1)=-G.*(mask(:)*(c.ts*s));
    end
    for k=1:4
        mask=zeros(4); mask(:,k)=1;
        JP(:,:,delayOffset+3+k)=-G.*(mask(:)*(c.ts*s));
    end
    for k=1:16
        JP(k,:,aOffset+k)=-G(k,:).*(p.a(k)*s)./den(k,:);
        JP(k,:,bOffset+k)=-G(k,:).*(p.b(k)*s.^2)./den(k,:);
    end
    J=reshape(JP,16*N,nParameters);
end
end
