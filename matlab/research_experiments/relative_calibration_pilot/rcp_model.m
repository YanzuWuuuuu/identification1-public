function [G,J,p] = rcp_model(x,c,s)
% Structured local chart, physical SI units. G is 16-by-N.
s=s(:).'; N=numel(s); wantJ=nargout>1;
A=reshape(x(1:9),3,3); b=x(10:12); row=x(13:15).';
v=A\b; w=row/A;
K=[A,b;row,row*v]; p.M=c.ms*c.U*K*c.V';
p.D=c.ds*reshape(x(16:31),4,4);
omega=2*pi*exp(x(32:35)); zeta=exp(x(36:39));
p.lambda=-zeta.*omega+1i*omega.*sqrt(1-zeta.^2);
p.R=complex(zeros(4,4,4));
P=p.M(:)*(1./s.^2)+p.D(:)*ones(1,N);
if wantJ, JP=complex(zeros(16,N,134)); end
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
for i=1:4
    ix=39+(i-1)*14+(1:14); order=c.order(:,i);
    left=ones(4,1); left(order(2:4))=x(ix(1:3))+1i*x(ix(4:6));
    right=c.rs(i)*(x(ix(7:10))+1i*x(ix(11:14))).';
    R=left*right; p.R(:,:,i)=R; lam=p.lambda(i);
    f=1./(s-lam); fc=1./(s-conj(lam));
    P=P+R(:)*f+conj(R(:))*fc;
    if wantJ
        dl=[lam, -omega(i)*zeta(i)-1i*omega(i)*zeta(i)^2/sqrt(1-zeta(i)^2)];
        JP(:,:,31+i)=R(:)*(dl(1)*f.^2)+conj(R(:))*(conj(dl(1))*fc.^2);
        JP(:,:,35+i)=R(:)*(dl(2)*f.^2)+conj(R(:))*(conj(dl(2))*fc.^2);
        for k=1:14
            if k<=6
                e=zeros(4,1); e(order(1+mod(k-1,3)+1))=1; % indices 2..4 in order
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
p.ty=[0;x(96:98)]*c.ts; p.tu=x(99:102)*c.ts;
p.a=reshape(exp(x(103:118))*c.as,4,4);
p.b=reshape(exp(x(119:134))*c.bs,4,4);
ell=p.ty+p.tu.'; den=1+p.a(:)*s+p.b(:)*s.^2;
H=exp(-ell(:)*s)./den; G=P.*H;
if wantJ
    JP=JP.*H;
    for j=2:4
        mask=zeros(4); mask(j,:)=1;
        JP(:,:,94+j)=-G.*(mask(:)*(c.ts*s));
    end
    for k=1:4
        mask=zeros(4); mask(:,k)=1;
        JP(:,:,98+k)=-G.*(mask(:)*(c.ts*s));
    end
    for k=1:16
        JP(k,:,102+k)=-G(k,:).*(p.a(k)*s)./den(k,:);
        JP(k,:,118+k)=-G(k,:).*(p.b(k)*s.^2)./den(k,:);
    end
    J=reshape(JP,16*N,134);
end
end
