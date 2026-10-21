function [p,dm,dr]=clv_decode(x,c)
A=reshape(x(1:9),3,3);b=x(10:12);r=x(13:15).'; v=A\b;w=r/A;
p.M=c.ms*c.U*[A,b;r,r*v]*c.V';p.D=zeros(4);
omega=2*pi*exp(x(c.fo+(1:c.K)));zeta=exp(x(c.zo+(1:c.K)));
p.lambda=-omega.*zeta+1i*omega.*sqrt(1-zeta.^2);p.omega=omega;p.zeta=zeta;
p.ty=[0;x(c.to+(1:3))]*c.ts;p.tu=x(c.to+(4:7))*c.ts;
if isfield(c,'additive_a')&&c.additive_a
    p.a=reshape(c.as*c.a_map*x(c.ao+(1:7)),4,4);
else,p.a=reshape(c.as*exp(x(c.ao+(1:16))),4,4);end
p.b=reshape(c.bs*exp(x(c.bo+(1:16))),4,4);
if nargout>1
    dm=zeros(4,4,15);dr=complex(zeros(4,4,c.K,14));
    for k=1:15
        d=zeros(4);
        if k<=9, E=zeros(3);E(k)=1;d(1:3,1:3)=E;d(4,4)=-w*E*v;
        elseif k<=12,e=zeros(3,1);e(k-9)=1;d(1:3,4)=e;d(4,4)=w*e;
        else,e=zeros(1,3);e(k-12)=1;d(4,1:3)=e;d(4,4)=e*v;end
        dm(:,:,k)=c.ms*c.U*d*c.V';
    end
end
for i=1:c.K
    ix=c.ro+14*(i-1)+(1:14); order=c.order(:,i);
    left=ones(4,1);left(order(2:4))=x(ix(1:3))+1i*x(ix(4:6));
    right=c.rs(i)*(x(ix(7:10))+1i*x(ix(11:14))).';p.R(:,:,i)=left*right;
    if nargout>1
        for k=1:14
            if k<=6
                e=zeros(4,1);e(order(2+mod(k-1,3)))=1;if k>3,e=1i*e;end
                dr(:,:,i,k)=e*right;
            else
                e=zeros(1,4);e(1+mod(k-7,4))=c.rs(i);if k>10,e=1i*e;end
                dr(:,:,i,k)=left*e;
            end
        end
    end
end
end
