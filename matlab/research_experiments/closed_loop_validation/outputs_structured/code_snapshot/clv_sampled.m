function [G,J,p]=clv_sampled(x,c,f,Ts)
wantJ=nargout>1;
if wantJ,[p,dm,dr]=clv_decode(x,c);else,p=clv_decode(x,c);end
z=exp(2i*pi*f(:).'*Ts);nf=numel(f);G=complex(zeros(16,nf));
if wantJ,J=complex(zeros(16,nf,c.n));else,J=[];end
ell=p.ty+p.tu.';
for edge=1:16
    [j,k]=ind2sub([4 4],edge);a=p.a(edge);b=p.b(edge);tau=ell(edge);
    [H,dH]=clv_kernel(0,0,a,b,tau,z,Ts,false,wantJ);G(edge,:)=p.M(edge)*H;
    if wantJ
        J(edge,:,1:15)=reshape(H.'*reshape(dm(j,k,:),1,15),1,nf,15);
        da=p.M(edge)*dH(:,:,3);db=p.M(edge)*dH(:,:,4);dt=p.M(edge)*dH(:,:,5);
    end
    for i=1:c.K
        [H,dH]=clv_kernel(p.omega(i),p.zeta(i),a,b,tau,z,Ts,true,wantJ);
        R=p.R(j,k,i);G(edge,:)=G(edge,:)+real(R)*H(1,:)+imag(R)*H(2,:);
        if wantJ
            for t=1:2
                ix=c.fo+i;if t==2,ix=c.zo+i;end
                J(edge,:,ix)=real(R)*dH(1,:,t)+imag(R)*dH(2,:,t);
            end
            ix=c.ro+14*(i-1)+(1:14);rho=reshape(dr(j,k,i,:),1,14);
            J(edge,:,ix)=reshape(H(1,:).'*real(rho)+H(2,:).'*imag(rho),1,nf,14);
            da=da+real(R)*dH(1,:,3)+imag(R)*dH(2,:,3);
            db=db+real(R)*dH(1,:,4)+imag(R)*dH(2,:,4);
            dt=dt+real(R)*dH(1,:,5)+imag(R)*dH(2,:,5);
        end
    end
    if wantJ
        if isfield(c,'additive_a')&&c.additive_a
            J(edge,:,c.ao+(1:7))=reshape(da.'*(c.as*c.a_map(edge,:)),1,nf,7);
        else,J(edge,:,c.ao+edge)=da*p.a(edge);end
        J(edge,:,c.bo+edge)=db;
        if j>1,J(edge,:,c.to+j-1)=dt*c.ts;end
        J(edge,:,c.to+3+k)=dt*c.ts;
    end
end
if wantJ,J=reshape(J,16*nf,c.n);end
end
