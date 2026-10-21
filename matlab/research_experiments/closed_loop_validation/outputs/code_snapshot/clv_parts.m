function [G,Gr,Gm]=clv_parts(p,f,Ts)
nf=numel(f);z=exp(2i*pi*f(:).'*Ts);ell=p.ty+p.tu.';
Gr=complex(zeros(16,nf));Gm=complex(zeros(16,nf,numel(p.lambda)));
for e=1:16
    H=clv_kernel(0,0,p.a(e),p.b(e),ell(e),z,Ts,false,false);Gr(e,:)=p.M(e)*H;
    for i=1:numel(p.lambda)
        w=abs(p.lambda(i));zz=-real(p.lambda(i))/w;
        H=clv_kernel(w,zz,p.a(e),p.b(e),ell(e),z,Ts,true,false);R=p.R(:,:,i);R=R(e);
        Gm(e,:,i)=real(R)*H(1,:)+imag(R)*H(2,:);
    end
end
G=Gr+sum(Gm,3);
end
