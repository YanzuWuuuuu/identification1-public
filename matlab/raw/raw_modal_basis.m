function [forced,initial]=raw_modal_basis(u,ts,frequency,damping,holdMethod)
% Time-domain convolution basis. No estimated FRF or output filtering.
assert(ismatrix(u)&&all(isfinite(u(:)))&&isreal(u)&&ts>0);
frequency=frequency(:); damping=damping(:); K=numel(frequency);
assert(numel(damping)==K && all(frequency>0) && all(damping>0 & damping<1));
N=size(u,1); ni=size(u,2); t=(0:N-1)'*ts;
forced=zeros(N,ni*(K+1)); initial=zeros(N,2+2*K);
initial(:,1:2)=[ones(N,1),t];
for k=0:K
    if k==0
        denominator=[1,0,0];
    else
        w=2*pi*frequency(k); denominator=[1,2*damping(k)*w,w^2];
        wd=w*sqrt(1-damping(k)^2);
        initial(:,2*k+(1:2))=exp(-damping(k)*w*t).*[cos(wd*t),sin(wd*t)];
    end
    discrete=c2d(tf(1,denominator),ts,holdMethod);
    [b,a]=tfdata(discrete,'v'); padded=zeros(1,numel(a));
    padded(end-numel(b)+1:end)=b;
    forced(:,k*ni+(1:ni))=filter(padded,a,u,[],1);
end
end
