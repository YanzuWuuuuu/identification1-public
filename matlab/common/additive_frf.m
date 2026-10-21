function [G,J]=additive_frf(beta,f,noni,K)
% beta=[vec(M);a1;a2;vec(B1);...], s in rad/s, f in Hz.
beta=beta(:); f=f(:);
assert(numel(beta)==noni+K*(noni+2),'Additive parameter size mismatch.');
assert(all(f>0 & isfinite(f)),'Evaluation frequencies must be positive.');
G=zeros(noni,numel(f)); J=zeros(2*noni*numel(f),numel(beta));
for k=1:numel(f)
    s=1i*2*pi*f(k); response=beta(1:noni)/s^2;
    derivative=zeros(noni,numel(beta)); derivative(:,1:noni)=eye(noni)/s^2;
    for j=1:K
        ix=noni+(j-1)*(noni+2)+(1:noni+2);
        A=1+beta(ix(1))*s+beta(ix(2))*s^2; P=beta(ix(3:end))/A;
        response=response+P;
        derivative(:,ix)=[-s*P/A,-s^2*P/A,eye(noni)/A];
    end
    G(:,k)=response;
    J((k-1)*2*noni+(1:2*noni),:)=[real(derivative);imag(derivative)];
end
end
