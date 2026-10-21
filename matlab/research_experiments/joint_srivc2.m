function out = joint_srivc2(f, G, fpri, zpri, opts)
% Joint full-band additive SRIVC (frequency domain), rigid body M/s^2 + K second-order modes.
%   mode k: (B0k + B1k*s)/(1 + a1k s + a2k s^2)   (B1 only if opts.numdeg=1: general damping)
% opts.weight  : 'none' | 'relative'(data-based 1/|G|)
% opts.solver  : 'moment' (normal-equation form Z'X, as in thesis Remark 6) | 'qr' (unsquared)
% opts.colscale: true/false ; opts.numdeg: 0/1 ; opts.maxit ; opts.tol
d=struct('weight','none','solver','moment','colscale',true,'numdeg',0,'maxit',80,'tol',1e-10);
fn=fieldnames(d); for i=1:numel(fn), if ~isfield(opts,fn{i}), opts.(fn{i})=d.(fn{i}); end, end
f=f(:); N=numel(f); s=1i*2*pi*f.'; noni=size(G,1); K=numel(fpri); nd=opts.numdeg;
nb=2+noni*(1+nd); p=noni+K*nb;
if strcmp(opts.weight,'relative'), W=1./max(abs(G),1e-3*max(abs(G),[],2)); else, W=ones(size(G)); end
a=[2*zpri(:)./(2*pi*fpri(:)), 1./(2*pi*fpri(:)).^2];
% initial numerators by fixed-pole LS
A=1+a(:,1)*s+a(:,2)*s.^2; Phi=1./s.^2; for k=1:K, Phi=[Phi;1./A(k,:)]; if nd, Phi=[Phi; s./A(k,:)]; end, end %#ok<AGROW>
th=zeros(noni,size(Phi,1));
for r=1:noni
  X=(W(r,:).*Phi).'; y=(W(r,:).*G(r,:)).'; Xr=[real(X);imag(X)]; sc=vecnorm(Xr);
  th(r,:)=((Xr./sc)\[real(y);imag(y)]).'./sc;
end
M=th(:,1); B0=zeros(noni,K); B1=zeros(noni,K);
for k=1:K, c=1+(k-1)*(1+nd)+1; B0(:,k)=th(:,c); if nd, B1(:,k)=th(:,c+1); end, end
condhist=zeros(opts.maxit,1);
for it=1:opts.maxit
  A=1+a(:,1)*s+a(:,2)*s.^2; P=zeros(noni,N,K);
  for k=1:K, P(:,:,k)=(B0(:,k)+B1(:,k).*s)./A(k,:); end
  pred=M./s.^2+sum(P,3);
  X=zeros(2*N*noni,p); Z=X; Y=zeros(2*N*noni,K+1); I=eye(noni);
  for n=1:N
    rows=(n-1)*2*noni+(1:2*noni); w=W(:,n); sn=s(n);
    x=zeros(noni,p); z=x; y=zeros(noni,K+1);
    x(:,1:noni)=diag(w)/sn^2; z(:,1:noni)=x(:,1:noni); y(:,1)=w.*(G(:,n)-sum(P(:,n,:),3));
    for k=1:K
      ik=noni+(k-1)*nb+(1:nb); Ak=A(k,n); Gt=G(:,n)-pred(:,n)+P(:,n,k);
      numcols=I/Ak; if nd, numcols=[numcols, sn*I/Ak]; end %#ok<AGROW>
      x(:,ik)=w.*[-sn*Gt/Ak, -sn^2*Gt/Ak, numcols];
      z(:,ik)=w.*[-sn*P(:,n,k)/Ak, -sn^2*P(:,n,k)/Ak, numcols];
      y(:,k+1)=w.*Gt/Ak;
    end
    X(rows,:)=[real(x);imag(x)]; Z(rows,:)=[real(z);imag(z)]; Y(rows,:)=[real(y);imag(y)];
  end
  if opts.colscale, sx=max(vecnorm(X),realmin); sz=max(vecnorm(Z),realmin); else, sx=ones(1,p); sz=sx; end
  Xs=X./sx; Zs=Z./sz;
  switch opts.solver
    case 'moment', Mom=Zs.'*Xs; condhist(it)=cond(Mom); sol=(Mom\(Zs.'*Y))./sx.';
    case 'qr', [Q,~]=qr(Zs,0); QX=Q.'*Xs; condhist(it)=cond(QX); sol=(QX\(Q.'*Y))./sx.';
  end
  Mn=sol(1:noni,1); an=a; B0n=B0; B1n=B1;
  for k=1:K
    b=sol(noni+(k-1)*nb+(1:nb),k+1);
    r=roots([b(2),b(1),1]); r=complex(-abs(real(r)),imag(r)); dd=real(poly(r)); dd=dd/dd(end);
    an(k,:)=[dd(2),dd(1)]; B0n(:,k)=b(3:2+noni); if nd, B1n(:,k)=b(3+noni:end); end
  end
  ch=norm([Mn-M; an(:)./a(:)-1])/max(norm([M;ones(2*K,1)]),eps);
  M=Mn; a=an; B0=B0n; B1=B1n;
  if ch<opts.tol, break; end
end
A=1+a(:,1)*s+a(:,2)*s.^2; pred=M./s.^2; for k=1:K, pred=pred+(B0(:,k)+B1(:,k).*s)./A(k,:); end
wn=1./sqrt(a(:,2));
out=struct('M',M,'B0',B0,'B1',B1,'a',a,'f_hz',wn/2/pi,'zeta',a(:,1).*wn/2,'pred',pred, ...
  'iterations',it,'cond',condhist(1:it),'cost',norm(G-pred,'fro')^2);
end
