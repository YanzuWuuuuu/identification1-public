function out = joint_fullband_srivc(f, G, fpri, zpri, opts)
% Joint full-band additive SRIVC: G(s) = M/s^2 + sum_k B_k/(1+a1k s+a2k s^2)
% (rigid body = marginally stable subsystem with known denominator s^2, thesis Sec. 3.5.1)
% f: Hz (column), G: noni x N complex FRF, fpri/zpri: prior freq (Hz)/damping
% opts.colscale (true/false), opts.weight ('relative'|'none'), opts.maxit
if nargin<5, opts=struct; end
if ~isfield(opts,'colscale'), opts.colscale=true; end
if ~isfield(opts,'weight'), opts.weight='relative'; end
if ~isfield(opts,'maxit'), opts.maxit=60; end
if ~isfield(opts,'Mfixed'), opts.Mfixed=[]; end
fixM=~isempty(opts.Mfixed);
f=f(:); N=numel(f); noni=size(G,1); K=numel(fpri); s=1i*2*pi*f.';
switch opts.weight
    case 'relative', W=1./max(abs(G),1e-3*max(abs(G),[],2));   % data-based (state of the art)
    otherwise, W=ones(size(G));                                 % 'none' and initial 'model_relative'
end
modelW=strcmp(opts.weight,'model_relative');
a=[2*zpri(:)./(2*pi*fpri(:)), 1./(2*pi*fpri(:)).^2];       % K x 2
% --- initial numerators: fixed-pole weighted linear LS ---
A=1+a(:,1)*s+a(:,2)*s.^2;                                  % K x N
Phi=[1./s.^2; 1./A];                                       % (K+1) x N basis
theta=zeros(noni,K+1);
for r=1:noni
    X=(W(r,:).*Phi).'; y=(W(r,:).*G(r,:)).';
    Xr=[real(X);imag(X)]; yr=[real(y);imag(y)]; sc=vecnorm(Xr);
    if ~opts.colscale, sc=ones(size(sc)); end
    theta(r,:)=((Xr./sc)\yr).'./sc;
end
M=theta(:,1); B=theta(:,2:end);
if fixM, M=opts.Mfixed(:); for r=1:noni, X=(W(r,:).*Phi(2:end,:)).'; y=(W(r,:).*(G(r,:)-M(r)./s.^2)).'; Xr=[real(X);imag(X)]; yr=[real(y);imag(y)]; sc=vecnorm(Xr); B(r,:)=((Xr./sc)\yr).'./sc; end, end
p=noni+K*(2+noni); hist=zeros(opts.maxit,1); condhist=hist; condraw=hist;
for it=1:opts.maxit
    A=1+a(:,1)*s+a(:,2)*s.^2; P=zeros(noni,N,K);
    for k=1:K, P(:,:,k)=B(:,k)./A(k,:); end
    pred=M./s.^2+sum(P,3);
    if modelW, W=1./max(abs(pred),1e-3*max(abs(pred),[],2)); end
    X=zeros(2*N*noni,p); Z=X; Y=zeros(2*N*noni,K+1);
    rowsC=@(v) [real(v);imag(v)];
    % rigid block (linear, regressor = instrument)
    Gt0=G-sum(P,3); idx0=1:noni;
    for n=1:N
        rows=(n-1)*2*noni+(1:2*noni); w=W(:,n);
        x=zeros(noni,p); z=x; y=zeros(noni,K+1);
        if ~fixM, x(:,idx0)=diag(w)/s(n)^2; z(:,idx0)=x(:,idx0); y(:,1)=w.*Gt0(:,n); else, y(:,1)=0; end
        for k=1:K
            ik=noni+(k-1)*(2+noni)+(1:2+noni);
            Gtk=G(:,n)-pred(:,n)+P(:,n,k); Ak=A(k,n);
            x(:,ik)=w.*[-s(n)*Gtk/Ak, -s(n)^2*Gtk/Ak, eye(noni)/Ak];
            z(:,ik)=w.*[-s(n)*P(:,n,k)/Ak, -s(n)^2*P(:,n,k)/Ak, eye(noni)/Ak];
            y(:,k+1)=w.*Gtk/Ak;
        end
        X(rows,:)=rowsC(x); Z(rows,:)=rowsC(z); Y(rows,:)=rowsC(y);
    end
    if opts.colscale, sx=max(vecnorm(X),realmin); sz=max(vecnorm(Z),realmin); else, sx=ones(1,p); sz=sx; end
    if fixM, X(:,idx0)=[]; Z(:,idx0)=[]; sx(idx0)=[]; sz(idx0)=[]; end
    Mom=(Z./sz).'*(X./sx); condhist(it)=cond(Mom); condraw(it)=cond(Z.'*X);
    sol=(Mom\((Z./sz).'*Y))./sx.';
    if fixM, sol=[zeros(noni,K+1);sol]; sol(idx0,1)=M; end
    Mn=sol(idx0,1); an=a; Bn=B;
    for k=1:K
        ik=noni+(k-1)*(2+noni)+(1:2+noni); b=sol(ik,k+1);
        r=roots([b(2),b(1),1]); r=complex(-abs(real(r)),imag(r)); d=real(poly(r)); d=d/d(end);
        an(k,:)=[d(2),d(1)]; Bn(:,k)=b(3:end);
    end
    ch=norm([Mn-M;an(:)./a(:)-1;(Bn(:)-B(:))/max(norm(B(:)),eps)])/max(norm([M;ones(2*K,1)]),eps);
    M=Mn; a=an; B=Bn;
    An=1+a(:,1)*s+a(:,2)*s.^2; pr=M./s.^2; for k=1:K, pr=pr+B(:,k)./An(k,:); end
    hist(it)=norm(W.*(G-pr),'fro');
    if ch<1e-10, break; end
end
% recompute prediction cleanly
A=1+a(:,1)*s+a(:,2)*s.^2; pred=M./s.^2; for k=1:K, pred=pred+B(:,k)./A(k,:); end
wn=1./sqrt(a(:,2)); out=struct('M',M,'B',B,'a',a,'f_hz',wn/2/pi,'zeta',a(:,1).*wn/2, ...
    'iterations',it,'cond_moment',condhist(1:it),'cond_raw',condraw(1:it),'rel_err',norm(G-pred,'fro')/norm(G,'fro'), ...
    'rel_err_weighted',norm(W.*(G-pred),'fro')/norm(W.*G,'fro'),'pred',pred);
end
