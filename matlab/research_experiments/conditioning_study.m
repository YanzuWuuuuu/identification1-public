function T = conditioning_study(fmodes, zeta, band, df, weight, param)
% Conditioning of one additive SRIVC/IV iteration at the true parameters (SISO-equivalent).
% Regressor blocks: rigid 1/s^2 ; mode i: [-s*Gi/Ai, -s^2*Gi/Ai, 1/Ai]  (Gi = r_i/Ai)
% param: 'dc' -> Ai = 1 + a1 s + a2 s^2 ; 'monic' -> Ai = s^2 + a1 s + a0 (then 1/Ai numerators)
% Returns cond of raw X, column-scaled X, and of the IV moment Z'X (Z=X at truth, noise-free),
% plus the numerical error of one IV solve: normal-equation LU vs QR (unsquared) vs SVD.
if nargin<6, param='dc'; end
f=(band(1):df:band(2))'; s=1i*2*pi*f; K=numel(fmodes); w0=2*pi*fmodes(:)';
rng(7); r=0.2+0.4*rand(1,K); m=0.8;
switch param
  case 'dc',    A=1+(2*zeta(:)'./w0).*s+(1./w0.^2).*s.^2; num=r./w0.^2;   % residue scaled by 1/w^2
  case 'monic', A=s.^2+(2*zeta(:)'.*w0).*s+w0.^2;          num=r;
end
Gi=num./A; G=m./s.^2+sum(Gi,2);
switch weight, case 'relative', W=1./abs(G); otherwise, W=ones(size(G)); end
X=zeros(numel(f),1+3*K); X(:,1)=1./s.^2;
for i=1:K
  X(:,1+3*(i-1)+(1:3))=[-s.*Gi(:,i)./A(:,i), -s.^2.*Gi(:,i)./A(:,i), 1./A(:,i)];
end
if strcmp(param,'monic')   % monic: derivative wrt a1,a0 and numerator
  for i=1:K, X(:,1+3*(i-1)+(1:3))=[-s.*Gi(:,i)./A(:,i), -Gi(:,i)./A(:,i), 1./A(:,i)]; end
end
X=W.*X; Xr=[real(X);imag(X)];
sc=vecnorm(Xr); Xs=Xr./sc;
% "true" parameters in this linearised IV form are such that X*theta = G-type target; build target exactly
theta=randn(size(Xr,2),1); y=Xs*theta;
t_ne=(Xs.'*Xs)\(Xs.'*y);                 % normal-equation/moment form (what the code does)
[Q,R]=qr(Xs,0); t_qr=R\(Q.'*y);          % orthogonal (unsquared) solve
T=struct('nf',numel(f),'cond_raw',cond(Xr),'cond_scaled',cond(Xs),'cond_moment_scaled',cond(Xs.'*Xs), ...
  'err_moment',norm(t_ne-theta)/norm(theta),'err_qr',norm(t_qr-theta)/norm(theta));
end
