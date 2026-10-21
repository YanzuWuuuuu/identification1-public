function T = theory_check(action, varargin)
% Numerical checks of the derived theory for delay in additive modal identification.
% 'bias'   : first-order bias prediction  df/f = -zeta*sin(w*tau)  (isolated modes, real residues)
% 'info'   : profiled information for tau, I_tt = ||(I-P_theta) dG/dtau||^2, at true parameters
%            -> predicted J(0)/J(tau*) = 1 + tau*^2 I_tt /(sigma^2 N noni), compare with simulations
here=fileparts(mfilename('fullpath')); root=fileparts(here);
old=pwd; cd(fullfile(root,'3')); no=4; ni=4; model; cd(old); %#ok<NASGU>
ft=[286.43,344.58,488.27,705.67,738.83]; zt=0.01*ones(1,5);
Rk=arrayfun(@(k) reshape(C2(:,k)*B2(k,:),[],1),1:5,'uni',0); noni=16;
switch action
 case 'bias'
  tau=varargin{1};
  T=struct('pred_df_pct',-100*zt.*sin(2*pi*ft*tau));
 case 'linbias'
  % first-order bias of the delay-ignoring estimator: dtheta = (J'J)^-1 J' t,  t = -tau*s.*G0
  [tau,nmodes,fmax]=varargin{:};
  ft=ft(1:nmodes); zt=zt(1:nmodes); Rk=Rk(1:nmodes);
  f=(5:2:fmax)'; s=1i*2*pi*f.'; N=numel(f); G0=vecM./s.^2;
  for k=1:nmodes, G0=G0+Rk{k}./(s.^2+2*zt(k)*2*pi*ft(k)*s+(2*pi*ft(k))^2); end
  a1=2*zt./(2*pi*ft); a2=1./(2*pi*ft).^2;
  J=zeros(noni*N,noni+nmodes*(2+noni));
  for n=1:N, J((n-1)*noni+(1:noni),1:noni)=eye(noni)/s(n)^2; end
  for k=1:nmodes
    Bk=Rk{k}*a2(k); ik=noni+(k-1)*(2+noni)+(1:2+noni);
    for n=1:N, A=1+a1(k)*s(n)+a2(k)*s(n)^2; P=Bk/A;
      J((n-1)*noni+(1:noni),ik)=[-s(n)*P/A, -s(n)^2*P/A, eye(noni)/A]; end
  end
  t=-tau*(s.*G0); t=t(:);
  Jr=[real(J);imag(J)]; tr=[real(t);imag(t)]; sc=vecnorm(Jr); [Q,R]=qr(Jr./sc,0); d=(R\(Q.'*tr))./sc.';
  da1=d(noni+(0:nmodes-1)*(2+noni)+1).'; da2=d(noni+(0:nmodes-1)*(2+noni)+2).';
  T=struct('df_pct',-100*da2./(2*a2),'dz_pct',100*(da1./a1-da2./(2*a2)));
 case 'info'
  [band,numdeg,tau,snr_db]=varargin{:};
  f=(5:2:800)'; s=1i*2*pi*f.'; G0=vecM./s.^2;
  for k=1:5, G0=G0+Rk{k}./(s.^2+2*zt(k)*2*pi*ft(k)*s+(2*pi*ft(k))^2); end
  sigma=10^(-snr_db/20)*median(abs(G0(:,f>=200).*exp(-s(f>=200)*tau)),'all');  % same as data generator
  keep=f>=band(1)&f<=band(2); f=f(keep); s=s(keep); N=numel(f); G0=G0(:,keep);
  % model Jacobian (DC-normalised modes) at true parameters, real stacking [Re;Im] per frequency
  a1=2*zt./(2*pi*ft); a2=1./(2*pi*ft).^2; B=cellfun(@(R,w) R/w^2, Rk, num2cell(2*pi*ft),'uni',0);
  cols={}; 
  cols{end+1}=kron(eye(noni),1)*0; %#ok<NASGU>
  J=[]; 
  Jm=zeros(noni*N,noni); for n=1:N, Jm((n-1)*noni+(1:noni),:)=eye(noni)/s(n)^2; end; J=[J,Jm];
  for k=1:5
    Jk=zeros(noni*N,2+noni*(1+numdeg));
    for n=1:N
      A=1+a1(k)*s(n)+a2(k)*s(n)^2; P=B{k}/A;
      blk=[-s(n)*P/A, -s(n)^2*P/A, eye(noni)/A];
      if numdeg, blk=[blk, s(n)*eye(noni)/A]; end
      Jk((n-1)*noni+(1:noni),:)=blk;
    end
    J=[J,Jk]; %#ok<AGROW>
  end
  dG=-(s.*G0); dG=dG(:);                       % d/dtau [exp(-s tau) G0] at tau=0 (unit-modulus factor irrelevant)
  Jr=[real(J);imag(J)]; dr=[real(dG);imag(dG)];
  Jr_s=Jr./vecnorm(Jr); [Q,~]=qr(Jr_s,0); res=dr-Q*(Q.'*dr);
  Itt=norm(res)^2;                              % profiled (Schur-complement) information
  ratio=1+tau^2*Itt/(sigma^2*N*noni);
  T=struct('Itt',Itt,'ratio_pred',ratio,'N',N,'sigma',sigma,'frac_info_removed',1-norm(res)^2/norm(dr)^2);
end
end
