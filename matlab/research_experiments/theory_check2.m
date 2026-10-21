function T = theory_check2(action, varargin)
% Further numerical checks for the derivation document.
% 'stationarity' : SRIVC fixed point == stationary point of weighted OE cost (Lemma 2)
% 'decomp'       : first-order damping/frequency bias decomposed by source of the delay effect
% 'info_fmin'    : profiled information I_tt versus lowest frequency line used
here=fileparts(mfilename('fullpath')); root=fileparts(here);
[vecM,C2,B2]=get_truth(); %#ok<ASGLU>
ft=[286.43,344.58,488.27,705.67,738.83]; zt=0.01*ones(1,5); K=5; noni=16;
Rk=arrayfun(@(k) reshape(C2(:,k)*B2(k,:),[],1),1:5,'uni',0);
a1=2*zt./(2*pi*ft); a2=1./(2*pi*ft).^2; B0=cell2mat(cellfun(@(R,w) R/w^2,Rk,num2cell(2*pi*ft),'uni',0));
switch action
 case 'stationarity'
  [tau,snr_db,seed,solver]=varargin{:};
  f=(5:2:800)'; s=1i*2*pi*f.'; G0=truth(s);
  sigma=10^(-snr_db/20)*median(abs(G0(:,f>=200)),'all'); rng(seed);
  G=G0.*exp(-s*tau)+sigma*(randn(size(G0))+1i*randn(size(G0)))/sqrt(2);
  Gc=G.*exp(s*tau);
  o=joint_srivc2(f,Gc,[280,355,480,690,765],0.02*ones(1,5),struct('solver',solver,'maxit',200,'tol',1e-13));
  J=jac(s,o.M,o.B0,o.a); e=Gc-o.pred; e=e(:);
  Jr=[real(J);imag(J)]; er=[real(e);imag(e)];
  g=Jr.'*er;                                  % gradient of 0.5*||e||^2 (up to sign)
  cosang=abs(g)./(vecnorm(Jr).'*norm(er));     % normalised: 0 = exactly stationary
  T=struct('max_normalised_gradient',max(cosang),'median_normalised_gradient',median(cosang),'iterations',o.iterations);
 case 'decomp'
  [tau,fmax]=varargin{:};
  f=(5:2:fmax)'; s=1i*2*pi*f.';
  J=jac(s,vecM,B0,[a1(:),a2(:)]); Jr=[real(J);imag(J)]; sc=vecnorm(Jr); [Q,R]=qr(Jr./sc,0);
  comps={vecM./s.^2}; names={'rigid'};
  for k=1:K, comps{end+1}=Rk{k}./(s.^2+2*zt(k)*2*pi*ft(k)*s+(2*pi*ft(k))^2); names{end+1}=sprintf('mode%d',k); end %#ok<AGROW>
  D=zeros(numel(comps),K); Fq=D;
  for c=1:numel(comps)
    t=-tau*(s.*comps{c}); t=t(:); d=(R\(Q.'*[real(t);imag(t)]))./sc.';
    da1=d(noni+(0:K-1)*(2+noni)+1).'; da2=d(noni+(0:K-1)*(2+noni)+2).';
    D(c,:)=100*(da1./a1-da2./(2*a2)); Fq(c,:)=-100*da2./(2*a2);
  end
  T=struct('names',{names},'dzeta_pct',D,'dfreq_pct',Fq,'total_dzeta',sum(D,1),'total_dfreq',sum(Fq,1));
 case 'info_fmin'
  [fmins,numdeg,tau,snr_db]=varargin{:};
  I=zeros(size(fmins));
  for q=1:numel(fmins)
    Tq=theory_check('info',[fmins(q) 800],numdeg,tau,snr_db); I(q)=Tq.Itt;
  end
  T=struct('fmin',fmins,'Itt',I);
end
  function G=truth(s)
    G=vecM./s.^2; for kk=1:K, G=G+Rk{kk}./(s.^2+2*zt(kk)*2*pi*ft(kk)*s+(2*pi*ft(kk))^2); end
  end
  function J=jac(s,M,Bm,a)
    N=numel(s); J=zeros(noni*N,noni+K*(2+noni));
    for n=1:N, J((n-1)*noni+(1:noni),1:noni)=eye(noni)/s(n)^2; end
    for kk=1:K, ik=noni+(kk-1)*(2+noni)+(1:2+noni);
      for n=1:N, A=1+a(kk,1)*s(n)+a(kk,2)*s(n)^2; P=Bm(:,kk)/A;
        J((n-1)*noni+(1:noni),ik)=[-s(n)*P/A,-s(n)^2*P/A,eye(noni)/A]; end
    end
  end
end
