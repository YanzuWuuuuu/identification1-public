function out = delay_feasibility(nMC, tau_true, snr_flex_db)
% Feasibility: unmodeled vs estimated loop delay in joint additive SRIVC (frequency domain).
% Delay is compensated EXACTLY in the frequency domain: G(w)*exp(+j*w*tau).
% Outer: 1-D golden-section search on tau; inner: joint additive SRIVC (rigid + 5 modes).
if nargin<1, nMC=3; end, if nargin<2, tau_true=100e-6; end, if nargin<3, snr_flex_db=10; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
old=pwd; cd(fullfile(root,'3')); no=4; ni=4; model; cd(old); %#ok<NASGU>
ft=[286.43,344.58,488.27,705.67,738.83]; zt=0.01*ones(1,5); K=5;
fpri=[280,355,480,690,765]; zpri=0.02*ones(1,5);
f=(5:2:800)'; s=1i*2*pi*f.'; G0=vecM./s.^2;
for k=1:K, R=C2(:,k)*B2(k,:); G0=G0+R(:)./(s.^2+2*zt(k)*2*pi*ft(k)*s+(2*pi*ft(k))^2); end
G0=G0.*exp(-s*tau_true);
sigma=10^(-snr_flex_db/20)*median(abs(G0(:,f>=200)),'all');
opts=struct('colscale',true,'weight','none','maxit',80);
fit=@(G,tau) joint_fullband_srivc(f,G.*exp(s*tau),fpri,zpri,opts);
cost=@(o,G,tau) norm(G.*exp(s*tau)-o.pred,'fro')^2;
res=zeros(nMC,2*(2*K)+2);
for m=1:nMC
  rng(2000+m); G=G0+sigma*(randn(size(G0))+1i*randn(size(G0)))/sqrt(2);
  o0=fit(G,0);                                    % delay ignored (current practice)
  % golden-section search on tau in [0, 400] us
  a=0; b=400e-6; g=(sqrt(5)-1)/2; c=b-g*(b-a); d=a+g*(b-a);
  fc=cost(fit(G,c),G,c); fd=cost(fit(G,d),G,d);
  for it=1:14
    if fc<fd, b=d; d=c; fd=fc; c=b-g*(b-a); fc=cost(fit(G,c),G,c);
    else, a=c; c=d; fc=fd; d=a+g*(b-a); fd=cost(fit(G,d),G,d); end
  end
  tau_hat=(a+b)/2; o1=fit(G,tau_hat);
  res(m,:)=[o0.f_hz(:)',o0.zeta(:)',o1.f_hz(:)',o1.zeta(:)',tau_hat,0];
end
fb=@(F) 100*(mean(F)-ft)./ft; zb=@(Z) 100*(mean(Z)-zt)./zt; zsd=@(Z) 100*std(Z)./zt;
fprintf('true tau = %g us, SNR(flex) = %d dB, %d runs\n',tau_true*1e6,snr_flex_db,nMC);
fprintf('delay ignored : freq bias%% %s | zeta bias%% %s\n',mat2str(fb(res(:,1:5)),2),mat2str(zb(res(:,6:10)),2));
fprintf('delay estimated: freq bias%% %s | zeta bias%% %s | zeta std%% %s\n',mat2str(fb(res(:,11:15)),2),mat2str(zb(res(:,16:20)),2),mat2str(zsd(res(:,16:20)),2));
fprintf('tau_hat = %s us\n',mat2str(1e6*res(:,21)',4));
% cost profile for the last realization
taus=(0:25:300)*1e-6; J=arrayfun(@(t) cost(fit(G,t),G,t),taus);
fprintf('cost profile J(tau)/min: %s  (tau = 0:25:300 us)\n',mat2str(J/min(J),3));
out=struct('res',res,'taus',taus,'J',J);
save(fullfile(here,sprintf('delay_feas_tau%gus_snr%d.mat',tau_true*1e6,snr_flex_db)),'out');
end
