function T = mc_realistic_noise(method, nMC, snr_flex_db, band, tau, fpri)
% Monte Carlo with ADDITIVE noise of constant level (as from output noise):
% low-frequency (rigid) SNR is high, flexible-band SNR = snr_flex_db.
% method: 'joint_uniform' | 'joint_relative' | 'split_uniform' | 'split_relative'
if nargin<3, snr_flex_db=20; end
if nargin<4, band=[5 800]; end
if nargin<5, tau=0; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
old=pwd; cd(fullfile(root,'3')); no=4; ni=4; model; cd(old); %#ok<NASGU>
ftrue=[286.43,344.58,488.27,705.67,738.83]; ztrue=0.01*ones(1,5); K=5;
if nargin<6 || isempty(fpri), fpri=[280,355,480,690,765]; end
zpri=0.02*ones(1,5);
f=(band(1):2:band(2))'; s=1i*2*pi*f.'; Mt=reshape(vecM,4,4); M0=vecM;
G0=zeros(16,numel(f));
for k=1:K, R=C2(:,k)*B2(k,:); G0=G0+R(:)./(s.^2+2*ztrue(k)*2*pi*ftrue(k)*s+(2*pi*ftrue(k))^2); end
Gflex0=G0; G0=G0+M0./s.^2;
G0=G0.*exp(-s*tau); Gflex0=Gflex0.*exp(-s*tau);
fb=f>=200; sigma=10^(-snr_flex_db/20)*median(abs(G0(:,fb)),'all');
fprintf('%s: band %g-%g Hz, sigma=%.2e, SNR@5Hz=%.0f dB, median SNR 200-800=%.0f dB\n',method,band, ...
    sigma,20*log10(median(abs(G0(:,1)))/sigma),20*log10(median(abs(G0(:,fb)),'all')/sigma));
res=nan(nMC,2*K+6);
for m=1:nMC
  rng(1000+m);
  G=G0+sigma*(randn(size(G0))+1i*randn(size(G0)))/sqrt(2);
  try
    switch method
      case {'joint_uniform','joint_relative','joint_model_relative'}
        w='none'; if contains(method,'relative'), w='relative'; end; if contains(method,'model'), w='model_relative'; end
        o=joint_fullband_srivc(f,G,fpri,zpri,struct('colscale',true,'weight',w,'maxit',80));
      case {'split_uniform','split_relative','split_model_relative'}
        % stage A: rigid band 5-150 Hz, G ~ M/s^2 + D (static compliance of flexible modes)
        ia=f<=150; Mh=zeros(16,1);
        for r=1:16
          X=[1./s(ia).^2; ones(1,nnz(ia))].'; y=G(r,ia).';
          Xr=[real(X);imag(X)]; sc=vecnorm(Xr); th=((Xr./sc)\[real(y);imag(y)])./sc.'; Mh(r)=th(1);
        end
        w='none'; if contains(method,'relative'), w='relative'; end; if contains(method,'model'), w='model_relative'; end
        ib=f>=200;
        o=joint_fullband_srivc(f(ib),G(:,ib),fpri,zpri,struct('colscale',true,'weight',w,'maxit',80,'Mfixed',Mh));
    end
    A=1+o.a(:,1)*s+o.a(:,2)*s.^2; P=o.M./s.^2; for k=1:K, P=P+o.B(:,k)./A(k,:); end
    eflex=norm(G0(:,fb)-P(:,fb),'fro')/norm(Gflex0(:,fb),'fro');
    res(m,:)=[o.f_hz(:)',o.zeta(:)',100*norm(o.M-M0)/norm(M0),100*eflex,max(o.cond_raw),max(o.cond_moment),o.iterations,1];
  catch ME
    fprintf('run %d failed: %s\n',m,ME.message); res(m,end)=0;
  end
end
ok=res(:,end)==1; R=res(ok,:);
fb_=100*(mean(R(:,1:K))-ftrue)./ftrue; fs_=100*std(R(:,1:K))./ftrue;
zb_=100*(mean(R(:,K+1:2*K))-ztrue)./ztrue; zs_=100*std(R(:,K+1:2*K))./ztrue;
fprintf('  ok %d/%d | freq bias%% %s | freq std%% %s\n',nnz(ok),nMC,mat2str(fb_,2),mat2str(fs_,2));
fprintf('  zeta bias%% %s | zeta std%% %s\n',mat2str(zb_,2),mat2str(zs_,2));
fprintf('  M err%% mean %.3f | flex-band FRF err%% mean %.2f | cond raw %.1e | cond scaled %.1e | iters %.1f\n', ...
    mean(R(:,2*K+1)),mean(R(:,2*K+2)),median(R(:,2*K+3)),median(R(:,2*K+4)),mean(R(:,2*K+5)));
T=struct('method',method,'res',res,'sigma',sigma,'band',band,'snr',snr_flex_db);
T.tau=tau; T.fpri=fpri;
save(fullfile(here,sprintf('mc_%s_snr%d_%d-%d_tau%gus_p%d.mat',method,snr_flex_db,band,tau*1e6,round(fpri(1)))),'T');
end
