function r = delay_step(action, varargin)
% Small resumable steps for the delay feasibility study (each call < 60 s).
% delay_step('setup', tau_true, snr_db, seed)  -> stores data in delay_state.mat
% delay_step('fit', tau)                      -> returns [cost, f(5), zeta(5)] for compensation tau
here=fileparts(mfilename('fullpath')); st=fullfile(here,'delay_state.mat');
switch action
  case 'setup'
    [tau_true,snr_db,seed]=varargin{:};
    root=fileparts(here); old=pwd; cd(fullfile(root,'3')); no=4; ni=4; model; cd(old); %#ok<NASGU>
    ft=[286.43,344.58,488.27,705.67,738.83]; zt=0.01*ones(1,5);
    f=(5:2:800)'; s=1i*2*pi*f.'; G0=vecM./s.^2;
    for k=1:5, R=C2(:,k)*B2(k,:); G0=G0+R(:)./(s.^2+2*zt(k)*2*pi*ft(k)*s+(2*pi*ft(k))^2); end
    G0=G0.*exp(-s*tau_true); sigma=10^(-snr_db/20)*median(abs(G0(:,f>=200)),'all');
    rng(seed); G=G0+sigma*(randn(size(G0))+1i*randn(size(G0)))/sqrt(2);
    save(st,'f','s','G','ft','zt','tau_true','snr_db','seed','vecM'); r=sigma;
  case 'fit'
    tau=varargin{1}; S=load(st);
    o=joint_fullband_srivc(S.f,S.G.*exp(S.s*tau),[280,355,480,690,765],0.02*ones(1,5), ...
        struct('colscale',true,'weight','none','maxit',80));
    J=norm(S.G.*exp(S.s*tau)-o.pred,'fro')^2;
    r=[J,o.f_hz(:)',o.zeta(:)',100*norm(o.M-S.vecM)/norm(S.vecM)];
end
end
