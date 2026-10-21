% Does the rigid/flexible magnitude gap block JOINT full-band (5-800 Hz) identification?
% Compares: no column scaling vs column scaling; relative weighting vs none;
% noise-free exact FRF, noisy FRF, and the existing LRM FRF.
here=fileparts(mfilename('fullpath')); root=fileparts(here); addpath(here);
cd(fullfile(root,'3')); no=4; ni=4; model; cd(here);
ftrue=[286.43,344.58,488.27,705.67,738.83]; ztrue=0.01*ones(1,5);
f_ode3=[286.573 344.945 490.39 718.958 755.501]; z_ode3=[0.0119 0.0132 0.0186 0.0322 0.0346];
fpri=[286,347.5,487.5,700,740]; zpri=0.02*ones(1,5);
f=(5:2:800)'; Mt=reshape(vecM,4,4);
Gex=zeros(16,numel(f));
for n=1:numel(f)
  sn=1i*2*pi*f(n); Gn=Mt/sn^2;
  for k=1:5, Gn=Gn+C2(:,k)*B2(k,:)/(sn^2+2*ztrue(k)*2*pi*ftrue(k)*sn+(2*pi*ftrue(k))^2); end
  Gex(:,n)=Gn(:);
end
rng(1); Gno=Gex.*(1+0.01*(randn(size(Gex))+1i*randn(size(Gex)))/sqrt(2));
addpath(fullfile(root,'common')); cfg=identification_config; [Gl,fl]=read_lrm_data(cfg.data_file,4,4); fl=fl(:);
[~,ix]=min(abs(fl-f.'),[],1); Glrm=Gl(:,ix);
cases={'exact',Gex,ftrue,ztrue;'noisy1pct',Gno,ftrue,ztrue;'LRM(vs ode3 truth)',Glrm,f_ode3,z_ode3};
cfgs={struct('colscale',false,'weight','none'),struct('colscale',false,'weight','relative'), ...
      struct('colscale',true,'weight','none'),struct('colscale',true,'weight','relative')};
names={'noscale/noW','noscale/relW','colscale/noW','colscale/relW'};
if ~exist('caseList','var'), caseList=1:size(cases,1); end
for c=caseList
  fprintf('\n=== data: %s (joint 5-800 Hz, rigid + 5 modes) ===\n',cases{c,1});
  for q=1:numel(cfgs)
    try
      o=joint_fullband_srivc(f,cases{c,2},fpri,zpri,cfgs{q});
      fe=100*abs(o.f_hz(:)'-cases{c,3})./cases{c,3}; ze=100*abs(o.zeta(:)'-cases{c,4})./cases{c,4};
      Me=100*norm(o.M-vecM)/norm(vecM);
      fprintf('%-13s it=%3d cond=%.1e | freq err%% %s | zeta err%% %s | M err %.3f%% | fit %.2e (w %.2e)\n', ...
        names{q},o.iterations,max(o.cond_moment),mat2str(fe,2),mat2str(ze,2),Me,o.rel_err,o.rel_err_weighted);
    catch ME
      fprintf('%-13s FAILED: %s\n',names{q},ME.message);
    end
  end
end
