function r = study2_step(action, varargin)
% Stepwise experiments (each call < 60 s).
% 'data'  : (sysname, tau, snr_db, seed)  -> builds G and saves study2_state.mat
%           sysname 'nominal' (5 modes) | 'repeated' (adds a 6th mode 0.05% from mode 4, zeta .002)
% 'fit'   : (tau_comp, opts) -> fits on G.*exp(s*tau_comp), returns struct
here=fileparts(mfilename('fullpath')); st=fullfile(here,'study2_state.mat');
switch action
 case 'data'
  [sysname,tau,snr_db,seed]=varargin{:};
  root=fileparts(here); old=pwd; cd(fullfile(root,'3')); no=4; ni=4; model; cd(old); %#ok<NASGU>
  ft=[286.43,344.58,488.27,705.67,738.83]; zt=0.01*ones(1,5); Rk=arrayfun(@(k) reshape(C2(:,k)*B2(k,:),[],1),1:5,'uni',0);
  if strcmp(sysname,'no5'), ft=ft(1:4); zt=zt(1:4); Rk=Rk(1:4); end
  if strcmp(sysname,'repeated')
    ft=[ft, 705.67*1.0005]; zt=[zt 0.01]; zt(4)=0.002; zt(6)=0.002; rng(99); Rk{6}=reshape(randn(4,1)*randn(1,4)*0.2,[],1);
  end
  f=(5:2:800)'; s=1i*2*pi*f.'; G0=vecM./s.^2;
  for k=1:numel(ft), G0=G0+Rk{k}./(s.^2+2*zt(k)*2*pi*ft(k)*s+(2*pi*ft(k))^2); end
  G0=G0.*exp(-s*tau); sigma=10^(-snr_db/20)*median(abs(G0(:,f>=200)),'all');
  rng(seed); G=G0+sigma*(randn(size(G0))+1i*randn(size(G0)))/sqrt(2);
  if isinf(snr_db), G=G0; end
  save(st,'f','s','G','G0','ft','zt','tau','vecM','sysname'); r=sigma;
 case 'fit'
  [tc,opts]=varargin{:}; S=load(st);
  fpri=[280,355,480,690,765]; if strcmp(S.sysname,'repeated'), fpri=[280,355,480,700,765,712]; end
  if strcmp(S.sysname,'no5'), fpri=[280,355,480,690]; end
  o=joint_srivc2(S.f,S.G.*exp(S.s*tc),fpri,0.02*ones(1,numel(fpri)),opts);
  r=struct('f',o.f_hz','z',o.zeta','Merr',100*norm(o.M-S.vecM)/norm(S.vecM),'cost',o.cost, ...
    'cond',max(o.cond),'it',o.iterations,'ft',S.ft,'zt',S.zt);
  fid=fopen(fullfile(here,'study2_log.txt'),'a');
  fprintf(fid,'%s | sys=%s tau=%g | tc=%.2fus solver=%s scale=%d numdeg=%d | it=%d cond=%.2e cost=%.6e M=%.4f%% | f=%s | z=%s\n', ...
    datestr(now,'HH:MM:SS'),S.sysname,S.tau,tc*1e6,opts.solver,opts.colscale,getfield_def(opts,'numdeg',0),o.iterations,max(o.cond),o.cost,r.Merr,mat2str(r.f,7),mat2str(r.z,4));
  fclose(fid);
end
end

function v=getfield_def(o,f,d)
if isfield(o,f), v=o.(f); else, v=d; end
end
