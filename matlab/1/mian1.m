% Stage 1: fixed-pole, weighted linear least-squares initialization.
stage_dir=fileparts(mfilename('fullpath'));
addpath(stage_dir,fullfile(fileparts(stage_dir),'common'));
cfg=identification_config;
if ~exist('plot_enabled','var'), plot_enabled=true; end
[G,f]=read_lrm_data(cfg.data_file,cfg.no,cfg.ni);
[frf_select,f_select]=select_frequency_data(G,cfg.initial_frequency_hz,f);
noni=cfg.no*cfg.ni; isubsystem=cfg.K;
den_para_init=cell(1,isubsystem);
for k=1:isubsystem
    omega=2*pi*cfg.frequency_prior_hz(k);
    den_para_init{k}=[1/omega^2,2*cfg.damping_prior/omega];
end
weight=frf_weights(frf_select,cfg.weight_floor_relative);
[num,belta_init]=SRIVC_init(f_select,frf_select,den_para_init,weight,[],noni,[],[],isubsystem);
beta=zeros(noni+isubsystem*(noni+2),1); beta(1:noni)=belta_init(1:noni);
for k=1:isubsystem
    a=den_para_init{k}; ix=noni+(k-1)*(noni+2)+(1:noni+2);
    beta(ix)=[a(2);a(1);belta_init(noni+(k-1)*noni+(1:noni))];
end
frf_est=additive_frf(beta,f_select,noni,isubsystem);
stage_result=struct('cfg',cfg,'belta_init',belta_init,'den_para_init',{den_para_init}, ...
    'f_select',f_select,'frf_select',frf_select,'frf_est',frf_est,'isubsystem',isubsystem);
save_stage_result(cfg,1,stage_result);
fprintf('Stage 1: %d channels, %d flexible modes, %d frequency samples.\n',noni,isubsystem,numel(f_select));
if plot_enabled, plot_frf_comparison(f_select,frf_select,frf_est,cfg.no,cfg.ni,'Stage 1: fixed-pole initialization'); end
