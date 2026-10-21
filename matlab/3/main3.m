% Stage 3: four-mode additive RIV fit; fifth plant mode remains unmodeled.
stage_dir=fileparts(mfilename('fullpath'));
addpath(stage_dir,fullfile(fileparts(stage_dir),'common'));
cfg=identification_config;
if ~exist('plot_enabled','var'), plot_enabled=true; end
initial=load_stage_result(cfg,1); rigid=load_stage_result(cfg,2);
noni=cfg.no*cfg.ni; K=cfg.K;
[frf_select,f_select]=select_frequency_data(rigid.flexible_frf,cfg.flexible_frequency_hz,rigid.f);
belta0=zeros((noni+2)*K,K);
for k=1:K
    ix=(k-1)*(noni+2)+(1:noni+2);
    a=initial.den_para_init{k};
    belta0(ix,k)=[a(2);a(1);initial.belta_init(noni+(k-1)*noni+(1:noni))];
end
belta0(1:noni+2,1)=rigid.first_mode;
weight=frf_weights(frf_select,cfg.weight_floor_relative);
[num,belta_est,info]=SRIVC2(f_select,frf_select,weight,belta0,cfg.tolerance,cfg.max_iterations,K,noni,cfg.frequency_bounds_hz);
belta_riv=belta_est;
packed=zeros(K*(noni+2),1);
for k=1:K
    ix=(k-1)*(noni+2)+(1:noni+2); packed(ix)=belta_est(ix,k);
end
[packed,refinement]=refine_flexible_fit(packed,f_select,frf_select,weight,K, ...
    cfg.frequency_bounds_hz,cfg.max_iterations,cfg.tolerance);
for k=1:K
    ix=(k-1)*(noni+2)+(1:noni+2); belta_est(ix,k)=packed(ix);
end
belta_vecest=zeros(noni+K*(noni+2),1); belta_vecest(1:noni)=rigid.vecM_est;
for k=1:K
    ix=(k-1)*(noni+2)+(1:noni+2);
    belta_vecest(noni+ix)=belta_est(ix,k);
end
% Evaluate the ESTIMATED parameters; no reference-model variables enter here.
frf_est=additive_frf(belta_vecest,f_select,noni,K)-rigid.vecM_est./((1i*2*pi*f_select.').^2);
[validation_G,validation_f]=select_frequency_data(rigid.G,cfg.validation_frequency_hz,rigid.f);
validation_est=additive_frf(belta_vecest,validation_f,noni,K);
% Local curvature of the measured-FRF LS objective. This is an information
% metric, not an estimated statistical covariance from repeated experiments.
metric_f=unique([cfg.rigid_frequency_hz;cfg.flexible_frequency_hz]);
[metric_G,metric_f]=select_frequency_data(rigid.G,metric_f,rigid.f);
[~,J]=additive_frf(belta_vecest,metric_f,noni,K);
Jw=apply_sqrt_weight(J,frf_weights(metric_G,cfg.weight_floor_relative));
scales=max(vecnorm(Jw),realmin); [~,R]=qr(Jw./scales,0);
information_metric=R.*scales;
metrics=frf_metrics(validation_f,validation_G,validation_est,cfg.weight_floor_relative);
stage_result=struct('cfg',cfg,'belta_vecest',belta_vecest,'belta_est',belta_est, ...
    'isubsystem',K,'f_select',f_select,'frf_select',frf_select,'frf_est',frf_est, ...
    'validation_f',validation_f,'validation_G',validation_G,'validation_est',validation_est, ...
    'information_metric',information_metric,'metric_type','local_FRF_LS_information', ...
    'info',info,'refinement',refinement,'belta_riv',belta_riv,'metrics',metrics);
save_stage_result(cfg,3,stage_result);
fprintf('Stage 3: %d flexible modes, %d iterations, %s.\n',K,num,info.reason);
fprintf('  Weighted LS refinement: %d iterations, exit flag %d.\n',refinement.iterations,refinement.exitflag);
if plot_enabled, plot_frf_comparison(f_select,frf_select,frf_est,cfg.no,cfg.ni,'Stage 3: four identified flexible modes'); end
