% Stage 2: identify the rigid contribution without using reference parameters.
stage_dir=fileparts(mfilename('fullpath'));
addpath(stage_dir,fullfile(fileparts(stage_dir),'common'));
cfg=identification_config;
if ~exist('plot_enabled','var'), plot_enabled=true; end
initial=load_stage_result(cfg,1);
[G,f]=read_lrm_data(cfg.data_file,cfg.no,cfg.ni);
[frf_select,f_select]=select_frequency_data(G,cfg.rigid_frequency_hz,f);
noni=cfg.no*cfg.ni;
a=initial.den_para_init{1}; M0=initial.belta_init(1:noni);
B0=initial.belta_init(noni+(1:noni));
belta0=[a(2);a(1);M0;a(2)*M0;a(1)*M0+B0];
weight=frf_weights(frf_select,cfg.weight_floor_relative);
[num,belta_est,info]=SRIVC1(f_select,frf_select,weight,belta0,noni,cfg.tolerance,cfg.max_iterations,1);
vecM_est=belta_est(3:noni+2);
first_mode=[belta_est(1:2);belta_est(2*noni+3:3*noni+2)-belta_est(2)*vecM_est];
valid=f>0 & all(isfinite(G),1).'; f=f(valid); G=G(:,valid);
rigid_frf=vecM_est./((1i*2*pi*f.').^2);
flexible_frf=G-rigid_frf;
s=1i*2*pi*f_select.';
frf_est=(belta_est(3:noni+2)+belta_est(noni+3:2*noni+2)*s+ ...
    belta_est(2*noni+3:3*noni+2)*s.^2)./(s.^2.*(1+belta_est(1)*s+belta_est(2)*s.^2));
stage_result=struct('cfg',cfg,'vecM_est',vecM_est,'first_mode',first_mode, ...
    'belta_est',belta_est,'f',f,'G',G,'flexible_frf',flexible_frf, ...
    'f_select',f_select,'frf_select',frf_select,'frf_est',frf_est,'info',info);
save_stage_result(cfg,2,stage_result);
fprintf('Stage 2: %d iterations, %s.\n',num,info.reason);
if plot_enabled, plot_frf_comparison(f_select,frf_select,frf_est,cfg.no,cfg.ni,'Stage 2: rigid and first flexible mode'); end
