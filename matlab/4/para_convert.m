% Stage 4: three rigid directions and four rank-one flexible residues.
stage_dir=fileparts(mfilename('fullpath'));
addpath(stage_dir,fullfile(fileparts(stage_dir),'common'));
cfg=identification_config;
if ~exist('plot_enabled','var'), plot_enabled=true; end
additive=load_stage_result(cfg,3);
[rho0,initialization]=modal_initialization(additive.belta_vecest,cfg.no,cfg.ni,cfg.K,cfg.nrbm);
[num,rou_est,info]=Gauss_Newton(additive.belta_vecest,rho0,0.01,10,200, ...
    cfg.tolerance,cfg.no,cfg.ni,cfg.K,cfg.nrbm,additive.information_metric);
belta_modal=f_rho_construct(rou_est,cfg.no,cfg.ni,cfg.K,cfg.nrbm);
validation_est=additive_frf(belta_modal,additive.validation_f,cfg.no*cfg.ni,cfg.K);
metrics=frf_metrics(additive.validation_f,additive.validation_G,validation_est,cfg.weight_floor_relative);
stage_result=struct('cfg',cfg,'rho0',rho0,'rou_est',rou_est,'belta_modal',belta_modal, ...
    'initialization',initialization,'info',info,'metrics',metrics, ...
    'validation_f',additive.validation_f,'validation_G',additive.validation_G, ...
    'validation_est',validation_est,'metric_type',additive.metric_type);
save_stage_result(cfg,4,stage_result);
fprintf('Stage 4: %d rigid directions, %d flexible modes, %d iterations, %s.\n',cfg.nrbm,cfg.K,num,info.reason);
if plot_enabled
    fig=plot_frf_comparison(additive.validation_f,additive.validation_G,validation_est, ...
        cfg.no,cfg.ni,'Stage 4: four-mode model, validation includes omitted fifth mode');
    exportgraphics(fig,fullfile(cfg.result_dir,'modal_validation.png'),'Resolution',160);
end
