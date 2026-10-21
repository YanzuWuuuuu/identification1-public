function cfg=identification_config
% Four flexible modes are identified; plant mode 5 remains unmodeled.
cfg.root=fileparts(fileparts(mfilename('fullpath')));
cfg.no=4; cfg.ni=4; cfg.nrbm=3; cfg.K=4;
cfg.frequency_prior_hz=[286,347.5,487.5,700];
cfg.damping_prior=0.02;
cfg.initial_frequency_hz=(40:2:720).';
cfg.rigid_frequency_hz=(5:1:310).';
cfg.flexible_frequency_hz=(200:2:720).';
cfg.validation_frequency_hz=(5:1:800).';
cfg.frequency_bounds_hz=[cfg.frequency_prior_hz(:)*0.9, ...
    min(cfg.frequency_prior_hz(:)*1.1,720)];
cfg.weight_floor_relative=1e-3;
cfg.tolerance=1e-8; cfg.max_iterations=100;
cfg.result_dir=fullfile(cfg.root,'results');
% Stage 0 can export a new dataset here. Until then use the supplied FRF.
cfg.data_file=fullfile(cfg.root,'0','LRM_data.mat');
if ~isfile(cfg.data_file), cfg.data_file=fullfile(cfg.root,'1','LRM_data.mat'); end
source=dir(cfg.data_file);
assert(~isempty(source),'LRM data file is missing.');
cfg.source_bytes=source.bytes; cfg.source_modified=source.datenum;
end
