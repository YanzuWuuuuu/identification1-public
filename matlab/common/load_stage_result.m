function stage=load_stage_result(cfg,number)
filename=fullfile(cfg.result_dir,sprintf('stage%d.mat',number));
assert(isfile(filename),'Run stage %d first: missing %s',number,filename);
stage=load(filename);
assert(isfield(stage,'cfg')&&isequal(stage.cfg,cfg), ...
    'identification:StaleStage', ...
    'Configuration or source data changed. Rerun identification from stage 1.');
end
