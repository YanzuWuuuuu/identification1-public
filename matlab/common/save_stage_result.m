function save_stage_result(cfg,number,stage)
if ~isfolder(cfg.result_dir), mkdir(cfg.result_dir); end
save(fullfile(cfg.result_dir,sprintf('stage%d.mat',number)),'-struct','stage');
end
