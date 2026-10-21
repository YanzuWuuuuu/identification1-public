function result=run_identification(makePlots,rebuildFrf)
% Run stages 1-4; optionally rebuild stage-0 FRF from saved time signals.
% Usage: result=run_identification(true,true) to rebuild FRF and show plots.
% This does not rerun the Simulink simulation.
if nargin<1, makePlots=true; end
if nargin<2, rebuildFrf=false; end
root=fileparts(mfilename('fullpath'));
previousPath=path; previousFolder=pwd;
cleanup=onCleanup(@()restore_environment(previousPath,previousFolder));
addpath(fullfile(root,'common'));
if rebuildFrf
    addpath(fullfile(root,'0')); estimate_frf_from_saved_data;
end
entries={'1/mian1.m','2/main2.m','3/main3.m','4/para_convert.m'};
result=struct;
for k=1:numel(entries)
    stage=run_stage(fullfile(root,entries{k}),makePlots);
    result.(sprintf('stage%d',k))=stage;
end
cfg=result.stage4.cfg; rows=zeros(cfg.K,5);
for k=1:cfg.K
    ix=cfg.no*cfg.ni+(k-1)*(cfg.no*cfg.ni+2)+(1:2);
    a=result.stage3.belta_vecest(ix); b=result.stage4.belta_modal(ix);
    rows(k,:)=[k,1/(2*pi*sqrt(a(2))),a(1)/(2*sqrt(a(2))), ...
        1/(2*pi*sqrt(b(2))),b(1)/(2*sqrt(b(2)))];
end
result.modes=array2table(rows,'VariableNames',{'Mode','AdditiveFrequencyHz', ...
    'AdditiveDamping','ModalFrequencyHz','ModalDamping'});
writetable(result.modes,fullfile(cfg.result_dir,'identified_modes.csv'));
disp(result.modes);
fprintf('Modal weighted FRF error, 200-720 Hz: %.3f%%; omitted-mode band >720-800 Hz: %.3f%%.\n', ...
    100*result.stage4.metrics.flexible_band.relative_weighted_error, ...
    100*result.stage4.metrics.unmodeled_band.relative_weighted_error);
if any(result.stage3.refinement.frequency_bound_active)
    fprintf('Note: additive mode(s) %s reached their frequency bounds; inspect truncation bias.\n', ...
        mat2str(find(result.stage3.refinement.frequency_bound_active).'));
end
save(fullfile(cfg.result_dir,'identification_result.mat'),'result');
end

function stage=run_stage(filename,plot_enabled) %#ok<INUSD>
run(filename);
stage=stage_result;
end

function restore_environment(previousPath,previousFolder)
path(previousPath); cd(previousFolder);
end
