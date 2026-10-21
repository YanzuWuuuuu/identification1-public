function summary=clv_run_batch(solver,structured)
if nargin<1,solver='constrained-lsq';end
if nargin<2,structured=false;end
here=fileparts(mfilename('fullpath'));sourceDir=fullfile(here,'outputs');outdir=sourceDir;addpath(here);
if structured,outdir=fullfile(here,'outputs_structured');end
if ~isfolder(outdir),mkdir(outdir);end
diary(fullfile(outdir,'batch.log'));cleanDiary=onCleanup(@()diary('off'));
S=load(fullfile(sourceDir,'data_distinct.mat'));cfg=S.cfg;cfg.solver=solver;cfg.additive_first_coefficient=structured;
manifest=struct('matlab_version',version,'solver',solver,'config',cfg, ...
    'run_command',sprintf('clv_run_batch(''%s'',%d)',solver,structured), ...
    'scope','Five-mode D=0 digital closed-loop validation. No truth enters the estimator.', ...
    'selection','Minimum training cost among stable, feasible, converged starts. If none, a best feasible diagnostic candidate is explicitly labelled.', ...
    'pilot_distinct_seed',2026101401,'training_distinct_seed',2026101402, ...
    'pilot_common_seed',2026101404,'training_common_seed',2026101405, ...
    'validation_seed',2026101406,'replicate_seeds',[2026101411 2026101412]);
writejson(fullfile(outdir,'manifest.json'),manifest);
codeDir=fullfile(outdir,'code_snapshot');if ~isfolder(codeDir),mkdir(codeDir);end
copyfile(fullfile(here,'*.m'),codeDir);
pool=gcp('nocreate');createdPool=isempty(pool);if createdPool,pool=parpool('Processes',4);end
cleanupPool=onCleanup(@()closepool(pool,createdPool)); %#ok<NASGU>
testCfg=cfg;testCfg.frequencies_hz=cfg.validation_frequencies_hz;testCfg.save_raw=true;
testDistinct=clv_data(S.sys,S.K,testCfg,8,2026101406,cfg.sensor_noise_std);
save(fullfile(outdir,'test_distinct.mat'),'testDistinct','testCfg','-v7.3');
common=S.p;common.filter_y=20e-6*ones(4,1);common.filter_u=20e-6*ones(4,1);
common.a=40e-6*ones(4);common.b=(20e-6)^2*ones(4);
[sysC,KC,infoC,plantC]=clv_system(common,S.nom,cfg);
pilotC=clv_data(sysC,KC,cfg,32,2026101404,cfg.sensor_noise_std);
[initialC,initInfoC]=clv_initial(pilotC,KC,cfg);
dataC=clv_data(sysC,KC,cfg,32,2026101405,cfg.sensor_noise_std);
testCommon=clv_data(sysC,KC,testCfg,8,2026101407,cfg.sensor_noise_std);
save(fullfile(outdir,'data_common.mat'),'common','sysC','KC','infoC','plantC','pilotC','initialC','initInfoC','dataC','testCommon','cfg','-v7.3');
cases={'distinct_N32','distinct_N8','common_N32','distinct_repeat2_N32','distinct_repeat3_N32'};
summary=struct([]);
for ic=1:numel(cases)
    label=cases{ic};truth=S.p;nominal=S.nom;pilot=S.pilot;initial=S.initial;K=S.K;data=S.data;test=testDistinct;starts=16;
    if ic==2,data=truncate(data,8);
    elseif ic==3,truth=common;pilot=pilotC;initial=initialC;K=KC;data=dataC;test=testCommon;
    elseif ic>=4
        data=clv_data(S.sys,S.K,cfg,32,2026101411+(ic-4),cfg.sensor_noise_std);starts=4;
    end
    caseDir=fullfile(outdir,label);if ~isfolder(caseDir),mkdir(caseDir);end
    save(fullfile(caseDir,'data.mat'),'data','pilot','initial','cfg','-v7.3');
    fprintf('BEGIN %s with %d starts\n',label,starts);
    parfor r=1:starts
        path=fullfile(caseDir,sprintf('start_%02d.mat',r));
        if ~isfile(path)
            try
                fit=clv_fit(data,pilot,K,cfg,initial,r,cfg.max_iterations);
                try
                    [~,~,stability]=clv_system(fit.estimate,nominal,cfg);
                    fit.closed_stable=true;fit.closed_spectral_radius=stability.closed_spectral_radius;fit.stability_error='';
                catch e
                    fit.closed_stable=false;fit.closed_spectral_radius=NaN;fit.stability_error=e.message;
                end
                result=struct('fit',fit);saveone(path,result);
                compact=rmfield(fit,{'estimate','x','x0','chart','sigma'});
            catch e
                result=struct('failure',getReport(e,'extended','hyperlinks','off'));saveone(path,result);
                compact=struct('start_index',r,'failure',result.failure);
            end
            writejson(fullfile(caseDir,sprintf('start_%02d.json',r)),compact);
        end
    end
    rows=struct([]);best=Inf;bestAny=Inf;selected=[];diagnostic=[];successful=0;feasible=0;failures=0;
    for r=1:starts
        res=load(fullfile(caseDir,sprintf('start_%02d.mat',r)));
        if ~isfield(res,'fit'),failures=failures+1;continue;end
        fit=res.fit;compact=rmfield(fit,{'estimate','x','x0','chart','sigma'});
        if isempty(rows),rows=compact;else,rows(end+1)=compact;end %#ok<AGROW>
        valid=fit.constraint_violation<=1e-8&&fit.closed_stable;
        feasible=feasible+valid;successful=successful+(valid&&fit.solver_pass);
        if valid&&fit.cost<bestAny,bestAny=fit.cost;diagnostic=fit;end
        if valid&&fit.solver_pass&&fit.cost<best,best=fit.cost;selected=fit;end
    end
    accepted=~isempty(selected);if ~accepted,selected=diagnostic;end
    entry=struct('case',label,'groups',data.groups,'starts',starts,'strict_converged_starts',successful, ...
        'stable_feasible_starts',feasible,'exception_starts',failures,'accepted_solver_result',accepted, ...
        'time_command_identity_relative_error',data.time_command_identity_relative_error, ...
        'input_noise_correlation',data.input_sensor_noise_complex_correlation, ...
        'selected_start',NaN,'cost',NaN,'selected_first_order_optimality',NaN,'metrics',struct(),'failure','');
    if ~isempty(selected)
        metrics=clv_metrics(truth,nominal,selected,K,cfg,test);
        entry.selected_start=selected.start_index;entry.cost=selected.cost;
        entry.selected_first_order_optimality=selected.first_order_optimality;
        entry.metrics=metrics;
        save(fullfile(caseDir,'selected.mat'),'selected','metrics','entry','rows','-v7.3');
    else,entry.failure='No stable feasible fit; no calibration claim.';end
    writejson(fullfile(caseDir,'summary.json'),entry);
    if isempty(summary),summary=entry;else,summary(end+1)=entry;end %#ok<AGROW>
    writejson(fullfile(outdir,'summary.json'),summary);
    fprintf('END %s: strict %d/%d, selected accepted=%d\n',label,successful,starts,accepted);
end
save(fullfile(outdir,'summary.mat'),'summary','manifest');
end
function d=truncate(d,N)
d.Y=d.Y(:,:,1:N);d.U=d.U(:,:,1:N);d.V=d.V(:,:,1:N);d.groups=N;d.mean=mean(d.Y,3);
end
function saveone(path,value)
save(path,'-struct','value');
end
function writejson(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
function closepool(pool,created)
if created&&isvalid(pool),delete(pool);end
end
