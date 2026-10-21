function result=run_raw_identification(makePlots,compareLrm)
% Feasibility experiment: direct raw input/output identification, four modes.
% run_raw_identification(true,false) never estimates or loads an LRM FRF.
if nargin<1, makePlots=true; end
if nargin<2, compareLrm=true; end
root=fileparts(mfilename('fullpath')); oldPath=path;
cleanup=onCleanup(@()path(oldPath));
addpath(fullfile(root,'common'),fullfile(root,'raw'));
out=fullfile(root,'results_raw'); if ~isfolder(out), mkdir(out); end
ts=.0002; ni=4; K=4; nrbm=3;
frequency0=[286,347.5,487.5,700]; damping0=repmat(.02,1,K);
bounds=[frequency0(:)*.9,min(frequency0(:)*1.1,720)];
u=[]; y=[];
for k=1:ni
    input=load(fullfile(root,'0',sprintf('u%d.mat',k)));
    output=load(fullfile(root,'0',sprintf('y%d.mat',k)));
    u(:,k)=input.(sprintf('u%d',k))(:); %#ok<AGROW>
    y(:,k)=output.(sprintf('y%d',k))(:); %#ok<AGROW>
end
assert(isequal(size(u),size(y))&&all(isfinite([u(:);y(:)])));
ntrain=floor(.7*size(u,1)); train=1:ntrain;
fprintf('Raw experiment: %d training and %d held-out samples, %.0f Hz.\n',ntrain,size(u,1)-ntrain,1/ts);
models=struct;
models.raw_time=fit_raw_modal(u(train,:),y(train,:),ts,frequency0,damping0,bounds,'foh',100);
fprintf('Direct time fit complete; applying modal rank constraints.\n');
models.raw_time_structured=fit_raw_structured(u(train,:),y(train,:),models.raw_time,nrbm,100);
models.raw_band=fit_raw_modal(u(train,:),y(train,:),ts,frequency0,damping0,bounds,'foh',100,[0,720]);
models.raw_band_structured=fit_raw_structured(u(train,:),y(train,:),models.raw_band,nrbm,100);
models.raw_zoh=fit_raw_modal(u(train,:),y(train,:),ts,frequency0,damping0,bounds,'zoh',100);
if compareLrm
    fprintf('Retraining the LRM reference on training samples only.\n');
    addpath(fullfile(root,'1'),fullfile(root,'2'),fullfile(root,'3'),fullfile(root,'4'));
    models.lrm_reference=fit_lrm_training_reference(u(train,:),y(train,:),ts);
end
names=fieldnames(models); rows=zeros(numel(names)*K,4); labels=strings(numel(names)*K,1);
errors=zeros(numel(names),6); evaluation=struct;
for j=1:numel(names)
    fit=models.(names{j}); evaluation.(names{j})=evaluate_raw_fit(u,y,fit,ntrain);
    errors(j,:)=100*evaluation.(names{j}).metrics;
    ix=(j-1)*K+(1:K); labels(ix)=string(names{j});
    rows(ix,:)=[(1:K)',fit.frequency_hz,fit.damping,abs(fit.frequency_hz-bounds(:,2))<1e-4];
end
modeTable=addvars(array2table(rows,'VariableNames',{'Mode','FrequencyHz','Damping','UpperBoundActive'}),labels,'Before',1,'NewVariableNames','Method');
metricTable=addvars(array2table(errors,'VariableNames',{'TrainTimePercent','HoldoutTimePercent', ...
    'Holdout5to200Percent','Holdout200to720Percent','HoldoutOver720to800Percent','Holdout200to800Percent'}), ...
    string(names),'Before',1,'NewVariableNames','Method');
result=struct('models',models,'evaluation',evaluation,'modes',modeTable,'metrics',metricTable, ...
    'sample_time',ts,'training_samples',ntrain,'sample_count',size(u,1), ...
    'frequency_prior_hz',frequency0,'frequency_bounds_hz',bounds, ...
    'nrbm',nrbm,'flexible_modes',K,'input',u,'output',y);
writetable(modeTable,fullfile(out,'mode_comparison.csv'));
writetable(metricTable,fullfile(out,'holdout_metrics.csv'));
save(fullfile(out,'raw_identification_result.mat'),'result');
disp(modeTable); disp(metricTable);
if makePlots, plot_raw_comparison(result,out); end
end

function plot_raw_comparison(result,out)
names={'raw_time_structured','raw_band_structured'};
labels={'Direct time, rank constrained','Direct band, rank constrained'};
if isfield(result.models,'lrm_reference')
    names{end+1}='lrm_reference'; labels{end+1}='LRM, same training data';
end
frequency=zeros(4,numel(names)); damping=frequency; errors=zeros(3,numel(names));
for k=1:numel(names)
    frequency(:,k)=result.models.(names{k}).frequency_hz;
    damping(:,k)=result.models.(names{k}).damping;
    errors(:,k)=100*result.evaluation.(names{k}).metrics(3:5).';
end
% Reference frequencies are displayed only after fitting; never fit inputs.
reference=[286.43;344.58;488.27;705.67];
fig=figure('Color','w','Position',[60,60,1400,900]); tiledlayout(2,2,'TileSpacing','compact');
nexttile; plot(1:4,frequency-reference,'o-','LineWidth',1.5); grid on;
xticks(1:4); xlabel('Mode'); ylabel('Frequency minus source nominal (Hz)');
title('Four identified modes; fifth remains unmodeled'); legend(labels,'Location','northwest');
nexttile; plot(1:4,damping,'o-','LineWidth',1.5); hold on; yline(.01,'k--','Source nominal');
grid on; xticks(1:4); xlabel('Mode'); ylabel('Damping ratio'); title('Mode 4 reaches the 720 Hz fit bound');
nexttile; bar(errors); set(gca,'YScale','log'); grid on;
xticklabels({'5-200 Hz','200-720 Hz','>720-800 Hz'}); ylabel('Held-out output error (%)');
title('Same unused 30% record; fixed Hann-window band scores');
nexttile; t=(0:result.sample_count-1)'*result.sample_time;
ix=find(t>=t(end)-.15);
plot(t(ix),1e6*result.output(ix,1),'k','LineWidth',1); hold on;
for k=1:2
    p=result.evaluation.(names{k}).prediction; plot(t(ix),1e6*p(ix,1),'LineWidth',1);
end
grid on; xlabel('Time (s)'); ylabel('Output 1 (micro-units)');
title('Held-out free simulation; no output correction');
legend([{'Recorded output'},labels(1:2)],'Location','best');
exportgraphics(fig,fullfile(out,'raw_vs_lrm.png'),'Resolution',160);
end
