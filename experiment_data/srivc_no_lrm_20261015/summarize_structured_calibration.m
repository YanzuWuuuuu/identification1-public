function summary=summarize_structured_calibration()
% Summarize frozen calibration runs; no fitting or validation-based selection.
% Original physical metrics remain in the individual evaluator outputs.
base=fileparts(mfilename('fullpath'));
out=fullfile(base,'structured_calibration_20261016');
tags={'zm0p5','z0','z0p5','z1','z0_lowhalf30','z0_lowhalf15'};
calibrations=cell(size(tags));evaluations=cell(size(tags));
for k=1:numel(tags)
    cp=fullfile(out,['calibration_' tags{k} '.mat']);
    ep=fullfile(out,['calibration_' tags{k} '_evaluation.mat']);
    assert(isfile(cp) && isfile(ep),'Missing calibration/evaluation for %s. Run its evaluator first.',tags{k});
    a=load(cp,'calib');b=load(ep,'evaluation');
    calibrations{k}=a.calib;evaluations{k}=b.evaluation;
    assert(abs(a.calib.position_mm-b.evaluation.position_mm)<1e-12,'Position mismatch for %s.',tags{k});
    row=makeRow(tags{k},a.calib,b.evaluation);
    if k==1,rows=repmat(row,numel(tags),1);end
    rows(k)=row;
end
assert(isequal([rows(1:4).position_mm],[-.5 0 .5 1]),'Unexpected normal-position ordering.');
normal=struct2table(rows(1:4));
% Preserve all original per-band, per-block values without recalculation.
allMetrics=table;
for k=1:numel(tags)
    t=evaluations{k}.metrics;
    t=addvars(t,repmat(string(tags{k}),height(t),1),repmat(calibrations{k}.position_mm,height(t),1), ...
        repmat(rows(k).low_lrm_half_window,height(t),1),'Before',1, ...
        'NewVariableNames',{'tag','position_mm','low_lrm_half_window'});
    if k==1,allMetrics=t;else,allMetrics=[allMetrics;t];end %#ok<AGROW>
end

% Sensitivity is an ordered diagnostic, not a search for the best window.
sensitivityIndices=[2 5 6];refCal=calibrations{2};refFit=refCal.fits.joint;
[refDy,refDu,refC0]=relativeDelays(refFit.delay_seconds);
for j=1:numel(sensitivityIndices)
    k=sensitivityIndices(j);s=rows(k);fit=calibrations{k}.fits.joint;
    [dy,du,c0]=relativeDelays(fit.delay_seconds);
    s.output_null_angle_to_half45_deg=directionAngle(fit.Ny,refFit.Ny);
    s.input_null_angle_to_half45_deg=directionAngle(fit.Nu,refFit.Nu);
    for c=2:4
        s.(sprintf('dy%d_minus_half45_us',c))=1e6*(dy(c)-refDy(c));
        s.(sprintf('du%d_minus_half45_us',c))=1e6*(du(c)-refDu(c));
    end
    s.c0_minus_half45_us=1e6*(c0-refC0);
    s.relative_delay_change_l2_us=1e6*norm([dy(2:4)-refDy(2:4);du(2:4).'-refDu(2:4).']);
    s.relative_delay_change_max_us=1e6*max(abs([dy(2:4)-refDy(2:4);du(2:4).'-refDu(2:4).']));
    s.delay_matrix_change_frobenius_us=1e6*norm(fit.delay_seconds-refFit.delay_seconds,'fro');
    s.interpretation='Front-end sensitivity only; frozen weights/model/high poles/start rule; no winner selected';
    if j==1,sensitivityRows=repmat(s,numel(sensitivityIndices),1);end
    sensitivityRows(j)=s;
end
sensitivity=struct2table(sensitivityRows);
assert(isequal(sensitivity.low_lrm_half_window(:).',[45 30 15]),'Expected 45/30/15 window diagnostics.');

files={};p=fullfile(out,'structured_calibration_positions_summary.csv');writetable(normal,p);files{end+1}=p;
p=fullfile(out,'structured_calibration_frontend_sensitivity.csv');writetable(sensitivity,p);files{end+1}=p;
p=fullfile(out,'structured_calibration_all_original_metrics.csv');writetable(allMetrics,p);files{end+1}=p;

colors=[.58 .58 .58;.03 .27 .52;.50 .18 .59];
blocks={'logical_inputs_to_N','N_to_logical_outputs','N_to_N'};
names={'Logical inputs to nullspace output','Nullspace input to logical outputs','Nullspace input to nullspace output'};
ratios=zeros(4,3,3);increment=zeros(4,3);
for k=1:4
    for b=1:3
        rmsValues=blockRms(evaluations{k}.metrics,'3-25',blocks{b});
        assert(rmsValues(1)>0,'Zero B0 reference for position %g, block %s.',normal.position_mm(k),blocks{b});
        ratios(k,:,b)=rmsValues/rmsValues(1);
        increment(k,b)=20*log10(max(rmsValues(2),realmin)/max(rmsValues(3),realmin));
    end
end
fig=figure('Visible','off','Color','w','Position',[80 80 1550 500]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for b=1:3
    ax=nexttile(tl);bars=bar(ax,normal.position_mm,ratios(:,:,b),'grouped');
    for m=1:3,bars(m).FaceColor=colors(m,:);bars(m).EdgeColor='none';end
    hold(ax,'on');yline(ax,1,':','Color',[.4 .4 .4],'HandleVisibility','off');
    style(ax);xticks(ax,normal.position_mm);xlabel(ax,'Position z (mm)');
    ylabel(ax,'Validation RMS / same-position B0 RMS');
    title(ax,names{b},'FontWeight','normal');ylim(ax,[0 max(1.05,1.1*max(ratios(:,:,b),[],'all'))]);
    if b==1,legendHandles=bars;end
end
lg=legend(legendHandles,{'B0: nominal static','B1: estimated static','B4: same static + delays'}, ...
    'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
title(tl,'3-25 Hz validation responses at four normal positions','FontWeight','normal');
subtitle(tl,{'Ratio display only: B0 = 1 separately at each position and for each block.', ...
    'Absolute recorded-coordinate RMS remains in the CSV; these responses include spring dynamics.'},'FontSize',10);
files=[files,exportFigure(fig,out,'structured_calibration_low_frequency_ratios')];clear cleanup;

fig=figure('Visible','off','Color','w','Position',[80 80 1550 500]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for b=1:3
    ax=nexttile(tl);plot(ax,normal.position_mm,increment(:,b),'o-','Color',colors(3,:), ...
        'MarkerFaceColor',colors(3,:),'MarkerSize',5,'LineWidth',1.0);hold(ax,'on');
    yline(ax,0,'--','Color',[.4 .4 .4]);style(ax);xticks(ax,normal.position_mm);
    lim=max(.1,1.15*max(abs(increment(:,b))));ylim(ax,[-lim lim]);
    xlabel(ax,'Position z (mm)');ylabel(ax,'20 log_{10}(B1 RMS / B4 RMS) (dB)');
    title(ax,names{b},'FontWeight','normal');
end
title(tl,'Increment from relative-delay alignment, 3-25 Hz','FontWeight','normal');
subtitle(tl,{'B1 and B4 share exactly the same estimated null vectors; positive values mean lower observed response.', ...
    'Zero is marked; each panel spans at least +/-0.1 dB. Numerical changes do not establish hardware significance.'},'FontSize',10);
files=[files,exportFigure(fig,out,'structured_calibration_delay_increment')];clear cleanup;

matPath=fullfile(out,'structured_calibration_summary.mat');files{end+1}=matPath;
summary=struct('positions',normal,'frontend_sensitivity',sensitivity,'original_metrics',allMetrics, ...
    'tags',{tags},'files',{files},'low_frequency_ratios',ratios,'delay_increment_db',increment, ...
    'low_frequency_blocks',{blocks}, ...
    'notes',{{ ...
    'All values use saved second-half validation FRFs; no model is selected using these summaries.', ...
    'output_null means nominal logical inputs to the null output; input_null means null input to nominal logical outputs.', ...
    'RMS columns retain the evaluator recorded-coordinate units and original frequency integration.', ...
    'c0 = ell(1,1); dy_j = ell(j,1)-ell(1,1); du_k = ell(1,k)-ell(1,1). These are effective model delays, not certified hardware delays.', ...
    'The common delay after causal alignment is reported separately from c0.', ...
    '260-310 Hz ratios are observed total responses relative to B0, not isolated flexible-mode retention.', ...
    'grad_inf is an unconstrained full gradient and includes frozen delay coordinates for the fixed candidate; KKT is the stationarity diagnostic.', ...
    '45/30/15 window changes are front-end sensitivity diagnostics only; weights/model/high poles/start rule remain frozen.', ...
    'Normal-position plots exclude the anomalous z=-1 mm record; no re-estimation or digital compensation is implemented.'}});
save(matPath,'summary','-v7.3');disp(normal);disp(sensitivity);
end

function row=makeRow(tag,calib,evaluation)
fit=calib.fits.joint;fixed=calib.fits.fixed;
row=struct('tag',tag,'position_mm',calib.position_mm, ...
    'low_lrm_half_window',optionalField(calib,'low_lrm_half_window',45), ...
    'joint_provisional',evaluation.joint_candidate_provisional);
blocks={'logical_inputs_to_N','N_to_logical_outputs','N_to_N'};
short={'output_null','input_null','NN'};methods={'B0','B1','B4'};
for b=1:3
    values=blockRms(evaluation.metrics,'3-25',blocks{b});
    for m=1:3,row.(sprintf('low_%s_%s_rms',short{b},methods{m}))=values(m);end
    row.(['low_' short{b} '_B0_to_B1_reduction_db'])=20*log10(max(values(1),realmin)/max(values(2),realmin));
    row.(['low_' short{b} '_B1_to_B4_reduction_db'])=20*log10(max(values(2),realmin)/max(values(3),realmin));
end
for b=[1 3]
    values=blockRms(evaluation.metrics,'260-310',blocks{b});
    for m=1:3
        row.(sprintf('target_%s_%s_rms',short{b},methods{m}))=values(m);
        row.(sprintf('target_%s_%s_over_B0',short{b},methods{m}))=values(m)/max(values(1),realmin);
    end
end
row.output_null_angle_to_nominal_deg=directionAngle(fit.Ny,calib.nominal.ny);
row.input_null_angle_to_nominal_deg=directionAngle(fit.Nu,calib.nominal.nu);
row.fixed_to_joint_output_angle_deg=directionAngle(fixed.Ny,fit.Ny);
row.fixed_to_joint_input_angle_deg=directionAngle(fixed.Nu,fit.Nu);
[dy,du,c0]=relativeDelays(fit.delay_seconds);
for j=2:4,row.(sprintf('dy%d_us',j))=1e6*dy(j);row.(sprintf('du%d_us',j))=1e6*du(j);end
row.c0_us=1e6*c0;
row.common_delay_after_alignment_us=1e6*evaluation.common_delay_after_alignment_seconds;
row.joint_exitflag=fit.exitflag;row.fixed_exitflag=fixed.exitflag;
row.joint_KKT=optionalField(fit,'kkt_firstorderopt',NaN);
row.fixed_KKT=optionalField(fixed,'kkt_firstorderopt',NaN);
row.joint_unconstrained_full_gradient=fit.grad_inf;
row.fixed_unconstrained_full_gradient=fixed.grad_inf;
row.analytic_gradient_relative_check=optionalField(calib,'gradient_check',NaN);
row.joint_feasible=optionalField(fit,'feasible',false);
row.joint_spatial_rank_valid=optionalField(fit,'spatial_rank_valid',false);
sigma=optionalField(fit,'normalized_spatial_sigma3_over1',[NaN NaN]);
assert(numel(sigma)==2,'Expected two normalized spatial rank diagnostics.');
row.output_spatial_sigma3_over1=sigma(1);row.input_spatial_sigma3_over1=sigma(2);
agreement=optionalField(calib,'start_agreement_deg',[NaN NaN]);
assert(numel(agreement)==2,'Expected output/input start-agreement angles.');
row.two_start_output_angle_deg=agreement(1);row.two_start_input_angle_deg=agreement(2);
conditions=optionalField(calib,'high_scaled_QR_condition',NaN);
row.high_scaled_QR_condition_max=max(conditions(:));
row.training_weighted_error_percent=fit.train_weighted_error_pct;
fitTable=evaluation.fit_metrics;index=strcmp(string(fitTable.candidate),'joint');
assert(nnz(index)==1,'Expected exactly one joint fit-metrics row.');
row.validation_integrated_relative_error_percent=fitTable.validation_integrated_relative_error_percent(index);
end

function values=blockRms(metrics,band,block)
values=zeros(1,3);methods={'B0','B1','B4'};
for m=1:3
    mask=strcmp(string(metrics.band_hz),band) & strcmp(string(metrics.block),block) & strcmp(string(metrics.method),methods{m});
    assert(nnz(mask)==1,'Expected exactly one metric for %s/%s/%s.',band,block,methods{m});
    values(m)=metrics.validation_response_rms(mask);
end
assert(all(isfinite(values)) && all(values>=0),'Invalid saved RMS values.');
end

function [dy,du,c0]=relativeDelays(ell)
assert(isequal(size(ell),[4 4]) && all(isfinite(ell(:))),'Invalid delay matrix.');
c0=ell(1,1);dy=ell(:,1)-c0;du=ell(1,:)-c0;
assert(norm(ell-(c0+dy+du),'fro')<1e-10,'Delay matrix is not input/output additive.');
end

function value=directionAngle(a,b)
a=a(:);b=b(:);assert(norm(a)>0 && norm(b)>0,'Zero direction.');
value=acosd(min(1,max(0,abs(a'*b)/(norm(a)*norm(b)))));
end

function value=optionalField(source,name,fallback)
if isfield(source,name),value=source.(name);else,value=fallback;end
end

function style(ax)
set(ax,'FontName','Times New Roman','FontSize',10,'Box','on','LineWidth',.65,'GridAlpha',.18);
grid(ax,'on');
end

function files=exportFigure(fig,out,stem)
png=fullfile(out,[stem '.png']);pdf=fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);exportgraphics(fig,pdf,'ContentType','vector');files={png,pdf};
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
