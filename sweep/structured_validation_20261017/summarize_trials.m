function summary=summarize_trials
out=fileparts(mfilename('fullpath'));tags={'zm0p6','zm0p3','z0','z0p5','z1'};z=[-.6 -.3 0 .5 1];
rows={};calibrationAll=table;predictionAll=table;
for p=1:5
    e=evaluate_structured_trial(tags{p});a=e.calibration_metrics;b=e.fit_metrics;
    a.position_mm=repmat(z(p),height(a),1);b.position_mm=repmat(z(p),height(b),1);
    calibrationAll=[calibrationAll;a];predictionAll=[predictionAll;b]; %#ok<AGROW>
    source=load(fullfile(out,['trial_' tags{p} '.mat']));r=source.result;
    assert(r.delay.nonnegative_total_delay,'Negative total delays: physical interpretation rejected.');
    fixedError=@(name,block) b.error_pct(strcmp(b.model,name)&strcmp(b.block,block)&b.fmin_Hz==25&b.fmax_Hz==800);
    change=@(cal,block) a.change_dB_vs_nominal(strcmp(a.calibrator,cal)&strcmp(a.block,block)&a.fmin_Hz==25&a.fmax_Hz==80);
    retain=@(cal) a.ratio_to_nominal(strcmp(a.calibrator,cal)&strcmp(a.block,'N_to_N')&a.fmin_Hz==260&a.fmax_Hz==310);
    rows(end+1,:)={z(p),fixedError('SRIVC_free','physical'),fixedError('Structured_updated_delay','physical'), ...
        fixedError('SRIVC_free','N_to_N'),fixedError('Structured_updated_delay','N_to_N'), ...
        change('Updated_static','N_to_logical'),change('Updated_static_and_delay','N_to_logical'),change('Direct_band_static','N_to_logical'), ...
        change('Updated_static','logical_to_N'),change('Updated_static_and_delay','logical_to_N'),change('Direct_band_static','logical_to_N'), ...
        retain('Updated_static_and_delay'),retain('Direct_band_static'), ...
        min(r.delay.delay_matrix_seconds(:))*1e3,max(r.delay.delay_matrix_seconds(:))*1e3,r.delay.optimization.firstorderopt}; %#ok<AGROW>
end
summary=cell2table(rows,'VariableNames',{'z_mm','baseline_physical_error_pct','structured_physical_error_pct', ...
    'baseline_NN_error_pct','structured_NN_error_pct','N_to_logical_static_change_dB','N_to_logical_dynamic_change_dB','N_to_logical_direct_change_dB', ...
    'logical_to_N_static_change_dB','logical_to_N_dynamic_change_dB','logical_to_N_direct_change_dB', ...
    'NN_first_resonance_band_ratio_structured','NN_first_resonance_band_ratio_direct', ...
    'min_effective_total_delay_ms','max_effective_total_delay_ms','firstorderopt'});
writetable(summary,fullfile(out,'five_position_summary.csv'));
writetable(calibrationAll,fullfile(out,'all_calibration_metrics.csv'));
writetable(predictionAll,fullfile(out,'all_prediction_metrics.csv'));
save(fullfile(out,'five_position_summary.mat'),'summary','calibrationAll','predictionAll');
disp(summary);
fig=figure('Visible','off','Color','w','Position',[100 100 1400 780]);tl=tiledlayout(fig,2,2,'Padding','compact');
cols=[.1 .6 .5;.05 .28 .65;.55 .3 .7];
for block=1:2
    ax=nexttile(tl);if block==1,V=[summary.N_to_logical_static_change_dB,summary.N_to_logical_dynamic_change_dB,summary.N_to_logical_direct_change_dB];ttl='N input to logical outputs';
    else,V=[summary.logical_to_N_static_change_dB,summary.logical_to_N_dynamic_change_dB,summary.logical_to_N_direct_change_dB];ttl='Logical inputs to N output';end
    h=bar(ax,z,V);for k=1:3,h(k).FaceColor=cols(k,:);end
    title(ax,[ttl ': 25-80 Hz']);ylabel(ax,'Change from nominal geometry (dB)');xlabel(ax,'z position (mm)');grid(ax,'on');
    legend(ax,{'Structured static','Structured + delays','Direct band static'},'Location','best','FontSize',9);
end
ax=nexttile(tl);plot(ax,z,summary.baseline_NN_error_pct,'o--','Color',[.55 .3 .7],'LineWidth',1.4);hold(ax,'on');
plot(ax,z,summary.structured_NN_error_pct,'o-','Color',[.05 .28 .65],'LineWidth',1.4);grid(ax,'on');
title(ax,'Fixed geometric N to N prediction: 25-800 Hz');xlabel(ax,'z position (mm)');ylabel(ax,'Complex FRF error (%)');legend(ax,{'Free SRIVC','Structure + delays'},'Location','best');
ax=nexttile(tl);plot(ax,z,100*summary.NN_first_resonance_band_ratio_structured,'o-','Color',[.05 .28 .65],'LineWidth',1.4);hold(ax,'on');
plot(ax,z,100*summary.NN_first_resonance_band_ratio_direct,'s-.','Color',[.55 .3 .7],'LineWidth',1.4);grid(ax,'on');yline(ax,100,':k');
title(ax,'Measured N to N band response: 260-310 Hz');xlabel(ax,'z position (mm)');ylabel(ax,'RMS ratio to nominal (%)');
legend(ax,{'Structure + delays','Direct band static'},'Location','best');
title(tl,{'Five-position check on measured FRFs','Training-derived directions; same-record held-out development check, not a hardware decoupling test'},'FontSize',13);
exportgraphics(fig,fullfile(out,'five_position_check.png'),'Resolution',170);exportgraphics(fig,fullfile(out,'five_position_check.pdf'),'ContentType','vector');close(fig);
plot_calibration_check('z0');
end
