function audit=audit_retest_sensitivity
% Re-estimate on the other half for sensitivity only; never replace training fit.
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
addpath(fullfile(root,'experiment_data','srivc_no_lrm_20261015'));
tags={'m0p6','m0p3','0','0p5','1'};rows={};audit=cell(1,5);
opts=struct('fixed_modes',[],'max_iterations',300,'score_tolerance',1e-6,'verbose',false);
for k=1:5
    s=load(fullfile(out,['model_z' tags{k} '.mat']),'result');r=s.result;fit=r.selected.fit;
    phase=exp(2i*pi*r.f.*r.delay_matrix_seconds(:).');
    v=frf_srivc_fit(r.f,r.validationFRF.*phase,fit.frequency_hz,fit.damping,fit.bounds_hz,r.weights,opts);
    audit{k}=struct('position_mm',r.position_mm,'fit',v,'scope','Conditional split sensitivity; not confidence limits and not used for selection');
    for m=1:numel(fit.frequency_hz)
        rows(end+1,:)={r.position_mm,m,fit.frequency_hz(m),v.frequency_hz(m), ...
            v.frequency_hz(m)-fit.frequency_hz(m),fit.damping(m),v.damping(m),v.damping(m)-fit.damping(m),v.converged,v.score}; %#ok<AGROW>
    end
    fprintf('Split refit z%g score%.3g converged%d\n',r.position_mm,v.score,v.converged);
end
T=cell2table(rows,'VariableNames',{'position_mm','component','train_frequency_hz','check_frequency_hz','frequency_difference_hz', ...
    'train_damping','check_damping','damping_difference','check_converged','check_score'});
save(fullfile(out,'split_sensitivity.mat'),'audit','-v7.3');writetable(T,fullfile(out,'split_parameter_sensitivity.csv'));
end
