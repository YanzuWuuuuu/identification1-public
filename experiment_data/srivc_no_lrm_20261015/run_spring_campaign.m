function outputs=run_spring_campaign
% Reproduce this bounded campaign; original measurements/baselines stay read-only.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'spring_identification_20261016');
prepare_spring_data;
run_spring_identification;
check_spring_sensitivity;
balance_spring_fit;
final=refine_spring_delays;
hi=final.f>=25;W=final.weights;H=final.trainingFRF;
upper=100*norm((final.selected.H(hi,:)-H(hi,:)).*sqrt(W(hi,:)),'fro')/norm(H(hi,:).*sqrt(W(hi,:)),'fro');
assert(final.selected.fit.converged && final.selected_matches_optimizer && upper<=final.upper_training_gate_percent, ...
    'Final model does not meet the declared training gate.');
copyfile(fullfile(out,'spring_selected.mat'),fullfile(out,'spring_unbalanced_selected.mat'));
result=final;save(fullfile(out,'spring_selected.mat'),'result','-v7.3');
fit=final.selected.fit;phase=exp(1i*2*pi*final.f.*final.delay_matrix_seconds(:).');
options=struct('fixed_modes',4:numel(fit.frequency_hz),'max_iterations',250,'score_tolerance',1e-6,'verbose',false);
spring_val_refit=frf_srivc_fit(final.f,final.validationFRF.*phase,fit.frequency_hz,fit.damping,fit.bounds_hz,final.weights,options);
save(fullfile(out,'spring_validation_refit.mat'),'spring_val_refit');
outputs.verification=summarize_spring_identification;
spring_checks=plot_spring_identification;save(fullfile(out,'spring_plot_checks.mat'),'spring_checks');
outputs.plots=spring_checks;outputs.spatial=audit_spring_mode_shapes;
end
