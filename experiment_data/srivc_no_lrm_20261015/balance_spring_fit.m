function campaign=balance_spring_fit
% Retain inverse-|H| weights within each band, with an explicit band multiplier.
% Choose the largest low-band multiplier meeting a TRAINING upper-band gate.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'spring_identification_20261016');
S=load(fullfile(out,'spring_selected.mat'),'result');parent=S.result;
S=load(fullfile(base,'error_repair_20261016','final_channel_delay_srivc.mat'),'result');old=S.result;
f=parent.f;H=parent.trainingFRF;V=parent.validationFRF;W0=parent.weights;lo=f<25;hi=~lo;seed=parent.selected.fit;
phase=exp(1i*2*pi*f.*parent.delay_matrix_seconds(:).');target=H.*phase;
oldUpper=100*norm((old.selected.H-H(hi,:)).*sqrt(W0(hi,:)),'fro')/norm(H(hi,:).*sqrt(W0(hi,:)),'fro');
gate=1.10*oldUpper;
models=struct('alpha',{},'fit',{},'H',{},'weights',{},'cost',{},'upper_train_percent',{});rows={};
for alpha=[1 .3 .1 .03 .01 .003]
    W=W0;W(lo,:)=alpha*W(lo,:);
    opts=struct('fixed_modes',[],'max_iterations',250,'score_tolerance',1e-6,'verbose',false);
    fit=frf_srivc_fit(f,target,seed.frequency_hz,seed.damping,seed.bounds_hz,W,opts);pred=fit.prediction./phase;
    cost=norm((pred-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
    upper=100*norm((pred(hi,:)-H(hi,:)).*sqrt(W0(hi,:)),'fro')/norm(H(hi,:).*sqrt(W0(hi,:)),'fro');
    models(end+1)=struct('alpha',alpha,'fit',fit,'H',pred,'weights',W,'cost',cost,'upper_train_percent',upper); %#ok<AGROW>
    rows(end+1,:)={alpha,fit.converged,fit.score,100*sqrt(cost),upper, ...
        100*norm(pred(lo,:)-V(lo,:),'fro')/norm(V(lo,:),'fro'), ...
        100*norm(pred(hi,:)-V(hi,:),'fro')/norm(V(hi,:),'fro')}; %#ok<AGROW>
    fprintf('Spring alpha %g: score %.3g upper train %.3f gate %.3f\n',alpha,fit.score,upper,gate);
end
summary=cell2table(rows,'VariableNames',{'low_band_multiplier','converged','score','train_weighted_percent', ...
    'upper_weighted_train_percent','low_validation_percent','upper_validation_percent'});
campaign=struct('models',models,'summary',summary,'upper_gate_percent',gate,'old_upper_train_percent',oldUpper, ...
    'low_weighted_energy_fraction',norm(H(lo,:).*sqrt(W0(lo,:)),'fro')^2/norm(H.*sqrt(W0),'fro')^2);
save(fullfile(out,'spring_balance.mat'),'campaign','-v7.3');writetable(summary,fullfile(out,'spring_balance.csv'));disp(summary);
valid=find(arrayfun(@(m)m.fit.converged && m.upper_train_percent<=gate,models));
if isempty(valid),warning('No converged model retains upper band within the predeclared training gate.');return;end
m=models(valid(1));result=parent;result.selected.fit=m.fit;result.selected.H=m.H;result.selected.cost=m.cost;
result.selected.tag=sprintf('band_balanced_three_springs_alpha_%g',m.alpha);result.weights=m.weights;
result.band_balance_alpha=m.alpha;result.upper_training_gate_percent=gate;
result.selection='Largest tested low-band multiplier with converged SRIVC and upper-band weighted training error <=1.10 times the prior fit. No validation-based choice.';
save(fullfile(out,'spring_balanced_selected.mat'),'result','-v7.3');
end
