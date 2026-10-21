function campaign=run_background_order_check
% One additional below-band rational background tests weak cross-channel tails.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
S=load(fullfile(out,'mode_repair_selected.mat'),'result');parent=S.result;
S=load(fullfile(base,'delay_identification_20261016/common_delay_srivc.mat'),'result');C=S.result;
f=C.f;H=C.trainingFRF;W=C.weights;phase=exp(1i*2*pi*f.*parent.delay_matrix_seconds(:).');
target=H.*phase;seed=parent.selected.fit;
attempts=struct('tag',{},'fit',{},'H',{},'cost',{});
opts=struct('fixed_modes',[],'max_iterations',450,'score_tolerance',1e-6,'verbose',false);
for initial=[5 10 21]
    fn=[seed.frequency_hz;initial];z=[seed.damping;.04];bounds=[seed.bounds_hz;1 24];
    fit=frf_srivc_fit(f,target,fn,z,bounds,W,opts);pred=fit.prediction./phase;
    cost=norm((pred-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
    attempts(end+1)=struct('tag',sprintf('extra_background_start_%g',initial),'fit',fit,'H',pred,'cost',cost); %#ok<AGROW>
    fprintf('Extra bg %g: train %.6f%%, score %.3g, converged %d, reason %s, bg %.4g/%.4g Hz\n', ...
        initial,100*sqrt(cost),fit.score,fit.converged,fit.reason,fit.frequency_hz(1),fit.frequency_hz(end));
    save(fullfile(out,'background_order_checkpoint.mat'),'attempts','-v7.3');
end
rows={};
for k=1:numel(attempts)
    m=attempts(k);ix=f<100;ch=[4 7];
    rows(k,:)={m.tag,m.fit.converged,m.fit.score,m.fit.reason,100*sqrt(m.cost), ...
        100*norm(m.H-C.welchValidation,'fro')/norm(C.welchValidation,'fro'), ...
        100*norm(m.H-C.lrmValidation,'fro')/norm(C.lrmValidation,'fro'), ...
        100*norm(m.H(ix,ch)-C.lrmValidation(ix,ch),'fro')/norm(C.lrmValidation(ix,ch),'fro')};
end
summary=cell2table(rows,'VariableNames',{'attempt','converged','score','reason','weighted_train_percent', ...
    'welch_validation_percent','lrm_validation_percent','weak_cross_lrm_25_100_percent'});
campaign=struct('attempts',attempts,'summary',summary,'f',f,'interpretation', ...
    'Extra background is a below-band nuisance component, not a certified mechanical mode. Delays and W frozen.');
save(fullfile(out,'background_order_campaign.mat'),'campaign','-v7.3');writetable(summary,fullfile(out,'background_order_summary.csv'));
disp(summary);valid=find(arrayfun(@(q)q.fit.converged,attempts));
if isempty(valid),fprintf('No converged extra-background candidate: retain parent.\n');return;end
[~,b]=min([attempts(valid).cost]);chosen=attempts(valid(b));
result=parent;result.selected.fit=chosen.fit;result.selected.H=chosen.H;result.selected.cost=chosen.cost;
result.label='two_low_backgrounds_plus_residual_modes';result.interpretation=campaign.interpretation;
save(fullfile(out,'background_order_selected.mat'),'result','-v7.3');
end
