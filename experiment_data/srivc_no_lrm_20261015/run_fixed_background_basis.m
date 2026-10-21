function campaign=run_fixed_background_basis
% Fix below-band pole locations to avoid poorly informed denominator updates.
% Their real numerator matrices remain estimated. All in-band poles use SRIVC.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
S=load(fullfile(out,'mode_repair_selected.mat'),'result');parent=S.result;
S=load(fullfile(base,'delay_identification_20261016/common_delay_srivc.mat'),'result');C=S.result;
f=C.f;H=C.trainingFRF;W=C.weights;phase=exp(1i*2*pi*f.*parent.delay_matrix_seconds(:).');
target=H.*phase;seed=parent.selected.fit;models=struct('tag',{},'fit',{},'H',{},'cost',{});
for anchor=[5 10 21]
    fn=[seed.frequency_hz;anchor];z=[seed.damping;.04];bounds=[seed.bounds_hz;1 24];
    opts=struct('fixed_modes',[1 numel(fn)],'max_iterations',250,'score_tolerance',1e-6,'verbose',false);
    fit=frf_srivc_fit(f,target,fn,z,bounds,W,opts);pred=fit.prediction./phase;
    cost=norm((pred-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
    models(end+1)=struct('tag',sprintf('fixed_background_anchor_%g',anchor),'fit',fit,'H',pred,'cost',cost); %#ok<AGROW>
    fprintf('Fixed bg %g: train %.6f%%, score %.3g, converged %d, reason %s\n',anchor,100*sqrt(cost),fit.score,fit.converged,fit.reason);
end
rows=cell(numel(models),7);
for k=1:numel(models)
    m=models(k);ix=f<100;ch=[4 7];rows(k,:)={m.tag,m.fit.converged,m.fit.score,100*sqrt(m.cost), ...
        100*norm(m.H-C.welchValidation,'fro')/norm(C.welchValidation,'fro'), ...
        100*norm(m.H-C.lrmValidation,'fro')/norm(C.lrmValidation,'fro'), ...
        100*norm(m.H(ix,ch)-C.lrmValidation(ix,ch),'fro')/norm(C.lrmValidation(ix,ch),'fro')};
end
summary=cell2table(rows,'VariableNames',{'model','converged','score','weighted_train_percent', ...
    'welch_validation_percent','lrm_validation_percent','weak_cross_lrm_25_100_percent'});
campaign=struct('models',models,'summary',summary,'f',f,'selection','Lowest training cost among converged candidates', ...
    'interpretation','Fixed below-band rational bases; fitted numerator matrices, not identified gravity-spring frequencies.');
save(fullfile(out,'fixed_background_campaign.mat'),'campaign','-v7.3');writetable(summary,fullfile(out,'fixed_background_summary.csv'));disp(summary);
valid=find(arrayfun(@(q)q.fit.converged,models));assert(~isempty(valid),'All fixed-background trials failed.');
[~,b]=min([models(valid).cost]);m=models(valid(b));result=parent;
result.selected.fit=m.fit;result.selected.H=m.H;result.selected.cost=m.cost;result.selected.gradient=[];
result.label='fixed_background_basis_and_residual_modes';result.interpretation=campaign.interpretation;
result.parent_file=fullfile(out,'mode_repair_selected.mat');
result.selection=campaign.selection;
save(fullfile(out,'fixed_background_selected.mat'),'result','-v7.3');
end
