function campaign=run_broad_residual_check
% Test one damped residual component for coherent 180--320 Hz mismatch.
% It is a phenomenological model term, not a newly certified mechanical mode.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
S=load(fullfile(out,'repaired_channel_delay_srivc.mat'),'result');parent=S.result;
S=load(fullfile(base,'delay_identification_20261016/common_delay_srivc.mat'),'result');C=S.result;
f=C.f;H=C.trainingFRF;W=C.weights;phase=exp(1i*2*pi*f.*parent.delay_matrix_seconds(:).');
target=H.*phase;seed=parent.selected.fit;starts=[250 .2;270 .07;320 .2];
models=struct('tag',{},'fit',{},'H',{},'cost',{});
for k=1:size(starts,1)
    fn=[seed.frequency_hz;starts(k,1)];z=[seed.damping;starts(k,2)];bounds=[seed.bounds_hz;150 450];
    opts=struct('fixed_modes',seed.fixed_modes,'max_iterations',450,'score_tolerance',1e-6,'verbose',false);
    fit=frf_srivc_fit(f,target,fn,z,bounds,W,opts);pred=fit.prediction./phase;
    cost=norm((pred-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
    models(end+1)=struct('tag',sprintf('broad_start_%g',starts(k,1)),'fit',fit,'H',pred,'cost',cost); %#ok<AGROW>
    fprintf('%s: train %.5f%%, score %.3g, converged %d, reason %s, added %.4g Hz damping %.4g\n', ...
        models(end).tag,100*sqrt(cost),fit.score,fit.converged,fit.reason,fit.frequency_hz(end),fit.damping(end));
    save(fullfile(out,'broad_residual_checkpoint.mat'),'models','-v7.3');
end
rows=cell(numel(models),10);
for k=1:numel(models)
    m=models(k);local=zeros(1,2);cent=local;
    for j=1:2,ch=1+5*(j-1);ib=f>=180&f<=260;fb=f(ib);[~,ci]=min(abs(H(ib,ch)));ic=ib&abs(f-fb(ci))<=5;
        local(j)=100*norm(m.H(ib,ch)-C.lrmValidation(ib,ch))/norm(C.lrmValidation(ib,ch));
        cent(j)=100*norm(m.H(ic,ch)-C.lrmValidation(ic,ch))/norm(C.lrmValidation(ic,ch));
    end
    rows(k,:)={m.tag,m.fit.converged,m.fit.score,100*sqrt(m.cost), ...
        100*norm(m.H-C.welchValidation,'fro')/norm(C.welchValidation,'fro'), ...
        100*norm(m.H-C.lrmValidation,'fro')/norm(C.lrmValidation,'fro'),local(1),local(2),cent(1),cent(2)};
end
summary=cell2table(rows,'VariableNames',{'model','converged','score','weighted_train_percent', ...
    'welch_validation_percent','lrm_validation_percent','notch1_band','notch2_band','notch1_center','notch2_center'});
campaign=struct('models',models,'summary',summary,'f',f,'selection','Lowest training cost among converged candidates', ...
    'interpretation','One added damped residual component; physical attribution unresolved; fixed W, delays, below-band poles.');
save(fullfile(out,'broad_residual_campaign.mat'),'campaign','-v7.3');writetable(summary,fullfile(out,'broad_residual_summary.csv'));disp(summary);
valid=find(arrayfun(@(q)q.fit.converged,models));if isempty(valid),return;end
[~,b]=min([models(valid).cost]);m=models(valid(b));result=struct('f',f,'selected', ...
    struct('fit',m.fit,'H',m.H,'cost',m.cost,'coordinates_ms',parent.selected.coordinates_ms,'gradient',[]), ...
    'weights',W,'delay_design',parent.delay_design,'delay_matrix_seconds',parent.delay_matrix_seconds, ...
    'label','broad_residual_augmented_srivc','parent_file',fullfile(out,'repaired_channel_delay_srivc.mat'), ...
    'selection',campaign.selection,'interpretation',campaign.interpretation);
save(fullfile(out,'broad_residual_selected.mat'),'result','-v7.3');
end
