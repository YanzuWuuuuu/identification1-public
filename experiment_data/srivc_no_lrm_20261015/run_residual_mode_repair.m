function campaign=run_residual_mode_repair
% Target residual structures around 124/128 Hz and the repeated 772-Hz peak.
% All cases share the same training set, inverse-magnitude W and seven delays.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
S=load(fullfile(out,'background_and_delays.mat'),'result');parent=S.result;
S=load(fullfile(base,'delay_identification_20261016/common_delay_srivc.mat'),'result');C=S.result;
f=C.f;H=C.trainingFRF;W=C.weights;phase=exp(1i*2*pi*f.*parent.delay_matrix_seconds(:).');
target=H.*phase;seed=parent.selected.fit;
models=struct('tag',{},'fit',{},'H',{},'cost',{},'converged',{});
models(1)=entry('background_and_delays',seed,parent.selected.H);
attempts=struct('tag',{},'fit',{},'H',{},'cost',{},'converged',{});
tags={'split_124_128','add_772','split_124_128_add_772'};
opts=struct('fixed_modes',[],'max_iterations',450,'score_tolerance',1e-6,'verbose',false);
for k=1:3
    fn=seed.frequency_hz;z=seed.damping;bounds=seed.bounds_hz;
    if k==1 || k==3
        fn(2)=124;z(2)=.004;bounds(2,:)=[118 126];
        fn(end+1)=128.7;z(end+1)=.004;bounds(end+1,:)=[126 133];
    end
    if k==2 || k==3
        fn(end+1)=772;z(end+1)=.003;bounds(end+1,:)=[755 790];
    end
    startIndices=[];
    for start=1:2
        damping=z;if start==2,damping(2:end)=.6*damping(2:end);end
        fit=frf_srivc_fit(f,target,fn,damping,bounds,W,opts);
        m=entry(sprintf('%s_start%d',tags{k},start),fit,fit.prediction./phase);
        attempts(end+1)=m;startIndices(end+1)=numel(attempts); %#ok<AGROW>
        fprintf('%s: score %.3g, converged %d, train %.5f%%, reason %s\n', ...
            m.tag,fit.score,fit.converged,100*sqrt(m.cost),fit.reason);
        save(fullfile(out,'mode_repair_checkpoint.mat'),'attempts','-v7.3');
    end
    valid=startIndices([attempts(startIndices).converged]);
    if isempty(valid),valid=startIndices;end
    [~,b]=min([attempts(valid).cost]);m=attempts(valid(b));m.tag=tags{k};models(end+1)=m; %#ok<AGROW>
end
rows={};
for k=1:numel(models)
    m=models(k);rows(k,:)={m.tag,m.converged,m.fit.score,numel(m.fit.frequency_hz), ...
        100*sqrt(m.cost),err(m.H,C.welchValidation,true(size(f))),err(m.H,C.lrmValidation,true(size(f))), ...
        err(m.H,C.lrmValidation,f>=100&f<170),err(m.H,C.lrmValidation,f>=180&f<260), ...
        err(m.H,C.lrmValidation,f>=735)};
end
summary=cell2table(rows,'VariableNames',{'model','converged','score','second_order_components', ...
    'weighted_train_percent','welch_validation_percent','lrm_validation_percent', ...
    'lrm_100_170_percent','lrm_180_260_percent','lrm_735_800_percent'});
campaign=struct('models',models,'attempts',attempts,'summary',summary,'f',f, ...
    'fixed_conditions','Same LRM training, inverse W, seven input/output delays; backgrounds/poles free in SRIVC', ...
    'selection','Training cost among converged candidates; exploratory reused validation');
save(fullfile(out,'mode_repair_campaign.mat'),'campaign','-v7.3');
writetable(summary,fullfile(out,'mode_repair_summary.csv'));disp(summary);
valid=find([models.converged]);[~,b]=min([models(valid).cost]);m=models(valid(b));
result=struct('f',f,'selected',struct('fit',m.fit,'H',m.H,'cost',m.cost, ...
    'coordinates_ms',parent.selected.coordinates_ms,'gradient',[]), ...
    'weights',W,'delay_design',parent.delay_design,'delay_matrix_seconds',parent.delay_matrix_seconds, ...
    'label',m.tag,'parent_file',fullfile(out,'background_and_delays.mat'), ...
    'selection',campaign.selection,'method','Genuine additive SRIVC; shared-pole residual-structure repair; seven additive delays retained');
save(fullfile(out,'mode_repair_selected.mat'),'result','-v7.3');

    function m=entry(tag,fit,prediction)
        cost=norm((prediction-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
        m=struct('tag',tag,'fit',fit,'H',prediction,'cost',cost,'converged',fit.converged);
    end
    function e=err(a,b,ix)
        e=100*norm(a(ix,:)-b(ix,:),'fro')/norm(b(ix,:),'fro');
    end
end
