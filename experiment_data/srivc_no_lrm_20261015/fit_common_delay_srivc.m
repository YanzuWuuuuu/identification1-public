function result=fit_common_delay_srivc(resume)
% Shared effective delay: outer training-cost profile, inner genuine SRIVC.
% All inner modal parameters remain adjustable, except the fixed background.
if nargin<1,resume=false;end
base=fileparts(mfilename('fullpath'));out=fullfile(base,'delay_identification_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'frf_preprocessing_comparison.mat'),'comparison');C=S.comparison;
f=C.f;H=C.lrmTrain;floorValue=max(1e-3*max(abs(H),[],1),eps);
W=1./max(abs(H),floorValue);W=W/median(W(:));
seed=C.models(strcmp({C.models.tag},'lrm_inverse_magnitude')).fit;
opts=struct('fixed_modes',1,'max_iterations',200,'score_tolerance',1e-6,'verbose',false);
% Warm pole initialization is fixed across candidates; no held-out tuning.
frequency0=seed.frequency_hz;damping0=seed.damping;bounds=seed.bounds_hz;
normalizer=norm(H.*sqrt(W),'fro')^2;
cache=struct('tau_us',{},'fit',{},'cost',{},'H',{},'gradient',{});
if resume
    saved=load(fullfile(out,'common_delay_checkpoint.mat'),'cache','normalizer');
    assert(abs(saved.normalizer-normalizer)<1e-12*normalizer && ...
        isequal(saved.cache(1).fit.quadratic_weights,W),'Checkpoint does not match current weights.');
    cache=saved.cache;
end
grid=(0:100:1000)';values=zeros(size(grid));
for k=1:numel(grid),values(k)=objective(grid(k));end
[~,best]=min(values);lo=grid(max(best-1,1));hi=grid(min(best+1,numel(grid)));
if lo<hi
    bestTau=fminbnd(@objective,lo,hi,optimset('TolX',.01,'Display','off','MaxFunEvals',35));
    objective(bestTau);
end
costs=[cache.cost];[~,chosen]=min(costs);selected=cache(chosen);
assert(isfinite(selected.cost)&&selected.fit.converged,'No converged delay profile solution.');
trace=zeros(numel(cache),8);
for k=1:numel(cache)
    q=cache(k);ridge=max([0;q.fit.history(2:end,5)]);
    trace(k,:)=[q.tau_us,q.cost,100*q.fit.relative_error,q.fit.converged, ...
        q.fit.score,q.fit.iterations,q.gradient,ridge];
end
trace=array2table(trace,'VariableNames',{'tau_us','normalized_training_cost','training_weighted_error_percent', ...
    'inner_converged','inner_score','inner_iterations','training_profile_gradient_per_us','max_ridge'});
result=struct('f',f,'trainingFRF',H,'welchValidation',C.welchValidation,'lrmValidation',C.lrmValidation, ...
    'weights',W,'baseline',seed,'baselineH',seed.prediction,'selected',selected,'trace',trace, ...
    'coarse_grid_us',grid,'coarse_costs',values,'profile_cache',cache, ...
    'delay_matrix_seconds',repmat(selected.tau_us*1e-6,4,4), ...
    'method','Fixed inverse-magnitude training weights; outer common-delay profile; inner additive SRIVC; selected using training cost only', ...
    'interpretation','Effective phase parameter; not certified hardware delay or unbiased closed-loop estimate');
save(fullfile(out,'common_delay_srivc.mat'),'result','-v7.3');
writetable(trace,fullfile(out,'common_delay_profile.csv'));
fprintf('SELECTED common delay %.6f us; inner score %.3g; weighted train %.4f%%; baseline %.4f%%\n', ...
    selected.tau_us,selected.fit.score,100*selected.fit.relative_error,100*seed.relative_error);
disp(trace);

    function cost=objective(tauUs)
        old=find(abs([cache.tau_us]-tauUs)<1e-9,1);
        if ~isempty(old),cost=cache(old).cost;return;end
        phase=exp(1i*2*pi*f*tauUs*1e-6);corrected=H.*phase;
        fit=frf_srivc_fit(f,corrected,frequency0,damping0,bounds,W,opts);
        h=fit.prediction./phase;res=corrected-fit.prediction;
        rawCost=norm((H-h).*sqrt(W),'fro')^2/normalizer;
        assert(abs(rawCost-fit.relative_error^2)<1e-10,'Delay compensation cost identity failed.');
        grad=2*real(sum(conj(res).*W.*(1i*2*pi*f*1e-6.*corrected),'all'))/normalizer;
        cost=rawCost;
        if ~fit.converged,cost=Inf;end
        cache(end+1)=struct('tau_us',tauUs,'fit',fit,'cost',cost,'H',h,'gradient',grad);
        fprintf('tau=%8.3f us; train=%7.4f%%; inner converged=%d score=%.3g iter=%d\n', ...
            tauUs,100*sqrt(rawCost),fit.converged,fit.score,fit.iterations);
        save(fullfile(out,'common_delay_checkpoint.mat'),'cache','grid','normalizer','-v7.3');
    end
end
