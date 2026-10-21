function campaign=run_background_repair
% Isolate fixed sub-band background error. Keep data, inverse weight, delays.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'delay_identification_20261016/common_delay_srivc.mat'),'result');C=S.result;
S=load(fullfile(base,'delay_identification_20261016/channel_delay_srivc.mat'),'result');R=S.result;
f=C.f;H=C.trainingFRF;W=C.weights;ell=R.delay_matrix_seconds;
phase=exp(1i*2*pi*f.*ell(:).');target=H.*phase;
seed=R.selected.fit;fn=seed.frequency_hz;z=seed.damping;bounds=seed.bounds_hz;
models=struct('tag',{},'fit',{},'H',{},'delay',{},'cost',{});
models(1)=entry('original_fixed15',seed,R.selected.H);
% Only two low-frequency pole coordinates change in this diagnostic profile.
starts=[15 .15;10 .05];profiles=struct('x',{},'cost',{},'exitflag',{},'output',{});
for k=1:size(starts,1)
    x0=log(starts(k,:));
    [x,cost,flag,output]=fminsearch(@profileCost,x0,optimset('Display','off', ...
        'MaxIter',180,'MaxFunEvals',350,'TolX',1e-7,'TolFun',1e-11));
    profiles(k)=struct('x',x,'cost',cost,'exitflag',flag,'output',output);
end
[~,best]=min([profiles.cost]);bg=exp(profiles(best).x);
fixedFn=fn;fixedZ=z;fixedFn(1)=bg(1);fixedZ(1)=bg(2);
opts=struct('fixed_modes',1:numel(fn),'max_iterations',1,'score_tolerance',1e-6,'verbose',false);
fixed=frf_srivc_fit(f,target,fixedFn,fixedZ,bounds,W,opts);
models(end+1)=entry('profile_background_fixed_other_poles',fixed,fixed.prediction./phase);
% Genuine SRIVC: release all shared poles, including sub-band background.
bounds(1,:)=[1 24];opts.fixed_modes=[];opts.max_iterations=250;
free=frf_srivc_fit(f,target,fixedFn,fixedZ,bounds,W,opts);
models(end+1)=entry('free_background_srivc',free,free.prediction./phase);
rows={};
for k=1:numel(models)
    m=models(k);rows(k,:)={m.tag,m.fit.converged,m.fit.score,m.fit.frequency_hz(1),m.fit.damping(1), ...
        100*sqrt(m.cost),100*norm(m.H-C.welchValidation,'fro')/norm(C.welchValidation,'fro'), ...
        100*norm(m.H-C.lrmValidation,'fro')/norm(C.lrmValidation,'fro')};
end
summary=cell2table(rows,'VariableNames',{'model','converged','score','background_hz','background_damping', ...
    'weighted_train_percent','welch_validation_percent','lrm_validation_percent'});
campaign=struct('models',models,'summary',summary,'background_profiles',profiles,'f',f, ...
    'fixed_conditions','Original LRM training, W inverse magnitude, additive input/output delays frozen', ...
    'selection','Training cost only; validation reused for exploratory diagnostics, not a fresh test set');
save(fullfile(out,'background_repair_campaign.mat'),'campaign','-v7.3');
writetable(summary,fullfile(out,'background_repair_summary.csv'));
eligible=find(arrayfun(@(q)q.fit.converged,models));[~,index]=min([models(eligible).cost]);m=models(eligible(index));
result=R;result.selected.fit=m.fit;result.selected.H=m.H;result.selected.cost=m.cost;
result.selected.gradient=[];
result=rmfield(result,{'optimizer_coordinates_ms','optimizer_cost','optimizer_exitflag','optimizer_output','trace_cache'});
result.method='Background profile followed by constrained SRIVC; input/output delays held fixed';
result.parent_file=fullfile(base,'delay_identification_20261016/channel_delay_srivc.mat');
result.label=m.tag;result.parent='channel_delay_srivc';result.selection=campaign.selection;
save(fullfile(out,'background_repair_selected.mat'),'result','-v7.3');
disp(summary);

    function m=entry(tag,fit,prediction)
        cost=norm((prediction-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
        m=struct('tag',tag,'fit',fit,'H',prediction,'delay',ell,'cost',cost);
    end
    function cost=profileCost(x)
        bgTrial=exp(x);
        if bgTrial(1)<1 || bgTrial(1)>24 || bgTrial(2)<.00003 || bgTrial(2)>.95
            cost=1e3+sum(x.^2);return;
        end
        trialFn=fn;trialZ=z;trialFn(1)=bgTrial(1);trialZ(1)=bgTrial(2);
        q=1i*f./trialFn.';A=1+2*q.*trialZ.'+q.^2;
        B=complex(zeros(numel(f),2*numel(fn)+1));
        for n=1:numel(fn),B(:,2*n-1)=1./A(:,n);B(:,2*n)=q(:,n)./A(:,n);end
        B(:,end)=1;coef=zeros(size(B,2),16);
        for ch=1:16
            X=B.*sqrt(W(:,ch));y=target(:,ch).*sqrt(W(:,ch));
            X=[real(X);imag(X)];y=[real(y);imag(y)];scale=max(vecnorm(X),realmin);
            coef(:,ch)=((X./scale)\y)./scale.';
        end
        cost=norm((B*coef-target).*sqrt(W),'fro')^2/norm(target.*sqrt(W),'fro')^2;
    end
end
