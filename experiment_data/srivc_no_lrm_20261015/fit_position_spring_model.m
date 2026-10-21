function result=fit_position_spring_model(positionIndex,refineDelay,resumePrevious)
% Fit one position with the same 15-component family and weighting contract.
if nargin<2,refineDelay=true;end
if nargin<3,resumePrevious=false;end
base=fileparts(mfilename('fullpath'));out=fullfile(base,'position_spring_identification_20261016');
S=load(fullfile(out,'position_data.mat'),'campaign');D=S.campaign.datasets{positionIndex};
S=load(fullfile(base,'spring_identification_20261016','spring_selected.mat'),'result');parent=S.result;
seed=parent.selected.fit;f=D.f;H=D.trainingFRF;W=D.weights;E=parent.delay_design;
bounds=[4 24.5;4 24.5;4 24.5;110 140;125 165;260 315;315 352;350 380;475 515;600 700;600 700;695 725;110 140;740 795;150 450];
opts=struct('fixed_modes',[],'max_iterations',300,'score_tolerance',1e-6,'verbose',false);
start=parent.selected.coordinates_ms;
previousPath=fullfile(out,['model_z' D.tag '.mat']);
if resumePrevious && isfile(previousPath)
    previous=load(previousPath,'result');seed=previous.result.selected.fit;start=previous.result.selected.coordinates_ms;
end
phase=exp(1i*2*pi*f.*(E*start*1e-3).');norm2=norm(H.*sqrt(W),'fro')^2;
initials=struct('fit',{},'H',{},'cost',{},'coordinates_ms',{},'gradient',{});
for trial=1:2
    fn=seed.frequency_hz;z=seed.damping;
    if trial==2
        % Local peak seeds are taken from TRAINING FRFs only. Bounds/order
        % remain identical; local singular-value peaks do not label physics.
        targets=[4 5 6 7 8 9 12 14];windows=[112 135;133 160;270 295;318 350;352 375;482 505;698 715;750 790];
        curve=sqrt(sum(abs(H).^2,2));
        for kk=1:numel(targets)
            ix=find(f>=windows(kk,1)&f<=windows(kk,2));[~,jj]=max(curve(ix));fn(targets(kk))=f(ix(jj));
        end
        fn(13)=min(139,max(111,fn(4)+3.7));
    end
    fit=frf_srivc_fit(f,H.*phase,fn,z,bounds,W,opts);prediction=fit.prediction./phase;
    cost=norm((prediction-H).*sqrt(W),'fro')^2/norm2;
    initials(end+1)=struct('fit',fit,'H',prediction,'cost',cost,'coordinates_ms',start,'gradient',[]); %#ok<AGROW>
    fprintf('z=%g start %d: weighted train %.5f%% score %.3g converged %d (%s)\n',D.position_mm,trial,100*sqrt(cost),fit.score,fit.converged,fit.reason);
    save(fullfile(out,['initial_z' D.tag '.mat']),'initials','bounds','-v7.3');
end
valid=find(arrayfun(@(v)v.fit.converged,initials));
if isempty(valid),[~,j]=min([initials.cost]);else,[~,jj]=min([initials(valid).cost]);j=valid(jj);end
chosen=initials(j);seed=chosen.fit;cache=struct('coordinates_ms',{},'fit',{},'H',{},'cost',{},'gradient',{});
outerStatus='not_run';outerFlag=NaN;outerOutput=struct;optimizerX=start;outerMessage='';
if refineDelay && chosen.fit.converged
    options=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'Display','none', ...
        'MaxIterations',65,'MaxFunctionEvaluations',150,'OptimalityTolerance',1e-7,'StepTolerance',1e-7,'ConstraintTolerance',1e-9);
    try
        [optimizerX,~,outerFlag,outerOutput]=fmincon(@objective,start,[E;-E],[ones(16,1);zeros(16,1)], ...
            [],[],[0;-ones(6,1)],[1;ones(6,1)],[],options);
        objective(optimizerX);outerStatus='completed';
    catch err
        outerStatus='interrupted';outerMessage=err.message;fprintf('z=%g delay refinement stopped: %s\n',D.position_mm,err.message);
    end
    if ~isempty(cache),[~,j]=min([cache.cost]);chosen=cache(j);end
end
result=D;result.selected=chosen;result.delay_design=E;result.delay_matrix_seconds=reshape(E*chosen.coordinates_ms,4,4)*1e-3;
result.coordinates_ms=chosen.coordinates_ms;result.common_bounds_hz=bounds;result.initial_trials=initials;
result.outer_status=outerStatus;result.optimizer_exitflag=outerFlag;result.optimizer_output=outerOutput;
result.optimizer_coordinates_ms=optimizerX;result.outer_message=outerMessage;
result.selected_matches_optimizer=norm(chosen.coordinates_ms-optimizerX,inf)<1e-10;
result.inner_converged=chosen.fit.converged;
result.delay_optimization_complete=strcmp(outerStatus,'completed') && outerFlag>0 && result.selected_matches_optimizer;
if chosen.fit.converged && result.delay_optimization_complete,result.status='converged';else,result.status='nonconverged';end
result.selection='Minimum training cost among stationary inner fits when available; no validation-based choice. Same family/bounds/weight rule across positions.';
result.method='LRM half-window 45; all 15 additive second-order components by SRIVC; seven additive delays profiled when inner convergence permits';
result.validation_scope='Exploratory first/second halves of one record at each position; not independent repeat experiments';
save(fullfile(out,['model_z' D.tag '.mat']),'result','-v7.3');
trace=cache;save(fullfile(out,['profile_z' D.tag '.mat']),'trace','-v7.3');
lo=f<25;hi=~lo;
fprintf('RESULT z=%g: %s outer %g, LRM val low %.5f%% upper %.5f%%\n',D.position_mm,result.status,outerFlag, ...
    100*norm(chosen.H(lo,:)-D.validationFRF(lo,:),'fro')/norm(D.validationFRF(lo,:),'fro'), ...
    100*norm(chosen.H(hi,:)-D.validationFRF(hi,:),'fro')/norm(D.validationFRF(hi,:),'fro'));
disp([chosen.fit.frequency_hz,chosen.fit.damping]);
    function [value,gradient]=objective(x)
        for cacheIndex=1:numel(cache)
            if norm(cache(cacheIndex).coordinates_ms-x,inf)<1e-12,value=cache(cacheIndex).cost;gradient=cache(cacheIndex).gradient;return;end
        end
        ph=exp(1i*2*pi*f.*(E*x*1e-3).');corrected=H.*ph;
        ff=frf_srivc_fit(f,corrected,seed.frequency_hz,seed.damping,bounds,W,opts);
        if ~ff.converged
            retry=opts;retry.max_iterations=500;
            ff=frf_srivc_fit(f,corrected,ff.frequency_hz,ff.damping,bounds,W,retry);
        end
        if ~ff.converged,error('Nonstationary inner SRIVC: score %.3g, %s',ff.score,ff.reason);end
        residual=corrected-ff.prediction;value=norm(residual.*sqrt(W),'fro')^2/norm2;
        gg=2*real(sum(conj(residual).*W.*(1i*2*pi*f*1e-3.*corrected),1)).'/norm2;gradient=E.'*gg;
        cache(end+1)=struct('coordinates_ms',x,'fit',ff,'H',ff.prediction./ph,'cost',value,'gradient',gradient);
        if mod(numel(cache),15)==1,fprintf('z=%g delay evaluation %d train %.5f%% g %.3g\n',D.position_mm,numel(cache),100*sqrt(value),norm(gradient,inf));end
        % Save compact progress; complete fits are retained once at completion.
        progress=struct('coordinates_ms',{cache.coordinates_ms},'cost',{cache.cost},'gradient',{cache.gradient});
        save(fullfile(out,['progress_z' D.tag '.mat']),'progress');
    end
end
