function result=refine_spring_delays
% Seven additive input/output delay coordinates; three spring poles stay free.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'spring_identification_20261016');
S=load(fullfile(out,'spring_balanced_selected.mat'),'result');parent=S.result;
f=parent.f;H=parent.trainingFRF;W=parent.weights;E=parent.delay_design;seed=parent.selected.fit;
norm2=norm(H.*sqrt(W),'fro')^2;start=parent.coordinates_ms;
opts=struct('fixed_modes',[],'max_iterations',250,'score_tolerance',1e-6,'verbose',false);
cache=struct('coordinates_ms',{},'fit',{},'H',{},'cost',{},'gradient',{});
options=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'Display','none', ...
    'MaxIterations',60,'MaxFunctionEvaluations',140,'OptimalityTolerance',1e-7,'StepTolerance',1e-7,'ConstraintTolerance',1e-9);
[x,value,flag,output]=fmincon(@objective,start,[E;-E],[ones(16,1);zeros(16,1)], ...
    [],[],[0;-ones(6,1)],[1;ones(6,1)],[],options);
objective(x);[~,i]=min([cache.cost]);chosen=cache(i);result=parent;result.selected=chosen;
result.delay_matrix_seconds=reshape(E*chosen.coordinates_ms,4,4)*1e-3;result.coordinates_ms=chosen.coordinates_ms;
result.optimizer_coordinates_ms=x;result.optimizer_cost=value;result.optimizer_exitflag=flag;result.optimizer_output=output;
result.selected_matches_optimizer=norm(x-chosen.coordinates_ms,inf)<1e-10;result.trace_cache=cache;
result.method='Three free low-frequency poles and upper-band modal components by SRIVC; seven additive delay coordinates profiled by SQP; inverse-magnitude band-balanced weights';
save(fullfile(out,'spring_delay_refined.mat'),'result','-v7.3');
lo=f<25;hi=~lo;V=parent.validationFRF;
fprintf('Spring delay: exit %d inner %.3g outer %.3g low val %.5f upper val %.5f\n',flag,chosen.fit.score,norm(chosen.gradient,inf), ...
    100*norm(chosen.H(lo,:)-V(lo,:),'fro')/norm(V(lo,:),'fro'),100*norm(chosen.H(hi,:)-V(hi,:),'fro')/norm(V(hi,:),'fro'));
disp([chosen.fit.frequency_hz(1:3),chosen.fit.damping(1:3)]);disp(result.delay_matrix_seconds*1e6);
    function [value,gradient]=objective(x)
        for j=1:numel(cache)
            if norm(cache(j).coordinates_ms-x,inf)<1e-12,value=cache(j).cost;gradient=cache(j).gradient;return;end
        end
        phase=exp(1i*2*pi*f.*(E*x*1e-3).');corrected=H.*phase;
        fit=frf_srivc_fit(f,corrected,seed.frequency_hz,seed.damping,seed.bounds_hz,W,opts);
        if ~fit.converged,error('Spring profile has nonstationary inner fit: %.4g',fit.score);end
        residual=corrected-fit.prediction;value=norm(residual.*sqrt(W),'fro')^2/norm2;
        g=2*real(sum(conj(residual).*W.*(1i*2*pi*f*1e-3.*corrected),1)).'/norm2;gradient=E.'*g;
        cache(end+1)=struct('coordinates_ms',x,'fit',fit,'H',fit.prediction./phase,'cost',value,'gradient',gradient);
        if mod(numel(cache),10)==1,fprintf('Spring delay profile %d: weighted train %.5f%% g %.3g\n',numel(cache),100*sqrt(value),norm(gradient,inf));end
        save(fullfile(out,'spring_delay_checkpoint.mat'),'cache','-v7.3');
    end
end
