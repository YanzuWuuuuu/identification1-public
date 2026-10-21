function result=fit_channel_delay_srivc
% Seven delay coordinates: ell_jk=t0+dy_j+du_k, dy_1=du_1=0.
% The outer SQP uses a stationary SRIVC branch and its analytic envelope
% gradient; it is not a claim of a globally minimized profile likelihood.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'delay_identification_20261016');
S=load(fullfile(out,'common_delay_srivc.mat'),'result');C=S.result;
f=C.f;H=C.trainingFRF;W=C.weights;normalizer=norm(H.*sqrt(W),'fro')^2;
seed=C.selected.fit;fn=seed.frequency_hz;z=seed.damping;bounds=seed.bounds_hz;
opts=struct('fixed_modes',1,'max_iterations',160,'score_tolerance',1e-6,'verbose',false);
E=ones(16,7);E(:,2:end)=0;
for u=1:4
    for y=1:4
        ch=y+4*(u-1);if y>1,E(ch,y)=1;end;if u>1,E(ch,3+u)=1;end
    end
end
start=[C.selected.tau_us/1000;zeros(6,1)];
cache=struct('coordinates_ms',{},'fit',{},'H',{},'cost',{},'gradient',{});
options=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true, ...
    'Display','iter','MaxIterations',65,'MaxFunctionEvaluations',130, ...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-7,'ConstraintTolerance',1e-9);
[x,cost,exitflag,output]=fmincon(@objective,start,[E;-E],[ones(16,1);zeros(16,1)], ...
    [],[],[-0; -ones(6,1)],[1;ones(6,1)],[],options);
objective(x);valid=find(arrayfun(@(q)q.fit.converged && all(E*q.coordinates_ms>=-1e-8) && all(E*q.coordinates_ms<=1+1e-8),cache));
[~,k]=min([cache(valid).cost]);selected=cache(valid(k));
result=struct('f',f,'selected',selected,'delay_design',E,'delay_matrix_seconds', ...
    reshape(E*selected.coordinates_ms,4,4)*1e-3,'weights',W,'trace_cache',cache, ...
    'optimizer_coordinates_ms',x,'optimizer_cost',cost,'optimizer_exitflag',exitflag,'optimizer_output',output, ...
    'method','Input/output additive effective delays; fixed inverse-magnitude W; inner additive SRIVC; outer constrained profile SQP', ...
    'gauge','ell_jk=t0+dy_j+du_k, dy1=du1=0; 0<=ell_jk<=1 ms; relative dy,du may be negative', ...
    'selection','Lowest training cost among feasible converged inner fits; no validation-based selection');
save(fullfile(out,'channel_delay_srivc.mat'),'result','-v7.3');
trace=zeros(numel(cache),12);
for j=1:numel(cache)
    q=cache(j);trace(j,:)=[q.coordinates_ms.',q.cost,q.fit.converged,q.fit.score,q.fit.iterations,norm(q.gradient,inf)];
end
writetable(array2table(trace,'VariableNames',{'t0_ms','dy2_ms','dy3_ms','dy4_ms','du2_ms','du3_ms','du4_ms', ...
    'normalized_training_cost','inner_converged','inner_score','inner_iterations','outer_gradient_infinity'}),fullfile(out,'channel_delay_profile.csv'));
fprintf('Selected additive delays (microseconds), outer exit=%d:\n',exitflag);disp(1e6*result.delay_matrix_seconds);
fprintf('Training error %.6f%%, inner score %.3g, profile gradient %.3g\n',100*sqrt(selected.cost),selected.fit.score,norm(selected.gradient,inf));

    function [value,gradient]=objective(coordinates)
        for old=1:numel(cache)
            if norm(cache(old).coordinates_ms-coordinates,inf)<1e-12
                value=cache(old).cost;gradient=cache(old).gradient;return;
            end
        end
        ell=(E*coordinates).'*1e-3;phase=exp(1i*2*pi*f.*ell);corrected=H.*phase;
        fit=frf_srivc_fit(f,corrected,fn,z,bounds,W,opts);
        if ~fit.converged
            retry=opts;retry.max_iterations=400;
            fit=frf_srivc_fit(f,corrected,fit.frequency_hz,fit.damping,bounds,W,retry);
        end
        assert(fit.converged,'Inner SRIVC did not reach stationarity; profile gradient is not certified.');
        pred=fit.prediction./phase;error=corrected-fit.prediction;
        value=norm((H-pred).*sqrt(W),'fro')^2/normalizer;
        channelGradient=2*real(sum(conj(error).*W.*(1i*2*pi*f*1e-3.*corrected),1)).'/normalizer;
        gradient=E.'*channelGradient;
        cache(end+1)=struct('coordinates_ms',coordinates,'fit',fit,'H',pred,'cost',value,'gradient',gradient);
        fprintf('Profile %d: train %.5f%%, inner %.2g, g %.2g\n',numel(cache),100*sqrt(value),fit.score,norm(gradient,inf));
        save(fullfile(out,'channel_delay_checkpoint.mat'),'cache','E','start','-v7.3');
    end
end
