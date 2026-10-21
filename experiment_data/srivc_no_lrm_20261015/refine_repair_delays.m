function result=refine_repair_delays(inputFile,outputStem,maxIterations)
% Refine seven additive delays, preserving the parent model's fixed poles.
if nargin<3,maxIterations=65;end
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
S=load(inputFile,'result');parent=S.result;
S=load(fullfile(base,'delay_identification_20261016/common_delay_srivc.mat'),'result');C=S.result;
f=C.f;H=C.trainingFRF;W=parent.weights;E=parent.delay_design;
seed=parent.selected.fit;norm2=norm(H.*sqrt(W),'fro')^2;
opts=struct('fixed_modes',seed.fixed_modes,'max_iterations',220,'score_tolerance',1e-6,'verbose',false);
cache=struct('coordinates_ms',{},'fit',{},'H',{},'cost',{},'gradient',{});
start=parent.selected.coordinates_ms;
options=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'Display','none', ...
    'MaxIterations',maxIterations,'MaxFunctionEvaluations',2*maxIterations+10, ...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-7,'ConstraintTolerance',1e-9);
[x,value,flag,output]=fmincon(@objective,start,[E;-E],[ones(16,1);zeros(16,1)], ...
    [],[],[0;-ones(6,1)],[1;ones(6,1)],[],options);
objective(x);valid=find(arrayfun(@(q)q.fit.converged && all(E*q.coordinates_ms>=-1e-9) ...
    && all(E*q.coordinates_ms<=1+1e-9),cache));
assert(~isempty(valid),'No feasible converged profile candidate.');
[~,b]=min([cache(valid).cost]);chosen=cache(valid(b));
result=struct('f',f,'selected',chosen,'weights',W,'delay_design',E, ...
    'delay_matrix_seconds',reshape(E*chosen.coordinates_ms,4,4)*1e-3, ...
    'label',outputStem,'optimizer_coordinates_ms',x,'optimizer_cost',value, ...
    'optimizer_exitflag',flag,'optimizer_output',output,'trace_cache',cache, ...
    'parent_file',inputFile,'selection','Training objective only. Reused validation is exploratory.', ...
    'method','Parent fixed-pole constraints retained; inner additive SRIVC; seven additive-delay outer SQP', ...
    'selected_matches_optimizer',norm(chosen.coordinates_ms-x,inf)<1e-10);
save(fullfile(out,[outputStem '.mat']),'result','-v7.3');
fprintf('%s: outer %d; training %.6f%%; Welch %.6f%%; LRM %.6f%%; inner %.3g; outer g %.3g\n', ...
    outputStem,flag,100*sqrt(chosen.cost),100*norm(chosen.H-C.welchValidation,'fro')/norm(C.welchValidation,'fro'), ...
    100*norm(chosen.H-C.lrmValidation,'fro')/norm(C.lrmValidation,'fro'),chosen.fit.score,norm(chosen.gradient,inf));
disp(result.delay_matrix_seconds*1e6);

    function [value,gradient]=objective(coordinates)
        for old=1:numel(cache)
            if norm(cache(old).coordinates_ms-coordinates,inf)<1e-12
                value=cache(old).cost;gradient=cache(old).gradient;return;
            end
        end
        phase=exp(1i*2*pi*f.*(E*coordinates*1e-3).');corrected=H.*phase;
        fit=frf_srivc_fit(f,corrected,seed.frequency_hz,seed.damping,seed.bounds_hz,W,opts);
        if ~fit.converged
            retry=opts;retry.max_iterations=450;
            fit=frf_srivc_fit(f,corrected,fit.frequency_hz,fit.damping,seed.bounds_hz,W,retry);
        end
        if ~fit.converged
            fprintf('Nonstationary trial: %s, score %.3g\n',fit.reason,fit.score);
            error('Repair inner profile did not converge; no envelope gradient used.');
        end
        prediction=fit.prediction./phase;residual=corrected-fit.prediction;
        value=norm(residual.*sqrt(W),'fro')^2/norm2;
        g=2*real(sum(conj(residual).*W.*(1i*2*pi*f*1e-3.*corrected),1)).'/norm2;
        gradient=E.'*g;
        cache(end+1)=struct('coordinates_ms',coordinates,'fit',fit,'H',prediction,'cost',value,'gradient',gradient);
        if mod(numel(cache),10)==1,fprintf('%s profile %d, train %.5f%%, g %.3g\n',outputStem,numel(cache),100*sqrt(value),norm(gradient,inf));end
        save(fullfile(out,[outputStem '_checkpoint.mat']),'cache','start','-v7.3');
    end
end
