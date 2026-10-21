function check_fixed_pole_weights
% Isolate weights from nonlinear optimizer behavior: hold ALL poles fixed.
out=fileparts(mfilename('fullpath'));
S=load(fullfile(out,'frf_preprocessing_comparison.mat'),'comparison');C=S.comparison;
R=load(fullfile(out,'srivc_z0_mid_fixed_full.mat'),'result');p=R.result.fit;
opts=struct('fixed_modes',1:numel(p.frequency_hz),'verbose',false,'max_iterations',0);
rows={};
for e=1:2
    if e==1,target=C.welchTrain;source='welch';else,target=C.lrmTrain;source='lrm';end
    for w=1:2
        if w==1,W=ones(size(target));weight='uniform';
        else,W=1./max(abs(target),max(1e-3*max(abs(target),[],1),eps));weight='inverse_magnitude';end
        fit=frf_srivc_fit(C.f,target,p.frequency_hz,p.damping,p.bounds_hz,W,opts);
        assert(fit.score<1e-6 && fit.iterations==0,'Fixed-pole linear optimum was not attained.');
        for j=1:2
            ch=1+5*(j-1);band=C.f>=180 & C.f<=260;
            ref=C.welchValidation(band,ch);h=fit.prediction(band,ch);
            rows(end+1,:)={source,weight,j,100*norm(h-ref)/norm(ref),fit.score}; %#ok<AGROW>
        end
    end
end
T=cell2table(rows,'VariableNames',{'frf_source','quadratic_weight','diagonal_channel','notch_band_complex_error_percent','linear_correlation_score'});
writetable(T,fullfile(out,'frf_fixed_pole_weight_check.csv'));disp(T);
end
