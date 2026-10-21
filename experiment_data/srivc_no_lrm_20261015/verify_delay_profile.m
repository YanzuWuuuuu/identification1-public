function audit=verify_delay_profile
% Numerical audit of additive delay coordinates and stationary-profile gradient.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'delay_identification_20261016');
S=load(fullfile(out,'common_delay_srivc.mat'),'result');C=S.result;
S=load(fullfile(out,'channel_delay_srivc.mat'),'result');R=S.result;
f=C.f;H=C.trainingFRF;W=C.weights;E=R.delay_design;
normalizer=norm(H.*sqrt(W),'fro')^2;
seed=C.selected.fit;
opts=struct('fixed_modes',1,'max_iterations',400,'score_tolerance',1e-7,'verbose',false);
% A nonstationary outer point avoids comparing two almost-zero gradients.
x=[C.selected.tau_us/1000;.02;-.01;.015;-.015;.01;-.02];
[cost,g,fit,identityError]=profile(x);
step=1e-4;fd=zeros(7,1);scores=zeros(7,2);
for k=1:7
    xp=x;xm=x;xp(k)=xp(k)+step;xm(k)=xm(k)-step;
    [cp,~,fp]=profile(xp);[cm,~,fm]=profile(xm);
    fd(k)=(cp-cm)/(2*step);scores(k,:)=[fp.score fm.score];
end
relativeError=norm(fd-g)/max(norm(g),realmin);
L=R.delay_matrix_seconds;
crossDifference=L-L(:,1)-L(1,:)+L(1,1);
coordinatesError=norm(L(:)-1e-3*E*R.selected.coordinates_ms);
% The physical input/output gauge leaves every end-to-end delay unchanged.
ty=[0;R.selected.coordinates_ms(2:4)];
tu=R.selected.coordinates_ms(1)+[0;R.selected.coordinates_ms(5:7)];
gaugeShift=.031;
gaugeError=norm((ty+gaugeShift)+(tu-gaugeShift).'-(ty+tu.'),'fro');
audit=struct('coordinate_rank',rank(E),'delay_matrix_reconstruction_error',coordinatesError, ...
    'additive_cross_difference_max_seconds',max(abs(crossDifference),[],'all'), ...
    'physical_delay_gauge_error_ms',gaugeError,'test_coordinates_ms',x, ...
    'test_cost',cost,'central_difference_step_ms',step,'analytic_gradient',g, ...
    'finite_difference_gradient',fd,'gradient_relative_error',relativeError, ...
    'phase_compensation_cost_relative_error',identityError,'center_inner_score',fit.score, ...
    'maximum_perturbation_inner_score',max(scores,[],'all'));
assert(rank(E)==7 && coordinatesError<1e-15 && max(abs(crossDifference),[],'all')<1e-15);
assert(gaugeError<1e-12 && identityError<1e-12);
assert(relativeError<1e-3,'Profile gradient disagrees with central differences.');
save(fullfile(out,'delay_profile_audit.mat'),'audit');
writetable(table((1:7).',g,fd,'VariableNames',{'coordinate','analytic_gradient_per_ms','finite_difference_per_ms'}), ...
    fullfile(out,'delay_profile_gradient_audit.csv'));
disp(audit);

    function [value,gradient,inner,identity]=profile(coordinates)
        phase=exp(1i*2*pi*f.*(1e-3*E*coordinates).');
        corrected=H.*phase;
        inner=frf_srivc_fit(f,corrected,seed.frequency_hz,seed.damping,seed.bounds_hz,W,opts);
        assert(inner.converged,'Audit inner fit did not converge: reason=%s, score=%.4g, iterations=%d.', ...
            inner.reason,inner.score,inner.iterations);
        residual=corrected-inner.prediction;
        value=norm(residual.*sqrt(W),'fro')^2/normalizer;
        originalCost=norm((H-inner.prediction./phase).*sqrt(W),'fro')^2/normalizer;
        identity=abs(originalCost-value)/max(value,realmin);
        channelGradient=2*real(sum(conj(residual).*W.*(1i*2*pi*f*1e-3.*corrected),1)).'/normalizer;
        gradient=E.'*channelGradient;
    end
end
