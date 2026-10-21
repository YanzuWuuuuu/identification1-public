function verify_frf_srivc
f=(25:1:800).';s=1i*2*pi*f;fn=[180 430];damping=[.015 .007];
H=repmat([.001 -.002 .003],numel(f),1);
for k=1:2
    wn=2*pi*fn(k);numerator=wn^2*[1 -.4 .2]*(1+.1*k)+s*wn*[.03 -.02 .04];
    H=H+numerator./(s.^2+2*damping(k)*wn*s+wn^2);
end
W=1./max(abs(H),1e-3*max(abs(H),[],1));
opts=struct('fixed_modes',[],'debug',true,'verbose',false,'max_iterations',80);
fit=frf_srivc_fit(f,H,[175 435],[.02 .01],[160 200;410 450],W,opts);
assert(fit.converged && max(abs(fit.frequency_hz-fn.'))<1e-5);
assert(fit.debug.schur_full_update_relative_error<1e-7);
assert(fit.debug.instrument_derivative_relative_error<1e-6);
save(fullfile(fileparts(mfilename('fullpath')),'frf_synthetic_verification.mat'),'fit');
disp(fit.debug);fprintf('Synthetic FRF recovery: %.9f and %.9f Hz; weighted relative error %.3g\n',fit.frequency_hz,fit.relative_error);
end
