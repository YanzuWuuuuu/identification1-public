function verify_raw_srivc
% Known continuous-frequency 3-input/2-output system, generated independently.
rs=RandStream('mt19937ar','Seed',427);f=linspace(30,650,1200).';s=2i*pi*f;
U=(randn(rs,numel(f),3)+1i*randn(rs,numel(f),3))/sqrt(2);
ft=[180,430];zt=[.015,.007];Y=zeros(numel(f),2);
for m=1:2
    B0=randn(rs,3,2);B1=.05*randn(rs,3,2);
    q=s/(2*pi*ft(m));den=1+2*zt(m)*q+q.^2;
    boundary=randn(rs,1,2)*.03+q.*randn(rs,1,2)*.01;
    Y=Y+(U*B0+q.*(U*B1)+boundary)./den;
end
Y=Y+U*(randn(rs,3,2)*.02)+ones(numel(f),1)*[.01,-.02];
opts=struct('max_iterations',100,'weight',ones(size(f)),'debug',true,'verbose',true);
fit=raw_srivc_fit(f,U,Y,ft.*[1.025,.98],[.03,.02],[150 210;390 470],opts);
fprintf('Known-system recovered Hz: %s; zeta: %s\n',mat2str(fit.frequency_hz.',9),mat2str(fit.damping.',8));
disp(fit.debug);
assert(max(abs(fit.frequency_hz(:)-ft(:)))<1e-4);
assert(max(abs(fit.damping(:)-zt(:)))<1e-6);
assert(fit.relative_error<1e-7);
assert(fit.debug.schur_denominator_relative_error<1e-7);
assert(fit.debug.instrument_derivative_relative_error<1e-6);
assert(fit.debug.initial_regressor_instrument_gap>.01);
save(fullfile(fileparts(mfilename('fullpath')),'synthetic_verification.mat'),'fit','ft','zt');
fprintf('PASS: pole recovery, damping recovery, IV Schur solve, instrument derivative, X distinct from Z.\n');
end
