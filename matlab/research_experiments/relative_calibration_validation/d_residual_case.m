function results=d_residual_case(scenario,seed,maxIterations)
%D_RESIDUAL_CASE Five actual modes, optionally three synthetic high modes.
% There is no true constant feedthrough. D is only a fitted residual term.
if nargin<2, seed=1; end
if nargin<3, maxIterations=400; end
folder=fileparts(mfilename('fullpath')); outdir=fullfile(folder,'d_residual_outputs');
if ~isfolder(outdir), mkdir(outdir); end
data=load(fullfile(folder,'outputs','five_mode_control_expanded.mat'),'fullTruth','nominal','f','fv');
core=data.fullTruth; core.D=zeros(4); nominal=data.nominal;
assert(max(abs(abs(core.lambda(:))/(2*pi)-[286.43;344.58;488.27;705.67;738.83]))<1e-7);
switch scenario
    case 'original', extraFrequency=[];
    case 'near', extraFrequency=[950 1200 1500];
    case 'mid', extraFrequency=[1600 2000 2500];
    case 'far', extraFrequency=[3200 4000 5000];
    otherwise, error('Unknown scenario.');
end
% These are synthetic residues; they are not extracted FEA modes.
amplitudes=zeros(5,1);
for k=1:5
    amplitudes(k)=norm(real(core.R(:,:,k)*(core.lambda(k)-conj(core.lambda(k)))),'fro');
end
amplitude=median(amplitudes); A=zeros(4,4,3);
rng(6142026,'twister');
for k=1:3
    left=randn(4,1); right=randn(4,1);
    A(:,:,k)=amplitude*(left/norm(left))*(right/norm(right)).';
end
truth=core;
for k=1:numel(extraFrequency)
    w=2*pi*extraFrequency(k); lambda=-.01*w+1i*w*sqrt(1-.01^2);
    truth.lambda(5+k,1)=lambda;
    truth.R(:,:,5+k)=A(:,:,k)/(lambda-conj(lambda));
end
f=data.f; fv=data.fv; G=response(truth,f); N=8; noiseLevel=.003;
noiseStd=noiseLevel*sqrt(abs(G).^2+1e-14);
rng(20261014+seed,'twister');
raw=G+noiseStd.*(randn(16,numel(f),N)+1i*randn(16,numel(f),N))/sqrt(2);
Y=mean(raw,3); variance=sum(abs(raw-Y).^2,3)/(N-1);
variance=max(variance,median(variance(:))*1e-8); sigmaMean=sqrt(variance/N);
options=struct('num_modes',5,'frequency_band',[200 800],'sigma',sigmaMean);
tic; [initial,info]=init_five_from_frf(Y,f,options); initialization_seconds=toc;
diagnostic=tail_diagnostic(core,truth,fv);
results=struct([]);
for includeD=[false true]
    fit=fit_residual_d_from_frf(Y,f,sigmaMean,initial,'expanded',includeD,maxIterations);
    metrics=d_residual_metrics(core,fit.estimate,truth,nominal,fv);
    metrics.scenario=scenario; metrics.seed=seed; metrics.include_D=includeD;
    metrics.true_constant_D_is_zero=true; metrics.true_modes=numel(truth.lambda);
    metrics.fitted_modes=5; metrics.synthetic_extra_frequencies_hz=extraFrequency;
    metrics.extra_numerator_frobenius_norm=amplitude;
    metrics.noise_level=noiseLevel; metrics.independent_records=N;
    metrics.initial_frequency_hz=abs(initial.lambda)/(2*pi);
    metrics.estimated_frequency_hz=abs(fit.estimate.lambda)/(2*pi);
    metrics.initialization_seconds=initialization_seconds;
    metrics.tail_diagnostic=diagnostic;
    drop=intersect(fieldnames(fit),{'estimate','x','x0','c','fitted','lb','ub','active_parameter_indices','x_free','x0_free','x_expanded','lb_free','ub_free','gradient_free'});
    metrics.solver=rmfield(fit,drop);
    metrics.estimated_D_frobenius_norm=norm(fit.estimate.D,'fro');
    metrics.Dfit_minus_D0_frobenius_norm=norm(fit.estimate.D-diagnostic.D0,'fro');
    metrics.question='Does a fitted constant approximate three actual unmodelled high-frequency modes?';
    metrics.intervention='Same noisy data, estimated weights, five-mode initializer and expanded channel model; only 16 constant-D parameters are removed or fitted.';
    metrics.comparability='Paired within scenario and seed; higher-frequency scenarios share numerator matrices but tails become both flatter and weaker.';
    metrics.evidence_scope='Synthetic sensitivity experiment; added modes are not claimed as actual platform or FEA modes. Frequencies below 800 Hz only are fitted.';
    tag=sprintf('%s_seed%02d_D%d',scenario,seed,includeD);
    if maxIterations~=400, tag=[tag '_budget' num2str(maxIterations)]; end
    save(fullfile(outdir,[tag '.mat']),'core','truth','nominal','A','f','fv','raw','Y','sigmaMean','initial','info','fit','metrics');
    compact=metrics; if isfield(compact,'curves'), compact=rmfield(compact,'curves'); end
    fid=fopen(fullfile(outdir,[tag '.json']),'w','n','UTF-8'); assert(fid>=0);
    fprintf(fid,'%s',jsonencode(compact,'PrettyPrint',true)); fclose(fid);
    if isempty(results), results=compact; else, results(end+1)=compact; end %#ok<AGROW>
    fprintf('%s: cost %.5g, frequency %.4g%%, damping %.4g%%, dynamic %.4g%%, leak %.3g, flag %d, it %d\n', ...
        tag,fit.mean_whitened_residual,100*max(abs(metrics.frequency_relative_error)), ...
        100*max(abs(metrics.damping_relative_error)),100*max(metrics.relative.output_max_error,metrics.relative.input_max_error), ...
        metrics.baseline.joint.leak_to_flex_amplitude_ratio,fit.exitflag,fit.iterations);
end
end

function G=response(p,f)
s=2i*pi*f(:).'; P=mechanical(p,s);
ell=p.ty+p.tu.'; G=P.*exp(-ell(:)*s)./(1+p.a(:)*s+p.b(:)*s.^2);
end

function P=mechanical(p,s)
P=p.M(:)*(1./s.^2)+p.D(:)*ones(size(s));
for k=1:numel(p.lambda)
    R=p.R(:,:,k); P=P+R(:)*(1./(s-p.lambda(k)))+conj(R(:))*(1./(s-conj(p.lambda(k))));
end
end

function d=tail_diagnostic(core,truth,f)
f=f(:); sel=f>=200; fb=f(sel); s=2i*pi*fb.';
tail=zeros(16,numel(fb)); D0=zeros(4);
for k=6:numel(truth.lambda)
    R=truth.R(:,:,k); lambda=truth.lambda(k);
    tail=tail+R(:)*(1./(s-lambda))+conj(R(:))*(1./(s-conj(lambda)));
    D0=D0+real(R*(lambda-conj(lambda)))/abs(lambda)^2;
end
F=mechanical(core,s)-core.M(:)*(1./s.^2);
integral=@(X)trapz(fb,sum(abs(X).^2,1).');
E=integral(tail); EF=integral(F);
best=real(trapz(fb,tail,2)/(fb(end)-fb(1)));
d.D0=D0; d.best_real_constant=reshape(best,4,4);
d.frequency_band_hz=[fb(1) fb(end)];
d.tail_to_target_flexible_l2=sqrt(E/EF);
d.D0_sum_individual_norms=0;
for k=6:numel(truth.lambda)
    d.D0_sum_individual_norms=d.D0_sum_individual_norms+norm(truth.R(:,:,k)*(truth.lambda(k)-conj(truth.lambda(k))),'fro')/abs(truth.lambda(k))^2;
end
if E>0
    d.D0_relative_tail_l2=sqrt(integral(tail-D0(:))/E);
    d.best_real_constant_relative_tail_l2=sqrt(integral(tail-best)/E);
    d.tail_relative_error_at_upper_edge=norm(tail(:,end)-D0(:))/norm(tail(:,end));
else
    d.D0_relative_tail_l2=NaN; d.best_real_constant_relative_tail_l2=NaN;
    d.tail_relative_error_at_upper_edge=NaN;
end
d.note='Deterministic mechanical tail only, without channel filters. Best constant is an unweighted real L2 projection over the stated flexible band, used only for evaluation and never initialization.';
end
