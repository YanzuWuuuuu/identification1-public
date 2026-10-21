function result=validation_case(scenario,seed,noiseLevel,N,methods,starts,extraPole,maxIterations)
% A single paired run. Generation truth never enters initialization/fit.
if nargin<3, noiseLevel=.003; end
if nargin<4, N=8; end
if nargin<5, methods={'expanded','physical'}; end
if nargin<6, starts=1; end
if nargin<7, extraPole=0; end
if nargin<8, maxIterations=200; end
folder=fileparts(mfilename('fullpath')); pilot=fullfile(fileparts(folder),'relative_calibration_pilot'); addpath(pilot);
outdir=fullfile(folder,'outputs'); if ~isfolder(outdir), mkdir(outdir); end
loaded=load(fullfile(pilot,'outputs','common_noise0.mat'),'truth','nominal');
truth=loaded.truth; nominal=loaded.nominal;
switch scenario
    case 'common'
        truth.filter_y=25e-6*ones(4,1); truth.filter_u=truth.filter_y;
        truth.ty=zeros(4,1); truth.tu=12e-6*ones(4,1);
    case 'repeated'
        truth.filter_y=[20;20;35;35]*1e-6; truth.filter_u=[20;35;20;35]*1e-6;
        truth.ty=[0;8;17;5]*1e-6; truth.tu=[4;17;10;25]*1e-6;
    case 'near_repeated'
        truth.filter_y=[25;25.05;25.1;25.15]*1e-6; truth.filter_u=[25.15;25.1;25.05;25]*1e-6;
        truth.ty=[0;8;17;5]*1e-6; truth.tu=[4;17;10;25]*1e-6;
    otherwise
        truth.filter_y=[20;31;9;25]*1e-6; truth.filter_u=[17;11;23;35]*1e-6;
        truth.ty=[0;8;17;5]*1e-6; truth.tu=[4;17;10;25]*1e-6;
end
truth.a=truth.filter_y+truth.filter_u.'; truth.b=truth.filter_y*truth.filter_u.';
f=unique([logspace(log10(2),log10(180),70),200:2:800]).';
fv=unique([logspace(log10(2.1),log10(180.1),100),200.7:1.1:799.5]).';
Gt=response(truth,f,extraPole); Gv=response(truth,fv,extraPole);
sigma=noiseLevel*sqrt(abs(Gt).^2+(1e-7)^2);
rng(20261013+seed,'twister'); raw=Gt+sigma.*(randn(16,numel(f),N)+1i*randn(16,numel(f),N))/sqrt(2);
Y=mean(raw,3); estimatedVariance=sum(abs(raw-Y).^2,3)/(N-1);
% This floor depends only on measured repeat variance; no true noise weight.
estimatedVariance=max(estimatedVariance,median(estimatedVariance(:))*1e-8);
sigmaMean=sqrt(estimatedVariance/N);
options=struct('num_modes',4,'frequency_band',[200 800],'sigma',sigmaMean);
tic; [initial,info]=init_from_frf(Y,f,options); initialization_seconds=toc;
result=struct([]);
for mm=1:numel(methods)
    method=methods{mm}; candidates=cell(starts,1);
    for start=1:starts, candidates{start}=fit_from_frf(Y,f,sigmaMean,initial,method,start,maxIterations); end
    costs=cellfun(@(v)v.cost,candidates); [~,best]=min(costs); fit=candidates{best};
    metrics=rcp_metrics(truth,fit.estimate,nominal,fv);
    if extraPole>0
        % Metrics reference the projected four-mode core. Added unmodelled
        % channel response requires separately recomputing actual leakage.
        metrics=extra_metrics(metrics,truth,fit.estimate,fv,extraPole);
    end
    ev=response_est(fit.estimate,fv); metrics.validation_relative_frf_error=norm(ev-Gv,'fro')/norm(Gv,'fro');
    metrics.frequency_relative_error=(abs(fit.estimate.lambda)/(2*pi)-truth.freq(:))./truth.freq(:);
    metrics.damping_relative_error=(-real(fit.estimate.lambda)./abs(fit.estimate.lambda)-.01)/.01;
    metrics.scenario=scenario; metrics.seed=seed; metrics.noise_level=noiseLevel; metrics.independent_records=N;
    metrics.method=method; metrics.extra_output_pole_seconds=extraPole;
    if extraPole<0
        metrics.mismatch.kind='fifth_mode_present_four_modes_fitted';
        metrics.mismatch.additional_true_frequency_hz=738.83;
        metrics.mismatch.fitted_frequency_hz=abs(fit.estimate.lambda)/(2*pi);
        metrics.mismatch.actual_rigid_to_target_four_l2=metrics.baseline.joint.leak_to_flex_amplitude_ratio;
        metrics.mismatch.note='extra_output_pole_seconds=-1 is a scenario marker, not a negative filter time. True object is legacy five modes plus D; fitted model has four modes. Data initialization may select the fifth peak instead of the fourth; additional_true_frequency_hz does not label the omitted fitted mode. Task pass is intentionally undefined.';
    elseif extraPole>0
        metrics.mismatch.kind='additional_output_lowpass';
    end
    diagnostics=rmfield(fit,{'estimate','x','x0','c','fitted','lb','ub'});
    metrics.solver=diagnostics; metrics.initialization_seconds=initialization_seconds;
    metrics.initial_frequency_hz=abs(initial.lambda)/(2*pi); metrics.candidate_costs=costs;
    metrics.task_pass=max(abs(metrics.frequency_relative_error))<.01 && ...
        max([metrics.relative.output_max_error,metrics.relative.input_max_error])<.05 && ...
        metrics.baseline.joint.leak_to_flex_amplitude_ratio<.01 && ...
        all(arrayfun(@(z)abs(z.joint.retention_amplitude_ratio-1)<.1,metrics.mode));
    if extraPole~=0
        metrics.core_task_pass=metrics.task_pass; metrics.task_pass=NaN;
        metrics.task_pass_note='Not defined for misspecified objects; inspect actual mismatch metrics separately.';
    end
    tag=sprintf('%s_seed%02d_noise%.4g_N%d_%s_s%d_extra%.4g',scenario,seed,noiseLevel,N,method,starts,extraPole);
    if maxIterations~=200, tag=[tag '_budget' num2str(maxIterations)]; end
    save(fullfile(outdir,[tag '.mat']),'truth','nominal','f','fv','raw','Y','sigmaMean','initial','info','fit','metrics','candidates');
    compact=rmfield(metrics,'curves'); fid=fopen(fullfile(outdir,[tag '.json']),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(compact,'PrettyPrint',true)); fclose(fid);
    if isempty(result), result=compact; else, result(end+1)=compact; end %#ok<AGROW>
    fprintf('%s: cost %.4g, it %d, dyn %.3g%%, leak/flex %.3g, pass %d\n', ...
        tag,fit.mean_whitened_residual,fit.iterations,100*max(metrics.relative.output_max_error,metrics.relative.input_max_error), ...
        metrics.baseline.joint.leak_to_flex_amplitude_ratio,metrics.task_pass);
end
end

function G=response(p,f,extra)
G=response_est(p,f);
if extra<0
    s=2i*pi*f(:).'; w=2*pi*738.83;
    A=[-.9896;1.0118;1.0857;-.9832]*[-.2543,.1529,-.0934,-.0085];
    ell=p.ty+p.tu.'; H=exp(-ell(:)*s)./(1+p.a(:)*s+p.b(:)*s.^2);
    G=G+H.*(A(:)*(1./(s.^2+2*.01*w*s+w^2)));
end
if extra>0
    s=2i*pi*f(:).'; ex=1./(1+extra*[1;1.4;.7;1.1]*s);
    G=G.*repmat(ex,4,1);
end
end

function G=response_est(p,f)
s=2i*pi*f(:).'; P=p.M(:)*(1./s.^2)+p.D(:)*ones(1,numel(f));
for i=1:4
    R=p.R(:,:,i); P=P+R(:)*(1./(s-p.lambda(i)))+conj(R(:))*(1./(s-conj(p.lambda(i))));
end
ell=p.ty+p.tu.'; G=P.*exp(-ell(:)*s)./(1+p.a(:)*s+p.b(:)*s.^2);
end

function out=extra_metrics(out,t,e,f,extra)
% Additive to the pilot evaluator: evaluate the actual extra output filters.
[U,~,V]=svd(e.M); ny=U(:,4); nu=V(:,4);
dy=e.ty-e.ty(1); du=e.tu-e.tu(1); Ly=max([0;dy])+1e-5; Lu=max([0;du])+1e-5;
rig=zeros(numel(f),1); flex=rig;
for i=1:numel(f)
    s=2i*pi*f(i); d=1+e.a*s+e.b*s^2;
    cy=exp(-s*(Ly-dy)).*d(:,1)/d(1,1); cu=exp(-s*(Lu-du)).*d(1,:).'/d(1,1);
    hy=exp(-s*t.ty)./(1+s*t.filter_y)./(1+extra*[1;1.4;.7;1.1]*s);
    hu=exp(-s*t.tu)./(1+s*t.filter_u); F=zeros(4);
    for j=1:4, F=F+t.R(:,:,j)/(s-t.lambda(j))+conj(t.R(:,:,j))/(s-conj(t.lambda(j))); end
    left=(ny.*cy.*hy).'; right=nu.*cu.*hu;
    rig(i)=abs(left*t.M*right/s^2); flex(i)=abs(left*F*right);
end
out.mismatch.actual_rigid_to_flexible_l2=sqrt(trapz(f,rig.^2)/trapz(f,flex.^2));
out.mismatch.note='Original baseline metrics refer to the core model; use this actual-object ratio for mismatch.';
end
