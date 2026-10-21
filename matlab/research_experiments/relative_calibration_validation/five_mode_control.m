function results=five_mode_control(maxIterations)
%FIVE_MODE_CONTROL Fit five modes to the identical saved full-five-mode data.
% Truth is used only after fitting for evaluation; initial parameters are
% generated from Y/f/sigmaMean by init_five_from_frf. First four mechanical
% modes remain the target for leakage and retention metrics.
if nargin<1, maxIterations=200; end
folder=fileparts(mfilename('fullpath'));
pilot=fullfile(fileparts(folder),'relative_calibration_pilot'); addpath(pilot);
outdir=fullfile(folder,'outputs');
sourceName='distinct_seed01_noise0.003_N8_expanded_s1_extra-1.mat';
sourcePath=fullfile(outdir,sourceName);
source=load(sourcePath,'Y','f','sigmaMean','fv','truth','nominal','fit','metrics');
Y=source.Y; f=source.f; sigmaMean=source.sigmaMean; fv=source.fv;
options=struct('num_modes',5,'frequency_band',[200 800],'sigma',sigmaMean);
tic; [initial,info]=init_five_from_frf(Y,f,options); initialization_seconds=toc;
fprintf('Five-mode data initialization:'); fprintf(' %.4f',abs(initial.lambda)/(2*pi)); fprintf(' Hz\n');
results=struct([]);
for methodCell={'expanded','physical'}
    method=methodCell{1};
    fit=fit_general_from_frf(Y,f,sigmaMean,initial,method,1,maxIterations);
    % Evaluation starts here. The added fifth pole never enters init or fit.
    truth=source.truth; nominal=source.nominal;
    fullTruth=append_fifth_mode(truth);
    metrics=rcp_metrics(truth,fit.estimate,nominal,fv);
    metrics=rmfield(metrics,'fit');
    fullResponse=response(fullTruth,fv); estimatedResponse=response(fit.estimate,fv);
    E=estimatedResponse-fullResponse;
    metrics.validation_relative_frf_error=norm(E,'fro')/norm(fullResponse,'fro');
    metrics.full_five_mode_frf_l2_error=sqrt(trapz(fv,sum(abs(E).^2,1).')/ ...
        trapz(fv,sum(abs(fullResponse).^2,1).'));
    referenceFrequency=abs(fullTruth.lambda(:))/(2*pi);
    referenceDamping=-real(fullTruth.lambda(:))./abs(fullTruth.lambda(:));
    metrics.frequency_relative_error=(abs(fit.estimate.lambda(:))/(2*pi)-referenceFrequency)./referenceFrequency;
    metrics.damping_relative_error=(-real(fit.estimate.lambda(:))./abs(fit.estimate.lambda(:))-referenceDamping)./referenceDamping;
    metrics.reference_frequencies_hz=referenceFrequency;
    metrics.estimated_frequencies_hz=abs(fit.estimate.lambda(:))/(2*pi);
    metrics.estimated_damping=-real(fit.estimate.lambda(:))./abs(fit.estimate.lambda(:));
    metrics.scenario='full_five_mode_same_data_control'; metrics.method=method;
    metrics.seed=source.metrics.seed; metrics.noise_level=source.metrics.noise_level;
    metrics.independent_records=source.metrics.independent_records;
    metrics.fitted_mode_count=5; metrics.target_mode_indices=1:4;
    metrics.source_data_file=sourcePath;
    metrics.initial_frequency_hz=abs(initial.lambda)/(2*pi);
    metrics.initialization_seconds=initialization_seconds;
    metrics.solver=rmfield(fit,{'estimate','x','x0','c','fitted','lb','ub'});
    metrics.task_pass=max(abs(metrics.frequency_relative_error))<.01 && ...
        max(metrics.relative.output_max_error,metrics.relative.input_max_error)<.05 && ...
        metrics.baseline.joint.leak_to_flex_amplitude_ratio<.01 && ...
        all(arrayfun(@(z)abs(z.joint.retention_amplitude_ratio-1)<.1,metrics.mode));
    metrics.task_pass_definition='All five frequency errors <1%, composite dynamics <5%, rigid/target-first-four L2 <1%, first-four residue amplitude deviations <10%; solver status is separate.';
    metrics.scope='Identical noisy full-five-mode data and estimated weights as the saved four-mode experiment. The fifth mode is fitted explicitly. Calibration targets remain the first four modes; D and mode five are excluded from the target-flexibility metric.';
    metrics.four_mode_comparison=four_mode_comparison(source,outdir,method,Y,f,sigmaMean);
    tag=['five_mode_control_' method];
    if maxIterations~=200, tag=[tag '_budget' num2str(maxIterations)]; end
    save(fullfile(outdir,[tag '.mat']),'Y','f','fv','sigmaMean','initial','info', ...
        'fit','metrics','truth','fullTruth','nominal','sourcePath');
    compact=rmfield(metrics,'curves');
    fid=fopen(fullfile(outdir,[tag '.json']),'w','n','UTF-8'); assert(fid>=0);
    fprintf(fid,'%s',jsonencode(compact,'PrettyPrint',true)); fclose(fid);
    if isempty(results), results=compact; else, results(end+1)=compact; end %#ok<AGROW>
    fprintf('%s: cost %.4g, flag %d, it %d, dynamic error %.3g%%, first-four leak/flex %.3g, task pass %d\n', ...
        tag,fit.mean_whitened_residual,fit.exitflag,fit.iterations, ...
        100*max(metrics.relative.output_max_error,metrics.relative.input_max_error), ...
        metrics.baseline.joint.leak_to_flex_amplitude_ratio,metrics.task_pass);
end
end

function full=append_fifth_mode(core)
% Matches exactly the extra<0 evaluation generator in validation_case.m.
full=core; w=2*pi*738.83; zeta=.01;
lambda=-zeta*w+1i*w*sqrt(1-zeta^2);
numerator=[-.9896;1.0118;1.0857;-.9832]*[-.2543,.1529,-.0934,-.0085];
full.lambda=[core.lambda(:);lambda];
full.R(:,:,5)=numerator/(lambda-conj(lambda));
if isfield(core,'freq'), full.freq=[core.freq(:);738.83]; end
end

function G=response(p,f)
s=2i*pi*f(:).'; P=p.M(:)*(1./s.^2)+p.D(:)*ones(1,numel(f));
for i=1:numel(p.lambda)
    R=p.R(:,:,i); P=P+R(:)*(1./(s-p.lambda(i)))+conj(R(:))*(1./(s-conj(p.lambda(i))));
end
ell=p.ty+p.tu.'; G=P.*exp(-ell(:)*s)./(1+p.a(:)*s+p.b(:)*s.^2);
end

function comparison=four_mode_comparison(source,outdir,method,Y,f,sigmaMean)
name=sprintf('distinct_seed01_noise0.003_N8_%s_s1_extra-1.mat',method);
comparison.file=fullfile(outdir,name); comparison.available=false;
if strcmp(method,'expanded')
    ref=source;
elseif isfile(comparison.file)
    ref=load(comparison.file,'Y','f','sigmaMean','fit','metrics');
else
    comparison.note='Same-method four-mode fit is not saved; no paired numerical comparison.';
    return;
end
assert(isequaln(ref.Y,Y) && isequaln(ref.f,f) && isequaln(ref.sigmaMean,sigmaMean), ...
    'Four- and five-mode comparisons must use identical data and weights.');
comparison.available=true; comparison.identical_data_and_weights=true;
comparison.fitted_mode_count=4;
comparison.mean_whitened_residual=ref.fit.mean_whitened_residual;
comparison.exitflag=ref.fit.exitflag; comparison.iterations=ref.fit.iterations;
comparison.relative_dynamic_max_error=max(ref.metrics.relative.output_max_error,ref.metrics.relative.input_max_error);
comparison.actual_rigid_to_target_four_l2=ref.metrics.baseline.joint.leak_to_flex_amplitude_ratio;
comparison.validation_relative_frf_error=ref.metrics.validation_relative_frf_error;
comparison.note='Four-mode fitting of the same five-mode object is misspecified; no matched-model task pass is claimed.';
end
