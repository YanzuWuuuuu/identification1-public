function fit=fit_raw_modal(u,y,ts,frequency0,damping0,bounds,holdMethod,maxIterations,band)
% Direct raw time-domain output-error fit, with variable projection.
% Residue matrices and initial-condition responses are linear nuisance terms.
% Only frequency/damping are nonlinear. No LRM, FFT ratios, or true parameters.
assert(size(u,1)==size(y,1)&&all(isfinite(y(:)))&&isreal(y));
if nargin<9, band=[]; end
assert(isempty(band)||(numel(band)==2 && band(1)==0 && band(2)>0 && band(2)<.5/ts), ...
    'Keep DC in a fitted residual band so initial offset and velocity remain identifiable.');
K=numel(frequency0); ni=size(u,2); no=size(y,2);
assert(isequal(size(bounds),[K,2]));
scale=max(sqrt(mean(y.^2,1)),realmin);
target=y./scale;
weightedTarget=raw_residual_transform(target,ts,band);
normalizer=max(norm(weightedTarget,'fro'),realmin);
x0=log([frequency0(:);damping0(:)]);
lower=log([bounds(:,1);repmat(.001,K,1)]);
upper=log([bounds(:,2);repmat(.15,K,1)]);
options=optimoptions('lsqnonlin','Display','off','MaxIterations',maxIterations, ...
    'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-5, ...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1e-12,'StepTolerance',1e-10);
initialResidual=objective(x0);
[x,cost,~,flag,output]=lsqnonlin(@objective,x0,lower,upper,options);
[~,coeff,prediction]=objective(x);
coeff=coeff.*scale;
fit=struct('frequency_hz',exp(x(1:K)),'damping',exp(x(K+1:end)), ...
    'residues',reshape(coeff(1:ni*(K+1),:).',no,ni,K+1), ...
    'initial_coefficients',coeff(ni*(K+1)+1:end,:), ...
    'prediction',prediction.*scale,'hold_method',holdMethod,'sample_time',ts, ...
    'exitflag',flag,'output',output,'initial_cost',sum(initialResidual.^2), ...
    'cost',cost,'frequency_bounds_hz',bounds,'training_samples',size(u,1), ...
    'method','direct_raw_time_domain_variable_projection','residual_band_hz',band);

    function [residual,coeff,prediction]=objective(x)
        [forced,initial]=raw_modal_basis(u,ts,exp(x(1:K)),exp(x(K+1:end)),holdMethod);
        design=[forced,initial];
        weightedDesign=raw_residual_transform(design,ts,band);
        coeff=scaled_lsq(weightedDesign,weightedTarget);
        prediction=design*coeff;
        residual=(weightedDesign*coeff-weightedTarget)/normalizer; residual=residual(:);
    end
end
