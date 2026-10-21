function [beta,info]=refine_flexible_fit(beta0,f,G,weight,K,bounds,maxIterations,tolerance)
% Bounded weighted output-error refinement after RIV initialization.
% This is explicitly reported separately from the paper's RIV iterations.
assert(exist('lsqnonlin','file')==2,'Weighted refinement requires Optimization Toolbox.');
noni=size(G,1); beta0=beta0(:); scales=zeros(size(beta0));
lower=-inf(size(beta0)); upper=inf(size(beta0));
for k=1:K
    ix=(k-1)*(noni+2)+(1:noni+2);
    scales(ix(1:2))=max(abs(beta0(ix(1:2))),realmin);
    scales(ix(3:end))=max(norm(beta0(ix(3:end)))/sqrt(noni),eps);
    lower(ix(1))=eps;
    lower(ix(2))=1/(2*pi*bounds(k,2))^2;
    upper(ix(2))=1/(2*pi*bounds(k,1))^2;
end
normalizer=max(norm(apply_sqrt_weight(stack_complex(G),weight)),eps);
x0=min(max(beta0,lower),upper)./scales;
initialResidual=objective(x0);
options=optimoptions('lsqnonlin','Display','off','SpecifyObjectiveGradient',true, ...
    'MaxIterations',maxIterations,'OptimalityTolerance',tolerance, ...
    'StepTolerance',tolerance*0.01,'FunctionTolerance',tolerance*0.01);
[x,cost,~,exitflag,output]=lsqnonlin(@objective,x0,lower./scales,upper./scales,options);
beta=x.*scales;
frequencies=1./(2*pi*sqrt(beta(2:noni+2:end)));
info=struct('method','bounded_weighted_output_error_after_RIV', ...
    'converged',exitflag>0,'exitflag',exitflag,'iterations',output.iterations, ...
    'initial_cost',sum(initialResidual.^2),'final_cost',cost, ...
    'first_order_optimality',output.firstorderopt,'message',output.message, ...
    'frequency_bound_active',abs(frequencies-bounds(:,1))<1e-5 | abs(frequencies-bounds(:,2))<1e-5);

    function [r,J]=objective(x)
        [prediction,derivative]=additive_frf([zeros(noni,1);x.*scales],f,noni,K);
        r=apply_sqrt_weight(stack_complex(prediction-G),weight)/normalizer;
        J=apply_sqrt_weight(derivative(:,noni+1:end).*scales.',weight)/normalizer;
    end
end
