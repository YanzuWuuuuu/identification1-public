function [num,rou_est,info] = Gauss_Newton(belta,rho0,u,update_factor,Max_iter,e,no,ni,isubsystem,nrbm,metric)
% Scaled damped Gauss-Newton, objective ||metric*(f(rho)-belta)||^2.
if nargin<10 || isempty(nrbm), nrbm=1; end
belta=belta(:); rho=rho0(:); noni=no*ni;
if nargin<11 || isempty(metric)
    scales=zeros(size(belta));
    scales(1:noni)=max(norm(belta(1:noni))/sqrt(noni),eps);
    for k=1:isubsystem
        ix=noni+(k-1)*(noni+2)+(1:noni+2);
        scales(ix(1:2))=max(abs(belta(ix(1:2))),eps);
        scales(ix(3:end))=max(norm(belta(ix(3:end)))/sqrt(noni),eps);
    end
    metric=diag(1./scales);
end
assert(size(metric,2)==numel(belta),'Projection weight size mismatch.');
mu=max(u,1e-12); factor=max(update_factor,2);
[prediction,J]=f_rho_construct(rho,no,ni,isubsystem,nrbm);
r=metric*(prediction-belta); cost=r.'*r;
history=cost; converged=false; reason='maximum_iterations'; num=0;
for iteration=1:Max_iter
    A=metric*J;
    colnorm=vecnorm(A);
    scale=1./max(colnorm,max(max(colnorm)*1e-12,realmin));
    As=A.*scale;
    if norm(As.'*r,inf)<=e*max(norm(r),realmin)
        converged=true; reason='gradient'; break;
    end
    accepted=false;
    for attempt=1:25
        % QR on the augmented problem avoids normal-equation inverses.
        stepScaled=[As;sqrt(mu)*eye(numel(rho))]\[-r;zeros(numel(rho),1)];
        step=scale.'.*stepScaled;
        candidate=balance_factors(rho+step,no,ni,isubsystem,nrbm);
        poleIndex=nrbm*(no+ni)+1:(no+ni+2):numel(rho);
        stable=all(candidate(poleIndex)<0)&all(candidate(poleIndex+1)>0);
        if stable && all(isfinite(candidate))
            candidatePrediction=f_rho_construct(candidate,no,ni,isubsystem,nrbm);
            candidateResidual=metric*(candidatePrediction-belta);
            candidateCost=candidateResidual.'*candidateResidual;
        else
            candidateCost=inf;
        end
        if candidateCost<cost, accepted=true; break; end
        mu=mu*factor;
    end
    num=iteration;
    if ~accepted, reason='no_descent'; break; end
    relativeStep=norm(stepScaled)/max(norm(rho./scale.'),realmin);
    relativeDecrease=(cost-candidateCost)/max(cost,realmin);
    rho=candidate; r=candidateResidual; cost=candidateCost;
    history(end+1,1)=cost; %#ok<AGROW>
    [prediction,J]=f_rho_construct(rho,no,ni,isubsystem,nrbm);
    mu=max(mu/factor,1e-12);
    if relativeStep<=e || relativeDecrease<=e
        converged=true; reason='step_or_cost'; break;
    end
end
rou_est=rho;
info=struct('converged',converged,'reason',reason,'iterations',num, ...
    'cost_history',history(:),'weighted_residual_norm',sqrt(cost), ...
    'relative_parameter_residual',norm(prediction-belta)/max(norm(belta),eps));
end

function rho=balance_factors(rho,no,ni,K,nrbm)
% Equal norms fix scalar factor ambiguity without changing any residue.
starts=[(0:nrbm-1)*(no+ni),nrbm*(no+ni)+(0:K-1)*(no+ni+2)+2];
for first=starts
    il=first+(1:no); ir=first+no+(1:ni);
    nl=norm(rho(il)); nr=norm(rho(ir));
    if nl>0 && nr>0
        alpha=sqrt(nr/nl); rho(il)=rho(il)*alpha; rho(ir)=rho(ir)/alpha;
    end
end
end
