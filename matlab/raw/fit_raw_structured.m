function fit=fit_raw_structured(u,y,initialFit,nrbm,maxIterations)
% Direct raw fit with rank-nrbm rigid matrix and rank-one modal residues.
% Rigid input subspace comes only from the unconstrained raw-data fit.
K=numel(initialFit.frequency_hz); ni=size(u,2); no=size(y,2);
ts=initialFit.sample_time; band=initialFit.residual_band_hz;
holdMethod=initialFit.hold_method; bounds=initialFit.frequency_bounds_hz;
[~,~,rigidRight]=svd(initialFit.residues(:,:,1),'econ');
rigidRight=rigidRight(:,1:nrbm); right0=zeros(ni,K);
for k=1:K
    [~,~,V]=svd(initialFit.residues(:,:,k+1),'econ'); right0(:,k)=V(:,1);
end
scale=max(sqrt(mean(y.^2,1)),realmin); target=y./scale;
weightedTarget=raw_residual_transform(target,ts,band);
normalizer=max(norm(weightedTarget,'fro'),realmin);
x0=[log(initialFit.frequency_hz(:));log(initialFit.damping(:));right0(:)];
lower=[log(bounds(:,1));log(repmat(.001,K,1));-inf(ni*K,1)];
upper=[log(bounds(:,2));log(repmat(.15,K,1));inf(ni*K,1)];
options=optimoptions('lsqnonlin','Display','off','MaxIterations',maxIterations, ...
    'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-5, ...
    'OptimalityTolerance',1e-10,'FunctionTolerance',1e-12,'StepTolerance',1e-10);
initialResidual=objective(x0);
[x,cost,~,flag,output]=lsqnonlin(@objective,x0,lower,upper,options);
[~,coeff,prediction,right]=objective(x); coeff=coeff.*scale;
residues=zeros(no,ni,K+1); residues(:,:,1)=coeff(1:nrbm,:).'*rigidRight.';
for k=1:K, residues(:,:,k+1)=coeff(nrbm+k,:).'*right(:,k).'; end
fit=initialFit;
fit.frequency_hz=exp(x(1:K)); fit.damping=exp(x(K+1:2*K));
fit.residues=residues; fit.initial_coefficients=coeff(nrbm+K+1:end,:);
fit.prediction=prediction.*scale; fit.exitflag=flag; fit.output=output;
fit.cost=cost; fit.initial_cost=sum(initialResidual.^2);
fit.method='direct_raw_structured_output_error'; fit.rigid_input_subspace=rigidRight;

    function [residual,coeff,prediction,right]=objective(x)
        [forced,initial]=raw_modal_basis(u,ts,exp(x(1:K)),exp(x(K+1:2*K)),holdMethod);
        right=reshape(x(2*K+1:end),ni,K);
        right=right./max(vecnorm(right),realmin);
        design=zeros(size(u,1),nrbm+K+size(initial,2));
        design(:,1:nrbm)=forced(:,1:ni)*rigidRight;
        for mode=1:K
            design(:,nrbm+mode)=forced(:,mode*ni+(1:ni))*right(:,mode);
        end
        design(:,nrbm+K+1:end)=initial;
        weightedDesign=raw_residual_transform(design,ts,band);
        coeff=scaled_lsq(weightedDesign,weightedTarget);
        prediction=design*coeff;
        residual=(weightedDesign*coeff-weightedTarget)/normalizer; residual=residual(:);
    end
end
