function calib=run_coupled_support_calibration(positionIndex)
% Nested support-model diagnostic: add three symmetric damping couplings.
% No raw data, previous model, or controller is changed. The original C
% anchors are retained, so this is a restricted nonproportional extension,
% not a parameterization of every mass-normalized three-DOF mechanical model.
% SRIVC initializes the predecessor; this stage is profiled output-error SQP.
if nargin<1,positionIndex=3;end
base=fileparts(mfilename('fullpath'));
out=fullfile(base,'coupled_support_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'position_spring_identification_20261016','position_campaign.mat'),'campaign');
source=S.campaign.results{positionIndex};
assert(source.position_mm~=-1,'Exclude the anomalous z=-1 record.');
tag=positionTag(source.position_mm);
O=load(fullfile(base,'structured_calibration_20261016',['calibration_' tag '.mat']),'calib');
old=O.calib;oldFit=old.fits.joint;
assert(numel(oldFit.x)==34 && isequal(size(old.high_basis,2),25));
f=old.f(:);H=old.H_training;Hvalidation=old.H_validation;
W=old.initial_free.weights;E=old.initial_free.delay_design;
assert(isequal(size(E),[16 7]) && isequal(size(W),size(H)));
assert(all(isfinite(H(:))) && all(isfinite(W(:))) && all(W(:)>0));
nf=numel(f);sqrtW=sqrt(W);normalizer=norm(H.*sqrtW,'fro');
highBasis=old.high_basis;anchors=old.anchors(:);gains=old.amplitude_scales(:);
nominal=old.nominal;freeRows=cell(3,1);
for m=1:3,freeRows{m}=setdiff(1:4,anchors(m));end
qrHigh=cell(16,1);highConditions=zeros(16,1);
for channel=1:16
    bw=highBasis.*sqrtW(:,channel);br=[real(bw);imag(bw)];scale=max(vecnorm(br),realmin);
    [qq,rr]=qr(br./scale,0);
    qrHigh{channel}=struct('Q',qq,'R',rr,'scale',scale);highConditions(channel)=cond(rr);
end
% First 27 parameters match the predecessor exactly. The additional rho
% coordinates precede its unchanged seven separable delay coordinates.
x0=[oldFit.x(1:27);zeros(3,1);oldFit.x(28:34)];
lb=[old.bounds.lower(1:27);-.999*ones(3,1);old.bounds.lower(28:34)];
ub=[old.bounds.upper(1:27); .999*ones(3,1);old.bounds.upper(28:34)];
linearA=[zeros(16,30),E;zeros(16,30),-E];linearB=[ones(16,1);zeros(16,1)];
positiveMargin=1e-5;
calib=struct('position_mm',source.position_mm,'tag',tag,'f',f,'H_training',H, ...
    'H_validation',Hvalidation,'nominal',nominal,'initial_free',old.initial_free, ...
    'source_data',old.source_data,'source_geometry',old.source_geometry, ...
    'source_predecessor',fullfile(base,'structured_calibration_20261016',['calibration_' tag '.mat']), ...
    'method','Coupled-support profiled output-error SQP, initialized from the preceding SRIVC-based structured fit', ...
    'model_scope','C*(s^2*I+s*D+K)^(-1)*K*B; diagonal K, positive symmetric D; three retained C anchors restrict the nonproportional family', ...
    'weight_contract',old.weight_contract,'validation_scope',old.validation_scope, ...
    'anchors',anchors,'amplitude_scales',gains,'high_basis',highBasis, ...
    'high_scaled_QR_condition',highConditions,'positive_margin',positiveMargin, ...
    'output_dir',out,'parameter_layout','1:27 old low coordinates; 28:30 rho12,rho13,rho23; 31:37 old delay coordinates (ms)');

% An exact nesting check separates a model change from a changed data cost.
[v0,~,rhoZeroFit]=objective(x0);
nesting.low_relative_error=norm(rhoZeroFit.H_low_undelayed-oldFit.H_low_undelayed,'fro')/max(norm(oldFit.H_low_undelayed,'fro'),realmin);
nesting.full_relative_error=norm(rhoZeroFit.H-oldFit.H,'fro')/max(norm(oldFit.H,'fro'),realmin);
nesting.objective_absolute_error=abs(v0-oldFit.train_objective);
nesting.previous_objective=oldFit.train_objective;nesting.new_rho_zero_objective=v0;
assert(nesting.low_relative_error<1e-9 && nesting.full_relative_error<1e-9, ...
    'rho=0 must reproduce the predecessor FRF before any optimization.');
assert(nesting.objective_absolute_error<1e-10+1e-8*abs(oldFit.train_objective), ...
    'rho=0 changed the profiled training objective.');
calib.rho_zero_equivalence=nesting;

% Check away from rho=0, where damping cross derivatives are nontrivial.
test=x0;test(1:27)=test(1:27)+.003*sin((1:27).');test(28:30)=[.12;-.08;.10];
[~,ga]=objective(test);gn=zeros(37,1);
for coordinate=1:37
    step=2e-6*max(1,abs(test(coordinate)));xp=test;xm=test;
    xp(coordinate)=xp(coordinate)+step;xm(coordinate)=xm(coordinate)-step;
    gn(coordinate)=(objective(xp)-objective(xm))/(2*step);
end
groups={setdiff(1:27,[8 9 17 18 26 27]),[8 9 17 18 26 27],28:30,31:37};
groupErrors=zeros(1,4);
for group=1:4
    indices=groups{group};groupErrors(group)=norm(ga(indices)-gn(indices))/max(norm(gn(indices)),1e-10);
end
[cc,~,gc]=dampingConstraint(test);gcn=zeros(37,numel(cc));
for coordinate=28:30
    step=2e-6;xp=test;xm=test;xp(coordinate)=xp(coordinate)+step;xm(coordinate)=xm(coordinate)-step;
    gcn(coordinate,:)=(dampingConstraint(xp)-dampingConstraint(xm)).'/(2*step);
end
checks=struct('relative_profile_gradient',norm(ga-gn)/max(norm(gn),1e-12), ...
    'group_relative_errors',groupErrors,'group_names',{{'spatial','frequency_and_damping','rho','delay'}}, ...
    'constraint_gradient',norm(gc-gcn,'fro')/max(norm(gcn,'fro'),1e-12), ...
    'test_rho',test(28:30));
assert(checks.relative_profile_gradient<2e-5 && all(groupErrors<5e-4) && checks.constraint_gradient<1e-7, ...
    'Coupled-support derivative check failed: profile %.3g, groups %s, constraint %.3g.', ...
    checks.relative_profile_gradient,mat2str(groupErrors,4),checks.constraint_gradient);
calib.gradient_check=checks;
fprintf('z=%g coupled nesting %.3g, gradient %.3g, rho gradient %.3g, initial cost %.6g\n', ...
    source.position_mm,nesting.full_relative_error,checks.relative_profile_gradient,groupErrors(3),v0);

logEntries={};stageName='fixed_delay';stageTimer=tic;
options=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true,'Display','none','MaxIterations',350, ...
    'MaxFunctionEvaluations',1800,'OptimalityTolerance',1e-7, ...
    'StepTolerance',1e-9,'ConstraintTolerance',1e-9,'OutputFcn',@progress);
[xf,~,flag,info]=fmincon(@fixedObjective,x0(1:30),[],[],[],[],lb(1:30),ub(1:30),@dampingConstraint,options);
xFixed=[xf;x0(31:37)];calib.fits.fixed=assemble(xFixed,flag,info);
save(fullfile(out,['calibration_' tag '_progress.mat']),'calib','logEntries','-v7.3');
fprintf('z=%g coupled fixed: weighted %.5f%% exit %d KKT %.3g elapsed %.1fs\n', ...
    source.position_mm,calib.fits.fixed.train_weighted_error_pct,flag,info.firstorderopt,toc(stageTimer));

% One prescribed start only: a bounded model diagnostic, not a global search.
stageName='joint_delay';stageTimer=tic;
[xx,~,flag,info]=fmincon(@objective,xFixed,linearA,linearB,[],[],lb,ub,@dampingConstraint,options);
calib.fits.joint=assemble(xx,flag,info);calib.optimization_log=logEntries;
calib.bounds=struct('lower',lb,'upper',ub,'linearA',linearA,'linearB',linearB);
calib.accepted=calib.fits.joint.exitflag>0 && calib.fits.joint.feasible ...
    && calib.fits.joint.spatial_rank_valid && calib.fits.joint.stable ...
    && calib.fits.joint.kkt_firstorderopt<=1e-6;
if calib.accepted
    calib.selection_status='Single prescribed coupled-model start passed feasibility, stability, rank and KKT checks';
else
    calib.selection_status='Single coupled-model diagnostic retained with unmet acceptance checks';
end
calib.predecessor_summary=struct('objective',oldFit.train_objective,'Ny',oldFit.Ny,'Nu',oldFit.Nu, ...
    'frequency_hz',oldFit.frequency_hz,'damping',oldFit.damping,'coordinates_ms',oldFit.coordinates_ms);
calib.saved_at=char(datetime('now'));
save(fullfile(out,['calibration_' tag '.mat']),'calib','-v7.3');
fprintf('SAVED coupled %s accepted=%d weighted %.5f%% KKT %.3g rho [%g %g %g]\n', ...
    tag,calib.accepted,calib.fits.joint.train_weighted_error_pct,info.firstorderopt,xx(28:30));

    function [value,g]=fixedObjective(x)
        [value,allGradient]=objective([x;x0(31:37)]);g=allGradient(1:30);
    end

    function [value,g,detail]=objective(x)
        wantGradient=nargout>=2;
        [low,meta,jacobian]=lowModel(x,wantGradient);
        phase=exp(1i*2*pi*f.*(E*x(31:37)*1e-3).');corrected=H.*phase;
        residual=complex(zeros(nf,16));coeff=zeros(size(highBasis,2),16);value=0;
        for channel=1:16
            target=(corrected(:,channel)-low(:,channel)).*sqrtW(:,channel)/normalizer;
            targetReal=[real(target);imag(target)];factor=qrHigh{channel};
            projection=factor.Q.'*targetReal;ep=targetReal-factor.Q*projection;
            residual(:,channel)=ep(1:nf)+1i*ep(nf+1:end);value=value+sum(ep.^2);
            if nargout>=3,coeff(:,channel)=(factor.R\(projection*normalizer))./factor.scale.';end
        end
        g=zeros(37,1);
        if wantGradient
            weightedResidual=residual.*sqrtW/normalizer;
            g(1:30)=-2*real(reshape(jacobian,nf*16,30)'*weightedResidual(:));
            delayDerivative=(1i*2*pi*f*1e-3).*corrected;
            g(31:37)=E.'*(2*real(sum(conj(weightedResidual).*delayDerivative,1))).';
        end
        if nargout>=3
            high=highBasis*coeff;detail=meta;detail.H=(low+high)./phase;
            detail.H_low_undelayed=low;detail.H_high_undelayed=high;
            detail.high_coefficients=coeff;
            detail.corrected_training_residual=residual*normalizer./sqrtW;
        end
    end

    function [low,meta,jacobian]=lowModel(x,wantGradient)
        C=zeros(4,3);B=zeros(3,4);fn=zeros(3,1);ze=zeros(3,1);
        for mode=1:3
            offset=(mode-1)*9;C(anchors(mode),mode)=1;
            C(freeRows{mode},mode)=x(offset+(1:3));
            B(mode,:)=gains(mode)*x(offset+(4:7)).';
            fn(mode)=exp(x(offset+8));ze(mode)=exp(x(offset+9));
        end
        omega=2*pi*fn;K=diag(omega.^2);diagonalD=2*ze.*omega;
        rho=x(28:30);R=[1 rho(1) rho(2);rho(1) 1 rho(3);rho(2) rho(3) 1];
        rootD=sqrt(diagonalD);D=(rootD*rootD.').*R;
        low=complex(zeros(nf,16));jacobian=[];
        if wantGradient
            jacobian=complex(zeros(nf,16,30));dK=zeros(3,3,9);dD=zeros(3,3,9);
            for mode=1:3
                e=zeros(3);e(mode,mode)=1;
                dK(:,:,2*mode-1)=2*K(mode,mode)*e;
                dD(:,:,2*mode-1)=.5*(e*D+D*e);dD(:,:,2*mode)=dD(:,:,2*mode-1);
            end
            pairs=[1 2;1 3;2 3];
            for pair=1:3
                ii=pairs(pair,1);jj=pairs(pair,2);dD(ii,jj,6+pair)=rootD(ii)*rootD(jj);
                dD(jj,ii,6+pair)=dD(ii,jj,6+pair);
            end
            dynamicIndices=[8 9 17 18 26 27 28 29 30];
        end
        for frequency=1:nf
            s=1i*2*pi*f(frequency);inverseZ=(s^2*eye(3)+s*D+K)\eye(3);
            F=inverseZ*K;response=C*F*B;low(frequency,:)=response(:).';
            if wantGradient
                FB=F*B;CF=C*F;
                for mode=1:3
                    offset=(mode-1)*9;
                    for row=1:3
                        derivative=zeros(4);derivative(freeRows{mode}(row),:)=FB(mode,:);
                        jacobian(frequency,:,offset+row)=derivative(:).';
                    end
                    for input=1:4
                        derivative=zeros(4);derivative(:,input)=CF(:,mode)*gains(mode);
                        jacobian(frequency,:,offset+3+input)=derivative(:).';
                    end
                end
                for dynamic=1:9
                    dF=inverseZ*(dK(:,:,dynamic)-(dK(:,:,dynamic)+s*dD(:,:,dynamic))*F);
                    derivative=C*dF*B;jacobian(frequency,:,dynamicIndices(dynamic))=derivative(:).';
                end
            end
        end
        meta=struct('C_low',C,'B_low',B,'K',K,'D',D,'R',R,'rho',rho, ...
            'frequency_hz',fn,'damping',ze, ...
            'frequency_parameter_note','sqrt(diag(K))/(2*pi); damping entries are diagonal coordinates, not coupled modal damping ratios');
    end

    function [c,ceq,gc,gceq]=dampingConstraint(x)
        r=x(28:30);detR=1+2*prod(r)-sum(r.^2);
        c=[r.^2-(1-positiveMargin);positiveMargin-detR];ceq=[];
        if nargout>=3
            gc=zeros(numel(x),4);gc(28:30,1:3)=diag(2*r);
            gc(28:30,4)=2*[r(1)-r(2)*r(3);r(2)-r(1)*r(3);r(3)-r(1)*r(2)];gceq=[];
        end
    end

    function fit=assemble(x,flag,info)
        [value,g,fit]=objective(x);[uy,sy]=svd(fit.C_low);[uu,su]=svd(fit.B_low.');
        ny=uy(:,end);nu=uu(:,end);
        if ny.'*nominal.ny<0,ny=-ny;end
        if nu.'*nominal.nu<0,nu=-nu;end
        fit.Ny=ny;fit.Nu=nu;fit.coordinates_ms=x(31:37);
        fit.delay_seconds=reshape(E*x(31:37),4,4)*1e-3;
        fit.exitflag=flag;fit.optimizer_output=info;fit.x=x;fit.train_objective=value;
        fit.kkt_firstorderopt=info.firstorderopt;fit.train_weighted_error_pct=100*sqrt(value);fit.grad_inf=norm(g,inf);
        fit.spatial_singular_y=diag(sy);fit.spatial_singular_u=diag(su);
        fit.spatial_rank_valid=rank(fit.C_low)==3 && rank(fit.B_low)==3;
        cy=svd(fit.C_low./vecnorm(fit.C_low,2,1));cu=svd(fit.B_low./vecnorm(fit.B_low,2,2));
        fit.normalized_spatial_sigma3_over1=[cy(3)/cy(1),cu(3)/cu(1)];
        fit.correlation_eigenvalues=eig(fit.R);fit.damping_eigenvalues=eig(fit.D);
        fit.support_poles=eig([zeros(3),eye(3);-fit.K,-fit.D]);fit.stable=all(real(fit.support_poles)<0);
        positivePoles=fit.support_poles(imag(fit.support_poles)>1e-8);
        [~,order]=sort(abs(positivePoles));positivePoles=positivePoles(order);
        fit.coupled_modal_frequency_hz=abs(positivePoles)/(2*pi);
        fit.coupled_modal_damping=-real(positivePoles)./abs(positivePoles);
        fit.damping_constraint_values=dampingConstraint(x);
        fit.bound_violation=max([0;lb-x;x-ub]);fit.linear_violation=max([0;linearA*x-linearB]);
        fit.constraint_violation=max([0;fit.damping_constraint_values;fit.bound_violation;fit.linear_violation]);
        fit.feasible=fit.constraint_violation<=1e-8 && min(fit.correlation_eigenvalues)>0;
        fit.model='Restricted coupled-support model with symmetric positive damping + fixed upper poles/free numerators/D';
        fit.low_null_check=[norm(ny.'*fit.C_low),norm(fit.B_low*nu)];
        fit.training_relative_error_pct=100*norm(fit.H-H,'fro')/norm(H,'fro');
        fit.validation_relative_error_pct=100*norm(fit.H-Hvalidation,'fro')/norm(Hvalidation,'fro');
    end

    function stop=progress(~,optimValues,state)
        stop=false;
        if strcmp(state,'iter') && mod(optimValues.iteration,25)==0
            fprintf('z=%g coupled %s iter %d weighted %.5f%% opt %.3g sec %.1f\n', ...
                source.position_mm,stageName,optimValues.iteration,100*sqrt(optimValues.fval),optimValues.firstorderopt,toc(stageTimer));
            logEntries(end+1,:)={stageName,optimValues.iteration,optimValues.fval,optimValues.firstorderopt,toc(stageTimer)};
            progressLog=logEntries;save(fullfile(out,['optimization_' tag '_log.mat']),'progressLog');
        end
    end
end

function tag=positionTag(z)
if z==0,tag='z0';elseif z==-.5,tag='zm0p5';elseif z==.5,tag='z0p5';elseif z==1,tag='z1';else,error('Unsupported position');end
end
