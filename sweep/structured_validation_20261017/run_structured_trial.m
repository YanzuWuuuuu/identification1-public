function result=run_structured_trial(tag,strict)
% Fixed-SRIVC-pole structural refinement; no raw data or baseline is changed.
% Low block: common real 3D spaces, sum Re(residues)=0. Necessary relaxation,
% not a claim of a complete three-DOF second-order mechanical realization.
if nargin<1,tag='z0';end
if nargin<2,strict=false;end
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
source=load(fullfile(root,'sweep','srivc_retest_20261016',['model_' tag '.mat']));
src=source.result;fit=src.selected.fit;f=src.f(:);nf=numel(f);nm=numel(fit.reference_hz);
assert(nm==15);H=src.trainingFRF;W=src.weights;sw=sqrt(W);
normalizer=norm(H.*sw,'fro');s=2i*pi*f;
lambda=zeros(nm,1);residues=complex(zeros(4,4,nm));
for m=1:nm
    a=fit.a(m,:);q=(-a(1)+sqrt(complex(a(1)^2-4*a(2))))/(2*a(2));
    omega=2*pi*fit.reference_hz(m);lambda(m)=omega*q;
    residues(:,:,m)=reshape(omega*(fit.coefficients(2*m-1,:)+q*fit.coefficients(2*m,:))/(a(1)+2*a(2)*q),4,4);
end
kr=1./(s-lambda.')+1./(s-conj(lambda).');
ki=1i*(1./(s-lambda.')-1./(s-conj(lambda).'));
reconstructed=kr*reshape(real(residues),16,nm).'+ki*reshape(imag(residues),16,nm).'+fit.coefficients(end,:);
conversionError=norm(reconstructed-fit.prediction,'fro')/norm(fit.prediction,'fro');
assert(conversionError<1e-10,'Complex residue conversion failed.');
left=[];right=[];
for m=1:3
    left=[left,real(residues(:,:,m)),imag(residues(:,:,m))]; %#ok<AGROW>
    right=[right,real(residues(:,:,m)).',imag(residues(:,:,m)).']; %#ok<AGROW>
end
[uy,sy]=svd(left,'econ');[uu,su]=svd(right,'econ');
Ly0=uy(:,1:3);Lu0=uu(:,1:3);ny0=uy(:,4);nu0=uu(:,4);
lowPart=residues(:,:,1:3);lowGain=norm(lowPart(:));
klow=[kr(:,1)-kr(:,3),kr(:,2)-kr(:,3),ki(:,1:3)]*lowGain;
highGain=zeros(nm,1);anchors=zeros(nm,1);freeRows=cell(nm,1);
% 6 graph coordinates; 45 low coefficients; 12*14 high factors; 16 D; 7 delays.
nBase=6+45+12*14+16;nPar=nBase+7;x0=zeros(nPar,1);
for m=4:nm
    [us,ss,vs]=svd(residues(:,:,m));[~,anchor]=max(abs(us(:,1)));
    anchors(m)=anchor;freeRows{m}=setdiff(1:4,anchor);
    highGain(m)=ss(1,1);c=us(:,1)/us(anchor,1);b=us(anchor,1)*conj(vs(:,1));
    off=51+(m-4)*14;x0(off+(1:14))=[real(c(freeRows{m}));imag(c(freeRows{m}));real(b);imag(b)];
end
dScale=max(abs(fit.coefficients(end,:)).',median(abs(fit.coefficients(end,:)))*.1);
dScale=max(dScale,1e-15);dIndices=51+12*14+(1:16);
x0(dIndices)=fit.coefficients(end,:).'./dScale;
E=zeros(16,7);
for input=1:4
    for output=1:4
        c=output+4*(input-1);E(c,input)=1;
        if output>1,E(c,4+output-1)=1;end
    end
end
tau=src.delay_matrix_seconds*1e3;
x0(nBase+(1:7))=[tau(1,:).';tau(2:4,1)-tau(1,1)];
assert(norm(reshape(E*x0(nBase+(1:7)),4,4)-tau,'fro')<1e-9);
% Initialize the 45 low coefficients by one weighted linear least-squares solve.
[r0,j0]=residual(x0);idx=7:51;x0(idx)=-j0(:,idx)\r0;
lb=-100*ones(nPar,1);ub=100*ones(nPar,1);lb(1:6)=-2;ub(1:6)=2;
lb(nBase+(1:7))=[zeros(4,1);-.6*ones(3,1)];
ub(nBase+(1:7))=[1.2*ones(4,1);.6*ones(3,1)];
assert(all(x0>=lb & x0<=ub));
% Directional derivative checks span both structure and delay coordinates.
[r0,j0]=residual(x0);rng(20261017);checks=zeros(3,1);
for k=1:3
    v=randn(nPar,1);v=v/norm(v);h=2e-6;
    fd=(residual(x0+h*v)-residual(x0-h*v))/(2*h);
    checks(k)=norm(fd-j0*v)/max(norm(fd),1e-12);
end
assert(max(checks)<2e-5,'Analytic Jacobian check failed.');
result=struct('tag',tag,'position_mm',src.position_mm,'f',f,'H_training',H, ...
    'H_validation',src.validationFRF,'weights',W,'baseline_H',src.selected.H, ...
    'baseline_delay_seconds',src.delay_matrix_seconds,'lambda',lambda, ...
    'source_file',src.source_file,'model_scope','Common support-space necessary relaxation; high complex rank-one residues; SRIVC poles fixed', ...
    'data_scope','Previously inspected within-record held-out development check; closed-loop U/Y FRF, not unbiased IV', ...
    'residue_conversion_error',conversionError,'directional_Jacobian_errors',checks, ...
    'initial_support_singular_values_y',diag(sy),'initial_support_singular_values_u',diag(su), ...
    'x_initial',x0,'lower_bounds',lb,'upper_bounds',ub,'delay_design',E);
result.baseline_weighted_error_pct=100*norm((src.selected.H-H).*sw,'fro')/normalizer;
log=[];stage='fixed';timer=tic;
ftol=1e-10;stol=1e-9;if strict,ftol=1e-14;stol=1e-12;end
options=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','none', ...
    'MaxIterations',120,'MaxFunctionEvaluations',500,'FunctionTolerance',ftol, ...
    'StepTolerance',stol,'OptimalityTolerance',1e-7,'OutputFcn',@progress);
result.solver_configuration=struct('FunctionTolerance',ftol,'StepTolerance',stol,'OptimalityTolerance',1e-7,'MaxIterations',120,'strict_recheck',strict);
result.code_snapshot=fileread([mfilename('fullpath') '.m']);
fprintf('%s conversion %.3g Jacobian %.3g initial %.5g%% baseline %.5g%%\n',tag,conversionError,max(checks),100*norm(r0),result.baseline_weighted_error_pct);
[xx,~,~,flag,info]=lsqnonlin(@fixedResidual,x0(1:nBase),lb(1:nBase),ub(1:nBase),options);
xf=[xx;x0(nBase+(1:7))];result.fixed=assemble(xf,flag,info);
save(fullfile(out,['trial_' tag '_progress.mat']),'result','-v7.3');
stage='delay';timer=tic;
[xx,~,~,flag,info]=lsqnonlin(@residual,xf,lb,ub,options);
result.delay=assemble(xx,flag,info);result.log=log;
result.saved_at=char(datetime('now'));
save(fullfile(out,['trial_' tag '.mat']),'result','-v7.3');
fprintf('SAVED %s fixed %.5g%% (exit %d), delay %.5g%% (exit %d)\n',tag, ...
    result.fixed.train_weighted_error_pct,result.fixed.exitflag,result.delay.train_weighted_error_pct,result.delay.exitflag);

    function [r,j]=fixedResidual(x)
        [r,jj]=residual([x;x0(nBase+(1:7))]);j=jj(:,1:nBase);
    end
    function [r,j,meta]=residual(x)
        wantJ=nargout>=2;L=Ly0+ny0*x(1:3).';R=Lu0+nu0*x(4:6).';
        map=kron(R,L);A=reshape(x(7:51),9,5);core=klow*A.';
        low=core*map.';pred=low;res=complex(zeros(4,4,nm));
        res(:,:,1)=lowGain*L*(reshape(A(:,1),3,3)+1i*reshape(A(:,3),3,3))*R.';
        res(:,:,2)=lowGain*L*(reshape(A(:,2),3,3)+1i*reshape(A(:,4),3,3))*R.';
        res(:,:,3)=lowGain*L*(-reshape(A(:,1)+A(:,2),3,3)+1i*reshape(A(:,5),3,3))*R.';
        if wantJ
            jc=complex(zeros(nf*16,nPar));
            for k=1:3
                el=zeros(4,3);el(:,k)=ny0;er=zeros(4,3);er(:,k)=nu0;
                tmp=core*kron(R,el).';jc(:,k)=tmp(:);
                tmp=core*kron(er,L).';jc(:,3+k)=tmp(:);
            end
            for k=1:5
                for aidx=1:9
                    tmp=klow(:,k)*map(:,aidx).';jc(:,6+(k-1)*9+aidx)=tmp(:);
                end
            end
        end
        for m=4:nm
            off=51+(m-4)*14;c=complex(zeros(4,1));c(anchors(m))=1;
            c(freeRows{m})=x(off+(1:3))+1i*x(off+(4:6));
            b=x(off+(7:10))+1i*x(off+(11:14));
            rr=highGain(m)*c*b.';res(:,:,m)=rr;
            pred=pred+kr(:,m)*real(rr(:)).'+ki(:,m)*imag(rr(:)).';
            if wantJ
                for p=1:14
                    dc=complex(zeros(4,1));db=complex(zeros(4,1));
                    if p<=3,dc(freeRows{m}(p))=1;
                    elseif p<=6,dc(freeRows{m}(p-3))=1i;
                    elseif p<=10,db(p-6)=1;
                    else,db(p-10)=1i;end
                    dr=highGain(m)*(dc*b.'+c*db.');
                    tmp=kr(:,m)*real(dr(:)).'+ki(:,m)*imag(dr(:)).';jc(:,off+p)=tmp(:);
                end
            end
        end
        direct=(x(dIndices).*dScale).';pred=pred+direct;
        delayms=E*x(nBase+(1:7));phase=exp(-s*(delayms.'*1e-3));
        prediction=pred.*phase;err=(prediction-H).*sw/normalizer;r=[real(err(:));imag(err(:))];
        if wantJ
            for p=1:16
                tmp=zeros(nf,16);tmp(:,p)=dScale(p);jc(:,dIndices(p))=tmp(:);
            end
            jc=jc.*phase(:);
            for p=1:7
                tmp=-s.*prediction.*(E(:,p).'*1e-3);jc(:,nBase+p)=tmp(:);
            end
            jc=jc.*(sw(:)/normalizer);j=[real(jc);imag(jc)];
        end
        if nargout>=3
            meta=struct('H',prediction,'H_low',low.*phase,'residues',res,'direct',direct, ...
                'ny',(ny0-Ly0*x(1:3))/sqrt(1+norm(x(1:3))^2), ...
                'nu',(nu0-Lu0*x(4:6))/sqrt(1+norm(x(4:6))^2), ...
                'delay_matrix_seconds',reshape(delayms*1e-3,4,4),'support_L',L,'support_R',R);
        end
    end
    function model=assemble(x,flag,info)
        [r,j,model]=residual(x);model.x=x;model.exitflag=flag;model.optimization=info;
        model.train_weighted_error_pct=100*norm(r);model.gradient_inf=norm(j.'*r,inf);
        model.iterations=info.iterations;model.parameter_count=numel(x);
        model.active_parameter_count=nBase+7*double(strcmp(stage,'delay'));
        model.support_sum_real_residue=norm(sum(real(model.residues(:,:,1:3)),3),'fro');
        model.max_input_output_compensation_gain=1;
        model.nonnegative_total_delay=all(model.delay_matrix_seconds(:)>=-1e-12);
        fprintf('%s %s train %.6g%% iter %d exit %d optimality %.3g time %.1fs\n',tag,stage, ...
            model.train_weighted_error_pct,info.iterations,flag,info.firstorderopt,toc(timer));
    end
    function stop=progress(~,optim,state)
        stop=false;
        if strcmp(state,'iter') && mod(optim.iteration,10)==0
            fprintf('%s %s iter %d weighted %.6g%% optimality %.3g time %.1fs\n', ...
                tag,stage,optim.iteration,100*sqrt(optim.resnorm),optim.firstorderopt,toc(timer));
            log=[log;double(strcmp(stage,'delay')),optim.iteration,optim.resnorm,optim.firstorderopt,toc(timer)]; %#ok<AGROW>
        end
    end
end
