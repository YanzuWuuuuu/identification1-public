function calib=run_structured_calibration(positionIndex,lowHalfWindow)
% Exploratory hardware calibration: shared low-frequency spaces and delays.
% SRIVC supplies initialization. Final estimation is constrained output-error
% least squares with analytic envelope derivatives, not another SRIVC step.
% The three support modes use real rank-one numerators (proportional damping).
% Upper-band poles remain at their training SRIVC estimates; their free real
% numerators and D are re-estimated as linear nuisance parameters each step.
if nargin<1,positionIndex=3;end
if nargin<2,lowHalfWindow=45;end
base=fileparts(mfilename('fullpath'));out=fullfile(base,'structured_calibration_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'position_spring_identification_20261016','position_campaign.mat'),'campaign');
source=S.campaign.results{positionIndex};assert(source.position_mm~=-1,'Exclude the anomalous record.');
P=load(fullfile(base,'static_decoupling_20261016','decoupled_campaign.mat'),'pack');
tag=positionTag(source.position_mm);ft=source.selected.fit;
if lowHalfWindow~=45
    assert(positionIndex==3 && ismember(lowHalfWindow,[15 30]),'Only prespecified z0 low-window sensitivity is supported.');
    L=load(fullfile(base,'spring_identification_20261016','spring_frf_data.mat'),'data');
    index=find([L.data.low.half_window]==lowHalfWindow,1);assert(~isempty(index));
    lowMask=source.f<25;lowF=L.data.low_f>=3;
    assert(max(abs(source.f(lowMask)-L.data.low_f(lowF)))<1e-6,'Low frequency grids differ.');
    source.trainingFRF(lowMask,:)=L.data.low(index).train(lowF,:);
    source.validationFRF(lowMask,:)=L.data.low(index).validation(lowF,:);
    tag=[tag '_lowhalf' num2str(lowHalfWindow)];
end
f=source.f(:);H=source.trainingFRF;W=source.weights;
assert(all(isfinite(H(:))) && all(W(:)>0));
nf=numel(f);nm=numel(ft.reference_hz);E=source.delay_design;
assert(isequal(size(E),[16 7]) && nm==15 && ft.direct);
normalizer=norm(H.*sqrt(W),'fro');sqrtW=sqrt(W);
q=1i*f./ft.reference_hz(:).';A=1+q.*ft.a(:,1).'+q.^2.*ft.a(:,2).';
highBasis=complex(zeros(nf,2*(nm-3)+1));
for m=4:nm,highBasis(:,2*(m-4)+1)=1./A(:,m);highBasis(:,2*(m-4)+2)=q(:,m)./A(:,m);end
highBasis(:,end)=1;
qrHigh=cell(16,1);highConditions=zeros(16,1);
for c=1:16
    bw=highBasis.*sqrtW(:,c);br=[real(bw);imag(bw)];scale=max(vecnorm(br),realmin);
    [Q,R]=qr(br./scale,0);qrHigh{c}=struct('Q',Q,'R',R,'scale',scale);highConditions(c)=cond(R);
end

x0=zeros(34,1);lb=-20*ones(34,1);ub=20*ones(34,1);anchors=zeros(3,1);freeRows=cell(3,1);gains=zeros(3,1);
Yspace=[];Uspace=[];modeRank=zeros(3,1);
for m=1:3
    C0=reshape(ft.coefficients(2*m-1,:),4,4);C1=reshape(ft.coefficients(2*m,:),4,4);
    [U,sv,V]=svd(C0);cm=U(:,1);[~,anchor]=max(abs(cm));
    if cm(anchor)<0,cm=-cm;V(:,1)=-V(:,1);end
    anchors(m)=anchor;freeRows{m}=setdiff(1:4,anchor);gains(m)=sv(1,1)*cm(anchor);
    cm=cm/cm(anchor);ind=(m-1)*9;
    x0(ind+(1:3))=cm(freeRows{m});x0(ind+(4:7))=V(:,1);
    x0(ind+8)=log(ft.frequency_hz(m));x0(ind+9)=log(ft.damping(m));
    lb(ind+8)=log(4);ub(ind+8)=log(24.5);lb(ind+9)=log(.00003);ub(ind+9)=log(.95);
    poles=roots([ft.a(m,2),ft.a(m,1),1]);pole=poles(imag(poles)>0);
    residue=2*pi*ft.reference_hz(m)*(C0+pole*C1)/(ft.a(m,1)+2*ft.a(m,2)*pole);
    Yspace=[Yspace,real(residue),imag(residue)]; %#ok<AGROW>
    Uspace=[Uspace,real(residue).',imag(residue).']; %#ok<AGROW>
    ss=svd(residue);modeRank(m)=ss(2)/ss(1);
end
x0(28:34)=source.coordinates_ms(:);lb(28:34)=[0;-ones(6,1)];ub(28:34)=ones(7,1);
[Ybasis,Sy]=svd(Yspace);[Ubasis,Su]=svd(Uspace);sy=diag(Sy);su=diag(Su);
nominal=struct('ny',P.pack.ny,'nu',P.pack.nu, ...
    'Ty',P.pack.views(2).Ty(1:3,:),'Tu',P.pack.views(2).Tu(:,1:3));
calib=struct('position_mm',source.position_mm,'tag',tag,'f',f,'H_training',H, ...
    'H_validation',source.validationFRF,'nominal',nominal,'initial_free',source, ...
    'source_geometry',P.pack.geometry_source,'source_data',source.source_file, ...
    'method','SRIVC initialization + structured support-mode variable-projection output-error + seven delays', ...
    'model_scope','Three real rank-one support modes; fixed high poles and free high numerators/D; no channel filters', ...
    'weight_contract','Frozen original training weights, inverse magnitude including original low-band coefficient 0.03', ...
    'validation_scope','Same-record second-half LRM, total closed-loop commands; not unbiased plant certification', ...
    'precheck',struct('singular_values_y',sy,'singular_values_u',su, ...
    'sigma4_over3_y',sy(4)/sy(3),'sigma4_over3_u',su(4)/su(3), ...
    'sigma3_over1_y',sy(3)/sy(1),'sigma3_over1_u',su(3)/su(1), ...
    'pole_sigma2_over1',modeRank,'unconstrained_ny',Ybasis(:,end),'unconstrained_nu',Ubasis(:,end)), ...
    'low_lrm_half_window',lowHalfWindow,'anchors',anchors,'amplitude_scales',gains,'high_basis',highBasis, ...
    'high_scaled_QR_condition',highConditions,'output_dir',out);

% Verify the profiled analytic derivative before optimizing real data.
[v0,g0]=objective(x0);gfd=zeros(size(x0));
for j=1:numel(x0)
    h=2e-6*max(1,abs(x0(j)));xp=x0;xm=x0;xp(j)=xp(j)+h;xm(j)=xm(j)-h;
    gfd(j)=(objective(xp)-objective(xm))/(2*h);
end
calib.gradient_check=norm(gfd-g0)/max(norm(gfd),1e-12);
assert(calib.gradient_check<2e-5,'Analytic profile gradient failed: %.3g',calib.gradient_check);
fprintf('z=%g precheck sy4/3 %.4f su4/3 %.4f, gradient error %.3g, initial cost %.5g\n', ...
    source.position_mm,sy(4)/sy(3),su(4)/su(3),calib.gradient_check,v0);

logEntries={};stageName='fixed';stageTimer=tic;
options=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true, ...
    'Display','none','MaxIterations',350,'MaxFunctionEvaluations',1800, ...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-9,'ConstraintTolerance',1e-9,'OutputFcn',@progress);
fixedStart=x0;
[xf,~,flag,info]=fmincon(@fixedObjective,x0(1:27),[],[],[],[],lb(1:27),ub(1:27),[],options);
xFixed=[xf;x0(28:34)];calib.fits.fixed=assemble(xFixed,flag,info);
save(fullfile(out,['calibration_' tag '_progress.mat']),'calib','logEntries','-v7.3');
fprintf('z=%g fixed-delay: weighted %.5f%% exit %d opt %.3g elapsed %.1fs\n', ...
    source.position_mm,calib.fits.fixed.train_weighted_error_pct,flag,info.firstorderopt,toc(stageTimer));

% Two prespecified delay starts. Selection only by training cost among
% feasible, converged candidates; if neither converges the result is flagged.
linearA=[zeros(16,27),E;zeros(16,27),-E];linearB=[ones(16,1);zeros(16,1)];
starts=[xFixed,xFixed];starts(28:34,2)=[.3;zeros(6,1)];
trials=cell(1,2);
for trial=1:2
    stageName=sprintf('joint_%d',trial);stageTimer=tic;
    [xx,~,flag,info]=fmincon(@objective,starts(:,trial),linearA,linearB,[],[],lb,ub,[],options);
    trials{trial}=assemble(xx,flag,info);
    fprintf('z=%g %s: weighted %.5f%% exit %d opt %.3g elapsed %.1fs\n', ...
        source.position_mm,stageName,trials{trial}.train_weighted_error_pct,flag,info.firstorderopt,toc(stageTimer));
    calib.trials=trials;save(fullfile(out,['calibration_' tag '_progress.mat']),'calib','logEntries','-v7.3');
end
valid=find(cellfun(@(a)a.exitflag>0 && a.feasible && a.spatial_rank_valid && a.kkt_firstorderopt<=1e-6,trials));
if isempty(valid),valid=1:2;calib.selection_status='No fully accepted candidate; best training candidate retained as diagnostic';
else,calib.selection_status='Minimum training cost among feasible positive-exit spatial-rank-valid candidates with KKT <= 1e-6';end
costs=cellfun(@(a)a.train_objective,trials(valid));[~,ii]=min(costs);calib.selected_trial=valid(ii);
calib.fits.joint=trials{valid(ii)};calib.trials=trials;calib.optimization_log=logEntries;
calib.start_agreement_deg=[acosd(min(1,abs(trials{1}.Ny.'*trials{2}.Ny))), ...
    acosd(min(1,abs(trials{1}.Nu.'*trials{2}.Nu)))];
calib.bounds=struct('lower',lb,'upper',ub,'linearA',linearA,'linearB',linearB);
calib.saved_at=char(datetime('now'));
save(fullfile(out,['calibration_' tag '.mat']),'calib','-v7.3');
fprintf('SAVED %s, selected start %d, %s\n',tag,calib.selected_trial,calib.selection_status);

    function [value,g]=fixedObjective(xx)
        [value,gg]=objective([xx;fixedStart(28:34)]);g=gg(1:27);
    end

    function [value,g,detail]=objective(x)
        [low,C,B,fn,ze,basis,df,dz]=lowModel(x);
        phase=exp(1i*2*pi*f.*(E*x(28:34)*1e-3).');corrected=H.*phase;
        value=0;g=zeros(34,1);coeff=zeros(size(highBasis,2),16);
        residual=complex(zeros(nf,16));
        for c=1:16
            target=(corrected(:,c)-low(:,c)).*sqrtW(:,c)/normalizer;
            targetReal=[real(target);imag(target)];qrC=qrHigh{c};
            projection=qrC.Q.'*targetReal;ep=targetReal-qrC.Q*projection;
            er=ep(1:nf)+1i*ep(nf+1:end);value=value+sum(abs(er).^2);
            dweight=conj(er).*sqrtW(:,c)/normalizer;
            y=mod(c-1,4)+1;u=floor((c-1)/4)+1;
            for m=1:3
                ind=(m-1)*9;row=find(freeRows{m}==y,1);
                if ~isempty(row),g(ind+row)=g(ind+row)-2*real(sum(dweight.*basis(:,m)))*B(m,u);end
                g(ind+3+u)=g(ind+3+u)-2*real(sum(dweight.*basis(:,m)))*C(y,m)*gains(m);
                product=C(y,m)*B(m,u);
                g(ind+8)=g(ind+8)-2*real(sum(dweight.*df(:,m)))*product;
                g(ind+9)=g(ind+9)-2*real(sum(dweight.*dz(:,m)))*product;
            end
            dg=2*real(sum(dweight.*(1i*2*pi*f*1e-3).*corrected(:,c)));
            g(28:34)=g(28:34)+E(c,:).'*dg;
            if nargout>=3
                coeff(:,c)=(qrC.R\(projection*normalizer))./qrC.scale.';
                residual(:,c)=er*normalizer./sqrtW(:,c);
            end
        end
        if nargout>=3
            high=highBasis*coeff;prediction=(low+high)./phase;
            detail=struct('H',prediction,'H_low_undelayed',low,'H_high_undelayed',high, ...
                'C_low',C,'B_low',B,'frequency_hz',fn,'damping',ze,'high_coefficients',coeff, ...
                'corrected_training_residual',residual);
        end
    end

    function [low,C,B,fn,ze,bas,df,dz]=lowModel(x)
        C=zeros(4,3);B=zeros(3,4);fn=zeros(3,1);ze=fn;
        for m=1:3
            ind=(m-1)*9;C(anchors(m),m)=1;C(freeRows{m},m)=x(ind+(1:3));
            B(m,:)=gains(m)*x(ind+(4:7)).';fn(m)=exp(x(ind+8));ze(m)=exp(x(ind+9));
        end
        qq=1i*f./fn.';aa=1+2*qq.*ze.'+qq.^2;bas=1./aa;
        df=(2*qq.*ze.'+2*qq.^2)./aa.^2;dz=-2*qq.*ze.'./aa.^2;
        low=complex(zeros(nf,16));
        for m=1:3,co=C(:,m)*B(m,:);low=low+bas(:,m)*co(:).';end
    end

    function fit=assemble(x,flag,info)
        [value,g,fit]=objective(x);[Uy,SyLow]=svd(fit.C_low);[Uu,SuLow]=svd(fit.B_low.');
        ny=Uy(:,end);nu=Uu(:,end);
        if ny.'*nominal.ny<0,ny=-ny;end
        if nu.'*nominal.nu<0,nu=-nu;end
        fit.Ny=ny;fit.Nu=nu;fit.coordinates_ms=x(28:34);fit.delay_seconds=reshape(E*x(28:34),4,4)*1e-3;
        fit.exitflag=flag;fit.optimizer_output=info;fit.x=x;fit.train_objective=value;
        fit.kkt_firstorderopt=info.firstorderopt;
        fit.train_weighted_error_pct=100*sqrt(value);fit.grad_inf=norm(g,inf);
        fit.spatial_singular_y=diag(SyLow);fit.spatial_singular_u=diag(SuLow);
        fit.spatial_rank_valid=rank(fit.C_low)==3 && rank(fit.B_low)==3;
        cy=svd(fit.C_low./vecnorm(fit.C_low,2,1));cu=svd(fit.B_low./vecnorm(fit.B_low,2,2));
        fit.normalized_spatial_sigma3_over1=[cy(3)/cy(1),cu(3)/cu(1)];
        fit.feasible=all(E*x(28:34)>=-1e-8 & E*x(28:34)<=1+1e-8);
        fit.model='Three proportional-support modes + free upper-mode nuisance numerators and D';
        fit.low_null_check=[norm(ny.'*fit.C_low),norm(fit.B_low*nu)];
        fit.training_relative_error_pct=100*norm(fit.H-H,'fro')/norm(H,'fro');
    end

    function stop=progress(~,optimValues,state)
        stop=false;
        if strcmp(state,'iter') && mod(optimValues.iteration,25)==0
            msg=sprintf('z=%g %s iter %d weighted %.5f%% opt %.3g sec %.1f',source.position_mm, ...
                stageName,optimValues.iteration,100*sqrt(optimValues.fval),optimValues.firstorderopt,toc(stageTimer));
            fprintf('%s\n',msg);logEntries(end+1,:)={stageName,optimValues.iteration,optimValues.fval,optimValues.firstorderopt,toc(stageTimer)};
            progressLog=logEntries;save(fullfile(out,['optimization_' tag '_log.mat']),'progressLog');
        end
    end
end

function tag=positionTag(z)
if z==0,tag='z0';elseif z==-.5,tag='zm0p5';elseif z==.5,tag='z0p5';elseif z==1,tag='z1';else,error('Unsupported position');end
end
