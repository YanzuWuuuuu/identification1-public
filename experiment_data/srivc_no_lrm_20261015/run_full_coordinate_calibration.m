function result=run_full_coordinate_calibration
% Fixed z=0, first plate mode; same-record exploratory coordinate experiment.
% Baseline files are read only. No online controller or FRF is re-identified.
base=fileparts(mfilename('fullpath'));
out=fullfile(base,'full_coordinate_calibration_20261016');
if ~isfolder(out),mkdir(out);end
source=fullfile(base,'structured_calibration_20261016','calibration_z0.mat');
tmp=load(source,'calib');c=tmp.calib;j=c.fits.joint;ft=c.initial_free.selected.fit;
f=c.f(:);nf=numel(f);nm=numel(ft.reference_hz);
assert(c.position_mm==0 && nm==15 && all(diff(f)>0));
target=find(ft.frequency_hz>275 & ft.frequency_hz<290);
assert(isscalar(target),'Prespecified first plate mode is ambiguous.');
Pphys=c.nominal.Ty;Qphys=c.nominal.Tu;
scaleY=vecnorm(Pphys,2,2);scaleU=vecnorm(Qphys,2,1);
P=Pphys./scaleY;Q=Qphys./scaleU;
ny0=unit(c.nominal.ny);nu0=unit(c.nominal.nu);
ny=aligned(j.Ny,ny0);nu=aligned(j.Nu,nu0);
ell=j.delay_seconds;dy=ell(:,1)-ell(1,1);du=ell(1,:)-ell(1,1);
assert(norm(ell-(ell(1,1)+dy+du),'fro')<1e-12);
extraY=max([0;dy])-dy;extraU=max([0 du])-du;
assert(min([extraY(:);extraU(:)])>-1e-14);
phaseY=exp(-2i*pi*f*extraY.');phaseU=exp(-2i*pi*f*extraU);
commonDelay=ell(1,1)+max([0;dy])+max([0 du]);

% Audit all upper residues; only the prespecified target designs the maps.
residues=cell(nm,1);realModes=cell(nm,1);spatialRows={};
for m=4:nm
    rows=2*(m-4)+(1:2);
    C0=reshape(j.high_coefficients(rows(1),:),4,4);
    C1=reshape(j.high_coefficients(rows(2),:),4,4);
    R=modalResidue(C0,C1,ft.a(m,:),ft.reference_hz(m));
    shape=realRankOne(R);residues{m}=R;realModes{m}=shape;
    sv=svd(R);sy=svd([real(R),imag(R)]);su=svd([real(R).',imag(R).']);
    spatialRows(end+1,:)={m,ft.frequency_hz(m),sv(2)/sv(1),sy(2)/sy(1), ...
        su(2)/su(1),shape.error,abs(ny.'*shape.c),abs(shape.b.'*nu),m==target}; %#ok<AGROW>
end
spatialChecks=cell2table(spatialRows,'VariableNames',{'mode','frequency_hz', ...
    'complex_sigma2_over1','real_output_sigma2_over1','real_input_sigma2_over1', ...
    'best_common_real_rank1_relative_error','eta_output','eta_input','selected_target'});
shape=realModes{target};cf=shape.c;bf=shape.b;
alphaY=ny.'*cf;alphaU=bf.'*nu;
assert(abs(alphaY)>1e-8 && abs(alphaU)>1e-8,'Target projection is numerically zero.');
deltaP=-(P*cf)/alphaY*ny.';deltaQ=-nu*(bf.'*Q)/alphaU;
Pfull=P+deltaP;Qfull=Q+deltaQ;
gainCap=1.25;
lambdaY=maxFraction(P,deltaP,gainCap,'y',ny);
lambdaU=maxFraction(Q,deltaQ,gainCap,'u',nu);
Pbounded=P+lambdaY*deltaP;Qbounded=Q+lambdaU*deltaQ;
fprintf('Target %.6f Hz: real-rank1 error %.4f%%, projections y/u %.6f %.6f\n', ...
    ft.frequency_hz(target),100*shape.error,abs(alphaY),abs(alphaU));
fprintf('Gain cap %.2f gives fractions y/u %.8f %.8f\n',gainCap,lambdaY,lambdaU);

% Physical scaling is kept; only the fixed nominal scaling is used in metrics.
specs={ 'B0','Nominal static',P,Q,ny0,nu0,false; ...
        'B1','Estimated null spaces',P,Q,ny,nu,false; ...
        'B4','Estimated null spaces + delay alignment',P,Q,ny,nu,true; ...
        'F1','Full target-aware static coordinates',Pfull,Qfull,ny,nu,false; ...
        'F2','Full coordinates + delay alignment',Pfull,Qfull,ny,nu,true; ...
        'F3','Gain-limited coordinates + delay alignment',Pbounded,Qbounded,ny,nu,true};
rawTrain=reshape(c.H_training.',4,4,[]);rawVal=reshape(c.H_validation.',4,4,[]);
rawModel=reshape(j.H.',4,4,[]);
idx=2*(target-4)+(1:2);
targetFlat=c.high_basis(:,idx)*j.high_coefficients(idx,:);
targetDelayed=targetFlat.*exp(-2i*pi*f*ell(:).');
rawTarget=reshape(targetDelayed.',4,4,[]);
lowFlat=j.H_low_undelayed.*exp(-2i*pi*f*ell(:).');
rawLow=reshape(lowFlat.',4,4,[]);
methods=struct([]);gainRows={};checks=struct;
for m=1:size(specs,1)
    pm=specs{m,3};qm=specs{m,4};nym=specs{m,5};num=specs{m,6};
    Ty=[pm;nym.'];Tu=[qm,num];useDelay=specs{m,7};
    assert(all(isfinite(Ty(:))) && all(isfinite(Tu(:))) && rank(Ty)==4 && rank(Tu)==4);
    methods(m).code=specs{m,1};methods(m).name=specs{m,2}; %#ok<AGROW>
    methods(m).Ty_static=Ty;methods(m).Tu_static=Tu;
    methods(m).Ty_physical=[pm.*scaleY;nym.'];
    methods(m).Tu_physical=[qm.*scaleU,num];
    methods(m).use_delay=useDelay;
    if useDelay,py=phaseY;pu=phaseU;else,py=ones(nf,4);pu=ones(nf,4);end
    methods(m).H_training=transform(rawTrain,Ty,Tu,py,pu);
    methods(m).H_validation=transform(rawVal,Ty,Tu,py,pu);
    methods(m).H_model=transform(rawModel,Ty,Tu,py,pu);
    methods(m).H_target_diagnostic=transform(rawTarget,Ty,Tu,py,pu);
    methods(m).H_support_diagnostic=transform(rawLow,Ty,Tu,py,pu);
    rowG=vecnorm(pm,2,2)./vecnorm(P,2,2);
    colG=vecnorm(qm,2,1)./vecnorm(Q,2,1);
    gainRows(end+1,:)={specs{m,1},rowG(1),rowG(2),rowG(3),colG(1),colG(2),colG(3), ...
        norm(pm,2)/norm(P,2),norm(qm,2)/norm(Q,2),norm(Ty,2)/norm([P;ny.'],2), ...
        norm(Tu,2)/norm([Q,nu],2),cond(Ty),cond(Tu), ...
        norm(pm*j.C_low-P*j.C_low,'fro')/norm(P*j.C_low,'fro'), ...
        norm(j.B_low*qm-j.B_low*Q,'fro')/norm(j.B_low*Q,'fro')}; %#ok<AGROW>
    writematrix(methods(m).Ty_physical,fullfile(out,[specs{m,1} '_Ty_output.csv']));
    writematrix(methods(m).Tu_physical,fullfile(out,[specs{m,1} '_Tu_input.csv']));
end
gainMetrics=cell2table(gainRows,'VariableNames',{'method','output_Rx_gain','output_Ry_gain', ...
    'output_Z_gain','input_Rx_gain','input_Ry_gain','input_Z_gain','logical_output_operator_gain', ...
    'logical_input_operator_gain','full_output_operator_gain','full_input_operator_gain', ...
    'output_coordinate_condition','input_coordinate_condition','estimated_support_output_invariance', ...
    'estimated_support_input_invariance'});

% Algebra checks do not count as experimental decoupling evidence.
checks.target_direction_annihilation=[norm(Pfull*cf),norm(bf.'*Qfull)];
checks.support_null=[norm(ny.'*j.C_low)/norm(j.C_low,'fro'),norm(j.B_low*nu)/norm(j.B_low,'fro')];
checks.support_invariance=max(gainMetrics{4:6,{'estimated_support_output_invariance','estimated_support_input_invariance'}},[],'all');
checks.NN_same_static=relative(methods(4).H_validation(4,4,:),methods(2).H_validation(4,4,:));
checks.NN_same_delayed=max(relative(methods(5).H_validation(4,4,:),methods(3).H_validation(4,4,:)), ...
    relative(methods(6).H_validation(4,4,:),methods(3).H_validation(4,4,:)));
perfectTarget=cf*bf.';checks.ideal_rank1_cross= ...
    norm(Pfull*perfectTarget*nu)+norm(ny.'*perfectTarget*Qfull);
checks.aligned_target=relative(methods(3).H_target_diagnostic, ...
    transform(reshape(targetFlat.',4,4,[]),[P;ny.'],[Q,nu], ...
    repmat(exp(-2i*pi*f*commonDelay),1,4),ones(nf,4)));
assert(max([checks.target_direction_annihilation,checks.support_null,checks.support_invariance, ...
    checks.NN_same_static,checks.NN_same_delayed,checks.ideal_rank1_cross,checks.aligned_target])<1e-9);
assert(all(gainMetrics{6,2:11}<=gainCap+1e-10));

% Reproduce the old B0/B1/B4 values independently.
oldPath=fullfile(base,'structured_calibration_20261016','calibration_z0_evaluation.mat');
checks.baseline_saved_file=oldPath;
if isfile(oldPath)
    old=load(oldPath);vn=fieldnames(old);old=old.(vn{1});
    if isfield(old,'methods')
        checks.baseline_reproduction=zeros(3,1);
        for m=1:3,checks.baseline_reproduction(m)=relative(methods(m).H_validation,old.methods(m).H_validation);end
        assert(max(checks.baseline_reproduction)<1e-10,'Old baseline changed.');
    end
end

bands=[3 25;25 80;100 170;260 310;320 380;450 800;25 800];
labels={'3-25','25-80','100-170','260-310','320-380','450-800','25-800'};
blockNames={'logical_to_N','N_to_logical','N_to_N'};rows={};diagnosticRows={};
for b=1:size(bands,1)
    ix=f>=bands(b,1) & f<=bands(b,2)+1e-9;fb=f(ix);assert(numel(fb)>2);
    for block=1:3
        baseline=rmsBand(fb,blockData(methods(1).H_validation,block,ix));
        baselineTrain=rmsBand(fb,blockData(methods(1).H_training,block,ix));
        prior=rmsBand(fb,blockData(methods(3).H_validation,block,ix));
        for m=1:numel(methods)
            vv=blockData(methods(m).H_validation,block,ix);tt=blockData(methods(m).H_training,block,ix);
            hh=blockData(methods(m).H_model,block,ix);
            val=rmsBand(fb,vv);train=rmsBand(fb,tt);
            rows(end+1,:)={methods(m).code,labels{b},blockNames{block},val,train, ...
                val/baseline,train/baselineTrain,val/prior, ...
                rmsBand(fb,hh-vv)/max(val,realmin)}; %#ok<AGROW>
            targetValue=blockData(methods(m).H_target_diagnostic,block,ix);
            lowValue=blockData(methods(m).H_support_diagnostic,block,ix);
            diagnosticRows(end+1,:)={methods(m).code,labels{b},blockNames{block}, ...
                rmsBand(fb,targetValue),rmsBand(fb,lowValue),rmsBand(fb,hh-targetValue-lowValue), ...
                'Fitted component only; not measured modal ground truth'}; %#ok<AGROW>
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'method','band_hz','block','check_rms','train_rms', ...
    'check_over_B0','train_over_B0','check_over_B4','model_check_complex_error'});
componentMetrics=cell2table(diagnosticRows,'VariableNames',{'method','band_hz','block', ...
    'model_target_rms','model_support_rms','model_other_rms','scope'});

% Conditional sensitivity only: freeze training poles/delay/weights, refit
% free real numerators on each half; never use check directions for design.
refs=ft.reference_hz(:);aa=ft.a;refs(1:3)=j.frequency_hz;aa(1:3,:)=[2*j.damping,ones(3,1)];
q=1i*f./refs.';den=1+q.*aa(:,1).'+q.^2.*aa(:,2).';
Phi=complex(zeros(nf,2*nm+1));
for m=1:nm,Phi(:,2*m-1)=1./den(:,m);Phi(:,2*m)=q(:,m)./den(:,m);end
Phi(:,end)=1;unDelay=exp(2i*pi*f*ell(:).');
[coefT,condT]=refitNumerators(Phi,c.H_training.*unDelay,c.initial_free.weights);
[coefV,condV]=refitNumerators(Phi,c.H_validation.*unDelay,c.initial_free.weights);
sensT=spaceDiagnostic(coefT,aa,refs,target);sensV=spaceDiagnostic(coefV,aa,refs,target);
sensitivity=struct('training',sensT,'check',sensV,'design_condition_train',condT, ...
    'design_condition_check',condV,'free_output_null_angle_deg',angleBetween(sensT.ny,sensV.ny), ...
    'free_input_null_angle_deg',angleBetween(sensT.nu,sensV.nu), ...
    'target_output_angle_deg',angleBetween(sensT.target.c,sensV.target.c), ...
    'target_input_angle_deg',angleBetween(sensT.target.b,sensV.target.b), ...
    'scope','Conditional fixed-pole/delay sensitivity; not confidence intervals or independent blind test');

result=struct('source',source,'position_mm',0,'f',f,'target_index',target, ...
    'target_frequency_hz',ft.frequency_hz(target),'methods',methods,'metrics',metrics, ...
    'gain_metrics',gainMetrics,'spatial_checks',spatialChecks,'component_metrics',componentMetrics, ...
    'checks',checks,'sensitivity',sensitivity,'target_shape',shape,'residues',{residues}, ...
    'target_projections',[alphaY alphaU],'bounded_fractions',[lambdaY lambdaU], ...
    'gain_cap',gainCap,'extra_output_delay_seconds',extraY,'extra_input_delay_seconds',extraU, ...
    'aligned_common_delay_seconds',commonDelay,'fixed_output_scales',scaleY, ...
    'fixed_input_scales',scaleU,'matlab_version',version,'created',char(datetime('now')), ...
    'source_code',fileread([mfilename('fullpath') '.m']),'out_dir',out, ...
    'scope','Offline exploratory transforms of closed-loop-total-u FRFs. No actual controller and no new independent measurements.');
save(fullfile(out,'result.mat'),'result','-v7.3');
writetable(metrics,fullfile(out,'metrics.csv'));writetable(gainMetrics,fullfile(out,'gain_metrics.csv'));
writetable(spatialChecks,fullfile(out,'spatial_checks.csv'));
writetable(componentMetrics,fullfile(out,'component_metrics.csv'));
disp(gainMetrics);disp(metrics(ismember(metrics.band_hz,{'25-80','260-310'}) & ~strcmp(metrics.block,'N_to_N'),:));
disp(sensitivity);fprintf('Saved %s\n',fullfile(out,'result.mat'));
end

function shape=realRankOne(R)
rn=norm(R,'fro');A=R/rn;grid=linspace(0,pi,365);score=zeros(size(grid));
for i=1:numel(grid),score(i)=phaseScore(grid(i),A);end
inds=find(score(2:end-1)>=score(1:end-2) & score(2:end-1)>=score(3:end))+1;
candidate=[0 pi];values=[score(1) score(end)];
if score(1)>=score(2) && score(1)>=score(end-1)
    [ph,val]=fminbnd(@(p)-phaseScore(p,A),-grid(2),grid(2),optimset('TolX',1e-11));
    candidate(end+1)=mod(ph,pi);values(end+1)=-val;
end
for i=inds
    [ph,val]=fminbnd(@(p)-phaseScore(p,A),grid(i-1),grid(i+1),optimset('TolX',1e-11));
    candidate(end+1)=ph;values(end+1)=-val; %#ok<AGROW>
end
[~,best]=max(values);phase=candidate(best);
[U,~,V]=svd(real(exp(-1i*phase)*A));cv=U(:,1);bv=V(:,1);
[~,ii]=max(abs(cv));cv=cv*sign(cv(ii));[~,ii]=max(abs(bv));bv=bv*sign(bv(ii));
gamma=cv.'*R*bv;approx=gamma*cv*bv.';
shape=struct('c',cv,'b',bv,'gamma',gamma,'error',norm(R-approx,'fro')/rn,'phase',phase, ...
    'stationarity_y',norm((eye(4)-cv*cv.')*real((R*bv)*(R*bv)')*cv)/rn^2, ...
    'stationarity_u',norm((eye(4)-bv*bv.')*real((R.'*cv)*(R.'*cv)')*bv)/rn^2);
end
function value=phaseScore(p,A)
s=svd(real(exp(-1i*p)*A));value=s(1)^2;
end
function R=modalResidue(C0,C1,a,reference)
p=roots([a(2),a(1),1]);p=p(imag(p)>0);assert(isscalar(p));
R=2*pi*reference*(C0+p*C1)/(a(1)+2*a(2)*p);
end
function l=maxFraction(A,D,cap,side,anchor)
if gainValid(A+D,A,cap,side,anchor),l=1;return;end
lo=0;hi=1;
for k=1:70
    mid=(lo+hi)/2;if gainValid(A+mid*D,A,cap,side,anchor),lo=mid;else,hi=mid;end
end
l=lo;
end
function ok=gainValid(A,B,cap,side,anchor)
if side=='y'
    v=vecnorm(A,2,2)./vecnorm(B,2,2);fullA=[A;anchor.'];fullB=[B;anchor.'];
else
    v=vecnorm(A,2,1)./vecnorm(B,2,1);fullA=[A,anchor];fullB=[B,anchor];
end
ok=max(v)<=cap && norm(A,2)<=cap*norm(B,2) && norm(fullA,2)<=cap*norm(fullB,2);
end
function out=transform(H,Ty,Tu,py,pu)
out=complex(zeros(size(H)));
for k=1:size(H,3),out(:,:,k)=Ty*diag(py(k,:))*H(:,:,k)*diag(pu(k,:))*Tu;end
end
function v=unit(v)
v=v(:)/norm(v);
end
function v=aligned(v,ref)
v=unit(v);if v.'*ref<0,v=-v;end
end
function r=relative(a,b)
r=norm(a(:)-b(:))/max(norm(b(:)),realmin);
end
function x=blockData(H,b,ix)
if b==1,z=H(4,1:3,ix);elseif b==2,z=H(1:3,4,ix);else,z=H(4,4,ix);end
x=reshape(z,[],nnz(ix)).';
end
function v=rmsBand(f,x)
v=sqrt(trapz(f,sum(abs(x).^2,2))/(f(end)-f(1)));
end
function [co,conditions]=refitNumerators(Phi,H,W)
co=zeros(size(Phi,2),16);conditions=zeros(16,1);
for k=1:16
    A=Phi.*sqrt(W(:,k));v=H(:,k).*sqrt(W(:,k));A=[real(A);imag(A)];v=[real(v);imag(v)];
    scale=vecnorm(A);assert(all(scale>0));As=A./scale;conditions(k)=cond(As);
    co(:,k)=(As\v)./scale.';
end
end
function s=spaceDiagnostic(co,a,refs,target)
Y=[];U=[];
for m=1:3
    R=modalResidue(reshape(co(2*m-1,:),4,4),reshape(co(2*m,:),4,4),a(m,:),refs(m));
    Y=[Y,real(R),imag(R)];U=[U,real(R).',imag(R).']; %#ok<AGROW>
end
[uy,sy]=svd(Y);[uu,su]=svd(U);sy=diag(sy);su=diag(su);
R=modalResidue(reshape(co(2*target-1,:),4,4),reshape(co(2*target,:),4,4),a(target,:),refs(target));
s=struct('sy',sy,'su',su,'sigma4_over3_y',sy(4)/sy(3),'sigma4_over3_u',su(4)/su(3), ...
    'ny',uy(:,end),'nu',uu(:,end),'target',realRankOne(R));
end
function a=angleBetween(x,y)
a=acosd(min(1,abs(x(:).'*y(:))/norm(x)/norm(y)));
end
