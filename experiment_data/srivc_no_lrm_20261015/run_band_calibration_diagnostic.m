function campaign=run_band_calibration_diagnostic
% Auxiliary direct-target diagnostic; NOT the structured modal estimator.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'band_calibration_20261016');
if ~isfolder(out),mkdir(out);end
tags={'zm0p5','z0','z0p5','z1'};records=cell(1,4);results=cell(1,4);
for p=1:4
    a=load(fullfile(base,'structured_calibration_20261016',['calibration_' tags{p} '.mat']),'calib');
    records{p}=a.calib;
end
config=struct('training_band',[25 80],'delay_bound_ms',.25, ...
    'bands',[3 25;25 80;100 170;260 310;25 800], ...
    'method','Profiled real common null vectors and separable relative delays', ...
    'scope','Auxiliary direct-target diagnostic, not modal SRIVC or full paper algorithm', ...
    'split','Previously inspected same-record first/second halves; no fresh blind test');
for p=1:4
    results{p}=fitRecord(records{p},config.training_band,config.delay_bound_ms);
    fprintf('%s direct diagnostic: output %.5g -> %.5g, input %.5g -> %.5g\n', ...
        tags{p},results{p}.y.static.cost,results{p}.y.dynamic.cost,results{p}.u.static.cost,results{p}.u.dynamic.cost);
end
sens=cell(1,3);sens{1}=fitRecord(records{2},[25 60],.25);
sens{2}=fitRecord(records{2},[35 80],.25);sens{3}=fitRecord(records{2},[25 80],1);
campaign=struct('config',config,'records',{records},'results',{results},'sensitivity',{sens}, ...
    'created',char(datetime('now')),'matlab_version',version);
rows={};directionRows={};delayRows={};
for p=1:4
    c=records{p};r=results{p};
    methods=makeMethods(c,r);
    [rr,methods]=evaluateMethods(c,methods,config.bands,'own_position');rows=[rows;rr]; %#ok<AGROW>
    campaign.results{p}.methods=methods;
    for s=1:2
        names={'y','u'};side=names{s};v=r.(side);
        directionRows(end+1,:)={c.position_mm,side,angleDegrees(v.static.n,v.dynamic.n), ...
            v.dynamic.cost,v.static.cost,v.dynamic.exitflag,v.dynamic.firstorderopt,v.dynamic.at_bound,v.gradient_error}; %#ok<AGROW>
        for j=1:4
            delayRows(end+1,:)={c.position_mm,side,j,v.static.n(j),v.dynamic.n(j), ...
                v.dynamic.delay_ms(j)*1000,(max(v.dynamic.delay_ms)-v.dynamic.delay_ms(j))*1000}; %#ok<AGROW>
        end
    end
    if p~=2
        transferred=makeMethods(c,results{2});
        [rr,transferred]=evaluateMethods(c,transferred,config.bands,'z0_transfer');rows=[rows;rr]; %#ok<AGROW>
        campaign.results{p}.transfer_methods=transferred;
    end
end
for k=1:3
    label=sprintf('sensitivity_%d',k);mm=makeMethods(records{2},sens{k});
    [rr,mm]=evaluateMethods(records{2},mm,config.bands,label);rows=[rows;rr]; %#ok<AGROW>
    campaign.sensitivity{k}.methods=mm;
end
campaign.metrics=cell2table(rows,'VariableNames',{'position_mm','scope','method','band_lo','band_hi', ...
    'block','rms','B0_rms','over_B0','original_FRF_rms','percent_original_FRF'});
campaign.direction_metrics=cell2table(directionRows,'VariableNames',{'position_mm','side','static_dynamic_angle_deg', ...
    'dynamic_training_cost','static_training_cost','exitflag','firstorderopt','at_bound','gradient_error'});
campaign.parameters=cell2table(delayRows,'VariableNames',{'position_mm','side','channel','static_null','dynamic_null', ...
    'effective_relative_delay_us','causal_added_delay_us'});
writetable(campaign.metrics,fullfile(out,'band_calibration_metrics.csv'));
writetable(campaign.direction_metrics,fullfile(out,'band_calibration_directions.csv'));
writetable(campaign.parameters,fullfile(out,'band_calibration_parameters.csv'));
save(fullfile(out,'band_calibration_campaign.mat'),'campaign','-v7.3');
disp(campaign.direction_metrics);
fprintf('Saved %s\n',out);
end

function r=fitRecord(c,band,bound)
f=c.f(:);mask=f>=band(1)-1e-6 & f<=band(2)+1e-6;ff=f(mask);
H=reshape(c.H_training(mask,:).',4,4,[]);w=quadrature(ff);
normalizer=sum(w.*reshape(sum(abs(H).^2,[1 2]),[],1));
ell=c.fits.joint.delay_seconds;dy=(ell(:,1)-ell(1,1))*1000;du=(ell(1,:)-ell(1,1)).'*1000;
r=struct('position_mm',c.position_mm,'band',band,'bound_ms',bound);
r.y=fitSide(H,ff,w/normalizer,1,c.nominal.ny,dy,bound);
r.u=fitSide(H,ff,w/normalizer,-1,c.nominal.nu,du,bound);
end

function fit=fitSide(H,f,w,sgn,reference,oldDelay,bound)
nf=numel(f);A=complex(zeros(4,4,nf));
for k=1:nf
    if sgn==1,A(:,:,k)=H(:,:,k)*H(:,:,k)';else,A(:,:,k)=H(:,:,k)'*H(:,:,k);end
end
[v0,~,n0,Q0]=profile(zeros(3,1));
nom=reference(:)/norm(reference);nomCost=nom.'*Q0*nom;
assert(v0<=nomCost+1e-12,'Static spectral solution exceeds nominal training objective.');
static=struct('n',align(n0,reference),'delay_ms',zeros(4,1),'cost',v0,'eigenvalues',eig(Q0));
xp=[.013;-.021;.034];[~,ga]=profile(xp);gf=zeros(3,1);
for j=1:3,h=1e-5;x1=xp;x2=xp;x1(j)=x1(j)+h;x2(j)=x2(j)-h;gf(j)=(profile(x1)-profile(x2))/(2*h);end
gradientError=norm(ga-gf)/max([norm(ga),norm(gf),1e-9]);
assert(gradientError<5e-5,'Delay profile derivative check failed.');
opts=optimoptions('fmincon','Algorithm','sqp','Display','none','SpecifyObjectiveGradient',true, ...
    'MaxIterations',120,'MaxFunctionEvaluations',400,'OptimalityTolerance',1e-10, ...
    'StepTolerance',1e-11,'ConstraintTolerance',1e-10);
starts=[zeros(3,1),max(-bound,min(bound,oldDelay(2:4))),max(-bound,min(bound,-oldDelay(2:4)))];
trials=cell(1,3);
for t=1:3
    [x,value,flag,output]=fmincon(@profile,starts(:,t),[],[],[],[],-bound*ones(3,1),bound*ones(3,1),[],opts);
    [~,~,n,Q]=profile(x);
    trials{t}=struct('n',align(n,reference),'delay_ms',[0;x],'cost',value,'exitflag',flag, ...
        'firstorderopt',output.firstorderopt,'eigenvalues',eig(Q),'at_bound',any(abs(x)>.999*bound));
end
good=find(cellfun(@(v)v.exitflag>0 && v.firstorderopt<1e-7,trials));
if isempty(good),good=1:3;end
[~,j]=min(cellfun(@(v)v.cost,trials(good)));selected=trials{good(j)};
assert(selected.cost<=v0+1e-10,'Dynamic fit exceeds zero-delay feasible objective.');
fit=struct('static',static,'dynamic',selected,'trials',{trials},'gradient_error',gradientError, ...
    'nominal_training_cost',nomCost,'selection','Minimum training cost among converged candidates if available');
    function [value,g,n,Q]=profile(x)
        phase=exp(sgn*1i*2*pi*f(:).'*1e-3.*[0;x]);
        Ar=A.*reshape(phase,4,1,nf).*reshape(conj(phase),1,4,nf);
        Q=real(sum(Ar.*reshape(w,1,1,nf),3));Q=(Q+Q.')/2;
        [V,L]=eig(Q,'vector');[value,ix]=min(L);n=V(:,ix);
        if nargout>1
            g=zeros(3,1);
            for jj=2:4
                an=reshape(sum(Ar(jj,:,:).*reshape(n,1,4,1),2),[],1);
                g(jj-1)=2*real(sum(w.*(sgn*1i*2*pi*f*1e-3).*an))*n(jj);
            end
        end
    end
end

function methods=makeMethods(c,r)
ny=c.nominal.ny;nu=c.nominal.nu;ell=c.fits.joint.delay_seconds;
dy=ell(:,1)-ell(1,1);du=(ell(1,:)-ell(1,1)).';z=zeros(4,1);
codes={'B0','B1','B4','D0','D1','D2'};
ys={ny,c.fits.joint.Ny,c.fits.joint.Ny,r.y.static.n,r.y.dynamic.n,r.y.dynamic.n};
us={nu,c.fits.joint.Nu,c.fits.joint.Nu,r.u.static.n,r.u.dynamic.n,r.u.dynamic.n};
dys={z,z,dy,z,z,r.y.dynamic.delay_ms*1e-3};dus={z,z,du,z,z,r.u.dynamic.delay_ms*1e-3};
methods=struct([]);
for m=1:6
    item=struct('code',codes{m},'ny',ys{m},'nu',us{m},'dy',dys{m},'du',dus{m});
    if m==1,methods=repmat(item,1,6);end
    methods(m)=item;
end
end

function [rows,methods]=evaluateMethods(c,methods,bands,scope)
f=c.f(:);N=numel(f);H=reshape(c.H_validation.',4,4,[]);
P=c.nominal.Ty./vecnorm(c.nominal.Ty,2,2);Q=c.nominal.Tu./vecnorm(c.nominal.Tu,2,1);
rows={};values=zeros(6,size(bands,1),5);
for m=1:6
    v=methods(m);ny=v.ny(:);nu=v.nu(:);
    assert(abs(norm(ny)-1)<1e-10 && abs(norm(nu)-1)<1e-10);
    py=exp(-1i*2*pi*f.*(max(v.dy)-v.dy).');pu=exp(-1i*2*pi*f.*(max(v.du)-v.du).');
    assert(all(max(v.dy)-v.dy>=-1e-14) && all(max(v.du)-v.du>=-1e-14));
    left=[P;ny.'];right=[Q,nu];Z=complex(zeros(4,4,N));fullY=complex(zeros(N,4));fullU=fullY;
    for k=1:N
        corrected=diag(py(k,:))*H(:,:,k)*diag(pu(k,:));Z(:,:,k)=left*corrected*right;
        fullY(k,:)=ny.'*corrected;fullU(k,:)=(corrected*nu).';
    end
    blocks={reshape(Z(4,1:3,:),3,[]).',reshape(Z(1:3,4,:),3,[]).',reshape(Z(4,4,:),[],1),fullY,fullU};
    names={'logical_to_N','N_to_logical','N_to_N','all_inputs_to_N','N_to_all_outputs'};
    methods(m).H_validation=Z;methods(m).full_output_null=fullY;methods(m).full_input_null=fullU;
    methods(m).Ty_static=left;methods(m).Tu_static=right;
    for b=1:size(bands,1)
        mask=f>=bands(b,1)-1e-6 & f<=bands(b,2)+1e-6;
        original=rmsFreq(f(mask),c.H_validation(mask,:));
        for h=1:5
            val=rmsFreq(f(mask),blocks{h}(mask,:));values(m,b,h)=val;
            rows(end+1,:)={c.position_mm,scope,v.code,bands(b,1),bands(b,2),names{h}, ...
                val,values(1,b,h),val/values(1,b,h),original,100*val/original}; %#ok<AGROW>
        end
    end
end
assert(norm(methods(5).Ty_static-methods(6).Ty_static,'fro')==0);
assert(norm(methods(5).Tu_static-methods(6).Tu_static,'fro')==0);
end

function w=quadrature(f)
d=diff(f);w=[d(1);d(1:end-1)+d(2:end);d(end)]/2;
end
function value=rmsFreq(f,z)
value=sqrt(sum(quadrature(f).*sum(abs(z).^2,2))/(f(end)-f(1)));
end
function n=align(n,reference)
if n.'*reference<0,n=-n;end
end
function angle=angleDegrees(a,b)
angle=acosd(min(1,abs(a.'*b)/(norm(a)*norm(b))));
end
