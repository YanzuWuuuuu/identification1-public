function comparison=compare_frf_preprocessing
% Same transfer model/seeds, two FRF estimators, and three fixed objectives.
% Original records and previous results are read only. No OE refinement.
out=fileparts(mfilename('fullpath'));base=fileparts(out);
P=load(fullfile(out,'srivc_z0_mid_fixed_full.mat'),'result');r=P.result;
S=load(fullfile(base,'analysis_20261015','auto_identify_motor1234_z0_spectral.mat'),'train','validation');
select=S.train.f>=25 & S.train.f<=800;f=S.train.f(select);
Ht=reshape(S.train.H(:,:,select),16,[]).';Hv=reshape(S.validation.H(:,:,select),16,[]).';
cache=fullfile(out,'local_rational_comparison_data.mat');
if isfile(cache)
    L=load(cache,'localTrain','localValidation','localTrainChecks','localValidationChecks');
else
    [L.localTrain,L.localTrainChecks]=localRationalAt(r.training_f,r.training_U,r.training_Y,f,r);
    [L.localValidation,L.localValidationChecks]=localRationalAt(r.validation.f,r.validation.U,r.validation.Y,f,r);
    save(cache,'-struct','L');
end
frequency0=[15 125 137 281 335 359 493 678 685 705];
damping0=[.15 .01 .01 .02 .02 .02 .02 .005 .005 .005];
bounds=[5 24;115 132;132 145;260 300;320 348;348 375;475 510;665 682;682 695;695 725];
opts=struct('fixed_modes',1,'max_iterations',100,'score_tolerance',1e-6,'verbose',true);
sources={'welch','lrm'};weights={'uniform','inverse_magnitude','frequency4'};
models=struct('tag',{},'fit',{},'H',{},'weight',{});rows={};notches={};
% Include raw SRIVC without changing its original saved coefficients.
hraw=rawTransfer(r,f);row=metrics('raw_UY',hraw,Hv,L.localValidation,f,r,[]);
rows(end+1,:)=row;notches=[notches;notchMetrics('raw_UY',hraw,Hv,f)];
for e=1:2
    if e==1,target=Ht;else,target=L.localTrain;end
    for a=1:3
        tag=[sources{e} '_' weights{a}];
        switch weights{a}
            case 'uniform',W=ones(size(target));
            case 'inverse_magnitude'
                % Same objective weight as common/frf_weights.m, rho=1e-3.
                floorValue=max(1e-3*max(abs(target),[],1),eps);
                W=1./max(abs(target),floorValue);
            case 'frequency4',W=repmat((f/100).^4,1,16);
        end
        W=W/median(W(:)); % One global scalar, leaves the objective unchanged.
        fprintf('\nCASE %s\n',tag);
        fit=frf_srivc_fit(f,target,frequency0,damping0,bounds,W,opts);
        h=frfTransfer(fit,f);
        models(end+1)=struct('tag',tag,'fit',fit,'H',h,'weight',W); %#ok<AGROW>
        rows(end+1,:)=metrics(tag,h,Hv,L.localValidation,f,r,fit); %#ok<AGROW>
        notches=[notches;notchMetrics(tag,h,Hv,f)]; %#ok<AGROW>
    end
end
names={'method','converged','iterations','instrument_score','welch_complex_error_25_800_percent', ...
    'welch_complex_error_200_735_percent','lrm_complex_error_25_800_percent', ...
    'conditional_raw_output_error_percent','forced_raw_output_error_percent'};
summary=cell2table(rows,'VariableNames',names);
notchTable=cell2table(notches,'VariableNames',{'method','output','input','band_low_hz','band_high_hz', ...
    'reference_minimum_hz','model_minimum_hz','reference_minimum_db','model_minimum_db','band_complex_error_percent'});
comparison=struct('f',f,'welchTrain',Ht,'welchValidation',Hv, ...
    'lrmTrain',L.localTrain,'lrmValidation',L.localValidation, ...
    'lrmTrainChecks',L.localTrainChecks,'lrmValidationChecks',L.localValidationChecks, ...
    'rawH',hraw,'models',models,'summary',summary,'notches',notchTable, ...
    'definition','Same shared-pole transfer model and seeds; training-only weights; second-half validation', ...
    'lrm','Local complex equation-error rational fit; R=2, denominator degree=2, half-window=45 DFT lines; no closed-loop bias correction');
writetable(summary,fullfile(out,'frf_preprocessing_comparison.csv'));
writetable(notchTable,fullfile(out,'frf_preprocessing_notches.csv'));
save(fullfile(out,'frf_preprocessing_comparison.mat'),'comparison','-v7.3');
disp(summary);disp(notchTable);
fprintf('LRM scaled local design condition median/max: train %.4g/%.4g, validation %.4g/%.4g\n', ...
    median(L.localTrainChecks.condition(:)),max(L.localTrainChecks.condition(:)), ...
    median(L.localValidationChecks.condition(:)),max(L.localValidationChecks.condition(:)));
end

function [H,checks]=localRationalAt(f,U,Y,centers,r)
% Same local equation model as matlab/common/lrm_miso_core.m; centers only.
half=45;R=2;nu=size(U,2);ny=size(Y,2);p=(nu+1)*(R+1)+2;
H=complex(zeros(numel(centers),ny*nu));condition=zeros(numel(centers),ny);
ranks=condition;df=median(diff(f));centerError=zeros(size(centers));
for k=1:numel(centers)
    [~,center]=min(abs(f-centers(k)));centerError(k)=f(center)-centers(k);
    first=max(1,min(center-half,numel(f)-2*half));ix=first:first+2*half;
    x=(f(ix)-centers(k))/(half*df);
    for j=1:ny
        X=complex(zeros(numel(ix),p));output=Y(ix,j);
        for power=0:R
            X(:,power*nu+(1:nu))=U(ix,:).*x.^power;
            X(:,nu*(R+1)+power+1)=x.^power;
        end
        X(:,end-1:end)=[-x.*output,-x.^2.*output];
        scale=max(vecnorm(X),realmin);Xs=X./scale;singular=svd(Xs,'econ');
        condition(k,j)=singular(1)/singular(end);ranks(k,j)=sum(singular>1e-12*singular(1));
        theta=lsqminnorm(Xs,output,1e-12)./scale.';
        for u=1:nu,H(k,j+(u-1)*ny)=theta(u)*r.output_scale/r.input_scale(u);end
    end
end
checks=struct('condition',condition,'rank',ranks,'parameter_count',p, ...
    'window_points',2*half+1,'window_span_hz',2*half*df,'center_grid_error_hz',centerError);
assert(all(isfinite(H(:))) && all(ranks(:)==p),'LRM nonfinite result or deficient local design.');
end

function row=metrics(tag,h,Hv,Lv,f,r,fit)
b=f>=200&f<=735;
err=100*norm(h-Hv,'fro')/norm(Hv,'fro');errBand=100*norm(h(b,:)-Hv(b,:),'fro')/norm(Hv(b,:),'fro');
errL=100*norm(h-Lv,'fro')/norm(Lv,'fro');
if isempty(fit),hf=rawTransfer(r,r.validation.f);basefit=r.fit;
else,hf=frfTransfer(fit,r.validation.f);basefit=fit;end
forced=complex(zeros(size(r.validation.Y)));
for u=1:4
    forced=forced+hf(:,(u-1)*4+(1:4)).*(r.validation.U(:,u)*r.input_scale(u)/r.output_scale);
end
% Same validation boundary family for all transfer models, with their own poles.
[B,~,~,~,~,bc]=raw_srivc_design(r.validation.f,r.validation.U,basefit.a,basefit.reference_hz,true);
L=B(:,bc);w=(r.validation.f/100).^2;A=[real(L.*w);imag(L.*w)];
target=r.validation.Y-forced;yy=[real(target.*w);imag(target.*w)];sc=max(vecnorm(A),realmin);
boundary=((A./sc)\yy)./sc.';prediction=forced+L*boundary;den=norm(r.validation.Y.*w,'fro');
er=100*norm((prediction-r.validation.Y).*w,'fro')/den;
ef=100*norm((forced-r.validation.Y).*w,'fro')/den;
row={tag,basefit.converged,basefit.iterations,basefit.score,err,errBand,errL,er,ef};
end

function rows=notchMetrics(tag,h,Hv,f)
rows=cell(2,10);ii=f>=180&f<=260;fb=f(ii);channels=[1 6];
for j=1:2
    ref=Hv(ii,channels(j));model=h(ii,channels(j));[rv,ri]=min(abs(ref));[mv,mi]=min(abs(model));
    rows(j,:)={tag,j,j,180,260,fb(ri),fb(mi),20*log10(rv),20*log10(mv),100*norm(ref-model)/norm(ref)};
end
end

function h=frfTransfer(fit,f)
q=1i*f(:)./fit.reference_hz(:).';A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';
B=complex(zeros(numel(f),2*numel(fit.reference_hz)+1));
for m=1:numel(fit.reference_hz),B(:,2*m-1:2*m)=[ones(numel(f),1),q(:,m)]./A(:,m);end
B(:,end)=1;h=B*fit.coefficients;
end

function h=rawTransfer(r,f)
fit=r.fit;q=1i*f(:)./fit.reference_hz(:).';A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';
nm=numel(fit.reference_hz);h=complex(zeros(numel(f),16));nb=10;
for u=1:4
    value=complex(zeros(numel(f),4));
    for m=1:nm
        ix=(m-1)*nb+2*u-1;
        value=value+(ones(numel(f),1)*fit.coefficients(ix,:)+q(:,m)*fit.coefficients(ix+1,:))./A(:,m);
    end
    value=value+fit.coefficients(nm*nb+u,:);
    h(:,(u-1)*4+(1:4))=value*r.output_scale/r.input_scale(u);
end
end
