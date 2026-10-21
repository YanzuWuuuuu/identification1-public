function ref=run_recovered_reference_identification
% Pilot: known reconstructed reference -> [U;Y] local transfers, then
% fixed-pole/delay modal numerator fit to aggregated transfer moments.
% No single-bin DD' inverse; no assertion of finite-sample unbiasedness.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'full_coordinate_calibration_20261016');
S=load(fullfile(out,'recovered_reference_z0_candidate.mat'));T=load(fullfile(base,'structured_calibration_20261016','calibration_z0.mat'),'calib');
c=T.calib;j=c.fits.joint;f=c.f(:);nf=numel(f);ft=c.initial_free.selected.fit;nm=numel(ft.reference_hz);
it=c.initial_free.training_indices(1):c.initial_free.training_indices(2);
iv=c.initial_free.validation_indices(1):c.initial_free.validation_indices(2);
assert(numel(it)==numel(iv));n=numel(it);fs=S.fs;fr=(0:n-1)'*fs/n;
us=c.initial_free.input_scale;ys=c.initial_free.output_scale;ds=std(S.D_candidate(it,:));
vs=std(S.V_logic(it,:));assert(all(vs>0));
assert(norm(S.U(:,4)-S.D_candidate(:,4))<1e-10,'Fourth actuator is not direct external command.');
sets={it,iv};local=cell(1,2);
for split=1:2
    ix=sets{split};D=fft((S.D_candidate(ix,:)-mean(S.D_candidate(ix,:)))./ds)/n;
    % Regress the recorded feedback commands, then add the exact feedthrough.
    % Fitting total U directly makes its almost-identity rational denominator
    % poorly determined at high frequencies. No actual K is assumed here.
    Z=[(S.V_logic(ix,:)-mean(S.V_logic(ix,:)))./vs,zeros(n,1), ...
        (S.Y(ix,:)-mean(S.Y(ix,:)))/ys];Z=fft(Z)/n;
    [H,diagnostics]=localReference(fr,D,Z,f);
    Tvd=H(:,1:3,:);Tyd=H(:,5:8,:);Tud=complex(zeros(nf,4,4));
    feedbackMap=(S.Q_allocation.*vs)./us.';
    for k=1:nf,Tud(k,:,:)=diag(ds./us)+feedbackMap*squeeze(Tvd(k,:,:));end
    G=complex(zeros(nf,16));condU=zeros(nf,1);
    for k=1:nf
        U=squeeze(Tud(k,:,:));Y=squeeze(Tyd(k,:,:));condU(k)=cond(U);
        Hg=(Y/U)*ys./us;G(k,:)=Hg(:).';
    end
    assert(all(isfinite(G(:))));
    local{split}=struct('Tud',Tud,'Tyd',Tyd,'Tvd',Tvd,'FRF',G,'condition_Tud',condU,'diagnostics',diagnostics);
    fprintf('Reference split %d: local condition median/max %.3g %.3g; Tud median/max %.3g %.3g\n', ...
        split,median(diagnostics.condition(:)),max(diagnostics.condition(:)),median(condU),max(condU));
end

% Fixed pole, delay and order isolate the change of the statistical entrance.
refs=ft.reference_hz(:);a=ft.a;refs(1:3)=j.frequency_hz;a(1:3,:)=[2*j.damping,ones(3,1)];
q=1i*f./refs.';den=1+q.*a(:,1).'+q.^2.*a(:,2).';Phi=complex(zeros(nf,2*nm+1));
for m=1:nm,Phi(:,2*m-1)=1./den(:,m);Phi(:,2*m)=q(:,m)./den(:,m);end
Phi(:,end)=1;ell=j.delay_seconds;
% Same frozen original inverse-magnitude weight for both FRF refits.
weights=c.initial_free.weights;
directCo=fitFRF(Phi,c.H_training,weights,ell,f);
referenceCo=fitFRF(Phi,local{1}.FRF,weights,ell,f);

% Additional fit to Tyd-G*Tud. The right-hand side consists of local
% aggregates over 91 DFT bins, not one noisy Fourier outer product.
% Weights are fixed from the old training model's predicted Tyd.
Tu=local{1}.Tud;Ty=local{1}.Tyd;np=size(Phi,2);momentCo=zeros(np,16);momentCondition=zeros(4,1);
oldGn=reshape(j.H.',4,4,[]).*reshape(us,1,4,1)/ys;
pred=complex(zeros(nf,4,4));
for k=1:nf,pred(k,:,:)=oldGn(:,:,k)*squeeze(Tu(k,:,:));end
df=c.initial_free.frequency_quadrature(:);beta=ones(nf,1);beta(f<25)=.03;
momentWeights=zeros(nf,4);
for y=1:4
    magn=sqrt(sum(abs(squeeze(pred(:,y,:))).^2,2));
    floorValue=1e-3*max(magn);momentWeights(:,y)=df.*beta./max(magn,floorValue);
    X=complex(zeros(nf*4,np*4));rhs=reshape(squeeze(Ty(:,y,:)),[],1);
    for u=1:4
        phase=exp(-2i*pi*f*ell(y,u));pu=Phi.*phase;
        for r=1:4,X((r-1)*nf+(1:nf),(u-1)*np+(1:np))=pu.*Tu(:,u,r);end
    end
    w=sqrt(repmat(momentWeights(:,y),4,1));Ar=X.*w;br=rhs.*w;
    Ar=[real(Ar);imag(Ar)];br=[real(br);imag(br)];sc=vecnorm(Ar);
    assert(all(isfinite(sc)) && all(sc>0) && all(isfinite(br)));
    As=Ar./sc;assert(rank(As)==size(As,2),'Moment design is not full column rank.');
    momentCondition(y)=cond(As);coef=(As\br)./sc.';
    for u=1:4,momentCo(:,y+(u-1)*4)=coef((u-1)*np+(1:np))*ys/us(u);end
end
codes={'previous_structured','fixed_poles_direct_FRF','fixed_poles_reference_FRF','fixed_poles_reference_moments'};
models=cell(1,4);models{1}=j.H;
cos={directCo,referenceCo,momentCo};
for m=1:3,models{m+1}=(Phi*cos{m}).*exp(-2i*pi*f*ell(:).');end
bands=[3 25;25 80;100 170;260 310;25 800];names={'3-25','25-80','100-170','260-310','25-800'};
rows={};frontRows={};
for b=1:size(bands,1)
    ix=f>=bands(b,1) & f<=bands(b,2)+1e-9;fb=f(ix);
    for split=1:2
        direct=c.H_training;if split==2,direct=c.H_validation;end
        frontRows(end+1,:)={split,names{b},rmsBand(fb,local{split}.FRF(ix,:)-direct(ix,:))/rmsBand(fb,direct(ix,:)), ...
            median(local{split}.condition_Tud(ix)),max(local{split}.condition_Tud(ix))}; %#ok<AGROW>
    end
    for m=1:4
        Gn=reshape(models{m}.',4,4,[]).*reshape(us,1,4,1)/ys;
        values=zeros(1,2);
        for split=1:2
            residual=complex(zeros(nnz(ix),4,4)); kk=find(ix);
            for z=1:numel(kk),k=kk(z);residual(z,:,:)=squeeze(local{split}.Tyd(k,:,:))-Gn(:,:,k)*squeeze(local{split}.Tud(k,:,:));end
            rr=reshape(residual,nnz(ix),[]);yy=reshape(local{split}.Tyd(ix,:,:),nnz(ix),[]);
            values(split)=rmsBand(fb,rr)/rmsBand(fb,yy);
        end
        rows(end+1,:)={codes{m},names{b},values(1),values(2), ...
            rmsBand(fb,models{m}(ix,:)-local{2}.FRF(ix,:))/rmsBand(fb,local{2}.FRF(ix,:)), ...
            rmsBand(fb,models{m}(ix,:)-c.H_validation(ix,:))/rmsBand(fb,c.H_validation(ix,:))}; %#ok<AGROW>
    end
end
metrics=cell2table(rows,'VariableNames',{'model','band_hz','train_reference_transfer_error', ...
    'check_reference_transfer_error','check_reference_FRF_error','check_old_FRF_error'});
frontend=cell2table(frontRows,'VariableNames',{'split','band_hz','reference_vs_direct_difference', ...
    'median_reference_input_condition','max_reference_input_condition'});
ref=struct('f',f,'local',{local},'models',{models},'model_codes',{codes},'metrics',metrics,'frontend',frontend, ...
    'direct_coefficients',directCo,'reference_coefficients',referenceCo,'moment_coefficients',momentCo, ...
    'moment_condition',momentCondition,'moment_weights',momentWeights,'reference_hz',refs,'a',a,'delay',ell, ...
    'input_scale',us,'output_scale',ys,'reference_scale',ds,'logic_scale',vs, ...
    'feedback_allocation',S.Q_allocation,'training_indices',[it(1),it(end)], ...
    'check_indices',[iv(1),iv(end)],'source_code',fileread([mfilename('fullpath') '.m']), ...
    'scope','Reconstructed d; local 3 Hz rational approximation; same-record check reused; fixed poles/delays and free numerators; no asymptotic or closed-loop control claim');
save(fullfile(out,'reference_identification.mat'),'ref','-v7.3');
writetable(metrics,fullfile(out,'reference_identification_metrics.csv'));
writetable(frontend,fullfile(out,'reference_frontend_comparison.csv'));
disp(frontend);disp(metrics);fprintf('Moment condition numbers: ');fprintf('%.4g ',momentCondition);fprintf('\n');
end

function [H,C]=localReference(f,D,Z,centers)
half=45;df=median(diff(f));p=17;H=complex(zeros(numel(centers),8,4));conditions=zeros(numel(centers),7);ranks=conditions;
outputs=[1:3,5:8];
for k=1:numel(centers)
    [~,center]=min(abs(f-centers(k)));imin=2;imax=floor(numel(f)/2);
    if centers(k)>=25,imin=find(f>=25,1);imax=find(f<=800,1,'last');end
    first=max(imin,min(center-half,imax-2*half));ix=first:first+2*half;
    x=(f(ix)-centers(k))/(half*df);
    for j=1:numel(outputs)
        y=Z(ix,outputs(j));X=complex(zeros(numel(ix),p));
        for power=0:2,X(:,power*4+(1:4))=D(ix,:).*x.^power;X(:,13+power)=x.^power;end
        X(:,16:17)=[-x.*y,-x.^2.*y];sc=max(vecnorm(X),realmin);Xs=X./sc;sv=svd(Xs,'econ');
        conditions(k,j)=sv(1)/max(sv(end),realmin);ranks(k,j)=sum(sv>1e-12*sv(1));
        theta=lsqminnorm(Xs,y,1e-12)./sc.';
        assert(ranks(k,j)==p && all(isfinite(theta)),'Reference local rank/finite check failed.');
        H(k,outputs(j),:)=theta(1:4);
    end
    if mod(k,400)==0,fprintf('Reference local regressions %d/%d\n',k,numel(centers));end
end
C=struct('condition',conditions,'rank',ranks,'parameter_count',p,'window_points',91,'window_span_hz',90*df, ...
    'known_direct_fourth_input_row',true);
end
function co=fitFRF(Phi,H,W,ell,f)
co=zeros(size(Phi,2),16);H=H.*exp(2i*pi*f*ell(:).');
for ch=1:16
    A=Phi.*sqrt(W(:,ch));b=H(:,ch).*sqrt(W(:,ch));A=[real(A);imag(A)];b=[real(b);imag(b)];sc=vecnorm(A);
    co(:,ch)=((A./sc)\b)./sc.';
end
end
function v=rmsBand(f,x)
v=sqrt(trapz(f,sum(abs(x).^2,2))/(f(end)-f(1)));
end
