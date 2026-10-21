function diagnostics=diagnose_antiresonance_structure
% Focused causal checks beyond weights. Delayed profiles are DIAGNOSTICS,
% fixed-pole LS plus one searched delay, not relabelled as SRIVC solutions.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'weight_diagnosis_20261016');
S=load(fullfile(out,'weight_sweep.mat'),'audit');A=S.audit;f=A.f;
ix=find(strcmp({A.models.tag},'lrm_p2_floor0.0001'));p=A.models(ix).fit;
H=A.lrmTrain;W=1./max(abs(H),max(1e-4*max(abs(H),[],1),eps)).^2;W=W/median(W(:));
sw=sqrt(W);models=struct('tag',{},'H',{},'fit',{},'delay_seconds',{});rows={};
models(1)=struct('tag','baseline_p2_floor1e-4','H',p.prediction,'fit',p,'delay_seconds',0);
rows(1,:)=measure(models(1),A,W);
% Keep the same model structure/initial seeds; only the fixed background changes.
fn=[15 125 137 281 335 359 493 678 685 705];z=[.15 .01 .01 .02 .02 .02 .02 .005 .005 .005];
bounds=p.bounds_hz;opts=struct('fixed_modes',1,'max_iterations',180,'score_tolerance',1e-6,'verbose',false);
for low=[10 20]
    fn(1)=low;fit=frf_srivc_fit(f,H,fn,z,bounds,W,opts);
    model=struct('tag',sprintf('background_%gHz',low),'H',fit.prediction,'fit',fit,'delay_seconds',0);
    models(end+1)=model;rows(end+1,:)=measure(model,A,W); %#ok<AGROW>
end
% A second initial point checks numerical repeatability for the squared weight.
fn=[15 125 137 281 335 359 493 678 685 705];
fn(2:end)=fn(2:end)+(.1*(bounds(2:end,2)-bounds(2:end,1)).*(-1).^(2:10).').';
z(2:end)=1.5*z(2:end);fit=frf_srivc_fit(f,H,fn,z,bounds,W,opts);
models(end+1)=struct('tag','squared_weight_alternate_start','H',fit.prediction,'fit',fit,'delay_seconds',0);
rows(end+1,:)=measure(models(end),A,W);
% Fix ALL baseline poles: exact linear minimum under each p/rho.
B=basis(p,f);fixedRows={};
for power=1:2
    for rho=[1e-3 1e-4 0]
        ww=1./max(abs(H),max(rho*max(abs(H),[],1),realmin)).^power;ww=ww/median(ww(:));
        fit0=linearProfile(B,H,sqrt(ww),f,0);
        model=struct('tag',sprintf('fixed_p%d_floor%g',power,rho),'H',fit0.prediction,'fit',fit0,'delay_seconds',0);
        fixedRows(end+1,:)=measure(model,A,ww); %#ok<AGROW>
    end
end
% Multiplicative delay cannot change |H| directly; it can change the best
% fitted real numerators when amplitude and phase are fitted jointly.
grid=(0:10:1000)'*1e-6;profile=zeros(size(grid));
for k=1:numel(grid),fit0=linearProfile(B,H,sw,f,grid(k));profile(k)=fit0.weighted_cost;end
[~,best]=min(profile);lo=grid(max(best-1,1));hi=grid(min(best+1,numel(grid)));
objective=@(t)profileCost(B,H,sw,f,t);
if lo<hi,tau=fminbnd(objective,lo,hi,optimset('TolX',1e-10));else,tau=grid(best);end
candidates=[0;tau;grid(best);grid(end)];values=arrayfun(objective,candidates);[~,b]=min(values);tau=candidates(b);
fit0=linearProfile(B,H,sw,f,tau);
models(end+1)=struct('tag','fixed_poles_profile_common_delay','H',fit0.prediction,'fit',fit0,'delay_seconds',tau);
rows(end+1,:)=measure(models(end),A,W);
names={'tag','fit_gate','training_weighted_error_percent','welch_error_25_800_percent', ...
    'welch_error_200_735_percent','notch1_complex_error_percent','notch2_complex_error_percent', ...
    'notch1_min_db','notch2_min_db','delay_microseconds'};
summary=cell2table(rows,'VariableNames',names);fixed=cell2table(fixedRows,'VariableNames',names);
diagnostics=struct('models',models,'summary',summary,'fixed_poles',fixed, ...
    'delay_grid_seconds',grid,'delay_training_cost',profile,'f',f, ...
    'delay_definition','All baseline poles fixed; real linear numerators refitted; one common delay chosen ONLY by training weighted cost; not SRIVC nor a hardware delay estimate');
save(fullfile(out,'structure_diagnostics.mat'),'diagnostics','-v7.3');
writetable(summary,fullfile(out,'structure_diagnostics.csv'));writetable(fixed,fullfile(out,'fixed_poles_weight_sweep.csv'));
disp(summary);disp(fixed);
end

function row=measure(m,A,W)
f=A.f;h=m.H;iv=f>=200&f<=735;ib=f>=180&f<=260;
errors=zeros(1,2);deep=errors;
for j=1:2
    ch=1+5*(j-1);ref=A.welchValidation(ib,ch);pred=h(ib,ch);
    errors(j)=100*norm(pred-ref)/norm(ref);deep(j)=20*log10(min(abs(pred)));
end
row={m.tag,m.fit.converged,100*norm((h-A.lrmTrain).*sqrt(W),'fro')/norm(A.lrmTrain.*sqrt(W),'fro'), ...
    100*norm(h-A.welchValidation,'fro')/norm(A.welchValidation,'fro'), ...
    100*norm(h(iv,:)-A.welchValidation(iv,:),'fro')/norm(A.welchValidation(iv,:),'fro'),errors(1),errors(2),deep(1),deep(2),1e6*m.delay_seconds};
end

function cost=profileCost(B,H,sw,f,tau)
q=linearProfile(B,H,sw,f,tau);cost=q.weighted_cost;
end

function fit=linearProfile(B,H,sw,f,tau)
rot=exp(1i*2*pi*f*tau);C=zeros(size(B,2),size(H,2));
for ch=1:size(H,2)
    X=[real(B.*sw(:,ch));imag(B.*sw(:,ch))];y=H(:,ch).*rot.*sw(:,ch);y=[real(y);imag(y)];
    scale=max(vecnorm(X),realmin);C(:,ch)=((X./scale)\y)./scale.';
end
prediction=(B*C)./rot;cost=norm((prediction-H).*sw,'fro')^2;
fit=struct('coefficients',C,'prediction',prediction,'weighted_cost',cost,'converged',true, ...
    'method','Fixed-pole real linear least squares with profiled common delay','tau',tau);
end

function B=basis(p,f)
q=1i*f./p.reference_hz(:).';den=1+q.*p.a(:,1).'+q.^2.*p.a(:,2).';
B=complex(zeros(numel(f),2*size(p.a,1)+1));
for m=1:size(p.a,1),B(:,2*m-1)=1./den(:,m);B(:,2*m)=q(:,m)./den(:,m);end
B(:,end)=1;
end
