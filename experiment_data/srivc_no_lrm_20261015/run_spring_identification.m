function campaign=run_spring_identification
% Fit explicitly observed suspension poles, replacing the old extrapolation bases.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'spring_identification_20261016');
S=load(fullfile(out,'spring_frf_data.mat'),'data');data=S.data;
S=load(fullfile(base,'error_repair_20261016','final_channel_delay_srivc.mat'),'result');old=S.result;seed=old.selected.fit;
ix=data.low_f>=3;f=[data.low_f(ix);data.upper.f];H=[data.low(3).train(ix,:);data.upper.lrmTrain];
V=[data.low(3).validation(ix,:);data.upper.lrmValidation];
% Preserve the old anti-notch floor and inverse-magnitude rule; quadrature
% prevents the finer low-frequency grid receiving five-fold artificial weight.
floorValue=max(1e-3*max(abs(data.upper.lrmTrain),[],1),eps);
df=[.1*ones(nnz(ix),1);.5*ones(numel(data.upper.f),1)];
W=df./max(abs(H),floorValue);W=W/median(W(:));
phase=exp(1i*2*pi*f.*old.delay_matrix_seconds(:).');target=H.*phase;
hi=setdiff(1:numel(seed.frequency_hz),seed.fixed_modes);
models=struct('tag',{},'fit',{},'H',{},'cost',{});
starts={[11;16],[10.6;16.5],[11;12.5;16],[10.5;13;16.5],[10.8;11.7;16],[10;11.3;16.5]};
for k=1:numel(starts)
    lo=starts{k};nl=numel(lo);fn=[lo;seed.frequency_hz(hi)];z=[.05*ones(nl,1);seed.damping(hi)];
    if nl==2,lb=[7 14;14 21];elseif k<=4,lb=[7 12;12 14.5;14.5 21];else,lb=[8 14;8 14;14 21];end
    bounds=[lb;seed.bounds_hz(hi,:)];
    opts=struct('fixed_modes',[],'max_iterations',260,'score_tolerance',1e-6,'verbose',false);
    fit=frf_srivc_fit(f,target,fn,z,bounds,W,opts);pred=fit.prediction./phase;
    cost=norm((pred-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
    tag=sprintf('spring_%d_start_%d',nl,k);models(end+1)=struct('tag',tag,'fit',fit,'H',pred,'cost',cost); %#ok<AGROW>
    fprintf('%s: train %.5f%%, score %.3g converged %d reason %s\n',tag,100*sqrt(cost),fit.score,fit.converged,fit.reason);
    disp([fit.frequency_hz(1:nl),fit.damping(1:nl)]);
    save(fullfile(out,'spring_fit_checkpoint.mat'),'models','f','W','-v7.3');
end
oldH=evaluateFit(seed,f).*exp(-1i*2*pi*f.*old.delay_matrix_seconds(:).');
rows=cell(numel(models)+1,9);rows(1,:)=row('old_extrapolation',oldH,seed);
for k=1:numel(models),rows(k+1,:)=row(models(k).tag,models(k).H,models(k).fit);end
summary=cell2table(rows,'VariableNames',{'model','converged','score','train_weighted_percent','validation_3_25_percent', ...
    'validation_25_50_percent','validation_25_800_percent','validation_100_170_percent','validation_200_800_percent'});
campaign=struct('f',f,'models',models,'summary',summary,'trainingFRF',H,'validationFRF',V,'weights',W,'frequency_quadrature',df, ...
    'reference_window',45,'delay_matrix_seconds',old.delay_matrix_seconds,'oldH',oldH);
save(fullfile(out,'spring_campaign.mat'),'campaign','-v7.3');writetable(summary,fullfile(out,'spring_candidate_summary.csv'));disp(summary);
valid=find(arrayfun(@(m)m.fit.converged,models));
if isempty(valid),warning('No spring fit meets score tolerance; candidates retained, none selected.');return;end
[~,j]=min([models(valid).cost]);chosen=models(valid(j));
result=struct('f',f,'trainingFRF',H,'validationFRF',V,'weights',W,'frequency_quadrature',df, ...
    'selected',chosen,'delay_matrix_seconds',old.delay_matrix_seconds,'delay_design',old.delay_design, ...
    'coordinates_ms',old.selected.coordinates_ms,'training_low_hz',3, ...
    'selection','Lowest training cost among converged candidates; validation reused exploratory split', ...
    'method','All spring and upper-band poles estimated by additive SRIVC; seven previous additive delays retained for controlled comparison');
save(fullfile(out,'spring_selected.mat'),'result','-v7.3');

    function rr=row(tag,p,fit)
        bands=[3 25;25 50;25 800;100 170;200 800];err=zeros(1,5);
        for b=1:5,jj=f>=bands(b,1)&f<bands(b,2);err(b)=100*norm(p(jj,:)-V(jj,:),'fro')/norm(V(jj,:),'fro');end
        rr=[{tag,fit.converged,fit.score,100*norm((p-H).*sqrt(W),'fro')/norm(H.*sqrt(W),'fro')},num2cell(err)];
    end
end

function H=evaluateFit(fit,f)
q=1i*f(:)./fit.reference_hz(:).';A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';
B=complex(zeros(numel(f),2*numel(fit.reference_hz)+1));
for m=1:numel(fit.reference_hz),B(:,2*m-1:2*m)=[ones(numel(f),1),q(:,m)]./A(:,m);end
B(:,end)=1;H=B*fit.coefficients;
end
