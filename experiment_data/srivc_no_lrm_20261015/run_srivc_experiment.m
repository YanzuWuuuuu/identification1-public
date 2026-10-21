function result=run_srivc_experiment(positionMm,maxIterations,variant)
% Raw experimental spectra -> additive SRIVC, without any LRM or OE refinement.
if nargin<1,positionMm=0;end
if nargin<2,maxIterations=100;end
if nargin<3,variant='base';end
out=fileparts(mfilename('fullpath'));base=fileparts(out);
S=load(fullfile(base,sprintf('auto_identify_motor1234_z%g.mat',positionMm)));
t=S.u_Motor1(:,1);fs=1/median(diff(t));U=zeros(numel(t),4);Y=U;
for k=1:4
    U(:,k)=S.(sprintf('u_Motor%d',k))(:,2);
    Y(:,k)=S.(sprintf('pre_y_%d__0_',k-1))(:,2);
end
first=round(fs)+1;n=floor((numel(t)-first+1)/2);
it=first:first+n-1;iv=first+n:first+2*n-1;
uscale=std(U(it,:),0,1);yscale=1e-6;
Ut=fft((U(it,:)-mean(U(it,:)))./uscale,n)/n;
Yt=fft((Y(it,:)-mean(Y(it,:)))/yscale,n)/n;
Uv=fft((U(iv,:)-mean(U(iv,:)))./uscale,n)/n;
Yv=fft((Y(iv,:)-mean(Y(iv,:)))/yscale,n)/n;
f=(0:n-1)'*fs/n;all=find(f>=25 & f<=800);train=all(1:3:end);val=all;
% Rough resonance seeds, not nonlinear output-error fitted parameters.
frequency0=[15,281,335,359,493,678,685,705];
damping0=[.15,.02,.02,.02,.02,.005,.005,.005];
bounds=[5 24;260 300;320 348;348 375;475 510;665 682;682 695;695 725];
opts=struct('max_iterations',maxIterations,'verbose',true,'direct',true);
if startsWith(variant,'extended')
    frequency0=[15,125,137,281,335,359,493,678,685,705,772];
    damping0=[.10,.01,.01,.02,.02,.02,.02,.005,.005,.005,.01];
    bounds=[5 24;115 132;132 145;260 300;320 348;348 375;475 510;665 682;682 695;695 725;750 795];
end
if startsWith(variant,'mid')
    frequency0=[15,125,137,281,335,359,493,678,685,705];
    damping0=[.15,.01,.01,.02,.02,.02,.02,.005,.005,.005];
    bounds=[5 24;115 132;132 145;260 300;320 348;348 375;475 510;665 682;682 695;695 725];
end
if contains(variant,'fixed')
    opts.fixed_modes=1;
end
if contains(variant,'low10'),frequency0(1)=10;end
if contains(variant,'low20'),frequency0(1)=20;end
if contains(variant,'full'),train=all;end
if contains(variant,'shifted')
    offsets=.1*(bounds(:,2)-bounds(:,1)).*(-1).^(1:numel(frequency0)).';
    frequency0(2:end)=frequency0(2:end)+offsets(2:end).';
    damping0(2:end)=1.5*damping0(2:end);
end
fit=raw_srivc_fit(f(train),Ut(train,:),Yt(train,:),frequency0,damping0,bounds,opts);
[B,~,~,~,pc,bc]=raw_srivc_design(f(val),Uv(val,:),fit.a,fit.reference_hz,fit.direct);
forced=B(:,pc)*fit.coefficients(pc,:);w=(f(val)/100).^2;
L=B(:,bc);Ar=[real(L.*w);imag(L.*w)];target=Yv(val,:)-forced;
tr=[real(target.*w);imag(target.*w)];scale=max(vecnorm(Ar),realmin);
boundary=((Ar./scale)\tr)./scale.';predicted=forced+L*boundary;
bands=[25 100;100 200;200 315;315 385;385 550;550 650;650 735;735 800;25 800];
rows=zeros(size(bands,1),4);
for j=1:size(bands,1)
    id=f(val)>=bands(j,1)&f(val)<=bands(j,2);
    den=norm(Yv(val(id),:).*w(id),'fro');
    rows(j,:)=[bands(j,:),100*norm((predicted(id,:)-Yv(val(id),:)).*w(id),'fro')/den, ...
        100*norm((forced(id,:)-Yv(val(id),:)).*w(id),'fro')/den];
end
validation=struct('f',f(val),'U',Uv(val,:),'Y',Yv(val,:),'prediction',predicted, ...
    'forced_prediction',forced,'boundary_coefficients',boundary, ...
    'conditional_error',norm((predicted-Yv(val,:)).*w,'fro')/norm(Yv(val,:).*w,'fro'));
result=struct('fit',fit,'validation',validation,'position_mm',positionMm, ...
    'input_scale',uscale,'output_scale',yscale,'sample_rate_hz',fs, ...
    'segment_seconds',n/fs,'training_frequency_step_hz',(train(2)-train(1))*fs/n, ...
    'training_f',f(train),'training_U',Ut(train,:),'training_Y',Yt(train,:), ...
    'source_file',sprintf('auto_identify_motor1234_z%g.mat',positionMm));
tag=sprintf('z%g',positionMm);
if ~strcmp(variant,'base'),tag=[tag '_' variant];end
result.variant=variant;
save(fullfile(out,['srivc_' tag '.mat']),'result','-v7.3');
modalTable=table((1:numel(fit.frequency_hz)).',fit.frequency_hz,100*fit.damping, ...
    fit.initial_frequency_hz,fit.bounds_hz(:,1),fit.bounds_hz(:,2), ...
    'VariableNames',{'component','frequency_hz','damping_percent','initial_frequency_hz','lower_bound_hz','upper_bound_hz'});
writetable(modalTable,fullfile(out,['modes_' tag '.csv']));
writetable(array2table(rows,'VariableNames',{'fmin_hz','fmax_hz','conditional_error_percent','forced_only_error_percent'}),fullfile(out,['validation_' tag '.csv']));
writetable(array2table(fit.history,'VariableNames',{'iteration','weighted_error','instrument_score','accepted_fraction','ridge','moment_condition','regressor_instrument_gap'}),fullfile(out,['iteration_' tag '.csv']));
fprintf('FINAL: %s, converged=%d, iterations=%d, train %.4f%%, validation %.4f%%, score %.3g\n', ...
    fit.reason,fit.converged,fit.iterations,100*fit.relative_error,100*validation.conditional_error,fit.score);
disp(modalTable);disp(rows);
end
