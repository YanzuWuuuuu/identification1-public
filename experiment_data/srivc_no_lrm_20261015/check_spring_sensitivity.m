function campaign=check_spring_sensitivity
% Bounded checks of low-FRF preprocessing and one additional low-band pole.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'spring_identification_20261016');
S=load(fullfile(out,'spring_selected.mat'),'result');parent=S.result;
S=load(fullfile(out,'spring_frf_data.mat'),'data');D=S.data;
f=parent.f;seed=parent.selected.fit;nl=nnz(f<25);phase=exp(1i*2*pi*f.*parent.delay_matrix_seconds(:).');
models=struct('tag',{},'fit',{},'H',{},'trainingFRF',{},'weights',{},'cost',{});rows={};
for k=1:5
    H=parent.trainingFRF;fn=seed.frequency_hz;z=seed.damping;bounds=seed.bounds_hz;
    if k==1,H(1:nl,:)=D.low(2).train(D.low_f>=3,:);tag='lrm_half30';
    elseif k==2,H(1:nl,:)=interp1(D.welchTrain.f,reshape(D.welchTrain.H,16,[]).',f(1:nl));tag='welch_low_interpolated';
    else
        extras=[8 14 21];extra=extras(k-2);fn=[fn;extra];z=[z;.08];bounds=[bounds;4 24.5];tag=sprintf('four_low_start_%g',extra);
    end
    fl=max(1e-3*max(abs(D.upper.lrmTrain),[],1),eps);W=parent.frequency_quadrature./max(abs(H),fl);W=W/median(W(:));
    opts=struct('fixed_modes',[],'max_iterations',300,'score_tolerance',1e-6,'verbose',false);
    fit=frf_srivc_fit(f,H.*phase,fn,z,bounds,W,opts);p=fit.prediction./phase;cost=norm((p-H).*sqrt(W),'fro')^2/norm(H.*sqrt(W),'fro')^2;
    models(end+1)=struct('tag',tag,'fit',fit,'H',p,'trainingFRF',H,'weights',W,'cost',cost); %#ok<AGROW>
    lo=f<25;hi=~lo;
    rows(end+1,:)={tag,fit.converged,fit.score,100*sqrt(cost), ...
        100*norm(p(lo,:)-parent.validationFRF(lo,:),'fro')/norm(parent.validationFRF(lo,:),'fro'), ...
        100*norm(p(hi,:)-parent.validationFRF(hi,:),'fro')/norm(parent.validationFRF(hi,:),'fro')}; %#ok<AGROW>
    fprintf('%s: score %.3g converged %d train %.3f\n',tag,fit.score,fit.converged,100*sqrt(cost));
    disp([fit.frequency_hz(fit.frequency_hz<25),fit.damping(fit.frequency_hz<25)]);
end
summary=cell2table(rows,'VariableNames',{'model','converged','score','weighted_train_percent','validation_lrm45_low_percent','validation_upper_percent'});
campaign=struct('f',f,'models',models,'summary',summary);save(fullfile(out,'spring_sensitivity.mat'),'campaign','-v7.3');
writetable(summary,fullfile(out,'spring_sensitivity.csv'));disp(summary);
end
