function result=fit_retest_position(index)
% Actual additive SRIVC: all denominator coefficients are updated by IV.
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
old=fullfile(root,'experiment_data','srivc_no_lrm_20261015');addpath(old);
s=load(fullfile(out,'data_campaign.mat'),'campaign');D=s.campaign.records{index};
s=load(fullfile(old,'spring_identification_20261016','spring_selected.mat'),'result');parent=s.result;
f=D.f;H=D.trainingFRF;W=D.weights;delay=parent.delay_matrix_seconds;
phase=exp(2i*pi*f.*delay(:).');seed=parent.selected.fit;
boundsOld=[4 25;4 25;4 25;85 155;110 190;230 325;290 365;330 410;440 550;620 700;650 700;690 735;90 165;735 799;150 450];
fnOld=seed.frequency_hz;zOld=seed.damping;
% Training-only alternative: replace the old broad 264 Hz auxiliary term
% by a newly visible 650--675 Hz feature, keeping the same number of terms.
fnNew=[10.8;12.3;16.2;124;129.5;136.5;281.5;339.5;360.5;493;669;678.5;681.5;705.5;778];
boundsNew=[4 25;4 25;4 25;85 145;100 155;115 190;230 325;290 365;330 410;440 550;620 676;672 686;676 698;690 735;735 799];
zNew=[.05;.05;.04;.006;.005;.006;.004;.003;.003;.002;.008;.002;.002;.002;.015];
windows=[8 11.8;11 14;14.5 19;115 126;126 132;132 150;260 310;315 350;350 385;470 520;650 674;675 680;680 690;695 718;750 795];
curve=sqrt(sum(abs(H).^2,2)).*f.^2;
for m=[1 3 5:15]
    ii=find(f>=windows(m,1)&f<=windows(m,2));[~,j]=max(curve(ii));fnNew(m)=f(ii(j));
end
fnNew=min(max(fnNew,boundsNew(:,1)+.01),boundsNew(:,2)-.01);
% Separate actual high-frequency peaks; fixed overlapping windows can seed
% identical modes at some positions and make the numerator design singular.
hi=f>=650&f<=698;[~,peaks,~,prom]=findpeaks(log10(curve(hi)),f(hi),'MinPeakProminence',.02,'MinPeakDistance',1);
if numel(peaks)>=3
    [~,jj]=sort(prom,'descend');three=sort(peaks(jj(1:3)));
elseif numel(peaks)==2
    three=sort([665;peaks]);
else
    three=[665;678;683];
end
fnNew(11:13)=three;zNew(11:13)=.0015;
boundsNew(11:13,:)=[620,(three(1)+three(2))/2; (three(1)+three(2))/2,(three(2)+three(3))/2; (three(2)+three(3))/2,700];
opts=struct('fixed_modes',[],'max_iterations',220,'score_tolerance',1e-6,'verbose',false,'direct',true);
keep=setdiff(1:15,4);z14=zNew(keep);z14(end)=.0008;
starts={fnOld,fnNew,fnNew(keep),[fnNew(keep);264]};
damps={zOld,zNew,z14,[z14;.06]};
bounds={boundsOld,boundsNew,boundsNew(keep,:),[boundsNew(keep,:);200 320]};
trials=struct('fit',{},'H',{},'cost',{},'elapsed_seconds',{},'start_index',{});
diary(fullfile(out,['fit_z' D.tag '.log']));cleanup=onCleanup(@()diary('off'));
for k=1:numel(starts)
    tic;fit=frf_srivc_fit(f,H.*phase,starts{k},damps{k},bounds{k},W,opts);
    pred=fit.prediction./phase;
    trials(k)=struct('fit',fit,'H',pred,'cost',fit.weighted_cost,'elapsed_seconds',toc,'start_index',k); %#ok<AGROW>
    fprintf('z=%g trial%d SRIVC train%.4f%% score%.3g converged%d (%s), %dit, %.1fs\n',D.position_mm,k,100*fit.relative_error,fit.score,fit.converged,fit.reason,fit.iterations,trials(k).elapsed_seconds);
    disp([fit.frequency_hz,fit.damping]);
    save(fullfile(out,['trials_z' D.tag '.mat']),'trials','-v7.3');
end
% Choose training cost, with strict stationarity preferred when available.
valid=find(arrayfun(@(x)x.fit.converged,trials));
if isempty(valid),[~,pick]=min([trials.cost]);else,[~,j]=min([trials(valid).cost]);pick=valid(j);end
result=D;result.selected=trials(pick);result.trials=trials;result.delay_matrix_seconds=delay;
result.inner_converged=result.selected.fit.converged;
if result.inner_converged,result.status='converged';else,result.status='nonconverged';end
result.delay_status='fixed_old_z0_for_comparison';result.delay_source=fullfile(old,'spring_identification_20261016','spring_selected.mat');
result.selection='Training cost among strict stationary candidates if any; otherwise best training cost with nonconverged flag';
result.validation_scope='First/second halves of one new record, no independent repeat or certified closed-loop unbiasedness';
result.source_code=fileread([mfilename('fullpath') '.m']);result.srivc_source=fileread(fullfile(old,'frf_srivc_fit.m'));
save(fullfile(out,['model_z' D.tag '.mat']),'result','-v7.3');
end
