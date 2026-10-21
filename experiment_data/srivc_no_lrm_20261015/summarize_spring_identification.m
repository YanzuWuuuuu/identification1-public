function verification=summarize_spring_identification
% Common-band assessment, reconstruction checks, and final compact model.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'spring_identification_20261016');
S=load(fullfile(out,'spring_selected.mat'),'result');result=S.result;
S=load(fullfile(base,'error_repair_20261016','final_channel_delay_srivc.mat'),'result');old=S.result;
S=load(fullfile(out,'spring_frf_data.mat'),'data');data=S.data;
S=load(fullfile(out,'spring_validation_refit.mat'),'spring_val_refit');vf=S.spring_val_refit;
nm=numel(result.selected.fit.frequency_hz);
tableModes=table((1:nm)',result.selected.fit.frequency_hz,result.selected.fit.damping, ...
    'VariableNames',{'component','natural_frequency_hz','damping_ratio'});
writetable(tableModes,fullfile(out,'final_modal_parameters.csv'));
lowTable=table(result.selected.fit.frequency_hz(1:3),vf.frequency_hz(1:3), ...
    result.selected.fit.damping(1:3),vf.damping(1:3), ...
    'VariableNames',{'training_frequency_hz','second_half_refit_frequency_hz','training_damping','second_half_refit_damping'});
writetable(lowTable,fullfile(out,'spring_split_parameter_check.csv'));
models={old,result};tags={'old_25plus_extrapolated','joint_spring_model'};
bands=[3 25;3 8;8 13;13 18;18 25;25 50;25 800;100 170;200 800];rows={};
for estimator=1:2
    if estimator==1
        f=result.f;ref=result.validationFRF;name='LRM';
    else
        % Preserve the exact upper-band bins used in the previous comparison;
        % legacy selection omitted the nominal 25 Hz line due to fs rounding.
        fw=data.welchValidation.f;hw=reshape(data.welchValidation.H,16,[]).';
        iw=round(fw,6)>=3 & round(fw,6)<25;
        f=[fw(iw);data.upper.f];ref=[hw(iw,:);data.upper.welchValidation];name='Welch';
    end
    for m=1:2
        pred=evaluate(models{m},f);
        for b=1:size(bands,1)
            fq=round(f,6);ix=fq>=bands(b,1)&fq<bands(b,2);
            if bands(b,2)==800,ix=fq>=bands(b,1)&fq<=800;end
            err=100*norm(pred(ix,:)-ref(ix,:),'fro')/norm(ref(ix,:),'fro');
            rows(end+1,:)={tags{m},name,bands(b,1),bands(b,2),nnz(ix),err}; %#ok<AGROW>
        end
    end
end
bandTable=cell2table(rows,'VariableNames',{'model','validation_reference','band_start_hz','band_end_hz','frequency_count','relative_complex_error_percent'});
writetable(bandTable,fullfile(out,'final_band_metrics.csv'));
% Same original training-LRM-defined notch masks as the prior campaign.
f=data.upper.f;H=data.upper.lrmValidation;itr=f>=180&f<=260;fb=f(itr);notches={};
for m=1:2
    pred=evaluate(models{m},f);
    for ch=[1 6]
        [~,mi]=min(abs(data.upper.lrmTrain(itr,ch)));fc=fb(mi);im=itr&abs(f-fc)<=5;
        notches(end+1,:)={tags{m},ch,fc, ...
            100*norm(pred(itr,ch)-H(itr,ch))/norm(H(itr,ch)), ...
            100*norm(pred(im,ch)-H(im,ch))/norm(H(im,ch))}; %#ok<AGROW>
    end
end
notchTable=cell2table(notches,'VariableNames',{'model','flattened_channel','training_notch_center_hz','band_180_260_error_percent','center_pm5_error_percent'});
writetable(notchTable,fullfile(out,'final_notch_comparison.csv'));
reconstruction=norm(evaluate(result,result.f)-result.selected.H,'fro')/norm(result.selected.H,'fro');
ell=result.delay_matrix_seconds;additivity=max(abs(ell-ell(:,1)-ell(1,:)+ell(1,1)),[],'all');
% Check the outer profile gradient using fresh inner SRIVC solutions on two
% directions. Coordinates are milliseconds: 2e-4 is a 0.2 microsecond step.
E=result.delay_design;x=result.selected.coordinates_ms;fit=result.selected.fit;W=result.weights;
op=struct('fixed_modes',[],'max_iterations',250,'score_tolerance',1e-7,'verbose',false);
dirs=zeros(7,2);dirs(1,1)=1;dirs([2 5],2)=[1;-1]/sqrt(2);step=2e-4;fd=zeros(2,1);score=zeros(2,2);
for d=1:2
    costs=zeros(2,1);
    for signIndex=1:2
        sgn=2*signIndex-3;phase=exp(1i*2*pi*result.f.*(E*(x+sgn*step*dirs(:,d))*1e-3).');
        ff=frf_srivc_fit(result.f,result.trainingFRF.*phase,fit.frequency_hz,fit.damping,fit.bounds_hz,W,op);
        assert(ff.converged,'Profile gradient check has a nonstationary inner fit.');score(d,signIndex)=ff.score;
        costs(signIndex)=norm((result.trainingFRF.*phase-ff.prediction).*sqrt(W),'fro')^2/norm(result.trainingFRF.*sqrt(W),'fro')^2;
    end
    fd(d)=(costs(2)-costs(1))/(2*step);
end
analytic=dirs'*result.selected.gradient;gradientError=norm(fd-analytic,inf);
verification=struct('reconstruction_relative_error',reconstruction,'delay_additivity_seconds',additivity, ...
    'inner_score',fit.score,'outer_gradient_inf',norm(result.selected.gradient,inf), ...
    'outer_exitflag',result.optimizer_exitflag,'selected_matches_optimizer',result.selected_matches_optimizer, ...
    'profile_finite_difference',fd,'profile_analytic_directions',analytic,'profile_gradient_absolute_error',gradientError, ...
    'gradient_check_inner_scores',score,'low_split_refit_converged',vf.converged, ...
    'bands',bandTable,'notches',notchTable,'low_parameters',lowTable);
assert(reconstruction<1e-10 && additivity<1e-12 && fit.converged && gradientError<1e-6,'Final model verification failed.');
save(fullfile(out,'final_verification.mat'),'verification');
% Full search history remains in spring_delay_refined.mat.
if isfield(result,'trace_cache'),result=rmfield(result,'trace_cache');end
result.audit_file=fullfile(out,'final_verification.mat');save(fullfile(out,'spring_selected.mat'),'result','-v7.3');
disp(lowTable);disp(bandTable);disp(notchTable);
fprintf('Reconstruction %.3g; delay-additivity %.3g s; profile gradient error %.3g\n',reconstruction,additivity,gradientError);
end

function H=evaluate(r,f)
fit=r.selected.fit;f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';B=complex(zeros(numel(f),2*nm+1));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
B(:,end)=1;H=(B*fit.coefficients).*exp(-1i*2*pi*f.*r.delay_matrix_seconds(:).');
end
