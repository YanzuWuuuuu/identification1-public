function decompose_antiresonance_error
% Exact amplitude/phase error identity; no new fitted parameters.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'weight_diagnosis_20261016');
S=load(fullfile(out,'weight_sweep.mat'),'audit');A=S.audit;
S=load(fullfile(out,'structure_diagnostics.mat'),'diagnostics');D=S.diagnostics;
tags={'lrm_p1_floor0.001','lrm_p1_floor0','lrm_p2_floor0.0001'};models=struct('tag',{},'H',{});
for k=1:3
    ii=find(strcmp({A.models.tag},tags{k}));models(k)=struct('tag',tags{k},'H',A.models(ii).H);
end
ii=find(strcmp({D.models.tag},'fixed_poles_profile_common_delay'));
models(4)=struct('tag',D.models(ii).tag,'H',D.models(ii).H);
rows={};repeat={};f=A.f;band=f>=180 & f<=260;
for channel=[1 6]
    ff=f(band);[~,im]=min(abs(A.lrmTrain(band,channel)));center=ff(im);
    bottom=band & abs(f-center)<=5;shoulders=band & ~bottom;
    sets={band,bottom,shoulders};names={'whole_180_260','bottom_plusminus5Hz','remaining_shoulders'};
    for refname={'welch','lrm'}
        if strcmp(refname{1},'welch'),reference=A.welchValidation(:,channel);training=A.welchTrain(:,channel);
        else,reference=A.lrmValidation(:,channel);training=A.lrmTrain(:,channel);end
        for region=1:3
            mask=sets{region};hr=reference(mask);normalizer=sum(abs(hr).^2);
            for k=1:numel(models)
                hm=models(k).H(mask,channel);amp=sum((abs(hm)-abs(hr)).^2);
                phase=sum(2*abs(hm).*abs(hr).*(1-cos(angle(hm.*conj(hr)))));
                total=sum(abs(hm-hr).^2);
                assert(abs(amp+phase-total)<1e-10*max(total,realmin));
                rows(end+1,:)={models(k).tag,refname{1},channel,names{region},center, ...
                    100*sqrt(amp/normalizer),100*sqrt(phase/normalizer),100*sqrt(total/normalizer), ...
                    100*amp/max(total,realmin)}; %#ok<AGROW>
            end
        end
        [a,ia]=min(abs(training(band)));[b,ib]=min(abs(reference(band)));
        repeat(end+1,:)={refname{1},channel,ff(ia),20*log10(a),ff(ib),20*log10(b), ...
            100*norm(training(band)-reference(band))/norm(reference(band))}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'tag','reference','channel','region','training_notch_center_hz', ...
    'amplitude_error_rms_percent','phase_error_rms_percent','complex_error_percent','amplitude_fraction_of_squared_error_percent'});
R=cell2table(repeat,'VariableNames',{'estimator','channel','training_min_hz','training_min_db', ...
    'validation_min_hz','validation_min_db','split_complex_difference_percent'});
writetable(T,fullfile(out,'exact_error_decomposition.csv'));writetable(R,fullfile(out,'notch_reference_repeatability.csv'));
disp(T(strcmp(T.region,'whole_180_260') & strcmp(T.reference,'lrm'),:));disp(R);
end
