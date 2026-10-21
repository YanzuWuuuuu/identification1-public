function report=analyze_repair_residuals(repairResultPath,beforebase)
% Read-only model diagnostics; this function NEVER identifies another model.
% analyze_repair_residuals(repairResultPath)
% analyze_repair_residuals(repairResultPath,beforebase)
% beforebase: optional MAT path containing result, or a result structure.
% Defaults to the original seven-coordinate channel-delay SRIVC result.
% All comparisons use the original common f grid and training-only p=1 W.
% Reused half-record validation and estimator differences are NOT noise truth.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
if ~isfolder(out),mkdir(out);end
if nargin<1 || isempty(repairResultPath)
    repairResultPath=fullfile(out,'background_repair_selected.mat');
end
if nargin<2 || isempty(beforebase)
    beforebase=fullfile(base,'delay_identification_20261016','channel_delay_srivc.mat');
end
Q=load(fullfile(base,'delay_identification_20261016','common_delay_srivc.mat'),'result');common=Q.result;
Q=load(fullfile(base,'frf_preprocessing_comparison.mat'),'comparison');data=Q.comparison;
[before,beforeSource]=readResult(beforebase);
[repair,repairSource]=readResult(repairResultPath);
f=common.f(:);W=common.weights;nf=numel(f);nc=16;
assert(isequal(f,data.f(:)),'Preprocessing and common-model grids differ.');
assert(isequal(size(W),[nf,nc]) && isreal(W) && all(isfinite(W(:))) && all(W(:)>0), ...
    'Expected positive finite frequency-by-16 diagnostic quadratic weights.');
refs=struct('name',{'welch_train','welch_validation','lrm_train','lrm_validation'}, ...
    'estimator',{'Welch_H1','Welch_H1','LRM','LRM'}, ...
    'split',{'train','validation','train','validation'}, ...
    'H',{data.welchTrain,data.welchValidation,data.lrmTrain,data.lrmValidation});
for r=1:numel(refs)
    assert(isequal(size(refs(r).H),[nf,nc]) && all(isfinite(refs(r).H(:))), ...
        'Reference FRF layout or finite-value check failed.');
end
assert(relativeDifference(refs(3).H,common.trainingFRF)<1e-12 && ...
    relativeDifference(refs(2).H,common.welchValidation)<1e-12 && ...
    relativeDifference(refs(4).H,common.lrmValidation)<1e-12, ...
    'References must be identical to the saved training/validation data.');
models=struct('stage',{'before','repair'},'result',{before,repair},'H',{[],[]});
reconstruction=zeros(2,1);
for m=1:2
    q=models(m).result;
    assert(isequal(q.f(:),f),'Model and reference frequency grids differ.');
    models(m).H=transfer(q.selected.fit,f,q.delay_matrix_seconds);
    reconstruction(m)=relativeDifference(models(m).H,q.selected.H);
    assert(reconstruction(m)<1e-10,'Saved prediction disagrees with coefficients or delay layout.');
end
[~,stem]=fileparts(repairSource);
if isempty(stem),stem='repair_struct';end
stem=regexprep(stem,'[^A-Za-z0-9_-]','_');stem=stem(1:min(60,numel(stem)));
prefix=['diagnostics_' stem];

% The eight sub-bands are disjoint; 25-800 is an additional overall row.
bands=[25 800;25 50;50 100;100 180;180 260;260 550;550 650;650 735;735 800];
bandNames={'full_25_800','low_25_50','low_50_100','100_180', ...
    'antiresonance_180_260','260_550','550_650','650_735','upper_735_800'};
assert(all(sum(bandMasks(f,bands),2)==2), ...
    'Each frequency must belong to the full band and exactly one sub-band.');
masks=bandMasks(f,bands);metricRows={};
for m=1:2
    for r=1:4
        H=refs(r).H;R=models(m).H-H;
        amplitudeSquared=(abs(models(m).H)-abs(H)).^2;
        phaseSquared=2*abs(models(m).H).*abs(H).* ...
            (1-cos(angle(models(m).H.*conj(H))));
        splitIdentity=norm(abs(R).^2-amplitudeSquared-phaseSquared,'fro')/ ...
            max(norm(abs(R).^2,'fro'),realmin);
        assert(splitIdentity<1e-10,'Exact amplitude/phase squared-error identity failed.');
        totalResidual=sum(abs(R).^2,'all');totalWeighted=sum(W.*abs(R).^2,'all');
        totalReference=sum(abs(H).^2,'all');totalWeightedReference=sum(W.*abs(H).^2,'all');
        for b=1:size(bands,1)
            use=masks(:,b);
            for ch=0:nc
                if ch==0,cols=1:nc;y=0;u=0;else,cols=ch;y=mod(ch-1,4)+1;u=ceil(ch/4);end
                e=R(use,cols);h=H(use,cols);w=W(use,cols);
                ee=sum(abs(e).^2,'all');hh=sum(abs(h).^2,'all');
                we=sum(w.*abs(e).^2,'all');wh=sum(w.*abs(h).^2,'all');
                ae=sum(amplitudeSquared(use,cols),'all');pe=sum(phaseSquared(use,cols),'all');
                metricRows(end+1,:)={models(m).stage,refs(r).estimator,refs(r).split, ...
                    ch,y,u,bandNames{b},bands(b,1),bands(b,2),nnz(use),numel(e), ...
                    relativeRms(ee,hh),relativeRms(we,wh),relativeRms(ae,hh),relativeRms(pe,hh), ...
                    ee,hh,we,wh,percent(ee,totalResidual),percent(we,totalWeighted), ...
                    percent(hh,totalReference),percent(wh,totalWeightedReference)}; %#ok<AGROW>
            end
        end
    end
end
metricNames={'model','reference_estimator','split','channel','output','input','band', ...
    'fmin_hz','fmax_hz','frequency_count','complex_value_count', ...
    'relative_complex_error_percent','relative_weighted_complex_error_percent', ...
    'amplitude_component_rms_percent','phase_component_rms_percent', ...
    'absolute_residual_energy_recorded_units','reference_energy_recorded_units', ...
    'weighted_residual_loss','weighted_reference_energy', ...
    'residual_share_of_model_full_band_percent','weighted_loss_share_of_model_full_band_percent', ...
    'reference_share_of_full_band_percent','weighted_reference_share_of_full_band_percent'};
metrics=cell2table(metricRows,'VariableNames',metricNames);

% These discrepancies are observable reproducibility/estimator-agreement
% diagnostics, NOT a measurement-noise estimate or independent confidence CI.
pairs=[1 2;3 4;3 1;4 2];
pairNames={'Welch_train_vs_validation','LRM_train_vs_validation', ...
    'train_LRM_vs_Welch','validation_LRM_vs_Welch'};
agreementRows={};ratioRows={};weakChannels=[4 7]; % u1->y4 and u2->y3.
ratioFloorRelative=1e-3;
for p=1:4
    A=refs(pairs(p,1)).H;B=refs(pairs(p,2)).H;
    for b=1:size(bands,1)
        use=masks(:,b);
        for ch=0:nc
            if ch==0,cols=1:nc;y=0;u=0;else,cols=ch;y=mod(ch-1,4)+1;u=ceil(ch/4);end
            aa=A(use,cols);bb=B(use,cols);dd=aa-bb;
            de=sum(abs(dd).^2,'all');be=sum(abs(bb).^2,'all');
            se=sum(abs((aa+bb)/2).^2,'all');
            agreementRows(end+1,:)={pairNames{p},refs(pairs(p,1)).name,refs(pairs(p,2)).name, ...
                ch,y,u,bandNames{b},bands(b,1),bands(b,2),nnz(use), ...
                relativeRms(de,be),relativeRms(de,se),de,be}; %#ok<AGROW>
        end
        for ch=weakChannels
            floorValue=max(ratioFloorRelative*max(abs(B(:,ch))),realmin);
            valid=use & abs(B(:,ch))>=floorValue & abs(A(:,ch))>0;
            aa=A(valid,ch);bb=B(valid,ch);ratio=aa./bb;
            if isempty(ratio)
                gain=NaN;gainPhase=NaN;medianDb=NaN;lowDb=NaN;highDb=NaN;phaseMean=NaN;phaseSpread=NaN;
            else
                % A scalar complex LS ratio avoids silently dividing by notches.
                gainComplex=sum(conj(bb).*aa)/sum(abs(bb).^2);
                gain=abs(gainComplex);gainPhase=angle(gainComplex)*180/pi;
                db=20*log10(abs(ratio));medianDb=median(db);lowDb=quantileLocal(db,.1);highDb=quantileLocal(db,.9);
                phaseMean=angle(sum(ratio./abs(ratio)));
                phaseSpread=sqrt(mean(angle(ratio*exp(-1i*phaseMean)).^2))*180/pi;
                phaseMean=phaseMean*180/pi;
            end
            ratioRows(end+1,:)={pairNames{p},refs(pairs(p,1)).name,refs(pairs(p,2)).name, ...
                ch,mod(ch-1,4)+1,ceil(ch/4),bandNames{b},bands(b,1),bands(b,2), ...
                nnz(valid),nnz(use)-nnz(valid),floorValue,gain,gainPhase,medianDb,lowDb,highDb,phaseMean,phaseSpread}; %#ok<AGROW>
        end
    end
end
agreement=cell2table(agreementRows,'VariableNames',{'comparison','numerator_reference','denominator_reference', ...
    'channel','output','input','band','fmin_hz','fmax_hz','frequency_count', ...
    'difference_relative_to_denominator_percent','symmetric_relative_difference_percent', ...
    'absolute_difference_energy_recorded_units','denominator_energy_recorded_units'});
weakRatios=cell2table(ratioRows,'VariableNames',{'comparison','numerator_reference','denominator_reference', ...
    'channel','output','input','band','fmin_hz','fmax_hz','valid_ratio_count','excluded_near_zero_count', ...
    'denominator_floor_recorded_units','complex_LS_ratio_magnitude','complex_LS_ratio_phase_deg', ...
    'point_ratio_median_db','point_ratio_p10_db','point_ratio_p90_db', ...
    'point_ratio_circular_mean_phase_deg','point_ratio_circular_phase_rms_deg'});

% Prespecified 772-Hz candidate: inspect references AND before/after residuals.
% No smoothing, no pole fitting, and no automatic mode addition is performed.
% A local peak and its ratio to the 3-8-Hz side bands are merely diagnostics.
peakRows={};peakNumbers=NaN(3,4,nc+1,5);
sources={'experimental_FRF','before_residual','repair_residual'};
for s=1:3
    for r=1:4
        if s==1,V=refs(r).H;else,V=models(s-1).H-refs(r).H;end
        for ch=0:nc
            if ch==0,amp=sqrt(sum(abs(V).^2,2));y=0;u=0;else,amp=abs(V(:,ch));y=mod(ch-1,4)+1;u=ceil(ch/4);end
            search=f>=765 & f<=779;ids=find(search);[height,k]=max(amp(search));ix=ids(k);fp=f(ix);
            side=f>=730 & f<=800 & abs(f-fp)>=3 & abs(f-fp)<=8;
            sideMedian=median(amp(side));prominence=20*log10(max(height,realmin)/max(sideMedian,realmin));
            localMaximum=ix>1 && ix<nf && amp(ix)>=amp(ix-1) && amp(ix)>=amp(ix+1) && ...
                (amp(ix)>amp(ix-1) || amp(ix)>amp(ix+1));
            atSearchBoundary=k==1 || k==numel(ids);
            peakRows(end+1,:)={sources{s},refs(r).estimator,refs(r).split,ch,y,u,fp,height, ...
                20*log10(max(height,realmin)),sideMedian,prominence,localMaximum,atSearchBoundary}; %#ok<AGROW>
            peakNumbers(s,r,ch+1,:)=reshape([fp,prominence,localMaximum,atSearchBoundary,height],1,1,1,5);
        end
    end
end
peaks=cell2table(peakRows,'VariableNames',{'source','reference_estimator','split','channel','output','input', ...
    'candidate_peak_hz','candidate_amplitude_recorded_units','candidate_amplitude_db', ...
    'sideband_median_amplitude_recorded_units','peak_to_sideband_db','is_discrete_local_maximum', ...
    'is_search_boundary'});
stabilityRows={};
for s=1:3
    for ch=0:nc
        v=reshape(peakNumbers(s,:,ch+1,:),4,5);spread=max(v(:,1))-min(v(:,1));
        heuristic=spread<=1.5 && all(v(:,2)>=3) && all(v(:,3)>0) && ~any(v(:,4)>0) && ...
            all(v(:,1)>=768 & v(:,1)<=776);
        stabilityRows(end+1,:)={sources{s},ch,mod(ch-1,4)+1,ceil(ch/4), ...
            v(1,1),v(2,1),v(3,1),v(4,1),abs(v(1,1)-v(2,1)),abs(v(3,1)-v(4,1)), ...
            spread,min(v(:,2)),all(v(:,3)>0),any(v(:,4)>0),heuristic}; %#ok<AGROW>
        if ch==0,stabilityRows{end,3}=0;end
    end
end
peakStability=cell2table(stabilityRows,'VariableNames',{'source','channel','output','input', ...
    'welch_train_peak_hz','welch_validation_peak_hz','lrm_train_peak_hz','lrm_validation_peak_hz', ...
    'welch_split_peak_difference_hz','lrm_split_peak_difference_hz','all_estimates_peak_spread_hz', ...
    'minimum_peak_to_sideband_db','all_discrete_local_maxima','any_search_boundary', ...
    'repeated_candidate_heuristic_not_a_mode_test'});

files=struct('channel_band_metrics',fullfile(out,[prefix '_channel_band_metrics.csv']), ...
    'aggregate_band_metrics',fullfile(out,[prefix '_aggregate_band_metrics.csv']), ...
    'reference_agreement',fullfile(out,[prefix '_reference_agreement.csv']), ...
    'weak_channel_ratios',fullfile(out,[prefix '_weak_channel_ratios.csv']), ...
    'high_frequency_candidates',fullfile(out,[prefix '_high_frequency_candidates.csv']), ...
    'high_frequency_repetition',fullfile(out,[prefix '_high_frequency_repetition.csv']), ...
    'mat',fullfile(out,[prefix '.mat']));
writetable(metrics(metrics.channel>0,:),files.channel_band_metrics);
writetable(metrics(metrics.channel==0,:),files.aggregate_band_metrics);
writetable(agreement,files.reference_agreement);writetable(weakRatios,files.weak_channel_ratios);
writetable(peaks,files.high_frequency_candidates);writetable(peakStability,files.high_frequency_repetition);
% Keep reconstructed responses and fit metadata, not large optimizer caches.
reportModels=struct('stage',{},'H',{},'frequency_hz',{},'damping',{}, ...
    'delay_matrix_seconds',{},'inner_converged',{},'inner_score',{});
for m=1:2
    p=models(m).result.selected.fit;
    reportModels(m)=struct('stage',models(m).stage,'H',models(m).H, ...
        'frequency_hz',p.frequency_hz,'damping',p.damping, ...
        'delay_matrix_seconds',models(m).result.delay_matrix_seconds, ...
        'inner_converged',p.converged,'inner_score',p.score);
end
report=struct('f',f,'weights',W,'references',refs,'models',reportModels, ...
    'channel_and_band_metrics',metrics,'reference_agreement',agreement,'weak_channel_ratios',weakRatios, ...
    'high_frequency_candidates',peaks,'high_frequency_repetition',peakStability, ...
    'before_source',beforeSource,'repair_source',repairSource,'reconstruction_relative_errors',reconstruction, ...
    'files',files,'bands_hz',bands,'band_names',{bandNames}, ...
    'weight_definition','Original fixed common-result inverse-magnitude quadratic weights, identical for both models and both splits', ...
    'energy_definition','Sums of squared recorded-unit complex FRF differences; equal frequency grid, no dB subtraction', ...
    'contribution_definition','Percentages use each model/reference/split full-band all-channel loss; sub-bands are disjoint', ...
    'reference_difference_caveat','Split and Welch/LRM discrepancies are observed agreement proxies, not true measurement noise or confidence intervals', ...
    'weak_ratio_definition','A/B after excluding denominator amplitudes below 1e-3 of that denominator channel full-band maximum', ...
    'candidate_peak_definition','Search 765-779 Hz; side bands 3-8 Hz from each candidate; raw amplitude, no smoothing', ...
    'peak_heuristic_definition','Across Welch/LRM and both halves: spread<=1.5 Hz, all local maxima, no search-boundary peaks, all prominence>=3 dB, all peaks in 768-776 Hz; not a significance test', ...
    'validation_caveat','Both halves come from the same acquisition; validation has been reused during exploratory diagnosis', ...
    'created_utc',char(datetime('now','TimeZone','UTC')));
save(files.mat,'report','-v7.3');
disp(metrics(metrics.channel==0 & strcmp(metrics.split,'validation') & strcmp(metrics.band,'full_25_800'), ...
    {'model','reference_estimator','relative_complex_error_percent','relative_weighted_complex_error_percent'}));
disp(peakStability(peakStability.channel==0,:));
fprintf('Saved read-only repair diagnostics: %s\n',files.mat);
end

function [r,source]=readResult(value)
if isstruct(value)
    r=value;source='provided_result_struct';
else
    source=char(value);s=load(source,'result');
    assert(isfield(s,'result'),'MAT file must contain result.');r=s.result;
end
assert(isfield(r,'f') && isfield(r,'selected') && isfield(r,'delay_matrix_seconds') && isfield(r.selected,'fit') && ...
    isfield(r.selected,'H') && isequal(size(r.delay_matrix_seconds),[4,4]), ...
    'Result must contain f, selected.fit/H, and 4-by-4 delay_matrix_seconds.');
end

function H=transfer(p,f,ell)
q=1i*f(:)./p.reference_hz(:).';A=1+q.*p.a(:,1).'+q.^2.*p.a(:,2).';
nm=size(p.a,1);hasDirect=isfield(p,'direct') && p.direct;
B=complex(zeros(numel(f),2*nm+double(hasDirect)));
for k=1:nm,B(:,2*k-1)=1./A(:,k);B(:,2*k)=q(:,k)./A(:,k);end
if hasDirect,B(:,end)=1;end
assert(size(p.coefficients,1)==size(B,2) && size(p.coefficients,2)==16,'Unexpected FRF coefficient layout.');
H=(B*p.coefficients).*exp(-1i*2*pi*f(:).*ell(:).');
assert(all(isfinite(H(:))),'Nonfinite reconstructed transfer.');
end

function masks=bandMasks(f,bands)
masks=false(numel(f),size(bands,1));tol=1e-7;
for b=1:size(bands,1)
    if bands(b,2)==800
        masks(:,b)=f>=bands(b,1)-tol & f<=bands(b,2)+tol;
    else
        masks(:,b)=f>=bands(b,1)-tol & f<bands(b,2)-tol;
    end
end
end

function value=relativeRms(numerator,denominator)
if denominator>0,value=100*sqrt(max(0,numerator)/denominator);else,value=NaN;end
end

function value=percent(numerator,denominator)
if denominator>0,value=100*numerator/denominator;else,value=NaN;end
end

function value=relativeDifference(a,b)
assert(isequal(size(a),size(b)),'Compared matrices must have identical dimensions.');
value=norm(a-b,'fro')/max(norm(b,'fro'),realmin);
end

function q=quantileLocal(x,p)
x=sort(x(:));n=numel(x);
if n==0,q=NaN;return;end
t=1+(n-1)*p;lo=floor(t);hi=ceil(t);q=x(lo)+(t-lo)*(x(hi)-x(lo));
end
