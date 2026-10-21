function checks=plot_srivc_error_repair(repairResultPath,referenceEstimator)
% Compare all 16 FRFs and their errors against the SAME validation reference.
% Optional result: f, selected.fit/H, delay_matrix_seconds, label.
% referenceEstimator: 'welch' (default) or 'lrm'; proxy uses the same method.
% Residuals are evaluated on the original bins, never from display curves.
if nargin<2 || isempty(referenceEstimator),referenceEstimator='welch';end
referenceEstimator=lower(char(string(referenceEstimator)));
assert(any(strcmp(referenceEstimator,{'welch','lrm'})),'referenceEstimator must be welch or lrm.');
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
if ~isfolder(out),mkdir(out);end
D=load(fullfile(base,'delay_identification_20261016','common_delay_srivc.mat'),'result');common=D.result;
D=load(fullfile(base,'delay_identification_20261016','channel_delay_srivc.mat'),'result');before=D.result;
D=load(fullfile(base,'frf_preprocessing_comparison.mat'),'comparison');original=D.comparison;
f=common.f(:);
if strcmp(referenceEstimator,'welch')
    Hreference=common.welchValidation;Htrain=original.welchTrain;
    originalReference=original.welchValidation;referenceLabel='Welch H1';referenceKey='Welch';
else
    Hreference=common.lrmValidation;Htrain=original.lrmTrain;
    originalReference=original.lrmValidation;referenceLabel='LRM';referenceKey='LRM';
end
assert(isequal(f,original.f(:)),'Train/validation reference grids differ.');
assert(norm(Hreference-originalReference,'fro')<=1e-12*max(norm(Hreference,'fro'),realmin), ...
    'The repeatability proxy must use the same validation reference and estimator.');
referenceLegend=sprintf('Experimental FRF (validation %s)',referenceLabel);
assert(isequal(size(Hreference),[numel(f),16]) && isequal(size(Htrain),size(Hreference)), ...
    'Expected frequency-by-16 FRFs with output index varying fastest.');
assert(all(isfinite(Hreference(:))) && all(isfinite(Htrain(:))),'Nonfinite experimental FRF data.');
hasRepair=nargin>0 && ~isempty(repairResultPath);
models={before};names={'Existing channel-delay SRIVC'};prefix='baseline';repairLabel='';
if hasRepair
    R=load(repairResultPath,'result');after=R.result;models{2}=after;
    if isfield(after,'label'),repairLabel=char(string(after.label));else,repairLabel='Revised SRIVC model';end
    names={'Before repair','After repair'};
    [~,stem]=fileparts(repairResultPath);stem=regexprep(stem,'[^A-Za-z0-9_-]','_');
    prefix=['repair_' stem(1:min(numel(stem),70))];
end
if strcmp(referenceEstimator,'lrm'),prefix=['lrm_' prefix];end
nmodels=numel(models);blue=[.03 .25 .45];orange=[.78 .34 .08];gray=[.65 .65 .65];
if hasRepair,colors={orange,blue};else,colors={blue};end

% Dense evaluation is used only for FRF display, not for numerical errors.
fd=logspace(log10(min(f)),log10(max(f)),8000).';
for m=1:nmodels
    fit=models{m}.selected.fit;
    for j=1:numel(fit.frequency_hz)
        fc=fit.frequency_hz(j);span=max(30*fc*fit.damping(j),1);
        low=max(min(f),fc-span);high=min(max(f),fc+span);
        if high>low,fd=[fd;linspace(low,high,1200).'];end %#ok<AGROW>
    end
end
fd=unique([min(f);max(f);fd]);
Hmodel=cell(1,nmodels);Hdense=Hmodel;complexResidual=Hmodel;magnitudeError=Hmodel;phaseError=Hmodel;
reconstructionError=zeros(nmodels,1);invalidPhaseBins=zeros(nmodels,1);
for m=1:nmodels
    r=models{m};ell=r.delay_matrix_seconds(:).';
    Hsaved=transfer(r.selected.fit,r.f(:),ell);
    assert(isequal(size(Hsaved),size(r.selected.H)),'Saved prediction dimensions disagree.');
    reconstructionError(m)=norm(Hsaved-r.selected.H,'fro')/max(norm(r.selected.H,'fro'),realmin);
    assert(reconstructionError(m)<1e-10,'Coefficient/delay reconstruction disagrees with saved selected.H.');
    Hmodel{m}=transfer(r.selected.fit,f,ell);Hdense{m}=transfer(r.selected.fit,fd,ell);
    assert(all(isfinite(Hmodel{m}(:))) && all(isfinite(Hdense{m}(:))),'Nonfinite evaluated model FRF.');
    complexResidual{m}=Hmodel{m}-Hreference;
    valid=abs(Hreference)>0 & abs(Hmodel{m})>0;
    magnitudeError{m}=NaN(size(Hreference));phaseError{m}=NaN(size(Hreference));
    magnitudeError{m}(valid)=20*log10(abs(Hmodel{m}(valid))./abs(Hreference(valid)));
    % Principal phase error, deliberately NOT unwrapped or branch-aligned.
    phaseError{m}(valid)=angle(Hmodel{m}(valid).*conj(Hreference(valid)))*180/pi;
    invalidPhaseBins(m)=nnz(~valid);
    assert(all(isfinite(magnitudeError{m}(valid))) && all(isfinite(phaseError{m}(valid))), ...
        'Unexpected nonfinite error at nonzero FRF bins.');
end
repeatability=Htrain-Hreference; % Same chosen estimator, separate 30 s halves.
files={};

for kind={'complex_residual','magnitude_error','phase_error','frf_magnitude','frf_phase','frf_phase_independent'}
    name=kind{1};fig=figure('Visible','off','Color','w','Position',[45 45 1460 1040]);
    cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');handles=gobjects(0);
            if strcmp(name,'complex_residual')
                handles(end+1)=semilogx(ax,f,toDb(repeatability(:,ch)),':','Color',gray,'LineWidth',.9); %#ok<AGROW>
                for m=1:nmodels
                    handles(end+1)=semilogx(ax,f,toDb(complexResidual{m}(:,ch)), ...
                        'Color',colors{m},'LineWidth',.85); %#ok<AGROW>
                end
            elseif strcmp(name,'magnitude_error') || strcmp(name,'phase_error')
                yline(ax,0,':','Color',gray,'HandleVisibility','off');
                for m=1:nmodels
                    if strcmp(name,'magnitude_error'),v=magnitudeError{m}(:,ch);else,v=phaseError{m}(:,ch);end
                    handles(end+1)=semilogx(ax,f,v,'Color',colors{m},'LineWidth',.85); %#ok<AGROW>
                end
                if strcmp(name,'phase_error'),ylim(ax,[-180 180]);yticks(ax,[-180 -90 0 90 180]);end
            elseif strcmp(name,'frf_magnitude')
                handles(end+1)=semilogx(ax,f,toDb(Hreference(:,ch)),'Color',gray,'LineWidth',.7); %#ok<AGROW>
                for m=1:nmodels
                    handles(end+1)=semilogx(ax,fd,toDb(Hdense{m}(:,ch)),'Color',colors{m},'LineWidth',1.05); %#ok<AGROW>
                end
            else
                phaseDense=zeros(numel(fd),nmodels);
                for m=1:nmodels,phaseDense(:,m)=unwrap(angle(Hdense{m}(:,ch)))*180/pi;end
                if strcmp(name,'frf_phase')
                    anchor=phaseDense(:,end);
                    for m=1:nmodels-1
                        phaseDense(:,m)=phaseDense(:,m)+360*round((anchor(1)-phaseDense(1,m))/360);
                    end
                    referencePhase=angle(Hreference(:,ch))*180/pi;
                    atReference=interp1(fd,anchor,f,'linear');
                    referencePhase=referencePhase+360*round((atReference-referencePhase)/360);
                else
                    referencePhase=unwrap(angle(Hreference(:,ch)))*180/pi;
                end
                handles(end+1)=semilogx(ax,f,referencePhase,'.','Color',gray,'MarkerSize',3); %#ok<AGROW>
                for m=1:nmodels
                    handles(end+1)=semilogx(ax,fd,phaseDense(:,m),'Color',colors{m},'LineWidth',1.05); %#ok<AGROW>
                end
            end
            styleAxes(ax,y,u,[min(f) max(f)]);
            if y==1 && u==1,legendHandles=handles;end
        end
    end
    switch name
        case 'complex_residual'
            labels=[{sprintf('Train-validation %s difference (repeatability proxy)',referenceLabel)},names];
            axisLabel='Complex residual magnitude (dB, recorded output / motor-command unit)';
            mainTitle='Per-channel complex FRF residuals';
            note={'Model residual: H_{model}-H_{validation}, calculated at original bins before magnitude conversion.', ...
                sprintf('Gray proxy compares first/second-half %s FRFs; it includes drift and estimation variability, not only measurement noise.',referenceLabel)};
        case 'magnitude_error'
            labels=names;axisLabel='Magnitude error (dB)';mainTitle='Per-channel FRF magnitude error';
            note='20 log10(|H_{model}| / |H_{validation}|), at original bins. Zero indicates equal magnitudes.';
        case 'phase_error'
            labels=names;axisLabel='Principal phase error (deg)';mainTitle='Per-channel FRF phase error';
            note='angle(H_{model} conj(H_{validation})) at original bins. No unwrap, branch alignment, or smoothing.';
        case 'frf_magnitude'
            labels=[{referenceLegend},names];
            axisLabel='FRF magnitude (dB, recorded output / motor-command unit)';mainTitle='Full measured-boundary FRF comparison';
            note='Dense model evaluation includes all per-input/per-output additive delay phases. Experimental FRF remains an estimate.';
        case 'frf_phase'
            labels=[{referenceLegend},names];axisLabel='Unwrapped phase (deg)';mainTitle='Full measured-boundary FRF phase';
            note={'Gray samples select the primary model branch using integer 360-deg shifts. Each earlier model gets at most one constant whole-curve shift.', ...
                'Display only: reported phase errors use untouched complex FRFs. Independent-unwrap diagnostic is also saved.'};
        otherwise
            labels=[{referenceLegend},names];axisLabel='Independently unwrapped phase (deg)';mainTitle='Phase display diagnostic';
            note='Every experimental/model curve is unwrapped independently; no branch alignment or smoothing.';
    end
    lg=legend(legendHandles,labels,'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
    xlabel(tl,'Frequency (Hz)','FontSize',12);ylabel(tl,axisLabel,'FontSize',12);
    title(tl,sprintf('%s (%s reference)',mainTitle,referenceLabel),'FontWeight','normal');subtitle(tl,note,'FontSize',10);
    path=fullfile(out,[prefix '_' name '.png']);exportgraphics(fig,path,'Resolution',200);files{end+1}=path; %#ok<AGROW>
    if strcmp(name,'frf_magnitude') || strcmp(name,'frf_phase')
        path=fullfile(out,[prefix '_' name '.pdf']);exportgraphics(fig,path,'ContentType','vector');files{end+1}=path; %#ok<AGROW>
    end
    clear cleanup;
end

% Preserve frequency-domain errors for numerical diagnosis, without changing
% training or validation estimates. Each model is judged against one reference.
bands=[25 800;180 260;200 735;650 735];rows={};
for m=1:nmodels
    for y=1:4
        for u=1:4
            ch=y+4*(u-1);
            for b=1:size(bands,1)
                ii=f>=bands(b,1) & f<=bands(b,2);ref=Hreference(ii,ch);den=norm(ref);
                relative=100*norm(complexResidual{m}(ii,ch))/max(den,realmin);
                proxy=100*norm(repeatability(ii,ch))/max(den,realmin);
                me=magnitudeError{m}(ii,ch);pe=phaseError{m}(ii,ch);
                rows(end+1,:)={names{m},y,u,bands(b,1),bands(b,2),nnz(ii),relative, ...
                    sqrt(mean(me.^2,'omitnan')),sqrt(mean(pe.^2,'omitnan')),proxy}; %#ok<AGROW>
            end
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'model','output_channel','input_channel','fmin_hz','fmax_hz', ...
    'frequency_count',sprintf('relative_complex_difference_to_%s_percent',referenceKey),'rms_magnitude_error_db', ...
    'rms_principal_phase_error_deg','train_validation_proxy_relative_difference_percent'});
writetable(metrics,fullfile(out,[prefix '_per_channel_errors.csv']));
checks=struct('has_repair',hasRepair,'repair_label',repairLabel,'output_directory',out, ...
    'reconstruction_relative_errors',reconstructionError,'dense_frequency_count',numel(fd), ...
    'undefined_phase_or_magnitude_bin_count',invalidPhaseBins,'model_names',{names}, ...
    'reference_estimator',referenceEstimator,'reference',sprintf('Same second-half %s for all model errors',referenceLabel), ...
    'proxy',sprintf('First-half minus second-half %s FRF; repeatability proxy, not a measurement-noise truth',referenceLabel), ...
    'phase_error','Principal angle of model*conj(reference), no unwrap', ...
    'channel_order','column = output + 4*(input-1)','files',{files});
save(fullfile(out,[prefix '_error_curves.mat']),'f','Hreference','Hmodel','complexResidual', ...
    'magnitudeError','phaseError','repeatability','metrics','checks','referenceEstimator','-v7.3');
disp(checks);
end

function H=transfer(fit,f,delaySeconds)
f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if fit.direct,B(:,end)=1;end
P=B*fit.coefficients;assert(size(P,2)==16 && numel(delaySeconds)==16,'Expected 16 delayed transfer functions.');
H=exp(-1i*2*pi*f.*delaySeconds(:).').*P;
end

function styleAxes(ax,y,u,range)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on','LineWidth',.65, ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');xlim(ax,range);xticks(ax,[25 100 300 800]);
title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
if y<4,xticklabels(ax,{});end
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
