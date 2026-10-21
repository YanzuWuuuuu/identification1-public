function checks=plot_spring_identification
% Display a saved joint spring/flexible-mode fit; never run identification.
% Original complex FRFs and the 30 s record split are preserved unchanged.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'spring_identification_20261016');
dataPath=fullfile(out,'spring_frf_data.mat');resultPath=fullfile(out,'spring_selected.mat');
oldPath=fullfile(base,'error_repair_20261016','final_channel_delay_srivc.mat');
for p={dataPath,resultPath,oldPath}
    if ~isfile(p{1})
        error('plot_spring_identification:NotReady', ...
            'Required saved data/model is not ready: %s. No figures generated.',p{1});
    end
end
S=load(dataPath,'data');data=S.data;S=load(resultPath,'result');joint=S.result;
S=load(oldPath,'result');old=S.result;
required={'f','selected','delay_matrix_seconds','training_low_hz','validationFRF'};
assert(all(isfield(joint,required)),'spring_selected.mat does not contain the agreed result fields.');
assert(isfield(joint.selected,'fit') && isfield(joint.selected,'H'),'Selected model is not ready.');
k=find([data.low.half_window]==45);assert(numel(k)==1,'Expected one low-frequency half-window=45 estimate.');
fr=[data.low_f(:);data.upper.f(:)];
Htrain=[data.low(k).train;data.upper.lrmTrain];Href=[data.low(k).validation;data.upper.lrmValidation];
assert(all(diff(fr)>0) && isequal(size(Htrain),size(Href)) && isequal(size(Href),[numel(fr),16]), ...
    'Combined LRM reference must be sorted Nfrequency-by-16.');
fitLow=joint.training_low_hz;assert(isscalar(fitLow) && isfinite(fitLow) && fitLow>0 && fitLow<25);
% Check that the saved model used the same half-window-45 validation reference.
[present,idx]=ismembertol(joint.f(:),fr,1e-10,'DataScale',1);
assert(all(present),'Saved fitted frequency bins do not match the FRF data.');
referenceAgreement=norm(joint.validationFRF-Href(idx,:),'fro')/max(norm(joint.validationFRF,'fro'),realmin);
assert(referenceAgreement<1e-10,'The saved validation reference is not the selected half-window-45 FRF.');
models={old,joint};names={'Old model (below 25 Hz: extrapolation)','Joint low-frequency model'};
orange=[.78 .34 .08];blue=[.03 .25 .45];gray=[.65 .65 .65];
fd=logspace(log10(fitLow),log10(800),8000).';
for m=1:2
    fit=models{m}.selected.fit;
    for j=1:numel(fit.frequency_hz)
        fc=fit.frequency_hz(j);span=max(30*fc*fit.damping(j),1);
        lo=max(fitLow,fc-span);hi=min(800,fc+span);
        if hi>lo,fd=[fd;linspace(lo,hi,1200).'];end %#ok<AGROW>
    end
end
fd=unique([fitLow;25;50;800;fd]);Hmodel=cell(2,1);Hdense=Hmodel;reconstruction=zeros(2,1);
for m=1:2
    r=models{m};ell=r.delay_matrix_seconds(:).';
    saved=transfer(r.selected.fit,r.f(:),ell);
    reconstruction(m)=norm(saved-r.selected.H,'fro')/max(norm(r.selected.H,'fro'),realmin);
    assert(reconstruction(m)<1e-10,'Model reconstruction differs from saved selected.H.');
    Hmodel{m}=transfer(r.selected.fit,fr,ell);Hdense{m}=transfer(r.selected.fit,fd,ell);
    assert(all(isfinite(Hdense{m}(:))) && all(isfinite(Hmodel{m}(:))),'Nonfinite model FRF.');
end
used=fr>=fitLow & fr<=800;assert(all(isfinite(Href(used,:)),'all'),'Nonfinite LRM reference inside fitted band.');
files={};
specs={ ...
    'spring_low_magnitude','magnitude',[fitLow 50]; ...
    'spring_low_phase','phase',[fitLow 50]; ...
    'spring_full_magnitude','magnitude',[fitLow 800]; ...
    'spring_upper_band_comparison','magnitude',[25 800]; ...
    'spring_complex_residual','residual',[fitLow 800]};
for s=1:size(specs,1)
    name=specs{s,1};kind=specs{s,2};range=specs{s,3};
    ir=fr>=range(1) & fr<=range(2);id=fd>=range(1) & fd<=range(2);
    fig=figure('Visible','off','Color','w','Position',[40 40 1460 1040]);
    cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');
            if strcmp(kind,'phase')
                phases=zeros(nnz(id),2);
                for m=1:2,phases(:,m)=unwrap(angle(Hdense{m}(id,ch)))*180/pi;end
                phases(:,1)=phases(:,1)+360*round((phases(1,2)-phases(1,1))/360);
                experimental=unwrap(angle(Href(ir,ch)))*180/pi;
                % One constant whole-curve branch shift, never pointwise
                % alignment of experimental phases to the fitted model.
                experimental=experimental+360*round((phases(1,2)-experimental(1))/360);
                h1=semilogx(ax,fr(ir),experimental,'.','Color',gray,'MarkerSize',3);
                h2=semilogx(ax,fd(id),phases(:,1),'--','Color',orange,'LineWidth',1);
                h3=semilogx(ax,fd(id),phases(:,2),'Color',blue,'LineWidth',1.1);
            elseif strcmp(kind,'residual')
                h1=semilogx(ax,fr(ir),toDb(Href(ir,ch)),'Color',gray,'LineWidth',.7);
                h2=semilogx(ax,fr(ir),toDb(Hmodel{1}(ir,ch)-Href(ir,ch)),'--','Color',orange,'LineWidth',.85);
                h3=semilogx(ax,fr(ir),toDb(Hmodel{2}(ir,ch)-Href(ir,ch)),'Color',blue,'LineWidth',.85);
            else
                h1=semilogx(ax,fr(ir),toDb(Href(ir,ch)),'Color',gray,'LineWidth',.7);
                h2=semilogx(ax,fd(id),toDb(Hdense{1}(id,ch)),'--','Color',orange,'LineWidth',1);
                h3=semilogx(ax,fd(id),toDb(Hdense{2}(id,ch)),'Color',blue,'LineWidth',1.1);
            end
            matrixAxes(ax,y,u,range);
            if range(1)<25 && range(2)>25,xline(ax,25,':','Color',[.75 .75 .75],'HandleVisibility','off');end
            if y==1 && u==1,handles=[h1 h2 h3];end
        end
    end
    if strcmp(kind,'residual')
        labels={'Validation LRM magnitude','Old-model complex residual','Joint-model complex residual'};
        axisLabel='Magnitude (dB, recorded output / motor-command unit)';
        mainTitle='Per-channel complex residuals to the same validation FRF';
        note={sprintf('New fit begins at %g Hz. Old-model behavior below 25 Hz is extrapolation.',fitLow), ...
            'Residuals are complex differences at identical original bins. The validation record is the second half of the same experiment.'};
    elseif strcmp(kind,'phase')
        labels=[{'Validation LRM'},names];axisLabel='Unwrapped phase (deg)';mainTitle='Low-frequency phase including spring dynamics';
        note={sprintf('New fit begins at %g Hz; old model below 25 Hz is extrapolation. Individual channel delays are included.',fitLow), ...
            'Each curve is unwrapped independently; only one constant whole-curve 360-deg branch shift is allowed. No smoothing.'};
    else
        labels=[{'Validation LRM'},names];axisLabel='Magnitude (dB, recorded output / motor-command unit)';
        if range(1)>=25
            mainTitle='Retained response above 25 Hz';
        elseif range(2)<=50
            mainTitle='Low-frequency FRF including spring dynamics';
        else
            mainTitle='Joint low-frequency and flexible-response model';
        end
        note={sprintf('New fit begins at %g Hz. Vertical line: old 25 Hz fitting limit; lower old-model curve is extrapolation.',fitLow), ...
            'Low-frequency reference: half-window-45 LRM. First/second 30 s are the fitting/validation split of one experiment.'};
        if range(1)>=25
            note={'Both models are evaluated against the same upper-band validation LRM; this panel checks preservation of the earlier fit.', ...
                'First/second 30 s are the fitting/validation split of one experiment; improvement is not assumed.'};
        end
    end
    lg=legend(handles,labels,'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
    xlabel(tl,'Frequency (Hz)','FontSize',12);ylabel(tl,axisLabel,'FontSize',12);
    title(tl,mainTitle,'FontWeight','normal');subtitle(tl,note,'FontSize',10);
    files=[files,exportFigure(fig,out,name)]; %#ok<AGROW>
    clear cleanup;
end

% A compact quality check covers the unfit range down to 0.5 Hz separately.
% Singular values describe the FRF, whereas normalized input-spectrum
% condition numbers diagnose excitation. These are not interchangeable.
q=fr>=.5 & fr<=50;fq=fr(q);sv=NaN(nnz(q),4);split=NaN(nnz(q),1);ii=find(q);
for j=1:numel(ii)
    ht=Htrain(ii(j),:);hv=Href(ii(j),:);
    if all(isfinite(hv)),sv(j,:)=svd(reshape(hv,4,4)).';end
    if all(isfinite([ht hv])) && norm(hv)>0,split(j)=100*norm(ht-hv)/norm(hv);end
end
wt=data.welchTrain;wv=data.welchValidation;
assert(all(isfield(wt,{'f','condition','coherence'})) && all(isfield(wv,{'f','condition','coherence'})), ...
    'Complete Welch diagnostic structures are required.');
it=wt.f>=.5 & wt.f<=50;iv=wv.f>=.5 & wv.f<=50;
fig=figure('Visible','off','Color','w','Position',[60 60 1330 880]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);semilogx(ax,fq,20*log10(max(sv,realmin)),'LineWidth',1.1);
qualityAxes(ax,fitLow);ylabel(ax,'FRF singular value (dB, recorded units)');
title(ax,'Validation LRM: four FRF singular values','FontWeight','normal');
legend(ax,{'sigma_1','sigma_2','sigma_3','sigma_4'},'Location','best','Box','off');
ax=nexttile(tl);loglog(ax,wt.f(it),wt.condition(it),'--','Color',orange,'LineWidth',1.05);hold(ax,'on');
loglog(ax,wv.f(iv),wv.condition(iv),'Color',blue,'LineWidth',1.1);
qualityAxes(ax,fitLow);ylabel(ax,'Normalized input-spectrum condition number');
title(ax,'Welch input-spectrum conditioning','FontWeight','normal');
legend(ax,{'First 30 s','Second 30 s'},'Location','best','Box','off');
ax=nexttile(tl);semilogx(ax,fq,split,'Color',blue,'LineWidth',1.05);
qualityAxes(ax,fitLow);ylabel(ax,'Relative complex FRF difference (%)');
title(ax,'LRM split-record repeatability proxy','FontWeight','normal');
ax=nexttile(tl);semilogx(ax,wv.f(iv),wv.coherence(:,iv).','LineWidth',1.05);
qualityAxes(ax,fitLow);ylim(ax,[0 1.03]);ylabel(ax,'Multiple coherence');
title(ax,'Validation Welch: four output coherences','FontWeight','normal');
legend(ax,{'y_1','y_2','y_3','y_4'},'Location','best','Box','off');
title(tl,'Low-frequency data quality and split-record consistency','FontWeight','normal');
subtitle(tl,{sprintf('Diagnostics extend to 0.5 Hz; the joint model starts at %g Hz (vertical line).',fitLow), ...
    'Two halves of one closed-loop record: repeatability is not independent testing, and coherence does not establish an unbiased plant estimate.'},'FontSize',10);
files=[files,exportFigure(fig,out,'spring_low_frequency_quality')];clear cleanup;

% Comparable same-estimator errors quantify low-band improvement and whether
% fitting the spring dynamics changes the already fitted upper band.
bands=[fitLow 25;25 800;fitLow 800;100 150];rows={};aggregate=zeros(2,size(bands,1));
for m=1:2
    for b=1:size(bands,1)
        mask=fr>=bands(b,1) & fr<=bands(b,2);
        if b==1,mask=mask & fr<25;end
        aggregate(m,b)=100*norm(Hmodel{m}(mask,:)-Href(mask,:),'fro')/norm(Href(mask,:),'fro');
        for y=1:4
            for u=1:4
                ch=y+4*(u-1);relative=100*norm(Hmodel{m}(mask,ch)-Href(mask,ch))/max(norm(Href(mask,ch)),realmin);
                rows(end+1,:)={names{m},y,u,bands(b,1),bands(b,2),nnz(mask),relative}; %#ok<AGROW>
            end
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'model','output_channel','input_channel','fmin_hz','fmax_hz', ...
    'frequency_count','relative_complex_difference_to_validation_LRM_percent'});
csvPath=fullfile(out,'spring_band_comparison.csv');writetable(metrics,csvPath);files{end+1}=csvPath;
checks=struct('files',{files},'reference_half_window',45,'training_low_hz',fitLow, ...
    'reference_agreement_relative_error',referenceAgreement,'prediction_reconstruction_relative_errors',reconstruction, ...
    'dense_frequency_count',numel(fd),'evaluation_bands_hz',bands,'aggregate_relative_errors_percent',aggregate, ...
    'aggregate_row_order','old model, joint model','lowest_diagnostic_frequency_hz',min(fq), ...
    'low_frequency_nonfinite_singular_value_rows',nnz(any(~isfinite(sv),2)), ...
    'validation_scope','First/second halves of one closed-loop record, not independent experimental runs', ...
    'old_model_low_frequency_scope','Below 25 Hz is extrapolation', ...
    'units','Recorded output per motor-command unit; no unverified physical input conversion');
disp(checks);
end

function H=transfer(fit,f,delaySeconds)
f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if fit.direct,B(:,end)=1;end
P=B*fit.coefficients;assert(size(P,2)==16 && numel(delaySeconds)==16,'Expected 16 delayed FRFs.');
H=exp(-1i*2*pi*f.*delaySeconds(:).').*P;
end

function files=exportFigure(fig,out,stem)
png=fullfile(out,[stem '.png']);pdf=fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);exportgraphics(fig,pdf,'ContentType','vector');files={png,pdf};
end

function matrixAxes(ax,y,u,range)
style(ax);xlim(ax,range);
if range(2)<=50,ticks=[3 5 10 25 50];else,ticks=[3 10 25 100 300 800];end
ticks=ticks(ticks>=range(1) & ticks<=range(2));xticks(ax,ticks);
title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
if y<4,xticklabels(ax,{});end
end

function qualityAxes(ax,fitLow)
style(ax);xlim(ax,[.5 50]);xticks(ax,[.5 1 3 10 25 50]);xlabel(ax,'Frequency (Hz)');
xline(ax,fitLow,':','Color',[.55 .55 .55],'HandleVisibility','off');
end

function style(ax)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on','LineWidth',.65, ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
