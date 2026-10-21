function checks=plot_srivc_results(tag)
% Publication-style display of raw-spectrum additive SRIVC experiment.
% Gray curves are ordinary Welch references, not LRM or unbiased plant truth.
% No plot processing changes the complex data used to calculate errors.
if nargin<1,tag='z0_mid_fixed_full';end
out=fileparts(mfilename('fullpath'));base=fileparts(out);
S=load(fullfile(out,['srivc_' tag '.mat']),'result');r=S.result;
W=load(fullfile(base,'analysis_20261015', ...
    sprintf('auto_identify_motor1234_z%g_spectral.mat',r.position_mm)),'validation');
ref=W.validation;fit=r.fit;ny=fit.n_outputs;nu=fit.n_inputs;
blue=[.03 .24 .43];gray=[.64 .64 .64];orange=[.72 .31 .12];green=[.08 .32 .16];

% Dense logarithmic grid augmented around every narrow fitted resonance.
f=logspace(log10(25),log10(800),8000).';
for m=1:numel(fit.frequency_hz)
    fm=fit.frequency_hz(m);span=max(30*fm*fit.damping(m),1);
    f=[f;linspace(max(25,fm-span),min(800,fm+span),1200).']; %#ok<AGROW>
end
f=unique(f(f>=25 & f<=800));H=transferMatrix(r,f);
idx=ref.f>=25 & ref.f<=800;fr=ref.f(idx);Hr=ref.H(:,:,idx);

for kind={'magnitude','phase','phase_independent'}
    name=kind{1};phase=~strcmp(name,'magnitude');
    fig=figure('Visible','off','Color','w','Position',[50 50 1440 1020]);
    cleanup=onCleanup(@()closeOwned(fig));
    tl=tiledlayout(fig,ny,nu,'TileSpacing','compact','Padding','compact');
    for y=1:ny
        for u=1:nu
            ax=nexttile(tl);hold(ax,'on');hm=squeeze(H(y,u,:));hr=squeeze(Hr(y,u,:));
            if phase
                ym=unwrap(angle(hm))*180/pi;
                if strcmp(name,'phase')
                    referencePhase=interp1(f,ym,fr,'linear');
                    yr=angle(hr)*180/pi;
                    % Integer 360-degree branch display only; no smoothing.
                    yr=yr+360*round((referencePhase-yr)/360);
                else
                    yr=unwrap(angle(hr))*180/pi;
                end
                p1=semilogx(ax,fr,yr,'.','Color',gray,'MarkerSize',3);
            else
                ym=20*log10(max(abs(hm),realmin));
                yr=20*log10(max(abs(hr),realmin));
                p1=semilogx(ax,fr,yr,'Color',gray,'LineWidth',.65);
            end
            p2=semilogx(ax,f,ym,'Color',blue,'LineWidth',1.15);
            styleAxes(ax);xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);
            title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
            if y<ny,xticklabels(ax,{});end
            if y==1 && u==1,legendHandles=[p1 p2];end
        end
    end
    lg=legend(legendHandles,{'Experimental FRF estimate (Welch H1)','SRIVC model FRF'}, ...
        'Orientation','horizontal','Box','off','FontSize',11);
    lg.Layout.Tile='north';
    xlabel(tl,'Frequency (Hz)','FontSize',12);
    if phase
        ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
        if strcmp(name,'phase')
            subtitle(tl,'FRF phase: gray points shifted by integer multiples of 360 deg onto the model branch; no smoothing','FontSize',10);
        else
            subtitle(tl,'FRF phase diagnostic: experimental and model phase unwrapped independently; no branch alignment','FontSize',10);
        end
    else
        ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);
        subtitle(tl,sprintf('z = %g mm; FRF comparison over 25-800 Hz; experimental FRF is a Welch H1 estimate',r.position_mm),'FontSize',10);
    end
    exportgraphics(fig,fullfile(out,['srivc_' tag '_' name '.png']),'Resolution',200);
    if ~strcmp(name,'phase_independent')
        exportgraphics(fig,fullfile(out,['srivc_' tag '_' name '.pdf']),'ContentType','vector');
    end
    clear cleanup;
end

% Complex FRF residual at the exact experimental FRF frequencies.
% H1 is an experimental estimate, not a known plant or unbiased truth.
% Subtract complex responses BEFORE taking magnitude or logarithms.
HmodelRef=transferMatrix(r,fr);Hresidual=Hr-HmodelRef;
chosen=[1 1;3 2;4 4]; % Fixed (input, output) pairs, not selected by fit error.
fig=figure('Visible','off','Color','w','Position',[60 60 1400 850]);
cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for c=1:size(chosen,1)
    u=chosen(c,1);y=chosen(c,2);hr=squeeze(Hr(y,u,:));
    hm=squeeze(H(y,u,:));he=squeeze(Hresidual(y,u,:));
    ax=nexttile(tl,c);hold(ax,'on');
    p1=semilogx(ax,fr,20*log10(max(abs(hr),realmin)),'Color',gray,'LineWidth',.85);
    p2=semilogx(ax,f,20*log10(max(abs(hm),realmin)),'Color',blue,'LineWidth',1.2);
    if c==1
        p3=plot(ax,NaN,NaN,'Color',green,'LineWidth',1);
        legendHandles=[p1 p2 p3];
    end
    styleAxes(ax);xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);
    title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
    ylabel(ax,'FRF magnitude (dB)');
    ax=nexttile(tl,3+c);hold(ax,'on');
    semilogx(ax,fr,20*log10(max(abs(hr),realmin)),'Color',gray,'LineWidth',.85);
    semilogx(ax,fr,20*log10(max(abs(he),realmin)),'Color',green,'LineWidth',.85);
    styleAxes(ax);xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);
    xlabel(ax,'Frequency (Hz)');ylabel(ax,'FRF / residual magnitude (dB)');
end
lg=legend(legendHandles,{'Experimental FRF estimate (Welch H1)','SRIVC model FRF', ...
    'Complex FRF residual |H_{exp}-H_{model}|'},'Orientation','horizontal','Box','off','FontSize',10);
lg.Layout.Tile='north';
subtitle(tl,{'Top: experimental and identified FRFs. Bottom: experimental FRF and complex FRF residual.', ...
    'Residual computed at identical frequencies; dB referenced to recorded output / motor-command unit.'},'FontSize',10);
exportgraphics(fig,fullfile(out,['srivc_' tag '_frf_residual_examples.png']),'Resolution',200);
exportgraphics(fig,fullfile(out,['srivc_' tag '_frf_residual_examples.pdf']),'ContentType','vector');
clear cleanup;

% Full matrix prevents the illustrative three-channel selection hiding errors.
fig=figure('Visible','off','Color','w','Position',[50 50 1440 1020]);
cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,ny,nu,'TileSpacing','compact','Padding','compact');
for y=1:ny
    for u=1:nu
        ax=nexttile(tl);hold(ax,'on');hr=squeeze(Hr(y,u,:));he=squeeze(Hresidual(y,u,:));
        p1=semilogx(ax,fr,20*log10(max(abs(hr),realmin)),'Color',gray,'LineWidth',.7);
        p2=semilogx(ax,fr,20*log10(max(abs(he),realmin)),'Color',green,'LineWidth',.8);
        styleAxes(ax);xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);
        title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
        if y<ny,xticklabels(ax,{});end
        if y==1 && u==1,legendHandles=[p1 p2];end
    end
end
lg=legend(legendHandles,{'Experimental FRF estimate (Welch H1)', ...
    'Complex FRF residual |H_{exp}-H_{model}|'},'Orientation','horizontal','Box','off','FontSize',11);
lg.Layout.Tile='north';xlabel(tl,'Frequency (Hz)','FontSize',12);
ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);
subtitle(tl,'All 16 channels: complex difference at the experimental FRF frequencies; no dB subtraction','FontSize',10);
exportgraphics(fig,fullfile(out,['srivc_' tag '_frf_residual_matrix.png']),'Resolution',200);
clear cleanup;

% Unweighted relative complex differences against Welch H1, separately by band.
frfBands=[25 800;200 735];frfRows=zeros(ny*nu*size(frfBands,1),8);row=0;
for y=1:ny
    for u=1:nu
        hr=squeeze(Hr(y,u,:));he=squeeze(Hresidual(y,u,:));
        for b=1:size(frfBands,1)
            inBand=fr>=frfBands(b,1) & fr<=frfBands(b,2);
            valid=inBand & isfinite(hr) & isfinite(he);
            denominator=norm(hr(valid));relative=NaN;
            if denominator>0,relative=100*norm(he(valid))/denominator;end
            row=row+1;
            frfRows(row,:)=[y,u,frfBands(b,:),nnz(valid),nnz(inBand)-nnz(valid),relative,norm(he(valid))];
        end
    end
end
frfTable=array2table(frfRows,'VariableNames',{'output_channel','input_channel', ...
    'fmin_hz','fmax_hz','finite_frequency_count','excluded_nonfinite_count', ...
    'relative_complex_difference_to_Welch_H1_percent','complex_residual_l2_recorded_units'});
writetable(frfTable,fullfile(out,['srivc_' tag '_frf_reference_errors.csv']));
frfAggregate=zeros(size(frfBands,1),6);
for b=1:size(frfBands,1)
    inBand=fr>=frfBands(b,1) & fr<=frfBands(b,2);
    hb=Hr(:,:,inBand);eb=Hresidual(:,:,inBand);valid=isfinite(hb)&isfinite(eb);
    perChannel=frfTable.relative_complex_difference_to_Welch_H1_percent(frfTable.fmin_hz==frfBands(b,1));
    frfAggregate(b,:)=[frfBands(b,:),100*norm(eb(valid))/norm(hb(valid)), ...
        min(perChannel),median(perChannel),max(perChannel)];
end
frfSummary=array2table(frfAggregate,'VariableNames',{'fmin_hz','fmax_hz', ...
    'matrix_complex_difference_percent','minimum_channel_percent','median_channel_percent','maximum_channel_percent'});
writetable(frfSummary,fullfile(out,['srivc_' tag '_frf_reference_summary.csv']));

% Raw output spectra. Error is calculated before binning, using complex data.
v=r.validation;fv=v.f(:);Y=v.Y;P=v.prediction;R=Y-P;
bands=[25 200;200 385;385 550;650 800];errors=zeros(4,1);
fig=figure('Visible','off','Color','w','Position',[60 60 1340 830]);
cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for b=1:4
    ax=nexttile(tl);hold(ax,'on');ii=fv>=bands(b,1) & fv<=bands(b,2);
    ww=(fv(ii)/100).^2;
    errors(b)=100*norm(R(ii,:).*ww,'fro')/norm(Y(ii,:).*ww,'fro');
    % Bin squared magnitudes of all four outputs, not their complex average.
    [fb,amp]=binPowers(fv(ii),[sum(abs(Y(ii,:)).^2,2), ...
        sum(abs(P(ii,:)).^2,2),sum(abs(R(ii,:)).^2,2)],.5);
    amp=amp*r.output_scale;
    p1=semilogx(ax,fb,20*log10(max(amp(:,1),realmin)),'Color',gray,'LineWidth',1.2);
    p2=semilogx(ax,fb,20*log10(max(amp(:,2),realmin)),'Color',blue,'LineWidth',1.1);
    p3=semilogx(ax,fb,20*log10(max(amp(:,3),realmin)),'Color',orange,'LineWidth',.8);
    styleAxes(ax);xlim(ax,bands(b,:));
    title(ax,sprintf('%g-%g Hz: weighted complex error %.2f%%',bands(b,1),bands(b,2),errors(b)),'FontWeight','normal');
    xlabel(ax,'Frequency (Hz)');ylabel(ax,'Output spectral amplitude (dB, recorded units)');
    if b==1,legendHandles=[p1 p2 p3];end
end
lg=legend(legendHandles,{'Recorded output','SRIVC prediction + fitted boundary terms','Complex residual'}, ...
    'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
subtitle(tl,{'Second 30 s: transfer coefficients frozen; boundary/transient terms refitted', ...
    'Display: 0.5 Hz power bins across four outputs. Error: original complex bins weighted by f^2.'},'FontSize',10);
exportgraphics(fig,fullfile(out,['srivc_' tag '_raw_output_validation.png']),'Resolution',180);
clear cleanup;

% Iteration trace shows the actual IV score, not only a visual fit criterion.
h=fit.history;fig=figure('Visible','off','Color','w','Position',[80 80 1190 480]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);semilogy(ax,h(:,1),100*h(:,2),'o-','Color',blue,'LineWidth',1.3,'MarkerSize',3);
grid(ax,'on');xlabel(ax,'SRIVC iteration');ylabel(ax,'Weighted training error (%)');
title(ax,'Raw input-output spectral fit','FontWeight','normal');
ax=nexttile(tl);semilogy(ax,h(:,1),h(:,3),'o-','Color',blue,'LineWidth',1.3,'MarkerSize',3);hold(ax,'on');
yline(ax,1e-6,'--','Score threshold 10^{-6}','Color',gray,'LabelHorizontalAlignment','left');
grid(ax,'on');xlabel(ax,'SRIVC iteration');ylabel(ax,'Normalized instrument-correlation score');
title(ax,sprintf('Final score %.3g',fit.score),'FontWeight','normal');
title(tl,sprintf('Additive SRIVC: %s; %d iterations',strrep(fit.reason,'_',' '),fit.iterations),'FontWeight','normal');
subtitle(tl,'Stabilized IV updates with measured regressors and model-generated instruments; no LRM','FontSize',10);
exportgraphics(fig,fullfile(out,['srivc_' tag '_convergence.png']),'Resolution',180);
clear cleanup;

% Independently reconstruct forced spectra from the plotted transfer matrix.
Hv=transferMatrix(r,fv);reconstructed=complex(zeros(size(Y)));
for k=1:nu
    reconstructed=reconstructed+squeeze(Hv(:,k,:)).'.*(v.U(:,k)*r.input_scale(k));
end
forcedPhysical=v.forced_prediction*r.output_scale;
checks=struct('tag',tag,'dense_grid_count',numel(f), ...
    'all_transfer_values_finite',all(isfinite(H(:))), ...
    'transfer_reconstruction_relative_error',norm(reconstructed-forcedPhysical,'fro')/max(norm(forcedPhysical,'fro'),realmin), ...
    'display_validation_bands_hz',bands,'validation_weighted_error_percent',errors, ...
    'frf_comparison','Complex Hexp-Hmodel at identical Welch H1 frequencies; not dB subtraction or true-plant error', ...
    'frf_reference_error_table',frfTable,'fixed_example_channels_input_output',chosen, ...
    'frf_reference_summary',frfSummary, ...
    'gray_reference','Ordinary 2-second Welch; second half; display only, not unbiased truth', ...
    'phase_display','Measured principal phases shifted by integer 360 deg to model branch; independent-unwrap diagnostic also exported', ...
    'validation','Transfer frozen, record-specific boundary terms fitted on second half; no smoothing in error calculation', ...
    'fixed_background_modes',fit.fixed_modes);
assert(checks.all_transfer_values_finite,'Nonfinite plotted transfer function.');
assert(checks.transfer_reconstruction_relative_error<1e-10,'Transfer plotting coefficient layout or physical scaling mismatch.');
save(fullfile(out,['srivc_' tag '_plot_checks.mat']),'checks');
disp(checks);
end

function H=transferMatrix(r,f)
fit=r.fit;f=f(:);nm=numel(fit.reference_hz);nu=fit.n_inputs;ny=fit.n_outputs;
q=1i*f./fit.reference_hz(:).';A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';
H=complex(zeros(ny,nu,numel(f)));nb=2*(nu+1);
for u=1:nu
    values=complex(zeros(numel(f),ny));
    for m=1:nm
        ix=(m-1)*nb+2*u-1;
        values=values+(ones(numel(f),1)*fit.coefficients(ix,:)+q(:,m)*fit.coefficients(ix+1,:))./A(:,m);
    end
    if fit.direct,values=values+fit.coefficients(nm*nb+u,:);end
    values=values*(r.output_scale/r.input_scale(u));
    H(:,u,:)=reshape(values.',ny,1,[]);
end
end

function [fbin,amplitude]=binPowers(f,powerValues,width)
group=floor((f-min(f))/width)+1;groups=unique(group);
fbin=zeros(numel(groups),1);amplitude=zeros(numel(groups),size(powerValues,2));
for k=1:numel(groups)
    ii=group==groups(k);fbin(k)=mean(f(ii));amplitude(k,:)=sqrt(mean(powerValues(ii,:),1));
end
end

function styleAxes(ax)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10, ...
    'LineWidth',.65,'Box','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
