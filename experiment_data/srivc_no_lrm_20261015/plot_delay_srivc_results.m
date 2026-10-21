function checks=plot_delay_srivc_results
% Plot the complete measured-boundary FRF, including all channel delays.
% Main model: ell_jk=t0+dy_j+du_k. A common delay is only a comparison.
% The fitted rational P(s) alone is never presented as the measured FRF.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'delay_identification_20261016');
S=load(fullfile(out,'common_delay_srivc.mat'),'result');common=S.result;
S=load(fullfile(out,'channel_delay_srivc.mat'),'result');channel=S.result;
fr=common.f(:);assert(isequal(fr,channel.f(:)),'Saved comparison grids differ.');
Hexp=common.welchValidation;assert(size(Hexp,2)==16,'Expected 16 FRF channels.');
ell=channel.delay_matrix_seconds(:).';ellCommon=common.delay_matrix_seconds(:).';
assert(numel(ell)==16 && numel(ellCommon)==16,'Delay matrix must contain all 16 responses.');
blue=[.03 .25 .45];gray=[.65 .65 .65];green=[.07 .32 .16];
orange=[.79 .34 .08];purple=[.55 .16 .56];

% Dense model evaluation retains narrow resonances rather than interpolating
% sampled FRF values. Add resonance neighborhoods for each comparison model.
f=logspace(log10(25),log10(800),8000).';
fits={channel.selected.fit,common.selected.fit,common.baseline};
for j=1:numel(fits)
    fit=fits{j};
    for m=1:numel(fit.frequency_hz)
        center=fit.frequency_hz(m);span=max(30*center*fit.damping(m),1);
        low=max(25,center-span);high=min(800,center+span);
        if high>low,f=[f;linspace(low,high,1200).'];end %#ok<AGROW>
    end
end
f=unique([25;180;260;800;f]);H=delayedTransfer(channel.selected.fit,f,ell);
Href=delayedTransfer(channel.selected.fit,fr,ell);
Hcommon=delayedTransfer(common.selected.fit,f,ellCommon);
HcommonRef=delayedTransfer(common.selected.fit,fr,ellCommon);
Hbaseline=delayedTransfer(common.baseline,f,zeros(1,16));
HbaselineRef=delayedTransfer(common.baseline,fr,zeros(1,16));
channelError=norm(Href-channel.selected.H,'fro')/max(norm(channel.selected.H,'fro'),realmin);
commonError=norm(HcommonRef-common.selected.H,'fro')/max(norm(common.selected.H,'fro'),realmin);
baselineError=norm(HbaselineRef-common.baselineH,'fro')/max(norm(common.baselineH,'fro'),realmin);
assert(all(isfinite(H(:))) && all(isfinite(Href(:))),'Nonfinite delay-augmented FRF.');
assert(channelError<1e-10 && commonError<1e-10 && baselineError<1e-10, ...
    'Plotted FRF disagrees with saved model: check coefficient layout or delay phase.');
assert(isequal(size(Href),size(Hexp)),'Model and experimental FRFs must use exactly the same grid.');
Hresidual=Hexp-Href; % Complex subtraction BEFORE taking abs or dB.

for kind={'magnitude','phase','phase_independent'}
    name=kind{1};isPhase=~strcmp(name,'magnitude');
    fig=figure('Visible','off','Color','w','Position',[45 45 1450 1040]);
    cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');
            if isPhase
                model=unwrap(angle(H(:,ch)))*180/pi;
                if strcmp(name,'phase')
                    experimental=angle(Hexp(:,ch))*180/pi;
                    branch=interp1(f,model,fr,'linear');
                    experimental=experimental+360*round((branch-experimental)/360);
                else
                    experimental=unwrap(angle(Hexp(:,ch)))*180/pi;
                end
                h1=semilogx(ax,fr,experimental,'.','Color',gray,'MarkerSize',3);
            else
                model=toDb(H(:,ch));experimental=toDb(Hexp(:,ch));
                h1=semilogx(ax,fr,experimental,'Color',gray,'LineWidth',.7);
            end
            h2=semilogx(ax,f,model,'Color',blue,'LineWidth',1.05);
            matrixAxes(ax,y,u);
            if y==1 && u==1,handles=[h1 h2];end
        end
    end
    lg=legend(handles,{'Experimental FRF','Delay-augmented SRIVC'}, ...
        'Orientation','horizontal','Box','off','FontSize',11);lg.Layout.Tile='north';
    xlabel(tl,'Frequency (Hz)','FontSize',12);
    if isPhase
        ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
        if strcmp(name,'phase')
            subtitle(tl,{'Per-input/per-output additive delays are included in the plotted model.', ...
                'Gray phase points select the model branch by integer 360-deg shifts; no smoothing. Independent unwrap is saved separately.'},'FontSize',10);
        else
            subtitle(tl,{'Diagnostic: experimental and model phase unwrapped independently; no phase-branch alignment.', ...
                'The complete FRF includes the individual input/output additive delays.'},'FontSize',10);
        end
    else
        ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);
        subtitle(tl,{'Per-input/per-output additive effective delays; common delay is only a comparison.', ...
            'Experimental FRF: second-half Welch H1. Effective delays are not certified hardware delays.'},'FontSize',10);
    end
    exportgraphics(fig,fullfile(out,['delay_srivc_' name '.png']),'Resolution',200);
    if ~strcmp(name,'phase_independent')
        exportgraphics(fig,fullfile(out,['delay_srivc_' name '.pdf']),'ContentType','vector');
    end
    clear cleanup;
end

% All 16 residuals use precisely the original experimental frequency bins.
fig=figure('Visible','off','Color','w','Position',[45 45 1450 1040]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
for y=1:4
    for u=1:4
        ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');
        h1=semilogx(ax,fr,toDb(Hexp(:,ch)),'Color',gray,'LineWidth',.7);
        h2=semilogx(ax,fr,toDb(Hresidual(:,ch)),'Color',green,'LineWidth',.8);
        matrixAxes(ax,y,u);
        if y==1 && u==1,handles=[h1 h2];end
    end
end
lg=legend(handles,{'Experimental FRF','Complex FRF residual |H_{exp}-H_{model}|'}, ...
    'Orientation','horizontal','Box','off','FontSize',11);lg.Layout.Tile='north';
xlabel(tl,'Frequency (Hz)','FontSize',12);ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);
subtitle(tl,{'Full per-input/per-output delay model evaluated at identical experimental frequencies.', ...
    'Complex FRFs are subtracted before magnitude conversion; Welch H1 is an estimate, not a known plant.'},'FontSize',10);
exportgraphics(fig,fullfile(out,'delay_srivc_residual_matrix.png'),'Resolution',200);
exportgraphics(fig,fullfile(out,'delay_srivc_residual_matrix.pdf'),'ContentType','vector');clear cleanup;

% Representative antiresonances, selected by fixed channel indices.
local=f>=180 & f<=260;localRef=fr>=180 & fr<=260;fl=f(local);flr=fr(localRef);
cols=[1 6];fig=figure('Visible','off','Color','w','Position',[60 60 1340 870]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    for col=1:2
        ch=cols(col);ax=nexttile(tl);hold(ax,'on');
        hr=Hexp(localRef,ch);hm=[Hbaseline(local,ch),Hcommon(local,ch),H(local,ch)];
        if row==1
            yr=toDb(hr);ym=toDb(hm);
        else
            yr=unwrap(angle(hr))*180/pi;ym=unwrap(angle(hm),[],1)*180/pi;
            for m=1:3,ym(:,m)=ym(:,m)+360*round((yr(1)-ym(1,m))/360);end
        end
        h1=semilogx(ax,flr,yr,'Color',gray,'LineWidth',1.1);
        h2=semilogx(ax,fl,ym(:,1),'Color',orange,'LineWidth',1.15);
        h3=semilogx(ax,fl,ym(:,2),'--','Color',purple,'LineWidth',1.15);
        h4=semilogx(ax,fl,ym(:,3),'Color',blue,'LineWidth',1.3);
        style(ax);xlim(ax,[180 260]);xticks(ax,[180 200 220 240 260]);
        title(ax,sprintf('$u_%d\\rightarrow y_%d$',col,col),'Interpreter','latex','FontWeight','normal');
        xlabel(ax,'Frequency (Hz)');
        if row==1,ylabel(ax,'FRF magnitude (dB)');else,ylabel(ax,'Unwrapped phase (deg)');end
        if row==1 && col==1,handles=[h1 h2 h3 h4];end
    end
end
lg=legend(handles,{'Experimental FRF','Without delay','Common delay','Input/output additive delays'}, ...
    'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
title(tl,'Antiresonance comparison: channel-dependent delays','FontWeight','normal');
subtitle(tl,{'All curves include their fitted delay phases. Model selection used training data only.', ...
    'Local phase: independent unwrap, then one constant 360-deg branch shift per whole curve; no smoothing.'},'FontSize',10);
exportgraphics(fig,fullfile(out,'delay_srivc_antiresonance.png'),'Resolution',200);
exportgraphics(fig,fullfile(out,'delay_srivc_antiresonance.pdf'),'ContentType','vector');clear cleanup;

checks=struct('main_model','per-input/per-output additive delays', ...
    'dense_grid_count',numel(f),'dense_values_finite',all(isfinite(H(:))), ...
    'channel_prediction_reconstruction_relative_error',channelError, ...
    'common_prediction_reconstruction_relative_error',commonError, ...
    'baseline_prediction_reconstruction_relative_error',baselineError, ...
    'delay_matrix_microseconds',channel.delay_matrix_seconds*1e6, ...
    'welch_reference_complex_relative_error_percent',100*norm(Hresidual,'fro')/norm(Hexp,'fro'), ...
    'main_phase_display','Model unwrap; gray samples select display branch by integer 360-degree shifts; independent diagnostic saved', ...
    'notch_phase_display','Independent local unwrap with one constant integer 360-degree offset per model curve', ...
    'residual_definition','Hexperimental-Hmodel at exactly the same frequencies, before abs or log10', ...
    'model_definition','H_jk(s)=exp(-s*ell_jk)*P_jk(s); ell_jk=t0+dy_j+du_k', ...
    'physical_scaling','none; coefficients already in recorded output/motor-command units');
disp(checks);
end

function H=delayedTransfer(fit,f,delaySeconds)
f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';
B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if fit.direct,B(:,end)=1;end
P=B*fit.coefficients;
assert(size(P,2)==numel(delaySeconds),'Delay/channel mapping mismatch.');
H=exp(-1i*2*pi*f.*delaySeconds(:).').*P;
end

function matrixAxes(ax,y,u)
style(ax);xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);
title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
if y<4,xticklabels(ax,{});end
end

function style(ax)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on', ...
    'LineWidth',.65,'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
