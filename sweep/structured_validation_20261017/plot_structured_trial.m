function checks=plot_structured_trial(tag)
% Publication-style comparison of saved estimates; no fitting or data changes.
% Phase curves are independently unwrapped. Each whole curve receives one
% constant 360-degree shift placing its first plotted phase in [-270,90).
% This branch rule never uses a model to change measured points.
if nargin<1,tag='z0';end
tag=char(string(tag));out=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(out));source=fullfile(out,['trial_' tag '.mat']);
loaded=load(source,'result');r=loaded.result;
assert(all(isfield(r,{'f','H_training','H_validation','baseline_H','fixed','delay'})), ...
    'Saved result is missing required fields.');
f=r.f(:);nf=numel(f);assert(all(isfinite(f))&&all(diff(f)>0));
arrays={r.H_training,r.H_validation,r.baseline_H,r.fixed.H,r.delay.H};
for k=1:numel(arrays)
    assert(isequal(size(arrays{k}),[nf,16])&&all(isfinite(arrays{k}(:))), ...
        'Each saved FRF must be finite, frequency-by-16.');
end
Ty=readmatrix(fullfile(root,'sweep','srivc_retest_20261016','geometric_Ty_full.csv'));
Tu=readmatrix(fullfile(root,'sweep','srivc_retest_20261016','geometric_Tu_full.csv'));
assert(isequal(size(Ty),[4,4])&&isequal(size(Tu),[4,4]) ...
    &&all(isfinite([Ty(:);Tu(:)])));
keep=f>=25 & f<=800;assert(nnz(keep)>20);fp=f(keep);
physical=struct('data',r.H_validation(keep,:),'baseline',r.baseline_H(keep,:), ...
    'fixed',r.fixed.H(keep,:),'delay',r.delay.H(keep,:));
logical=struct;
names=fieldnames(physical);
for k=1:numel(names),logical.(names{k})=transform(physical.(names{k}),Ty,Tu);end
% The logical axes are fixed external geometric coordinates, not estimated
% nullspaces selected separately for the competing models.
labelsPhysicalIn={'u_1','u_2','u_3','u_4'};
labelsPhysicalOut={'y_1','y_2','y_3','y_4'};
labelsLogical={'R_x','R_y','Z','N'};
files={};
for phase=[false,true]
    kind='magnitude';if phase,kind='phase';end
    files=[files,drawGrid(fp,physical,1:4,labelsPhysicalIn,labelsPhysicalOut, ...
        phase,'Physical motor-to-sensor channels',out,[tag '_physical_' kind],r)]; %#ok<AGROW>
    files=[files,drawGrid(fp,logical,1:4,labelsLogical,labelsLogical, ...
        phase,'Fixed geometric coordinates, including nullspace',out,[tag '_full_logical_' kind],r)]; %#ok<AGROW>
    files=[files,drawGrid(fp,logical,1:3,labelsLogical,labelsLogical, ...
        phase,'Fixed geometric logical axes',out,[tag '_logical_' kind],r)]; %#ok<AGROW>
end
files=[files,drawNullFocus(fp,logical,out,[tag '_nullspace_resonance_residual'],r)];
checks=struct('source',source,'files',{files},'figure_count',numel(files)/2, ...
    'static_Ty',Ty,'static_Tu',Tu,'plotted_band_hz',[fp(1),fp(end)], ...
    'phase_display','Independent unwrap; one constant 360-degree shift per whole curve puts its first plotted sample in [-270,90). No pointwise model-based branch selection.', ...
    'frequency_display','Saved FRF points only; no interpolation, smoothing, or synthetic dense evaluation.', ...
    'validation_scope','Previously inspected same-record held-out development FRF, not independent repeat data or a certified unbiased estimate.', ...
    'algorithm_scope','Baseline additive SRIVC versus fixed-pole structural refinement; delay stage updates additive channel delays with SRIVC poles fixed.', ...
    'logical_scope','One unchanged geometric input/output coordinate system applied to data and every model; no newly identified calibration applied.', ...
    'saved_at',char(datetime('now')));
save(fullfile(out,['plot_checks_' tag '.mat']),'checks');
fprintf('Saved %d comparison figures for %s.\n',checks.figure_count,tag);
end

function L=transform(H,Ty,Tu)
L=complex(zeros(size(H)));
for k=1:size(H,1)
    g=Ty*reshape(H(k,:),4,4)*Tu;L(k,:)=g(:).';
end
end

function files=drawGrid(f,h,axesToShow,inputLabels,outputLabels,isPhase,caption,out,stem,r)
n=numel(axesToShow);fig=figure('Visible','off','Color','white', ...
    'Units','pixels','Position',[100,100,1440,1080]);
closer=onCleanup(@()closeFigure(fig)); %#ok<NASGU>
tl=tiledlayout(fig,n,n,'TileSpacing','compact','Padding','compact');
palette=colors;firstHandles=[];
for rr=1:n
    for cc=1:n
        output=axesToShow(rr);input=axesToShow(cc);ch=output+4*(input-1);
        ax=nexttile(tl);hold(ax,'on');
        p1=semilogx(ax,f,displayValues(h.data(:,ch),isPhase),'-', ...
            'Color',palette.data,'LineWidth',.65);
        p2=semilogx(ax,f,displayValues(h.baseline(:,ch),isPhase),'--', ...
            'Color',palette.baseline,'LineWidth',1.0);
        p3=semilogx(ax,f,displayValues(h.delay(:,ch),isPhase),'-', ...
            'Color',palette.delay,'LineWidth',1.15);
        if isempty(firstHandles),firstHandles=[p1,p2,p3];end
        styleAxes(ax,f,isPhase);
        title(ax,[inputLabels{input} ' \rightarrow ' outputLabels{output}], ...
            'FontWeight','normal','Interpreter','tex','FontSize',11);
        if rr<n,ax.XTickLabel=[];end
    end
end
lg=legend(firstHandles,{'Same-record held-out development FRF','Additive SRIVC baseline', ...
    'Fixed-pole structure + channel-delay refinement'}, ...
    'Location','northoutside','Orientation','horizontal','Box','off','FontSize',10);
lg.Layout.Tile='north';
xlabel(tl,'Frequency (Hz)','FontSize',12);
if isPhase,ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
else,ylabel(tl,'Magnitude (dB, recorded-unit ratio)','FontSize',12);end
title(tl,sprintf('%s | z = %g mm',caption,r.position_mm), ...
    'FontWeight','normal','FontSize',15,'Interpreter','none');
if isPhase
    subtitle(tl,'Independent phase unwrap; one constant branch shift per curve. Same-record development check.', ...
        'FontSize',10,'Interpreter','none');
else
    subtitle(tl,'Measured FRF is an estimate, not ground truth. Structural refinement keeps the SRIVC poles fixed.', ...
        'FontSize',10,'Interpreter','none');
end
files=exportPair(fig,out,stem);
end

function files=drawNullFocus(f,h,out,stem,r)
ch=16;p=colors;fig=figure('Visible','off','Color','white', ...
    'Units','pixels','Position',[100,100,1440,920]);
closer=onCleanup(@()closeFigure(fig)); %#ok<NASGU>
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
modelNames={'baseline','fixed','delay'};
modelLabels={'Additive SRIVC','Structure, fixed delays','Structure + channel delays'};
palette=[p.baseline;p.fixed;p.delay];styles={'--','-.','-'};
first=[];
for panel=1:3
    ax=nexttile(tl);hold(ax,'on');
    q=semilogx(ax,f,displayValues(h.data(:,ch),panel==2),'-', ...
        'Color',p.data,'LineWidth',.85);handles=q;
    for k=1:3
        q=semilogx(ax,f,displayValues(h.(modelNames{k})(:,ch),panel==2),styles{k}, ...
            'Color',palette(k,:),'LineWidth',1.15);handles(end+1)=q; %#ok<AGROW>
    end
    styleAxes(ax,f,panel==2);
    if panel==1,title(ax,'N \rightarrow N: full band');ylabel(ax,'Magnitude (dB)');first=handles;
    elseif panel==2,title(ax,'N \rightarrow N: independently unwrapped phase');ylabel(ax,'Phase (deg)');
    else,title(ax,'N \rightarrow N: first plate-mode region');xlim(ax,[250,315]);xticks(ax,[250,275,300,315]);ylabel(ax,'Magnitude (dB)');end
    xlabel(ax,'Frequency (Hz)');
end
for panel=1:3
    ax=nexttile(tl);hold(ax,'on');
    for k=1:3
        e=h.(modelNames{k})(:,ch)-h.data(:,ch);
        if panel==1
            v=displayValues(e,false);
        elseif panel==2
            % Complex relative error, not a magnitude-only difference.
            v=100*abs(e)./max(abs(h.data(:,ch)),realmin);
        else
            v=displayValues(e,false);
        end
        semilogx(ax,f,v,styles{k},'Color',palette(k,:),'LineWidth',1.05);
    end
    styleAxes(ax,f,false);xlabel(ax,'Frequency (Hz)');
    if panel==1
        title(ax,'Absolute complex residual |G_{model}-G_{FRF}|');ylabel(ax,'Residual magnitude (dB)');
    elseif panel==2
        title(ax,'Complex relative residual (can diverge near zeros)');
        ylabel(ax,'100 |G_{model}-G_{FRF}| / |G_{FRF}| (%)');set(ax,'YScale','log');
    else
        title(ax,'First plate-mode region: complex residual');ylabel(ax,'Residual magnitude (dB)');
        xlim(ax,[250,315]);xticks(ax,[250,275,300,315]);
    end
end
lg=legend(first,[{'Same-record held-out development FRF'},modelLabels], ...
    'Location','northoutside','Orientation','horizontal','Box','off','FontSize',10);
lg.Layout.Tile='north';
title(tl,sprintf('Fixed geometric nullspace | z = %g mm',r.position_mm), ...
    'FontWeight','normal','FontSize',15);
subtitle(tl,{'No newly estimated coordinate maps are applied to these curves.', ...
    'Fixed-pole structural refinement; residuals are relative to an estimated same-record development FRF.'}, ...
    'FontSize',10,'Interpreter','none');
files=exportPair(fig,out,stem);
end

function v=displayValues(z,isPhase)
if isPhase
    v=unwrap(angle(z))*180/pi;
    % Same own-first-point branch convention for data and each model.
    v=v-360*floor((v(1)+270)/360);
else
    v=20*log10(max(abs(z),realmin));
end
end

function styleAxes(ax,f,isPhase)
set(ax,'XScale','log','FontName','Arial','FontSize',9,'LineWidth',.65, ...
    'Box','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',.16,'MinorGridAlpha',.08);
grid(ax,'on');xlim(ax,[f(1),f(end)]);xticks(ax,[25,50,100,200,400,800]);
if isPhase
    ax.YAxis.Exponent=0;
end
end

function c=colors
c=struct('data',[.64,.64,.64],'baseline',[.60,.47,.73], ...
    'fixed',[.10,.48,.43],'delay',[.025,.24,.47]);
end

function files=exportPair(fig,out,stem)
files={fullfile(out,[stem '.png']),fullfile(out,[stem '.pdf'])};
exportgraphics(fig,files{1},'Resolution',200,'BackgroundColor','white');
exportgraphics(fig,files{2},'ContentType','vector','BackgroundColor','white');
end

function closeFigure(fig)
if isgraphics(fig),close(fig);end
end
