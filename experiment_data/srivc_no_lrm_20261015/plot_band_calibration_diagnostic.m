function files=plot_band_calibration_diagnostic
% Plot saved held-out FRFs only; no fitting, retuning, or model selection.
% D0/D1/D2 are auxiliary direct-target diagnostics, not the modal estimator.
base=fileparts(mfilename('fullpath'));
out=fullfile(base,'band_calibration_20261016');
source=fullfile(out,'band_calibration_campaign.mat');
assert(isfile(source),'Run run_band_calibration_diagnostic first.');
S=load(source,'campaign');campaign=S.campaign;T=campaign.metrics;
positions=cellfun(@(v)v.position_mm,campaign.records);
assert(numel(positions)==4 && all(abs(positions-[-.5 0 .5 1])<1e-12), ...
    'Unexpected position ordering.');
assert(isequal(campaign.config.training_band,[25 80]),'Unexpected primary training band.');
z0=find(abs(positions)<1e-12);assert(isscalar(z0));
f=campaign.records{z0}.f(:);methods=campaign.results{z0}.methods;
assert(all(diff(f)>0) && all(isfinite(f)),'Expected finite increasing frequencies.');
colors=struct('B0',[.57 .57 .57],'B1',[.48 .20 .59], ...
    'B4',[.48 .20 .59],'D0',[.04 .28 .55],'D1',[.04 .28 .55],'D2',[.12 .48 .26]);
codes={'B0','B4','D0','D2'};
labels={'B0: nominal static','B4: prior modal pilot + delay', ...
    'D0: band-optimal static','D2: band-profiled static + delay'};
blocks={'logical_to_N','N_to_logical'};
blockTitles={'Logical inputs to nullspace output','Nullspace input to logical outputs'};
files={};

% 1. Primary 25-80 Hz comparison. Ratios use each fixed same-position B0.
fig=newFigure([1440 560]);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
barCodes={'B1','D0','D2'};
for b=1:2
    ax=nexttile(tl);values=zeros(4,3);
    for p=1:4
        for m=1:3
            values(p,m)=metric(T,positions(p),'own_position',barCodes{m},blocks{b},'over_B0');
        end
    end
    h=bar(ax,positions,values,'grouped');hold(ax,'on');
    for m=1:3,h(m).FaceColor=colors.(barCodes{m});h(m).EdgeColor='none';end
    yline(ax,1,':','Color',colors.B0,'LineWidth',1.2,'HandleVisibility','off');
    style(ax);xticks(ax,positions);xlabel(ax,'Position z (mm)');
    ylabel(ax,'25-80 Hz held-out response RMS / B0');
    title(ax,blockTitles{b},'FontWeight','normal');
    ylim(ax,[0 max(1.05,1.12*max(values,[],'all'))]);
    if b==1,legendHandles=h;end
end
lg=legend(legendHandles,{'B1: prior modal-pilot static','D0: band-optimal static', ...
    'D2: band-profiled static + delay'},'Orientation','horizontal','Box','off');
lg.Layout.Tile='north';
title(tl,'25-80 Hz held-out LRM: both cross-response blocks','FontWeight','normal');
subtitle(tl,{'D0/D2 training band: 25-80 Hz. Auxiliary diagnostic, not the structured modal method.', ...
    'Recorded command coordinates; B0 = 1 separately for each block and position. Lower does not establish zero rigid leakage.'},'FontSize',10);
files=[files,exportPair(fig,out,'band_diagnostic_crossblock_ratios')];clear cleanup;

% 2. Freeze z=0 directions/delays and compare with separately fitted positions.
fig=newFigure([1440 560]);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
transferCodes={'D0','D0','D2','D2'};scopes={'own_position','z0_transfer','own_position','z0_transfer'};
lineStyles={'-','--','-','--'};markers={'o','s','o','s'};
for b=1:2
    ax=nexttile(tl);hold(ax,'on');h=gobjects(1,4);v=zeros(4,4);
    for m=1:4
        for p=1:4
            scope=scopes{m};if p==z0,scope='own_position';end
            v(p,m)=metric(T,positions(p),scope,transferCodes{m},blocks{b},'over_B0');
        end
        h(m)=plot(ax,positions,v(:,m),'LineStyle',lineStyles{m}, ...
            'Marker',markers{m},'Color',colors.(transferCodes{m}), ...
            'LineWidth',1.3,'MarkerSize',6,'MarkerFaceColor','w');
    end
    yline(ax,1,':','Color',colors.B0,'HandleVisibility','off');style(ax);
    xticks(ax,positions);xlabel(ax,'Evaluation position z (mm)');
    ylabel(ax,'25-80 Hz held-out response RMS / same-position B0');
    title(ax,blockTitles{b},'FontWeight','normal');
    ylim(ax,[0 max(1.05,1.12*max(v,[],'all'))]);
    if b==1,legendHandles=h;end
end
lg=legend(legendHandles,{'D0: fit each position','D0: reuse z=0', ...
    'D2: fit each position','D2: reuse z=0'},'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
title(tl,'Transfer of fixed z=0 calibration to other positions','FontWeight','normal');
subtitle(tl,{'Held-out LRM; all fits use the first-half 25-80 Hz band. At z=0, own and transfer values coincide.', ...
    'Recorded command coordinates; auxiliary direct-target diagnostic. Cross-position checks are not independent repeat experiments.'},'FontSize',10);
files=[files,exportPair(fig,out,'band_diagnostic_z0_transfer')];clear cleanup;

% 3. A vector block has a norm, not a meaningful aggregate phase.
fig=newFigure([1180 1040]);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
curveTitles=[blockTitles,{'Nullspace input to nullspace output'}];
for b=1:3
    ax=nexttile(tl);hold(ax,'on');h=gobjects(1,4);
    for m=1:4
        Z=response(methods,codes{m});
        if b==1,z=reshape(Z(4,1:3,:),3,[]).';
        elseif b==2,z=reshape(Z(1:3,4,:),3,[]).';
        else,z=reshape(Z(4,4,:),[],1);end
        h(m)=semilogx(ax,f,db(sqrt(sum(abs(z).^2,2))), ...
            'Color',colors.(codes{m}),'LineWidth',lineWidth(codes{m}));
    end
    style(ax);set(ax,'XScale','log');markTrainingBand(ax);
    xlim(ax,[f(1) f(end)]);xticks(ax,[3 10 25 80 200 500 800]);
    ylabel(ax,'Block magnitude (dB re 1 recorded-unit ratio)');
    xlabel(ax,'Frequency (Hz)');title(ax,curveTitles{b},'FontWeight','normal');
    if b==1,legendHandles=h;end
end
lg=legend(legendHandles,labels,'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
title(tl,'z=0: transformed held-out LRM responses','FontWeight','normal');
subtitle(tl,{'Cross blocks show the Euclidean norm across three logical reference channels; N to N shows its magnitude.', ...
    'Dotted vertical lines: 25-80 Hz training band. Command units are uncalibrated. Diagnostic, not a modal fit.'},'FontSize',10);
files=[files,exportPair(fig,out,'band_diagnostic_z0_crossblock_spectra')];clear cleanup;

% 4. Same static vectors in D1/D2 isolate the additional delay alignment.
fig=newFigure([1440 560]);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
increments=zeros(4,2);
for b=1:2
    ax=nexttile(tl);hold(ax,'on');
    for p=1:4
        a=metric(T,positions(p),'own_position','D1',blocks{b},'rms');
        d=metric(T,positions(p),'own_position','D2',blocks{b},'rms');
        increments(p,b)=20*log10(max(a,realmin)/max(d,realmin));
    end
    bar(ax,positions,increments(:,b),.50,'FaceColor',colors.D2,'EdgeColor','none');
    yline(ax,0,'--','Color',[.35 .35 .35]);style(ax);xticks(ax,positions);
    xlabel(ax,'Position z (mm)');ylabel(ax,'20 log_{10}(D1 RMS / D2 RMS) (dB)');
    title(ax,blockTitles{b},'FontWeight','normal');
    lo=min([0;increments(:,b)]);hi=max([0;increments(:,b)]);
    pad=max(.05,.12*max(hi-lo,.1));ylim(ax,[lo-pad hi+pad]);
end
title(tl,'Increment from delay alignment alone, 25-80 Hz','FontWeight','normal');
subtitle(tl,{'D1 and D2 use exactly the same profiled static vectors; only D2 applies the causal extra delays.', ...
    'Positive values mean lower held-out LRM cross-response. Training band: 25-80 Hz; diagnostic, not a hardware-delay measurement.'},'FontSize',10);
files=[files,exportPair(fig,out,'band_diagnostic_delay_increment')];clear cleanup;

% 5. Display the observed target-band response without calling it retention.
fig=newFigure([1440 560]);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for b=1:2
    ax=nexttile(tl);hold(ax,'on');h=gobjects(1,4);
    mask=true(size(f));if b==2,mask=f>=260-1e-6 & f<=310+1e-6;end
    for m=1:4
        Z=response(methods,codes{m});z=reshape(Z(4,4,:),[],1);
        h(m)=plot(ax,f(mask),db(z(mask)),'Color',colors.(codes{m}),'LineWidth',lineWidth(codes{m}));
    end
    style(ax);xlabel(ax,'Frequency (Hz)');ylabel(ax,'N to N magnitude (dB re 1 recorded-unit ratio)');
    if b==1,set(ax,'XScale','log');xlim(ax,[f(1) f(end)]);markTrainingBand(ax);title(ax,'Full measured frequency range','FontWeight','normal');
    else,xlim(ax,[260 310]);title(ax,'260-310 Hz: first plate-mode neighborhood','FontWeight','normal');end
    if b==1,legendHandles=h;end
end
lg=legend(legendHandles,labels,'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
title(tl,'z=0: nullspace response near the target mode','FontWeight','normal');
subtitle(tl,{'Held-out LRM, recorded command coordinates. This is a total observed response, not an isolated flexible-mode preservation ratio.', ...
    'D0/D2 training band: 25-80 Hz. Auxiliary direct-target diagnostic; no additional tuning in the target band.'},'FontSize',10);
files=[files,exportPair(fig,out,'band_diagnostic_z0_NN_magnitude')];clear cleanup;

% 6. Phase is shown only for the scalar N-to-N response, independently unwrap.
fig=newFigure([1440 560]);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for b=1:2
    ax=nexttile(tl);hold(ax,'on');h=gobjects(1,4);
    mask=true(size(f));if b==2,mask=f>=260-1e-6 & f<=310+1e-6;end
    for m=1:4
        Z=response(methods,codes{m});z=reshape(Z(4,4,:),[],1);phase=unwrap(angle(z))*180/pi;
        h(m)=plot(ax,f(mask),phase(mask),'Color',colors.(codes{m}),'LineWidth',lineWidth(codes{m}));
    end
    style(ax);xlabel(ax,'Frequency (Hz)');ylabel(ax,'N to N independently unwrapped phase (deg)');
    if b==1,set(ax,'XScale','log');xlim(ax,[f(1) f(end)]);markTrainingBand(ax);title(ax,'Full measured frequency range','FontWeight','normal');
    else,xlim(ax,[260 310]);title(ax,'260-310 Hz: first plate-mode neighborhood','FontWeight','normal');end
    if b==1,legendHandles=h;end
end
lg=legend(legendHandles,labels,'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
title(tl,'z=0: scalar nullspace phase, held-out LRM','FontWeight','normal');
subtitle(tl,{'Each full frequency trace is unwrapped independently from its own first principal-value sample; no branch alignment to another method.', ...
    'Low-magnitude phase may be noisy. D2 includes its causal common delay. Training band: 25-80 Hz; auxiliary diagnostic.'},'FontSize',10);
files=[files,exportPair(fig,out,'band_diagnostic_z0_NN_phase')];clear cleanup;

% 7-8. Complete D0 matrix, with B0 reference, in fixed recorded coordinates.
axisNames={'Rx ref.','Ry ref.','Z ref.','N'};
B0=response(methods,'B0');D0=response(methods,'D0');
assert(isequal(size(B0),size(D0)) && size(D0,1)==4 && size(D0,2)==4);
for kind=1:2
    fig=newFigure([1470 1190]);cleanup=onCleanup(@()closeOwned(fig));
    tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ax=nexttile(tl);hold(ax,'on');
            a=reshape(B0(y,u,:),[],1);b=reshape(D0(y,u,:),[],1);
            if kind==1,aa=db(a);bb=db(b);else,aa=unwrap(angle(a))*180/pi;bb=unwrap(angle(b))*180/pi;end
            h1=semilogx(ax,f,aa,'Color',colors.B0,'LineWidth',.75);
            h2=semilogx(ax,f,bb,'Color',colors.D0,'LineWidth',1.0);
            style(ax);set(ax,'XScale','log');xlim(ax,[f(1) f(end)]);xticks(ax,[10 100 800]);
            markTrainingBand(ax);title(ax,sprintf('%s to %s',axisNames{u},axisNames{y}),'FontWeight','normal');
            if y<4,ax.XTickLabel=[];end
            if u==1
                if kind==1,ylabel(ax,'Magnitude (dB)');else,ylabel(ax,'Unwrapped phase (deg)');end
            end
            if y==4,xlabel(ax,'Frequency (Hz)');end
        end
    end
    lg=legend([h1 h2],{'B0: nominal static / held-out LRM','D0: band-optimal static / held-out LRM'}, ...
        'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
    if kind==1
        title(tl,'z=0: complete 4-by-4 static-transformed response magnitudes','FontWeight','normal');
        subtitle(tl,{'Rows: output reference axes; columns: input reference axes. dB relative to one recorded-coordinate ratio; command units uncalibrated.', ...
            'Dotted lines: 25-80 Hz training band. D0 is an auxiliary approximate-nullspace diagnostic, not a modal identification model.'},'FontSize',10);
        stem='band_diagnostic_z0_D0_full_magnitude';
    else
        title(tl,'z=0: complete 4-by-4 static-transformed response phases','FontWeight','normal');
        subtitle(tl,{'Each scalar FRF and each method is unwrapped independently on the full frequency grid, without model-based branch alignment.', ...
            'Held-out LRM; noisy low-magnitude phase is retained. Training band: 25-80 Hz. D0 is an auxiliary diagnostic, not the modal estimator.'},'FontSize',10);
        stem='band_diagnostic_z0_D0_full_phase';
    end
    files=[files,exportPair(fig,out,stem)];clear cleanup;
end

% 9. Compact paper figure: the same frozen D0 transfer metrics, no new fit.
fig=newFigure([1160 430]);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
paperTitles={'(a) Logical inputs to N','(b) N to logical outputs'};
for b=1:2
    ax=nexttile(tl);hold(ax,'on');values=zeros(4,2);
    for p=1:4
        values(p,1)=metric(T,positions(p),'own_position','D0',blocks{b},'over_B0');
        scope='z0_transfer';if p==z0,scope='own_position';end
        values(p,2)=metric(T,positions(p),scope,'D0',blocks{b},'over_B0');
    end
    h1=plot(ax,positions,values(:,1),'-o','Color',colors.D0,'LineWidth',1.8, ...
        'MarkerSize',8,'MarkerFaceColor','w');
    h2=plot(ax,positions,values(:,2),'--s','Color',colors.D0,'LineWidth',1.8, ...
        'MarkerSize',8,'MarkerFaceColor','w');
    yline(ax,1,':','Color',colors.B0,'LineWidth',1.1,'HandleVisibility','off');
    style(ax);set(ax,'FontSize',17,'XMinorGrid','off','YMinorGrid','off');
    xticks(ax,positions);xlabel(ax,'Position z (mm)','FontSize',18);
    ylabel(ax,'Response RMS / B0','FontSize',18);
    title(ax,paperTitles{b},'FontWeight','normal','FontSize',18);
    ylim(ax,[0 max(1.05,1.12*max(values,[],'all'))]);
    if b==1,legendHandles=[h1 h2];end
end
lg=legend(legendHandles,{'Fit each position','Reuse z = 0'}, ...
    'Orientation','horizontal','Box','off','FontSize',17);
lg.Layout.Tile='north';
files=[files,exportPair(fig,out,'band_diagnostic_z0_transfer_paper')];clear cleanup;

manifestPath=fullfile(out,'band_calibration_figure_manifest.mat');
files{end+1}=manifestPath;
manifest=struct('files',{files},'source',source,'positions_mm',positions, ...
    'training_band',[25 80],'delay_increment_db',increments, ...
    'notes',{{'All plotted curves are transformed held-out LRM frequency responses; none is a fitted model curve.', ...
    'D0/D1/D2 are auxiliary direct-target diagnostics, not replacements for the structured modal estimator.', ...
    'B1 means the previous modal pilot static directions; B4 adds the previous effective relative delays.', ...
    'Phase is unwrapped independently for each scalar trace; no phase is assigned to vector norms.', ...
    'Recorded sensor/command units remain uncalibrated; logical axes are fixed normalized reference directions.', ...
    'Reusing z=0 at other positions keeps its D0/D2 null vectors and relative delays unchanged.', ...
    'The previously inspected second halves are held-out checks, not fresh blind or independent repeat experiments.'}});
save(manifestPath,'manifest');
assert(all(cellfun(@isfile,files)),'A plotted output file is missing.');
fprintf('Saved %d PNG/PDF figure files and manifest to %s\n',numel(files)-1,out);
end

function v=metric(T,pos,scope,code,block,column)
keep=abs(T.position_mm-pos)<1e-12 & string(T.scope)==string(scope) & ...
    string(T.method)==string(code) & abs(T.band_lo-25)<1e-12 & abs(T.band_hi-80)<1e-12 & ...
    string(T.block)==string(block);
assert(nnz(keep)==1,'Expected one metric for z=%g, %s, %s, %s.',pos,scope,code,block);
values=T.(column);v=values(keep);assert(isscalar(v) && isfinite(v),'Nonfinite or nonscalar plotted metric.');
end
function Z=response(methods,code)
index=find(strcmp({methods.code},code));assert(isscalar(index),'Missing or repeated method %s.',code);
Z=methods(index).H_validation;assert(all(isfinite(Z(:))),'Nonfinite method response.');
end
function fig=newFigure(sz)
fig=figure('Visible','off','Color','w','Position',[60 60 sz]);
end
function style(ax)
set(ax,'Color','w','FontName','Times New Roman','FontSize',11,'LineWidth',.65, ...
    'Box','on','TickDir','out','XMinorGrid','on','YMinorGrid','on', ...
    'GridColor',[.76 .76 .76],'MinorGridColor',[.86 .86 .86], ...
    'GridAlpha',.45,'MinorGridAlpha',.35,'Layer','top');grid(ax,'on');
end
function markTrainingBand(ax)
xline(ax,25,':','Color',[.72 .72 .72],'LineWidth',.8,'HandleVisibility','off');
xline(ax,80,':','Color',[.72 .72 .72],'LineWidth',.8,'HandleVisibility','off');
end
function value=db(z)
value=20*log10(max(abs(z),realmin));
end
function width=lineWidth(code)
if strcmp(code,'B0'),width=.8;elseif strcmp(code,'B4'),width=.95;else,width=1.15;end
end
function paths=exportPair(fig,out,stem)
drawnow;paths={fullfile(out,[stem '.png']),fullfile(out,[stem '.pdf'])};
exportgraphics(fig,paths{1},'Resolution',220,'BackgroundColor','white');
exportgraphics(fig,paths{2},'ContentType','vector','BackgroundColor','white');
end
function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
