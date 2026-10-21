function checks=plot_antiresonance_weight_diagnosis
% Focused plots of weighting and model-structure diagnostics.
% These plots do not fit, smooth, or modify the saved complex FRF estimates.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'weight_diagnosis_20261016');
S=load(fullfile(out,'weight_sweep.mat'),'audit');A=S.audit;
T=load(fullfile(out,'structure_diagnostics.mat'),'diagnostics');D=T.diagnostics;
f=A.f(:);local=f>=180 & f<=260;fl=f(local);channels=[1 6];
assert(nnz(local)>2,'The requested local band is not covered.');
assert(isequal(D.f(:),f),'Weight and structure diagnostics must share the frequency grid.');
p1=selectModel(A.models,'lrm_p1_floor0.001');
p1zero=selectModel(A.models,'lrm_p1_floor0');
p2=selectModel(A.models,'lrm_p2_floor0.0001');
delayed=selectModel(D.models,'fixed_poles_profile_common_delay');
blue=[.03 .27 .48];orange=[.80 .35 .08];purple=[.55 .16 .56];
gray=[.66 .66 .66];black=[.18 .18 .18];

% Main comparison: show only the two weight exponents to avoid obscuring the
% antiresonance flanks with almost coincident clipped/unclipped p=1 curves.
fig=figure('Visible','off','Color','w','Position',[50 50 1330 860]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    for col=1:2
        ch=channels(col);H=[A.welchValidation(local,ch),A.lrmValidation(local,ch), ...
            p1.H(local,ch),p2.H(local,ch)];
        values=displayValues(H,row==2);ax=nexttile(tl);hold(ax,'on');
        h1=semilogx(ax,fl,values(:,1),'Color',gray,'LineWidth',1.0);
        h2=semilogx(ax,fl,values(:,2),'.','Color',black,'MarkerSize',4);
        h3=semilogx(ax,fl,values(:,3),'Color',blue,'LineWidth',1.3);
        h4=semilogx(ax,fl,values(:,4),'Color',orange,'LineWidth',1.3);
        labelAxes(ax,row,col);
        if row==1 && col==1,handles=[h1 h2 h3 h4];end
    end
end
lg=legend(handles,{'Validation FRF: Welch H1','Validation FRF: LRM', ...
    'SRIVC: inverse magnitude','SRIVC: inverse squared magnitude'}, ...
    'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
title(tl,'Weight choice changes the antiresonance and its flanks','FontWeight','normal');
subtitle(tl,{sprintf('Quadratic weights: exponent 1, floor 0.001 (converged=%s); exponent 2, floor 0.0001 (converged=%s).', ...
    yesNo(p1.fit.converged),yesNo(p2.fit.converged)), ...
    'Phase is locally unwrapped; one constant multiple of 360 deg per whole curve sets the first-point display branch.'},'FontSize',10);
exportgraphics(fig,fullfile(out,'weights_notch_comparison.png'),'Resolution',200);clear cleanup;

% Isolate clipping while retaining the same exponent and training FRF source.
floorDifferences=zeros(2,1);floorMaxDb=zeros(2,1);
fig=figure('Visible','off','Color','w','Position',[60 60 1290 500]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for col=1:2
    ch=channels(col);ax=nexttile(tl);hold(ax,'on');
    h1=semilogx(ax,fl,toDb(A.lrmValidation(local,ch)),'.','Color',black,'MarkerSize',4);
    h2=semilogx(ax,fl,toDb(p1.H(local,ch)),'Color',blue,'LineWidth',1.5);
    h3=semilogx(ax,fl,toDb(p1zero.H(local,ch)),'--','Color',orange,'LineWidth',1.15);
    floorDifferences(col)=100*norm(p1.H(local,ch)-p1zero.H(local,ch))/norm(p1.H(local,ch));
    floorMaxDb(col)=max(abs(toDb(p1.H(local,ch))-toDb(p1zero.H(local,ch))));
    labelAxes(ax,1,col);
    title(ax,sprintf('u_%d to y_%d: complex difference %.3g%%',col,col,floorDifferences(col)),'Interpreter','tex','FontWeight','normal');
    if col==1,handles=[h1 h2 h3];end
end
lg=legend(handles,{'Validation FRF: LRM','Inverse magnitude, floor 0.001', ...
    'Inverse magnitude, no floor'},'Orientation','horizontal','Box','off','FontSize',10);
lg.Layout.Tile='north';title(tl,'Clipping check with the same inverse-magnitude weight','FontWeight','normal');
subtitle(tl,'Curves are saved SRIVC fits; differences are calculated from complex responses before any plotting.','FontSize',10);
exportgraphics(fig,fullfile(out,'weights_floor_check.png'),'Resolution',200);clear cleanup;

% A profiled common delay is a model-structure diagnostic, not a new physical
% delay estimate. All poles remain fixed while real linear residues are fit.
fig=figure('Visible','off','Color','w','Position',[70 70 1330 880]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    for col=1:2
        ch=channels(col);H=[A.welchValidation(local,ch),A.lrmValidation(local,ch), ...
            p2.H(local,ch),delayed.H(local,ch)];
        values=displayValues(H,row==2);ax=nexttile(tl);hold(ax,'on');
        h1=semilogx(ax,fl,values(:,1),'Color',gray,'LineWidth',1.0);
        h2=semilogx(ax,fl,values(:,2),'.','Color',black,'MarkerSize',4);
        h3=semilogx(ax,fl,values(:,3),'Color',orange,'LineWidth',1.3);
        h4=semilogx(ax,fl,values(:,4),'Color',purple,'LineWidth',1.3);
        labelAxes(ax,row,col);
        if row==1 && col==1,handles=[h1 h2 h3 h4];end
    end
end
lg=legend(handles,{'Validation FRF: Welch H1','Validation FRF: LRM', ...
    'Squared-weight SRIVC baseline','Fixed-pole LS + common-delay profile'}, ...
    'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
title(tl,'Does an additional common phase term explain the remaining mismatch?','FontWeight','normal');
subtitle(tl,{sprintf('Diagnostic only: fixed poles, refitted real numerators, one profiled common delay %.3f microseconds.',1e6*delayed.delay_seconds), ...
    'This is not an SRIVC delay estimate or a hardware-delay calibration.', ...
    'Phase display uses one constant 360-deg branch offset per curve; complex responses are unchanged.'},'FontSize',10);
exportgraphics(fig,fullfile(out,'structure_phase_check.png'),'Resolution',200);clear cleanup;

checks=struct('channels',[1 1;2 2],'frequency_band_hz',[180 260], ...
    'floor_change_complex_difference_percent',floorDifferences,'floor_change_max_magnitude_difference_db',floorMaxDb, ...
    'profiled_diagnostic_delay_microseconds',1e6*delayed.delay_seconds, ...
    'delay_interpretation','Fixed-pole LS + common delay profile only; not SRIVC or a physical hardware delay estimate', ...
    'phase_display','Local independent unwrap; constant whole-curve integer 360-degree shift chosen at first local point', ...
    'smoothing','none');
disp(checks);
end

function m=selectModel(models,tag)
ix=find(strcmp({models.tag},tag));assert(isscalar(ix),'Expected exactly one model: %s',tag);m=models(ix);
end

function values=displayValues(H,phase)
if ~phase,values=toDb(H);return;end
values=unwrap(angle(H),[],1)*180/pi;
for c=2:size(values,2)
    values(:,c)=values(:,c)+360*round((values(1,1)-values(1,c))/360);
end
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function labelAxes(ax,row,col)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on', ...
    'LineWidth',.65,'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');xlim(ax,[180 260]);xticks(ax,[180 200 220 240 260]);
title(ax,sprintf('$u_%d\\rightarrow y_%d$',col,col),'Interpreter','latex','FontWeight','normal');
xlabel(ax,'Frequency (Hz)');
if row==1,ylabel(ax,'FRF magnitude (dB)');else,ylabel(ax,'Unwrapped phase (deg)');end
end

function out=yesNo(value)
if value,out='yes';else,out='no';end
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
