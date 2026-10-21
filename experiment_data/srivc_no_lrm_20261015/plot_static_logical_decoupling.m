function checks=plot_static_logical_decoupling
% Plot saved static logical-axis transforms, without temporal compensation.
% Models use the supplied dense grid; validation stays on original LRM bins.
% Phase display permits independent unwrap and whole-curve 360-deg offsets.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'static_decoupling_20261016');
source=fullfile(out,'decoupled_campaign.mat');
if ~isfile(source)
    error('plot_static_logical_decoupling:NotReady','Saved static-decoupling pack is not ready: %s',source);
end
S=load(source,'pack');pack=S.pack;
assert(all(isfield(pack,{'positions_mm','normal_indices','axis_labels','views','records'})), ...
    'The saved pack lacks required fields.');
normal=reshape(pack.normal_indices,1,[]);positions=reshape(pack.positions_mm,1,[]);
assert(numel(normal)==4 && numel(unique(normal))==4,'Expected four normal-position indices.');
assert(numel(pack.axis_labels)==4 && numel(pack.records)==numel(positions),'Axis/record dimensions disagree.');
assert(max(abs(positions(normal)-[-.5 0 .5 1]))<1e-10,'Normal-position order must be -0.5,0,0.5,1 mm.');
mass=find(strcmp({pack.views.tag},'massbalanced'));geo=find(strcmp({pack.views.tag},'geometric'));
assert(numel(mass)==1 && numel(geo)==1,'Expected massbalanced and geometric views.');
zero=find(positions==0,1);assert(~isempty(zero) && ismember(zero,normal),'The zero-position record is required.');
palette=[.03 .27 .52;.02 .48 .43;.39 .24 .65;.70 .16 .52];
finiteChecks=true(2,4);nearZeroCounts=cell(2,4);exactZeroCounts=cell(2,4);
viewIndices=[mass geo];
for vi=1:2
    v=viewIndices(vi);
    assert(isequal(size(pack.views(v).Ty),[4 4]) && isequal(size(pack.views(v).Tu),[4 4]), ...
        'Static transforms must be ordinary 4-by-4 matrices.');
    assert(all(isfinite(pack.views(v).Ty(:))) && all(isfinite(pack.views(v).Tu(:))),'Nonfinite static matrices.');
    for k=1:4
        r=pack.records{normal(k)};
        assert(~r.anomalous,'normal_indices includes an anomalous record.');
        assert(abs(r.position_mm-positions(normal(k)))<1e-10,'Record position disagrees with pack.positions_mm.');
        f=r.f(:);fd=r.f_dense(:);a=r.views(v);
        assert(all(isfinite(f)) && all(isfinite(fd)) && all(diff(f)>0) && all(diff(fd)>0), ...
            'Frequency grids must be finite and strictly increasing.');
        assert(f(1)<=3+1e-6 && f(end)>=800-1e-6 && fd(1)<=3+1e-6 && fd(end)>=800-1e-6, ...
            'This plot requires measured/model coverage of 3-800 Hz; it will not extrapolate.');
        assert(isequal(size(a.H_model),[4 4 numel(f)]) && isequal(size(a.H_validation),[4 4 numel(f)]) ...
            && isequal(size(a.H_training),[4 4 numel(f)]) && isequal(size(a.H_dense),[4 4 numel(fd)]), ...
            'Expected output-by-input-by-frequency transformed FRFs.');
        finiteChecks(vi,k)=all(isfinite(a.H_model(:))) && all(isfinite(a.H_validation(:))) ...
            && all(isfinite(a.H_training(:))) && all(isfinite(a.H_dense(:)));
        assert(finiteChecks(vi,k),'A normal-position transformed FRF is nonfinite.');
        amplitude=abs(a.H_dense);scale=max(amplitude,[],3);
        nearZeroCounts{vi,k}=sum(amplitude<=1e-10*max(scale,realmin),3);
        exactZeroCounts{vi,k}=sum(amplitude==0,3);
    end
end

files={};
% Four figures: logical 3-by-3 position overlays for both static conventions.
for v=[mass geo]
    for phase=[false true]
        kind=kindName(phase);stem=sprintf('%s_logical_overlay_%s',pack.views(v).tag,kind);
        files=[files,drawGrid(pack,v,normal,3,phase,false,palette,out,stem)]; %#ok<AGROW>
    end
end
% Two figures: full mass-balanced map, including all nullspace cross channels.
for phase=[false true]
    stem=['massbalanced_full_overlay_' kindName(phase)];
    files=[files,drawGrid(pack,mass,normal,4,phase,false,palette,out,stem)]; %#ok<AGROW>
end
% Four figures: measured/identified logical axes at z=0 for both views.
for v=[mass geo]
    for phase=[false true]
        stem=sprintf('%s_z0_logical_validation_%s',pack.views(v).tag,kindName(phase));
        files=[files,drawGrid(pack,v,zero,3,phase,true,palette,out,stem)]; %#ok<AGROW>
    end
end
% One figure: a dedicated nullspace response, with a z=0 validation reference.
files=[files,drawNullspace(pack,mass,normal,zero,palette,out)];
% Two figures: measured/identified full 4-by-4 map at z=0.
for phase=[false true]
    stem=['massbalanced_z0_full_validation_' kindName(phase)];
    files=[files,drawGrid(pack,mass,zero,4,phase,true,palette,out,stem)]; %#ok<AGROW>
end
checks=struct('files',{files},'figure_count',numel(files)/2,'normal_positions_mm',positions(normal), ...
    'excluded_positions_mm',positions(setdiff(1:numel(positions),normal)), ...
    'all_normal_record_values_finite',all(finiteChecks(:)), ...
    'view_order',{{'massbalanced','geometric'}},'finite_by_view_and_position',finiteChecks, ...
    'dense_near_zero_counts_relative_1e10',{nearZeroCounts},'dense_exact_zero_counts',{exactZeroCounts}, ...
    'phase_display','Independent unwrap; constant whole-curve integer 360-degree offsets only; no pointwise alignment', ...
    'near_zero_phase','Phase near zero magnitude can be poorly defined; no such curve values were suppressed or replaced', ...
    'operations','Plots saved static matrix transforms only; no FIR, inverse filter, or delay compensation is applied', ...
    'frequency_band_hz',[3 800],'validation','Original LRM bins, unchanged; dense supplied model grid', ...
    'plot_log_zero','Only log(0) uses the machine-positive floor for numeric display; no signal-level floor or smoothing');
assert(checks.figure_count==13,'Unexpected number of output figures.');
disp(checks);
end

function files=drawGrid(pack,v,indices,naxes,phase,showValidation,palette,out,stem)
normal=reshape(pack.normal_indices,1,[]);indices=reshape(indices,1,[]);
zero=find(pack.positions_mm==0,1);anchor=zero;
if ~ismember(anchor,indices),anchor=indices(1);end
if naxes==3,window=[50 50 1210 970];else,window=[40 40 1460 1040];end
fig=figure('Visible','off','Color','w','Position',window);cleanup=onCleanup(@()closeOwned(fig));
tl=tiledlayout(fig,naxes,naxes,'TileSpacing','compact','Padding','compact');
gray=[.65 .65 .65];modelBlue=[.03 .25 .45];labels={};
for y=1:naxes
    for u=1:naxes
        ax=nexttile(tl);hold(ax,'on');handles=gobjects(0);
        anchorRecord=pack.records{anchor};ad=anchorRecord.f_dense(:);ia=inBand(ad);
        if phase
            av=unwrap(angle(response(anchorRecord.views(v).H_dense,y,u,ia)))*180/pi;branchStart=av(1);
        end
        if showValidation
            r=pack.records{indices(1)};f=r.f(:);iv=inBand(f);
            h=response(r.views(v).H_validation,y,u,iv);
            if phase
                value=unwrap(angle(h))*180/pi;value=value+360*round((branchStart-value(1))/360);
            else,value=toDb(h);end
            handles(end+1)=semilogx(ax,f(iv),value,'Color',gray,'LineWidth',.75); %#ok<AGROW>
            labels={'Validation LRM','Identified model'};
        end
        for k=1:numel(indices)
            index=indices(k);r=pack.records{index};f=r.f_dense(:);id=inBand(f);
            h=response(r.views(v).H_dense,y,u,id);
            if phase
                value=unwrap(angle(h))*180/pi;value=value+360*round((branchStart-value(1))/360);
            else,value=toDb(h);end
            if showValidation,color=modelBlue;else,color=palette(find(normal==index,1),:);end
            lineStyle='-';width=1.05;
            if ~strcmp(char(string(r.status)),'converged'),lineStyle=':';width=.8;end
            handles(end+1)=semilogx(ax,f(id),value,'Color',color,'LineStyle',lineStyle,'LineWidth',width); %#ok<AGROW>
            if ~showValidation
                labels{k}=sprintf('z = %g mm',r.position_mm); %#ok<AGROW>
                if strcmp(lineStyle,':'),labels{k}=[labels{k} ' (nonconverged)'];end
            elseif strcmp(lineStyle,':'),labels{2}='Identified model (nonconverged)';end
        end
        style(ax);title(ax,[pack.axis_labels{u} ' \rightarrow ' pack.axis_labels{y}], ...
            'Interpreter','tex','FontWeight','normal');
        if y<naxes,xticklabels(ax,{});end
        if y==1 && u==1,legendHandles=handles;legendLabels=labels;end
    end
end
lg=legend(legendHandles,legendLabels,'Orientation','horizontal','Box','off','FontSize',11);lg.Layout.Tile='north';
xlabel(tl,'Frequency (Hz)','FontSize',12);
if phase,ylabel(tl,'Unwrapped phase (deg)','FontSize',12);else,ylabel(tl,'Magnitude (dB)','FontSize',12);end
if showValidation
    title(tl,[viewName(pack.views(v).tag) ': z = 0 mm, static only'],'FontWeight','normal');
else
    title(tl,[viewName(pack.views(v).tag) ': static only'],'FontWeight','normal');
    subtitle(tl,'z = -1 mm excluded (anomalous record)','FontSize',10);
end
files=exportFigure(fig,out,stem);clear cleanup;
end

function files=drawNullspace(pack,v,normal,zero,palette,out)
fig=figure('Visible','off','Color','w','Position',[60 60 1130 850]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
gray=[.65 .65 .65];anchor=pack.records{zero};fanchor=anchor.f_dense(:);ia=inBand(fanchor);
for row=1:2
    phase=row==2;ax=nexttile(tl);hold(ax,'on');
    if phase
        value=unwrap(angle(response(anchor.views(v).H_dense,4,4,ia)))*180/pi;branchStart=value(1);
    end
    fr=anchor.f(:);ir=inBand(fr);href=response(anchor.views(v).H_validation,4,4,ir);
    if phase
        value=unwrap(angle(href))*180/pi;value=value+360*round((branchStart-value(1))/360);
    else,value=toDb(href);end
    handles=semilogx(ax,fr(ir),value,'Color',gray,'LineWidth',.8);labels={'Validation LRM, z = 0 mm'};
    for k=1:numel(normal)
        r=pack.records{normal(k)};f=r.f_dense(:);ii=inBand(f);h=response(r.views(v).H_dense,4,4,ii);
        if phase
            value=unwrap(angle(h))*180/pi;value=value+360*round((branchStart-value(1))/360);
        else,value=toDb(h);end
        lineStyle='-';width=1.15;label=sprintf('z = %g mm',r.position_mm);
        if ~strcmp(char(string(r.status)),'converged'),lineStyle=':';width=.8;label=[label ' (nonconverged)'];end
        handles(end+1)=semilogx(ax,f(ii),value,'Color',palette(k,:),'LineStyle',lineStyle,'LineWidth',width); %#ok<AGROW>
        labels{end+1}=label; %#ok<AGROW>
    end
    style(ax);title(ax,[pack.axis_labels{4} ' \rightarrow ' pack.axis_labels{4}], ...
        'Interpreter','tex','FontWeight','normal');
    if phase,ylabel(ax,'Unwrapped phase (deg)');xlabel(ax,'Frequency (Hz)');else,ylabel(ax,'Magnitude (dB)');end
    if row==1,legendHandles=handles;legendLabels=labels;end
end
lg=legend(legendHandles,legendLabels,'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
title(tl,'Mass balanced nullspace: static only','FontWeight','normal');
subtitle(tl,'z = -1 mm excluded (anomalous record)','FontSize',10);
files=exportFigure(fig,out,'massbalanced_nullspace_overlay');clear cleanup;
end

function name=viewName(tag)
if strcmp(tag,'massbalanced'),name='Mass balanced';elseif strcmp(tag,'geometric'),name='Geometric';else,name=char(string(tag));end
end

function name=kindName(phase)
if phase,name='phase';else,name='magnitude';end
end

function value=response(H,y,u,mask)
value=reshape(H(y,u,mask),[],1);
end

function mask=inBand(f)
mask=f>=3-1e-8 & f<=800+1e-8;
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function style(ax)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on','LineWidth',.65, ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');xlim(ax,[3 800]);xticks(ax,[3 10 25 100 300 800]);
end

function files=exportFigure(fig,out,stem)
png=fullfile(out,[stem '.png']);pdf=fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);exportgraphics(fig,pdf,'ContentType','vector');files={png,pdf};
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
