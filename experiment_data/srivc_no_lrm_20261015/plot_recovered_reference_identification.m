function files = plot_recovered_reference_identification(resultFile)
%PLOT_RECOVERED_REFERENCE_IDENTIFICATION Plot the reference-assisted pilot.
% Neither the recovered-reference FRF nor its fixed-pole fit is ground truth.
% The record's second half has been inspected before: this is a same-record
% check, not an independent experiment or a claim of unbiased estimation.

base = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resultFile)
    resultFile = fullfile(base,'full_coordinate_calibration_20261016', ...
        'reference_identification.mat');
end
saved = load(resultFile,'ref');
r = saved.ref;
fAll = r.f(:);
assert(all(isfinite(fAll) & fAll>0) && all(diff(fAll)>0), ...
    'Frequencies must be positive, finite, and strictly increasing.');
assert(size(r.local{2}.FRF,1)==numel(fAll) && size(r.local{2}.FRF,2)==16);
assert(isequal(size(r.models{4}),size(r.local{2}.FRF)));
assert(strcmp(r.model_codes{4},'fixed_poles_reference_moments'));
out = fileparts(resultFile);
use = fAll>=25 & fAll<=800+1e-9;
f = fAll(use);
assert(numel(f)>1,'No usable 25-800 Hz frequency range.');
comparison = {r.local{2}.FRF(use,:),r.models{4}(use,:)};
files = {};

for showPhase = 0:1
    fig = figure('Visible','off','Color','w','Position',[80 80 1480 1160]);
    cleanup = onCleanup(@() close(fig));
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for iy = 1:4
        for iu = 1:4
            ax = nexttile(tl);hold(ax,'on');h = gobjects(1,2);
            ch = iy+4*(iu-1);
            for source = 1:2
                v = comparison{source}(:,ch);
                if showPhase
                    y = independentUnwrap(v);
                else
                    y = db(abs(v));
                end
                if source==1
                    h(source) = plot(ax,f,y,'Color',[.60 .60 .60],'LineWidth',.75);
                else
                    h(source) = plot(ax,f,y,'Color',[.04 .29 .55],'LineWidth',1.10);
                end
            end
            style(ax);set(ax,'XScale','log');
            xlim(ax,[25 800]);xticks(ax,[25 50 100 200 400 800]);
            title(ax,sprintf('u_%d to y_%d',iu,iy),'FontWeight','normal');
            if iy==4,xlabel(ax,'Frequency (Hz)');end
            if iu==1
                if showPhase,ylabel(ax,'Unwrapped phase (deg)');
                else,ylabel(ax,'Magnitude (dB)');end
            end
            if iy==1 && iu==1,legendHandles=h;end
        end
    end
    lg = legend(legendHandles,{'Recovered-reference FRF (second half)', ...
        'Fixed-pole reference-moment fit (first half)'}, ...
        'Orientation','horizontal','Box','off');
    lg.Layout.Tile = 'north';
    title(tl,'Recovered-reference FRF and fixed-pole modal fit','FontWeight','normal');
    if showPhase
        subtitle(tl,{'z = 0; exploratory same-record check. The recovered-reference FRF is an estimate, not ground truth.', ...
            'Independent phase unwrapping from 25 Hz; no model-based branch alignment or removed samples.'},'FontSize',10);
        stem = 'reference_FFR_phase';
    else
        subtitle(tl,{'z = 0; exploratory same-record check. The recovered-reference FRF is an estimate, not ground truth.', ...
            'Physical input/output channels in recorded units; pole locations and channel delays are frozen.'},'FontSize',10);
        stem = 'reference_FFR_magnitude';
    end
    files = [files,exportPair(fig,out,stem)]; %#ok<AGROW>
    clear cleanup;
end

% Equal-parameter-count fits are compared on the same reference-transfer
% discrepancy. Their training objectives differ; this is a pilot comparison.
T = r.metrics;
bands = {'25-80','100-170','260-310','25-800'};
codes = {'fixed_poles_direct_FRF','fixed_poles_reference_FRF', ...
    'fixed_poles_reference_moments'};
colors = [.36 .49 .62;.47 .20 .61;.06 .45 .28];
values = zeros(4,3);
fig = figure('Visible','off','Color','w','Position',[100 100 1200 900]);
cleanup = onCleanup(@() close(fig));
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for band = 1:4
    ax = nexttile(tl);hold(ax,'on');
    for method = 1:3
        useRow = strcmp(T.band_hz,bands{band}) & strcmp(T.model,codes{method});
        assert(nnz(useRow)==1,'Missing or duplicated metric for %s, %s.',codes{method},bands{band});
        values(band,method) = 100*T.check_reference_transfer_error(useRow);
    end
    bh = bar(ax,1:3,values(band,:),.58,'FaceColor','flat','EdgeColor','none');
    bh.CData = colors;
    style(ax);set(ax,'XMinorGrid','off','YMinorGrid','off');
    xticks(ax,1:3);xticklabels(ax,{'Direct FRF','Reference FRF','Reference moments'});
    ylabel(ax,'Relative reference-transfer error (%)');
    title(ax,[bands{band} ' Hz'],'FontWeight','normal');
    top = max(values(band,:));ylim(ax,[0 1.25*top]);
    for method = 1:3
        text(ax,method,values(band,method)+.035*top, ...
            sprintf('%.2f%%',values(band,method)), ...
            'HorizontalAlignment','center','FontName','Arial','FontSize',12);
    end
end
title(tl,'Fixed poles, delays and parameter count: reference-assisted pilot','FontWeight','normal');
subtitle(tl,{'Same-record check: relative RMS of T_{yd} - G T_{ud}, evaluated against the same second-half reference transfers.', ...
    'Direct FRF and reference FRF use the same frozen weights. Reference moments change the objective and weights.', ...
    'The total gain is not attributable solely to the FRF frontend; this is not proof of unbiasedness or successful decoupling.'},'FontSize',10);
files = [files,exportPair(fig,out,'reference_transfer_error_comparison')];
clear cleanup;
end

function y = db(v)
y = 20*log10(max(v,realmin('double')));
y(~isfinite(v)) = NaN;
end

function y = independentUnwrap(v)
y = nan(size(v));
valid = isfinite(real(v)) & isfinite(imag(v));
edges = diff([false;valid(:);false]);
starts = find(edges==1);stops = find(edges==-1)-1;
for k = 1:numel(starts)
    ix = starts(k):stops(k);
    y(ix) = unwrap(angle(v(ix)))*180/pi;
end
end

function style(ax)
set(ax,'FontName','Arial','FontSize',10,'Box','on', ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.15,'MinorGridAlpha',.07);
grid(ax,'on');
end

function files = exportPair(fig,out,stem)
png = fullfile(out,[stem '.png']);pdf = fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);
exportgraphics(fig,pdf,'ContentType','vector');
files = {png,pdf};
end
