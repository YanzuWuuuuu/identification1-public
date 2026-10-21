function files = plot_reference_calibration_tradeoff(resultFile)
%PLOT_REFERENCE_CALIBRATION_TRADEOFF Compare calibration on a fixed frontend.
% Every response metric uses the old direct FRF of the same check record.
% Gains are operator-norm proxies, not measured actuation/noise RMS.
% No invariance of the full logical 3-by-3 response is claimed: completing
% the modal coordinates changes its flexible contribution. The preserved
% quantity in the construction is the fitted support-space mapping.

base = fileparts(mfilename('fullpath'));
if nargin<1 || isempty(resultFile)
    resultFile = fullfile(base,'full_coordinate_calibration_20261016', ...
        'reference_calibration_evaluation.mat');
end
S = load(resultFile,'e');e = S.e;T = e.metrics;G = e.gains;
out = fileparts(resultFile);
codes = {'B4','R1','R2','R3'};
bands = {'25-80','260-310'};
blocks = {'logical_to_N','N_to_logical'};
values = zeros(4,2,2);
gain = zeros(4,1);
for m=1:4
    for band=1:2
        for block=1:2
            ix = strcmp(T.frontend,'old_direct_FRF') & strcmp(T.method,codes{m}) ...
                & strcmp(T.band_hz,bands{band}) & strcmp(T.block,blocks{block});
            assert(nnz(ix)==1,'Expected one calibration metric for each method/band/block.');
            values(m,block,band) = T.over_old_B4(ix);
        end
    end
    ix = strcmp(G.method,codes{m});
    assert(nnz(ix)==1,'Expected one operator-gain row for %s.',codes{m});
    gain(m) = G.full_input_operator_vs_R1(ix);
end
assert(all(isfinite(values(:)) & values(:)>0) && all(isfinite(gain) & gain>0));
assert(max(abs(reshape(values(1,:,:),[],1)-1))<1e-10, ...
    'Old B4 must be the common cross-response baseline.');
colors = [.04 .29 .55;.06 .45 .28];
fig = figure('Visible','off','Color','w','Position',[80 80 1540 690]);
cleanup = onCleanup(@() close(fig));
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
% The two response panels deliberately share a log range; narrow zooming
% would visually exaggerate the near-unity target-band differences.
upper = 10^ceil(log10(1.35*max(values(:))));
lower = min(.5,min(values(:))/.8);
for band=1:2
    ax = nexttile(tl);hold(ax,'on');
    h = bar(ax,1:4,values(:,:,band),'grouped','BaseValue',lower,'EdgeColor','none');
    for block=1:2,h(block).FaceColor=colors(block,:);end
    style(ax);set(ax,'YScale','log');ylim(ax,[lower upper]);
    xticks(ax,1:4);xticklabels(ax,codes);
    yline(ax,1,':','B4 = 1','HandleVisibility','off', ...
        'LabelHorizontalAlignment','left','Color',[.25 .25 .25]);
    ylabel(ax,'Cross-response RMS / old B4');
    title(ax,[bands{band} ' Hz | fixed old-FRF check'],'FontWeight','normal');
    for block=1:2
        for m=1:4
            text(ax,h(block).XEndPoints(m),values(m,block,band)*1.13, ...
                sprintf('%.3g',values(m,block,band)), ...
                'HorizontalAlignment','center','FontName','Arial','FontSize',10);
        end
    end
    if band==1,legendHandles=h;end
end
ax = nexttile(tl);hold(ax,'on');
h = bar(ax,1:4,gain,.58,'FaceColor',[.47 .20 .61], ...
    'EdgeColor','none','BaseValue',.5);
style(ax);set(ax,'YScale','log');
ylim(ax,[.5 10^ceil(log10(1.35*max(gain)))]);
xticks(ax,1:4);xticklabels(ax,codes);
yline(ax,1,':','R1 = 1','HandleVisibility','off', ...
    'LabelHorizontalAlignment','left','Color',[.25 .25 .25]);
yline(ax,1.25,'--','Budget = 1.25','HandleVisibility','off', ...
    'LabelHorizontalAlignment','left','Color',[.50 .50 .50]);
title(ax,'Complete input-transform operator norm','FontWeight','normal');
ylabel(ax,'Input operator norm / R1');
for m=1:4
    text(ax,h.XEndPoints(m),gain(m)*1.13,sprintf('%.3g',gain(m)), ...
        'HorizontalAlignment','center','FontName','Arial','FontSize',10);
end
lg = legend(legendHandles,{'Logical inputs to N output','N input to logical outputs'}, ...
    'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
title(tl,'Reference-assisted calibration: response tradeoff and input gain','FontWeight','normal');
subtitle(tl,{'B4: previous calibration; R1: reference-fit null spaces; R2: complete target-mode coordinates; R3: gain-limited completion.', ...
    'Exploratory same-record check on one fixed FRF frontend. The input gain is a proxy, not measured hardware effort.', ...
    'Completing coordinates preserves the estimated support mapping, not the full logical response or closed-loop performance.'},'FontSize',10);
png = fullfile(out,'reference_calibration_tradeoff.png');
pdf = fullfile(out,'reference_calibration_tradeoff.pdf');
exportgraphics(fig,png,'Resolution',200);
exportgraphics(fig,pdf,'ContentType','vector');
files = {png,pdf};
clear cleanup;
end

function style(ax)
set(ax,'FontName','Arial','FontSize',11,'Box','on','YMinorGrid','on', ...
    'GridAlpha',.15,'MinorGridAlpha',.07);
grid(ax,'on');
end
