function files = plot_full_coordinate_calibration(resultFile)
%PLOT_FULL_COORDINATE_CALIBRATION Inspect complete-coordinate calibration.
% The saved responses are already transformed and share the evaluation
% coordinate scaling. Do not apply Ty_static/Tu_static a second time.
% This is an exploratory same-record check, not a closed-loop validation.
%
% Required result fields:
%   f: N-by-1 frequency in Hz
%   methods: struct array with code, name, H_validation, H_training, H_model
%   Each response has size 4-by-4-by-N, ordered [Rx,Ry,Z,N].
%
% Phase is unwrapped independently for the measured and model responses.
% No model-based phase alignment or point rejection is performed.

base = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resultFile)
    resultFile = fullfile(base,'full_coordinate_calibration_20261016','result.mat');
end
saved = load(resultFile,'result');
assert(isfield(saved,'result'),'Expected a variable named result in %s.',resultFile);
r = saved.result;
f = r.f(:);
assert(all(isfinite(f) & f > 0) && all(diff(f)>0), ...
    'Frequencies must be finite, positive, and strictly increasing.');
assert(isstruct(r.methods),'result.methods must be a struct array.');
out = fileparts(resultFile);
codes = {'B0','B1','B4','F1','F2','F3'};
idx = zeros(1,numel(codes));
for m = 1:numel(codes)
    match = arrayfun(@(v) strcmp(char(v.code),codes{m}),r.methods);
    assert(nnz(match)==1,'Expected exactly one method with code %s.',codes{m});
    idx(m) = find(match);
    checkResponse(r.methods(idx(m)).H_validation,numel(f),'H_validation');
    checkResponse(r.methods(idx(m)).H_model,numel(f),'H_model');
end
methods = r.methods(idx);

% Retain all six methods in the result; show the four decisive comparisons.
colors = [.56 .56 .56; .40 .62 .78; .04 .29 .55; ...
          .68 .48 .73; .47 .20 .61; .06 .45 .28];
lineStyles = {'-',':','-',':','--','-'};
lineWidths = [1.05 1.05 1.35 1.05 1.35 1.55];
plotIdx = [1 3 5 6];
legendNames = cell(1,numel(codes));
for m = 1:numel(codes)
    legendNames{m} = sprintf('%s: %s',codes{m},char(methods(m).name));
end
files = {};
targetHz = 281.24;
if isfield(r,'target_frequency_hz') && isscalar(r.target_frequency_hz)
    targetHz = r.target_frequency_hz;
end

fig = figure('Visible','off','Color','w','Position',[80 80 1540 960]);
cleanup = onCleanup(@() close(fig));
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
blockNames = {'Logical inputs to N output','N input to logical outputs','N input to N output'};
response = zeros(numel(f),3,numel(codes));
for m = 1:numel(codes)
    H = methods(m).H_validation;
    response(:,1,m) = reshape(sqrt(sum(abs(H(4,1:3,:)).^2,2)),[],1);
    response(:,2,m) = reshape(sqrt(sum(abs(H(1:3,4,:)).^2,1)),[],1);
    response(:,3,m) = reshape(abs(H(4,4,:)),[],1);
end
for row = 1:2
    for block = 1:3
        ax = nexttile(tl);
        hold(ax,'on');
        h = gobjects(1,numel(plotIdx));
        for curve = 1:numel(plotIdx)
            m = plotIdx(curve);
            h(curve) = plot(ax,f,db(response(:,block,m)), ...
                'Color',colors(m,:),'LineStyle',lineStyles{m}, ...
                'LineWidth',lineWidths(m));
        end
        style(ax);
        if row == 1
            set(ax,'XScale','log');
            xlim(ax,[max(3,f(1)),min(800,f(end))]);
            xticks(ax,[3 10 30 100 300 800]);
            title(ax,blockNames{block},'FontWeight','normal');
        else
            xlim(ax,[260 310]);
            xticks(ax,260:10:310);
            xline(ax,targetHz,':','HandleVisibility','off','Color',[.35 .35 .35]);
            title(ax,[blockNames{block} ' | target band'],'FontWeight','normal');
            % Fit the limits to the displayed band, not to out-of-band peaks.
            inBand = f>=260 & f<=310;
            shown = db(response(inBand,block,plotIdx));
            shown = shown(isfinite(shown));
            if ~isempty(shown)
                margin = max(2,.07*(max(shown)-min(shown)));
                ylim(ax,[min(shown)-margin,max(shown)+margin]);
            end
        end
        xlabel(ax,'Frequency (Hz)');
        if block < 3
            ylabel(ax,'Cross-response norm (dB re recorded coordinates)');
        else
            ylabel(ax,'Magnitude (dB re recorded coordinates)');
        end
        if row==1 && block==1
            legendHandles = h;
        end
    end
end
lg = legend(legendHandles,legendNames(plotIdx),'Orientation','horizontal', ...
    'NumColumns',2,'Box','off','Interpreter','none','FontSize',9);
lg.Layout.Tile = 'north';
title(tl,'z = 0: full-coordinate calibration and target-channel response','FontWeight','normal');
subtitle(tl,{'Exploratory same-record check: previously inspected second-half FRF; transformations fitted on the first half.', ...
    'Identical coordinate scaling for every method. Cross-response norms combine three logical channels; no closed-loop claim.'},'FontSize',10);
files = [files,exportPair(fig,out,'full_coordinate_cross_responses')];
clear cleanup;

% Detailed element-wise identification check in the gain-limited coordinates.
chosen = methods(strcmp(codes,'F3'));
coordinates = {'Rx','Ry','Z','N'};
comparison = {chosen.H_validation,chosen.H_model};
for showPhase = 0:1
    fig = figure('Visible','off','Color','w','Position',[80 80 1480 1160]);
    cleanup = onCleanup(@() close(fig));
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for iy = 1:4
        for iu = 1:4
            ax = nexttile(tl);
            hold(ax,'on');
            h = gobjects(1,2);
            for source = 1:2
                v = reshape(comparison{source}(iy,iu,:),[],1);
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
            style(ax);
            set(ax,'XScale','log');
            xlim(ax,[max(3,f(1)),min(800,f(end))]);
            xticks(ax,[3 10 30 100 300 800]);
            title(ax,[coordinates{iu} ' to ' coordinates{iy}],'FontWeight','normal');
            if iy==4
                xlabel(ax,'Frequency (Hz)');
            end
            if iu==1
                if showPhase
                    ylabel(ax,'Unwrapped phase (deg)');
                else
                    ylabel(ax,'Magnitude (dB re recorded coordinates)');
                end
            end
            if iy==1 && iu==1
                legendHandles = h;
            end
        end
    end
    lg = legend(legendHandles,{'Second-half FRF','Training-identified model'}, ...
        'Orientation','horizontal','Box','off');
    lg.Layout.Tile = 'north';
    title(tl,'z = 0: F3 gain-limited full-coordinate transformation','FontWeight','normal');
    if showPhase
        subtitle(tl,{'Exploratory same-record check; both responses use the same frozen F3 transformation.', ...
            'Each curve is unwrapped independently along frequency; no alignment to the model and no removed samples.'},'FontSize',10);
        stem = 'full_coordinate_F3_phase';
    else
        subtitle(tl,{'Exploratory same-record check; both responses use the same frozen F3 transformation.', ...
            'Recorded-coordinate response, not a verified physical m/N calibration. No closed-loop performance claim.'},'FontSize',10);
        stem = 'full_coordinate_F3_magnitude';
    end
    files = [files,exportPair(fig,out,stem)]; %#ok<AGROW>
    clear cleanup;
end
end

function checkResponse(H,nf,fieldName)
assert(isnumeric(H) && size(H,1)==4 && size(H,2)==4 && size(H,3)==nf, ...
    '%s must have size 4-by-4-by-numel(f).',fieldName);
end

function y = db(v)
y = 20*log10(max(v,realmin('double')));
y(~isfinite(v)) = NaN;
end

function y = independentUnwrap(v)
% A NaN must not invalidate all later samples; unwrap each finite segment.
% Segment starts retain their principal phase, with no external alignment.
y = nan(size(v));
valid = isfinite(real(v)) & isfinite(imag(v));
edges = diff([false;valid(:);false]);
starts = find(edges==1);
stops = find(edges==-1)-1;
for k = 1:numel(starts)
    segment = starts(k):stops(k);
    y(segment) = unwrap(angle(v(segment)))*180/pi;
end
end

function style(ax)
set(ax,'FontName','Arial','FontSize',10,'Box','on', ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.15,'MinorGridAlpha',.07);
grid(ax,'on');
end

function files = exportPair(fig,out,stem)
png = fullfile(out,[stem '.png']);
pdf = fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);
exportgraphics(fig,pdf,'ContentType','vector');
files = {png,pdf};
end
