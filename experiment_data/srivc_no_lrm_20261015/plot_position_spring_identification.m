function checks=plot_position_spring_identification
% Display all five independently saved position fits. No old-baseline curve.
% Original FRF bins are preserved; dense model lines are analytical evaluations.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'position_spring_identification_20261016');
source=fullfile(out,'position_campaign.mat');
if ~isfile(source)
    error('plot_position_spring_identification:NotReady','Position campaign is not saved yet: %s',source);
end
S=load(source,'campaign');assert(isfield(S,'campaign') && isfield(S.campaign,'results'),'Missing campaign.results.');
records=S.campaign.results;
assert(iscell(records) && numel(records)==5,'campaign.results must contain five position records.');
records=reshape(records,1,[]);
missing=find(cellfun(@isempty,records));
if ~isempty(missing)
    error('plot_position_spring_identification:IncompleteCampaign', ...
        'Position campaign has empty result cells %s. No missing curve has been fabricated.',mat2str(missing));
end
required={'position_mm','f','trainingFRF','validationFRF','selected','delay_matrix_seconds','status','record_warning','source_file'};
for k=1:5
    assert(all(isfield(records{k},required)),'Position result %d lacks required fields.',k);
    assert(all(isfield(records{k}.selected,{'fit','H'})),'Position result %d has no saved fit/prediction.',k);
end
positions=cellfun(@(r)r.position_mm,records);[positions,order]=sort(positions);records=records(order);
assert(max(abs(positions-[-1 -.5 0 .5 1]))<1e-10,'Expected exactly z=-1,-0.5,0,0.5,1 mm.');
palette=[.48 .48 .48;.03 .27 .52;.02 .48 .43;.39 .24 .65;.70 .16 .52];
anomalous=false(1,5);converged=false(1,5);labels=cell(1,5);lineStyles=cell(1,5);
lineWidths=ones(1,5)*1.1;denseF=cell(1,5);denseH=denseF;reconstruction=zeros(1,5);
warnings=cell(1,5);candidateRows={};
for k=1:5
    r=records{k};fit=r.selected.fit;f=r.f(:);
    assert(all(diff(f)>0) && f(1)<=3+1e-6 && f(end)>=800-1e-6, ...
        'Position %g mm does not cover the requested 3-800 Hz band.',r.position_mm);
    assert(isequal(size(r.validationFRF),[numel(f),16]) && all(isfinite(r.validationFRF(:))), ...
        'Position %g mm has nonfinite/mismatched validation FRFs.',r.position_mm);
    warnings{k}=char(string(r.record_warning));
    anomalous(k)=r.position_mm==-1 || strlength(string(r.record_warning))>0;
    converged(k)=strcmp(char(string(r.status)),'converged') && isfield(fit,'converged') && fit.converged;
    labels{k}=sprintf('z = %g mm',r.position_mm);lineStyles{k}='-';
    if anomalous(k)
        labels{k}=[labels{k} ' (anomalous record)'];lineStyles{k}='--';
    end
    if ~converged(k)
        labels{k}=[labels{k} ' (nonconverged)'];lineStyles{k}=':';lineWidths(k)=.8;
    end
    evalSaved=transfer(fit,f,r.delay_matrix_seconds(:).');
    assert(isequal(size(evalSaved),size(r.selected.H)),'Saved prediction shape mismatch.');
    reconstruction(k)=norm(evalSaved-r.selected.H,'fro')/max(norm(r.selected.H,'fro'),realmin);
    assert(reconstruction(k)<1e-10,'Position %g mm model differs from saved selected.H.',r.position_mm);
    fd=logspace(log10(3),log10(800),8000).';
    for m=1:numel(fit.frequency_hz)
        center=fit.frequency_hz(m);span=max(30*center*fit.damping(m),.15);
        lo=max(3,center-span);hi=min(800,center+span);
        if hi>lo,fd=[fd;linspace(lo,hi,600).'];end %#ok<AGROW>
        seed=NaN;
        if isfield(fit,'initial_frequency_hz'),seed=fit.initial_frequency_hz(m);end
        candidateRows(end+1,:)={r.position_mm,m,seed,fit.frequency_hz(m),100*fit.damping(m), ...
            converged(k),anomalous(k)}; %#ok<AGROW>
    end
    fd=unique([3;25;30;800;fd]);hd=transfer(fit,fd,r.delay_matrix_seconds(:).');
    assert(all(isfinite(hd(:))) && all(isfinite(evalSaved(:))),'Position %g mm has nonfinite model FRFs.',r.position_mm);
    denseF{k}=fd;denseH{k}=hd;
end

files={};groups={1:5};groupNames={'all_positions'};
% Always supply the four-record view so the anomalous record cannot hide
% small differences by expanding axes limits in the all-position view.
normal=find(~anomalous);
if ~isempty(normal),groups{end+1}=normal;groupNames{end+1}='nonanomalous_positions';end
for g=1:numel(groups)
    indices=groups{g};
    if g==1
        bands={[3 800],[25 800],[3 30]};bandNames={'full','above25','low'};
    else
        bands={[3 800]};bandNames={'full'};
    end
    anchor=find(positions==0 & converged & ~anomalous,1);
    if isempty(anchor) || ~ismember(anchor,indices),anchor=indices(1);end
    for b=1:numel(bands)
        range=bands{b};
        for phase=[false true]
            fig=figure('Visible','off','Color','w','Position',[40 40 1470 1060]);
            cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
            for y=1:4
                for u=1:4
                    ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');handles=gobjects(0);
                    if phase
                        ia=denseF{anchor}>=range(1) & denseF{anchor}<=range(2);
                        pa=unwrap(angle(denseH{anchor}(ia,ch)))*180/pi;branchStart=pa(1);
                    end
                    for k=indices
                        mask=denseF{k}>=range(1) & denseF{k}<=range(2);
                        if phase
                            values=unwrap(angle(denseH{k}(mask,ch)))*180/pi;
                            values=values+360*round((branchStart-values(1))/360);
                        else
                            values=toDb(denseH{k}(mask,ch));
                        end
                        handles(end+1)=semilogx(ax,denseF{k}(mask),values,'LineStyle',lineStyles{k}, ...
                            'Color',palette(k,:),'LineWidth',lineWidths(k)); %#ok<AGROW>
                    end
                    matrixAxes(ax,y,u,range);
                    if y==1 && u==1,legendHandles=handles;end
                end
            end
            lg=legend(legendHandles,labels(indices),'Orientation','horizontal','Box','off','FontSize',9);
            lg.Layout.Tile='north';xlabel(tl,'Frequency (Hz)','FontSize',12);
            if phase
                name='phase';ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
                phaseNote='Independent unwrap per position, then one constant whole-curve 360-deg shift at the first plotted frequency; no smoothing.';
            else
                name='magnitude';ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);
                phaseNote='Analytical model evaluation includes each position''s fitted input/output delay phases; no previous-baseline curve.';
            end
            if g==1,titleText='Five-position FRF model comparison';
            else,titleText='Position comparison excluding anomalous records';end
            title(tl,titleText,'FontWeight','normal');
            subtitle(tl,{phaseNote, ...
                'Gray dashed: anomalous record. Thin dotted: nonconverged/provisional fit. Curve presence does not establish physical modal identity.'},'FontSize',9);
            stem=sprintf('%s_%s_%s',groupNames{g},bandNames{b},name);
            files=[files,exportFigure(fig,out,stem)]; %#ok<AGROW>
            clear cleanup;
        end
    end
end

% Per-position validation comparisons use ORIGINAL reference frequency bins.
% No model is replaced by a curve from another position or from the old fit.
range=[25 800];referenceGray=[.67 .67 .67];
for k=1:5
    r=records{k};f=r.f(:);ir=f>=range(1) & f<=range(2);
    fd=denseF{k};id=fd>=range(1) & fd<=range(2);
    posTag=strrep(strrep(sprintf('%g',r.position_mm),'-','minus'),'.','p');
    for phase=[false true]
        fig=figure('Visible','off','Color','w','Position',[45 45 1460 1040]);
        cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
        for y=1:4
            for u=1:4
                ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');
                if phase
                    model=unwrap(angle(denseH{k}(id,ch)))*180/pi;
                    experimental=unwrap(angle(r.validationFRF(ir,ch)))*180/pi;
                    experimental=experimental+360*round((model(1)-experimental(1))/360);
                else
                    model=toDb(denseH{k}(id,ch));experimental=toDb(r.validationFRF(ir,ch));
                end
                h1=semilogx(ax,f(ir),experimental,'Color',referenceGray,'LineWidth',.7);
                h2=semilogx(ax,fd(id),model,'Color',palette(k,:),'LineStyle',lineStyles{k},'LineWidth',lineWidths(k));
                matrixAxes(ax,y,u,range);
                if y==1 && u==1,handles=[h1 h2];end
            end
        end
        lg=legend(handles,{'Validation FRF',labels{k}},'Orientation','horizontal','Box','off','FontSize',10);
        lg.Layout.Tile='north';xlabel(tl,'Frequency (Hz)','FontSize',12);
        if phase
            name='phase';ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
            note='Experimental/model phases unwrapped independently; only one constant 360-deg display offset is allowed. No smoothing.';
        else
            name='magnitude';ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);
            note='Validation FRF is shown at its original frequency bins; the fitted model is evaluated analytically on a dense grid.';
        end
        title(tl,sprintf('Position z = %g mm: same-position FRF validation',r.position_mm),'FontWeight','normal');
        subtitle(tl,{note,'The validation split belongs to the same experimental record; nonconverged fits and anomalous records remain explicitly flagged.'},'FontSize',9);
        files=[files,exportFigure(fig,out,sprintf('position_z%s_validation_%s',posTag,name))]; %#ok<AGROW>
        clear cleanup;
    end
end

% Component indices are model coordinates, not a claim of physical matching.
candidates=cell2table(candidateRows,'VariableNames',{'position_mm','candidate_component_index','initial_frequency_hz', ...
    'fitted_frequency_hz','damping_percent','fit_converged','record_anomalous'});
csv=fullfile(out,'position_modal_candidate_coordinates.csv');writetable(candidates,csv);files{end+1}=csv;
checks=struct('files',{files},'positions_mm',positions,'converged',converged,'anomalous_records',anomalous, ...
    'record_warnings',{warnings},'reconstruction_relative_errors',reconstruction, ...
    'dense_grid_counts',cellfun(@numel,denseF),'channel_order','output + 4*(input-1)', ...
    'phase_display','Independent unwrap and one constant whole-curve 360-degree offset only', ...
    'component_matching','Candidate component indices only; no physical mode matching has been established', ...
    'old_baseline_curves','none','smoothing','none');
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
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on','LineWidth',.65, ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');xlim(ax,range);
if range(2)<=30,ticks=[3 5 10 20 30];else,ticks=[3 10 25 100 300 800];end
xticks(ax,ticks(ticks>=range(1) & ticks<=range(2)));
title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
if y<4,xticklabels(ax,{});end
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
