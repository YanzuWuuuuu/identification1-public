function checks=plot_retest_results(mode)
% Plot the new SRIVC sweep campaign without changing its data or estimates.
% mode: 'overview' (6), 'per_position' (20), 'all' (26), 'phase_only' (12).
% Phase curves are independently unwrapped, then shifted by one whole-curve
% integer multiple of 360 deg so their own phase at 50 Hz is in [-270,90).
% Only the branch-selection anchor is linearly interpolated at 50 Hz.
% The rule does not reference another curve or track subsequent zero crossings.
if nargin<1,mode='overview';end
mode=char(string(mode));
assert(any(strcmp(mode,{'overview','per_position','all','phase_only'})),'Unknown plot mode.');
phaseOnly=strcmp(mode,'phase_only');
phaseModes=[false true];if phaseOnly,phaseModes=true;end
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
oldbase=fullfile(root,'experiment_data','srivc_no_lrm_20261015');
source=fullfile(out,'campaign.mat');
S=load(source,'campaign');records=S.campaign.records;
assert(iscell(records) && numel(records)==5,'Expected five campaign.records.');
records=reshape(records,1,[]);positions=cellfun(@(r)r.position_mm,records);
[positions,order]=sort(positions);records=records(order);
assert(max(abs(positions-[-.6 -.3 0 .5 1]))<1e-8,'Unexpected new position order.');
T=load(fullfile(oldbase,'static_decoupling_20261016','decoupled_campaign.mat'),'pack');
geo=find(strcmp({T.pack.views.tag},'geometric'));
assert(isscalar(geo) && geo==2,'Expected the existing geometric view at index 2.');
Ty=T.pack.views(geo).Ty;Tu=T.pack.views(geo).Tu;
assert(isequal(size(Ty),[4 4]) && isequal(size(Tu),[4 4]),'Expected full static 4-by-4 maps.');
assert(all(isfinite(Ty(:))) && all(isfinite(Tu(:))),'Static map is nonfinite.');
logicalLabels={'Rx','Ry','Z'};
if isfield(T.pack,'axis_labels')
    assert(isequal(string(T.pack.axis_labels(1:3)),string(logicalLabels)), ...
        'The saved logical axis order differs from Rx, Ry, Z.');
end
palette=[.04 .28 .57;.00 .47 .43;.39 .25 .66;.70 .15 .46;.23 .31 .35];
new=cell(1,5);reconstruction=zeros(1,5);status=cell(1,5);
for k=1:5
    new{k}=prepareRecord(records{k},Ty,Tu);
    reconstruction(k)=new{k}.reconstruction;status{k}=new{k}.status;
end
files={};
if any(strcmp(mode,{'overview','all','phase_only'}))
    for logical=[false true]
        for phase=phaseModes
            stem=sprintf('new_positions_%s_%s',viewTag(logical),kindTag(phase));
            files=[files,drawOverlay(new,logical,phase,palette,logicalLabels,out,stem)]; %#ok<AGROW>
        end
    end
    if ~phaseOnly
        O=load(fullfile(oldbase,'position_spring_identification_20261016','position_campaign.mat'),'campaign');
        oi=find(cellfun(@(r)abs(r.position_mm)<1e-10,O.campaign.results));
        ni=find(abs(positions)<1e-10);
        assert(isscalar(oi) && isscalar(ni),'Need one old and one new z=0 record.');
        old=prepareRecord(O.campaign.results{oi},Ty,Tu);
        for logical=[false true]
            stem=sprintf('z0_before_after_%s_magnitude',viewTag(logical));
            files=[files,drawComparison(new{ni},old,logical,logicalLabels,out,stem)]; %#ok<AGROW>
        end
    end
end
if any(strcmp(mode,{'per_position','all','phase_only'}))
    for k=1:5
        for logical=[false true]
            for phase=phaseModes
                stem=sprintf('position_%s_%s_%s',safeTag(records{k}),viewTag(logical),kindTag(phase));
                files=[files,drawValidation(new{k},logical,phase,logicalLabels,out,stem)]; %#ok<AGROW>
            end
        end
    end
end
checks=struct('mode',mode,'files',{files},'figure_count',numel(files)/2, ...
    'positions_mm',positions,'status',{status},'reconstruction_relative_error',reconstruction, ...
    'static_output_map',Ty,'static_input_map',Tu,'logical_axis_order',{logicalLabels}, ...
    'frequency_band_hz',[3 800],'source_file',source, ...
    'phase_display','Independent unwrap, followed by one constant whole-curve 360-deg shift placing that curve''s own interpolated phase at 50 Hz in [-270,90). No other-curve reference, pointwise tracking, or smoothing.', ...
    'phase_initial_anchor','Fixed 50 Hz; linear interpolation of this curve''s independently unwrapped phase only to select the integer branch offset. No plotted sample is interpolated or replaced.', ...
    'reference_display','Original validation FRF bins, with no interpolation, smoothing, clipping, or replacement.', ...
    'model_display','Analytical modal evaluation on a dense grid with additional points around every fitted pole.', ...
    'logical_operation','Same previous geometric static maps applied to every new and old FRF; no FIR or delay compensation.', ...
    'validation_scope','First/second segments of the same sweep; not an independent repeat experiment.');
save(fullfile(out,['plot_checks_' mode '.mat']),'checks');
fprintf('Created %d figures (%s); reconstruction max %.3g.\n',checks.figure_count,mode,max(reconstruction));
end

function r=prepareRecord(input,Ty,Tu)
required={'position_mm','f','trainingFRF','validationFRF','selected','delay_matrix_seconds','status'};
assert(all(isfield(input,required)),'Incomplete input record.');
f=input.f(:);fit=input.selected.fit;
assert(all(isfinite(f)) && all(diff(f)>0) && f(1)<=3+1e-6 && f(end)>=800-1e-6, ...
    'Record does not cover 3-800 Hz with increasing frequency bins.');
assert(isequal(size(input.validationFRF),[numel(f) 16]) && all(isfinite(input.validationFRF(:))), ...
    'Invalid measured FRF array.');
assert(isequal(size(input.delay_matrix_seconds),[4 4]) && all(isfinite(input.delay_matrix_seconds(:))), ...
    'Invalid delay matrix.');
saved=transfer(fit,f,input.delay_matrix_seconds);
assert(isequal(size(saved),size(input.selected.H)),'Saved model shape mismatch.');
err=norm(saved-input.selected.H,'fro')/max(norm(input.selected.H,'fro'),realmin);
assert(err<1e-9,'Dense model formula does not reconstruct selected.H: %.3g.',err);
fd=logspace(log10(3),log10(800),7000).';
for m=1:numel(fit.frequency_hz)
    center=fit.frequency_hz(m);width=max(center*fit.damping(m),.005);
    lo=max(3,center-25*width);hi=min(800,center+25*width);
    if hi>lo,fd=[fd;linspace(lo,hi,500).'];end %#ok<AGROW>
end
fd=unique([3;25;100;800;fd]);hd=transfer(fit,fd,input.delay_matrix_seconds);
assert(all(isfinite(hd(:))),'Nonfinite dense model.');
keep=f>=3-1e-8 & f<=800+1e-8;
r=struct('position_mm',input.position_mm,'f',f(keep),'fd',fd, ...
    'physical_check',input.validationFRF(keep,:),'physical_model',hd, ...
    'logical_check',logicalTransform(input.validationFRF(keep,:),Ty,Tu), ...
    'logical_model',logicalTransform(hd,Ty,Tu),'reconstruction',err, ...
    'status',char(string(input.status)),'converged',logical(fit.converged));
% The inner iteration flag is the relevant convergence criterion when delays
% are intentionally held fixed; the full status is retained in the subtitle.
end

function H=transfer(fit,f,delay)
f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';
B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if fit.direct,B(:,end)=1;end
H=(B*fit.coefficients).*exp(-1i*2*pi*f.*delay(:).');
assert(size(H,2)==16,'Expected sixteen physical responses.');
end

function L=logicalTransform(H,Ty,Tu)
L=complex(zeros(size(H,1),9));
for k=1:size(H,1)
    G=reshape(H(k,:),4,4);A=Ty(1:3,:)*G*Tu(:,1:3);L(k,:)=A(:).';
end
end

function files=drawValidation(r,isLogical,isPhase,axisLabels,out,stem)
n=4-double(isLogical);[fig,tl]=newFigure(n);cleanup=onCleanup(@()closeOwned(fig));
[model,reference]=getArrays(r,isLogical);blue=[.03 .28 .56];gray=[.65 .65 .65];
for y=1:n
    for u=1:n
        ch=y+n*(u-1);ax=nexttile(tl);hold(ax,'on');
        a=semilogx(ax,r.f,values(reference(:,ch),isPhase,r.f),'Color',gray,'LineWidth',.70);
        b=semilogx(ax,r.fd,values(model(:,ch),isPhase,r.fd),'Color',blue,'LineWidth',1.10, ...
            'LineStyle',modelStyle(r));
        axesStyle(ax,y,u,n,isLogical,axisLabels);
        if y==1 && u==1,handles=[a b];end
    end
end
label='SRIVC model';if ~r.converged,label=[label ' (not converged)'];end
setLegend(handles,{'Check-segment FRF',label});
setGlobalLabels(tl,isPhase);
title(tl,sprintf('%s: z = %g mm%s',viewName(isLogical),r.position_mm,statusTitle(r)), ...
    'FontWeight','normal','FontSize',14);
subtitle(tl,{sprintf('Closed-loop sweep; static geometric maps on logical axes. Status: %s.',r.status),phaseNote(isPhase)},'FontSize',10,'Interpreter','none');
files=exportFigure(fig,out,stem);clear cleanup;
end

function files=drawOverlay(records,isLogical,isPhase,palette,axisLabels,out,stem)
n=4-double(isLogical);[fig,tl]=newFigure(n);cleanup=onCleanup(@()closeOwned(fig));
labels=cell(1,numel(records));nonconverged=[];
for k=1:numel(records)
    labels{k}=sprintf('z = %g mm',records{k}.position_mm);
    if ~records{k}.converged,labels{k}=[labels{k} ' (not converged)'];nonconverged(end+1)=records{k}.position_mm;end %#ok<AGROW>
end
for y=1:n
    for u=1:n
        ch=y+n*(u-1);ax=nexttile(tl);hold(ax,'on');handles=gobjects(1,numel(records));
        for k=1:numel(records)
            r=records{k};model=getArrays(r,isLogical);
            handles(k)=semilogx(ax,r.fd,values(model(:,ch),isPhase,r.fd),'Color',palette(k,:), ...
                'LineWidth',1.05,'LineStyle',modelStyle(r));
        end
        axesStyle(ax,y,u,n,isLogical,axisLabels);
        if y==1 && u==1,legendHandles=handles;end
    end
end
setLegend(legendHandles,labels);setGlobalLabels(tl,isPhase);
note='';if ~isempty(nonconverged),note=sprintf(' | Not converged at z = %s mm',mat2str(nonconverged));end
title(tl,['New sweeps: ' viewName(isLogical) note],'FontWeight','normal','FontSize',14,'Interpreter','none');
subtitle(tl,{'Same static geometric maps for all positions; curves show separately fitted models.',phaseNote(isPhase)},'FontSize',10);
files=exportFigure(fig,out,stem);clear cleanup;
end

function files=drawComparison(new,old,isLogical,axisLabels,out,stem)
n=4-double(isLogical);[fig,tl]=newFigure(n);cleanup=onCleanup(@()closeOwned(fig));
[newModel,newCheck]=getArrays(new,isLogical);oldModel=getArrays(old,isLogical);
for y=1:n
    for u=1:n
        ch=y+n*(u-1);ax=nexttile(tl);hold(ax,'on');
        a=semilogx(ax,new.f,values(newCheck(:,ch),false),'Color',[.66 .66 .66],'LineWidth',.7);
        b=semilogx(ax,new.fd,values(newModel(:,ch),false),'Color',[.03 .28 .56], ...
            'LineWidth',1.15,'LineStyle',modelStyle(new));
        c=semilogx(ax,old.fd,values(oldModel(:,ch),false),'Color',[.50 .27 .66], ...
            'LineWidth',1.0,'LineStyle','--');
        axesStyle(ax,y,u,n,isLogical,axisLabels);
        if y==1 && u==1,handles=[a b c];end
    end
end
labels={'New check-segment FRF','New SRIVC model','Previous SRIVC model'};
if ~new.converged,labels{2}=[labels{2} ' (not converged)'];end
if ~old.converged,labels{3}=[labels{3} ' (not converged)'];end
setLegend(handles,labels);setGlobalLabels(tl,false);
note='';if ~new.converged || ~old.converged,note=' | Includes nonconverged fit';end
title(tl,['Before/after, z = 0 mm: ' viewName(isLogical) note],'FontWeight','normal','FontSize',14,'Interpreter','none');
subtitle(tl,{'Same position and the same static geometric maps; no temporal compensation.', ...
    'Changes include all differences between the two acquisitions; frequency shifts alone do not establish a physical mode identity.'},'FontSize',10);
files=exportFigure(fig,out,stem);clear cleanup;
end

function [fig,tl]=newFigure(n)
if n==4,window=[30 30 1540 1160];else,window=[40 40 1340 1040];end
fig=figure('Visible','off','Color','w','Position',window);
tl=tiledlayout(fig,n,n,'TileSpacing','compact','Padding','compact');
end

function [model,reference]=getArrays(r,isLogical)
if isLogical,model=r.logical_model;reference=r.logical_check;
else,model=r.physical_model;reference=r.physical_check;end
end

function result=values(H,isPhase,f)
if isPhase
    result=unwrap(angle(H))*180/pi;
    % Choose a common branch from this curve alone. Keep all later relative
    % winding and zero-crossing differences from independent unwrap intact.
    anchor=interp1(f(:),result,50,'linear');
    assert(isfinite(anchor),'The phase curve does not cover the 50 Hz anchor.');
    result=result-360*floor((anchor+270)/360);
else
    result=20*log10(max(abs(H),realmin));
end
end

function axesStyle(ax,y,u,n,isLogical,labels)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on', ...
    'LineWidth',.6,'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.09);
grid(ax,'on');xlim(ax,[3 800]);xticks(ax,[3 10 25 100 300 800]);
if isLogical,text=[labels{u} ' \rightarrow ' labels{y}];
else,text=sprintf('u_%d \\rightarrow y_%d',u,y);end
title(ax,text,'Interpreter','tex','FontWeight','normal','FontSize',12);
if y<n,xticklabels(ax,{});end
end

function setLegend(handles,labels)
lg=legend(handles,labels,'Orientation','horizontal','Box','off','FontSize',10,'Interpreter','none');
lg.Layout.Tile='north';
end

function setGlobalLabels(tl,isPhase)
xlabel(tl,'Frequency (Hz)','FontSize',12);
if isPhase,ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
else,ylabel(tl,'Magnitude (dB)','FontSize',12);end
end

function text=phaseNote(isPhase)
if isPhase,text='Independently unwrapped; common 50 Hz branch convention; no smoothing. Each curve uses its own interpolated 50 Hz phase in [-270,90) deg.';
else,text='Recorded coordinate units; check FRFs remain on their original bins. Dense lines are analytical model evaluations.';end
end

function text=statusTitle(r)
if r.converged,text='';else,text=' | SRIVC not converged';end
end

function style=modelStyle(r)
if r.converged,style='-';else,style=':';end
end

function tag=safeTag(r)
if isfield(r,'tag'),tag=char(string(r.tag));
else,tag=strrep(strrep(sprintf('z%g',r.position_mm),'-','minus'),'.','p');end
assert(isempty(regexp(tag,'[^a-zA-Z0-9_\-]','once')),'Unsafe plot filename tag.');
end

function name=viewTag(isLogical)
if isLogical,name='logical';else,name='physical';end
end

function name=viewName(isLogical)
if isLogical,name='logical axes, static geometry';else,name='motor-to-sensor FRFs';end
end

function name=kindTag(isPhase)
if isPhase,name='phase';else,name='magnitude';end
end

function files=exportFigure(fig,out,stem)
png=fullfile(out,[stem '.png']);pdf=fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);exportgraphics(fig,pdf,'ContentType','vector');
files={png,pdf};
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
