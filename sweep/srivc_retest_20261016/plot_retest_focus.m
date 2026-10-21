function checks=plot_retest_focus
% Additional 25-800 Hz identification panels and measured-FRF comparisons.
% Existing 3-800 Hz plots are retained. No data smoothing or interpolation.
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
oldbase=fullfile(root,'experiment_data','srivc_no_lrm_20261015');
S=load(fullfile(out,'campaign.mat'),'campaign');records=S.campaign.records;
ni=find(cellfun(@(r)abs(r.position_mm)<1e-10,records));
assert(isscalar(ni),'Expected exactly one new z=0 record.');new=records{ni};
S=load(fullfile(oldbase,'position_spring_identification_20261016','position_campaign.mat'),'campaign');
oi=find(cellfun(@(r)abs(r.position_mm)<1e-10,S.campaign.results));
assert(isscalar(oi),'Expected exactly one previous z=0 record.');old=S.campaign.results{oi};
S=load(fullfile(oldbase,'static_decoupling_20261016','decoupled_campaign.mat'),'pack');
gi=find(strcmp({S.pack.views.tag},'geometric'));assert(isscalar(gi));
Ty=S.pack.views(gi).Ty;Tu=S.pack.views(gi).Tu;
assert(isequal(size(Ty),[4 4]) && isequal(size(Tu),[4 4]));
assert(all(isfinite(Ty(:))) && all(isfinite(Tu(:))));
validate(new);validate(old);
fit=new.selected.fit;fd=logspace(log10(25),log10(800),6000).';
for m=1:numel(fit.frequency_hz)
    center=fit.frequency_hz(m);width=max(center*fit.damping(m),.005);
    lo=max(25,center-25*width);hi=min(800,center+25*width);
    if hi>lo,fd=[fd;linspace(lo,hi,500).'];end %#ok<AGROW>
end
fd=unique([25;100;800;fd]);hd=transfer(fit,fd,new.delay_matrix_seconds);
hs=transfer(fit,new.f,new.delay_matrix_seconds);
reconstruction=norm(hs-new.selected.H,'fro')/max(norm(new.selected.H,'fro'),realmin);
assert(reconstruction<1e-9 && all(isfinite(hd(:))),'Saved/dense prediction mismatch.');
files={};gray=[.64 .64 .64];blue=[.03 .28 .56];
fr=new.f(:);keep=fr>=25-1e-8 & fr<=800+1e-8;
converged=logical(fit.converged);style='-';status='';
if ~converged,style=':';status=' | SRIVC not converged';end
for phase=[false true]
    fig=figure('Visible','off','Color','w','Position',[35 35 1540 1160]);
    cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'Padding','compact','TileSpacing','compact');
    for y=1:4
        for u=1:4
            ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');
            a=semilogx(ax,fr(keep),value(new.validationFRF(keep,ch),phase,fr(keep)),'Color',gray,'LineWidth',.70);
            b=semilogx(ax,fd,value(hd(:,ch),phase,fd),'Color',blue,'LineWidth',1.10,'LineStyle',style);
            gridStyle(ax,true);xlim(ax,[25 800]);xticks(ax,[25 50 100 200 400 800]);
            title(ax,sprintf('u_%d \\rightarrow y_%d',u,y),'Interpreter','tex','FontWeight','normal','FontSize',12);
            if y<4,xticklabels(ax,{});end
            if y==1 && u==1,handles=[a b];end
        end
    end
    labels={'Check-segment FRF','SRIVC model'};if ~converged,labels{2}=[labels{2} ' (not converged)'];end
    lg=legend(handles,labels,'Orientation','horizontal','Box','off','FontSize',11,'Interpreter','none');lg.Layout.Tile='north';
    xlabel(tl,'Frequency (Hz)','FontSize',12);
    if phase
        ylabel(tl,'Unwrapped phase (deg)','FontSize',12);kind='phase';
        note='Independently unwrapped; common 50 Hz branch convention; no smoothing. Each curve uses its own interpolated 50 Hz phase in [-270,90) deg.';
    else
        ylabel(tl,'Magnitude (dB)','FontSize',12);kind='magnitude';
        note='Recorded coordinate units; original check FRF bins and analytical dense-grid model.';
    end
    title(tl,['New sweep, z = 0 mm: motor-to-sensor identification' status],'FontSize',14,'FontWeight','normal');
    subtitle(tl,{note,'Closed-loop total-input/output fit; first/second segments of the same sweep.'},'FontSize',10);
    files=[files,exportFigure(fig,out,['z0_physical_25_800_' kind])]; %#ok<AGROW>
    clear cleanup;
end

% Four direct comparisons between measured check-segment FRF estimates.
% Neither curve below is a modal prediction.
newLogical=transform(new.validationFRF,Ty,Tu);oldLogical=transform(old.validationFRF,Ty,Tu);
ranges=[100 170;320 380;260 300;470 510];
newCurves={newLogical(:,2),new.validationFRF(:,6),new.validationFRF(:,1),new.validationFRF(:,16)};
oldCurves={oldLogical(:,2),old.validationFRF(:,6),old.validationFRF(:,1),old.validationFRF(:,16)};
names={'Rx \rightarrow Ry (static geometric axes)','u_2 \rightarrow y_2','u_1 \rightarrow y_1','u_4 \rightarrow y_4'};
fig=figure('Visible','off','Color','w','Position',[55 55 1260 820]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
counts=zeros(4,2);
for k=1:4
    ax=nexttile(tl);hold(ax,'on');fn=new.f(:);fo=old.f(:);
    in=fn>=ranges(k,1)-1e-8 & fn<=ranges(k,2)+1e-8;
    io=fo>=ranges(k,1)-1e-8 & fo<=ranges(k,2)+1e-8;
    counts(k,:)=[nnz(io),nnz(in)];assert(all(counts(k,:)>2),'Insufficient FRF bins in a focus band.');
    a=plot(ax,fo(io),value(oldCurves{k}(io),false),'Color',gray,'LineWidth',1.15);
    b=plot(ax,fn(in),value(newCurves{k}(in),false),'Color',blue,'LineWidth',1.20);
    gridStyle(ax,false);xlim(ax,ranges(k,:));
    title(ax,names{k},'Interpreter','tex','FontWeight','normal','FontSize',13);
    xlabel(ax,'Frequency (Hz)');ylabel(ax,'Magnitude (dB)');
    if k==1,handles=[a b];end
end
lg=legend(handles,{'Previous check-segment FRF','New check-segment FRF'},'Orientation','horizontal', ...
    'Box','off','FontSize',11,'Interpreter','none');lg.Layout.Tile='north';
title(tl,'Before/after at z = 0 mm: check-segment FRFs','FontSize',15,'FontWeight','normal');
subtitle(tl,{'Both curves are nonparametric FRF estimates, not modal fits. Original bins; no smoothing or interpolation.', ...
    'The same static geometric matrices are used for Rx to Ry. A peak shift alone does not identify its physical source.'},'FontSize',10);
files=[files,exportFigure(fig,out,'z0_measured_FRF_peak_comparison')];clear cleanup;
checks=struct('files',{files},'figure_count',3,'reconstruction_relative_error',reconstruction, ...
    'converged',converged,'focus_bands_hz',ranges,'focus_original_bin_counts_old_new',counts, ...
    'focus_sources','Previous and new z=0 validationFRF fields; neither curve is a model.', ...
    'logical_maps','Same previous geometric Ty/Tu for both records.', ...
    'phase_display','Independent unwrap, then one constant whole-curve 360-deg shift placing that curve''s own interpolated phase at 50 Hz in [-270,90). No other-curve reference, pointwise tracking, or smoothing.', ...
    'phase_initial_anchor','Fixed 50 Hz; linear interpolation of this curve''s independently unwrapped phase only to select the integer branch offset. No plotted sample is interpolated or replaced.', ...
    'data_processing','No smoothing, interpolation, clipping, or replacement of FRF bins.');
save(fullfile(out,'plot_checks_focus.mat'),'checks');
fprintf('Created 3 additional focus figures; reconstruction %.3g.\n',reconstruction);
end

function validate(r)
f=r.f(:);assert(all(isfinite(f)) && all(diff(f)>0));
assert(isequal(size(r.validationFRF),[numel(f) 16]) && all(isfinite(r.validationFRF(:))));
assert(f(1)<=3+1e-6 && f(end)>=800-1e-6,'Insufficient frequency coverage.');
end

function H=transfer(fit,f,delay)
f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if fit.direct,B(:,end)=1;end
H=(B*fit.coefficients).*exp(-1i*2*pi*f.*delay(:).');
end

function L=transform(H,Ty,Tu)
L=complex(zeros(size(H,1),9));
for k=1:size(H,1),G=Ty(1:3,:)*reshape(H(k,:),4,4)*Tu(:,1:3);L(k,:)=G(:).';end
end

function x=value(H,phase,f)
if phase
    x=unwrap(angle(H))*180/pi;
    % Each curve chooses its own constant branch offset. Subsequent winding
    % differences remain unchanged, including those near response zeros.
    anchor=interp1(f(:),x,50,'linear');
    assert(isfinite(anchor),'The phase curve does not cover the 50 Hz anchor.');
    x=x-360*floor((anchor+270)/360);
else
    x=20*log10(max(abs(H),realmin));
end
end

function gridStyle(ax,logFrequency)
set(ax,'FontName','Times New Roman','FontSize',11,'Box','on','LineWidth',.65, ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.09);
if logFrequency,set(ax,'XScale','log');else,set(ax,'XScale','linear');end
grid(ax,'on');
end

function files=exportFigure(fig,out,stem)
png=fullfile(out,[stem '.png']);pdf=fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);exportgraphics(fig,pdf,'ContentType','vector');files={png,pdf};
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
