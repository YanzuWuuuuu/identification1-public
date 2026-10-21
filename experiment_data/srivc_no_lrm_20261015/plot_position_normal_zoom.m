function checks=plot_position_normal_zoom
% Two focused bands for the four nonanomalous positions. No old baseline.
% Read saved campaign only; never run fitting or modify the active plotter.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'position_spring_identification_20261016');
path=fullfile(out,'position_campaign.mat');
if ~isfile(path),error('plot_position_normal_zoom:NotReady','Campaign is not saved: %s',path);end
S=load(path,'campaign');records=reshape(S.campaign.results,1,[]);
assert(numel(records)==5 && ~any(cellfun(@isempty,records)),'Five saved position records are required.');
positions=cellfun(@(r)r.position_mm,records);[positions,order]=sort(positions);records=records(order);
assert(max(abs(positions-[-1 -.5 0 .5 1]))<1e-10,'Unexpected campaign positions.');
indices=find(positions~=-1);
for k=indices
    assert(strlength(string(records{k}.record_warning))==0,'Position %g mm is also flagged anomalous.',positions(k));
end
palette=[.48 .48 .48;.03 .27 .52;.02 .48 .43;.39 .24 .65;.70 .16 .52];
denseF=cell(1,5);denseH=denseF;labels=cell(1,5);styles=cell(1,5);widths=ones(1,5)*1.1;
converged=false(1,5);reconstruction=NaN(1,5);
for k=indices
    r=records{k};fit=r.selected.fit;ell=r.delay_matrix_seconds(:).';
    prediction=transfer(fit,r.f(:),ell);
    reconstruction(k)=norm(prediction-r.selected.H,'fro')/max(norm(r.selected.H,'fro'),realmin);
    assert(reconstruction(k)<1e-10,'Position %g mm does not reproduce saved selected.H.',positions(k));
    fd=logspace(log10(3),log10(800),8000).';
    for m=1:numel(fit.frequency_hz)
        fc=fit.frequency_hz(m);span=max(30*fc*fit.damping(m),.15);
        lo=max(3,fc-span);hi=min(800,fc+span);
        if hi>lo,fd=[fd;linspace(lo,hi,600).'];end %#ok<AGROW>
    end
    fd=unique([3;25;30;800;fd]);hd=transfer(fit,fd,ell);
    assert(all(isfinite(hd(:))) && all(isfinite(prediction(:))),'Nonfinite FRF at position %g mm.',positions(k));
    denseF{k}=fd;denseH{k}=hd;labels{k}=sprintf('z = %g mm',positions(k));styles{k}='-';
    converged(k)=strcmp(char(string(r.status)),'converged') && fit.converged;
    if ~converged(k),labels{k}=[labels{k} ' (nonconverged)'];styles{k}=':';widths(k)=.8;end
end
anchor=find(positions==0,1);files={};ranges={[25 800],[3 30]};bandNames={'above25','low'};
for b=1:2
    range=ranges{b};
    for phase=[false true]
        fig=figure('Visible','off','Color','w','Position',[40 40 1470 1040]);
        cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
        for y=1:4
            for u=1:4
                ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');handles=gobjects(0);
                if phase
                    ia=denseF{anchor}>=range(1) & denseF{anchor}<=range(2);
                    anchorPhase=unwrap(angle(denseH{anchor}(ia,ch)))*180/pi;
                end
                for k=indices
                    ii=denseF{k}>=range(1) & denseF{k}<=range(2);
                    if phase
                        value=unwrap(angle(denseH{k}(ii,ch)))*180/pi;
                        value=value+360*round((anchorPhase(1)-value(1))/360);
                    else
                        value=20*log10(max(abs(denseH{k}(ii,ch)),realmin));
                    end
                    handles(end+1)=semilogx(ax,denseF{k}(ii),value,'Color',palette(k,:), ...
                        'LineStyle',styles{k},'LineWidth',widths(k)); %#ok<AGROW>
                end
                set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on', ...
                    'LineWidth',.65,'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
                grid(ax,'on');xlim(ax,range);
                if b==1,xticks(ax,[25 100 300 800]);else,xticks(ax,[3 5 10 20 30]);end
                title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
                if y<4,xticklabels(ax,{});end
                if y==1 && u==1,legendHandles=handles;end
            end
        end
        lg=legend(legendHandles,labels(indices),'Orientation','horizontal','Box','off','FontSize',11);
        lg.Layout.Tile='north';xlabel(tl,'Frequency (Hz)','FontSize',12);
        if phase
            name='phase';ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
            note='Independent unwrap per position; only a constant whole-curve 360-deg shift at the first plotted frequency. No smoothing.';
        else
            name='magnitude';ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);
            note='Each position uses its own fitted modal model and input/output delays. Dense analytical evaluation; no old-baseline curve.';
        end
        title(tl,sprintf('Four nonanomalous positions: %g-%g Hz',range(1),range(2)),'FontWeight','normal');
        subtitle(tl,{note,'The anomalous z = -1 mm record is excluded from this focused comparison.'},'FontSize',10);
        stem=sprintf('nonanomalous_positions_%s_%s',bandNames{b},name);
        png=fullfile(out,[stem '.png']);pdf=fullfile(out,[stem '.pdf']);
        exportgraphics(fig,png,'Resolution',200);exportgraphics(fig,pdf,'ContentType','vector');
        files=[files,{png,pdf}]; %#ok<AGROW>
        clear cleanup;
    end
end
checks=struct('files',{files},'positions_mm',positions(indices),'converged',converged(indices), ...
    'reconstruction_relative_errors',reconstruction(indices),'ranges_hz',[25 800;3 30], ...
    'phase_display','Independent unwrap plus constant whole-curve integer 360-degree offsets only', ...
    'excluded_position_mm',-1,'smoothing','none');
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

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
