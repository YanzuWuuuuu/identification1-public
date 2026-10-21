function checks=plot_repair_notches(resultPath)
% Detailed antiresonance comparison with the same second-half LRM reference.
% Dense model curves include their fitted individual channel delay phases.
% All numerical error metrics and residual curves use original frequency bins.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'error_repair_20261016');
if nargin<1 || isempty(resultPath),resultPath=fullfile(out,'final_channel_delay_srivc.mat');end
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'delay_identification_20261016','common_delay_srivc.mat'),'result');common=S.result;
S=load(fullfile(base,'delay_identification_20261016','channel_delay_srivc.mat'),'result');before=S.result;
S=load(resultPath,'result');after=S.result;
f=common.f(:);band=f>=180 & f<=260;fr=f(band);
LRM=common.lrmValidation;Welch=common.welchValidation;train=common.trainingFRF;
assert(nnz(band)>2 && isequal(size(LRM),size(Welch),size(train)), ...
    'Expected corresponding LRM/Welch reference arrays.');
assert(size(LRM,2)==16 && size(LRM,1)==numel(f),'Expected frequency-by-16 FRF arrays.');
models={before,after};names={'Before','Revised'};cols=[1 6];
orange=[.79 .34 .08];blue=[.03 .25 .45];gray=[.40 .40 .40];lightgray=[.74 .74 .74];
colors={orange,blue};fd=unique([linspace(180,260,6000).';fr]);
hd=cell(2,1);hb=cell(2,1);reconstruction=zeros(2,1);
for m=1:2
    r=models{m};ell=r.delay_matrix_seconds(:).';
    prediction=transfer(r.selected.fit,r.f(:),ell);
    reconstruction(m)=norm(prediction-r.selected.H,'fro')/max(norm(r.selected.H,'fro'),realmin);
    assert(reconstruction(m)<1e-10,'Coefficient/delay mapping disagrees with saved prediction.');
    hd{m}=transfer(r.selected.fit,fd,ell);hb{m}=transfer(r.selected.fit,f,ell);
    assert(all(isfinite(hd{m}(:))) && all(isfinite(hb{m}(:))),'Nonfinite model FRF.');
end
assert(all(isfinite(LRM(:))) && all(isfinite(Welch(:))),'Nonfinite experimental reference FRF.');

fig=figure('Visible','off','Color','w','Position',[50 50 1350 1170]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
for row=1:3
    for col=1:2
        ch=cols(col);ax=nexttile(tl);hold(ax,'on');
        if row==1
            h1=semilogx(ax,fr,toDb(LRM(band,ch)),'Color',gray,'LineWidth',1.1);
            h2=semilogx(ax,fr,toDb(Welch(band,ch)),'--','Color',lightgray,'LineWidth',1.0);
            h3=semilogx(ax,fd,toDb(hd{1}(:,ch)),'Color',orange,'LineWidth',1.25);
            h4=semilogx(ax,fd,toDb(hd{2}(:,ch)),'Color',blue,'LineWidth',1.3);
            ylabel(ax,'FRF magnitude (dB, recorded units)');
            if col==1,handles=[h1 h2 h3 h4];end
        elseif row==2
            phases=unwrap(angle([LRM(band,ch),Welch(band,ch)]),[],1)*180/pi;
            phases(:,2)=phases(:,2)+360*round((phases(1,1)-phases(1,2))/360);
            semilogx(ax,fr,phases(:,1),'Color',gray,'LineWidth',1.1);
            semilogx(ax,fr,phases(:,2),'--','Color',lightgray,'LineWidth',1.0);
            for m=1:2
                pm=unwrap(angle(hd{m}(:,ch)))*180/pi;
                % Only one constant whole-curve 360-degree shift, never a
                % pointwise adjustment that might hide a local phase error.
                pm=pm+360*round((phases(1,1)-pm(1))/360);
                semilogx(ax,fd,pm,'Color',colors{m},'LineWidth',1.25);
            end
            ylabel(ax,'Unwrapped phase (deg)');
        else
            for m=1:2
                residual=hb{m}(band,ch)-LRM(band,ch);
                semilogx(ax,fr,toDb(residual),'Color',colors{m},'LineWidth',1.1);
            end
            ylabel(ax,'Complex residual to LRM (dB, recorded units)');
        end
        style(ax);xlim(ax,[180 260]);xticks(ax,[180 200 220 240 260]);
        title(ax,sprintf('$u_%d\\rightarrow y_%d$',col,col),'Interpreter','latex','FontWeight','normal');
        if row<3,xticklabels(ax,{});else,xlabel(ax,'Frequency (Hz)');end
    end
end
lg=legend(handles,{'Validation LRM','Validation Welch','Before','Revised'}, ...
    'Orientation','horizontal','Box','off','FontSize',11);lg.Layout.Tile='north';
title(tl,'Antiresonance detail before and after model revision','FontWeight','normal');
subtitle(tl,{'Top: FRF magnitude. Middle: locally continuous phase. Bottom: complex residual to the same validation LRM.', ...
    'Phase uses whole-curve constant 360-deg branch offsets only. No smoothing; residuals use original bins.', ...
    'Magnitude units are recorded output / motor-command units; physical input calibration is unconfirmed.'},'FontSize',10);
exportgraphics(fig,fullfile(out,'final_antiresonance_detail.png'),'Resolution',200);
exportgraphics(fig,fullfile(out,'final_antiresonance_detail.pdf'),'ContentType','vector');clear cleanup;

% Center/flank masks are chosen by the TRAINING LRM notch minimum, before
% seeing the revised model. The identical masks assess both models.
rows={};minima={};
for col=1:2
    ch=cols(col);[~,k]=min(abs(train(band,ch)));center=fr(k);
    centerMask=band & abs(f-center)<=5;flankMask=band & abs(f-center)>5;
    [lrmin,lridx]=min(abs(LRM(band,ch)));[wmin,widx]=min(abs(Welch(band,ch)));
    minima(end+1,:)={'Validation LRM',col,col,fr(lridx),20*log10(max(lrmin,realmin)),'original bins'}; %#ok<AGROW>
    minima(end+1,:)={'Validation Welch',col,col,fr(widx),20*log10(max(wmin,realmin)),'original bins'}; %#ok<AGROW>
    regions={band,centerMask,flankMask};
    for m=1:2
        [minimum,im]=min(abs(hd{m}(:,ch)));errors=zeros(1,3);
        for region=1:3
            mask=regions{region};ref=LRM(mask,ch);
            errors(region)=100*norm(hb{m}(mask,ch)-ref)/max(norm(ref),realmin);
        end
        rows(end+1,:)={names{m},col,col,center,fd(im),20*log10(max(minimum,realmin)), ...
            fr(lridx),20*log10(max(lrmin,realmin)),errors(1),errors(2),errors(3)}; %#ok<AGROW>
        minima(end+1,:)={names{m},col,col,fd(im),20*log10(max(minimum,realmin)),'dense model evaluation'}; %#ok<AGROW>
    end
end
metrics=cell2table(rows,'VariableNames',{'model','output_channel','input_channel','training_LRM_center_hz', ...
    'model_minimum_hz','model_minimum_db','validation_LRM_minimum_hz','validation_LRM_minimum_db', ...
    'LRM_band_180_260_complex_error_percent','LRM_center_pm5_complex_error_percent','LRM_remaining_flanks_complex_error_percent'});
minimumTable=cell2table(minima,'VariableNames',{'curve','output_channel','input_channel','minimum_hz','minimum_db','evaluation_grid'});
writetable(metrics,fullfile(out,'final_antiresonance_detail_metrics.csv'));
writetable(minimumTable,fullfile(out,'final_antiresonance_detail_minima.csv'));
checks=struct('reconstruction_relative_error',reconstruction,'channels',[1 1;2 2], ...
    'frequency_band_hz',[180 260],'dense_model_frequency_count',numel(fd), ...
    'center_definition','Training LRM minimum in 180-260 Hz; fixed +/-5 Hz center; remaining band is flanks', ...
    'reference','Second-half LRM for all complex residuals and reported errors; Welch is a separate display reference', ...
    'phase_display','Local independent unwrap followed by a constant whole-curve integer 360-degree shift', ...
    'units','Recorded output per motor command; motor physical unit unconfirmed', ...
    'model_equation','H_jk(s)=exp(-s*ell_jk)*P_jk(s)','metrics',metrics,'minima',minimumTable);
disp(metrics);disp(checks);
end

function H=transfer(fit,f,delaySeconds)
f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if fit.direct,B(:,end)=1;end
P=B*fit.coefficients;assert(size(P,2)==16 && numel(delaySeconds)==16,'Expected 16 delayed FRFs.');
H=exp(-1i*2*pi*f.*delaySeconds(:).').*P;
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function style(ax)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on','LineWidth',.65, ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
