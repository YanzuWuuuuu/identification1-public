function checks=plot_frf_preprocessing_comparison
% Display the FRF-entry/weight ablation without selecting successful cases.
% All comparison curves are saved estimates. No smoothing or refitting here.
out=fileparts(mfilename('fullpath'));
S=load(fullfile(out,'frf_preprocessing_comparison.mat'),'comparison');C=S.comparison;
f=C.f(:);assert(size(C.welchValidation,2)==16,'Expected 4-output-by-4-input FRFs.');
assert(all(diff(f)>0),'Frequency grid must be increasing.');
blue=[.03 .27 .48];purple=[.55 .16 .56];gray=[.64 .64 .64];darkgray=[.27 .27 .27];teal=[.04 .48 .47];orange=[.80 .35 .08];
uniformWelch=selectModel(C,'welch_uniform');inverseWelch=selectModel(C,'welch_inverse_magnitude');
uniformLRM=selectModel(C,'lrm_uniform');inverseLRM=selectModel(C,'lrm_inverse_magnitude');
cols=[1 6];channelInputs=[1 2];channelOutputs=[1 2];
status=sprintf('Converged: Welch uniform=%s, inverse-magnitude=%s; LRM uniform=%s, inverse-magnitude=%s', ...
    yesNo(uniformWelch.fit.converged),yesNo(inverseWelch.fit.converged), ...
    yesNo(uniformLRM.fit.converged),yesNo(inverseLRM.fit.converged));

% Rows isolate the training FRF source; within each row only weights change.
fig=figure('Visible','off','Color','w','Position',[50 50 1300 820]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    if row==1,uniform=uniformWelch;inverse=inverseWelch;source='Welch H1';
    else,uniform=uniformLRM;inverse=inverseLRM;source='LRM';end
    for col=1:2
        ch=cols(col);ax=nexttile(tl);hold(ax,'on');
        p1=semilogx(ax,f,toDb(C.welchValidation(:,ch)),'Color',gray,'LineWidth',1.15);
        p2=semilogx(ax,f,toDb(uniform.H(:,ch)),'Color',blue,'LineWidth',1.3);
        p3=semilogx(ax,f,toDb(inverse.H(:,ch)),'Color',purple,'LineWidth',1.3);
        style(ax);xlim(ax,[180 260]);xticks(ax,[180 200 220 240 260]);
        title(ax,sprintf('%s training FRF: u_%d to y_%d',source,channelInputs(col),channelOutputs(col)),'FontWeight','normal');
        xlabel(ax,'Frequency (Hz)');ylabel(ax,'FRF magnitude (dB)');
        if row==1 && col==1,legendHandles=[p1 p2 p3];end
    end
end
lg=legend(legendHandles,{'Same validation FRF (Welch H1)','SRIVC: uniform weight', ...
    'SRIVC: inverse-magnitude weight'},'Orientation','horizontal','Box','off','FontSize',10);
lg.Layout.Tile='north';title(tl,'Antiresonance: changing the FRF input and objective weight','FontWeight','normal');
subtitle(tl,status,'FontSize',10);
exportgraphics(fig,fullfile(out,'frf_preprocessing_antiresonance.png'),'Resolution',200);clear cleanup;

% Three identification routes on a fixed local window. Every phase curve is
% unwrapped locally. Only ONE constant 360-degree offset per whole model
% curve is chosen at the first local frequency, preserving all phase errors.
local=f>=180 & f<=260;flocal=f(local);
assert(nnz(local)>2,'The route comparison needs frequencies in 180-260 Hz.');
fig=figure('Visible','off','Color','w','Position',[55 55 1330 850]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    for col=1:2
        ch=cols(col);ax=nexttile(tl);hold(ax,'on');
        allH=[C.welchValidation(local,ch),C.rawH(local,ch), ...
            inverseWelch.H(local,ch),inverseLRM.H(local,ch)];
        if row==1
            values=toDb(allH);
        else
            values=unwrap(angle(allH),[],1)*180/pi;
            for curve=2:4
                values(:,curve)=values(:,curve)+360*round((values(1,1)-values(1,curve))/360);
            end
        end
        p1=semilogx(ax,flocal,values(:,1),'Color',gray,'LineWidth',1.15);
        p2=semilogx(ax,flocal,values(:,2),'Color',blue,'LineWidth',1.3);
        p3=semilogx(ax,flocal,values(:,3),'Color',orange,'LineWidth',1.2);
        p4=semilogx(ax,flocal,values(:,4),'Color',purple,'LineWidth',1.2);
        style(ax);xlim(ax,[180 260]);xticks(ax,[180 200 220 240 260]);
        title(ax,sprintf('$u_%d\\rightarrow y_%d$',channelInputs(col),channelOutputs(col)), ...
            'Interpreter','latex','FontWeight','normal');
        xlabel(ax,'Frequency (Hz)');
        if row==1,ylabel(ax,'FRF magnitude (dB)');else,ylabel(ax,'Unwrapped phase (deg)');end
        if row==1 && col==1,legendHandles=[p1 p2 p3 p4];end
    end
end
lg=legend(legendHandles,{'Experimental FRF (Welch H1)','Direct SRIVC', ...
    'Welch + weighted SRIVC','LRM + weighted SRIVC'}, ...
    'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
title(tl,'Three identification routes near the antiresonance','FontWeight','normal');
subtitle(tl,{'Weighted routes use inverse-magnitude weights; the gray FRF comes from the validation record.', ...
    'Phase: one constant multiple of 360 deg per curve, chosen at the first local frequency; complex responses are unchanged.'},'FontSize',10);
exportgraphics(fig,fullfile(out,'frf_preprocessing_routes.png'),'Resolution',200);clear cleanup;

% Full matrix retains the chosen LRM/inverse-magnitude case even if its score
% gate failed. Status is explicit in both the magnitude and phase figures.
for phase=[false true]
    fig=figure('Visible','off','Color','w','Position',[50 50 1450 1040]);
    cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ch=y+4*(u-1);ax=nexttile(tl);hold(ax,'on');
            if phase
                raw=unwrap(angle(C.rawH(:,ch)))*180/pi;
                candidate=unwrap(angle(inverseLRM.H(:,ch)))*180/pi;
                experimental=unwrap(angle(C.welchValidation(:,ch)))*180/pi;
                % Display-only integer branch selection, preserving each
                % complex FRF phase exactly. Do not smooth the gray trace.
                experimental=experimental+360*round((raw-experimental)/360);
                p1=semilogx(ax,f,experimental,'.','Color',gray,'MarkerSize',3);
            else
                raw=toDb(C.rawH(:,ch));candidate=toDb(inverseLRM.H(:,ch));
                experimental=toDb(C.welchValidation(:,ch));
                p1=semilogx(ax,f,experimental,'Color',gray,'LineWidth',.7);
            end
            p2=semilogx(ax,f,raw,'Color',blue,'LineWidth',1.05);
            p3=semilogx(ax,f,candidate,'Color',purple,'LineWidth',1.0);
            style(ax);xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);
            title(ax,sprintf('$u_%d\\rightarrow y_%d$',u,y),'Interpreter','latex','FontWeight','normal');
            if y<4,xticklabels(ax,{});end
            if y==1 && u==1,legendHandles=[p1 p2 p3];end
        end
    end
    lg=legend(legendHandles,{'Experimental FRF (validation Welch H1)', ...
        'Raw U/Y SRIVC','LRM FRF + inverse-magnitude SRIVC'}, ...
        'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
    xlabel(tl,'Frequency (Hz)','FontSize',12);
    title(tl,sprintf('LRM/inverse-magnitude fit: converged=%s, IV score=%.3g, iterations=%d', ...
        yesNo(inverseLRM.fit.converged),inverseLRM.fit.score,inverseLRM.fit.iterations),'FontWeight','normal');
    if phase
        ylabel(tl,'Unwrapped phase (deg)','FontSize',12);name='phase';
        subtitle(tl,{'Curves unwrapped independently; gray points shifted only by integer 360 deg to the raw-model branch.', ...
            'No smoothing or complex-data modification; colored models are not phase-aligned to each other.'},'FontSize',10);
    else
        ylabel(tl,'Magnitude (dB, recorded output / motor-command unit)','FontSize',12);name='magnitude';
        subtitle(tl,'All 16 channels shown regardless of the LRM fit convergence status. Welch H1 is an estimate, not plant truth.','FontSize',10);
    end
    exportgraphics(fig,fullfile(out,['frf_preprocessing_matrix_' name '.png']),'Resolution',200);clear cleanup;
end

% Check the nonparametric entries themselves, without any fitted model.
assert(isfield(C,'welchTrain') && isfield(C,'lrmTrain'),'Train FRFs are required for the entry diagnostic.');
fig=figure('Visible','off','Color','w','Position',[60 60 1320 870]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    for col=1:2
        ch=cols(col);ax=nexttile(tl);hold(ax,'on');
        allH=[C.welchTrain(:,ch),C.welchValidation(:,ch),C.lrmTrain(:,ch),C.lrmValidation(:,ch)];
        if row==1,values=toDb(allH);else,values=unwrap(angle(allH),[],1)*180/pi;end
        p1=semilogx(ax,f,values(:,1),'--','Color',darkgray,'LineWidth',.95);
        p2=semilogx(ax,f,values(:,2),'-','Color',gray,'LineWidth',1.15);
        p3=semilogx(ax,f,values(:,3),'--','Color',teal,'LineWidth',.95);
        p4=semilogx(ax,f,values(:,4),'-','Color',purple,'LineWidth',1.1);
        style(ax);xlim(ax,[180 300]);xticks(ax,[180 200 220 260 300]);
        title(ax,sprintf('u_%d to y_%d',channelInputs(col),channelOutputs(col)),'FontWeight','normal');
        xlabel(ax,'Frequency (Hz)');
        if row==1,ylabel(ax,'FRF magnitude (dB)');else,ylabel(ax,'Independently unwrapped phase (deg)');end
        if row==1 && col==1,legendHandles=[p1 p2 p3 p4];end
    end
end
lg=legend(legendHandles,{'Welch H1: train','Welch H1: validation','LRM: train','LRM: validation'}, ...
    'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
title(tl,'Nonparametric FRF inputs near the antiresonance and first broad resonance','FontWeight','normal');
subtitle(tl,'First/second 30 s records. No model curves; phase unwrapped independently without branch alignment.','FontSize',10);
exportgraphics(fig,fullfile(out,'frf_preprocessing_nonparametric_detail.png'),'Resolution',200);clear cleanup;

checks=struct('frequency_count',numel(f),'channel_mapping','column = output + 4*(input-1)', ...
    'antiresonance_channels',[1 1;2 2],'antiresonance_window_hz',[180 260], ...
    'nonparametric_window_hz',[180 300],'matrix_model','lrm_inverse_magnitude', ...
    'matrix_model_converged',inverseLRM.fit.converged,'matrix_model_score',inverseLRM.fit.score, ...
    'phase_display','Independent unwrap, with gray points shifted by integer 360 deg to raw-SRIVC branch only', ...
    'routes_phase_display','Independent local unwrap; one whole-curve integer 360-degree offset set at first local frequency against gray reference', ...
    'smoothing','none','weight_definition','inverse_magnitude is the quadratic weight W, not its square');
disp(checks);
end

function model=selectModel(C,tag)
ix=find(strcmp({C.models.tag},tag));assert(numel(ix)==1,'Expected one comparison model tagged %s.',tag);
model=C.models(ix);
assert(isequal(size(model.H),size(C.welchValidation)),'Model/reference FRF grids must match.');
end

function out=yesNo(value)
if value,out='yes';else,out='no';end
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function style(ax)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on', ...
    'LineWidth',.65,'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
