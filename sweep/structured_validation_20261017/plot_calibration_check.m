function plot_calibration_check(tag)
if nargin<1,tag='z0';end
out=fileparts(mfilename('fullpath'));s=load(fullfile(out,['evaluation_' tag '.mat']));e=s.evaluation;
f=e.f;Hs=e.calibrated_check_FRFs;colors=[.55 .55 .55;.1 .6 .5;.05 .28 .65;.55 .3 .7];
styles={'-','--','-','-.'};labels={'Nominal geometry','Structured, static only','Structured + relative delays','Direct band-static baseline'};
fig=figure('Visible','off','Color','w','Position',[100 100 1450 850]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
titles={'N input to logical outputs','Logical inputs to N output','Nullspace input to nullspace output'};
cols={13:15,[4 8 12],16};mag=cell(numel(Hs),3);
for block=1:3
    ax=nexttile(tl,block);hold(ax,'on');
    for c=1:numel(Hs)
        mag{c,block}=sqrt(sum(abs(Hs{c}(:,cols{block})).^2,2));
        semilogx(ax,f,20*log10(max(mag{c,block},1e-30)),styles{c},'Color',colors(c,:),'LineWidth',1.15);
    end
    set(ax,'XScale','log');grid(ax,'on');xlim(ax,[3 800]);title(ax,titles{block});
    xlabel(ax,'Frequency (Hz)');ylabel(ax,'Magnitude (dB, unit-norm ports)');
    xline(ax,25,':','Color',[.7 .7 .7]);xline(ax,80,':','Color',[.7 .7 .7]);
end
for block=1:2
    ax=nexttile(tl,3+block);hold(ax,'on');
    for c=2:numel(Hs)
        semilogx(ax,f,20*log10(max(mag{c,block},1e-30)./max(mag{1,block},1e-30)),styles{c},'Color',colors(c,:),'LineWidth',1.15);
    end
    set(ax,'XScale','log');grid(ax,'on');xlim(ax,[3 800]);ylim(ax,[-25 25]);yline(ax,0,':k');
    title(ax,[titles{block} ': change']);xlabel(ax,'Frequency (Hz)');ylabel(ax,'Change from nominal (dB; lower is less cross-response)');
end
ax=nexttile(tl,6);hold(ax,'on');
for c=1:numel(Hs),plot(ax,f,20*log10(max(mag{c,3},1e-30)),styles{c},'Color',colors(c,:),'LineWidth',1.3);end
grid(ax,'on');xlim(ax,[250 315]);title(ax,'N to N near first plate resonance');xlabel(ax,'Frequency (Hz)');ylabel(ax,'Magnitude (dB)');
lg=legend(tl.Children(end),labels,'Orientation','horizontal');lg.Layout.Tile='north';
title(tl,{['Offline calibration check | ' tag], ...
    'All curves transform the SAME held-out development FRF; calibrators use training data only.', ...
    'Unit-norm directions; pure causal delay alignment has unit magnitude. This is not a hardware closed-loop test.'},'FontSize',12);
exportgraphics(fig,fullfile(out,[tag '_calibration_check.png']),'Resolution',170);
exportgraphics(fig,fullfile(out,[tag '_calibration_check.pdf']),'ContentType','vector');close(fig);
end
