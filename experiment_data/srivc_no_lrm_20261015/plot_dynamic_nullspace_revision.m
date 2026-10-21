function files=plot_dynamic_nullspace_revision
% Scientific diagnostic plots, including the failed unguarded alternative.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'dynamic_nullspace_20261016');
S=load(fullfile(out,'dynamic_evaluation.mat'),'evaluation');e=S.evaluation;T=e.metrics;
pos=[-.5 0 .5 1];codes={'D0','p40_l8','p40_l8_g'};
names={'Previous band-static','Dynamic, no cross-band guards','Dynamic, with cross-band guards'};
colors=[.48 .22 .55;.12 .48 .32;.04 .27 .53];blocks={'logical_to_N','N_to_logical'};
blockNames={'Logical input to nullspace output','Nullspace input to logical outputs'};
files={};
fig=figure('Visible','off','Color','w','Position',[100 100 1440 1080]);cleanup=onCleanup(@()close(fig));
tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');bands=[3 25;25 80;100 170];
for row=1:3,for col=1:2
    ax=nexttile(tl);values=zeros(4,3);
    for p=1:4,for m=1:3
        ix=T.position_mm==pos(p) & strcmp(T.scope,'own_position') & strcmp(T.method,codes{m}) & strcmp(T.block,blocks{col}) & T.fmin==bands(row,1) & T.fmax==bands(row,2);
        assert(nnz(ix)==1);values(p,m)=T.over_B0(ix);
    end,end
    h=bar(ax,pos,values,'grouped');hold(ax,'on');
    for m=1:3,h(m).FaceColor=colors(m,:);h(m).EdgeColor='none';end
    yline(ax,1,':','Nominal = 1','HandleVisibility','off','LabelHorizontalAlignment','left');
    style(ax);xticks(ax,pos);xlabel(ax,'Position z (mm)');ylabel(ax,'Response RMS / nominal');
    title(ax,sprintf('%s | %d-%d Hz',blockNames{col},bands(row,1),bands(row,2)),'FontWeight','normal');
    ylim(ax,[0 1.1*max([values(:);1])]);if row==1 && col==1,hl=h;end
end,end
lg=legend(hl,names,'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
title(tl,'Dynamic correction: in-band gain does not imply full-band decoupling','FontWeight','normal');
subtitle(tl,'Second-half LRM, previously inspected; all fits use first-half data. Lower is better. No uncertainty bars are available.','FontSize',10);
files=[files,exportPair(fig,out,'dynamic_revision_crossband_tradeoff')];clear cleanup;

f=e.responses{1}.f;B0=findResponse(e,0,'B0');D0=findResponse(e,0,'D0');NEW=findResponse(e,0,'p40_l8_g');
codeNames={'Rx ref.','Ry ref.','Z ref.','N'};Zs={B0,D0,NEW};cc=[.65 .65 .65;colors(1,:);colors(3,:)];
for phasePlot=0:1
    fig=figure('Visible','off','Color','w','Position',[100 100 1450 1120]);cleanup=onCleanup(@()close(fig));
    tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for iy=1:4,for iu=1:4
        ax=nexttile(tl);hold(ax,'on');h=gobjects(1,3);
        for m=1:3
            v=reshape(Zs{m}(iy,iu,:),[],1);
            if phasePlot,y=unwrap(angle(v))*180/pi;else,y=20*log10(max(abs(v),realmin));end
            h(m)=semilogx(ax,f,y,'Color',cc(m,:),'LineWidth',1.05);
        end
        style(ax);set(ax,'XScale','log');xlim(ax,[3 800]);xticks(ax,[3 10 30 100 300 800]);
        title(ax,[codeNames{iu} ' to ' codeNames{iy}],'FontWeight','normal');
        if iy==4,xlabel(ax,'Frequency (Hz)');end
        if iu==1
            if phasePlot,ylabel(ax,'Unwrapped phase (deg)');else,ylabel(ax,'Magnitude (dB re recorded units)');end
        end
        if iy==1 && iu==1,hl=h;end
    end,end
    lg=legend(hl,{'Nominal static','Previous band-static','Guarded dynamic'},'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
    title(tl,'z=0: complete transformed measured response','FontWeight','normal');
    subtitle(tl,'Second-half LRM, not an isolated rigid/flexible response. Logical 3-by-3 block is unchanged by the new algorithm.','FontSize',10);
    if phasePlot,name='dynamic_revision_z0_phase';else,name='dynamic_revision_z0_magnitude';end
    files=[files,exportPair(fig,out,name)];clear cleanup;
end
end
function Z=findResponse(e,pos,method)
ix=cellfun(@(r)r.position_mm==pos && strcmp(r.method,method) && strcmp(r.scope,'own_position'),e.responses);
assert(nnz(ix)==1);r=e.responses{ix};Z=r.Z;
end
function style(ax)
set(ax,'FontName','Arial','FontSize',10,'Box','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',.15,'MinorGridAlpha',.08);grid(ax,'on');
end
function files=exportPair(fig,out,name)
png=fullfile(out,[name '.png']);pdf=fullfile(out,[name '.pdf']);
exportgraphics(fig,png,'Resolution',180);exportgraphics(fig,pdf,'ContentType','vector');files={png,pdf};
end
