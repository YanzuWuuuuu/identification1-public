function report = make_calibration_target_figure()
% Show the actual identification targets, using an already selected estimate.
here=fileparts(mfilename('fullpath'));
out=fullfile(here,'figures','bode');
S=load(fullfile(here,'outputs','data_distinct.mat'),'p','nom');
Q=load(fullfile(here,'outputs_structured','distinct_N32','selected.mat'),'selected','metrics');
V=load(fullfile(out,'bode_source.mat'),'f');
f=V.f(:); s=2i*pi*f.'; truth=S.p; estimate=Q.selected.estimate;
relative=complex(zeros(3,numel(f),2,2)); % channel, frequency, side, truth/fit
nullVectors=zeros(4,3,2); angles=zeros(2,2);
for side=1:2
    if side==1
        [Ut,~,~]=svd(truth.M);[Un,~,~]=svd(S.nom.M);[Uh,~,~]=svd(estimate.M);
    else
        [~,~,Ut]=svd(truth.M);[~,~,Un]=svd(S.nom.M);[~,~,Uh]=svd(estimate.M);
    end
    nt=Ut(:,4); nn=Un(:,4); nh=Uh(:,4);
    if nn'*nt<0,nn=-nn;end
    if nh'*nt<0,nh=-nh;end
    nullVectors(:,:,side)=[nt nn nh];
    angles(side,:)=[asind(min(1,norm(nn*nn'-nt*nt',2))),asind(min(1,norm(nh*nh'-nt*nt',2)))];
    for version=1:2
        if version==1, p=truth;else,p=estimate;end
        d=reshape(1+p.a(:)*s+p.b(:)*s.^2,4,4,[]);
        for j=2:4
            if side==1
                delay=p.ty(j)-p.ty(1);den=reshape(d(j,1,:),1,[]);
            else
                delay=p.tu(j)-p.tu(1);den=reshape(d(1,j,:),1,[]);
            end
            relative(j-1,:,side,version)=exp(-delay*s).*reshape(d(1,1,:),1,[])./den;
        end
    end
end
err=max(abs(relative(:,:,:,2)./relative(:,:,:,1)-1),[],[1 2]);
assert(abs(err(1)-Q.metrics.relative_output_error)<1e-10);
assert(abs(err(2)-Q.metrics.relative_input_error)<1e-10);
assert(max(abs(angles(:,2)-[Q.metrics.output_null_angle_degrees;Q.metrics.input_null_angle_degrees]))<1e-10);
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 8.3 5.2]);
closer=onCleanup(@()close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','loose');
colors=[.23 .38 .50;.63 .40 .29;.40 .53 .44];
names={'Output','Input'};
for side=1:2
    ax=nexttile(layout,(side-1)*3+1);hold(ax,'on');
    n=nullVectors(:,:,side);
    h2=plot(ax,1:4,abs(n(:,2)-n(:,1)),'^:','Color',colors(2,:),'LineWidth',1,'MarkerSize',5);
    h3=plot(ax,1:4,abs(n(:,3)-n(:,1)),'s-','Color',colors(1,:),'LineWidth',1,'MarkerSize',5);
    set(ax,'YScale','log');ylim(ax,[1e-6 2e-2]);yticks(ax,[1e-6 1e-4 1e-2]);
    title(ax,[char('a'+(side-1)*3) ') ' names{side} ' null-direction error'],'FontWeight','normal');
    xlabel(ax,'Channel');ylabel(ax,'Absolute component error');xticks(ax,1:4);xlim(ax,[.75 4.25]);
    if side==1,legend(ax,[h2 h3],{'Nominal','Identified'},'Location','southwest','FontSize',7,'Box','off');end
    for quantity=1:2
        ax=nexttile(layout,(side-1)*3+quantity+1);hold(ax,'on');
        for j=1:3
            zt=relative(j,:,side,1);zh=relative(j,:,side,2);
            if quantity==1,yt=20*log10(abs(zt));yh=20*log10(abs(zh));else,yt=unwrap(angle(zt))*180/pi;yh=unwrap(angle(zh))*180/pi;end
            plot(ax,f,yt,'-','Color',colors(j,:),'LineWidth',1.4);
            plot(ax,f,yh,'--','Color',colors(j,:),'LineWidth',1);
        end
        xlabel(ax,'Frequency (Hz)');xlim(ax,[20 800]);xticks(ax,[20 400 800]);
        if quantity==1,ylabel(ax,'Relative gain (dB)');else,ylabel(ax,'Relative phase (deg)');end
        title(ax,[char('a'+(side-1)*3+quantity) ') ' names{side} ' relative dynamics'],'FontWeight','normal');
        if side==1 && quantity==1
            lh=gobjects(1,5);
            for ch=1:3,lh(ch)=plot(ax,NaN,NaN,'-','Color',colors(ch,:),'LineWidth',1.5);end
            lh(4)=plot(ax,NaN,NaN,'-','Color',[.3 .3 .3],'LineWidth',1.2);
            lh(5)=plot(ax,NaN,NaN,'--','Color',[.3 .3 .3],'LineWidth',1);
            lg=legend(ax,lh,{'Channel 2 / 1','Channel 3 / 1','Channel 4 / 1','True dynamics','Identified dynamics'}, ...
                'Orientation','horizontal','FontSize',7.2,'Box','off');
            lg.Layout.Tile='north';
        end
    end
end
axesList=findall(fig,'Type','axes');
set(axesList,'FontName','Times New Roman','FontSize',8,'LineWidth',.55,'TickDir','out','Box','on','XGrid','on','YGrid','on','GridAlpha',.15);
exportgraphics(fig,fullfile(out,'calibration_targets.png'),'Resolution',300,'Padding',6,'Units','points');
exportgraphics(fig,fullfile(out,'calibration_targets.pdf'),'ContentType','vector','Padding',6,'Units','points');
save(fullfile(out,'calibration_targets_source.mat'),'f','relative','nullVectors','angles','-v7');
report=struct('source_case','distinct_N32','source','Existing selected estimate and same simulated truth; no optimization rerun', ...
    'nominal_output_null_angle_deg',angles(1,1),'estimated_output_null_angle_deg',angles(1,2), ...
    'nominal_input_null_angle_deg',angles(2,1),'estimated_input_null_angle_deg',angles(2,2), ...
    'max_relative_output_error_percent',100*err(1),'max_relative_input_error_percent',100*err(2), ...
    'display_convention','Unit null vectors sign-aligned only for display; all three non-reference channel ratios shown; true solid and estimated dashed');
fid=fopen(fullfile(out,'calibration_targets_report.json'),'w','n','UTF-8');
assert(fid>=0);fileCloser=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));disp(report);
end
