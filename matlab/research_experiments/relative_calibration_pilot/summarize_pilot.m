function diagnostics=summarize_pilot()
% Diagnose saved pilot outcomes without rerunning or selecting fits.
folder=fileparts(mfilename('fullpath')); outdir=fullfile(folder,'outputs');
names={'common','repeated','distinct'}; diagnostics=struct([]);
for k=1:3
    dat=load(fullfile(outdir,[names{k} '_noise1.mat']));
    t=dat.truth; m=dat.metrics; J=dat.J; res=dat.res;
    diagout.case=m.case;
    diagout.rank3_chart_reciprocal_condition=rcond(reshape(dat.xh(1:9),3,3));
    diagout.active_bound_indices=find(abs(dat.xh-dat.lb)<1e-7 | abs(dat.xh-dat.ub)<1e-7).';
    diagout.max_normalized_residual_gradient=max(abs(J'*res)./(vecnorm(J).'*norm(res)));
    e=dat.estimate;
    diagout.output_compensator_global_max_gain=max(arrayfun(@(j) ratioPeak(e.a(j,1),e.b(j,1),e.a(1,1),e.b(1,1)),1:4));
    diagout.input_compensator_global_max_gain=max(arrayfun(@(j) ratioPeak(e.a(1,j),e.b(1,j),e.a(1,1),e.b(1,1)),1:4));
    ratio=zeros(size(dat.fv));
    for n=1:numel(ratio)
        s=2i*pi*dat.fv(n); F=zeros(4);
        for j=1:4, F=F+t.R(:,:,j)/(s-t.lambda(j))+conj(t.R(:,:,j))/(s-conj(t.lambda(j))); end
        hy=exp(-s*t.ty)./(1+s*t.filter_y); hu=exp(-s*t.tu)./(1+s*t.filter_u);
        ratio(n)=norm(diag(hy)*(t.M/s^2)*diag(hu),'fro')/norm(diag(hy)*F*diag(hu),'fro');
    end
    [diagout.max_unprojected_rigid_to_flexible_frobenius_ratio, imax]=max(ratio);
    diagout.frequency_of_max_unprojected_ratio_hz=dat.fv(imax);
    diagout.joint_vs_static_rigid_reduction_db=m.baseline.joint.rigid_reduction_vs_nominal_db-m.baseline.static_only.rigid_reduction_vs_nominal_db;
    if isempty(diagnostics), diagnostics=diagout; else, diagnostics(end+1)=diagout; end %#ok<AGROW>
    allMetrics(k)=m; %#ok<AGROW>
end

function peak=ratioPeak(an,bn,ad,bd)
% Exact rational magnitude maximum over all real frequencies; delays have
% unit magnitude. x=(omega/1e4)^2 conditions the stationary polynomial.
w=1e4; nn=[(bn*w^2)^2,(an*w)^2-2*bn*w^2,1];
dd=[(bd*w^2)^2,(ad*w)^2-2*bd*w^2,1];
rr=roots([nn(1)*dd(2)-nn(2)*dd(1),2*(nn(1)-dd(1)),nn(2)-dd(2)]);
xx=real(rr(abs(imag(rr))<1e-9 & real(rr)>0));
val=[1;nn(1)/dd(1);polyval(nn,xx)./polyval(dd,xx)];
peak=sqrt(max(val));
end
fid=fopen(fullfile(outdir,'diagnostics.json'),'w','n','UTF-8'); clean=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(diagnostics,'PrettyPrint',true));
fig=figure('Color','w','Position',[100 100 1300 880],'Visible','off');
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
nexttile; cv=allMetrics(3).curves;
loglog(cv.frequency_hz,cv.rigid_amplitude(:,1),'Color',[.65 .65 .65],'LineWidth',1.2); hold on;
loglog(cv.frequency_hz,cv.rigid_amplitude(:,2),'Color',[.9 .45 .1],'LineWidth',1.2);
loglog(cv.frequency_hz,cv.rigid_amplitude(:,3),'Color',[.1 .4 .8],'LineWidth',1.5);
loglog(cv.frequency_hz,cv.target_flexible_amplitude(:,3),'k','LineWidth',1.2);
grid on; xlabel('Frequency (Hz)'); ylabel('Projected FRF magnitude');
title('Distinct filters: actual truth after calibration');
legend('Nominal: rigid','Static only: rigid','Joint: rigid','Joint: target flexibility','Location','southwest');
nexttile; hold on;
cols=lines(3);
for k=1:3
    cv=allMetrics(k).curves;
    semilogy(cv.frequency_hz,100*max([cv.output_relative_error,cv.input_relative_error],[],2),'Color',cols(k,:),'LineWidth',1.5);
end
set(gca,'YScale','log'); grid on; xlabel('Frequency (Hz)'); ylabel('Worst relative dynamic error (%)');
title('Estimated relative channel response'); legend(names,'Location','southeast');
nexttile;
ret=zeros(4,3);
for k=1:3, for j=1:4, ret(j,k)=100*(allMetrics(k).mode(j).joint.retention_amplitude_ratio-1); end, end
bar(1:4,ret); grid on; xlabel('Flexible mode'); ylabel('Residue amplitude deviation (%)');
title('True poles: common factor removed'); legend(names,'Location','best');
nexttile;
rr=zeros(3,3);
for k=1:3
    mm=allMetrics(k).baseline;
    rr(k,:)=[mm.nominal.leak_to_flex_amplitude_ratio,mm.static_only.leak_to_flex_amplitude_ratio,mm.joint.leak_to_flex_amplitude_ratio];
end
bar(rr); set(gca,'YScale','log','XTickLabel',names); grid on;
ylabel('Rigid / target-flexible band L2 ratio');
title('2.07-799.8 Hz, same noisy realization'); legend('Nominal','Static only','Joint','Location','best');
title(tl,'Local feasibility pilot: known-model synthetic FRFs, near-truth initialization');
exportgraphics(fig,fullfile(outdir,'pilot_summary.png'),'Resolution',160);
savefig(fig,fullfile(outdir,'pilot_summary.fig')); close(fig);
disp(jsonencode(diagnostics,'PrettyPrint',true));
end
