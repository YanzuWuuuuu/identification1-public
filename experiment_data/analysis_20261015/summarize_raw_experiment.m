function summarize_raw_experiment
out=fileparts(mfilename('fullpath'));base=fileparts(out);
Q=load(fullfile(out,'raw_local_models.mat'));
T=Q.modalTable;G=Q.groupTable;
main=T((T.group==1)|(T.group==2&T.component==2)|(T.group==3)|(T.group==4&T.component==2),:);
writetable(main,fullfile(out,'main_modal_frequencies.csv'));
positions=unique(main.z_mm).';
fig=figure('Visible','off','Color','w','Position',[50 50 1300 850]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
names={'A: broad resonance near 280 Hz','B: resonance near 360 Hz', ...
    'C: resonance near 493 Hz','D: resonance near 705 Hz'};
for m=1:4
    ax=nexttile(tl);id=m:4:height(main);
    plot(ax,main.z_mm(id),main.frequency_hz(id),'o-','LineWidth',1.6,'MarkerSize',7);hold(ax,'on');
    plot(ax,main.z_mm(id),main.second_half_frequency_hz(id),'s--','LineWidth',1.2,'MarkerSize',6);
    title(ax,names{m});xlabel(ax,'z position (mm)');ylabel(ax,'Fitted frequency (Hz)');
    xticks(ax,positions);grid(ax,'on');
    if m==1,legend(ax,'First 30 s','Second 30 s','Location','best');end
end
title(tl,'Direct fits to raw U/Y spectra; z = -1 mm excluded due to an abnormal motor command');
exportgraphics(fig,fullfile(out,'raw_position_comparison.png'),'Resolution',170);close(fig);

% Four-by-four transfer functions from the direct local output-error models.
idx=find(cellfun(@(v)v.position_mm==0,Q.allRecords),1);rec=Q.allRecords{idx};
D=load(fullfile(out,'auto_identify_motor1234_z0_spectral.mat'));
for phase=[false,true]
    fig=figure('Visible','off','Color','w','Position',[40 40 1450 1000]);
    tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ax=nexttile(tl);hold(ax,'on');ib=D.validation.f>=25 & D.validation.f<=800;
            f=D.validation.f(ib);h=squeeze(D.validation.H(y,u,ib));
            if phase,v=unwrap(angle(h))*180/pi;else,v=20*log10(abs(h));end
            semilogx(ax,f,v,'Color',[.67 .67 .67],'LineWidth',.65);
            for g=1:4
                model=rec.groups{g}.train;ff=linspace(model.window_hz(1),model.window_hz(2),650).';
                HH=localTransfer(model,ff);hh=squeeze(HH(y,u,:));
                if phase
                    vv=unwrap(angle(hh))*180/pi;vv=vv+360*round((interp1(f,v,ff(1))-vv(1))/360);
                else
                    vv=20*log10(abs(hh));
                end
                semilogx(ax,ff,vv,'Color',[.04 .30 .62],'LineWidth',1.3);
            end
            set(ax,'XScale','log');xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);grid(ax,'on');
            title(ax,sprintf('u_%d to y_%d',u,y));
            if y==4,xlabel(ax,'Frequency (Hz)');end
            if u==1
                if phase,ylabel(ax,'Unwrapped phase (deg)');else,ylabel(ax,'Magnitude (dB, recorded units)');end
            end
        end
    end
    title(tl,'z = 0 mm: raw-data local model (blue); preliminary Welch estimate (gray)');
    if phase,name='raw_z0_bode_phase.png';else,name='raw_z0_bode_magnitude.png';end
    exportgraphics(fig,fullfile(out,name),'Resolution',160);close(fig);
end

% Quantify the abnormal command and verify the sensor-unit inference.
fileList=dir(fullfile(base,'*.mat'));anomaly=struct;diagRows={};
fig=figure('Visible','off','Color','w','Position',[60 60 1300 720]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact');ax1=nexttile(tl);ax2=nexttile(tl);hold(ax2,'on');
for k=1:numel(fileList)
    S=load(fullfile(base,fileList(k).name));t=S.u_Motor3(:,1);t=t-t(1);u=S.u_Motor3(:,2);
    z=str2double(extractBetween(fileList(k).name,'_z','.mat'));p=polyfit(t,u,1);
    ymean=zeros(1,4);for j=1:4,ymean(j)=mean(S.(sprintf('pre_y_%d__0_',j-1))(:,2));end
    diagRows(end+1,:)={z,u(1),u(end),mean(u),std(u),p(1),mean(S.d_z(:,2)),ymean(1),ymean(2),ymean(3),ymean(4)}; %#ok<AGROW>
    take=1:50:numel(t);
    if z==-1
        plot(ax1,t(take),u(take),'Color',[.7 .17 .08],'LineWidth',1.2);
        anomaly=struct('file',fileList(k).name,'motor3_start',u(1),'motor3_end',u(end), ...
            'motor3_mean',mean(u),'motor3_slope_per_second',p(1),'linear_detrended_std',std(u-polyval(p,t)));
    else
        plot(ax2,t(take),u(take),'LineWidth',.75,'DisplayName',sprintf('z=%g mm',z));
    end
end
grid(ax1,'on');grid(ax2,'on');xlabel(ax1,'Elapsed time (s)');xlabel(ax2,'Elapsed time (s)');
ylabel(ax1,'Motor 3 command, recorded units');ylabel(ax2,'Motor 3 command, recorded units');
title(ax1,'z = -1 mm: large monotonic command drift');title(ax2,'Other four positions: comparison of the same channel');
legend(ax2,'Location','eastoutside');
exportgraphics(fig,fullfile(out,'motor3_anomaly.png'),'Resolution',170);close(fig);
diagnostics=sortrows(cell2table(diagRows,'VariableNames',{'z_mm','motor3_start','motor3_end','motor3_mean','motor3_std', ...
    'motor3_slope_per_s','d_z_mean','sensor1_mean','sensor2_mean','sensor3_mean','sensor4_mean'}),'z_mm');
writetable(diagnostics,fullfile(out,'record_diagnostics.csv'));
summary=struct('method','Local shared-pole output-error fits to raw U/Y DFT; closed-loop total input', ...
    'frequency_audit_band_hz',[25 1500],'local_fit_windows_hz',[240 315;315 385;455 535;650 735], ...
    'sensor_unit_inference','Known 0.5 mm position steps and sensor mean increments are consistent with sensors in metres; motor command unit remains unverified', ...
    'excluded_position_mm',-1,'anomaly',anomaly, ...
    'main_modes',table2struct(main),'all_local_modes',table2struct(T),'validation',table2struct(G), ...
    'validation_scope','Mechanical/transfer coefficients frozen; record-specific boundary terms fitted to the second half. Not an independent experimental run.', ...
    'failed_resolution_check','4 s Welch spectral inversion was ill-conditioned and its modal fits are rejected.', ...
    'limitations','Preliminary local models; no controller-based closed-loop consistency claim; no relative-delay or nullspace calibration result.');
fid=fopen(fullfile(out,'identification_summary.json'),'w');clean=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(summary,PrettyPrint=true),'char');clear clean;
disp(main);disp(diagnostics);
end

function H=localTransfer(model,f)
nm=numel(model.frequency_hz);z=exp(2i*pi*f/model.sample_rate_hz);
B=complex(zeros(numel(f),2*nm+4));
for m=1:nm
    lambda=2*pi*model.frequency_hz(m)*(-model.damping(m)+1i*sqrt(1-model.damping(m)^2));
    p=exp(lambda/model.sample_rate_hz);a=1./(z-p);b=1./(z-conj(p));
    B(:,2*m-1)=a+b;B(:,2*m)=1i*(a-b);
end
v=(f-model.center_hz)/model.scale_hz;
B(:,end-3:end)=[ones(size(v)),1i*ones(size(v)),v,1i*v];
H=complex(zeros(4,4,numel(f)));
for u=1:4
    a=B*model.transfer_coefficients((u-1)*size(B,2)+(1:size(B,2)),:);
    H(:,u,:)=reshape(a.',4,1,[]);
end
end
