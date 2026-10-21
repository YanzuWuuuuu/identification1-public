function audit_experimental_records(windowSeconds)
% Preliminary MIMO spectral audit. Original MAT files are read only.
if nargin<1,windowSeconds=2;end
base=fileparts(fileparts(mfilename('fullpath')));
out=fileparts(mfilename('fullpath'));
if windowSeconds~=2
    out=fullfile(out,sprintf('resolution_%gs',windowSeconds));
    if ~exist(out,'dir'),mkdir(out);end
end
files=dir(fullfile(base,'*.mat'));
yn={'pre_y_0__0_','pre_y_1__0_','pre_y_2__0_','pre_y_3__0_'};
un={'u_Motor1','u_Motor2','u_Motor3','u_Motor4'};
allResults=struct([]);
for a=1:numel(files)
    S=load(fullfile(base,files(a).name));
    t=S.(un{1})(:,1);dt=median(diff(t));fs=1/dt;
    U=zeros(numel(t),4);Y=U;timeError=0;
    for j=1:4
        U(:,j)=S.(un{j})(:,2);Y(:,j)=S.(yn{j})(:,2);
        timeError=max([timeError;abs(S.(un{j})(:,1)-t);abs(S.(yn{j})(:,1)-t)]);
    end
    assert(all(isfinite([t;U(:);Y(:)])) && all(diff(t)>0));
    assert(timeError<dt*1e-3 && max(abs(diff(t)-dt))<dt*1e-3);
    nwin=round(windowSeconds*fs);nfirst=round(fs)+1;nhalf=floor((numel(t)-nfirst+1)/2);
    itrain=nfirst:nfirst+nhalf-1;ival=nfirst+nhalf:nfirst+2*nhalf-1;
    train=spectralEstimate(U(itrain,:),Y(itrain,:),fs,nwin);
    validation=spectralEstimate(U(ival,:),Y(ival,:),fs,nwin);
    whole=spectralEstimate(U(nfirst:end,:),Y(nfirst:end,:),fs,nwin);
    f=whole.f;band=f>=25 & f<=1500;
    indicator=squeeze(sqrt(sum(abs(whole.H).^2,[1,2]))).*(2*pi*f).^2;
    db=20*log10(max(indicator,realmin));
    smoothDb=movmedian(db,3);
    [pks,loc,~,prom]=findpeaks(smoothDb(band),f(band),'MinPeakProminence',3,'MinPeakDistance',12);
    [~,ord]=sort(prom,'descend');ord=ord(1:min(14,numel(ord)));
    peaks=table(loc(ord),prom(ord),pks(ord),'VariableNames',{'frequency_hz','prominence_db','indicator_db'});
    peaks=sortrows(peaks,'frequency_hz');
    repeatedInput=norm(detrend(U(itrain,:),'constant')-detrend(U(ival,:),'constant'),'fro')/norm(detrend(U(itrain,:),'constant'),'fro');
    diffH=train.H(:,:,band)-validation.H(:,:,band);
    avgH=(train.H(:,:,band)+validation.H(:,:,band))/2;
    splitError=norm(diffH(:))/norm(avgH(:));
    R=struct('filename',files(a).name,'sample_rate_hz',fs,'sample_count',numel(t), ...
        'duration_s',t(end)-t(1),'start_time_s',t(1),'time_alignment_error_s',timeError, ...
        'input_mean',mean(U),'output_mean',mean(Y),'input_std',std(U),'output_std',std(Y), ...
        'split_input_difference',repeatedInput,'split_frf_relative_difference',splitError, ...
        'median_multiple_coherence',median(whole.coherence(:,band),2).', ...
        'median_input_condition',median(whole.condition(band)), ...
        'maximum_input_condition',max(whole.condition(band)), ...
        'peak_frequency_hz',peaks.frequency_hz.','peak_prominence_db',peaks.prominence_db.');
    if a==1, allResults=repmat(R,numel(files),1); end
    allResults(a)=R;
    stem=erase(files(a).name,'.mat');
    save(fullfile(out,[stem '_spectral.mat']),'train','validation','whole','R','peaks');
    writetable(peaks,fullfile(out,[stem '_peaks.csv']));
    fprintf('%s fs=%.3f N=%d median coherence=%s cond=%.3g splitFRF=%.3f%% repeatU=%.4g\n', ...
        R.filename,fs,numel(t),mat2str(R.median_multiple_coherence,3),R.median_input_condition,100*splitError,repeatedInput);
    disp(peaks);
end
save(fullfile(out,'spectral_audit.mat'),'allResults');
fid=fopen(fullfile(out,'spectral_audit.json'),'w');c=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(allResults,PrettyPrint=true),'char');clear c;
fig=figure('Visible','off','Color','w','Position',[60 60 1400 800]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact');
ax1=nexttile(tl);hold(ax1,'on');ax2=nexttile(tl);hold(ax2,'on');
for a=1:numel(files)
    Q=load(fullfile(out,[erase(files(a).name,'.mat') '_spectral.mat']));
    f=Q.whole.f;b=f>=25 & f<=1500;
    strength=squeeze(sqrt(sum(abs(Q.whole.H).^2,[1,2]))).*(2*pi*f).^2;
    plot(ax1,f(b),20*log10(strength(b)),'LineWidth',1.2,'DisplayName',erase(files(a).name,'auto_identify_motor1234_'));
    plot(ax2,f(b),median(Q.whole.coherence(:,b),1),'LineWidth',1.1);
end
set(ax1,'XScale','log');set(ax2,'XScale','log');
xlim(ax1,[25 1500]);xlim(ax2,[25 1500]);ylim(ax2,[0 1.02]);
grid(ax1,'on');grid(ax2,'on');legend(ax1,'Location','northwest','Interpreter','none');
ylabel(ax1,'20 log_{10}(omega^2 ||H||_F), raw units');
ylabel(ax2,'Median multiple coherence');xlabel(ax2,'Frequency (Hz)');
title(ax1,'Measured four-input/four-output responses: preliminary closed-loop spectral estimate');
exportgraphics(fig,fullfile(out,'spectral_overview.png'),'Resolution',180);close(fig);
end

function R=spectralEstimate(U,Y,fs,nwin)
step=floor(nwin/2);starts=1:step:size(U,1)-nwin+1;nseg=numel(starts);
nf=floor(nwin/2)+1;Suu=complex(zeros(4,4,nf));Syu=Suu;Syy=Suu;
w=hann(nwin,'periodic');
for b=starts
    uu=fft(detrend(U(b:b+nwin-1,:)).*w,nwin);
    yy=fft(detrend(Y(b:b+nwin-1,:)).*w,nwin);
    for j=1:4
        for k=1:4
            Suu(j,k,:)=Suu(j,k,:)+reshape(uu(1:nf,j).*conj(uu(1:nf,k)),1,1,nf);
            Syu(j,k,:)=Syu(j,k,:)+reshape(yy(1:nf,j).*conj(uu(1:nf,k)),1,1,nf);
            Syy(j,k,:)=Syy(j,k,:)+reshape(yy(1:nf,j).*conj(yy(1:nf,k)),1,1,nf);
        end
    end
end
fac=1/(nseg*fs*sum(w.^2));Suu=Suu*fac;Syu=Syu*fac;Syy=Syy*fac;
H=complex(zeros(4,4,nf));coherence=zeros(4,nf);condition=zeros(nf,1);
for k=1:nf
    A=(Suu(:,:,k)+Suu(:,:,k)')/2;
    scale=sqrt(max(real(diag(A)),realmin));C=A./(scale*scale');
    condition(k)=cond(C);
    if rcond(C)>1e-12
        H(:,:,k)=(Syu(:,:,k)./scale')/C./scale';
        explained=real(diag(H(:,:,k)*Syu(:,:,k)'));
        coherence(:,k)=min(1,max(0,explained./max(real(diag(Syy(:,:,k))),realmin)));
    else
        H(:,:,k)=NaN;coherence(:,k)=NaN;
    end
end
R=struct('f',(0:nf-1)'*fs/nwin,'H',H,'coherence',coherence,'condition',condition, ...
    'Suu',Suu,'Syu',Syu,'Syy',Syy,'segment_count',nseg,'window_samples',nwin, ...
    'window_seconds',nwin/fs,'frequency_resolution_hz',fs/nwin, ...
    'method','MIMO Welch H1; total closed-loop inputs; preliminary, not an unbiased plant estimator');
end
