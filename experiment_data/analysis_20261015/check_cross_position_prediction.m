function check_cross_position_prediction
% Fixed z=0 transfer functions applied to other positions, with the same
% record-specific boundary-term treatment as the local validation metric.
out=fileparts(mfilename('fullpath'));base=fileparts(out);
Q=load(fullfile(out,'raw_local_models.mat'));rows={};
idx=find(cellfun(@(v)v.position_mm==0,Q.allRecords),1);
reference=Q.allRecords{idx};
for ip=1:numel(Q.allRecords)
    rec=Q.allRecords{ip};
    S=load(fullfile(base,sprintf('auto_identify_motor1234_z%g.mat',rec.position_mm)));
    fs=rec.sample_rate_hz;start=round(fs)+1;
    n=floor((size(S.u_Motor1,1)-start+1)/2);U=zeros(n,4);
    for k=1:4
        uk=S.(sprintf('u_Motor%d',k));
        U(:,k)=uk(start+n:start+2*n-1,2);
    end
    UF=fft(U-mean(U),n)/n;f=(0:n-1)'*fs/n;
    for g=1:4
        a=reference.groups{g}.train;
        id=f>=a.window_hz(1)&f<=a.window_hz(2);
        ff=f(id);uu=UF(id,:);Y=rec.groups{g}.measured_validation;
        assert(size(Y,1)==numel(ff),'Frequency grid does not match.');
        basis=modalBasis(a,ff);Xg=complex(zeros(numel(ff),4*size(basis,2)));
        for k=1:4
            Xg(:,(k-1)*size(basis,2)+(1:size(basis,2)))=basis.*uu(:,k);
        end
        forced=Xg*a.transfer_coefficients;
        w=(ff/mean(ff)).^2;
        A=[real(basis.*w);imag(basis.*w)];B=[real((Y-forced).*w);imag((Y-forced).*w)];
        sc=vecnorm(A);c=(A./sc)\B;c=c./sc.';
        pred=forced+basis*c;
        err=norm((pred-Y).*w,'fro')/norm(Y.*w,'fro');
        own=rec.groups{g}.validation_error;
        rows(end+1,:)={rec.position_mm,g,100*own,100*err,err/own}; %#ok<AGROW>
    end
end
comparison=cell2table(rows,'VariableNames',{'z_mm','group','own_model_error_percent', ...
    'z0_model_error_percent','error_ratio'});
writetable(comparison,fullfile(out,'cross_position_prediction.csv'));disp(comparison);

% Direct validation output spectra, rather than a ratio-based FRF benchmark.
fig=figure('Visible','off','Color','w','Position',[50 50 1300 850]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for g=1:4
    a=reference.groups{g};f=a.f;Y=a.measured_validation;P=a.predicted_validation;
    ax=nexttile(tl);hold(ax,'on');
    % Pool squared DFT magnitude in 0.5 Hz bins for readability only.
    bins=floor((f-f(1))/.5);nb=max(bins)+1;ff=zeros(nb,1);power=zeros(nb,3);
    for j=1:nb
        k=bins==j-1;ff(j)=mean(f(k));
        power(j,:)=[sum(abs(Y(k,:)).^2,'all'),sum(abs(P(k,:)).^2,'all'), ...
            sum(abs(Y(k,:)-P(k,:)).^2,'all')];
    end
    plot(ax,ff,10*log10(max(power(:,1),realmin)),'Color',[.55 .55 .55],'LineWidth',1.5);
    plot(ax,ff,10*log10(max(power(:,2),realmin)),'Color',[.04 .30 .62],'LineWidth',1.0);
    plot(ax,ff,10*log10(max(power(:,3),realmin)),'Color',[.75 .24 .08],'LineWidth',.8);
    xlabel(ax,'Frequency (Hz)');ylabel(ax,'Pooled spectral power (dB, recorded units)');
    title(ax,sprintf('%g-%g Hz; conditional validation error %.2f%%',f(1),f(end),100*a.validation_error));
    grid(ax,'on');xlim(ax,[f(1) f(end)]);
    if g==1,legend(ax,'Measured outputs','Predicted outputs','Residual','Location','best');end
end
title(tl,'z = 0 mm, second 30 s: transfer parameters fixed; boundary terms refitted');
exportgraphics(fig,fullfile(out,'raw_z0_output_validation.png'),'Resolution',170);close(fig);
end

function B=modalBasis(a,f)
nm=numel(a.frequency_hz);z=exp(2i*pi*f/a.sample_rate_hz);
B=complex(zeros(numel(f),2*nm+4));
for j=1:nm
    lambda=2*pi*a.frequency_hz(j)*(-a.damping(j)+1i*sqrt(1-a.damping(j)^2));
    p=exp(lambda/a.sample_rate_hz);x=1./(z-p);y=1./(z-conj(p));
    B(:,2*j-1)=x+y;B(:,2*j)=1i*(x-y);
end
v=(f-a.center_hz)/a.scale_hz;
B(:,end-3:end)=[ones(size(v)),1i*ones(size(v)),v,1i*v];
end
