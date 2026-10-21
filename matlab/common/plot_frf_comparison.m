function fig=plot_frf_comparison(f,G,Gfit,no,ni,titleText)
fig=figure('Name',titleText,'Color','w','Position',[100,100,1400,900]);
layout=tiledlayout(no,ni,'TileSpacing','compact','Padding','compact');
for output=1:no
    for input=1:ni
        nexttile; k=(input-1)*no+output;
        semilogx(f,20*log10(max(abs(G(k,:)),realmin)),'Color',[.15,.6,.25]); hold on;
        semilogx(f,20*log10(max(abs(Gfit(k,:)),realmin)),'b');
        semilogx(f,20*log10(max(abs(G(k,:)-Gfit(k,:)),realmin)),'r');
        grid on; title(sprintf('y%d / u%d',output,input));
        if output==no, xlabel('Frequency (Hz)'); end
        if input==1, ylabel('Magnitude (dB)'); end
    end
end
title(layout,titleText);
legend({'LRM data','Estimated model','Complex residual'},'Location','best');
end
