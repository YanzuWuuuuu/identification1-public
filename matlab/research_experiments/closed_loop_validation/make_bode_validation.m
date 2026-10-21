function report = make_bode_validation()
% Plot existing held-out records; this function never reruns identification.
here = fileparts(mfilename('fullpath'));
priorPath = path;
restorePath = onCleanup(@()path(priorPath)); %#ok<NASGU>
addpath(here);
out = fullfile(here,'figures','bode');
if ~isfolder(out), mkdir(out); end
S = load(fullfile(here,'outputs','data_distinct.mat'),'p','K','cfg');
Q = load(fullfile(here,'outputs_structured','distinct_N32','selected.mat'),'selected','metrics');
V = load(fullfile(here,'outputs_structured','test_distinct.mat'),'testDistinct');
truth = S.p; estimate = Q.selected.estimate; K = S.K; Ts = S.cfg.Ts;
validation = V.testDistinct;
f = validation.f(:); nf = numel(f); nr = validation.groups;
assert(nr == 8 && nf == 195);
fDense = logspace(log10(20),log10(800),1400).';
for j = 1:numel(truth.lambda)
    center = abs(truth.lambda(j))/(2*pi);
    fDense = [fDense;linspace(.94*center,1.06*center,220).']; %#ok<AGROW>
end
fDense = unique([fDense;f]);
fDense = fDense(fDense >= 20 & fDense <= 800);
[Gtrue,Ttrue] = responses(truth,f,K,Ts);
[Gfit,Tfit] = responses(estimate,f,K,Ts);
[GtrueDense,TtrueDense] = responses(truth,fDense,K,Ts);
[GfitDense,TfitDense] = responses(estimate,fDense,K,Ts);
Trecords = complex(zeros(16,nf,nr));
conditionD = zeros(nf,1);
for k = 1:nf
    D = validation.D(:,:,k);
    conditionD(k) = cond(D);
    for r = 1:nr
        % D is deterministic external excitation, never a noisy input ratio.
        Trecords(:,k,r) = reshape(reshape(validation.Y(:,k,r),4,4)/D,16,1);
    end
end
Tmeasured = mean(Trecords,3);
noiseSEM = sqrt(sum(abs(Trecords-Tmeasured).^2,3)/(nr*(nr-1)));
measuredResidual = abs(Tmeasured-Tfit);
truthResidual = abs(Ttrue-Tfit);
% Unwrap each channel independently along increasing frequency. Neither the
% measured data nor its branch choices use the model as a phase reference.
phaseMeasuredDeg = unwrap(angle(Tmeasured),[],2)*180/pi;
phaseModelDeg = unwrap(angle(TfitDense),[],2)*180/pi;
phaseModelAtDataDeg = interp1(fDense,phaseModelDeg.',f).';
phaseBranchDifference = round((phaseMeasuredDeg-phaseModelAtDataDeg)/360);
phaseComplexInvariantError = max(abs(exp(1i*phaseMeasuredDeg*pi/180) ...
    -Tmeasured./abs(Tmeasured)),[],'all');
assert(phaseComplexInvariantError < 1e-12);
channelRelativeError = sqrt(sum(abs(Tfit-Ttrue).^2,2)./sum(abs(Ttrue).^2,2));
truthError = norm(Tfit-Ttrue,'fro')/norm(Ttrue,'fro');
assert(abs(truthError-Q.metrics.heldout_closed_relative_error) < 1e-10);
assert(max(conditionD) < 1.000001);

% Three durable figures, all in the actual 4-output by 4-input arrangement.
plot_matrix('closed_loop_magnitude','magnitude',fDense,TfitDense,[], ...
    f,Tmeasured,measuredResidual,noiseSEM,out);
plot_matrix('closed_loop_phase','phase',fDense,TfitDense,[], ...
    f,Tmeasured,[],[],out);
plot_matrix('sampled_plant_magnitude','plant',fDense,GfitDense,GtrueDense, ...
    f,[],abs(Gfit-Gtrue),[],out);

% Reproducible numeric export, with measured and modeled responses separated.
[outputIndex,inputIndex,frequencyIndex] = ndgrid(1:4,1:4,1:nf);
repeatedFrequency = reshape(repmat(f.',16,1),[],1);
source = table(outputIndex(:),inputIndex(:),repeatedFrequency, ...
    real(Tmeasured(:)),imag(Tmeasured(:)),real(Tfit(:)),imag(Tfit(:)), ...
    real(Ttrue(:)),imag(Ttrue(:)),noiseSEM(:), ...
    real(Gfit(:)),imag(Gfit(:)),real(Gtrue(:)),imag(Gtrue(:)), ...
    'VariableNames',{'output','input','frequency_hz', ...
    'T_measured_real','T_measured_imag','T_estimated_real','T_estimated_imag', ...
    'T_true_real','T_true_imag','T_mean_noise_complex_sd', ...
    'G_estimated_real','G_estimated_imag','G_true_real','G_true_imag'});
assert(all(repeatedFrequency == f(frequencyIndex(:))));
writetable(source,fullfile(out,'bode_source.csv'));
save(fullfile(out,'bode_source.mat'),'f','fDense','Trecords','Tmeasured', ...
    'Tfit','Ttrue','TfitDense','TtrueDense','Gfit','Gtrue','GfitDense', ...
    'GtrueDense','noiseSEM','channelRelativeError','phaseMeasuredDeg', ...
    'phaseModelDeg','phaseModelAtDataDeg','phaseBranchDifference','-v7');
report = struct('case','distinct_N32','training_groups',32,'validation_groups',nr, ...
    'frequency_hz',[min(f) max(f)],'validation_frequency_count',nf, ...
    'dense_frequency_count',numel(fDense),'closed_truth_relative_error',truthError, ...
    'closed_measured_relative_error',norm(Tfit-Tmeasured,'fro')/norm(Tmeasured,'fro'), ...
    'nonparametric_noise_relative_error',norm(Tmeasured-Ttrue,'fro')/norm(Ttrue,'fro'), ...
    'per_channel_truth_relative_error',reshape(channelRelativeError,4,4), ...
    'plant_truth_relative_error',norm(Gfit-Gtrue,'fro')/norm(Gtrue,'fro'), ...
    'max_excitation_matrix_condition',max(conditionD), ...
    'phase_definition','Unwrapped independently along increasing frequency for each measured and modeled channel, starting at its principal phase. No model-guided branch shifts, smoothing, or point removal.', ...
    'phase_complex_invariant_error',phaseComplexInvariantError, ...
    'phase_max_dense_model_step_deg',max(abs(diff(phaseModelDeg,1,2)),[],'all'), ...
    'phase_branch_disagreement_count',nnz(phaseBranchDifference), ...
    'noise_definition','Complex standard deviation of the mean, estimated from 8 validation repetitions; not a 95% confidence interval.', ...
    'measurement_definition','Mean(Y_n / D), where D is known deterministic external excitation; not Y/U and not LRM.', ...
    'scope','Existing synthetic held-out records, new plots only, no re-estimation or hardware data. Plant G includes channel dynamics and exact hold/sampling.');
fid = fopen(fullfile(out,'bode_report.json'),'w','n','UTF-8');
assert(fid>=0); closer = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
disp(report);
end

function [G,T] = responses(p,f,K,Ts)
G = clv_parts(p,f,Ts);
Kf = freqresp(K,2*pi*f);
T = complex(zeros(size(G)));
for k = 1:numel(f)
    Gk = reshape(G(:,k),4,4);
    T(:,k) = reshape((eye(4)+Gk*Kf(:,:,k))\Gk,16,1);
end
end

function plot_matrix(name,kind,fd,model,truth,f,measured,residual,noise,out)
% A custom MIMO matrix is more appropriate than a generic line template.
blue = [0.23 0.38 0.50];
gray = [0.54 0.54 0.54];
rust = [0.65 0.36 0.24];
noiseColor = [0.48 0.55 0.47];
fig = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[1 1 7.4 6.7],'Renderer','painters');
closer = onCleanup(@()close(fig)); %#ok<NASGU>
layout = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
handles = gobjects(0); labels = {};
isPhase = strcmp(kind,'phase'); isPlant = strcmp(kind,'plant');
if isPhase
    measuredPhase = unwrap(angle(measured),[],2)*180/pi;
    modelPhase = unwrap(angle(model),[],2)*180/pi;
    phaseValues = [measuredPhase(:);modelPhase(:)];
    phaseLimits = 90*[floor(min(phaseValues)/90),ceil(max(phaseValues)/90)]+[-15 15];
    phaseTicks = (180*ceil(phaseLimits(1)/180)):180:(180*floor(phaseLimits(2)/180));
else
    if isPlant
        allValues = [db(model(:));db(truth(:));db(residual(:))];
    else
        allValues = [db(model(:));db(measured(:));db(residual(:));db(noise(:))];
    end
    limits = [20*floor(min(allValues)/20),20*ceil(max(allValues)/20)];
end
for j = 1:4
    for k = 1:4
        e = j+4*(k-1);
        ax = nexttile(layout,(j-1)*4+k); hold(ax,'on');
        set(ax,'XScale','log','FontName','Times New Roman','FontSize',8.3, ...
            'LineWidth',.55,'TickDir','out','Box','on', ...
            'XColor',[.25 .25 .25],'YColor',[.25 .25 .25], ...
            'GridColor',[.82 .82 .82],'GridAlpha',.5, ...
            'XMinorGrid','on','YMinorGrid','off','MinorGridAlpha',.18);
        if isPhase
            hm = plot(ax,f,measuredPhase(e,:),'.','Color',gray,'MarkerSize',3.4);
            hp = plot(ax,fd,modelPhase(e,:),'-','Color',blue,'LineWidth',1.05);
            ylim(ax,phaseLimits); yticks(ax,phaseTicks);
            if j==1 && k==1, handles=[hm hp]; labels={'Held-out data (8-record mean)','Identified model'}; end
        elseif isPlant
            ht = plot(ax,fd,db(truth(e,:)),'-','Color',[.48 .48 .48],'LineWidth',1.45);
            hp = plot(ax,fd,db(model(e,:)),'--','Color',blue,'LineWidth',.95);
            hr = plot(ax,f,db(residual(e,:)),'-','Color',rust,'LineWidth',.7);
            ylim(ax,limits); yticks(ax,limits(1):40:limits(2));
            if j==1 && k==1, handles=[ht hp hr]; labels={'True sampled plant','Identified sampled plant','Complex model error'}; end
        else
            hm = plot(ax,f,db(measured(e,:)),'.','Color',gray,'MarkerSize',3.4);
            hp = plot(ax,fd,db(model(e,:)),'-','Color',blue,'LineWidth',1.05);
            hr = plot(ax,f,db(residual(e,:)),'-','Color',rust,'LineWidth',.6);
            hn = plot(ax,f,db(noise(e,:)),':','Color',noiseColor,'LineWidth',.8);
            ylim(ax,limits); yticks(ax,limits(1):40:limits(2));
            if j==1 && k==1, handles=[hm hp hr hn]; labels={'Held-out data (8-record mean)','Identified model','Data-model residual','Mean noise: complex SD'}; end
        end
        xlim(ax,[20 800]); xticks(ax,[20 100 500]); grid(ax,'on');
        if j<4, ax.XTickLabel=[]; end
        if k>1, ax.YTickLabel=[]; end
        if isPlant, inputSymbol='u'; else, inputSymbol='d'; end
        channelTitle = ['$' inputSymbol '_' num2str(k) '\rightarrow y_' num2str(j) '$'];
        title(ax,channelTitle, ...
            'Interpreter','latex','FontSize',9,'FontWeight','normal');
    end
end
xlabel(layout,'Frequency (Hz)','FontName','Times New Roman','FontSize',10);
if isPhase
    ylabel(layout,'Phase (deg), unwrapped','FontName','Times New Roman','FontSize',10);
else
    ylabel(layout,'Magnitude (dB re 1 m/N)','FontName','Times New Roman','FontSize',10);
end
lgd = legend(handles,labels,'Orientation','horizontal','Box','off','FontName','Times New Roman','FontSize',8.3);
lgd.Layout.Tile='north';
if numel(labels)==4, lgd.NumColumns=2; else, lgd.NumColumns=numel(labels); end
drawnow;
exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',300,'Padding',6,'Units','points');
exportgraphics(fig,fullfile(out,[name '.pdf']),'ContentType','vector','Padding',6,'Units','points');
end

function y=db(z)
y=20*log10(max(abs(z),realmin));
end

