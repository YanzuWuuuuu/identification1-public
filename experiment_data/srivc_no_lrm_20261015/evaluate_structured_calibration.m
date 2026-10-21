function evaluation=evaluate_structured_calibration(tag)
% Offline evaluation of saved calibration directions and relative delays.
% B0: nominal static. B1: joint-fit static directions. B4: SAME B1 directions
% plus additional nonnegative, unit-magnitude frequency-domain delays.
% This function never estimates parameters or selects a model on validation.
if nargin<1 || isempty(tag),tag='z0';end
tag=char(string(tag));assert(~isempty(regexp(tag,'^[A-Za-z0-9_-]+$','once')),'Invalid result tag.');
base=fileparts(mfilename('fullpath'));out=fullfile(base,'structured_calibration_20261016');
source=fullfile(out,['calibration_' tag '.mat']);
if ~isfile(source),error('evaluate_structured_calibration:NotReady','Saved calibration is not ready: %s',source);end
S=load(source,'calib');calib=S.calib;f=calib.f(:);n=numel(f);
assert(n>=2 && all(isfinite(f)) && all(diff(f)>0),'Expected a finite increasing frequency grid.');
assert(isequal(size(calib.H_training),[n 16]) && isequal(size(calib.H_validation),[n 16]), ...
    'FRFs must be Nfrequency-by-16, output index varying fastest.');
assert(all(isfinite(calib.H_training(:))) && all(isfinite(calib.H_validation(:))),'Nonfinite experimental FRFs.');
joint=calib.fits.joint;fixed=calib.fits.fixed;
assert(isequal(size(joint.H),[n 16]) && isequal(size(fixed.H),[n 16]),'Saved fitted FRF size mismatch.');
assert(all(isfinite(joint.H(:))) && all(isfinite(fixed.H(:))),'Nonfinite fitted responses.');
assert(isequal(size(calib.nominal.Ty),[3 4]) && isequal(size(calib.nominal.Tu),[4 3]), ...
    'Nominal logical-axis maps must be 3-by-4 and 4-by-3.');

% Freeze nominal axis directions and recorded-command input budgets for all
% three methods. Each output row/input column has Euclidean norm one.
Pnom=calib.nominal.Ty;Qnom=calib.nominal.Tu;
py=vecnorm(Pnom,2,2);qu=vecnorm(Qnom,2,1);
assert(all(py>0) && all(qu>0) && all(isfinite([py(:);qu(:)])),'Zero/nonfinite nominal axis direction.');
Pnom=Pnom./py;Qnom=Qnom./qu;
ny0=unit(calib.nominal.ny);nu0=unit(calib.nominal.nu);
ny1=alignedUnit(joint.Ny,ny0);nu1=alignedUnit(joint.Nu,nu0);
nyFixed=alignedUnit(fixed.Ny,ny0);nuFixed=alignedUnit(fixed.Nu,nu0);
assert(numel(ny0)==4 && numel(nu0)==4 && numel(ny1)==4 && numel(nu1)==4,'Null vectors must have four entries.');

% The optimizer coordinate ordering is deliberately not used here.
ell=joint.delay_seconds;assert(isequal(size(ell),[4 4]) && isreal(ell) && all(isfinite(ell(:))), ...
    'Expected a real, finite, four-by-four delay matrix in seconds.');
dy=ell(:,1)-ell(1,1);du=ell(1,:)-ell(1,1);
additiveResidual=ell-(ell(1,1)+dy+du);
assert(norm(additiveResidual,'fro')<1e-10,'Channel delays do not have the assumed input/output additive form.');
extraY=max([0;dy])-dy;extraU=max([0 du])-du;
assert(all(extraY>=-1e-14) && all(extraU>=-1e-14),'Additional alignment delays must be nonnegative.');
phaseY=exp(-1i*2*pi*f.*extraY.');phaseU=exp(-1i*2*pi*f.*extraU);
codes={'B0','B1','B4'};methodNames={'Nominal static','Estimated static','Estimated static + relative-delay alignment'};
nullYs={ny0,ny1,ny1};nullUs={nu0,nu1,nu1};methods=struct([]);
rawValidation=reshape(calib.H_validation.',4,4,[]);rawTraining=reshape(calib.H_training.',4,4,[]);
rawModel=reshape(joint.H.',4,4,[]);
hasLow=isfield(joint,'H_low_undelayed') && isequal(size(joint.H_low_undelayed),[n 16]);
if hasLow
    % Self-fitted low-frequency contribution is retained ONLY as a diagnostic.
    assert(all(isfinite(joint.H_low_undelayed(:))),'Nonfinite fitted low-frequency diagnostic.');
    lowDelayed=joint.H_low_undelayed.*exp(-1i*2*pi*f.*ell(:).');
    rawLow=reshape(lowDelayed.',4,4,[]);
end
outputRowNorms=zeros(n,3);inputColumnNorms=zeros(n,3);outputOperatorNorms=zeros(n,3);
inputOperatorNorms=zeros(n,3);unitarity=zeros(n,2);
for m=1:3
    left=[Pnom;nullYs{m}'];right=[Qnom,nullUs{m}];
    Hv=complex(zeros(4,4,n));Ht=Hv;Hm=Hv;Hl=Hv;
    for k=1:n
        if m==3,Cy=diag(phaseY(k,:));Cu=diag(phaseU(k,:));else,Cy=eye(4);Cu=eye(4);end
        L=left*Cy;R=Cu*right;
        Hv(:,:,k)=L*rawValidation(:,:,k)*R;
        Ht(:,:,k)=L*rawTraining(:,:,k)*R;
        Hm(:,:,k)=L*rawModel(:,:,k)*R;
        if hasLow,Hl(:,:,k)=L*rawLow(:,:,k)*R;end
        outputRowNorms(k,m)=norm(nullYs{m}'*Cy,2);
        inputColumnNorms(k,m)=norm(Cu*nullUs{m},2);
        outputOperatorNorms(k,m)=norm(L,2);inputOperatorNorms(k,m)=norm(R,2);
        if m==3,unitarity(k,:)=[norm(Cy'*Cy-eye(4),'fro'),norm(Cu'*Cu-eye(4),'fro')];end
    end
    entry=struct('code',codes{m},'description',methodNames{m},'Ny',nullYs{m},'Nu',nullUs{m}, ...
        'Ty_static',left,'Tu_static',right,'H_validation',Hv,'H_training',Ht,'H_model_joint',Hm);
    if hasLow,entry.H_low_model_diagnostic=Hl;end
    if m==1,methods=repmat(entry,1,3);end
    methods(m)=entry;
end
assert(max(abs(outputRowNorms(:)-1))<1e-10 && max(abs(inputColumnNorms(:)-1))<1e-10, ...
    'Unit input/output normalization failed.');
assert(max(unitarity,[],'all')<1e-10,'Pure-delay alignment is not unitary.');
assert(max(abs(outputOperatorNorms(:,3)-outputOperatorNorms(:,2)))<1e-10 && ...
    max(abs(inputOperatorNorms(:,3)-inputOperatorNorms(:,2)))<1e-10, ...
    'B4 changed the B1 induced norm despite pure unitary phase alignment.');

% Original-bin, actual-frequency integration. Block norms are sums across
% their three one-at-a-time unit logical excitations/observations, not averages.
bands=[3 25;25 80;100 170;260 310;25 800];
bandNames={'3-25','25-80','100-170','260-310','25-800'};
blockNames={'logical_inputs_to_N','N_to_logical_outputs','N_to_N'};
rows={};lowRows={};rmsValues=zeros(3,5,3);referenceRms=zeros(5,1);
for b=1:5
    mask=frequencyMask(f,bands(b,:));fb=f(mask);
    assert(nnz(mask)>=2,'Insufficient frequency bins in prescribed band %s.',bandNames{b});
    referenceRms(b)=integratedRms(fb,calib.H_validation(mask,:));
    assert(referenceRms(b)>0,'Zero reference energy in band %s.',bandNames{b});
    for m=1:3
        for block=1:3
            z=blockValues(methods(m).H_validation,block);z=z(mask,:);
            rmsValue=integratedRms(fb,z);rmsValues(m,b,block)=rmsValue;
            rows(end+1,:)={codes{m},bandNames{b},blockNames{block},bands(b,1),bands(b,2), ...
                fb(1),fb(end),nnz(mask),rmsValue,referenceRms(b),100*rmsValue/referenceRms(b)}; %#ok<AGROW>
            if hasLow
                diagnostic=blockValues(methods(m).H_low_model_diagnostic,block);
                lowRows(end+1,:)={codes{m},bandNames{b},blockNames{block},integratedRms(fb,diagnostic(mask,:)), ...
                    'Self-fitted low component; algebra diagnostic, not validation evidence'}; %#ok<AGROW>
            end
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'method','band_hz','block','requested_fmin_hz','requested_fmax_hz', ...
    'actual_fmin_hz','actual_fmax_hz','frequency_count','validation_response_rms', ...
    'fixed_original_validation_frobenius_rms','percent_of_fixed_original_validation_rms'});
pairRows={};
for b=1:5
    for block=1:3
        a=rmsValues(2,b,block);c=rmsValues(3,b,block);
        pairRows(end+1,:)={bandNames{b},blockNames{block},a,c,c/max(a,realmin), ...
            20*log10(max(a,realmin)/max(c,realmin))}; %#ok<AGROW>
    end
end
pairSummary=cell2table(pairRows,'VariableNames',{'band_hz','block','B1_rms','B4_rms','B4_over_B1_rms', ...
    'B1_to_B4_reduction_db_positive_means_less_response'});

% Fitting quality for fixed/joint candidates is reported, not used to select.
fitRows={};for name={'fixed','joint'}
    fit=calib.fits.(name{1});certificate=candidateCertificate(fit);
    fitRows(end+1,:)={name{1},fit.exitflag,fit.grad_inf,certificate.kkt_firstorderopt, ...
        certificate.feasible,certificate.spatial_rank_valid,certificate.passed, ...
        certificate.normalized_spatial_sigma3_over1(1),certificate.normalized_spatial_sigma3_over1(2), ...
        'Unconstrained full gradient; includes frozen delay coordinates for fixed fit; not the convergence gate',fit.train_objective, ...
        fit.train_weighted_error_pct,relativeIntegratedError(f,fit.H,calib.H_training), ...
        relativeIntegratedError(f,fit.H,calib.H_validation), ...
        directionAngle(fit.Ny,ny0),directionAngle(fit.Nu,nu0)}; %#ok<AGROW>
end
fitMetrics=cell2table(fitRows,'VariableNames',{'candidate','solver_exitflag','gradient_infinity', ...
    'kkt_firstorderopt','feasible','spatial_rank_valid','numerical_certificate_passed', ...
    'normalized_output_spatial_sigma3_over1','normalized_input_spatial_sigma3_over1','full_gradient_interpretation', ...
    'training_objective','reported_training_weighted_error_percent','training_integrated_relative_error_percent', ...
    'validation_integrated_relative_error_percent','output_null_angle_to_nominal_deg','input_null_angle_to_nominal_deg'});
directionMetrics=table( ...
    {'Nominal to fixed';'Nominal to joint';'Fixed to joint'}, ...
    [directionAngle(ny0,nyFixed);directionAngle(ny0,ny1);directionAngle(nyFixed,ny1)], ...
    [directionAngle(nu0,nuFixed);directionAngle(nu0,nu1);directionAngle(nuFixed,nu1)], ...
    'VariableNames',{'comparison','output_null_angle_deg','input_null_angle_deg'});
jointCertificate=candidateCertificate(joint);
provisional=~jointCertificate.passed;
statusText='';if provisional,statusText=' [provisional: numerical certificate incomplete]';end
files={};prefix=['calibration_' tag];
path=fullfile(out,[prefix '_metrics.csv']);writetable(metrics,path);files{end+1}=path;
path=fullfile(out,[prefix '_pair_summary.csv']);writetable(pairSummary,path);files{end+1}=path;
path=fullfile(out,[prefix '_fit_metrics.csv']);writetable(fitMetrics,path);files{end+1}=path;
path=fullfile(out,[prefix '_direction_angles.csv']);writetable(directionMetrics,path);files{end+1}=path;
if hasLow
    lowMetrics=cell2table(lowRows,'VariableNames',{'method','band_hz','block','model_low_rms','interpretation'});
    path=fullfile(out,[prefix '_model_low_diagnostic.csv']);writetable(lowMetrics,path);files{end+1}=path;
else,lowMetrics=table;end

% The primary figures all use the validation FRF, never the fitted low term.
colors=[.55 .55 .55;.03 .27 .52;.50 .18 .59];
legendLabels={'B0: nominal static','B1: estimated static','B4: same static + relative delays'};
fig=figure('Visible','off','Color','w','Position',[60 60 1250 1000]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
displayRange=[max(3,min(f)) min(800,max(f))];plotMask=frequencyMask(f,displayRange);
titles={'Nominal logical inputs to nullspace output','Nullspace input to nominal logical outputs','Nullspace input to nullspace output'};
for block=1:3
    ax=nexttile(tl);hold(ax,'on');handles=gobjects(0);
    for m=1:3
        z=blockValues(methods(m).H_validation,block);value=sqrt(sum(abs(z).^2,2));
        handles(end+1)=semilogx(ax,f(plotMask),toDb(value(plotMask)),'Color',colors(m,:),'LineWidth',1.15); %#ok<AGROW>
    end
    style(ax,displayRange);ylabel(ax,'Response magnitude (dB)');title(ax,titles{block},'FontWeight','normal');
    if block==3,xlabel(ax,'Frequency (Hz)');end
    if block==1,legendHandles=handles;end
end
lg=legend(legendHandles,legendLabels,'Orientation','horizontal','Box','off','FontSize',10);lg.Layout.Tile='north';
title(tl,sprintf('Validation cross-block responses: z = %g mm%s',calib.position_mm,statusText),'FontWeight','normal');
subtitle(tl,{'Identical nominal axes and unit-norm input/output directions; B1 and B4 share the same null vectors.', ...
    'Observed low-frequency residuals include spring dynamics; they are not an isolated true rigid-body contribution.'},'FontSize',10);
files=[files,exportFigure(fig,out,[prefix '_validation_blocks'])];clear cleanup;

files=[files,plotNN(f,methods,colors,legendLabels,displayRange,calib.position_mm,statusText,out,[prefix '_NN_response'])];
files=[files,plotNN(f,methods,colors,legendLabels,[260 310],calib.position_mm,statusText,out,[prefix '_target_band_260_310'])];

axisLabels={'Rx','Ry','Z','N'};B4validation=methods(3).H_validation;B4model=methods(3).H_model_joint;
for phase=[false true]
    fig=figure('Visible','off','Color','w','Position',[40 40 1460 1040]);
    cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ax=nexttile(tl);hold(ax,'on');hv=reshape(B4validation(y,u,plotMask),[],1);hm=reshape(B4model(y,u,plotMask),[],1);
            if phase
                vm=unwrap(angle(hm))*180/pi;vv=unwrap(angle(hv))*180/pi;vv=vv+360*round((vm(1)-vv(1))/360);
            else,vv=toDb(hv);vm=toDb(hm);end
            h1=semilogx(ax,f(plotMask),vv,'Color',[.65 .65 .65],'LineWidth',.7);
            h2=semilogx(ax,f(plotMask),vm,'Color',colors(3,:),'LineWidth',1.1);
            style(ax,displayRange);title(ax,[axisLabels{u} ' \rightarrow ' axisLabels{y}], ...
                'Interpreter','tex','FontWeight','normal');
            if y<4,xticklabels(ax,{});end
            if y==1 && u==1,handles=[h1 h2];end
        end
    end
    lg=legend(handles,{'B4 validation LRM','B4 fitted model'},'Orientation','horizontal','Box','off','FontSize',11);lg.Layout.Tile='north';
    xlabel(tl,'Frequency (Hz)','FontSize',12);
    if phase,name='phase';ylabel(tl,'Unwrapped phase (deg)','FontSize',12);
    else,name='magnitude';ylabel(tl,'Magnitude (dB, recorded normalized coordinates)','FontSize',12);end
    title(tl,sprintf('B4 aligned coordinates: z = %g mm%s',calib.position_mm,statusText),'FontWeight','normal');
    subtitle(tl,'Offline exact frequency-domain alignment; original frequency bins, no smoothing or fractional-delay implementation.','FontSize',10);
    files=[files,exportFigure(fig,out,[prefix '_B4_full_' name])]; %#ok<AGROW>
    clear cleanup;
end

evaluationPath=fullfile(out,[prefix '_evaluation.mat']);files{end+1}=evaluationPath;
evaluation=struct('position_mm',calib.position_mm,'tag',tag,'files',{files},'metrics',metrics,'pair_summary',pairSummary, ...
    'fit_metrics',fitMetrics,'direction_metrics',directionMetrics,'low_model_diagnostic',lowMetrics, ...
    'f',f,'methods',methods,'nominal_output_rows_unit_norm',Pnom,'nominal_input_columns_unit_norm',Qnom, ...
    'relative_output_delay_seconds',dy,'relative_input_delay_seconds',du, ...
    'additional_causal_output_delay_seconds',extraY,'additional_causal_input_delay_seconds',extraU, ...
    'common_delay_after_alignment_seconds',ell(1,1)+max([0;dy])+max([0 du]), ...
    'joint_candidate_provisional',provisional, ...
    'numerical_diagnostics',struct('joint_certificate',jointCertificate, ...
        'fixed_certificate',candidateCertificate(fixed), ...
        'start_agreement_deg',optionalField(calib,'start_agreement_deg',[NaN NaN]), ...
        'high_scaled_QR_condition',optionalField(calib,'high_scaled_QR_condition',NaN), ...
        'gate','Positive exit flag, feasible delay constraints, rank-three spatial factors, and finite KKT first-order optimality <= 1e-6', ...
        'full_gradient_interpretation','grad_inf is the unconstrained full gradient, including frozen delay coordinates in the fixed-delay fit; it is not a constrained stationarity or convergence test'), ...
    'observed_target_band',struct('band_hz',[260 310],'B1_NN_rms',rmsValues(2,4,3), ...
        'B4_NN_rms',rmsValues(3,4,3),'B4_over_B1_rms',rmsValues(3,4,3)/max(rmsValues(2,4,3),realmin), ...
        'interpretation','Observed validation response, not isolated flexible-mode retention'), ...
    'checks',struct('delay_additivity_residual_seconds',norm(additiveResidual,'fro'), ...
        'max_unitarity_error',max(unitarity,[],'all'),'maximum_null_output_row_norm_error',max(abs(outputRowNorms(:)-1)), ...
        'maximum_null_input_column_norm_error',max(abs(inputColumnNorms(:)-1)), ...
        'unit_white_sensor_noise_power_gain_range',[min(outputRowNorms(:).^2),max(outputRowNorms(:).^2)], ...
        'null_direction_max_imaginary_part',max(abs(imag([ny0;nu0;ny1;nu1]))), ...
        'B1_B4_output_operator_norm_difference',max(abs(outputOperatorNorms(:,3)-outputOperatorNorms(:,2))), ...
        'B1_B4_input_operator_norm_difference',max(abs(inputOperatorNorms(:,3)-inputOperatorNorms(:,2))), ...
        'B1_B4_output_null_difference',norm(methods(2).Ny-methods(3).Ny), ...
        'B1_B4_input_null_difference',norm(methods(2).Nu-methods(3).Nu)), ...
    'normalization','One frozen nominal output row/input column per axis, unit Euclidean norm; recorded-command metric Ru=I', ...
    'block_convention','Vector block energy sums the three unit logical-axis responses; no hidden per-method scaling', ...
    'reference_denominator','Integrated Frobenius RMS of the same original untransformed validation FRF in each band', ...
    'noise_statement','Algebraic unit-norm gain only: pure diagonal delays preserve B1 gain for independent diagonal sensor-noise spectra. Correlated noise may differ; no measured noise gain is claimed.', ...
    'implementation_scope','Exact frequency-domain offline causal alignment. No FIR, digital fractional-delay implementation, controller intervention, or new experiment.', ...
    'validation_scope','Second-half closed-loop LRM FRFs; no external excitation instrument, no claim of unbiased plant/noise estimation or independent experimental validation.');
save(evaluationPath,'evaluation','-v7.3');
disp(pairSummary);disp(fitMetrics);disp(evaluation.checks);
end

function certificate=candidateCertificate(fit)
% Missing numerical certification is provisional, not evidence of convergence.
kkt=optionalField(fit,'kkt_firstorderopt',NaN);
feasible=optionalField(fit,'feasible',false);
rankValid=optionalField(fit,'spatial_rank_valid',false);
feasible=isscalar(feasible) && isfinite(feasible) && feasible==1;
rankValid=isscalar(rankValid) && isfinite(rankValid) && rankValid==1;
sigma=optionalField(fit,'normalized_spatial_sigma3_over1',[NaN NaN]);
if numel(sigma)~=2,sigma=[NaN NaN];end
passed=isscalar(fit.exitflag) && isfinite(fit.exitflag) && fit.exitflag>0 && ...
    feasible && rankValid && isscalar(kkt) && isfinite(kkt) && kkt>=0 && kkt<=1e-6;
certificate=struct('kkt_firstorderopt',kkt,'feasible',feasible, ...
    'spatial_rank_valid',rankValid,'normalized_spatial_sigma3_over1',sigma(:).','passed',passed);
end

function value=optionalField(source,name,fallback)
if isfield(source,name),value=source.(name);else,value=fallback;end
end

function files=plotNN(f,methods,colors,labels,range,position,statusText,out,stem)
mask=frequencyMask(f,range);assert(nnz(mask)>=2,'Requested NN display band is unavailable.');
fm=f(mask);h=complex(zeros(nnz(mask),3));
for m=1:3,h(:,m)=reshape(methods(m).H_validation(4,4,mask),[],1);end
model=reshape(methods(3).H_model_joint(4,4,mask),[],1);
fig=figure('Visible','off','Color','w','Position',[60 60 1240 850]);
cleanup=onCleanup(@()closeOwned(fig));tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
for row=1:2
    ax=nexttile(tl);hold(ax,'on');handles=gobjects(0);
    if row==1,value=toDb(h);vm=toDb(model);
    else
        value=unwrap(angle(h),[],1)*180/pi;vm=unwrap(angle(model))*180/pi;
        for m=1:2,value(:,m)=value(:,m)+360*round((value(1,3)-value(1,m))/360);end
        vm=vm+360*round((value(1,3)-vm(1))/360);
    end
    for m=1:3,handles(end+1)=semilogx(ax,fm,value(:,m),'Color',colors(m,:),'LineWidth',1.2);end %#ok<AGROW>
    handles(end+1)=semilogx(ax,fm,vm,'--','Color',[.24 .07 .30],'LineWidth',1.0);
    style(ax,range);
    if row==1,ylabel(ax,'Observed response magnitude (dB)');legendHandles=handles;
    else,ylabel(ax,'Unwrapped phase (deg)');xlabel(ax,'Frequency (Hz)');end
end
lg=legend(legendHandles,[strcat(labels,' (validation)'),{'B4 fitted model'}], ...
    'Orientation','horizontal','Box','off','FontSize',9);lg.Layout.Tile='north';
title(tl,sprintf('N to N observed response: z = %g mm%s',position,statusText),'FontWeight','normal');
subtitle(tl,{'Curves contain complete measured responses, not isolated flexible-mode components.', ...
    'Independent phase unwrap; at most one constant whole-curve 360-deg offset per curve. B4 is exact offline phase alignment.'},'FontSize',10);
files=exportFigure(fig,out,stem);clear cleanup;
end

function z=blockValues(H,index)
n=size(H,3);
if index==1,z=reshape(H(4,1:3,:),3,n).';
elseif index==2,z=reshape(H(1:3,4,:),3,n).';
else,z=reshape(H(4,4,:),n,1);end
end

function value=integratedRms(f,z)
assert(numel(f)>=2 && f(end)>f(1),'RMS integration requires at least two distinct frequency points.');
value=sqrt(trapz(f,sum(abs(z).^2,2))/(f(end)-f(1)));
end

function value=relativeIntegratedError(f,model,reference)
value=100*integratedRms(f,model-reference)/max(integratedRms(f,reference),realmin);
end

function mask=frequencyMask(f,range)
mask=f>=range(1)-1e-8 & f<=range(2)+1e-8;
end

function v=unit(v)
v=v(:);assert(all(isfinite(v)) && norm(v)>0,'Null direction must be finite and nonzero.');v=v/norm(v);
end

function v=alignedUnit(v,reference)
v=unit(v);a=reference'*v;if abs(a)>0,v=v*exp(-1i*angle(a));end
end

function value=directionAngle(a,b)
a=unit(a);b=unit(b);value=acosd(min(1,max(0,abs(a'*b))));
end

function value=toDb(H)
value=20*log10(max(abs(H),realmin));
end

function style(ax,range)
set(ax,'XScale','log','FontName','Times New Roman','FontSize',10,'Box','on','LineWidth',.65, ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',.20,'MinorGridAlpha',.10);
grid(ax,'on');xlim(ax,range);
if range(1)>=260 && range(2)<=310,xticks(ax,[260 270 280 290 300 310]);
else,ticks=[3 10 25 80 170 310 800];xticks(ax,ticks(ticks>=range(1) & ticks<=range(2)));end
end

function files=exportFigure(fig,out,stem)
png=fullfile(out,[stem '.png']);pdf=fullfile(out,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',200);exportgraphics(fig,pdf,'ContentType','vector');files={png,pdf};
end

function closeOwned(fig)
if isgraphics(fig),close(fig);end
end
