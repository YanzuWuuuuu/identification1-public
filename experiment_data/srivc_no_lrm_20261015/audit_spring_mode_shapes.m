function audit = audit_spring_mode_shapes(resultPath,rawPath,outputDir)
% Read-only source audit of complex residues and recorded output coordinates.
% No fitting of poles/delays, no physical mode labels, and no input overwrites.
% Defaults can be rerun after spring_selected.mat has been updated.
base=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(resultPath)
    resultPath=fullfile(base,'spring_identification_20261016','spring_selected.mat');
end
if nargin<2 || isempty(rawPath)
    rawPath=fullfile(fileparts(base),'auto_identify_motor1234_z0.mat');
end
if nargin<3 || isempty(outputDir)
    outputDir=fullfile(fileparts(resultPath),'mode_shape_diagnostics');
end
S=load(resultPath,'result'); result=S.result; fit=result.selected.fit;
nm=numel(fit.reference_hz); C=fit.coefficients; a=fit.a;
assert(isequal(size(a),[nm 2]) && size(C,1)>=2*nm && size(C,2)==16, ...
    'Expected additive second-order model with 16 columns, y+4*(u-1).');
ell=result.delay_matrix_seconds;
assert(isequal(size(ell),[4 4]) && all(isfinite(ell(:))) && isreal(ell), ...
    'delay_matrix_seconds must be a real finite 4-by-4 matrix.');
assert(all(a(:,2)>0) && all(4*a(:,2)-a(:,1).^2>0), ...
    'This audit expects the underdamped components enforced by frf_srivc_fit.');
if ~exist(outputDir,'dir'),mkdir(outputDir);end

% Snapshot of matlab/0/Multiple_main_Control.m:78-85, not measured geometry.
C1=[.1714 .0861 1.0014;.1697 -.0756 .9986; ...
    -.1708 -.1108 1.0012;-.1673 .1034 .9989];
C2=[-.5297 .5215 .6477 .4121 -.9896; ...
    .5585 .5490 .5350 .0562 1.0118; ...
    -.6198 .4402 .6328 -.0175 1.0857; ...
    .6258 .4527 .6219 .3133 -.9832];
referenceHz=[286.43 344.58 488.27 705.67 738.83];
sensorXY=[-.095 .160;.095 .160;.095 -.160;-.095 -.160];
Cnom=[sensorXY(:,2),-sensorXY(:,1),ones(4,1)];
nNom=[1;-1;1;-1]/2; Q1=orth(C1);
recorded=recordedMap(rawPath,outputDir);

Rational=complex(zeros(4,4,nm)); Observed=Rational;
leftVectors=complex(zeros(4,nm,2)); rightVectors=leftVectors;
poles=complex(zeros(nm,1)); singularValues=zeros(4,nm,2);
metrics=cell(2*nm,20); shapes=cell(8*nm,10); row=0; vrow=0;
views={'rational_delay_removed','observed_delay_included'};
for m=1:nm
    wr=2*pi*fit.reference_hz(m);
    qp=(-a(m,1)+1i*sqrt(4*a(m,2)-a(m,1)^2))/(2*a(m,2));
    p=wr*qp; poles(m)=p;
    dAds=(a(m,1)+2*a(m,2)*qp)/wr;
    % q=s/wr, so this is dA/ds, not dA/dq. Constant D has no residue.
    R=reshape((C(2*m-1,:)+C(2*m,:)*qp)/dAds,4,4);
    Rational(:,:,m)=R;
    % The analytic factor is evaluated at the complex pole, not at i*omega.
    Observed(:,:,m)=R.*exp(-p*ell);
    for v=1:2
        if v==1,Rv=R;else,Rv=Observed(:,:,m);end
        [U,SV,V]=svd(Rv); sv=diag(SV); nr=norm(Rv,'fro');
        if nr>0
            [~,pivot]=max(abs(U(:,1))); rot=exp(-1i*angle(U(pivot,1)));
            u=U(:,1)*rot; vv=V(:,1)*rot;
            rankRatio=sv(2)/sv(1); firstEnergy=sv(1)^2/sum(sv.^2);
            nominalNull=norm(nNom'*Rv,'fro')/nr;
            c1Inside=norm(Q1'*Rv,'fro')/nr;
            c1Outside=norm(Rv-Q1*(Q1'*Rv),'fro')/nr;
            mac=abs(C2'*u).^2./(sum(abs(C2).^2,1)'*sum(abs(u).^2));
            if strcmp(recorded.status,'completed')
                recordInside=norm(recorded.row_basis'*Rv,'fro')/nr;
                recordNull=norm(recorded.null_basis'*Rv,'fro')/nr;
            else
                recordInside=NaN;recordNull=NaN;
            end
        else
            u=complex(nan(4,1)); vv=u; rankRatio=NaN;firstEnergy=NaN;
            nominalNull=NaN;c1Inside=NaN;c1Outside=NaN;mac=nan(5,1);
            recordInside=NaN;recordNull=NaN;
        end
        leftVectors(:,m,v)=u;rightVectors(:,m,v)=vv;singularValues(:,m,v)=sv;
        row=row+1;
        metrics(row,:)=[{m,views{v},fit.frequency_hz(m),fit.damping(m), ...
            real(p),imag(p)/(2*pi),nr,rankRatio,firstEnergy,nominalNull, ...
            c1Inside,c1Outside,recordInside,recordNull},num2cell(mac.'), ...
            {abs(1+a(m,1)*qp+a(m,2)*qp^2)}];
        for k=1:4
            vrow=vrow+1;
            shapes(vrow,:)={m,views{v},k,real(u(k)),imag(u(k)),abs(u(k)), ...
                angle(u(k))*180/pi,real(vv(k)),imag(vv(k)),abs(vv(k))};
        end
    end
end
metricNames={'component','coordinate_view','natural_frequency_hz','damping_ratio', ...
    'pole_real_per_s','damped_frequency_hz','residue_frobenius','sigma2_over_sigma1', ...
    'first_singular_energy_fraction','nominal_null_amplitude_fraction', ...
    'C1_columnspace_amplitude_fraction','C1_outside_amplitude_fraction', ...
    'recorded_map_rowspace_amplitude_fraction','recorded_map_null_amplitude_fraction', ...
    'MAC_C2_mode1','MAC_C2_mode2','MAC_C2_mode3','MAC_C2_mode4','MAC_C2_mode5', ...
    'denominator_root_residual'};
metricTable=cell2table(metrics,'VariableNames',metricNames);
shapeTable=cell2table(shapes,'VariableNames',{'component','coordinate_view','channel', ...
    'left_vector_real','left_vector_imag','left_vector_amplitude','left_vector_phase_deg', ...
    'right_vector_real','right_vector_imag','right_vector_amplitude'});
writetable(metricTable,fullfile(outputDir,'modal_residue_metrics.csv'));
writetable(shapeTable,fullfile(outputDir,'modal_singular_vectors.csv'));
source=dir(resultPath);
audit=struct('source_result',resultPath,'source_bytes',source.bytes, ...
    'source_modified',source.date,'created',datestr(now,31), ...
    'fit_converged',fit.converged,'fit_score',fit.score,'poles',poles, ...
    'delay_matrix_seconds',ell,'delay_additivity_residual_seconds', ...
    max(abs(ell-ell(:,1)-ell(1,:)+ell(1,1)),[],'all'), ...
    'rational_residues',Rational,'observed_residues',Observed, ...
    'left_vectors',leftVectors,'right_vectors',rightVectors, ...
    'singular_values',singularValues,'metrics',metricTable,'shapes',shapeTable, ...
    'nominal_sensor_xy',sensorXY,'nominal_Cr',Cnom,'nominal_null',nNom, ...
    'original_C1',C1,'original_C2',C2,'original_C2_frequency_hz',referenceHz, ...
    'recorded_mapping',recorded);
save(fullfile(outputDir,'mode_shape_audit.mat'),'audit','-v7.3');
writeNotes(audit,outputDir);
fprintf('Mode-shape screening written to %s\n',outputDir);
disp(metricTable(:,[1:4,8:15]));
end

function rec=recordedMap(rawPath,outputDir)
rec=struct('status','not_available','raw_path',rawPath);
yn={'pre_y_0__0_','pre_y_1__0_','pre_y_2__0_','pre_y_3__0_'};
zn={'y_0__0_','y_0__1_','y_0__2_'};
if ~isfile(rawPath),rec.reason='Raw z0 file not found.';return;end
metadata=whos('-file',rawPath);rec.variable_metadata=metadata;
if ~all(ismember([yn zn],{metadata.name}))
    rec.reason='Required seven output records are not all present.';return;
end
try
    S=load(rawPath,yn{:},zn{:});
catch err
    rec.reason=['MATLAB could not load the optional output records: ',err.message];
    return;
end
t=S.(yn{1})(:,1);n=numel(t);dt=median(diff(t));
Y=zeros(n,4);Z=zeros(n,3);maxTimeError=0;
assert(n>10 && all(isfinite(t)) && all(diff(t)>0),'Invalid raw timestamp record.');
names=[yn zn];
for k=1:7
    x=S.(names{k});
    assert(isequal(size(x),[n 2]) && all(isfinite(x(:))), ...
        'Expected synchronous finite N-by-2 output records.');
    maxTimeError=max(maxTimeError,max(abs(x(:,1)-t)));
    if k<=4,Y(:,k)=x(:,2);else,Z(:,k-4)=x(:,2);end
end
assert(maxTimeError<dt*1e-3,'Output timestamps differ; no interpolation is applied.');
% Exclude the same initial second used by the existing spectral analysis.
keep=find(t-t(1)>=1);Y=Y(keep,:);Z=Z(keep,:);n=size(Y,1);h=floor(n/2);
it=1:h;iv=h+1:n;
allMap=affineMap(Y,Z);trainMap=affineMap(Y(it,:),Z(it,:));
secondMap=affineMap(Y(iv,:),Z(iv,:));
trainError=relativeError(Z(it,:),Y(it,:)*trainMap.L'+trainMap.offset');
validationError=relativeError(Z(iv,:),Y(iv,:)*trainMap.L'+trainMap.offset');
allError=relativeError(Z,Y*allMap.L'+allMap.offset');
rowChange=vecnorm(secondMap.L-trainMap.L,2,2)./max(vecnorm(trainMap.L,2,2),realmin);
[~,sv,V]=svd(allMap.L);ss=diag(sv);
mapRank=sum(ss>max(size(allMap.L))*eps(max(ss)));
rec.status='completed';rec.physical_output_names=yn;rec.recorded_output_names=zn;
rec.mapping=allMap.L;rec.offset=allMap.offset;rec.mapping_rank=mapRank;
rec.row_basis=V(:,1:mapRank);rec.null_basis=V(:,mapRank+1:end);
rec.all_samples=allMap;rec.first_half=trainMap;rec.second_half=secondMap;
rec.train_relative_error=trainError;rec.validation_relative_error=validationError;
rec.all_relative_error=allError;rec.mapping_row_relative_change=rowChange;
rec.sample_count=n;rec.first_half_samples=h;rec.timestamp_error_seconds=maxTimeError;
rec.equation='Z = Y*L'' + ones(N,1)*offset''; variable names are not physical-axis labels.';
rec.interpretation='Row/null space of a recorded signal map, not an independently calibrated rigid-body subspace.';
rows=cell(9,7);maps={allMap,trainMap,secondMap};tags={'all','first_half','second_half'};
for k=1:3
    for j=1:3
        rows{3*(k-1)+j,1}=tags{k};rows{3*(k-1)+j,2}=zn{j};
        rows(3*(k-1)+j,3:6)=num2cell(maps{k}.L(j,:));
        rows{3*(k-1)+j,7}=maps{k}.offset(j);
    end
end
mapTable=cell2table(rows,'VariableNames',{'mapping_set','recorded_signal', ...
    'pre_y_0','pre_y_1','pre_y_2','pre_y_3','offset'});
writetable(mapTable,fullfile(outputDir,'recorded_output_affine_map.csv'));
residualTable=table(zn(:),trainError(:),validationError(:),allError(:),rowChange(:), ...
    'VariableNames',{'recorded_signal','first_half_relative_residual', ...
    'second_half_prediction_relative_residual','all_samples_relative_residual', ...
    'split_mapping_row_relative_change'});
writetable(residualTable,fullfile(outputDir,'recorded_output_mapping_residuals.csv'));
nullNames=[{'sensor_channel'},arrayfun(@(k)sprintf('null_%d',k), ...
    1:size(rec.null_basis,2),'UniformOutput',false)];
writetable(array2table([(1:4)',rec.null_basis],'VariableNames',nullNames), ...
    fullfile(outputDir,'recorded_output_mapping_nullspace.csv'));
rec.map_table=mapTable;rec.residual_table=residualTable;
end

function m=affineMap(Y,Z)
muY=mean(Y,1);muZ=mean(Z,1);X=Y-muY;T=Z-muZ;
scale=sqrt(mean(X.^2,1));
assert(all(scale>0),'A physical sensor record is constant; mapping is not identifiable.');
Xs=X./scale;[U,S,V]=svd(Xs,'econ');sv=diag(S);
r=sum(sv>max(size(Xs))*eps(max(sv)));
beta=(V(:,1:r)*((U(:,1:r)'*T)./sv(1:r)))./scale';
L=beta';offset=(muZ-muY*beta)';
m=struct('L',L,'offset',offset,'design_rank',r,'scaled_design_singular_values',sv, ...
    'scaled_design_condition',sv(1)/max(sv(end),realmin));
end

function err=relativeError(Z,pred)
% Centered denominator prevents large DC positions masking residual error.
err=vecnorm(pred-Z,2,1)./max(vecnorm(Z-mean(Z,1),2,1),realmin);
end

function writeNotes(audit,outputDir)
fid=fopen(fullfile(outputDir,'mode_shape_sources_and_limits.md'),'w','n','UTF-8');
assert(fid>=0,'Cannot write diagnostic notes.');cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Modal residue and recorded-coordinate screening\n\n');
fprintf(fid,'Source result: `%s` (modified %s, %d bytes).\n\n', ...
    audit.source_result,audit.source_modified,audit.source_bytes);
fprintf(fid,'Audit created: %s. Fit converged=%d, score=%.6g. No poles, delays, or source data were changed.\n\n', ...
    audit.created,audit.fit_converged,audit.fit_score);
fprintf(fid,'For q=s/(2*pi*fref), A=1+a1*q+a2*q^2, the positive-imaginary pole p has residue R=(b0+b1*q(p))/((a1+2*a2*q(p))/(2*pi*fref)). Rows/columns are outputs/inputs; flattened column index is y+4*(u-1). A constant direct term has no pole residue.\n\n');
fprintf(fid,'Two views are reported: `rational_delay_removed` and `observed_delay_included`. The latter is R_observed=R.*exp(-p*ell), evaluated at the complex pole. Delay additivity residual is %.6g s. Removing estimated delay is model dependent, not independently verified hardware calibration.\n\n', ...
    audit.delay_additivity_residual_seconds);
fprintf(fid,'## Spatial references\n\n');
fprintf(fid,'- Nominal sensor coordinates come from `matlab/anlaysis/test.m:6-9`: (-.095,.160), (.095,.160), (.095,-.160), (-.095,-.160) m. The nominal null vector is [1,-1,1,-1]/2.\n');
fprintf(fid,'- C1 and five C2 columns are copied from `matlab/0/Multiple_main_Control.m:78-85`; their nominal modal frequencies are 286.43, 344.58, 488.27, 705.67, 738.83 Hz. These are the old simulation references, not experimental ground truth.\n');
fprintf(fid,'- Channel order, sensor gains, signs, timing, and exact hardware geometry are not confirmed. All geometric fractions and MAC values are exploratory.\n\n');
fprintf(fid,'## Meaning of metrics\n\n');
fprintf(fid,'`sigma2_over_sigma1` and `first_singular_energy_fraction` describe algebraic residue rank concentration. They do not establish a plate/suspension origin. Projection columns are amplitude norm ratios, not squared energy fractions. MAC compares the leading complex left singular vector with each original C2 column and is phase invariant. Near-equal leading singular values make the individual leading vector ambiguous.\n\n');
fprintf(fid,'Four sensors give only one nominal non-rigid direction. Small null projection does not prove rigid motion: a flexible shape may have rigid-compatible values at these four points. Large projection may also reflect geometry, gain, phase, or model errors. No component receives an automatic physical label.\n\n');
fprintf(fid,'## Recorded signal mapping\n\n');
fprintf(fid,'Raw source: `%s`. Status: %s.\n\n',audit.recorded_mapping.raw_path,audit.recorded_mapping.status);
if strcmp(audit.recorded_mapping.status,'completed')
    rec=audit.recorded_mapping;
    fprintf(fid,'Only the four pre_y records and three y_0__k_ records were loaded; metadata lists all %d file variables. The first second was excluded. A static affine map Z=Y*L''+offset'' was fitted on the first half and predicted on the second half. An all-sample map and both split maps are also saved. No coordinate names or physical units are inferred from variable names.\n\n',numel(rec.variable_metadata));
    fprintf(fid,'All-sample map rank=%d; standardized sensor design rank=%d, condition=%.6g. A deficient/ill-conditioned sensor design limits uniqueness/stability of the recovered map. Second-half relative residuals are [%s], with each denominator the corresponding centered output norm.\n\n', ...
        rec.mapping_rank,rec.all_samples.design_rank,rec.all_samples.scaled_design_condition, ...
        strtrim(sprintf(' %.6g',rec.validation_relative_error)));
    fprintf(fid,'The recorded-map row/null spaces represent an existing signal transformation. They do not independently establish the true rigid-body geometry. Their directly consistent comparison is the `observed_delay_included` residue view, because the affine map was fitted to raw sensor records. The rational view is also tabulated only as a coordinate-sensitivity check.\n\n');
else
    fprintf(fid,'%s\n\n',audit.recorded_mapping.reason);
end
fprintf(fid,'The same acquisition and previously reused split are being diagnosed; this is not a new independent validation experiment. Current calculations do not identify which hardware part creates the 124/128/138 Hz features.\n');
end
