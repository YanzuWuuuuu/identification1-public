function pack=build_static_logical_decoupling
% Static coordinate transforms of saved experimental FRFs and SRIVC models.
% This does not refit, filter, compensate delays, or change raw observations.
base=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(base));
out=fullfile(base,'static_decoupling_20261016');
if ~isfolder(out),mkdir(out);end
source=fullfile(base,'position_spring_identification_20261016','position_campaign.mat');
geometrySource=fullfile(root,'matlab','0','Multiple_main_Control.m');
S=load(source,'campaign');records=S.campaign.results;
positions=cellfun(@(r)r.position_mm,records);
[positions,order]=sort(positions);records=records(order);
assert(max(abs(positions-[-1 -.5 0 .5 1]))<1e-10,'Unexpected position records.');

% Read numeric assignments only. Never execute the simulation main script.
txt=fileread(geometrySource);
B1=readMatrix(txt,'B1');C1=readMatrix(txt,'C1');
Mr=diag([readScalar(txt,'Jrx'),readScalar(txt,'Jry'),readScalar(txt,'mz')]);
assert(isequal(size(B1),[3 4]) && isequal(size(C1),[4 3]));
assert(rank(B1)==3 && rank(C1)==3 && all(diag(Mr)>0));
sqrtMr=diag(sqrt(diag(Mr)));
ny=null(C1.');nu=null(B1);
assert(isequal(size(ny),[4 1]) && isequal(size(nu),[4 1]));
[~,iy]=max(abs(ny));ny=ny*sign(ny(iy));
[~,iu]=max(abs(nu));nu=nu*sign(nu(iu));
Ty=pinv(C1);Tu=pinv(B1);
views(1)=struct('tag','massbalanced','title','Square-root mass/inertia balanced', ...
    'Ty',[sqrtMr*Ty;ny.'],'Tu',[Tu*sqrtMr,nu]);
views(2)=struct('tag','geometric','title','Geometric logical axes', ...
    'Ty',[Ty;ny.'],'Tu',[Tu,nu]);
R0=C1*(Mr\B1);
checks.geometry=struct('left_inverse',norm(Ty*C1-eye(3),'fro'), ...
    'right_inverse',norm(B1*Tu-eye(3),'fro'), ...
    'output_null',norm(ny.'*C1),'input_null',norm(B1*nu), ...
    'massbalanced_rigid',norm(views(1).Ty*R0*views(1).Tu-diag([1 1 1 0]),'fro'), ...
    'geometric_rigid',norm(views(2).Ty*R0*views(2).Tu-blkdiag(inv(Mr),0),'fro'));
assert(all(structfun(@(x)x<1e-10,checks.geometry)),'Static geometry algebra failed.');
checks.model_reconstruction=zeros(1,5);checks.coordinate_reconstruction=zeros(2,5);
checks.scaling_equivalence=zeros(1,5);
checks.nominal_geometry_source=geometrySource;
pack=struct('positions_mm',positions,'normal_indices',2:5,'axis_labels',{{'Rx','Ry','Z','N'}}, ...
    'views',views,'records',{cell(1,5)},'output_dir',out, ...
    'path',fullfile(out,'decoupled_campaign.mat'),'B1',B1,'C1',C1,'Mr',Mr, ...
    'ny',ny,'nu',nu,'source',source,'geometry_source',geometrySource, ...
    'operations','Static Ty*G*Tu only; all identified delays retained; no compensation or refitting', ...
    'channel_mapping','pre_y_0..3 and u_Motor1..4 in numeric order; assumed S1..4/A1..4', ...
    'units','Nominal coordinate scaling of recorded commands; actuator-force calibration not verified', ...
    'validation_scope','Second-half LRM from closed-loop total commands; exploratory, not certified unbiased plant FRFs');

errorRows={};bands=[3 25;25 800.0001;100 170];bandLabels={'3-25','25-800','100-170'};
for k=1:numel(records)
    r=records{k};f=r.f(:);fit=r.selected.fit;
    H=transfer(fit,f,r.delay_matrix_seconds);
    checks.model_reconstruction(k)=norm(H-r.selected.H,'fro')/max(norm(r.selected.H,'fro'),realmin);
    assert(checks.model_reconstruction(k)<1e-10,'Model evaluation differs from saved prediction.');
    fd=logspace(log10(3),log10(800),8000).';
    for m=1:numel(fit.frequency_hz)
        fc=fit.frequency_hz(m);span=max(30*fc*fit.damping(m),.15);
        lo=max(3,fc-span);hi=min(800,fc+span);
        if hi>lo,fd=[fd;linspace(lo,hi,600).'];end %#ok<AGROW>
    end
    fd=unique([3;25;30;800;fd]);Hd=transfer(fit,fd,r.delay_matrix_seconds);
    H=reshape(H.',4,4,[]);Hd=reshape(Hd.',4,4,[]);
    Hv=reshape(r.validationFRF.',4,4,[]);Ht=reshape(r.trainingFRF.',4,4,[]);
    assert(all(isfinite(H(:))) && all(isfinite(Hv(:))) && all(isfinite(Hd(:))));
    record=struct('position_mm',r.position_mm,'anomalous',r.position_mm==-1, ...
        'status',r.status,'f',f,'f_dense',fd,'source_file',r.source_file, ...
        'record_warning',r.record_warning,'inner_converged',r.inner_converged, ...
        'delay_matrix_seconds',r.delay_matrix_seconds);
    for v=1:2
        ey=views(v).Ty;eu=views(v).Tu;
        a=struct('H_model',transform(H,ey,eu),'H_validation',transform(Hv,ey,eu), ...
            'H_training',transform(Ht,ey,eu),'H_dense',transform(Hd,ey,eu));
        record.views(v)=a;
        recovered=transform(a.H_validation,inv(ey),inv(eu));
        checks.coordinate_reconstruction(v,k)=norm(recovered(:)-Hv(:))/max(norm(Hv(:)),realmin);
        assert(checks.coordinate_reconstruction(v,k)<1e-10,'FRF coordinate inversion failed.');
        for b=1:size(bands,1)
            mask=f>=bands(b,1) & f<bands(b,2);fb=f(mask);
            for y=1:4
                for u=1:4
                    val=reshape(a.H_validation(y,u,mask),[],1);
                    err=reshape(a.H_model(y,u,mask),[],1)-val;
                    rmsVal=sqrt(trapz(fb,abs(val).^2)/(fb(end)-fb(1)));
                    rmsErr=sqrt(trapz(fb,abs(err).^2)/(fb(end)-fb(1)));
                    rel=100*rmsErr/max(rmsVal,realmin);
                    errorRows(end+1,:)={r.position_mm,views(v).tag, ...
                        pack.axis_labels{u},pack.axis_labels{y},bandLabels{b}, ...
                        rel,rmsErr,rmsVal,record.anomalous}; %#ok<AGROW>
                end
            end
        end
    end
    scale4=blkdiag(sqrtMr,1);
    scaled=transform(record.views(2).H_validation,scale4,scale4);
    target=record.views(1).H_validation;
    checks.scaling_equivalence(k)=norm(scaled(:)-target(:))/max(norm(target(:)),realmin);
    assert(checks.scaling_equivalence(k)<1e-10,'Mass scaling differs from pure constant scaling.');
    pack.records{k}=record;
end
pack.checks=checks;
errors=cell2table(errorRows,'VariableNames',{'position_mm','view','input_axis','output_axis', ...
    'band_hz','complex_relative_error_pct','complex_error_rms','reference_rms','anomalous_record'});
writetable(errors,fullfile(out,'decoupled_validation_errors.csv'));
pack.validation_errors=errors;
save(pack.path,'pack','-v7.3');
save(fullfile(out,'static_decoupling_checks.mat'),'checks');
for v=1:2
    writematrix(views(v).Ty,fullfile(out,[views(v).tag '_Ty_output.csv']));
    writematrix(views(v).Tu,fullfile(out,[views(v).tag '_Tu_input.csv']));
end
writematrix(B1,fullfile(out,'source_B1.csv'));
writematrix(C1,fullfile(out,'source_C1.csv'));
writematrix(Mr,fullfile(out,'mass_inertia.csv'));
disp(checks);
end

function H=transfer(fit,f,delay)
f=f(:);nm=numel(fit.reference_hz);
q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';
B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm
    B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);
end
if fit.direct,B(:,end)=1;end
H=(B*fit.coefficients).*exp(-1i*2*pi*f.*delay(:).');
end

function result=transform(H,Ty,Tu)
result=pagemtimes(pagemtimes(Ty,H),Tu);
end

function A=readMatrix(txt,name)
token=regexp(txt,['(?m)^\s*' name '\s*=\s*\[([^\]]+)\]\s*;'],'tokens','once');
assert(~isempty(token),'Missing numeric assignment for %s.',name);
body=regexprep(token{1},'%[^\r\n]*','');rows=strsplit(body,';');A=[];
for k=1:numel(rows)
    values=sscanf(strrep(rows{k},',',' '),'%f').';
    assert(~isempty(values),'Empty matrix row in %s.',name);
    A=[A;values]; %#ok<AGROW>
end
end

function value=readScalar(txt,name)
token=regexp(txt,['(?m)^\s*' name '\s*=\s*([0-9.eE+\-]+)\s*;'],'tokens','once');
assert(~isempty(token),'Missing numeric scalar %s.',name);
value=str2double(token{1});assert(isfinite(value));
end
