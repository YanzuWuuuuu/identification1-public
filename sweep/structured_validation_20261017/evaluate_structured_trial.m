function evaluation=evaluate_structured_trial(tag)
if nargin<1,tag='z0';end
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
pack=load(fullfile(out,['trial_' tag '.mat']));r=pack.result;f=r.f;
base=fullfile(root,'sweep','srivc_retest_20261016');
Ty=readmatrix(fullfile(base,'geometric_Ty_full.csv'));Tu=readmatrix(fullfile(base,'geometric_Tu_full.csv'));
bands=[3 25;25 80;100 170;260 310;25 800;320 400;470 515];
names={'SRIVC_free','Structured_fixed_delay','Structured_updated_delay','Train_check_FRF_gap'};
Hs={r.baseline_H,r.fixed.H,r.delay.H,r.H_training};blocks={'physical','logical','N_to_N','N_to_logical','logical_to_N'};
checkLogical=transform(r.H_validation,Ty,Tu);rows={};
for model=1:numel(names)
    pred=Hs{model};predLogical=transform(pred,Ty,Tu);
    for band=1:size(bands,1)
        mask=f>=bands(band,1)&f<=bands(band,2);
        for block=1:5
            if block==1,target=r.H_validation;prediction=pred;
            else,target=selectBlock(checkLogical,block);prediction=selectBlock(predLogical,block);end
            t=rmsIntegral(target(mask,:),f(mask));e=rmsIntegral(prediction(mask,:)-target(mask,:),f(mask));
            rows(end+1,:)={names{model},bands(band,1),bands(band,2),blocks{block},100*e/max(t,realmin),e,t}; %#ok<AGROW>
        end
    end
end
fitTable=cell2table(rows,'VariableNames',{'model','fmin_Hz','fmax_Hz','block','error_pct','absolute_complex_RMS_error','reference_complex_RMS'});
writetable(fitTable,fullfile(out,['prediction_metrics_' tag '.csv']));
% Keep nominal logical directions fixed and individually normalize all ports.
Py=Ty(1:3,:)./vecnorm(Ty(1:3,:),2,2);Pu=Tu(:,1:3)./vecnorm(Tu(:,1:3),2,1);
nyNom=Ty(4,:).';nyNom=nyNom/norm(nyNom);nuNom=Tu(:,4);nuNom=nuNom/norm(nuNom);
chosen=r.delay;ny=chosen.ny;nu=chosen.nu;
ny=ny*sign(ny.'*nyNom);nu=nu*sign(nu.'*nuNom);
tau=chosen.delay_matrix_seconds;tauU=tau(1,:);tauY=tau(:,1)-tau(1,1);
addedU=max(tauU)-tauU;addedY=max(tauY)-tauY;
% Frequency-by-channel order is output-fast, MATLAB column-major.
added=addedY+addedU;phase=exp(-2i*pi*f*added(:).');
calNames={'Nominal_static','Updated_static','Updated_static_and_delay'};
calH={transform(r.H_validation,[Py;nyNom.'],[Pu,nuNom]), ...
    transform(r.H_validation,[Py;ny.'],[Pu,nu]), ...
    transform(r.H_validation.*phase,[Py;ny.'],[Pu,nu])};
% A stronger, elementary baseline: minimize measured TRAINING cross energy
% in the already-prespecified 25--80 Hz band. No validation data enters it.
trainBand=find(f>=25 & f<=80);fb=f(trainBand);quadrature=zeros(numel(fb),1);
quadrature(1)=(fb(2)-fb(1))/2;quadrature(end)=(fb(end)-fb(end-1))/2;
quadrature(2:end-1)=(fb(3:end)-fb(1:end-2))/2;
Ay=zeros(4);Au=zeros(4);
for p=1:numel(trainBand)
    g=reshape(r.H_training(trainBand(p),:),4,4);gy=g*Pu;gu=Py*g;
    Ay=Ay+quadrature(p)*real(gy*gy');Au=Au+quadrature(p)*real(gu'*gu);
end
[vy,dy]=eig((Ay+Ay.')/2,'vector');[~,iy]=min(dy);nyBand=vy(:,iy);nyBand=nyBand*sign(nyBand.'*nyNom);
[vu,du]=eig((Au+Au.')/2,'vector');[~,iu]=min(du);nuBand=vu(:,iu);nuBand=nuBand*sign(nuBand.'*nuNom);
calNames{4}='Direct_band_static';
calH{4}=transform(r.H_validation,[Py;nyBand.'],[Pu,nuBand]);
rows={};
for band=1:size(bands,1)
    mask=f>=bands(band,1)&f<=bands(band,2);physical=rmsIntegral(r.H_validation(mask,:),f(mask));
    for block=2:5
        reference=selectBlock(calH{1},block);referenceRMS=rmsIntegral(reference(mask,:),f(mask));
        for candidate=1:numel(calNames)
            response=selectBlock(calH{candidate},block);rr=rmsIntegral(response(mask,:),f(mask));
            rows(end+1,:)={calNames{candidate},bands(band,1),bands(band,2),blocks{block},rr, ...
                20*log10(max(rr,realmin)),rr/max(referenceRMS,realmin),20*log10(max(rr,realmin)/max(referenceRMS,realmin)),rr/physical}; %#ok<AGROW>
        end
    end
end
calTable=cell2table(rows,'VariableNames',{'calibrator','fmin_Hz','fmax_Hz','block','response_complex_RMS', ...
    'response_dB_recorded_units','ratio_to_nominal','change_dB_vs_nominal','ratio_to_full_physical_RMS'});
writetable(calTable,fullfile(out,['calibration_metrics_' tag '.csv']));
angles=[acosd(min(1,abs(ny.'*nyNom))),acosd(min(1,abs(nu.'*nuNom)))];
diagnostic=table({'fixed';'delay'},[r.fixed.train_weighted_error_pct;r.delay.train_weighted_error_pct], ...
    [r.fixed.exitflag;r.delay.exitflag],[r.fixed.optimization.firstorderopt;r.delay.optimization.firstorderopt], ...
    [r.fixed.iterations;r.delay.iterations],[235;242], ...
    'VariableNames',{'stage','training_weighted_error_pct','exitflag','firstorderopt','iterations','active_parameters'});
writetable(diagnostic,fullfile(out,['optimization_' tag '.csv']));
evaluation=struct('tag',tag,'f',f,'fit_metrics',fitTable,'calibration_metrics',calTable, ...
    'calibration_names',{calNames},'calibrated_check_FRFs',{calH},'ny',ny,'nu',nu, ...
    'nominal_ny',nyNom,'nominal_nu',nuNom,'null_angles_deg',angles,'nominal_logical_output',Py, ...
    'nominal_logical_input',Pu,'added_output_delay_s',addedY,'added_input_delay_s',addedU, ...
    'calibrator_scope','Offline check FRF; pure additional causal delays; each port unit norm; nominal logical directions fixed', ...
    'max_compensator_magnitude',1);
evaluation.band_static_ny=nyBand;evaluation.band_static_nu=nuBand;
evaluation.band_static_training_eigenvalues_y=sort(dy);evaluation.band_static_training_eigenvalues_u=sort(du);
save(fullfile(out,['evaluation_' tag '.mat']),'evaluation','-v7.3');
writematrix([Py;ny.'],fullfile(out,['Ty_unit_' tag '.csv']));
writematrix([Pu,nu],fullfile(out,['Tu_unit_' tag '.csv']));
fprintf('PREDICTION %s 25-800 Hz\n',tag);disp(fitTable(fitTable.fmin_Hz==25 & fitTable.fmax_Hz==800,1:5));
fprintf('CALIBRATION %s 25-80 Hz (unit norm, same measured FRF)\n',tag);
disp(calTable(calTable.fmin_Hz==25 & calTable.fmax_Hz==80,[1:4,7:8]));
fprintf('Null angles output/input [%.5g %.5g] deg\n',angles);
fprintf('Updated total delay ms\n');disp(tau*1e3);
end

function h=transform(H,Ty,Tu)
h=zeros(size(H,1),size(Ty,1)*size(Tu,2),'like',H);
for k=1:size(H,1),m=Ty*reshape(H(k,:),4,4)*Tu;h(k,:)=m(:).';end
end
function B=selectBlock(H,index)
switch index
    case 2,cols=reshape(1:16,4,4);cols=cols(1:3,1:3);B=H(:,cols(:));
    case 3,B=H(:,16);
    case 4,B=H(:,13:15);
    case 5,B=H(:,[4 8 12]);
end
end
function v=rmsIntegral(H,f)
v=sqrt(trapz(f,sum(abs(H).^2,2))/max(f(end)-f(1),eps));
end
