function result=raw_periodic_validation()
%RAW_PERIODIC_VALIDATION Periodic time samples -> FFT -> joint FRF fit.
% This is an exactly periodic steady-state simulation. There are no startup
% transients, input errors, aliasing, quantization, or noncoherent leakage.
% The four-mode+D generator is isolated from data initialization and fitting.
% Post-calibration tests have independent sample-domain sensor noise.
folder=fileparts(mfilename('fullpath'));
pilot=fullfile(fileparts(folder),'relative_calibration_pilot'); addpath(pilot);
outdir=fullfile(folder,'outputs'); if ~isfolder(outdir), mkdir(outdir); end
source=load(fullfile(pilot,'outputs','common_noise0.mat'),'truth','nominal');
truth=source.truth; nominal=source.nominal;
truth.filter_y=[20;31;9;25]*1e-6; truth.filter_u=[17;11;23;35]*1e-6;
truth.ty=[0;8;17;5]*1e-6; truth.tu=[4;17;10;25]*1e-6;
truth.a=truth.filter_y+truth.filter_u.'; truth.b=truth.filter_y*truth.filter_u.';
Fs=5000; periodSamples=2500; recordCount=8; forcePeakPerLine=.1;
f=[2:2:180,200:2:800].'; bins=round(f/Fs*periodSamples)+1;
assert(all(abs((bins-1)*Fs/periodSamples-f)<1e-12),'Noncoherent grid.');
lineCount=numel(f); complexFrfNoiseStd=1e-8;
% FFT amplitude convention: X_amp = (2/N)*FFT(x) at positive frequencies.
% Real white sample noise sigma_t gives complex amplitude RMS 2*sigma_t/sqrt(N)
% and each-quadrature SD sigma_t*sqrt(2/N). Thus the requested per-record
% complex FRF SD is reached by sigma_t = A*sigma_FRF*sqrt(N)/2.
sensorNoiseStd=forcePeakPerLine*complexFrfNoiseStd*sqrt(periodSamples)/2;
sampleTime=(0:periodSamples-1)/Fs;
[Gtrue,~,~,~]=truth_parts(truth,f);
rng(2026101401,'twister');
phase=2*pi*rand(4,lineCount); inputAmplitude=forcePeakPerLine*exp(1i*phase);
rawFrf=complex(zeros(16,lineCount,recordCount));
exampleTraining=struct;
for input=1:4
    ua=inputAmplitude(input,:);
    u=time_from_amplitudes(ua,bins,periodSamples);
    measuredU=positive_amplitudes(u,bins,periodSamples);
    rows=(1:4)+4*(input-1);
    ya=Gtrue(rows,:).*ua;
    cleanY=time_from_amplitudes(ya,bins,periodSamples);
    for r=1:recordCount
        measuredY=cleanY+sensorNoiseStd*randn(4,periodSamples);
        measuredYa=positive_amplitudes(measuredY,bins,periodSamples);
        rawFrf(rows,:,r)=measuredYa./measuredU;
        if input==1 && r==1
            exampleTraining.time_seconds=sampleTime;
            exampleTraining.input_1_force_N=u;
            exampleTraining.clean_outputs_m=cleanY;
            exampleTraining.measured_outputs_m=measuredY;
        end
    end
end
Y=mean(rawFrf,3);
repeatVariance=sum(abs(rawFrf-Y).^2,3)/(recordCount-1);
repeatVariance=max(repeatVariance,median(repeatVariance(:))*1e-8);
sigmaMean=sqrt(repeatVariance/recordCount);
options=struct('num_modes',4,'frequency_band',[200 800],'sigma',sigmaMean);
tic; [initial,info]=init_from_frf(Y,f,options); initializationSeconds=toc;
fprintf('Periodic raw-data initialization: %.2f %.2f %.2f %.2f Hz\n', ...
    abs(initial.lambda)/(2*pi));
fit=fit_from_frf(Y,f,sigmaMean,initial,'expanded',1);
metrics=rcp_metrics(truth,fit.estimate,nominal,f);

% New experiment: a virtual multisine is distributed through the estimated
% input compensators, the true plant generates four sensor records, and
% fresh sample-domain noise is added before output calibration/projection.
% The same fresh noise realizations are reused across baseline comparisons.
% Compensation is evaluated at exact DFT bins, equivalent to its periodic
% steady-state action; this does not verify digital startup implementation.
rng(2026101402,'twister');
testPhase=2*pi*rand(1,lineCount); commandAmplitude=forcePeakPerLine*exp(1i*testPhase);
testNoiseSamples=sensorNoiseStd*randn(4,periodSamples,recordCount);
testNoiseAmplitudes=complex(zeros(4,lineCount,recordCount));
for r=1:recordCount
    testNoiseAmplitudes(:,:,r)=positive_amplitudes(testNoiseSamples(:,:,r),bins,periodSamples);
end
[G,Gr,Gf,Gd]=truth_parts(truth,f);
Gmode=complex(zeros(16,lineCount,4));
for i=1:4, Gmode(:,:,i)=single_mode_response(truth,i,f); end
[ny0,nu0]=null_pair(truth.M); [nye,nue]=null_pair(fit.estimate.M);
[nyn,nun]=null_pair(nominal.M);
names={'nominal','static_only','joint','oracle'};
ny=[nyn,nye,nye,ny0]; nu=[nun,nue,nue,nu0];
[cyJoint,cuJoint]=compensators(fit.estimate,f);
[cyOracle,cuOracle]=compensators(truth,f);
exampleTest=struct; testSpectra=struct; testResults=struct;
for m=1:4
    cy=ones(4,lineCount); cu=cy;
    if m==3, cy=cyJoint; cu=cuJoint; end
    if m==4, cy=cyOracle; cu=cuOracle; end
    appliedInput=cu.*nu(:,m).*commandAmplitude;
    left=cy.*ny(:,m);
    rigid=complex(zeros(1,lineCount)); flex=rigid; background=rigid;
    rawOutput=complex(zeros(4,lineCount)); singleMode=complex(zeros(4,lineCount));
    for k=1:lineCount
        rigid(k)=left(:,k).'*reshape(Gr(:,k),4,4)*appliedInput(:,k);
        flex(k)=left(:,k).'*reshape(Gf(:,k),4,4)*appliedInput(:,k);
        background(k)=left(:,k).'*reshape(Gd(:,k),4,4)*appliedInput(:,k);
        rawOutput(:,k)=reshape(G(:,k),4,4)*appliedInput(:,k);
        for i=1:4
            singleMode(i,k)=left(:,k).'*reshape(Gmode(:,k,i),4,4)*appliedInput(:,k);
        end
    end
    projectedNoise=reshape(sum(left.*testNoiseAmplitudes,1),lineCount,recordCount).';
    cleanTestY=time_from_amplitudes(rawOutput,bins,periodSamples);
    finalMeasured=complex(zeros(recordCount,lineCount));
    for r=1:recordCount
        actualTestY=cleanTestY+testNoiseSamples(:,:,r);
        actualTestAmplitude=positive_amplitudes(actualTestY,bins,periodSamples);
        finalMeasured(r,:)=sum(left.*actualTestAmplitude,1);
    end
    decomposition=rigid+flex+background+projectedNoise;
    v.all_excited=band_metrics(true(size(f)),rigid,flex,background,projectedNoise);
    v.low_frequency=band_metrics(f<=180,rigid,flex,background,projectedNoise);
    v.flexible_band=band_metrics(f>=200,rigid,flex,background,projectedNoise);
    v.max_output_compensation_gain=max(abs(cy),[],'all');
    v.max_input_compensation_gain=max(abs(cu),[],'all');
    v.command_virtual_rms_N=sqrt(sum(abs(commandAmplitude).^2)/2);
    v.applied_input_rms_N=sqrt(sum(abs(appliedInput).^2,2)/2);
    v.sample_fft_decomposition_relative_error=norm(finalMeasured-decomposition,'fro')/max(norm(finalMeasured,'fro'),realmin);
    for i=1:4
        fc=abs(truth.lambda(i))/(2*pi);
        modeMask=abs(f-fc)<=max(12,.04*fc);
        v.mode_window(i).frequency_hz=fc;
        v.mode_window(i).band_hz=[min(f(modeMask)),max(f(modeMask))];
        v.mode_window(i).metrics=band_metrics(modeMask,rigid,flex,background,projectedNoise);
        v.mode_window(i).single_mode_rms_m=sqrt(sum(abs(singleMode(i,modeMask)).^2)/2);
        v.mode_window(i).other_modes_rms_m=sqrt(sum(abs(flex(modeMask)-singleMode(i,modeMask)).^2)/2);
        v.mode_window(i).single_mode_to_noise_db=20*log10(max(v.mode_window(i).single_mode_rms_m,realmin)/ ...
            max(v.mode_window(i).metrics.noise_rms_m,realmin));
    end
    testResults.(names{m})=v;
    testSpectra.(names{m})=struct('rigid_amplitude_m',rigid,'flexible_amplitude_m',flex, ...
        'D_background_amplitude_m',background,'noise_amplitude_m',projectedNoise, ...
        'measured_projected_amplitude_m',finalMeasured,'single_mode_amplitude_m',singleMode);
    if m==3
        exampleTest.time_seconds=sampleTime;
        exampleTest.applied_input_force_N=time_from_amplitudes(appliedInput,bins,periodSamples);
        exampleTest.clean_sensor_outputs_m=cleanTestY;
        exampleTest.measured_sensor_outputs_m=exampleTest.clean_sensor_outputs_m+testNoiseSamples(:,:,1);
        exampleTest.projected_rigid_m=time_from_amplitudes(rigid,bins,periodSamples);
        exampleTest.projected_flexible_m=time_from_amplitudes(flex,bins,periodSamples);
        exampleTest.projected_D_m=time_from_amplitudes(background,bins,periodSamples);
        exampleTest.projected_band_noise_m=time_from_amplitudes(projectedNoise(1,:),bins,periodSamples);
        exampleTest.projected_measured_excited_band_m=time_from_amplitudes(finalMeasured(1,:),bins,periodSamples);
    end
end
result.description='Periodic sample-domain sensor-noise simulation, full 4x4 FRF, data-derived initialization, one expanded-model fit.';
result.scope='Exact periodic steady state. No transient, noncoherent leakage, input noise, aliasing, quantization, or digital filter startup is tested.';
result.generation_model='Four physical modal pairs plus a real constant D, with distinct first-order channel filters and pure delays.';
result.sampling_Hz=Fs; result.samples_per_period=periodSamples;
result.frequency_spacing_Hz=Fs/periodSamples;
result.frequency_lines=lineCount; result.frequency_intervals_Hz=[2 180;200 800];
result.input_channel_count=4; result.periods_per_input=recordCount;
result.training_duration_seconds=4*recordCount*periodSamples/Fs;
result.independent_test_periods=recordCount;
result.force_peak_per_frequency_line_N=forcePeakPerLine;
result.sensor_white_noise_sample_std_m=sensorNoiseStd;
result.designed_per_record_complex_frf_std_m_per_N=complexFrfNoiseStd;
result.realized_mean_per_record_complex_frf_std_m_per_N=sqrt(mean(repeatVariance,'all'));
result.estimated_mean_frf_std_m_per_N=sqrt(mean(sigmaMean.^2,'all'));
result.data_mean_truth_relative_error=norm(Y-Gtrue,'fro')/norm(Gtrue,'fro');
result.initialization_seconds=initializationSeconds;
result.initial_frequency_hz=abs(initial.lambda)/(2*pi);
result.frequency_relative_error=(abs(fit.estimate.lambda)/(2*pi)-truth.freq(:))./truth.freq(:);
result.damping_relative_error=(-real(fit.estimate.lambda)./abs(fit.estimate.lambda)-.01)/.01;
result.solver=rmfield(fit,{'estimate','x','x0','c','fitted','lb','ub'});
result.identification_metrics=rmfield(metrics,'curves');
result.independent_test=testResults;
result.noise_metric_note='Noise RMS is measured from fresh time-domain white noise after the output compensator and projection, restricted to excited frequency lines; it is not inferred from training fit residuals.';
result.statistical_certificate=false;
tag='raw_periodic_distinct_expanded';
save(fullfile(outdir,[tag '.mat']),'result','truth','nominal','fit','initial','info', ...
    'f','Y','sigmaMean','rawFrf','metrics','exampleTraining','exampleTest','testSpectra');
fid=fopen(fullfile(outdir,[tag '.json']),'w','n','UTF-8'); assert(fid>=0);
clean=onCleanup(@() fclose(fid)); fprintf(fid,'%s',jsonencode(result,'PrettyPrint',true));
fprintf('Raw periodic validation: flag %d, relative dynamic error %.3g%%, flexible-band leakage/flex %.3g, SNR %.2f dB\n', ...
    fit.exitflag,100*max(metrics.relative.output_max_error,metrics.relative.input_max_error), ...
    testResults.joint.flexible_band.rigid_to_flexible_rms_ratio, ...
    testResults.joint.flexible_band.flexible_to_noise_db);
end

function y=time_from_amplitudes(a,bins,N)
a=reshape(a,[],numel(bins)); spectrum=complex(zeros(size(a,1),N));
spectrum(:,bins)=(N/2)*a;
spectrum(:,N-bins+2)=(N/2)*conj(a);
y=real(ifft(spectrum,[],2));
end

function a=positive_amplitudes(y,bins,N)
spectrum=fft(y,[],2); a=(2/N)*spectrum(:,bins);
end

function [G,Gr,Gf,Gd]=truth_parts(p,f)
s=2i*pi*f(:).'; N=numel(f);
rig=p.M(:)*(1./s.^2); bg=p.D(:)*ones(1,N);
flex=complex(zeros(16,N));
for i=1:4
    R=p.R(:,:,i); flex=flex+R(:)*(1./(s-p.lambda(i)))+conj(R(:))*(1./(s-conj(p.lambda(i))));
end
ell=p.ty+p.tu.'; H=exp(-ell(:)*s)./(1+p.a(:)*s+p.b(:)*s.^2);
Gr=rig.*H; Gf=flex.*H; Gd=bg.*H; G=Gr+Gf+Gd;
end

function G=single_mode_response(p,i,f)
s=2i*pi*f(:).'; R=p.R(:,:,i);
P=R(:)*(1./(s-p.lambda(i)))+conj(R(:))*(1./(s-conj(p.lambda(i))));
ell=p.ty+p.tu.'; H=exp(-ell(:)*s)./(1+p.a(:)*s+p.b(:)*s.^2);
G=P.*H;
end

function [ny,nu]=null_pair(M)
[U,~,V]=svd(M); ny=U(:,end); nu=V(:,end);
end

function [cy,cu]=compensators(p,f)
s=2i*pi*f(:).'; dy=p.ty-p.ty(1); du=p.tu-p.tu(1);
Ly=max([0;dy])+10e-6; Lu=max([0;du])+10e-6;
cy=complex(zeros(4,numel(f))); cu=cy;
for k=1:numel(f)
    d=1+p.a*s(k)+p.b*s(k)^2;
    cy(:,k)=exp(-s(k)*(Ly-dy)).*d(:,1)/d(1,1);
    cu(:,k)=exp(-s(k)*(Lu-du)).*d(1,:).'/d(1,1);
end
end

function out=band_metrics(mask,rigid,flex,background,noise)
mask=mask(:).';
out.frequency_line_count=sum(mask);
out.rigid_rms_m=sqrt(sum(abs(rigid(mask)).^2)/2);
out.flexible_rms_m=sqrt(sum(abs(flex(mask)).^2)/2);
out.D_background_rms_m=sqrt(sum(abs(background(mask)).^2)/2);
out.noise_rms_m=sqrt(mean(sum(abs(noise(:,mask)).^2,2)/2));
out.noise_rms_per_period_m=sqrt(sum(abs(noise(:,mask)).^2,2)/2);
out.coherently_averaged_noise_rms_m=sqrt(sum(abs(mean(noise(:,mask),1)).^2)/2);
out.rigid_to_flexible_rms_ratio=out.rigid_rms_m/max(out.flexible_rms_m,realmin);
out.flexible_to_noise_db=20*log10(max(out.flexible_rms_m,realmin)/max(out.noise_rms_m,realmin));
out.flexible_to_rigid_and_noise_db=20*log10(max(out.flexible_rms_m,realmin)/ ...
    max(hypot(out.rigid_rms_m,out.noise_rms_m),realmin));
out.flexible_to_averaged_noise_db=20*log10(max(out.flexible_rms_m,realmin)/ ...
    max(out.coherently_averaged_noise_rms_m,realmin));
end
