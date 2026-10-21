function met=audit_band_metrics(truth,nominal,fit,K,cfg,test)
estimate=fit.estimate;f=test.f(:);nf=numel(f);s=2i*pi*f.';
[G,Gr,Gm]=clv_parts(truth,f,cfg.Ts);[Gh,~,~]=clv_parts(estimate,f,cfg.Ts);
Kf=freqresp(K,2*pi*f);Ip=repmat(eye(4),1,1,nf);
T=pagemldivide(Ip+pagemtimes(reshape(G,4,4,nf),Kf),reshape(G,4,4,nf));
Th=pagemldivide(Ip+pagemtimes(reshape(Gh,4,4,nf),Kf),reshape(Gh,4,4,nf));
met.frequency_relative_error=(abs(estimate.lambda)-abs(truth.lambda))./abs(truth.lambda);
met.damping_relative_error=(-real(estimate.lambda)./abs(estimate.lambda))./(-real(truth.lambda)./abs(truth.lambda))-1;
met.estimated_frequencies_hz=abs(estimate.lambda)/(2*pi);
met.estimated_damping=-real(estimate.lambda)./abs(estimate.lambda);
met.rigid_matrix_relative_error=norm(estimate.M-truth.M,'fro')/norm(truth.M,'fro');
[Ut,~,Vt]=svd(truth.M);[Uh,~,Vh]=svd(estimate.M);[Un,~,Vn]=svd(nominal.M);
nyt=Ut(:,4);nut=Vt(:,4);nyh=Uh(:,4);nuh=Vh(:,4);nyn=Un(:,4);nun=Vn(:,4);
if nyh'*nyt<0,nyh=-nyh;end;if nuh'*nut<0,nuh=-nuh;end
met.output_null_angle_degrees=asind(min(1,norm(nyh*nyh'-nyt*nyt',2)));
met.input_null_angle_degrees=asind(min(1,norm(nuh*nuh'-nut*nut',2)));
met.heldout_plant_relative_error=norm(Gh-G,'fro')/norm(G,'fro');
met.heldout_closed_relative_error=norm(Th(:)-T(:))/norm(T(:));
for z=1:2
    if z==1,pp=truth;else,pp=estimate;end
    d=1+pp.a(:)*s+pp.b(:)*s.^2;d=reshape(d,4,4,nf);
    for k=1:4
        ry(k,:,z)=exp(-(pp.ty(k)-pp.ty(1))*s).*reshape(d(1,1,:)./d(k,1,:),1,nf); %#ok<AGROW>
        ru(k,:,z)=exp(-(pp.tu(k)-pp.tu(1))*s).*reshape(d(1,1,:)./d(1,k,:),1,nf); %#ok<AGROW>
    end
end
met.relative_output_error=max(abs(ry(:,:,2)./ry(:,:,1)-1),[],'all');
met.relative_input_error=max(abs(ru(:,:,2)./ru(:,:,1)-1),[],'all');
met.relative_output_curve=max(abs(ry(:,:,2)./ry(:,:,1)-1),[],1);
met.relative_input_curve=max(abs(ru(:,:,2)./ru(:,:,1)-1),[],1);
met.frequency_hz=f;
[cyh,cuh,cySys,~]=clv_compensators(estimate,f,cfg.Ts);
[cyt,cut,cyTrueSys,~]=clv_compensators(truth,f,cfg.Ts);
names={'nominal','static','joint','oracle'};NY=[nyn,nyh,nyh,nyt];NU=[nun,nuh,nuh,nut];
% Rigid component is propagated under the full object's actual clean u.
rigidTruth=truth;rigidTruth.R(:)=0;[~,~,~,rigidPlant]=clv_system(rigidTruth,nominal,cfg);
N=cfg.period_samples;keep=cfg.burn_samples+(1:N);bins=round(f/cfg.Fs*N)+1;
rigidY=cell(4,1);
for inp=1:4,rigidY{inp}=lsim(rigidPlant,test.clean_u{inp},test.time_full);end
for im=1:4
    cy=ones(4,nf);cu=cy;cs=repmat({ss([],[],[],1,cfg.Ts)},4,1);
    if im==3,cy=cyh;cu=cuh;cs=cySys;elseif im==4,cy=cyt;cu=cut;cs=cyTrueSys;end
    ny=NY(:,im);nu=NU(:,im);ffr=zeros(1,nf);fff=zeros(1,nf);single=complex(zeros(5,nf));
    for k=1:nf
        left=(ny.*cy(:,k)).';right=nu.*cu(:,k);right=right/norm(right);
        ffr(k)=left*reshape(Gr(:,k),4,4)*right;
        for j=1:5,single(j,k)=left*reshape(Gm(:,k,j),4,4)*right;end
        fff(k)=sum(single(:,k));
    end
    item.double_null_rigid_flex_ratio=sqrt(sum(abs(ffr).^2)/sum(abs(fff).^2));
    item.double_null_rigid_curve=abs(ffr);item.double_null_flexible_curve=abs(fff);
    item.single_mode_l2=sqrt(sum(abs(single).^2,2));
    item.max_output_compensator_gain=max(abs(cy),[],'all');
    item.max_input_compensator_gain=max(abs(cu),[],'all');
    Ramps=complex(zeros(4,nf));Famps=Ramps;Noise=complex(zeros(4,nf,test.groups));
    for inp=1:4
        rr=zeros(numel(test.time_full),1);ff=rr;
        for j=1:4
            rr=rr+ny(j)*lsim(cs{j},rigidY{inp}(:,j),test.time_full);
            ff=ff+ny(j)*lsim(cs{j},test.clean_full{inp}(:,j)-rigidY{inp}(:,j),test.time_full);
        end
        ar=(2/N)*fft(rr(keep));af=(2/N)*fft(ff(keep));Ramps(inp,:)=ar(bins).';Famps(inp,:)=af(bins).';
        for r=1:test.groups
            vv=test.raw{inp,r}-test.clean_full{inp};nn=zeros(size(rr));
            for j=1:4,nn=nn+ny(j)*lsim(cs{j},vv(:,j),test.time_full);end
            an=(2/N)*fft(nn(keep));Noise(inp,:,r)=an(bins).';
        end
    end
    item.general_rigid_curve_m=sqrt(mean(abs(Ramps).^2,1)/2);
    item.general_noise_curve_m=sqrt(mean(abs(Noise).^2,[1 3])/2);
    bands={[20 800],[20 100],[200 800],[260 310]};bn={'all','low','modal','first_resonance'};
    for ib=1:4
        mask=f>=bands{ib}(1)&f<=bands{ib}(2);
        v.rigid_rms_m=sqrt(mean(sum(abs(Ramps(:,mask)).^2,2))/2);
        v.flexible_rms_m=sqrt(mean(sum(abs(Famps(:,mask)).^2,2))/2);
        v.noise_rms_m=sqrt(mean(sum(abs(Noise(:,mask,:)).^2,2),'all')/2);
        v.flexible_to_rigid_noise_db=20*log10(v.flexible_rms_m/hypot(v.rigid_rms_m,v.noise_rms_m));
        v.rigid_to_noise_ratio=v.rigid_rms_m/v.noise_rms_m;
        item.output_only.(bn{ib})=v;
    end
    met.calibration.(names{im})=item;
end
oracle=met.calibration.oracle.single_mode_l2;
for im=1:4,met.calibration.(names{im}).mode_retention_vs_digital_oracle=met.calibration.(names{im}).single_mode_l2./oracle;end
met.scope='Heldout frequencies and fresh sample-domain noise. Output-only calibration leaves feedback and input identical. Double-null maps use sampled plant, causal digital compensators and unit applied input per frequency; they are not closed-loop scalar transfer functions.';
met.digital_implementation='Bilinear rational ratios plus integer delay and first-order Thiran remainder; 400 microsecond common delay on each side. Digital oracle includes discretization/hold residual, not exact continuous cancellation.';
end
