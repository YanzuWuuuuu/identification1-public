function result=run_information_audit()
% Local diagnostics only. No refitting or modification of prior experiment data.
here=fileparts(mfilename('fullpath')); parent=fullfile(fileparts(here),'closed_loop_validation');addpath(parent);
S=load(fullfile(parent,'outputs','data_distinct.mat'));
Q=load(fullfile(parent,'outputs_structured','distinct_N32','selected.mat'),'selected');
cfg=S.cfg;cfg.additive_first_coefficient=true; f=S.data.f;nf=numel(f);
[xt,ct]=clv_chart(S.p,cfg);[Gt,Jgt]=clv_sampled(xt,ct,f,cfg.Ts);
Gp=reshape(Gt,4,4,nf);Id=repmat(eye(4),1,1,nf);
% S1: the identical fitted mean and fixed fitting weights, only a coordinate shear.
fit=Q.selected;c=fit.chart;Kf=freqresp(S.K,2*pi*f);
[~,Jfit]=clv_mean(fit.x,c,f,cfg.Ts,Kf,S.data.D);
J=[real(Jfit);imag(Jfit)]./fit.sigma/sqrt(2*numel(S.data.mean));
T=eye(c.n);T(c.to+(1:7),c.ao+(1:7))=-(c.as/c.ts)*eye(7);
% x_delay = eta - (as/ts)*x_a; x_a and x_b are retained, same 125 dimensions.
Jnew=J*T; C=relative_jacobian(fit.x,c,[20;100;286.43;800]);
old=target_bound(J,C);new=target_bound(Jnew,C*T);
assert(abs(old.target_trace/new.target_trace-1)<1e-7,'Target covariance changed under invertible coordinates.');
[~,svmat,V]=svd(J,0);weak=cell(3,1);
for k=1:3
 v=V(:,end-k+1);da=c.as*c.a_map*v(c.ao+(1:7));dl=c.ts*c.a_map*v(c.to+(1:7));
 db=fit.estimate.b(:).*v(c.bo+(1:16));
 e=[norm(v(1:c.to))^2,norm(v(c.to+(1:7)))^2,norm(v(c.ao+(1:7)))^2,norm(v(c.bo+(1:16)))^2];
 weak{k}=struct('singular_value',svmat(end-k+1,end-k+1),'mechanical_delay_a_b_energy_fraction',e, ...
  'slope_cancellation_ratio',norm(da+dl)/max(norm(da)+norm(dl),realmin), ...
  'curvature_cancellation_ratio',norm(db-fit.estimate.a(:).*da)/max(norm(db)+norm(fit.estimate.a(:).*da),realmin), ...
  'a_vs_negative_delay_cosine',-da'*dl/max(norm(da)*norm(dl),realmin), ...
  'delta_a_seconds',da.','delta_ell_seconds',dl.','delta_b_seconds_squared',db.','right_singular_vector',v.'); %#ok<AGROW>
end
weak=vertcat(weak{:});
result.coordinates=struct('old_condition',cond(J),'eta_condition',cond(Jnew), ...
 'old_unit_column_condition',cond(J./vecnorm(J)), 'eta_unit_column_condition',cond(Jnew./vecnorm(Jnew)), ...
 'target_trace_old',old.target_trace,'target_trace_eta',new.target_trace, ...
 'target_trace_relative_difference',abs(old.target_trace/new.target_trace-1),'weak_directions',weak);
writejson(fullfile(here,'coordinates.json'),result.coordinates);
% Correct proportional model tangent at the real-numerator truth: imaginary
% residues with real rank-one left/right factors, poles still unknown.
realids=[1:ct.ro,ct.to+1:ct.n];
for i=1:ct.K,realids=[realids,ct.ro+14*(i-1)+[1:3,11:14]];end %#ok<AGROW>
realids=sort(realids);assert(numel(realids)==90);
assert(norm(real(S.p.R(:)))/norm(S.p.R(:))<1e-12,'Proportional tangent requires imaginary residues at truth.');
Ctr=relative_jacobian(xt,ct,[20;100;286.43;800]);
delayC=zeros(6,ct.n);delayC(1:3,ct.to+(1:3))=ct.ts*eye(3);
for j=2:4,delayC(2+j,ct.to+3+j)=ct.ts;delayC(2+j,ct.to+4)=-ct.ts;end
[Un,Sn,Vn]=svd(S.nom.M);ctrlMap=Vn(:,1:3)*diag(1./diag(Sn(1:3,1:3)))*Un(:,1:3)';
[Ap,Bp,Cp,Dp]=ssdata(S.plant);assert(norm(Dp,'fro')<1e-14);
sigmaReal=cfg.sensor_noise_std*sqrt(2/cfg.period_samples);
JopenComplex=reshape(pagemtimes(reshape(Jgt,4,4,nf,ct.n),reshape(S.data.D,4,4,nf,1)),16*nf,ct.n);
Jopen=[real(JopenComplex);imag(JopenComplex)]/sigmaReal;
result.open_loop_reference=information(Jopen,ct,Ctr,delayC,32);
result.open_loop_reference.note='Only a frequency-information reference; no open-loop time experiment on the double integrator.';
designs=[15 30 60 100 150 200];scan=cell(numel(designs),1);
fcross=logspace(log10(2),log10(260),500);Gcross=freqresp(S.plant,2*pi*fcross);
for ib=1:numel(designs)
 design=designs(ib);wc=2*pi*design;wd=2*pi*120;
 kc=c2d(ss(tf([2*.8*wc+wc^2/wd,wc^2],[1/wd,1])),cfg.Ts,'tustin');
 K=minreal(ctrlMap*kc,1e-10,false);[Ak,Bk,Ck,Dk]=ssdata(K);
 ev=eig([Ap-Bp*Dk*Cp,-Bp*Ck;Bk*Cp,Ak]);rho=max(abs(ev));
 row=struct('controller_design_hz',design,'closed_radius',rho,'stable',rho<1-1e-8, ...
 'first_downcross_hz',[],'rigid_excitation_20_100_relative',NaN,'null_excitation_20_100_relative',NaN, ...
 'rigid_excitation_20hz_relative',NaN,'rigid_excitation_40hz_relative',NaN,'total_applied_input_energy_relative',NaN, ...
 'full_covariance_cancellation_relative_error',NaN,'complex',struct(),'proportional',struct(), ...
 'shaped_d_energy_factor',NaN,'shaped_d_peak_gain',NaN,'shaped_mean_information_relative_error',NaN);
 Kcf=freqresp(K,2*pi*fcross);loop=pagemtimes(Kcf,Gcross);vals=zeros(3,numel(fcross));
 for k=1:numel(fcross),sv=svd(loop(:,:,k));vals(:,k)=sv(1:3);end
 crossings=nan(1,3);for j=1:3,ix=find(vals(j,1:end-1)>=1 & vals(j,2:end)<1,1);if ~isempty(ix),crossings(j)=fcross(ix);end;end
 row.first_downcross_hz=crossings;
 if row.stable
  Kf=freqresp(K,2*pi*f);Su=pagemldivide(Id+pagemtimes(Kf,Gp),Id);Sy=pagemldivide(Id+pagemtimes(Gp,Kf),Id);
  SuD=pagemtimes(Su,S.data.D);
  JW=reshape(pagemtimes(reshape(Jgt,4,4,nf,ct.n),reshape(SuD,4,4,nf,1)),16*nf,ct.n);
  JR=[real(JW);imag(JW)]/sigmaReal;
  row.complex=information(JR,ct,Ctr,delayC,32);
  row.proportional=information(JR(:,realids),restrict_chart(ct,realids),Ctr(:,realids),delayC(:,realids),32);
  % Explicit check: whiten the full Sy Sigma Sy* covariance at one frequency.
  kk=10;Jloc=reshape(Jgt((kk-1)*16+(1:16),:),4,4,ct.n);
  Jclosed=reshape(pagemtimes(pagemtimes(Sy(:,:,kk),Jloc),SuD(:,:,kk)),16,ct.n);
  Sigma=Sy(:,:,kk)*Sy(:,:,kk)';L=chol(Sigma,'lower');
  Jwh=reshape(pagemldivide(L,reshape(Jclosed,4,4,ct.n)),16,ct.n);
  Jideal=JW((kk-1)*16+(1:16),:);
  row.full_covariance_cancellation_relative_error=norm(real(Jwh'*Jwh)-real(Jideal'*Jideal),'fro')/norm(real(Jideal'*Jideal),'fro');
  assert(row.full_covariance_cancellation_relative_error<1e-10);
  mask=f<=100;Ur=Vn(:,1:3);nu=Vn(:,4);rp=zeros(nnz(mask),1);np=rp;
  inds=find(mask);for j=1:numel(inds),k=inds(j);rp(j)=norm(Ur'*SuD(:,:,k),'fro')^2;np(j)=norm(nu'*SuD(:,:,k),'fro')^2;end
  row.rigid_excitation_20_100_relative=sum(rp)/(nnz(mask)*3*cfg.force_peak_per_line^2);
  row.null_excitation_20_100_relative=sum(np)/(nnz(mask)*cfg.force_peak_per_line^2);
  row.rigid_excitation_20hz_relative=norm(Ur'*SuD(:,:,find(f==20,1)),'fro')^2/(3*cfg.force_peak_per_line^2);
  row.rigid_excitation_40hz_relative=norm(Ur'*SuD(:,:,find(f==40,1)),'fro')^2/(3*cfg.force_peak_per_line^2);
  row.total_applied_input_energy_relative=norm(SuD(:))^2/norm(S.data.D(:))^2;
  Dshape=pagemtimes(Id+pagemtimes(Kf,Gp),S.data.D);Ushape=pagemtimes(Su,Dshape);
  row.shaped_d_energy_factor=norm(Dshape(:))^2/norm(S.data.D(:))^2;
  row.shaped_d_peak_gain=max(abs(Dshape),[],'all')/cfg.force_peak_per_line;
  Jshape=reshape(pagemtimes(reshape(Jgt,4,4,nf,ct.n),reshape(Ushape,4,4,nf,1)),16*nf,ct.n);
  row.shaped_mean_information_relative_error=norm(Jshape-JopenComplex,'fro')/norm(JopenComplex,'fro');
  assert(row.shaped_mean_information_relative_error<1e-10);
 end
 scan{ib}=row;writejson(fullfile(here,'bandwidth_scan.json'),vertcat(scan{1:ib}));
 fprintf('design %g Hz: rho %.7f stable %d\n',design,rho,row.stable);
end
scan=vertcat(scan{:});result.scan=scan;result.scope='Local mean information at true model, stationary output-white-noise spectral approximation, no covariance-derivative information; finite-record diagonal estimator and empirical coverage are not validated by this calculation.';
save(fullfile(here,'information_audit.mat'),'result','J','Jnew','T','Ctr','Jopen','xt','ct','weak','-v7.3');
writejson(fullfile(here,'information_audit.json'),result);
end

function c=restrict_chart(c,ids)
c.delay_ids=find(ismember(ids,c.to+(1:7)));c.n=numel(ids);
end
function out=information(J,c,C,D,N)
if isfield(c,'delay_ids'),ids=c.delay_ids;else,ids=c.to+(1:7);end
nu=setdiff(1:size(J,2),ids);A=J(:,nu);A=A./vecnorm(A);
[U,S,~]=svd(A,0);s=diag(S);U=U(:,s>max(s)*1e-12);
assert(size(U,2)==numel(nu),'Numerically unresolved nuisance rank; do not report truncated information.');
B=J(:,ids);R=B-U*(U'*B);ss=svd(R,0);dd=target_bound(J,D);
[~,sp,vp]=svd(R,0);profileTrace=norm((D(:,ids)*vp)./diag(sp).','fro')^2;
assert(abs(profileTrace/dd.target_trace-1)<1e-5,'Profile information and full target covariance disagree.');
tt=target_bound(J,C);ds=svd(J./vecnorm(J),0);
out=struct('profile_delay_min_eigenvalue',min(ss)^2,'profile_delay_trace',norm(R,'fro')^2, ...
 'delay_residual_energy_fraction',norm(R,'fro')^2/norm(B,'fro')^2, ...
 'delay_target_rms_bound_us',sqrt(dd.target_trace/6/N)*1e6, ...
 'relative_response_rms_bound_percent',sqrt(tt.target_trace/size(C,1)/N)*100, ...
 'unit_column_condition',ds(1)/ds(end),'nuisance_rank',size(U,2),'nuisance_dimension',numel(nu), ...
 'minimum_relative_unit_column_singular_value',ds(end)/ds(1));
end
function out=target_bound(J,C)
scale=vecnorm(J);[~,S,V]=svd(J./scale,0);s=diag(S);
out.target_trace=norm(((C./scale)*V)./s.','fro')^2;
out.condition=s(1)/s(end);
end
function C=relative_jacobian(x,c,f)
p=clv_decode(x,c);s=2i*pi*f(:);C=complex(zeros(6*numel(f),c.n));
for k=1:numel(f)
 d=1+p.a(:)*s(k)+p.b(:)*s(k)^2;
 Jd=zeros(16,c.n);Jd(:,c.ao+(1:7))=s(k)*c.as*c.a_map;Jd(:,c.bo+(1:16))=diag(s(k)^2*p.b(:));
 for j=2:4
  row=6*(k-1)+j-1;C(row,:)=Jd(1,:)/d(1)-Jd(j,:)/d(j);C(row,c.to+j-1)=C(row,c.to+j-1)-s(k)*c.ts;
  row=6*(k-1)+3+j-1;e=1+4*(j-1);C(row,:)=Jd(1,:)/d(1)-Jd(e,:)/d(e);
  C(row,c.to+3+j)=C(row,c.to+3+j)-s(k)*c.ts;C(row,c.to+4)=C(row,c.to+4)+s(k)*c.ts;
 end
end
end
function writejson(path,data)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(data,'PrettyPrint',true));
end
