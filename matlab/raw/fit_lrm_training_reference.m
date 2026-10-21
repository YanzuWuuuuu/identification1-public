function fit=fit_lrm_training_reference(u,y,ts)
% Existing LRM + stages 1-4, retrained on the same raw training subset.
% Entire comparison runs in memory; original LRM and stage MAT files untouched.
cfg=identification_config; K=cfg.K; no=size(y,2); ni=size(u,2); n=no*ni;
N=size(u,1); f=(0:floor(N/2))'/(N*ts);
U=fft(u,[],1)/N; Y=fft(y,[],1)/N;
U=U(1:numel(f),:).'; Y=Y(1:numel(f),:).'; G=zeros(n,numel(f));
for output=1:no
    [~,~,~,P]=lrm_miso_core(Y(output,:).',U,9,1);
    G(output+(0:ni-1)*no,:)=P;
end
den=cell(1,K);
for k=1:K
    w=2*pi*cfg.frequency_prior_hz(k); den{k}=[1/w^2,2*cfg.damping_prior/w];
end
[data,grid]=select_frequency_data(G,cfg.initial_frequency_hz,f);
[~,numerators]=SRIVC_init(grid,data,den,frf_weights(data,cfg.weight_floor_relative),[],n,[],[],K);
[data,grid]=select_frequency_data(G,cfg.rigid_frequency_hz,f);
a=den{1}; M=numerators(1:n); B=numerators(n+(1:n));
[~,rigid,rigidInfo]=SRIVC1(grid,data,frf_weights(data,cfg.weight_floor_relative), ...
    [a(2);a(1);M;a(2)*M;a(1)*M+B],n,cfg.tolerance,cfg.max_iterations,1);
M=rigid(3:n+2); first=[rigid(1:2);rigid(2*n+3:3*n+2)-rigid(2)*M];
valid=f>0; f=f(valid); G=G(:,valid);
flexible=G-M./((1i*2*pi*f.').^2);
[data,grid]=select_frequency_data(flexible,cfg.flexible_frequency_hz,f);
initial=zeros((n+2)*K,K);
for k=1:K
    ix=(k-1)*(n+2)+(1:n+2); a=den{k};
    initial(ix,k)=[a(2);a(1);numerators(n+(k-1)*n+(1:n))];
end
initial(1:n+2,1)=first; weight=frf_weights(data,cfg.weight_floor_relative);
[~,estimate,rivInfo]=SRIVC2(grid,data,weight,initial,cfg.tolerance,cfg.max_iterations,K,n,cfg.frequency_bounds_hz);
packed=zeros(K*(n+2),1);
for k=1:K, ix=(k-1)*(n+2)+(1:n+2); packed(ix)=estimate(ix,k); end
[packed,refinement]=refine_flexible_fit(packed,grid,data,weight,K,cfg.frequency_bounds_hz,cfg.max_iterations,cfg.tolerance);
beta=[M;packed];
[metricData,metricGrid]=select_frequency_data(G,unique([cfg.rigid_frequency_hz;cfg.flexible_frequency_hz]),f);
[~,J]=additive_frf(beta,metricGrid,n,K);
Jw=apply_sqrt_weight(J,frf_weights(metricData,cfg.weight_floor_relative));
scales=max(vecnorm(Jw),realmin); [~,R]=qr(Jw./scales,0); metric=R.*scales;
rho0=modal_initialization(beta,no,ni,K,cfg.nrbm);
[~,rho,projection]=Gauss_Newton(beta,rho0,.01,10,200,cfg.tolerance,no,ni,K,cfg.nrbm,metric);
modal=f_rho_construct(rho,no,ni,K,cfg.nrbm);
frequency=zeros(K,1); damping=frequency; residues=zeros(no,ni,K+1);
residues(:,:,1)=reshape(modal(1:n),no,ni);
for k=1:K
    b=modal(n+(k-1)*(n+2)+(1:n+2));
    frequency(k)=1/(2*pi*sqrt(b(2))); damping(k)=b(1)/(2*sqrt(b(2)));
    residues(:,:,k+1)=reshape(b(3:end)/b(2),no,ni);
end
fit=struct('frequency_hz',frequency,'damping',damping,'residues',residues, ...
    'sample_time',ts,'hold_method','foh','method','LRM_retrained_same_training_subset', ...
    'rigid_info',rigidInfo,'riv_info',rivInfo,'refinement',refinement,'projection',projection, ...
    'training_samples',N,'training_frf',G,'training_frequency_hz',f);
end
