function [initial,info]=init_five_from_frf(G,f,options)
%INIT_FIVE_FROM_FRF Separate data-derived initializer for a five-mode comparison.
% No true poles, modal directions, rigid factors, or geometry are consumed.
% Public priors: number of modes, search band, damping range, and a common
% initial channel time constant. G is 4-by-4-by-N or 16-by-N; f is in Hz.
% Optional sigma is the complex-observation standard deviation, not a
% weight recomputed from the noisy magnitude. This function does no LRM.
if nargin<3, options=struct; end
if ~isfield(options,'num_modes'), options.num_modes=5; end
if ~isfield(options,'frequency_band'), options.frequency_band=[200 800]; end
if ~isfield(options,'low_frequency_max'), options.low_frequency_max=min(40,options.frequency_band(1)/5); end
if ~isfield(options,'initial_time_constant'), options.initial_time_constant=20e-6; end
if ~isfield(options,'damping_range'), options.damping_range=[.001 .08]; end
K=options.num_modes;
assert(isscalar(K) && K>=1 && K==round(K),'num_modes must be a positive integer.');
f=f(:); G=reshape(G,16,[]);
assert(numel(f)==size(G,2) && all(f>0),'Positive frequency grid and G must agree.');
assert(all(isfinite(G(:))) && all(isfinite(f)),'Finite observations are required.');
[f,order]=sort(f); G=G(:,order); s=1i*2*pi*f.'; N=numel(f);
if isfield(options,'sigma') && ~isempty(options.sigma)
    sigma=reshape(options.sigma,16,[]); sigma=sigma(:,order);
    assert(isequal(size(sigma),size(G)) && all(sigma(:)>0),'Invalid sigma.');
else
    sigma=ones(size(G));
end
T=options.initial_time_constant;
% All channels start with the same stable two-factor denominator. Removing
% that initial denominator only prepares a start; the measured G is retained
% unchanged for the subsequent likelihood fit.
den0=(1+T*s).^2; P0=G.*den0; sigmaP=sigma.*abs(den0);

% Near zero: s^2 G = M + A_1 s + A_2 s^2 + ... . Real-coefficient
% extrapolation accounts for the first-order channel phase lag; taking one
% low-frequency sample as M would bias the static null spaces.
low=find(f<=options.low_frequency_max);
if numel(low)<6, low=1:min(N,max(6,ceil(.08*N))); end
assert(numel(low)>=4,'Too few low-frequency observations.');
scale=2*pi*max(f(low)); q=s(low)/scale;
X=[ones(numel(low),1),q.',(q.^2).'];
Mvec=zeros(16,1); lag1=zeros(16,1);
for k=1:16
    z=(G(k,low).*s(low).^2).';
    sig=(sigma(k,low).*abs(s(low)).^2).';
    coef=realWls(X,z,sig);
    Mvec(k)=coef(1); lag1(k)=coef(2)/scale;
end
Mraw=reshape(Mvec,4,4); [U,S,V]=svd(Mraw); S(4,4)=0;
initial.M=U*S*V';

% Multiplication by omega^2 removes the gross displacement-frequency scale.
% Prominence, rather than the largest bins, distinguishes K separate modes.
residual=P0-initial.M(:)*(1./s.^2);
amplitude=sqrt(sum(abs(residual).^2,1)).*abs(s).^2;
band=find(f>=options.frequency_band(1) & f<=options.frequency_band(2));
assert(numel(band)>=20,'Insufficient points in the modal search band.');
y=log(max(amplitude,realmin)); candidates=[]; prominence=[];
for ii=2:numel(band)-1
    k=band(ii);
    if y(k)>=y(band(ii-1)) && y(k)>y(band(ii+1))
        width=max(18,.09*f(k));
        left=band(f(band)>=f(k)-width & f(band)<f(k));
        right=band(f(band)<=f(k)+width & f(band)>f(k));
        if isempty(left) || isempty(right), continue; end
        p=min(y(k)-min(y(left)),y(k)-min(y(right)));
        candidates(end+1)=k; prominence(end+1)=p; %#ok<AGROW>
    end
end
[~,ranked]=sort(prominence,'descend'); selected=[];
for it=ranked
    k=candidates(it);
    if isempty(selected) || all(abs(f(k)-f(selected))>max(12,.035*f(k)))
        selected(end+1)=k; %#ok<AGROW>
        if numel(selected)==options.num_modes, break; end
    end
end
assert(numel(selected)==options.num_modes, ...
    'Only %d separated modal peaks found; inspect band, grid, and data.',numel(selected));
[~,ord]=sort(f(selected)); selected=selected(ord);
lambda=zeros(K,1); localInfo=struct([]);
for i=1:K
    fp=f(selected(i)); width=max(22,.09*fp);
    if i>1, width=min(width,.43*(fp-f(selected(i-1)))); end
    if i<K, width=min(width,.43*(f(selected(i+1))-fp)); end
    local=band(abs(f(band)-fp)<=width);
    if numel(local)<10
        [~,near]=sort(abs(f(band)-fp)); local=sort(band(near(1:min(14,numel(band)))));
    end
    step=median(diff(f(local)));
    freqLimits=[max(options.frequency_band(1),fp-max(3*step,.04*fp)), ...
        min(options.frequency_band(2),fp+max(3*step,.04*fp))];
    objective=@(x) localCost(x,s(local),P0(:,local),sigmaP(:,local), ...
        freqLimits,options.damping_range);
    [xp,value,flag]=fminsearch(objective,log([fp;.012]), ...
        optimset('Display','off','MaxIter',180,'MaxFunEvals',420, ...
        'TolX',1e-9,'TolFun',1e-10));
    fn=exp(xp(1)); zeta=exp(xp(2)); omega=2*pi*fn;
    lambda(i)=-zeta*omega+1i*omega*sqrt(1-zeta^2);
    localInfo(i).peak_hz=fp;
    localInfo(i).frequency_hz=fn;
    localInfo(i).damping=zeta;
    localInfo(i).objective=value;
    localInfo(i).exitflag=flag;
    localInfo(i).frequency_limits=freqLimits;
    localInfo(i).points=numel(local);
end
initial.lambda=lambda;

% With poles fixed, residues and a real constant D are linear parameters.
% Keep the low-frequency rank-three estimate fixed in this initial solve.
X=ones(N,1);
for i=1:K
    a=1./(s-lambda(i)); b=1./(s-conj(lambda(i)));
    X=[X,(a+b).',(1i*(a-b)).']; %#ok<AGROW>
end
coef=zeros(size(X,2),16);
for k=1:16
    coef(:,k)=realWls(X,(P0(k,:)-initial.M(k)./s.^2).',sigmaP(k,:).');
end
initial.D=reshape(coef(1,:),4,4); initial.R=complex(zeros(4,4,K));
for i=1:K
    R=reshape(coef(2*i,:)+1i*coef(2*i+1,:),4,4);
    [U,S,V]=svd(R); initial.R(:,:,i)=U(:,1)*S(1,1)*V(:,1)';
end
initial.a=2*T*ones(4); initial.b=T^2*ones(4);
initial.ty=zeros(4,1); initial.tu=zeros(4,1);
info.description='Low-frequency real polynomial, spectral prominence, local shared-pole variable projection, global residue WLS.';
info.uses_true_parameters=false; info.num_modes=K;
info.frequency_band=options.frequency_band;
info.low_frequency_range=[min(f(low)),max(f(low))];
info.low_frequency_points=numel(low);
info.rigid_before_rank_projection=Mraw;
info.rigid_rank_projection_relative_change=norm(Mraw-initial.M,'fro')/max(norm(Mraw,'fro'),eps);
info.low_frequency_linear_coefficient=reshape(lag1,4,4);
info.candidate_peak_hz=f(candidates);
info.candidate_prominence=prominence(:);
info.selected_peak_hz=f(selected);
info.local_modal_fit=localInfo;
info.initial_time_constant=T;
info.initial_delay_seconds=0;
info.spectral_frequency_hz=f(band);
info.spectral_indicator=amplitude(band).';
end

function cost=localCost(x,s,Y,sigma,frequencyLimits,dampingLimits)
fn=exp(x(1)); zeta=exp(x(2));
if ~isfinite(fn) || ~isfinite(zeta) || fn<frequencyLimits(1) || fn>frequencyLimits(2) ...
        || zeta<dampingLimits(1) || zeta>dampingLimits(2)
    cost=1e12+1e6*sum(x.^2); return;
end
w=2*pi*fn; lam=-zeta*w+1i*w*sqrt(1-zeta^2);
a=1./(s-lam); b=1./(s-conj(lam));
v=(imag(s)-mean(imag(s)))/max(std(imag(s)),eps);
% A complex linear background approximates the distant modes and smooth
% channel response in this local positive-frequency window.
X=[(a+b).',(1i*(a-b)).',ones(numel(s),1),1i*ones(numel(s),1),v.',1i*v.'];
cost=0;
for k=1:16
    beta=realWls(X,Y(k,:).',sigma(k,:).');
    r=(X*beta-Y(k,:).')./sigma(k,:).';
    cost=cost+sum(abs(r).^2);
end
end

function beta=realWls(X,y,sigma)
sigma=max(sigma,max(sigma)*1e-12);
A=X./sigma; b=y./sigma; A=[real(A);imag(A)]; b=[real(b);imag(b)];
col=sqrt(sum(A.^2,1)); col=max(col,realmin);
beta=(A./col)\b; beta=beta./col.';
end

