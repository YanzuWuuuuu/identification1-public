function fit=frf_srivc_fit(f,H,frequency0,damping0,bounds,W,opts)
% Additive SRIVC on complex FRFs, with shared real denominator coefficients.
% W is the nonnegative QUADRATIC cost weight: the regression uses sqrt(W).
% H is Nfrequency-by-Nresponse. Each response has its own real numerator
% coefficients and optional real direct term, but no record-boundary terms.
% Updates retain the measured-regressor IV matrix Z'*X, not Z'*Z/Gauss-Newton.
if nargin<7,opts=struct;end
if ~isfield(opts,'max_iterations'),opts.max_iterations=100;end
if ~isfield(opts,'score_tolerance'),opts.score_tolerance=1e-6;end
if ~isfield(opts,'direct'),opts.direct=true;end
if ~isfield(opts,'verbose'),opts.verbose=true;end
if ~isfield(opts,'debug'),opts.debug=false;end
if ~isfield(opts,'fixed_modes'),opts.fixed_modes=1;end
f=f(:);frequency0=frequency0(:);nm=numel(frequency0);nc=size(H,2);nf=numel(f);
if isscalar(damping0),damping0=repmat(damping0,nm,1);end
if nargin<6 || isempty(W),W=ones(size(H));end
if isscalar(W),W=repmat(W,size(H));end
if isvector(W) && numel(W)==nf,W=repmat(W(:),1,nc);end
assert(size(H,1)==nf && isequal(size(W),size(H)),'H/W must have Nf rows and matching channels.');
assert(isreal(W) && all(isfinite(W(:))) && all(W(:)>=0),'W must be finite nonnegative quadratic weights.');
assert(all(isfinite(H(:))) && all(isfinite(f)) && all(f>=0),'H and f must be finite; f cannot be negative.');
assert(isequal(size(bounds),[nm 2]) && numel(damping0)==nm,'Bounds or damping dimensions do not match modes.');
assert(all(frequency0>0) && all(damping0>0),'Initial frequency and damping must be positive.');
assert(all(ismember(opts.fixed_modes,1:nm)),'fixed_modes must contain valid mode indices.');
sqrtW=sqrt(W);a=[2*damping0(:),ones(nm,1)];referenceHz=frequency0;
fixed=reshape([2*opts.fixed_modes-1;2*opts.fixed_modes],1,[]);free=setdiff(1:2*nm,fixed);
[B,~,~]=design(f,a,referenceHz,opts.direct);nl=size(B,2);
C=zeros(nl,nc);normalizer=norm(stack(H.*sqrtW),'fro');
assert(normalizer>0,'Weighted FRF data must have nonzero energy.');
for c=1:nc
    Br=stack(B.*sqrtW(:,c));target=stack(H(:,c).*sqrtW(:,c));
    sb=max(vecnorm(Br),realmin);[Q,R]=qr(Br./sb,0);
    assert(size(R,1)==nl && size(R,2)==nl,'Insufficient weighted frequency points for linear coefficients.');
    C(:,c)=(R\(Q.'*target))./sb.';
end
history=[];reason='maximum_iterations';converged=false;
previousAlpha=NaN;previousRidge=NaN;debug=struct;
for it=0:opts.max_iterations
    [B,A,q]=design(f,a,referenceHz,opts.direct);
    pred=B*C;res=H-pred;E=stack(res.*sqrtW);cost=norm(E,'fro')^2;
    linear=cell(nc,1);Z=cell(nc,1);X=cell(nc,1);
    moment=zeros(2*nm);rhs=zeros(2*nm,1);znorm=zeros(2*nm,1);
    gap2=0;zenergy=0;linearScore=0;linearConditions=zeros(nc,1);
    for c=1:nc
        Br=stack(B.*sqrtW(:,c));sb=max(vecnorm(Br),realmin);
        [Q,R]=qr(Br./sb,0);linear{c}=struct('B',Br,'Q',Q,'R',R,'scale',sb);
        linearConditions(c)=cond(R);
        z=complex(zeros(nf,2*nm));x=z;
        for m=1:nm
            ix=2*m-1:2*m;part=B(:,ix)*C(ix,c);
            filt=[q(:,m),q(:,m).^2]./A(:,m);
            z(:,ix)=-filt.*part;
            x(:,ix)=-filt.*(part+res(:,c));
        end
        Z{c}=stack(z.*sqrtW(:,c));X{c}=stack(x.*sqrtW(:,c));
        zp=Z{c}-Q*(Q.'*Z{c});xp=X{c}-Q*(Q.'*X{c});
        ep=E(:,c)-Q*(Q.'*E(:,c));
        moment=moment+zp.'*xp;rhs=rhs+zp.'*ep;
        znorm=znorm+sum(zp.^2,1).';
        linearScore=max(linearScore,max(abs((Br./sb).'*E(:,c))));
        gap2=gap2+norm(X{c}-Z{c},'fro')^2;zenergy=zenergy+norm(Z{c},'fro')^2;
    end
    sd=max(sqrt(znorm),realmin);M=moment./(sd*sd.');g=rhs./sd;
    score=max([abs(g(free));linearScore])/max(sqrt(cost),1e-14*normalizer);
    if isempty(free),condition=1;else,condition=cond(M(free,free));end
    history(end+1,:)=[it,sqrt(cost)/normalizer,score,previousAlpha,previousRidge,condition, ...
        sqrt(gap2/max(zenergy,realmin))]; %#ok<AGROW>
    if opts.verbose && (mod(it,5)==0 || score<opts.score_tolerance)
        fprintf('FRF SRIVC %3d: weighted error %.5f%%, score %.3g, cond %.3g, X-Z gap %.3g\n', ...
            it,100*sqrt(cost)/normalizer,score,condition,history(end,7));
    end
    if opts.debug && it==0
        debug.initial_regressor_instrument_gap=history(end,7);
        % Audit exact elimination on the constrained (free-denominator) system.
        nd=numel(free);FF=zeros(nd+nc*nl);rr=zeros(nd+nc*nl,1);
        for c=1:nc
            Br=linear{c}.B;zz=Z{c}(:,free);xx=X{c}(:,free);
            ix=nd+(c-1)*nl+(1:nl);
            FF(1:nd,1:nd)=FF(1:nd,1:nd)+zz.'*xx;
            FF(1:nd,ix)=zz.'*Br;FF(ix,1:nd)=Br.'*xx;FF(ix,ix)=Br.'*Br;
            rr(1:nd)=rr(1:nd)+zz.'*E(:,c);rr(ix)=Br.'*E(:,c);
        end
        ss=sqrt(max(abs(diag(FF)),realmin));FS=FF./(ss*ss.');
        dfull=(FS\(rr./ss))./ss;
        da=zeros(2*nm,1);
        if nd>0,da(free)=(M(free,free)\g(free))./sd(free);end
        dlin=zeros(nl,nc);
        for c=1:nc
            L=linear{c};dlin(:,c)=(L.R\(L.Q.'*(E(:,c)-X{c}*da)))./L.scale.';
        end
        dschur=[da(free);dlin(:)];
        debug.schur_full_update_relative_error=norm(dfull-dschur)/max(norm(dfull),realmin);
        debug.schur_denominator_relative_error=norm(dfull(1:nd)-da(free))/max(norm(da(free)),realmin);
        debug.full_iv_scaled_condition=cond(FS);
        debug.full_iv_relative_equation_residual=norm(FF*dschur-rr)/max(norm(rr),realmin);
        if nd>0,j=free(1);else,j=1;end
        m=ceil(j/2);k=mod(j-1,2)+1;ap=a;am=a;step=1e-6;
        ap(m,k)=ap(m,k)+step;am(m,k)=am(m,k)-step;
        bp=design(f,ap,referenceHz,opts.direct);bm=design(f,am,referenceHz,opts.direct);
        fd=stack(((bp-bm)*C/(2*step)).*sqrtW);zz=zeros(size(fd));
        for c=1:nc,zz(:,c)=Z{c}(:,j);end
        debug.instrument_derivative_relative_error=norm(fd-zz,'fro')/max(norm(zz,'fro'),realmin);
    end
    if score<opts.score_tolerance || sqrt(cost)/normalizer<1e-10
        converged=true;reason='instrument_correlation';break
    end
    if it==opts.max_iterations,break;end
    accepted=false;
    for ridge=[0,1e-8,1e-6,1e-4,1e-2,1,100]
        da=zeros(2*nm,1);
        da(free)=((M(free,free)+ridge*eye(numel(free)))\g(free))./sd(free);
        dc=zeros(size(C));
        for c=1:nc
            L=linear{c};dc(:,c)=(L.R\(L.Q.'*(E(:,c)-X{c}*da)))./L.scale.';
        end
        for backtrack=0:18
            alpha=2^-backtrack;an=a+alpha*reshape(da,2,[]).';cn=C+alpha*dc;
            if ~all(isfinite(an(:))) || any(an(:,2)<=0),continue;end
            hz=referenceHz./sqrt(an(:,2));zeta=an(:,1)./(2*sqrt(an(:,2)));
            stable=all(zeta>=.00003 & zeta<=.95) && all(hz>=bounds(:,1) & hz<=bounds(:,2));
            if ~stable,continue;end
            bn=design(f,an,referenceHz,opts.direct);rn=stack((H-bn*cn).*sqrtW);trial=norm(rn,'fro')^2;
            if isfinite(trial) && trial<cost,accepted=true;break;end
        end
        if accepted,break;end
    end
    if ~accepted,reason='no_descent';break;end
    a=an;C=cn;previousAlpha=alpha;previousRidge=ridge;
end
if opts.direct,directCoefficient=C(end,:);else,directCoefficient=zeros(1,nc);end
fit=struct('method','frequency_domain_additive_SRIVC_FRF', ...
    'a',a,'reference_hz',referenceHz,'frequency_hz',referenceHz./sqrt(a(:,2)), ...
    'damping',a(:,1)./(2*sqrt(a(:,2))),'coefficients',C, ...
    'numerator_coefficients',C(1:2*nm,:),'direct',opts.direct,'direct_coefficients',directCoefficient, ...
    'n_channels',nc,'prediction',pred,'relative_error',sqrt(cost)/normalizer, ...
    'weighted_cost',cost,'quadratic_weights',W,'converged',converged,'reason',reason, ...
    'score',score,'score_tolerance',opts.score_tolerance,'iterations',it,'history',history,'debug',debug, ...
    'initial_frequency_hz',frequency0,'bounds_hz',bounds,'fixed_modes',opts.fixed_modes, ...
    'linear_design_condition',linearConditions,'initialization','fixed-pole linear LS; no nonlinear OE optimization', ...
    'frequency_band_hz',[min(f),max(f)],'f',f, ...
    'coefficient_layout','Rows 2m-1/2m: b0_m/b1_m, with q_m=i*f/reference_hz(m); optional final row: real D');
end

function [B,A,q]=design(f,a,referenceHz,direct)
nm=size(a,1);q=1i*f./referenceHz(:).';A=1+q.*a(:,1).'+q.^2.*a(:,2).';
B=complex(zeros(numel(f),2*nm+double(direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if direct,B(:,end)=1;end
end

function r=stack(z)
r=[real(z);imag(z)];
end
