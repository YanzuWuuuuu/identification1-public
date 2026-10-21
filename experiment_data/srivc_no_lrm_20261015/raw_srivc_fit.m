function fit=raw_srivc_fit(f,U,Y,frequency0,damping0,bounds,opts)
% Additive frequency-domain SRIVC/RIV on raw input/output DFT coefficients.
% X uses measured partial outputs, Z uses model-generated partial outputs.
% The linear block is eliminated by QR (exact Schur complement, not a GN fit).
if nargin<7,opts=struct;end
if ~isfield(opts,'max_iterations'),opts.max_iterations=100;end
if ~isfield(opts,'score_tolerance'),opts.score_tolerance=1e-6;end
if ~isfield(opts,'direct'),opts.direct=true;end
if ~isfield(opts,'verbose'),opts.verbose=true;end
if ~isfield(opts,'weight'),opts.weight=(f(:)/100).^2;end
if ~isfield(opts,'debug'),opts.debug=false;end
if ~isfield(opts,'fixed_modes'),opts.fixed_modes=[];end
f=f(:);frequency0=frequency0(:);nm=numel(frequency0);ny=size(Y,2);
if isscalar(damping0),damping0=repmat(damping0,nm,1);end
w=opts.weight(:);a=[2*damping0(:),ones(nm,1)];referenceHz=frequency0;
fixed=reshape([2*opts.fixed_modes-1;2*opts.fixed_modes],1,[]);
free=setdiff(1:2*nm,fixed);
[B,~,~,~,~,~]=raw_srivc_design(f,U,a,referenceHz,opts.direct);
Br=stack(B.*w);Yr=stack(Y.*w);sb=max(vecnorm(Br),realmin);
[Q,R]=qr(Br./sb,0);C=(R\(Q.'*Yr))./sb.';
normalizer=norm(Yr,'fro');history=[];reason='maximum_iterations';
converged=false;previousAlpha=NaN;previousRidge=NaN;debug=struct;
for it=0:opts.max_iterations
    [B,A,q,blocks,plantCols,boundaryCols]=raw_srivc_design(f,U,a,referenceHz,opts.direct);
    pred=B*C;res=Y-pred;Br=stack(B.*w);E=stack(res.*w);
    cost=norm(E,'fro')^2;sb=max(vecnorm(Br),realmin);
    [Q,R]=qr(Br./sb,0);
    Z=cell(ny,1);X=Z;Zp=Z;Xp=Z;
    moment=zeros(2*nm);rhs=zeros(2*nm,1);znorm=zeros(2*nm,1);
    rawMomentGradient=zeros(2*nm,1);gap2=0;zenergy=0;
    for y=1:ny
        z=complex(zeros(numel(f),2*nm));x=z;
        for m=1:nm
            part=B(:,blocks{m})*C(blocks{m},y);
            filt=[q(:,m),q(:,m).^2]./A(:,m);
            z(:,2*m-1:2*m)=-filt.*part;
            x(:,2*m-1:2*m)=-filt.*(part+res(:,y));
        end
        Z{y}=stack(z.*w);X{y}=stack(x.*w);
        Zp{y}=Z{y}-Q*(Q.'*Z{y});Xp{y}=X{y}-Q*(Q.'*X{y});
        ep=E(:,y)-Q*(Q.'*E(:,y));
        moment=moment+Zp{y}.'*Xp{y};rhs=rhs+Zp{y}.'*ep;
        znorm=znorm+sum(Zp{y}.^2,1).';
        rawMomentGradient=rawMomentGradient+Z{y}.'*E(:,y);
        gap2=gap2+norm(X{y}-Z{y},'fro')^2;zenergy=zenergy+norm(Z{y},'fro')^2;
    end
    sd=max(sqrt(znorm),realmin);
    M=moment./(sd*sd.');g=rhs./sd;
    % Reduced score also checks the full linear correlation equations.
    linearScore=max(abs((Br./sb).'*E),[],'all');
    score=max([abs(g(free));linearScore])/max(sqrt(cost),1e-14*normalizer);
    history(end+1,:)=[it,sqrt(cost)/normalizer,score,previousAlpha,previousRidge,cond(M(free,free)),sqrt(gap2/max(zenergy,realmin))]; %#ok<AGROW>
    if opts.verbose && (mod(it,5)==0 || score<opts.score_tolerance)
        fprintf('SRIVC %3d: error %.5f%%, score %.3g, cond %.3g, X-Z gap %.3g\n',it,100*sqrt(cost)/normalizer,score,cond(M),history(end,7));
    end
    if opts.debug && it==0
        % Verify Schur elimination against the full IV moment system.
        FF=zeros(2*nm+ny*size(Br,2));rr=zeros(size(FF,1),1);
        for y=1:ny
            ix=2*nm+(y-1)*size(Br,2)+(1:size(Br,2));
            FF(1:2*nm,1:2*nm)=FF(1:2*nm,1:2*nm)+Z{y}.'*X{y};
            FF(1:2*nm,ix)=Z{y}.'*Br;FF(ix,1:2*nm)=Br.'*X{y};FF(ix,ix)=Br.'*Br;
            rr(1:2*nm)=rr(1:2*nm)+Z{y}.'*E(:,y);rr(ix)=Br.'*E(:,y);
        end
        ss=sqrt(max(abs(diag(FF)),realmin));dd=((FF./(ss*ss.'))\(rr./ss))./ss;
        ds=(M\g)./sd;
        debug.schur_denominator_relative_error=norm(dd(1:2*nm)-ds)/max(norm(ds),realmin);
        debug.initial_regressor_instrument_gap=history(end,7);
        % Central difference independently checks a model-generated instrument.
        ap=a;am=a;h=1e-6;ap(1,1)=ap(1,1)+h;am(1,1)=am(1,1)-h;
        bp=raw_srivc_design(f,U,ap,referenceHz,opts.direct);bm=raw_srivc_design(f,U,am,referenceHz,opts.direct);
        fd=stack(((bp-bm)*C/(2*h)).*w);
        zz=zeros(size(fd));for y=1:ny,zz(:,y)=Z{y}(:,1);end
        debug.instrument_derivative_relative_error=norm(fd-zz,'fro')/norm(zz,'fro');
    end
    if score<opts.score_tolerance || sqrt(cost)/normalizer<1e-10
        converged=true;reason='instrument_correlation';break
    end
    if it==opts.max_iterations,break;end
    accepted=false;
    % Stabilized IV solve. Every trial retains Z'X; no Z'Z/GN substitution.
    for ridge=[0,1e-8,1e-6,1e-4,1e-2,1,100]
        da=zeros(2*nm,1);
        da(free)=((M(free,free)+ridge*eye(numel(free)))\g(free))./sd(free);
        dc=zeros(size(C));
        for y=1:ny
            dc(:,y)=(R\(Q.'*(E(:,y)-X{y}*da)))./sb.';
        end
        for backtrack=0:18
            alpha=2^-backtrack;an=a+alpha*reshape(da,2,[]).';cn=C+alpha*dc;
            hz=referenceHz./sqrt(an(:,2));zeta=an(:,1)./(2*sqrt(an(:,2)));
            stable=all(isfinite(an(:)))&&all(an(:,2)>0)&&all(zeta>=.00003 & zeta<=.95) ...
                &&all(hz>=bounds(:,1) & hz<=bounds(:,2));
            if ~stable,continue;end
            bn=raw_srivc_design(f,U,an,referenceHz,opts.direct);rn=stack((Y-bn*cn).*w);
            cnorm=norm(rn,'fro')^2;
            if isfinite(cnorm) && cnorm<cost
                accepted=true;break
            end
        end
        if accepted,break;end
    end
    if ~accepted,reason='no_descent';break;end
    a=an;C=cn;previousAlpha=alpha;previousRidge=ridge;
end
fit=struct('method','frequency_domain_additive_SRIVC_raw_UY', ...
    'a',a,'reference_hz',referenceHz,'frequency_hz',referenceHz./sqrt(a(:,2)), ...
    'damping',a(:,1)./(2*sqrt(a(:,2))),'coefficients',C,'direct',opts.direct, ...
    'n_inputs',size(U,2),'n_outputs',ny,'plant_columns',plantCols,'boundary_columns',boundaryCols, ...
    'prediction',pred,'relative_error',sqrt(cost)/normalizer,'converged',converged, ...
    'reason',reason,'score',score,'iterations',it,'history',history,'debug',debug, ...
    'initial_frequency_hz',frequency0,'bounds_hz',bounds, ...
    'fixed_modes',opts.fixed_modes, ...
    'linear_design_condition',cond(R),'initialization','fixed-pole linear LS; no nonlinear OE optimization', ...
    'frequency_band_hz',[min(f),max(f)]);
end

function r=stack(z)
r=[real(z);imag(z)];
end
