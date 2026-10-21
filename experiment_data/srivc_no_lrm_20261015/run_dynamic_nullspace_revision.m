function result=run_dynamic_nullspace_revision(positionIndex,poleHz,richBasis,laguerreOrder,alignmentDelay,crossGuards)
% Direct cross-response calibration with anchored stable triangular maps.
% poleHz=0: constant correction; positive: second-order stable proper bases.
if nargin<1,positionIndex=2;end
if nargin<2,poleHz=40;end
if nargin<3,richBasis=false;end
if nargin<4,laguerreOrder=0;end
if nargin<5,alignmentDelay=0;end
if nargin<6,crossGuards=false;end
assert(alignmentDelay>=0,'Only nonnegative causal alignment delays allowed.');
assert(~(richBasis && laguerreOrder>0),'Use one basis extension at a time.');
base=fileparts(mfilename('fullpath'));out=fullfile(base,'dynamic_nullspace_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'band_calibration_20261016','band_calibration_campaign.mat'),'campaign');
c=S.campaign.records{positionIndex};f=c.f(:);H=reshape(c.H_training.',4,4,[]);
P=c.nominal.Ty./vecnorm(c.nominal.Ty,2,2);Q=c.nominal.Tu./vecnorm(c.nominal.Tu,2,1);
ny=c.nominal.ny(:);nu=c.nominal.nu(:);Ey=[P;ny.'];Eu=[Q,nu];
assert(norm(P*ny)<1e-10 && norm(Q.'*nu)<1e-10,'Nominal anchor is not orthogonal to logical directions.');
Z=complex(zeros(size(H)));for k=1:numel(f),Z(:,:,k)=Ey*H(:,:,k)*Eu;end
config=struct('pole_hz',poleHz,'rich_basis',richBasis,'laguerre_order',laguerreOrder,'training_band',[25 80],'protected_bands',[260 310;320 380;450 800], ...
    'max_port_gain',1.05,'target_relative_rms_budget',.10,'ridge',1e-6,'max_sweeps',8);
config.alignment_delay_seconds=alignmentDelay;
config.cross_guard_bands=zeros(0,2);if crossGuards,config.cross_guard_bands=[3 25;80 170];end
config.cross_guard_rms_ratio=1;
tag=sprintf('%s_p%g',c.tag,poleHz);if richBasis,tag=[tag '_r4'];end
if laguerreOrder>0,tag=[tag '_l' num2str(laguerreOrder)];end
if alignmentDelay>0,tag=[tag '_d' num2str(1e3*alignmentDelay) 'ms'];end
if crossGuards,tag=[tag '_g'];end
phi=basis(f,poleHz,richBasis,laguerreOrder);nb=size(phi,2);np=3*nb;
train=mask(f,config.training_band);ft=f(train);zt=Z(:,:,train);pt=phi(train,:);
w=quadrature(ft);bvec=reshape(zt(1:3,4,:),3,[]).';cvec=reshape(zt(4,1:3,:),3,[]).';
at=exp(-2i*pi*ft*alignmentDelay);bvec=bvec.*at;cvec=cvec.*at;
my=designCross(zt,pt,'y');mu=designCross(zt,pt,'u');
sy=sqrt(repmat(w,3,1)/sum(w.*sum(abs(cvec).^2,2)));
su=sqrt(repmat(w,3,1)/sum(w.*sum(abs(bvec).^2,2)));
My=my.*sy;vy=cvec(:).*sy;Mu=mu.*su;vu=bvec(:).*su;
hy=real(My'*My)+config.ridge*eye(np);gy=real(My'*vy);cy=real(vy'*vy);
hu=real(Mu'*Mu)+config.ridge*eye(np);gu=real(Mu'*vu);cu=real(vu'*vu);
guards=unique([0;logspace(-1,4.4,500).';f(train);config.protected_bands(:)]);
pg=basis(guards,poleHz,richBasis,laguerreOrder);gainY=designDirections(pg,P);gainU=designDirections(pg,Q.');
crossY={};crossU={};
for guardIndex=1:size(config.cross_guard_bands,1)
    ix=mask(f,config.cross_guard_bands(guardIndex,:));fb=f(ix);zb=Z(:,:,ix);pp=phi(ix,:);wb=quadrature(fb);
    by=reshape(zb(4,1:3,:),3,[]).';bu=reshape(zb(1:3,4,:),3,[]).';
    syb=sqrt(repmat(wb,3,1)/sum(wb.*sum(abs(by).^2,2)));sub=sqrt(repmat(wb,3,1)/sum(wb.*sum(abs(bu).^2,2)));
    by=by.*exp(-2i*pi*fb*alignmentDelay);bu=bu.*exp(-2i*pi*fb*alignmentDelay);
    crossY{guardIndex}=struct('M',designCross(zb,pp,'y').*syb,'v',by(:).*syb,'budget',config.cross_guard_rms_ratio); %#ok<AGROW>
    crossU{guardIndex}=struct('M',designCross(zb,pp,'u').*sub,'v',bu(:).*sub,'budget',config.cross_guard_rms_ratio); %#ok<AGROW>
end
protected=cell(1,size(config.protected_bands,1));
for b=1:numel(protected)
    ix=mask(f,config.protected_bands(b,:));fb=f(ix);zb=Z(:,:,ix);pp=phi(ix,:);wb=quadrature(fb);
    dd=reshape(zb(4,4,:),[],1);scale=sqrt(wb/sum(wb.*abs(dd).^2));
    protected{b}=struct('Z',zb,'phi',pp,'scale',scale,'reference',dd,'alignment',exp(-2i*pi*fb*alignmentDelay));
end
opts=optimoptions('fmincon','Algorithm','sqp','Display','none','SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true,'MaxIterations',250,'MaxFunctionEvaluations',1200, ...
    'OptimalityTolerance',1e-9,'ConstraintTolerance',1e-9,'StepTolerance',1e-11);
xy=zeros(np,1);xu=zeros(np,1);history={};previous=inf;
if (richBasis || laguerreOrder>0) && alignmentDelay==0 && ~crossGuards
    prior=load(fullfile(out,sprintf('%s_p%g.mat',c.tag,poleHz)),'result');
    xy(1:6)=prior.result.xy;xu(1:6)=prior.result.xu;
end
% Check derivatives at a small, nonzero configuration, including protection.
testY=.005*sin((1:np).');testU=.006*cos((1:np).');
dataY=constraintData('y',testU);dataU=constraintData('u',testY);
checks.y=derivativeCheck(testY,hy,gy,cy,dataY);checks.u=derivativeCheck(testU,hu,gu,cu,dataU);
assert(checks.y<1e-5 && checks.u<1e-5,'Analytic derivative check failed.');
for sweep=1:config.max_sweeps
    dataY=constraintData('y',xu);
    [xy,jy,ey,oy]=fmincon(@(x)objective(x,hy,gy,cy),xy,[],[],[],[],[],[],@(x)constraints(x,dataY),opts);
    dataU=constraintData('u',xy);
    [xu,ju,eu,ou]=fmincon(@(x)objective(x,hu,gu,cu),xu,[],[],[],[],[],[],@(x)constraints(x,dataU),opts);
    val=jy+ju;cvY=constraints(xy,constraintData('y',xu));cvU=constraints(xu,constraintData('u',xy));
    violation=max([0;cvY;cvU]);
    history(end+1,:)={sweep,jy,ju,val,ey,eu,oy.firstorderopt,ou.firstorderopt,violation}; %#ok<AGROW>
    fprintf('%s sweep %d J %.6g y/u %.6g %.6g flags %d %d violation %.3g\n',tag,sweep,val,jy,ju,ey,eu,violation);
    assert(val<=previous+1e-6,'Alternating objective increased.');
    if abs(previous-val)<1e-8 && violation<1e-7,break;end
    previous=val;
end
result=struct('tag',tag,'position_mm',c.position_mm,'config',config,'xy',xy,'xu',xu, ...
    'nominal',c.nominal,'gradient_checks',checks,'history',{history},'f',f,'H_training',c.H_training, ...
    'H_validation',c.H_validation,'source',c.source_data,'matlab_version',version, ...
    'method','Anchored stable triangular cross-response calibration; not physical channel-delay identification');
result.training_objective=[objective(xy,hy,gy,cy),objective(xu,hu,gu,cu)];
result.feasibility=max([0;constraints(xy,constraintData('y',xu));constraints(xu,constraintData('u',xy))]);
result.side_exitflags=[ey eu];result.side_optimality=[oy.firstorderopt ou.firstorderopt];
result.unconstrained_xy=-(hy\gy);result.unconstrained_xu=-(hu\gu);
result.unconstrained_training_objective=[objective(result.unconstrained_xy,hy,gy,cy),objective(result.unconstrained_xu,hu,gu,cu)];
result.real_design_condition=[cond([real(My);imag(My)]),cond([real(Mu);imag(Mu)])];
result.training_constraint_y=constraints(xy,constraintData('y',xu));
result.training_constraint_u=constraints(xu,constraintData('u',xy));
result.protected_training_relative_error=sqrt(max(0,result.training_constraint_y(end-2:end)+config.target_relative_rms_budget^2));
result.filters=makeFilters(poleHz,richBasis,laguerreOrder,xy,xu,P,Q,ny,nu);
result.filters.alignment_delay_seconds=alignmentDelay;
% Exact orthogonal-anchor identity, valid also with a pure-delay anchor.
result.filters.output_hinf=sqrt(1+norm(ss(result.filters.fy*P),inf)^2);
result.filters.input_hinf=sqrt(1+norm(ss(Q*result.filters.fu),inf)^2);
result.stable=all(real(pole(ss(result.filters.fy)))<0) && all(real(pole(ss(result.filters.fu)))<0);
if alignmentDelay>0
    % Keep causal branches separately: a*anchor + stable rational correction.
    % No delay inversion is implemented; do not claim a causal inverse map.
    result.filters=rmfield(result.filters,{'Ny','Nu'});
    result.filters.output_anchor=ny.';result.filters.input_anchor=nu;
    result.filters.output_correction=result.filters.fy*P;
    result.filters.input_correction=Q*result.filters.fu;
end
result.accepted=result.feasibility<1e-7 && all(result.side_exitflags>0) && result.stable ...
    && result.filters.output_hinf<=config.max_port_gain+1e-5 && result.filters.input_hinf<=config.max_port_gain+1e-5;
save(fullfile(out,[tag '.mat']),'result','-v7.3');
fprintf('SAVED %s accepted=%d Hinf y/u %.7f %.7f gradient %.3g %.3g\n',tag,result.accepted, ...
    result.filters.output_hinf,result.filters.input_hinf,checks.y,checks.u);

    function d=constraintData(side,other)
        if side=='y',d.gain=gainY;d.anchor=ny;else,d.gain=gainU;d.anchor=nu;end
        d.anchor_rows=exp(-2i*pi*guards*alignmentDelay)*d.anchor.';
        d.gamma=config.max_port_gain;
        if side=='y',d.target=crossY;else,d.target=crossU;end
        ncross=numel(d.target);
        for bi=1:numel(protected)
            pr=protected{bi};zz=pr.Z;otherF=pr.phi*reshape(other,3,nb).';nf=size(zz,3);
            constant=complex(zeros(nf,1));v=complex(zeros(nf,3));
            for kk=1:nf
                a=zz(1:3,1:3,kk);bb=zz(1:3,4,kk);cc=zz(4,1:3,kk);
                if side=='y'
                    u=otherF(kk,:).';constant(kk)=pr.alignment(kk)*cc*u;v(kk,:)=(pr.alignment(kk)*bb+a*u).';
                else
                    y=otherF(kk,:);constant(kk)=pr.alignment(kk)*y*bb;v(kk,:)=pr.alignment(kk)*cc+y*a;
                end
            end
            mt=complex(zeros(nf,np));
            for ll=1:nb,for aa=1:3,mt(:,aa+3*(ll-1))=pr.phi(:,ll).*v(:,aa);end,end
            d.target{ncross+bi}=struct('M',mt.*pr.scale,'v',constant.*pr.scale,'budget',config.target_relative_rms_budget);
        end
    end
end

function phi=basis(f,poleHz,richBasis,laguerreOrder)
if poleHz==0,phi=ones(numel(f),1);else,q=1i*f(:)/poleHz;phi=[1./(1+q).^2,q./(1+q).^2];end
if richBasis
    extra=phi(:,1)./(1+1i*f(:)/20);phi=[phi,extra,extra./(1+1i*f(:)/80)];
end
if laguerreOrder>0
    q=1i*f(:)/poleHz;lag0=1./(1+q);alpha=(q-1)./(q+1);
    for ell=2:laguerreOrder-1,phi(:,end+1)=lag0.*alpha.^ell;end %#ok<AGROW>
end
end
function M=designCross(Z,phi,side)
nf=size(Z,3);nb=size(phi,2);M=complex(zeros(nf*3,3*nb));
for l=1:nb
    for a=1:3
        if side=='y',v=reshape(Z(a,1:3,:),3,[]).';else,v=reshape(Z(1:3,a,:),3,[]).';end
        v=v.*phi(:,l);M(:,a+3*(l-1))=v(:);
    end
end
end
function M=designDirections(phi,B)
nf=size(phi,1);nb=size(phi,2);M=complex(zeros(nf,4,3*nb));
for l=1:nb,for a=1:3,M(:,:,a+3*(l-1))=phi(:,l)*B(a,:);end,end
end
function [v,g]=objective(x,H,h,c)
v=x.'*H*x+2*h.'*x+c;if nargout>1,g=2*(H*x+h);end
end
function [c,ceq,gc,gceq]=constraints(x,d)
ng=size(d.gain,1);np=numel(x);G=reshape(d.gain,ng*4,np);
z=reshape(G*x,ng,4)+d.anchor_rows;c=sum(abs(z).^2,2)-d.gamma^2;
gc=zeros(np,ng+numel(d.target));
for a=1:np,gc(a,1:ng)=2*real(sum(conj(z).*d.gain(:,:,a),2)).';end
for b=1:numel(d.target)
    t=d.target{b};r=t.M*x+t.v;c(end+1,1)=real(r'*r)-t.budget^2; %#ok<AGROW>
    gc(:,ng+b)=2*real(t.M'*r);
end
ceq=[];gceq=[];
end
function errorValue=derivativeCheck(x,H,h,c,d)
[~,g]=objective(x,H,h,c);[cv,~,gc]=constraints(x,d);gn=zeros(size(x));cn=zeros(numel(x),numel(cv));
for a=1:numel(x)
    step=1e-6;xp=x;xm=x;xp(a)=xp(a)+step;xm(a)=xm(a)-step;
    gn(a)=(objective(xp,H,h,c)-objective(xm,H,h,c))/(2*step);
    cn(a,:)=(constraints(xp,d)-constraints(xm,d)).'/(2*step);
end
errorValue=max(norm(g-gn)/max(norm(g),1e-8),norm(gc-cn,'fro')/max(norm(gc,'fro'),1e-8));
end
function F=makeFilters(poleHz,richBasis,laguerreOrder,xy,xu,P,Q,ny,nu)
if poleHz==0
    F.fy=tf(xy.');F.fu=tf(xu);F.Ny=tf(ny.'+xy.'*P);F.Nu=tf(nu+Q*xu);
elseif laguerreOrder>0
    p=2*pi*poleHz;L=laguerreOrder;Ag=-p*eye(L)+p*diag(ones(L-1,1),-1);Bg=[p;zeros(L-1,1)];Cg=zeros(L);
    Cg(1,2)=1;Cg(2,1:2)=[1 -1];
    for ell=2:L-1,for rr=0:ell,Cg(ell+1,rr+1)=nchoosek(ell,rr)*(-2)^rr;end,end
    yy=reshape(xy,3,L)*Cg;uu=reshape(xu,3,L)*Cg;
    F.fy=ss(kron(eye(3),Ag),kron(eye(3),Bg),reshape(yy.',1,[]),zeros(1,3));
    F.fu=ss(Ag,Bg,uu,zeros(3,1));
    F.Ny=ss(ny.')+F.fy*P;F.Nu=ss(nu)+Q*F.fu;
else
    p=2*pi*poleHz;den=[1 2*p p^2];bt={tf(p^2,den),tf([p 0],den)};
    if richBasis,bt{3}=bt{1}*tf(2*pi*20,[1 2*pi*20]);bt{4}=bt{3}*tf(2*pi*80,[1 2*pi*80]);end
    F.fy=tf(zeros(1,3));F.fu=tf(zeros(3,1));
    for l=1:numel(bt),F.fy=F.fy+bt{l}*xy(3*(l-1)+(1:3)).';F.fu=F.fu+bt{l}*xu(3*(l-1)+(1:3));end
    F.Ny=tf(ny.')+F.fy*P;F.Nu=tf(nu)+Q*F.fu;
end
end
function ix=mask(f,band)
ix=f>=band(1)-1e-6 & f<=band(2)+1e-6;
end
function w=quadrature(f)
df=diff(f);w=[df(1);df(1:end-1)+df(2:end);df(end)]/2;
end
