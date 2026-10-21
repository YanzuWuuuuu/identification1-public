function fit=fit_residual_d_from_frf(Y,f,sigma,initial,method,includeD,maxIterations)
%FIT_RESIDUAL_D_FROM_FRF Joint fit with an optional static residual matrix.
% Copies the chart, bounds, weights and single-start solver of
% fit_general_from_frf. If includeD=false, D is exactly zero and its 16
% coordinates are omitted from the optimization and objective Jacobian.
% fit.x/x0/lb/ub use the full METHOD coordinate vector. Free-coordinate
% diagnostics and the explicit index mapping are also returned.
if nargin<5 || isempty(method), method='expanded'; end
if nargin<6 || isempty(includeD), includeD=true; end
if nargin<7 || isempty(maxIterations), maxIterations=400; end
assert(any(strcmp(method,{'expanded','physical'})),'Unknown method.');
assert(isscalar(includeD) && any(includeD==[0,1]),'includeD must be a logical scalar.');
includeD=logical(includeD);
assert(isequal(size(Y),size(sigma)) && all(sigma(:)>0),'Invalid observation weights.');
[c,x]=chart(initial); [lb,ub]=domain(x,c);
K=c.num_modes; base=c.filter_offset;
if strcmp(method,'physical')
    x=[x(1:base);zeros(8,1)];
    lb=[lb(1:base);log(.1e-6/20e-6)*ones(8,1)];
    ub=[ub(1:base);log(150e-6/20e-6)*ones(8,1)];
end
x=min(max(x,lb+1e-9),ub-1e-9);
if includeD
    activeidx=(1:numel(x)).';
else
    x(16:31)=0;
    activeidx=[(1:15),(32:numel(x))].';
end
template=x; xfree=x(activeidx); lbf=lb(activeidx); ubf=ub(activeidx);
s=2i*pi*f(:).';
fun=@(z)free_objective(z,template,activeidx,c,s,Y,sigma,method);
[rr,J0]=fun(xfree);
direction=sin((1:numel(xfree)).'); direction=direction/norm(direction);
hh=1e-5; fd=(fun(xfree+hh*direction)-fun(xfree-hh*direction))/(2*hh);
je=norm(fd-J0*direction)/max(norm(fd),realmin);
assert(je<1e-4,'Residual-D model Jacobian check failed: %.3g.',je);
opts=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off', ...
    'MaxIterations',maxIterations,'MaxFunctionEvaluations',4*maxIterations, ...
    'FunctionTolerance',1e-10,'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-6,'ScaleProblem','jacobian');
tic; [zh,cost,r,flag,out,~,J]=lsqnonlin(fun,xfree,lbf,ubf,opts); sec=toc;
xh=template; xh(activeidx)=zh;
xx=expand(xh,c,method); [G,~,estimate]=rcv_model_general(xx,c,s);
if ~includeD
    assert(all(xh(16:31)==0) && all(estimate.D(:)==0),'Fixed-zero D changed unexpectedly.');
end
fit.estimate=estimate; fit.method=method; fit.x=xh; fit.x0=x; fit.c=c;
fit.cost=cost; fit.mean_whitened_residual=cost/numel(Y); fit.exitflag=flag;
fit.iterations=out.iterations; fit.elapsed_seconds=sec;
fit.first_order_optimality=out.firstorderopt;
atBound=find(abs(zh-lbf)<1e-7 | abs(zh-ubf)<1e-7);
fit.active_bounds=activeidx(atBound).'; % Full method-coordinate indices.
fit.active_bounds_full=fit.active_bounds;
fit.active_bounds_free=atBound.';      % Indices into x_free only.
fit.normalized_gradient=max(abs(J'*r)./(vecnorm(J).'*max(norm(r),realmin)));
fit.gradient_free=J'*r;
fit.gradient_coordinate_system='free optimization coordinates';
fit.jacobian_check=je; fit.chart_rcond=rcond(reshape(xh(1:9),3,3));
fit.start_index=1; fit.initial_cost=norm(rr)^2;
fit.iteration_budget=maxIterations; fit.num_modes=K;
fit.fitted=G; fit.lb=lb; fit.ub=ub;
fit.include_D=includeD;
fit.active_parameter_indices=activeidx.';
fit.fixed_parameter_indices=setdiff(1:numel(xh),activeidx.');
fit.num_free_parameters=numel(activeidx);
fit.num_full_method_parameters=numel(xh);
fit.x_free=zh; fit.x0_free=xfree; fit.lb_free=lbf; fit.ub_free=ubf;
fit.x_expanded=xx;
end

function [r,J]=free_objective(z,template,activeidx,c,s,Y,sigma,method)
x=template; x(activeidx)=z;
if nargout>1
    [r,Jfull]=objective(x,c,s,Y,sigma,method);
    J=Jfull(:,activeidx);
else
    r=objective(x,c,s,Y,sigma,method);
end
end

function [r,J]=objective(x,c,s,Y,sigma,method)
if nargout>1
    [xx,T]=expand(x,c,method); [G,Jc]=rcv_model_general(xx,c,s);
    if strcmp(method,'physical')
        b=c.filter_offset; Jc=[Jc(:,1:b),Jc(:,b+(1:32))*T];
    end
    Jc=Jc./sigma(:); J=[real(Jc);imag(Jc)];
else
    xx=expand(x,c,method); G=rcv_model_general(xx,c,s);
end
v=(G-Y)./sigma; r=[real(v(:));imag(v(:))];
end

function [xx,T]=expand(x,c,method)
if ~strcmp(method,'physical'), xx=x; T=[]; return; end
base=c.filter_offset;
ty=20e-6*exp(x(base+(1:4))); tu=20e-6*exp(x(base+(5:8)));
a=ty+tu.'; b=ty*tu.';
xx=[x(1:base);log(a(:)/c.as);log(b(:)/c.bs)];
if nargout>1
    T=zeros(32,8);
    for j=1:4
        for k=1:4
            idx=j+4*(k-1);
            T(idx,j)=ty(j)/a(j,k); T(idx,4+k)=tu(k)/a(j,k);
            T(16+idx,j)=1; T(16+idx,4+k)=1;
        end
    end
end
end

function [c,x]=chart(p)
K=numel(p.lambda); assert(K>=1 && size(p.R,3)==K,'Inconsistent modal initializer.');
c.num_modes=K; c.delay_offset=31+16*K; c.filter_offset=c.delay_offset+7;
[c.U,~,c.V]=svd(p.M); c.ms=norm(p.M,2); c.ds=1e-8;
c.ts=1e-5; c.as=5e-5; c.bs=1e-9; x=zeros(70+16*K,1);
core=c.U'*p.M*c.V/c.ms; x(1:9)=reshape(core(1:3,1:3),[],1);
x(10:12)=core(1:3,4); x(13:15)=core(4,1:3).';
x(16:31)=p.D(:)/c.ds;
x(31+(1:K))=log(abs(p.lambda(:))/(2*pi));
x(31+K+(1:K))=log(-real(p.lambda(:))./abs(p.lambda(:)));
c.order=zeros(4,K); c.rs=zeros(1,K);
for i=1:K
    R=p.R(:,:,i); [~,pivot]=max(vecnorm(R,2,2));
    order=[pivot,setdiff(1:4,pivot,'stable')].'; c.order(:,i)=order;
    right=R(pivot,:); c.rs(i)=norm(right);
    assert(c.rs(i)>0,'Modal initialization has a zero residue.');
    [~,col]=max(abs(right)); left=R(:,col)/right(col);
    v=left(order(2:4)); r=right.'/c.rs(i);
    ix=31+2*K+(i-1)*14+(1:14); x(ix)=[real(v);imag(v);real(r);imag(r)];
end
x(c.delay_offset+(1:3))=p.ty(2:4)/c.ts;
x(c.delay_offset+(4:7))=p.tu/c.ts;
x(c.filter_offset+(1:16))=log(p.a(:)/c.as);
x(c.filter_offset+(17:32))=log(p.b(:)/c.bs);
end

function [lb,ub]=domain(x,c)
K=c.num_modes; lb=-inf(size(x)); ub=inf(size(x));
lb(1:15)=x(1:15)-.2; ub(1:15)=x(1:15)+.2;
fi=31+(1:K); zi=31+K+(1:K); f=exp(x(fi));
lo=max(200,.9*f); hi=min(800,1.1*f);
if K>1
    mid=(f(1:K-1)+f(2:K))/2;
    lo(2:K)=max(lo(2:K),mid); hi(1:K-1)=min(hi(1:K-1),mid);
end
assert(all(lo<hi),'Data-derived frequency intervals overlap incorrectly.');
lb(fi)=log(lo); ub(fi)=log(hi); lb(zi)=log(.001); ub(zi)=log(.08);
di=c.delay_offset+(1:7); lb(di)=-100e-6/c.ts; ub(di)=150e-6/c.ts;
ai=c.filter_offset+(1:16); bi=c.filter_offset+(17:32);
lb(ai)=log(2e-6/c.as); ub(ai)=log(300e-6/c.as);
lb(bi)=log(1e-12/c.bs); ub(bi)=log((150e-6)^2/c.bs);
end
