function fit=fit_from_frf(Y,f,sigma,initial,method,startIndex,maxIterations)
% All charts, scales and frequency intervals are made from data-derived start.
if nargin<6, startIndex=1; end
if nargin<7, maxIterations=200; end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'relative_calibration_pilot'));
[c,x]=chart(initial); [lb,ub]=domain(x,c);
if strcmp(method,'physical')
    x=[x(1:102);zeros(8,1)]; lb=[lb(1:102);log(.1e-6/20e-6)*ones(8,1)];
    ub=[ub(1:102);log(150e-6/20e-6)*ones(8,1)];
end
if startIndex>1
    % Prespecified alternative starts, independent of truth and outcome.
    rng(891+startIndex,'twister');
    x(32:35)=x(32:35)+.003*randn(4,1);
    x(36:39)=x(36:39)+.15*randn(4,1);
    x(96:102)=x(96:102)+.8*randn(7,1);
    if strcmp(method,'physical'), x(103:110)=x(103:110)+.35*randn(8,1);
    else, x(103:118)=x(103:118)+.25*randn(16,1); x(119:134)=x(119:134)+.5*randn(16,1); end
end
x=min(max(x,lb+1e-9),ub-1e-9); s=2i*pi*f(:).';
fun=@(xx) objective(xx,c,s,Y,sigma,method);
[rr,J0]=fun(x); direction=sin((1:numel(x)).'); direction=direction/norm(direction);
hh=1e-5; fd=(fun(x+hh*direction)-fun(x-hh*direction))/(2*hh);
je=norm(fd-J0*direction)/norm(fd); assert(je<1e-4,'Jacobian check failed.');
opts=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off', ...
    'MaxIterations',maxIterations,'MaxFunctionEvaluations',4*maxIterations,'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10,'OptimalityTolerance',1e-6,'ScaleProblem','jacobian');
tic; [xh,cost,r,flag,out,~,J]=lsqnonlin(fun,x,lb,ub,opts); sec=toc;
xx=expand(xh,c,method); [G,~,estimate]=rcp_model(xx,c,s);
fit.estimate=estimate; fit.method=method; fit.x=xh; fit.x0=x; fit.c=c;
fit.cost=cost; fit.mean_whitened_residual=cost/numel(Y); fit.exitflag=flag;
fit.iterations=out.iterations; fit.elapsed_seconds=sec; fit.first_order_optimality=out.firstorderopt;
fit.active_bounds=find(abs(xh-lb)<1e-7 | abs(xh-ub)<1e-7).';
fit.normalized_gradient=max(abs(J'*r)./(vecnorm(J).'*max(norm(r),realmin)));
fit.jacobian_check=je; fit.chart_rcond=rcond(reshape(xh(1:9),3,3));
fit.start_index=startIndex; fit.initial_cost=norm(rr)^2;
fit.iteration_budget=maxIterations;
fit.fitted=G; fit.lb=lb; fit.ub=ub;
end

function [r,J]=objective(x,c,s,Y,sigma,method)
if nargout>1
    [xx,T]=expand(x,c,method); [G,Jc]=rcp_model(xx,c,s);
    if strcmp(method,'physical'), Jc=[Jc(:,1:102),Jc(:,103:134)*T]; end
    Jc=Jc./sigma(:); J=[real(Jc);imag(Jc)];
else
    xx=expand(x,c,method); G=rcp_model(xx,c,s);
end
v=(G-Y)./sigma; r=[real(v(:));imag(v(:))];
end

function [xx,T]=expand(x,c,method)
if ~strcmp(method,'physical'), xx=x; T=[]; return; end
ty=20e-6*exp(x(103:106)); tu=20e-6*exp(x(107:110));
a=ty+tu.'; b=ty*tu.';
xx=[x(1:102);log(a(:)/c.as);log(b(:)/c.bs)];
if nargout>1
    T=zeros(32,8);
    for j=1:4, for k=1:4
        idx=j+4*(k-1); T(idx,j)=ty(j)/a(j,k); T(idx,4+k)=tu(k)/a(j,k);
        T(16+idx,j)=1; T(16+idx,4+k)=1;
    end, end
end
end

function [c,x]=chart(p)
[c.U,~,c.V]=svd(p.M); c.ms=norm(p.M,2); c.ds=1e-8;
c.ts=1e-5; c.as=5e-5; c.bs=1e-9; x=zeros(134,1);
K=c.U'*p.M*c.V/c.ms; x(1:9)=reshape(K(1:3,1:3),[],1);
x(10:12)=K(1:3,4); x(13:15)=K(4,1:3).'; x(16:31)=p.D(:)/c.ds;
x(32:35)=log(abs(p.lambda)/(2*pi)); x(36:39)=log(-real(p.lambda)./abs(p.lambda));
for i=1:4
    R=p.R(:,:,i); [~,pivot]=max(vecnorm(R,2,2)); order=[pivot,setdiff(1:4,pivot,'stable')].';
    c.order(:,i)=order; right=R(pivot,:); c.rs(i)=norm(right);
    [~,col]=max(abs(right)); left=R(:,col)/right(col);
    v=left(order(2:4)); r=right.'/c.rs(i);
    ix=39+(i-1)*14+(1:14); x(ix)=[real(v);imag(v);real(r);imag(r)];
end
x(96:98)=p.ty(2:4)/c.ts; x(99:102)=p.tu/c.ts;
x(103:118)=log(p.a(:)/c.as); x(119:134)=log(p.b(:)/c.bs);
end

function [lb,ub]=domain(x,c)
lb=-inf(134,1); ub=inf(134,1);
lb(1:15)=x(1:15)-.2; ub(1:15)=x(1:15)+.2;
f=exp(x(32:35)); lo=max(200,.9*f); hi=min(800,1.1*f);
mid=(f(1:3)+f(2:4))/2; lo(2:4)=max(lo(2:4),mid); hi(1:3)=min(hi(1:3),mid);
lb(32:35)=log(lo); ub(32:35)=log(hi); lb(36:39)=log(.001); ub(36:39)=log(.08);
lb(96:102)=-100e-6/c.ts; ub(96:102)=150e-6/c.ts;
lb(103:118)=log(2e-6/c.as); ub(103:118)=log(300e-6/c.as);
lb(119:134)=log(1e-12/c.bs); ub(119:134)=log((150e-6)^2/c.bs);
end
