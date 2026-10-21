function fit=clv_fit(data,pilot,K,cfg,initial,startIndex,maxIterations)
if nargin<7,maxIterations=cfg.max_iterations;end
[x,c,lb,ub]=clv_chart(initial,cfg);c.amin=.1*min(svd(reshape(x(1:9),3,3)));
if startIndex>1
    rng(891+startIndex,'twister');
    base=x;found=false;
    for attempt=1:100
        x=base;x(c.fo+(1:c.K))=base(c.fo+(1:c.K))+.03*randn(c.K,1);
        x(c.zo+(1:c.K))=base(c.zo+(1:c.K))+.3*randn(c.K,1);
        ty=exp(log(5e-6)+rand(4,1)*log(80/5));tu=exp(log(5e-6)+rand(4,1)*log(80/5));
        aa=ty+tu.';bb=ty*tu.';
        if c.additive_a,x(c.ao+(1:7))=[aa(2:4,1)-aa(1,1);aa(1,:).']/c.as;
        else,x(c.ao+(1:16))=log(aa(:)/c.as);end
        x(c.bo+(1:16))=log(bb(:)/c.bs);
        dy=100e-6*rand(4,1);du=150e-6*rand(4,1);
        x(c.to+(1:3))=(dy(2:4)-dy(1))/c.ts;x(c.to+(4:7))=(du+dy(1))/c.ts;
        [viol,~]=clv_constraints(x,c,cfg);
        if all(x>lb&x<ub)&&max(viol)<0,found=true;break;end
    end
    assert(found,'No feasible multistart seed found.');
end
Y=data.Y;yr=[real(reshape(Y,[],data.groups));imag(reshape(Y,[],data.groups))];
pilotMean=[real(pilot.mean(:));imag(pilot.mean(:))];
% Per-record coordinate variance, frozen across every optimizer iteration.
variance=var(yr,0,2);floorScale=max(abs(pilotMean),cfg.sensor_noise_std);
sigma=sqrt(max(variance,1e-6*floorScale.^2));
if data.sample_noise_std==0
    % Noiseless numerical gate uses a fixed nonzero instrument scale.
    sigma(:)=cfg.sensor_noise_std*sqrt(2/cfg.period_samples);
end
ybar=mean(yr,2);Kf=freqresp(K,2*pi*data.f);m=numel(ybar);
cacheX=[];cacheQ=[];cacheG=[];cacheJ=[];cacheE=[];
Alin=[];blin=[];
if c.additive_a
    Alin=zeros(32,c.n);Alin(:,c.ao+(1:7))=[-c.a_map;c.a_map];
    blin=[-cfg.filter_a_bounds(1)*ones(16,1);cfg.filter_a_bounds(2)*ones(16,1)]/c.as;
end
fun=@objective;[q0,~]=fun(x);J0=cacheJ;
direction=sin((1:c.n).');direction=direction/norm(direction);hh=1e-5;
[yp,~]=clv_mean(x+hh*direction,c,data.f,cfg.Ts,Kf,data.D);
[ym,~]=clv_mean(x-hh*direction,c,data.f,cfg.Ts,Kf,data.D);
fd=([real(yp(:)-ym(:));imag(yp(:)-ym(:))]./sigma)/(2*hh*sqrt(m));
jacError=norm(fd-J0*direction)/max(norm(fd),realmin);assert(jacError<1e-4,'Closed Jacobian error %.3g',jacError);
opts=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
    'Display','off','MaxIterations',maxIterations,'MaxFunctionEvaluations',cfg.max_function_evaluations, ...
    'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-8,'StepTolerance',1e-10);
tic;
if isfield(cfg,'solver')&&strcmp(cfg.solver,'constrained-lsq')
    opts=optimoptions('lsqnonlin','Algorithm','interior-point','SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true,'Display','off','MaxIterations',maxIterations, ...
        'MaxFunctionEvaluations',cfg.max_function_evaluations,'OptimalityTolerance',1e-6, ...
        'ConstraintTolerance',1e-8,'StepTolerance',1e-10,'FunctionTolerance',1e-12);
    [xh,cost,~,flag,out]=lsqnonlin(@residual,x,lb,ub,Alin,blin,[],[],@(v)clv_constraints(v,c,cfg),opts);
else
    [xh,cost,flag,out]=fmincon(fun,x,Alin,blin,[],[],lb,ub,@(v)clv_constraints(v,c,cfg),opts);
end
elapsed=toc;
[~,gradient]=fun(xh);J=cacheJ;e=cacheE;[viol,~]=clv_constraints(xh,c,cfg);
estimate=clv_decode(xh,c);singular=svd(J,0);
fit=struct('estimate',estimate,'x',xh,'x0',x,'chart',c,'cost',cost,'initial_cost',q0, ...
    'exitflag',flag,'iterations',out.iterations,'function_count',out.funcCount,'elapsed_seconds',elapsed, ...
    'first_order_optimality',out.firstorderopt,'gradient_inf',norm(gradient,inf), ...
    'constraint_violation',max([0;viol;lb-xh;xh-ub]),'jacobian_check',jacError, ...
    'normalized_gradient',max(abs(J'*e)./(max(vecnorm(J).',realmin)*max(norm(e),realmin))), ...
    'jacobian_singular_max',singular(1),'jacobian_singular_min',singular(end), ...
    'active_bounds',find(abs(xh-lb)<1e-7|abs(xh-ub)<1e-7).', ...
    'start_index',startIndex,'max_iterations',maxIterations,'sigma',sigma,'solver_message',out.message);
fit.solver_pass=flag>0&&out.firstorderopt<=1e-6&&fit.constraint_violation<=1e-8;
fit.asymptotic_score_tolerance_pass=norm(gradient)<=min(1e-6,data.groups^(-3/4))&&isempty(fit.active_bounds);
    function [e0,J0]=residual(xx)
        objective(xx);e0=cacheE;J0=cacheJ;
    end
    function [q,g]=objective(xx)
        if isequal(xx,cacheX),q=cacheQ;g=cacheG;return;end
        [mu,jc]=clv_mean(xx,c,data.f,cfg.Ts,Kf,data.D);
        ee=([real(mu(:));imag(mu(:))]-ybar)./sigma/sqrt(m);
        JJ=[real(jc);imag(jc)]./sigma/sqrt(m);
        q=ee'*ee;g=2*JJ'*ee;
        cacheX=xx;cacheQ=q;cacheG=g;cacheJ=JJ;cacheE=ee;
    end
end
