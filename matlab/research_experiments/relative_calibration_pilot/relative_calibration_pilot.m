function summary=relative_calibration_pilot(which)
% Local identification feasibility pilot. No LRM or legacy fifth mode.
if nargin<1, which='all'; end
folder=fileparts(mfilename('fullpath')); outdir=fullfile(folder,'outputs');
if ~isfolder(outdir), mkdir(outdir); end
[base,c,nominal]=makeBase();
f=unique([logspace(log10(2),log10(800),100),linspace(2,800,160), ...
    reshape(base.freq.*(1+linspace(-.055,.055,27).'),1,[])]).';
fv=unique([logspace(log10(2.07),log10(798.7),150),linspace(2.13,799.8,420), ...
    reshape(base.freq.*(1+linspace(-.05,.05,42).'),1,[])]).';
names={'common','repeated','distinct'}; levels=[0,1]; summary=struct([]);
for ic=1:3
    truth=base;
    switch ic
        case 1
            truth.filter_y=25e-6*ones(4,1); truth.filter_u=truth.filter_y;
            truth.ty=zeros(4,1); truth.tu=12e-6*ones(4,1);
        case 2
            truth.filter_y=[20;20;35;35]*1e-6; truth.filter_u=[20;35;20;35]*1e-6;
            truth.ty=[0;8;17;5]*1e-6; truth.tu=[4;17;10;25]*1e-6;
        case 3
            truth.filter_y=[20;31;9;25]*1e-6; truth.filter_u=[17;11;23;35]*1e-6;
            truth.ty=[0;8;17;5]*1e-6; truth.tu=[4;17;10;25]*1e-6;
    end
    truth.a=truth.filter_y+truth.filter_u.';
    truth.b=truth.filter_y*truth.filter_u.';
    xt=encode(truth,c); [Gt,Jt]=rcp_model(xt,c,1i*2*pi*f);
    sigma=.001*sqrt(abs(Gt).^2+(1e-7)^2);
    [Gv]=rcp_model(xt,c,1i*2*pi*fv); sv=.001*sqrt(abs(Gv).^2+(1e-7)^2);
    rng(20261013,'twister'); noise=(randn(size(Gt))+1i*randn(size(Gt)))/sqrt(2);
    nv=(randn(size(Gv))+1i*randn(size(Gv)))/sqrt(2);
    rng(31013,'twister'); dx=zeros(134,1);
    dx(1:15)=.003*randn(15,1); dx(16:31)=.1*randn(16,1);
    dx(32:35)=.001*randn(4,1); dx(36:39)=.08*randn(4,1);
    dx(40:95)=.015*randn(56,1); dx(96:102)=.25*randn(7,1);
    dx(103:134)=.02*randn(32,1); x0=xt+dx;
    [lb,ub]=bounds(xt,c);
    fun=@(x) residual(x,c,1i*2*pi*f,Gt,sigma);
    [~,Ja]=fun(x0); rng(42); direction=randn(134,1); direction=direction/norm(direction);
    step=1e-5; Jfd=(fun(x0+step*direction)-fun(x0-step*direction))/(2*step);
    jacerr=norm(Jfd-Ja*direction)/norm(Jfd);
    assert(jacerr<1e-5,'Analytic Jacobian check failed: %g',jacerr);
    Js=[real(Jt./sigma(:));imag(Jt./sigma(:))];
    cs=vecnorm(Js); singular=svd(Js./cs,0);
    for level=levels
        tag=sprintf('%s_noise%d',names{ic},level);
        if ~strcmp(which,'all') && ~strcmp(which,tag), continue; end
        fprintf('Starting %s: %d frequencies, %d parameters\n',tag,numel(f),numel(xt));
        Y=Gt+level*sigma.*noise;
        opts=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off', ...
            'MaxIterations',150,'MaxFunctionEvaluations',450, ...
            'FunctionTolerance',1e-11,'StepTolerance',1e-11,'OptimalityTolerance',1e-6, ...
            'ScaleProblem','jacobian');
        tic;
        [xh,resnorm,res,flag,output,~,J]=lsqnonlin(@(x) residual(x,c,1i*2*pi*f,Y,sigma),x0,lb,ub,opts);
        elapsed=toc; [Gh,~,estimate]=rcp_model(xh,c,1i*2*pi*f);
        Ghv=rcp_model(xh,c,1i*2*pi*fv);
        metrics=rcp_metrics(truth,estimate,nominal,fv);
        metrics.case=tag; metrics.noise_multiplier=level;
        metrics.exitflag=flag; metrics.iterations=output.iterations;
        metrics.first_order_optimality=output.firstorderopt; metrics.elapsed_seconds=elapsed;
        metrics.train_relative_frf_error=norm(Gh-Gt,'fro')/norm(Gt,'fro');
        metrics.validation_relative_frf_error=norm(Ghv-Gv,'fro')/norm(Gv,'fro');
        metrics.validation_whitened_mse=mean(abs((Ghv-Gv-level*sv.*nv)./sv).^2,'all');
        metrics.train_whitened_mse=resnorm/numel(Gt);
        metrics.jacobian_directional_relative_error=jacerr;
        metrics.truth_column_normalized_jacobian_smin=singular(end);
        metrics.truth_column_normalized_jacobian_condition=singular(1)/singular(end);
        metrics.frequency_relative_error=(abs(estimate.lambda)./(2*pi)-truth.freq(:))./truth.freq(:);
        metrics.damping_relative_error=(-real(estimate.lambda)./abs(estimate.lambda)-.01)/.01;
        metrics.parameters_at_bounds=sum(abs(xh-lb)<1e-7 | abs(xh-ub)<1e-7);
        save(fullfile(outdir,[tag '.mat']),'truth','estimate','nominal','c','f','fv','sigma', ...
            'Y','xt','x0','xh','metrics','output','J','res','lb','ub');
        compact=metrics;
        if isfield(compact,'curves'), compact=rmfield(compact,'curves'); end
        writeJson(fullfile(outdir,[tag '.json']),compact);
        if isempty(summary), summary=compact; else, summary(end+1)=compact; end %#ok<AGROW>
        fprintf('%s done: flag %d, iter %d, validation MSE %.5g, %.1fs\n', ...
            tag,flag,output.iterations,metrics.validation_whitened_mse,elapsed);
    end
end
if ~isempty(summary), writeJson(fullfile(outdir,['summary_' which '.json']),summary); end
end

function [res,J]=residual(x,c,s,Y,sigma)
if nargout>1
    [G,Jc]=rcp_model(x,c,s); Jc=Jc./sigma(:); J=[real(Jc);imag(Jc)];
else
    G=rcp_model(x,c,s);
end
r=(G-Y)./sigma; res=[real(r(:));imag(r(:))];
end

function [p,c,nominal]=makeBase()
B=[.0972,.0990,-.1405,-.0024;.1004,-.1002,9.8523e-4,.0011;.9717,.9893,.9564,1.0781];
C=[.1714,.0861,1.0014;.1697,-.0756,.9986;-.1708,-.1108,1.0012;-.1673,.1034,.9989];
BF=[-.5405,.5613,.0015,-.0292;.0725,.0316,.6871,-.1560;.2364,.2200,.0601,-.8658;-.2834,-.3269,-.5538,-.0619];
CF=[-.5297,.5215,.6477,.4121;.5585,.5490,.5350,.0562;-.6198,.4402,.6328,-.0175;.6258,.4527,.6219,.3133];
p.M=C*diag([1/.0439,1/.0384,1/4.2809])*B;
p.D=1e-8*eye(4); p.freq=[286.43,344.58,488.27,705.67];
w=2*pi*p.freq(:); p.lambda=-.01*w+1i*w*sqrt(1-.01^2); p.R=complex(zeros(4,4,4));
for i=1:4, p.R(:,:,i)=CF(:,i)*BF(i,:)/(p.lambda(i)-conj(p.lambda(i))); end
rng(808,'twister');
Cn=C+.005*norm(C,'fro')/sqrt(numel(C))*randn(size(C));
Bn=B+.005*norm(B,'fro')/sqrt(numel(B))*randn(size(B));
nominal.M=Cn*diag([1/.0439,1/.0384,1/4.2809])*Bn;
[c.U,~,c.V]=svd(nominal.M); c.ms=norm(p.M,2); c.ds=1e-8;
c.ts=1e-5; c.as=5e-5; c.bs=1e-9;
for i=1:4
    [~,pivot]=max(abs(CF(:,i))); c.order(:,i)=[pivot,setdiff(1:4,pivot,'stable')].';
    c.rs(i)=norm(p.R(pivot,:,i));
end
end

function x=encode(p,c)
x=zeros(134,1); K=c.U'*p.M*c.V/c.ms;
x(1:9)=reshape(K(1:3,1:3),[],1); x(10:12)=K(1:3,4); x(13:15)=K(4,1:3).';
x(16:31)=p.D(:)/c.ds; x(32:35)=log(abs(p.lambda)/(2*pi)); x(36:39)=log(-real(p.lambda)./abs(p.lambda));
for i=1:4
    order=c.order(:,i); right=p.R(order(1),:,i); [~,col]=max(abs(right));
    left=p.R(:,col,i)/right(col); v=left(order(2:4)); r=right.'/c.rs(i);
    ix=39+(i-1)*14+(1:14); x(ix)=[real(v);imag(v);real(r);imag(r)];
end
x(96:98)=p.ty(2:4)/c.ts; x(99:102)=p.tu/c.ts;
x(103:118)=log(p.a(:)/c.as); x(119:134)=log(p.b(:)/c.bs);
end

function [lb,ub]=bounds(xt,c)
lb=-inf(134,1); ub=inf(134,1);
lb(1:15)=xt(1:15)-.2; ub(1:15)=xt(1:15)+.2;
lb(32:35)=xt(32:35)+log(.97); ub(32:35)=xt(32:35)+log(1.03);
lb(36:39)=log(.002); ub(36:39)=log(.05);
lb(96:102)=-50e-6/c.ts; ub(96:102)=100e-6/c.ts;
lb(103:118)=log(2e-6/c.as); ub(103:118)=log(200e-6/c.as);
lb(119:134)=log(1e-12/c.bs); ub(119:134)=log((150e-6)^2/c.bs);
end

function writeJson(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); clean=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
