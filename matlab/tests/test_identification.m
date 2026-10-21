function tests = test_identification
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.root = root;
testCase.TestData.oldPath = path;
testCase.TestData.oldVisible = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
addpath(root,fullfile(root,'0'),fullfile(root,'1'),fullfile(root,'2'),fullfile(root,'3'),fullfile(root,'4'));
if isfolder(fullfile(root,'common')), addpath(fullfile(root,'common')); end
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
set(groot,'DefaultFigureVisible',testCase.TestData.oldVisible);
end

function testAllRigidDirectionsAreReconstructed(testCase)
left = [1,0,0;0,2,0;0,0,3];
right = [1,2,0;0,1,0;0,0,1;1,0,1];
rho = zeros(21,1);
for k=1:3, rho((k-1)*7+(1:7))=[left(:,k);right(:,k)]; end
actual = f_rho_construct(rho,3,4,0,3);
verifyEqual(testCase,actual,reshape(left*right.',[],1),'AbsTol',1e-12);
end

function testNormalizedJacobianMatchesFiniteDifferences(testCase)
rho = [1;2;-1;0.5;1;2;0.3; -0.4;0.3;0.8;1;-0.2;0.6;-0.1; ...
       -12;1800;0.3;-0.2;0.1;0.4;0.2;0.7;-0.3];
[~,J] = f_rho_construct(rho,3,4,1,2);
numeric = zeros(size(J));
for k=1:numel(rho)
    h=1e-5*max(abs(rho(k)),1);
    plus=rho; minus=rho; plus(k)=plus(k)+h; minus(k)=minus(k)-h;
    numeric(:,k)=(f_rho_construct(plus,3,4,1,2)-f_rho_construct(minus,3,4,1,2))/(2*h);
end
columnError=vecnorm(J-numeric)./max(vecnorm(numeric),1e-15);
verifyLessThan(testCase,max(columnError),1e-6);
end

function testProjectionStopsAtAnExactSolution(testCase)
rho = [1;0.4;0.8;0.2;-2;50;0.3;0.4;0.7;0.2];
beta = f_rho_construct(rho,2,2,1);
[iterations,actual] = Gauss_Newton(beta,rho,0.01,10,8,1e-8,2,2,1);
verifyLessThan(testCase,iterations,3);
verifyEqual(testCase,f_rho_construct(actual,2,2,1),beta,'AbsTol',1e-11);
end

function testProjectionReducesResidualAndReportsConvergence(testCase)
rho = [1;0.4;0.8;0.2;-2;50;0.3;0.4;0.7;0.2];
beta = f_rho_construct(rho,2,2,1);
initial=rho.*[1.1;0.9;1;1.05;1.2;1.04;0.9;1.1;1.1;0.9];
[~,actual,info] = Gauss_Newton(beta,initial,0.01,10,200,1e-9,2,2,1);
verifyLessThan(testCase,norm(f_rho_construct(actual,2,2,1)-beta)/norm(beta),1e-6);
verifyTrue(testCase,info.converged);
verifyLessThanOrEqual(testCase,max(diff(info.cost_history)),1e-12);
end

function testFourModeRivRecovery(testCase)
f=linspace(5,180,400).'; noni=4; K=4;
freq=[30,65,105,145]; damping=[.02,.025,.015,.03];
residues=[1,.2,-.3,.7;.3,1,.4,-.2;.5,-.4,1,.3;.2,.6,-.1,1];
truth=zeros((noni+2)*K,K); initial=truth; G=zeros(noni,numel(f));
for k=1:K
    w=2*pi*freq(k); idx=(k-1)*(noni+2)+(1:noni+2);
    truth(idx,k)=[2*damping(k)/w;1/w^2;residues(:,k)/w^2];
    initial(idx,k)=truth(idx,k).*[1.15;1.04;1.1;0.9;1.1;0.9];
    s=1i*2*pi*f.';
    G=G+residues(:,k)./(s.^2+2*damping(k)*w*s+w^2);
end

[~,actual] = SRIVC2(f,G,1,initial,1e-9,100,K,noni);
fitted=zeros(size(G));
for k=1:K
    idx=(k-1)*(noni+2)+(1:noni+2); b=actual(idx,k);
    fitted=fitted+b(3:end)./(1+b(1)*s+b(2)*s.^2);
end
verifyLessThan(testCase,norm(fitted-G,'fro')/norm(G,'fro'),1e-5);
verifySize(testCase,actual,[(noni+2)*K,K]);
end

function testProjectionUsesTheSuppliedMetric(testCase)
rho=[1;0.4;0.8;0.2;-2;50;0.3;0.4;0.7;0.2];
beta=f_rho_construct(rho,2,2,1);
beta(4)=beta(4)+0.2;
metric=diag([1;3;2;5;1e3;1e4;1e4;1e4;1e4;1e4]);
[~,actual,info]=Gauss_Newton(beta,rho,0.01,10,200,1e-9,2,2,1,1,metric);
expected=norm(metric*(f_rho_construct(actual,2,2,1)-beta));
verifyEqual(testCase,info.weighted_residual_norm,expected,'RelTol',1e-8);
end

function testRigidRivPreservesComplexPhase(testCase)
f=linspace(5,60,150).'; s=1i*2*pi*f.'; noni=4;
M=[1;.2;-.1;.7]; R=[.2;-.3;.5;.1]; w=2*pi*40; a=[2*.02/w;1/w^2];
G=M./s.^2+R./(s.^2+2*.02*w*s+w^2);
initial=[a.*[1.1;1.04];M;M*a(1);M*a(2)+R/w^2];
[~,b] = SRIVC1(f,G,1,initial,noni,1e-9,100,1);
fitted=(b(3:6)+b(7:10)*s+b(11:14)*s.^2)./(s.^2.*(1+b(1)*s+b(2)*s.^2));
verifyLessThan(testCase,norm(fitted-G,'fro')/norm(G,'fro'),1e-5);
end

function testNumericalEntriesRejectInvalidFrequencyData(testCase)
for invalid=[0,NaN,Inf]
    f=[1;invalid;3]; G=ones(4,3);
    verifyError(testCase,@()SRIVC_init(f,G,{[1,1]},1,[],4,[],[],1), ...
        'identification:InvalidFrequencyData');
    verifyError(testCase,@()riv_frequency_fit(f,G,1,ones(6,1),0,0,1e-8,10), ...
        'identification:InvalidFrequencyData');
end
for G={ones(4,2),[ones(4,2),NaN(4,1)]}
    verifyError(testCase,@()SRIVC_init([1;2;3],G{1},{[1,1]},1,[],4,[],[],1), ...
        'identification:InvalidFrequencyData');
    verifyError(testCase,@()riv_frequency_fit([1;2;3],G{1},1,ones(6,1),0,0,1e-8,10), ...
        'identification:InvalidFrequencyData');
end
end

function testProjectionIsInvariantToOverallMetricScale(testCase)
rho=[1;0.4;0.8;0.2;-2;50;0.3;0.4;0.7;0.2];
beta=f_rho_construct(rho,2,2,1);
initial=rho.*[1.1;0.9;1;1.05;1.2;1.04;0.9;1.1;1.1;0.9];
metric=diag([1;3;2;5;1e3;1e4;1e4;1e4;1e4;1e4]);
[~,actual,info]=Gauss_Newton(beta,initial,0.01,10,200,1e-9,2,2,1,1,1e-9*metric);
verifyLessThan(testCase,norm(f_rho_construct(actual,2,2,1)-beta)/norm(beta),1e-6);
verifyTrue(testCase,info.converged);
end

function testLrmSeparatesTheCommonTransient(testCase)
stream=RandStream('mt19937ar','Seed',21); n=81; center=41;
U=randn(stream,4,n)+1i*randn(stream,4,n); r=((1:n)-center)/14;
B=[1+2i;-.3+.4i;.8-.2i;-.5-.6i]; slope=[.2;-.1;.3;-.2];
T=0.4+0.7i;
Y=(sum((B+slope*r).*U,1)+T+(.2-.1i)*r)./(1+(.1+.05i)*r+.03*r.^2);
[~,transient,~,G]=LRM_MISO(Y.',U,15,1);
verifyEqual(testCase,G(:,center),B,'AbsTol',1e-9);
verifyEqual(testCase,transient(:,center),repmat(T,4,1),'AbsTol',1e-9);
end

function testRefinementHonorsFrequencyBound(testCase)
f=linspace(10,60,180).'; s=1i*2*pi*f.';
residue=[1;.3;-.2;.5]; w=2*pi*45;
G=residue./(s.^2+2*.04*w*s+w^2);
w0=2*pi*30; initial=[2*.04/w0;1/w0^2;residue/w0^2];
[actual,info]=refine_flexible_fit(initial,f,G,1,1,[25,40],100,1e-8);
frequency=1/(2*pi*sqrt(actual(2)));
verifyLessThanOrEqual(testCase,frequency,40+1e-7);
verifyGreaterThan(testCase,frequency,39.9);
verifyTrue(testCase,info.converged);
verifyLessThan(testCase,info.final_cost,info.initial_cost);
end

function testFrequencySelectionUsesActualPositiveBins(testCase)
f=[0;0.5;1]; G=[NaN,1+2i,3+4i;NaN,5+6i,7+8i];
[actual,frequencies]=select_frequency_data(G,[.5,.9],f.');
verifyEqual(testCase,actual,G(:,2:3));
verifyEqual(testCase,frequencies,[.5;1]);
end

function testPipelineHasConsistentDataAndFourModes(testCase)
result=run_identification(false);
cfg=result.stage1.cfg; n=cfg.no*cfg.ni;
verifySize(testCase,result.stage1.belta_init,[n*(cfg.K+1),1]);
verifySize(testCase,result.stage3.belta_vecest,[n+cfg.K*(n+2),1]);
verifySize(testCase,result.stage4.rou_est,[cfg.nrbm*(cfg.no+cfg.ni)+cfg.K*(cfg.no+cfg.ni+2),1]);
rigid=result.stage2;
expected=rigid.G-rigid.vecM_est./((1i*2*pi*rigid.f.').^2);
verifyEqual(testCase,rigid.flexible_frf,expected,'AbsTol',1e-14);
stage=result.stage3;
expected=additive_frf(stage.belta_vecest,stage.f_select,n,cfg.K)- ...
    rigid.vecM_est./((1i*2*pi*stage.f_select.').^2);
verifyEqual(testCase,stage.frf_est,expected,'AbsTol',1e-14);
verifyTrue(testCase,stage.refinement.converged);
verifyTrue(testCase,result.stage4.info.converged);
beta=result.stage4.belta_modal;
singular=svd(reshape(beta(1:n),cfg.no,cfg.ni));
verifyLessThan(testCase,singular(cfg.nrbm+1)/singular(1),1e-12);
for k=1:cfg.K
    ix=n+(k-1)*(n+2)+(1:n+2);
    verifyGreaterThan(testCase,beta(ix(1:2)),[0;0]);
    singular=svd(reshape(beta(ix(3:end)),cfg.no,cfg.ni));
    verifyLessThan(testCase,singular(2)/singular(1),1e-12);
end
verifyEqual(testCase,height(result.modes),4);
verifyGreaterThan(testCase,max(result.stage4.validation_f),max(cfg.flexible_frequency_hz));
saved=load(fullfile(cfg.result_dir,'stage3.mat'));
verifyEqual(testCase,saved.belta_vecest,stage.belta_vecest);
staleCfg=cfg; staleCfg.source_bytes=staleCfg.source_bytes+1;
verifyError(testCase,@()load_stage_result(staleCfg,1),'identification:StaleStage');
end
