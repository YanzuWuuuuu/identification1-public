function tests=test_raw_identification
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath=path;
addpath(fullfile(root,'common'));
if isfolder(fullfile(root,'raw')), addpath(fullfile(root,'raw')); end
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testRawInputResponseMatchesIndependentSimulation(testCase)
stream=RandStream('mt19937ar','Seed',14);
u=randn(stream,1500,2); ts=.001; t=(0:size(u,1)-1)'*ts;
M=[.3,-.2;.1,.5;.4,.2]; R=[.5,.2;-.1,.3;.2,-.4];
hz=45; zeta=.035; w=2*pi*hz;
system=M*tf(1,[1,0,0])+R*tf(1,[1,2*zeta*w,w^2]);
expected=lsim(system,u,t,'foh');
[basis,initial]=raw_modal_basis(u,ts,hz,zeta,'foh');
forced=basis*[M.';R.'];
% FOH state conversion can change the initial state; nuisance states absorb it.
correction=initial*((initial\(expected-forced)));
verifyLessThan(testCase,norm(expected-forced-correction,'fro')/norm(expected,'fro'),1e-9);
end

function testRawFitRecoversPolesWithoutAnyFrf(testCase)
stream=RandStream('mt19937ar','Seed',24);
u=randn(stream,2000,2); ts=.001; t=(0:size(u,1)-1)'*ts;
hz=[35,83]; zeta=[.025,.04];
M=[.2,.1;-.3,.4;.1,.5];
system=M*tf(1,[1,0,0]);
residues={ [.5;.2;-.3]*[.6,-.2],[-.2;.4;.3]*[.2,.7] };
for k=1:2
    w=2*pi*hz(k);
    system=system+residues{k}*tf(1,[1,2*zeta(k)*w,w^2]);
end
y=lsim(system,u,t,'foh');
fit=fit_raw_modal(u,y,ts,hz.*[1.04,.96],[.03,.03], ...
    [30,40;75,90],'foh',80);
verifyEqual(testCase,fit.frequency_hz(:),hz(:),'AbsTol',1e-3);
verifyEqual(testCase,fit.damping(:),zeta(:),'AbsTol',1e-4);
verifyLessThan(testCase,norm(fit.prediction-y,'fro')/norm(y,'fro'),1e-7);
verifyTrue(testCase,fit.exitflag>0);
weighted=fit_raw_modal(u,y,ts,hz.*[1.04,.96],[.03,.03], ...
    [30,40;75,90],'foh',80,[0,100]);
verifyEqual(testCase,weighted.frequency_hz(:),hz(:),'AbsTol',1e-3);
verifyEqual(testCase,weighted.damping(:),zeta(:),'AbsTol',1e-4);
verifyLessThan(testCase,norm(weighted.prediction-y,'fro')/norm(y,'fro'),1e-7);
weighted.frequency_hz=weighted.frequency_hz.*[1.02;.98];
structured=fit_raw_structured(u,y,weighted,2,80);
verifyEqual(testCase,structured.frequency_hz(:),hz(:),'AbsTol',1e-3);
verifyEqual(testCase,structured.damping(:),zeta(:),'AbsTol',1e-4);
verifyLessThan(testCase,norm(structured.prediction-y,'fro')/norm(y,'fro'),1e-7);
for k=2:3
    singular=svd(structured.residues(:,:,k));
    verifyLessThan(testCase,singular(2)/singular(1),1e-12);
end
evaluation=evaluate_raw_fit(u,y,structured,1400);
altered=y; altered(1401:end,:)=altered(1401:end,:)+1;
other=evaluate_raw_fit(u,altered,structured,1400);
verifyEqual(testCase,other.prediction,evaluation.prediction,'AbsTol',1e-14);
verifyLessThan(testCase,evaluation.metrics(2),1e-5);
end
