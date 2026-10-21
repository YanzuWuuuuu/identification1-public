function rows=run_band_audit()
here=fileparts(mfilename('fullpath'));parent=fullfile(fileparts(here),'closed_loop_validation');addpath(parent);
S=load(fullfile(parent,'outputs','data_distinct.mat'),'p','nom','K','cfg');
V=load(fullfile(parent,'outputs_structured','test_distinct.mat'),'testDistinct');
X=load(fullfile(parent,'outputs_structured','data_common.mat'),'common','KC','testCommon');
allrows=jsondecode(fileread(fullfile(parent,'outputs_structured','summary.json')));rows=cell(numel(allrows),1);
for j=1:numel(allrows)
 e=allrows(j);Q=load(fullfile(parent,'outputs_structured',e.case,'selected.mat'),'selected');
 p=S.p;K=S.K;test=V.testDistinct;
 if strcmp(e.case,'common_N32'),p=X.common;K=X.KC;test=X.testCommon;end
 met=audit_band_metrics(p,S.nom,Q.selected,K,S.cfg,test);
 r=struct('case',e.case,'calibration',struct());
 for name={'static','joint','oracle'}
  n=name{1};r.calibration.(n)=met.calibration.(n).output_only;
  % Original saved broadband and low-band results must be reproduced.
  prior=e.metrics.calibration.(n).output_only;
  assert(abs(r.calibration.(n).low.rigid_rms_m/prior.low.rigid_rms_m-1)<1e-7);
 end
 rows{j}=r;fprintf('band audit %s complete\n',e.case);
end
rows=vertcat(rows{:});fid=fopen(fullfile(here,'resonance_validation.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(rows,'PrettyPrint',true));
end
