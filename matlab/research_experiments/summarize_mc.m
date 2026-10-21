function summarize_mc(pattern)
here=fileparts(mfilename('fullpath')); d=dir(fullfile(here,pattern));
ftrue=[286.43,344.58,488.27,705.67,738.83]; ztrue=0.01*ones(1,5); K=5;
for i=1:numel(d)
  S=load(fullfile(here,d(i).name)); R=S.T.res; R=R(R(:,end)==1,:); if isempty(R), continue; end
  if isfield(S.T,'ftrue'), ft=S.T.ftrue; zt=S.T.ztrue; else, ft=ftrue; zt=ztrue; end
  fe=100*sqrt(mean((R(:,1:K)-ft).^2))./ft; ze=100*sqrt(mean((R(:,K+1:2*K)-zt).^2))./zt;
  fprintf('%-46s n=%2d | f RMSE%% %s | zeta RMSE%% %s | M %.3f%% | flexFRF %.2f%% | cond raw %.0e scaled %.0f\n', ...
    d(i).name,size(R,1),mat2str(fe,2),mat2str(ze,2),mean(R(:,2*K+1)),mean(R(:,2*K+2)),median(R(:,2*K+3)),median(R(:,2*K+4)));
end
end
