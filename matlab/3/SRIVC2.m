function [num,belta,info] = SRIVC2(f,frf,Weight,belta0,e,M,K,noni,bounds)
% Frequency-domain additive RIV with scaled solves and stable trial poles.
if nargin<9, bounds=[]; end
assert(size(frf,1)==noni,'FRF channel count mismatch.');
[num,belta,info]=riv_frequency_fit(f,frf,Weight,belta0,zeros(1,K),zeros(1,K),e,M,bounds);
end
