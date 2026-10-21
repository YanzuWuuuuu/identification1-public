function [selected,f_actual,indices]=select_frequency_data(G,requested,f)
f=f(:); requested=requested(:);
assert(size(G,2)==numel(f),'FRF and frequency lengths differ.');
valid=f>0 & isfinite(f) & all(isfinite(G),1).';
fv=f(valid); positions=find(valid);
assert(~isempty(fv),'No positive finite frequency samples.');
assert(all(requested>=min(fv)&requested<=max(fv)),'Requested frequency is outside the valid data range.');
indices=interp1(fv,positions,requested,'nearest');
indices=unique(indices,'stable');
selected=G(:,indices); f_actual=f(indices);
end
