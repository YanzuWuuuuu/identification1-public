function [G,f]=read_lrm_data(filename,no,ni)
data=load(filename); f=data.f(:);
G=zeros(no*ni,numel(f));
for output=1:no
    channel=data.(sprintf('P%d',output));
    assert(isequal(size(channel),[ni,numel(f)]),'FRF data shape mismatch.');
    G(output:no:end,:)=channel;
end
assert(all(isfinite(f))&&all(diff(f)>0),'Frequency grid must be finite and increasing.');
end
