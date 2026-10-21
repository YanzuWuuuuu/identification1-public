function weighted=apply_sqrt_weight(A,weight)
if isscalar(weight)
    assert(isfinite(weight)&&weight>0,'Weight must be positive.');
    weighted=sqrt(weight)*A;
elseif isvector(weight)
    weight=weight(:);
    assert(numel(weight)==size(A,1)&&all(isfinite(weight)&weight>0),'Invalid diagonal weights.');
    weighted=sqrt(weight).*A;
else
    assert(isequal(size(weight),[size(A,1),size(A,1)]),'Weight dimensions differ.');
    weighted=chol(weight)*A;
end
end
