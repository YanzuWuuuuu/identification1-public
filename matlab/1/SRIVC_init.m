function [num,belta] = SRIVC_init(f,frf,den_para,Weight,~,noni,~,~,K)
% Fixed-pole initialization is one linear weighted LS solve, not an iteration.
% Legacy arguments belta0/e/Max_iter are retained for call compatibility.
assert(isvector(f) && isreal(f) && all(isfinite(f(:)) & f(:)>0) && ...
    isequal(size(frf),[noni,numel(f)]) && all(isfinite(frf(:))), ...
    'identification:InvalidFrequencyData', ...
    'Initialization requires finite positive frequencies and matching finite FRFs.');
N=numel(f); X=zeros(2*noni*N,noni*(K+1));
for k=1:N
    s=1i*2*pi*f(k); basis=eye(noni)/s^2;
    for mode=1:K
        A=1+den_para{mode}(2)*s+den_para{mode}(1)*s^2;
        basis=[basis,eye(noni)/A]; %#ok<AGROW>
    end
    X((k-1)*2*noni+(1:2*noni),:)=[real(basis);imag(basis)];
end
belta=scaled_lsq(apply_sqrt_weight(X,Weight),apply_sqrt_weight(stack_complex(frf),Weight));
num=1;
end
