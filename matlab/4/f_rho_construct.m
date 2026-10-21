function [f_rho,Jacobi] = f_rho_construct(rho,no,ni,isubsystem,nrbm)
% Map modal factors to [vec(M); a1; a2; vec(B); ...].
% Flexible denominators are 1+a1*s+a2*s^2. rho uses pole real/imag
% parts and real left/right factors of each monic modal residue.
if nargin<5, nrbm=1; end
noni=no*ni; block=no+ni; offset=nrbm*block;
assert(numel(rho)==offset+isubsystem*(block+2),'Modal parameter size mismatch.');
f_rho=zeros(noni+isubsystem*(noni+2),1);
Jacobi=zeros(numel(f_rho),numel(rho));
M=zeros(no,ni);
for k=1:nrbm
    cols=(k-1)*block+(1:block);
    left=rho(cols(1:no)); right=rho(cols(no+1:end));
    M=M+left*right.';
    Jacobi(1:noni,cols)=[kron(right,eye(no)),kron(eye(ni),left)];
end
f_rho(1:noni)=M(:);
for k=1:isubsystem
    cols=offset+(k-1)*(block+2)+(1:block+2);
    rows=noni+(k-1)*(noni+2)+(1:noni+2);
    x=rho(cols(1)); y=rho(cols(2)); d=x*x+y*y;
    assert(isfinite(d)&&d>0,'A flexible pole must have nonzero finite magnitude.');
    left=rho(cols(3:2+no)); right=rho(cols(3+no:end));
    residue=left*right.'; b=residue(:)/d;
    f_rho(rows)=[-2*x/d;1/d;b];
    % Include normalization derivatives of denominator AND residue.
    Jacobi(rows,cols(1:2))=[-2/d+4*x*x/d^2,4*x*y/d^2; ...
        -2*x/d^2,-2*y/d^2; -2*x*b/d,-2*y*b/d];
    Jacobi(rows(3:end),cols(3:end))= ...
        [kron(right,eye(no)),kron(eye(ni),left)]/d;
end
end
