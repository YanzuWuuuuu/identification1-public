function [rho,info]=modal_initialization(beta,no,ni,K,nrbm)
noni=no*ni; stride=no+ni;
assert(nrbm<=min(no,ni),'Too many rigid-body directions.');
M=reshape(beta(1:noni),no,ni); [U,S,V]=svd(M,'econ');
rho=zeros(nrbm*stride+K*(stride+2),1);
for k=1:nrbm
    rho((k-1)*stride+(1:stride))=sqrt(S(k,k))*[U(:,k);V(:,k)];
end
info.rigid_singular_values=diag(S);
info.rigid_truncation_relative=norm(M-U(:,1:nrbm)*S(1:nrbm,1:nrbm)*V(:,1:nrbm).','fro')/max(norm(M,'fro'),eps);
for k=1:K
    b=beta(noni+(k-1)*(noni+2)+(1:noni+2));
    x=-b(1)/(2*b(2)); ySquared=1/b(2)-x*x;
    assert(b(1)>0 && b(2)>0 && ySquared>0,'Expected stable underdamped flexible modes.');
    [U,S,V]=svd(reshape(b(3:end)/b(2),no,ni),'econ');
    ix=nrbm*stride+(k-1)*(stride+2)+(1:stride+2);
    rho(ix)=[x;sqrt(ySquared);sqrt(S(1,1))*U(:,1);sqrt(S(1,1))*V(:,1)];
    info.flexible_singular_values(:,k)=diag(S);
end
end
