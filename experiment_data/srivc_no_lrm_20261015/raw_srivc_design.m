function [B,A,q,blocks,plantColumns,boundaryColumns]=raw_srivc_design(f,U,a,referenceHz,direct)
% Continuous-frequency filters. No FRF estimate or LRM enters this function.
f=f(:);nm=size(a,1);nu=size(U,2);q=1i*f./referenceHz(:).';
A=1+q.*a(:,1).'+q.^2.*a(:,2).';
nb=2*(nu+1);extra=2+direct*nu;
B=complex(zeros(numel(f),nm*nb+extra));blocks=cell(nm,1);
plantColumns=[];boundaryColumns=[];
for m=1:nm
    ix=(m-1)*nb+(1:nb);blocks{m}=ix;
    for k=1:nu
        B(:,ix(2*k-1))=U(:,k)./A(:,m);
        B(:,ix(2*k))=q(:,m).*U(:,k)./A(:,m);
    end
    B(:,ix(end-1))=1./A(:,m);B(:,ix(end))=q(:,m)./A(:,m);
    plantColumns=[plantColumns,ix(1:2*nu)]; %#ok<AGROW>
    boundaryColumns=[boundaryColumns,ix(end-1:end)]; %#ok<AGROW>
end
if direct
    ix=nm*nb+(1:nu);B(:,ix)=U;plantColumns=[plantColumns,ix];
end
B(:,end-1:end)=[ones(numel(f),1),1i*f/100];
boundaryColumns=[boundaryColumns,size(B,2)-1,size(B,2)];
end
