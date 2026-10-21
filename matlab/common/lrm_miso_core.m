function [transientDb,Tr,theta_store,P_LRM,magnitudeDb,phaseDeg]=lrm_miso_core(Y,U,B,R)
% Local rational MISO fit: one common transient, numerator order R and
% denominator order two. Tr repeats the SAME output transient for each
% input for compatibility; its rows are not separate input transients.
Y=Y(:); [nu,n]=size(U); nw=B-1; width=2*nw+1;
assert(any(R==[1,2])&&nw>0,'Use R=1 or 2 and a positive window radius.');
parameters=(nu+1)*(R+1)+2;
assert(n>=width&&width>=parameters,'LRM window is too short for the model.');
P_LRM=zeros(nu,n); Tr=zeros(nu,n); theta_store=zeros(n,parameters+1);
for k=1:n
    first=max(1,min(k-nw,n-width+1)); locations=first:first+width-1;
    r=(locations.'-k)/nw; output=Y(locations); inputs=U(:,locations).';
    X=zeros(width,parameters);
    for degree=0:R
        X(:,degree*nu+(1:nu))=inputs.*r.^degree;
        X(:,nu*(R+1)+degree+1)=r.^degree;
    end
    X(:,end-1:end)=[-r.*output,-r.^2.*output];
    scale=max(vecnorm(X),realmin);
    theta=lsqminnorm(X./scale,output,1e-12)./scale.';
    P_LRM(:,k)=theta(1:nu);
    Tr(:,k)=theta(nu*(R+1)+1);
    theta_store(k,:)=[k,theta.'];
end
magnitudeDb=20*log10(max(abs(P_LRM),realmin));
phaseDeg=rad2deg(angle(P_LRM));
transientDb=20*log10(max(abs(Tr),realmin));
end
