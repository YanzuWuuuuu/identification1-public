function [cy,cu,sy,su]=clv_compensators(p,f,Ts)
% Causal digital implementation: bilinear proper rational ratio, integer
% delay plus first-order Thiran for the remaining fractional sample.
dy=p.ty-p.ty(1);du=p.tu-p.tu(1);L=400e-6;
cy=complex(zeros(4,numel(f)));cu=cy;sy=cell(4,1);su=cell(4,1);
for j=1:4
    sy{j}=make(p.b(j,1),p.a(j,1),p.b(1,1),p.a(1,1),L-dy(j),Ts);
    su{j}=make(p.b(1,j),p.a(1,j),p.b(1,1),p.a(1,1),L-du(j),Ts);
    cy(j,:)=reshape(freqresp(sy{j},2*pi*f),1,[]);cu(j,:)=reshape(freqresp(su{j},2*pi*f),1,[]);
end
end
function sd=make(b,a,b0,a0,tau,Ts)
assert(tau>=-1e-12,'Noncausal compensation request.');
sd=c2d(tf([b a 1],[b0 a0 1]),Ts,'tustin');
samples=max(tau,0)/Ts;integer=floor(samples);fraction=samples-integer;
if fraction>1e-8
    alpha=(1-fraction)/(1+fraction);sd=sd*tf([alpha 1],[1 alpha],Ts);
end
sd.InputDelay=integer;sd=ss(absorbDelay(sd));assert(isstable(sd));
end
