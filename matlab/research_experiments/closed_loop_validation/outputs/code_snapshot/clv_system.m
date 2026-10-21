function [sys,K,info,plant]=clv_system(p,nominal,cfg)
% Independent generator: continuous SISO cascades -> MATLAB c2d -> MIMO.
% The estimator does not call this function or receive these states.
wr=2*pi*1000;nModes=numel(p.lambda);n=2+2*(1+nModes);cells=cell(16,1);
ell=p.ty+p.tu.';
for edge=1:16
    a=p.a(edge);b=p.b(edge);A=zeros(n);B=zeros(n,1);C=zeros(1,n);
    A(1:2,1:2)=[0 wr;-1/(b*wr) -a/b];B(2)=1/(b*wr);
    for mode=0:nModes
        ix=3+2*mode+(0:1); A(ix(2),1)=1/wr;
        if mode==0,A(ix,ix)=[0 wr;0 0];C(ix)=[p.M(edge) 0];
        else
            lam=p.lambda(mode);w=abs(lam);zeta=-real(lam)/w; R=p.R(:,:,mode);R=R(edge);
            A(ix,ix)=[0 wr;-w*w/wr -2*zeta*w];C(ix)=[-2*real(R*conj(lam)),2*wr*real(R)];
        end
    end
    ct=ss(A,B,C,0);ct.InputDelay=ell(edge);
    cells{edge}=c2d(ct,cfg.Ts,'zoh');
end
full=append(cells{:});left=zeros(4,16);right=zeros(16,4);
for edge=1:16,[j,k]=ind2sub([4 4],edge);left(j,edge)=1;right(edge,k)=1;end
full=absorbDelay(left*full*right); full=prescale(full);
plant=minreal(full,1e-12,false);
checkf=[20 200 300 500 710 800];g0=freqresp(full,2*pi*checkf);g1=freqresp(plant,2*pi*checkf);
info.realization_reduction_relative_error=norm(g0(:)-g1(:))/norm(g0(:));
assert(info.realization_reduction_relative_error<1e-7, ...
    'Plant reduction changed transfer: %.6g, %d to %d states.', ...
    info.realization_reduction_relative_error,size(full.A,1),size(plant.A,1));
[Un,Sn,Vn]=svd(nominal.M);wc=2*pi*15;wd=2*pi*120;kp=wc^2;kd=2*.8*wc;
kc=c2d(ss(tf([kd+kp/wd,kp],[1/wd,1])),cfg.Ts,'tustin');
K=minreal(Vn(:,1:3)*diag(1./diag(Sn(1:3,1:3)))*kc*Un(:,1:3)',1e-10,false);
[Ap,Bp,Cp,Dp]=ssdata(plant);[Ak,Bk,Ck,Dk]=ssdata(K);
assert(norm(Dp,'fro')<1e-14,'Unexpected discrete plant direct feedthrough.');
Acl=[Ap-Bp*Dk*Cp,-Bp*Ck;Bk*Cp,Ak];
Bcl=[Bp,-Bp*Dk;zeros(size(Ak,1),4),Bk];
Ccl=[Cp,zeros(4,size(Ak,1));-Dk*Cp,-Ck];Dcl=[zeros(4),eye(4);eye(4),-Dk];
sys=minreal(ss(Acl,Bcl,Ccl,Dcl,cfg.Ts),1e-10,false);
info.plant_states=size(Ap,1);info.closed_states=size(sys.A,1);
info.closed_spectral_radius=max(abs(eig(sys.A)));info.controller_hz=15;
assert(info.closed_spectral_radius<1-1e-8,'Closed loop is not internally stable.');
end
