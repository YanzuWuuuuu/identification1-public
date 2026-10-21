function [initial,info,Gpilot]=clv_initial(data,K,cfg)
% This interface intentionally accepts no truth, source model or nominal M.
here=fileparts(mfilename('fullpath'));addpath(fullfile(fileparts(here),'relative_calibration_validation'));
f=data.f;nf=numel(f);Kf=freqresp(K,2*pi*f);Gpilot=complex(zeros(16,nf));
for k=1:nf
    T=reshape(data.mean(:,k),4,4)/data.D(:,:,k); Z=eye(4)-Kf(:,:,k)*T;
    assert(rcond(Z)>1e-8,'Pilot feedback inversion is ill-conditioned.');
    Gpilot(:,k)=reshape(T/Z,16,1);
end
% Removing the known fundamental ZOH factor is an initializer only.
% All alias effects remain in the exact final predictor.
s=2i*pi*f.';hold=(1-exp(-s*cfg.Ts))./(s*cfg.Ts);
[initial,info]=init_five_from_frf(Gpilot./hold,f,struct('num_modes',5,'frequency_band',[200 800]));
initial.D=zeros(4);initial.ty=zeros(4,1);initial.tu=10e-6*ones(4,1);
z=exp(s*cfg.Ts);Kmode=numel(initial.lambda);
X=complex(zeros(nf,1+2*Kmode));X(:,1)=clv_kernel(0,0,initial.a(1),initial.b(1),10e-6,z,cfg.Ts,false,false).';
for i=1:Kmode
    w=abs(initial.lambda(i));zz=-real(initial.lambda(i))/w;
    H=clv_kernel(w,zz,initial.a(1),initial.b(1),10e-6,z,cfg.Ts,true,false);
    X(:,2*i:2*i+1)=H.';
end
Ar=[real(X);imag(X)];scale=vecnorm(Ar);coef=(Ar./scale)\[real(Gpilot.');imag(Gpilot.')];coef=coef./scale.';
[U,S,V]=svd(reshape(coef(1,:),4,4));S(4,4)=0;initial.M=U*S*V';
for i=1:Kmode
    R=reshape(coef(2*i,:)+1i*coef(2*i+1,:),4,4);[U,S,V]=svd(R);
    initial.R(:,:,i)=U(:,1)*S(1,1)*V(:,1)';
end
info.final_linear_refit='Exact sampled basis, D=0, then rank projections';
info.hold_fundamental_removal='Initialization only, not final prediction';
end
