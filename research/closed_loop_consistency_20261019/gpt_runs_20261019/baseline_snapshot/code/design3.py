import numpy as np, sys
from cl2 import *
N=30000; k=np.arange(1,N//2); f=k*fs/N; zk=np.exp(2j*np.pi*k/N)
bt=truth_beta(); G0=resp(zk,bt); Kz=Ktf(zk); S=np.linalg.inv(np.eye(4)[None]+Kz@G0); F=S@Kz
Tu=np.linalg.pinv(Br); T=np.vstack([np.linalg.pinv(Tu), sla.null_space(Tu.T).T])
pd=dshape(zk); sd=float(sys.argv[1]); sv=float(sys.argv[2])
bands=[(3,10),(10,25),(25,80),(80,170)]
for a,b in bands:
    sl=(f>=a)&(f<b); G=G0[sl]
    Pv=sv**2*np.eye(4)[None]+SW**2*G@G.conj().transpose(0,2,1)
    TS=T[None]@S[sl]; TF=T[None]@F[sl]
    nu=np.real(np.einsum('kij,kjl,kil->ki',TF,Pv,TF.conj())).sum(0)
    su=(np.abs(TS)**2*pd[sl][:,None,None]).sum(2).sum(0)*sd**2
    print(a,b,'share % per coord',np.round(100*nu/(nu+su),2))
Ty=G0@S
for a,b in bands:
    sl=(f>=a)&(f<b); G=G0[sl]; Pv=sv**2*np.eye(4)[None]+SW**2*G@G.conj().transpose(0,2,1)
    Sy=np.eye(4)[None]-G@F[sl]; nv=np.real(np.trace(Sy@Pv@Sy.conj().transpose(0,2,1),axis1=1,axis2=2)).sum()
    print(a,b,'y SNR dB',round(10*np.log10(np.sum(np.abs(Ty[sl])**2*pd[sl][:,None,None])*sd**2/nv),1))
print('rho',closed_loop()[-1])
