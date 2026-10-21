import numpy as np, os
from cl2 import *
import predict2 as P2
N=int(os.environ['NN']); k,f,zk,band,full,bt,G0,Kz,wg=P2.setup(N); I4=np.eye(4)
sel=band; zs=zk[sel]; G=G0[sel]; Ks=Kz[sel]
S=np.linalg.inv(I4[None]+Ks@G); F=S@Ks
Pd=(0.2**2*dshape(zs))[:,None,None]; Pv=(3e-7)**2*I4[None]+SW**2*G@G.conj().transpose(0,2,1)
Phiu=Pd*S@S.conj().transpose(0,2,1)+F@Pv@F.conj().transpose(0,2,1); Pvu=-Pv@F.conj().transpose(0,2,1)
NPp=P2.NP; Gl=P2.gl_tensor(zs,bt)
for wname in ['W1','Wg']:
    W=np.ones((sel.sum(),4)) if wname=='W1' else wg[sel]
    b=np.real(np.einsum('kijp,kij->p',Gl.conj(),W[:,:,None]*Pvu))
    L=P2.psqrt(Phiu); Q=np.einsum('ki,kijp,kjl->kilp',np.sqrt(W),Gl,L).reshape(-1,NPp); M=np.real(Q.conj().T@Q)
    C=np.zeros((NPp,NPp)); h=float(os.environ.get("HH","1e-6"))
    for m in range(IW,ITR):
        bp=bt.copy(); bp[m]+=h; Glm=(P2.gl_tensor(zs,bp)-Gl)/h
        col=np.real(np.einsum('kijp,kij->p',Glm.conj(),W[:,:,None]*Pvu))
        C[:,m]=col; C[m,:]=col
    def solve(A,v):
        sc=1/np.sqrt(np.abs(np.diag(A))); return sc*np.linalg.solve(A*sc[:,None]*sc[None,:],v*sc)
    d1=solve(M,b); d2=solve(M-C,b)
    # noise-free M0
    Q0=np.einsum('ki,kijp,kjl->kilp',np.sqrt(W),Gl,np.sqrt(Pd)*S).reshape(-1,NPp); M0=np.real(Q0.conj().T@Q0)
    d0=solve(M0,b)
    pz=lambda d: np.round((np.exp(d[IZ:ITR])-1)*100,2)
    print(wname,' M^-1 b',pz(d1),' (M-C)^-1 b',pz(d2),' M0^-1 b (small-noise)',pz(d0))
