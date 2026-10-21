"""Solve the population (N->inf) IV equation h(beta)=0 for the direct method (Theorem 1, exact limit)."""
import numpy as np, os
from cl2 import *
import predict2 as P2
N=int(os.environ['NN']); k,f,zk,band,full,bt,G0,Kz,wg=P2.setup(N); I4=np.eye(4)
sd,sv=0.2,3e-7
sel=band; zs=zk[sel]; G0s=G0[sel]; Ks=Kz[sel]
S=np.linalg.inv(I4[None]+Ks@G0s); F=S@Ks
Pd=(sd**2*dshape(zs))[:,None,None]; Pv=sv**2*I4[None]+SW**2*G0s@G0s.conj().transpose(0,2,1)
Phiu=Pd*S@S.conj().transpose(0,2,1)+F@Pv@F.conj().transpose(0,2,1)
EUV=-F@Pv          # E[U V^H]
for wname in ['W1','Wg']:
    W=np.ones((sel.sum(),4)) if wname=='W1' else wg[sel]
    beta=bt[:P2.NP].copy(); full_b=np.zeros(NB)
    for it in range(300):
        full_b[:P2.NP]=beta
        Gl=P2.gl_tensor(zs,full_b)                         # K x 4 x 4 x P
        Gb=resp(zs,full_b); dG=G0s-Gb
        # h_l = sum Re tr( G_l^H W dG Phiu ) + sum Re tr( G_l^H W E[V U^H] )   (E[V U^H] = (E[U V^H])^H)
        WdG=W[:,:,None]*dG; X=WdG@Phiu + W[:,:,None]*EUV.conj().transpose(0,2,1)
        h=np.real(np.einsum('kijp,kij->p',Gl.conj(),X))
        L=psq=P2.psqrt(Phiu); Q=np.einsum('ki,kijp,kjl->kilp',np.sqrt(W),Gl,L).reshape(-1,P2.NP)
        M=np.real(Q.conj().T@Q); sc=1/np.sqrt(np.diag(M))
        st=sc*np.linalg.solve(M*sc[:,None]*sc[None,:],h*sc)
        mx=np.abs(st[IW:ITR]).max()
        if mx>0.2: st*=0.2/mx
        beta=beta+st
        if np.abs(st[IW:ITR]).max()<1e-10: break
    d=beta[IW:ITR]-bt[IW:ITR]
    print(wname,'iters',it,'last step',np.abs(st[IW:ITR]).max(),'exact limit f%',np.round((np.exp(d[:NM])-1)*100,3),' z%',np.round((np.exp(d[NM:])-1)*100,2))
