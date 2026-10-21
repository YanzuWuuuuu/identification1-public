import numpy as np, sys
exec(open('fd_cl.py').read().split("N=int")[0])
N=20000; k=np.arange(1,N//2); z=np.exp(2j*np.pi*k/N); f=k/N
A=zpoly(a_true,z); B=zpoly(b_true,z); C=zpoly(cn,z)/zpoly(cd,z); F=C/(1+C*B/A)
g=np.stack([-z**-1*B/A**2,-z**-2*B/A**2,z**-1/A,z**-2/A],1)
# second derivatives G_lm
zi=[z**-1,z**-2]; H=np.zeros((len(z),4,4),complex)
for i in range(2):
    for j in range(2):
        H[:,i,j]=2*zi[i]*zi[j]*B/A**3            # d2G/da_i da_j
    for j in range(2):
        H[:,i,2+j]=H[:,2+j,i]=-zi[i]*zi[j]/A**2   # d2G/da_i db_j
w=((f>0.02)&(f<0.12))*1.0
for sv in [0.15,0.5]:
    Phiu=np.abs(F)**2*(1+sv**2); Pvu=-sv**2*np.conj(F)        # E[V U*]
    b=np.real((g.conj()*(w*Pvu)[:,None]).sum(0)); M=np.real((g.conj()*(w*Phiu)[:,None]).T@g)
    M0=np.real((g.conj()*(w*np.abs(F)**2)[:,None]).T@g)
    Cm=np.real((H.conj()*(w*Pvu)[:,None,None]).sum(0))
    print(sv,'M^-1b',np.round(np.linalg.solve(M,b),5),'(M-C)^-1b',np.round(np.linalg.solve(M-Cm,b),5),'M0^-1b',np.round(np.linalg.solve(M0,b),5))
