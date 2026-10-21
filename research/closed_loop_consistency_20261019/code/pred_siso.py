import numpy as np, sys
exec(open('fd_cl.py').read().split("N=int")[0])
N=20000; sv=float(sys.argv[1]); nmc=int(sys.argv[2])
k=np.arange(1,N//2); om=2*np.pi*k/N; z=np.exp(1j*om); f=k/N
A=zpoly(a_true,z); B=zpoly(b_true,z); C=zpoly(cn,z)/zpoly(cd,z)
F=C/(1+C*B/A)
g=np.stack([-z**-1*B/A**2,-z**-2*B/A**2,z**-1/A,z**-2/A],1)
th0=np.r_[a_true[1:],b_true[1:]]
def predict(w):
    Phu=np.abs(F)**2*(1+sv**2)
    b=-(w[:,None]*np.real(g*F[:,None])*sv**2).sum(0)
    M=np.real((g.conj()*(w*Phu)[:,None]).T@g)
    return np.linalg.solve(M,b)
band=((f>0.02)&(f<0.12))*1.0
for name,w in [('full',np.ones_like(f)),('band',band)]:
    d=predict(w); est=[]
    for m in range(nmc):
        r,u,y=sim(N,sv);R=np.fft.fft(r)[k];U=np.fft.fft(u)[k];Y=np.fft.fft(y)[k]
        est.append(fd_srivc(R,U,Y,z,w,'dir',40))
    est=np.array(est)-th0
    print(name,'pred dtheta',np.array2string(d,precision=5),'\n     MC mean ',np.array2string(est.mean(0),precision=5),' se',np.array2string(est.std(0)/np.sqrt(nmc),precision=5))
