import numpy as np, sys, os
TR=os.environ.get('TR','1')=='1'
from scipy.signal import lfilter
rng=np.random.default_rng(1)
# true DT plant: lightly damped 2nd order, strictly proper (b0=0)
fs=1.0; w0=2*np.pi*0.05; zeta=0.05; r_=np.exp(-zeta*w0)
a_true=np.array([1,-2*r_*np.cos(w0*np.sqrt(1-zeta**2)),r_**2])
b_true=np.array([0,1.0,0.5])*0.01
# controller C = kc*(1-0.9 z^-1)/(1-z^-1)
cn=np.array([0.8,-0.75]); cd=np.array([1,-1.0])
pm=np.convolve
# closed loop: den = A*Cd + B*Cn
den=np.polyadd(pm(a_true,cd),pm(b_true,cn))
assert np.all(np.abs(np.roots(den))<1), np.abs(np.roots(den))
def sim(N,sv):
    r=rng.standard_normal(N); v=sv*rng.standard_normal(N)
    y=lfilter(pm(b_true,cn),den,r)+lfilter(pm(a_true,cd),den,v)
    u=lfilter(pm(a_true,cn),den,r)-lfilter(pm(a_true,cn),den,v)
    return r,u,y
def zpoly(c,z):  # c[0]+c[1]z^-1+...
    return sum(ci*z**(-i) for i,ci in enumerate(c))
def fd_srivc(R,U,Y,z,w,mode,it=30):
    a=np.array([1,-1.8,0.9]); b=np.array([0,0.01,0.0])
    for _ in range(it):
        A=zpoly(a,z); B=zpoly(b,z)
        Yf=Y/A; Uf=U/A
        if mode=='ref':   # CLSRIVC-like: noise-free input from reference, known C
            C=zpoly(cn,z)/zpoly(cd,z); Uh=C/(1+B/A*C)*R
        else: Uh=U
        Xh=B/A*Uh
        T=[z**0/A, z**-1/A, z**-2/A] if TR else []
        phi=np.stack([-z**-1*Yf,-z**-2*Yf,z**-1*Uf,z**-2*Uf]+T,1)
        zet=np.stack([-z**-1*Xh/A,-z**-2*Xh/A,z**-1*Uh/A,z**-2*Uh/A]+T,1)
        M=np.real((zet.conj()*w[:,None]).T@phi); h=np.real((zet.conj()*w[:,None]).T@Yf)
        th=np.linalg.solve(M,h); a=np.r_[1,th[:2]]; b=np.r_[0,th[2:4]]
    return np.r_[th[:4]]
N=int(sys.argv[1]) if len(sys.argv)>1 else 20000; sv=float(sys.argv[2]) if len(sys.argv)>2 else 0.5; nmc=int(sys.argv[3]) if len(sys.argv)>3 else 100
k=np.arange(1,N//2); om=2*np.pi*k/N; z=np.exp(1j*om); f=k/N
Gt=zpoly(b_true,z)/zpoly(a_true,z)
band=(f>0.02)&(f<0.12)
cases={'full,uniform':(np.ones_like(f),'dir'),
       'band 0.02-0.12,uniform':(band*1.0,'dir'),
       'full,W=1/|G|^2 (Schur-like)':(1/np.abs(Gt)**2/np.mean(1/np.abs(Gt)**2),'dir'),
       'band,reference-IV':(band*1.0,'ref'),
       'full,reference-IV, W=1/|G|^2':(1/np.abs(Gt)**2/np.mean(1/np.abs(Gt)**2),'ref')}
th0=np.r_[a_true[1:],b_true[1:]]
res={c:[] for c in cases}
for m in range(nmc):
    r,u,y=sim(N,sv)
    R=np.fft.fft(r)[k];U=np.fft.fft(u)[k];Y=np.fft.fft(y)[k]
    for c,(w,mode) in cases.items(): res[c].append(fd_srivc(R,U,Y,z,w,mode))
def modal(th):
    p=np.roots(np.r_[1,th[:2]])[0]; s=np.log(p); return abs(s)/2/np.pi, -s.real/abs(s)
f0,z0=modal(th0)
print(f"N={N} sv={sv} nmc={nmc}  true f={f0:.5f} zeta={z0:.4f}")
for c,v in res.items():
    v=np.array(v); fz=np.array([modal(t) for t in v])
    bf=(fz[:,0].mean()/f0-1)*100; sf=fz[:,0].std()/f0*100
    bz=(fz[:,1].mean()/z0-1)*100; sz=fz[:,1].std()/z0*100
    bb=(v[:,2:].mean(0)/th0[2:]-1)*100
    print(f"{c:32s} f bias {bf:+.3f}% (sd {sf:.3f}, se {sf/np.sqrt(nmc):.3f})  zeta bias {bz:+.2f}% (sd {sz:.2f}, se {sz/np.sqrt(nmc):.2f})  b1 bias {bb[0]:+.2f}%")
