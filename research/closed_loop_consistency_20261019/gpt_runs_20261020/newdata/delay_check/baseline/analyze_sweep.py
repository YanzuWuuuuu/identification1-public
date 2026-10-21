import numpy as np, glob, os, sys
os.environ.setdefault('NOFLEX','1')
from cl2 import Cr, Br
import scipy.linalg as sla
ny0 = sla.null_space(Cr.T)[:,0]; nu0 = sla.null_space(Br)[:,0]
def angles(lin):
    C0 = [lin[32*j:32*j+16].reshape(4,4) for j in range(3)]
    Uo,_,_ = np.linalg.svd(np.hstack(C0)); ny = Uo[:,3]
    _,_,Vt = np.linalg.svd(np.vstack(C0)); nu = Vt[3]
    a = lambda x,y: np.degrees(np.arccos(min(1,abs(x@y)/np.linalg.norm(x)/np.linalg.norm(y))))
    return a(ny,ny0), a(nu,nu0)
for tag, files in [('N=90k',['sw_090k.npz']),('N=300k',['sw_300k.npz']),('N=900k',['sw_900k_a.npz','sw_900k_b.npz'])]:
    fs_=[f for f in files if os.path.exists(f)]
    if not fs_: continue
    r=[np.load(f) for f in fs_]
    print('=====',tag)
    for m in r[0].files:
        v=np.vstack([x[m] for x in r]); ok=np.isfinite(v[:,0]); v=v[ok]; n=len(v)
        z=(np.exp(v[:,3:6])-1)*100; fr=(np.exp(v[:,0:3])-1)*100
        ang=np.array([angles(row[7:]) for row in v])
        print(f"{m:12s} n={n} fail={(~ok).sum()} maxit={int(v[:,6].max())}")
        print(f"   zeta% mean {np.round(z.mean(0),2)} se {np.round(z.std(0)/np.sqrt(n),2)} sd {np.round(z.std(0),2)} rmse {np.round(np.sqrt((z**2).mean(0)),2)}")
        print(f"   f%    mean {np.round(fr.mean(0),3)} sd {np.round(fr.std(0),3)}")
        print(f"   null-space angle deg (out,in): mean {np.round(ang.mean(0),3)} rms {np.round(np.sqrt((ang**2).mean(0)),3)}")
