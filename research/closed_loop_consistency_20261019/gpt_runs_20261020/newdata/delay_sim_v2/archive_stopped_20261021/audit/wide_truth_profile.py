"""One saved P20 record; truth/initial and fixed-q profile, no outer fit."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from qr_profile import *

data=dict(np.load(ROOT/"mc/data/P20_r000.npz"))
m=SimModel(550.)
rows=[]
for init in ("truth","mc"):
    b0=initial(m,init)
    pr=QRProfileIV(m,select(data,"W","M3"),"M3")
    A,g,s,e,z=pr.eq(b0)
    row=dict(initial=init,initial_score=s,initial_linear_score=float(max(abs(g[:m.nl])/np.sqrt(z[:m.nl]*e))))
    b,ok,_,n=pr.profile(b0[m.nl:],b0[:m.nl],maxinner=8)
    A,g,s,e,z=pr.eq(b)
    _,LL,_=m.full.compact(pr.c,b,pr.U,False)
    comp=LL*b[:m.nl].reshape(1,4,m.full.L)
    row.update(profile_success=ok,iterations=n,profile_full_score=s,
               profile_linear_score=float(max(abs(g[:m.nl])/np.sqrt(z[:m.nl]*e))),
               profile_linear_meta=pr.last_qr_profile,
               cancellation_factor=float(la.norm(np.sum(abs(comp),axis=2))/la.norm(comp.sum(2))),
               max_abs_linear=float(max(abs(b[:m.nl]))),
               beta=b)
    rows.append(row)
    js(ROOT/"audit/wide_truth_profile.json",rows)
    print(init,row["initial_score"],row["profile_linear_score"],row["profile_full_score"],ok,row["cancellation_factor"],flush=True)
