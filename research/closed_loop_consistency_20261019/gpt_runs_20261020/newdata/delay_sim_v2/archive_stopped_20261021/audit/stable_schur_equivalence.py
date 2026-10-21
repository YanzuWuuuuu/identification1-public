"""Independent equivalence checks only; does not run the new optimizer."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from stable_profile import *

data=dict(np.load(ROOT/"mc/data/P20_r000.npz"))
b=np.array(rd(ROOT/"audit/wide_truth_profile.json")[0]["beta"])
m=SimModel(550.)
pr=StableProfileIV(m,select(data,"W","M3"),"M3")
nl=m.nl
S,meta=projected_schur(pr,b)
A,g,score,ee,zz=pr.eq(b)
sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
As=A*sc[:,None]*sc[None,:]
block=la.solve(As[:nl,:nl],As[:nl,nl:],assume_a="gen")
Ss=As[nl:,nl:]-As[nl:,:nl]@block
Snormal=Ss/sc[nl:,None]/sc[None,nl:]
uh=pr.uh(b)
_,ZL,ZN=m.full.compact(pr.c,b,uh,True)
fz=[tall_factor(real_rows(ZL[:,o])) for o in range(4)]
yr=[real_rows(pr.Y[:,o]) for o in range(4)]
Znr=[real_rows(ZN[:,o]) for o in range(4)]

def frozen_g(q):
    bb=b.copy();bb[nl:]=q
    _,XL,_=m.full.compact(pr.c,bb,pr.U,False)
    for o in range(4):
        fx=tall_factor(real_rows(XL[:,o]))
        cross=fz[o]["Q"].T@fx["Q"]
        gamma=la.solve(cross,fz[o]["Q"].T@yr[o])
        bb[o*m.full.L:(o+1)*m.full.L]=coefficients(fx,gamma)
    yh=m.output(pr.c,bb,pr.U,False)[0]
    return sum(Znr[o].T@real_rows(pr.Y[:,o]-yh[:,o]) for o in range(4))

rows=[];rng=np.random.default_rng(92501)
for kind in ("frequency","damping","delay","mixed"):
    v=np.zeros(m.nq)
    sl=slice(0,8) if kind=="frequency" else slice(8,16) if kind=="damping" else slice(16,23) if kind=="delay" else slice(0,23)
    v[sl]=rng.normal(size=len(v[sl]));v/=la.norm(v)
    for h in (1e-4,1e-5):
        fd=(frozen_g(b[nl:]+h*v)-frozen_g(b[nl:]-h*v))/(2*h)
        target=-S@v
        rows.append(dict(direction=kind,h=h,relative=float(la.norm(fd-target)/la.norm(target)),
                         normalized_relative=float(la.norm((fd-target)/np.sqrt(zz[nl:]*ee))/la.norm(target/np.sqrt(zz[nl:]*ee)))))
result=dict(source="W P20 r000, fixed truth-q profiled point; original data and fixed instruments",
    full_linear_rank=meta["full_linear_rank"],
    direct_normal_schur_relative=float(la.norm(S-Snormal)/la.norm(S)),
    single_double_projection_relative=meta["double_vs_single_projection_relative"],
    max_projected_J_orthogonality=max(x["projected_J_orthogonality"] for x in meta["outputs"]),
    max_projected_Z_orthogonality=max(x["projected_Z_orthogonality"] for x in meta["outputs"]),
    derivative_checks=rows,
    note="Finite differences freeze Z_L and Z_N. This verifies the stated GN approximation, not the full updated-instrument moment derivative.")
js(ROOT/"audit/stable_schur_equivalence.json",result)
print(result,flush=True)
