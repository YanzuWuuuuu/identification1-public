"""Read-only inspection of retained ranks in the current outer Schur step."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from qr_profile import *

rows=[]
for band in ("L","H","W"):
    P=20 if band=="W" else 10
    data=dict(np.load(ROOT/"mc/data"/f"P{P:02d}_r000.npz"))
    m=SimModel(BANDS[band][1]);nl=m.nl
    if band=="W":
        rr=rd(ROOT/"audit/wide_truth_profile.json")
        points=[("truth_profile",np.array(rr[0]["beta"])),("initial_profile",np.array(rr[1]["beta"]))]
    else:
        rr=rd(ROOT/"mc_qr_final/raw"/f"P10_r000_{band}_M3.json")
        points=[("endpoint",np.array(rr["endpoint"]["beta"]))]
    for point,b in points:
        pr=QRProfileIV(m,select(data,band,"M3"),"M3")
        A,g,sc,ee,zz=pr.eq(b)
        scale=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
        As=A*scale[:,None]*scale[None,:]
        x,res,rank,sv=la.lstsq(As[:nl,:nl],As[:nl,nl:],cond=1e-13,lapack_driver="gelsy")
        S=As[nl:,nl:]-As[nl:,:nl]@x
        ss=la.svdvals(S)
        aa=la.svdvals(As[:nl,:nl])
        st,_,srank,_=la.lstsq(S,g[nl:]*scale[nl:],cond=1e-12,lapack_driver="gelsy")
        st*=scale[nl:]
        cap=np.r_[np.full(16,.04),np.full(7,.4)]
        fac=min(1.,1/max(np.max(abs(st)/cap),1e-300))
        pos=st>0;neg=st<0
        if pos.any():fac=min(fac,float(np.min(.99*(m.ub[pos]-b[nl:][pos])/st[pos])))
        if neg.any():fac=min(fac,float(np.min(.99*(m.lb[neg]-b[nl:][neg])/st[neg])))
        _,X,_=m.full.compact(pr.c,b,pr.U,False)
        fx=[tall_factor(real_rows(X[:,o])) for o in range(4)]
        result=dict(band=band,point=point,full_score=sc,linear_columns=nl,
            tall_X_rank=sum(z["rank"] for z in fx),schur_linear_rank=int(rank),
            schur_linear_scaled_condition=float(aa[0]/max(aa[-1],1e-300)),
            schur_nonlinear_rank=int(srank),nonlinear_columns=23,
            schur_condition=float(ss[0]/max(ss[-1],1e-300)),
            schur_linear_rhs_relative=float(la.norm(As[:nl,:nl]@x-As[:nl,nl:])/la.norm(As[:nl,nl:])),
            raw_outer_step=st,step_factor=fac,actual_outer_step=st*fac)
        rows.append(result);js(ROOT/"audit/schur_rank_check.json",rows)
        print(band,point,"rank",result["tall_X_rank"],rank,srank,"fac",fac,flush=True)
