"""Inspect a saved, not-yet-accepted stable-profile endpoint; no fitting."""
from implementation_check import *

p=ROOT/"stable_profile_pilot/v2_P10_r000_H_M3/eval_040.json"
r=rd(p);b=np.array(r["beta"])
data=dict(np.load(ROOT/"mc/data/P10_r000.npz"))
m=SimModel(550.);pr=QRProfileIV(m,select(data,"H","M3"),"M3")
A,g,score,ee,zz=pr.eq(b);ss=abs(g)/np.sqrt(zz*ee)
y,XL,_=m.full.compact(pr.c,b,pr.U)
comp=XL*b[:m.nl].reshape(1,4,m.full.L)
cancel=float(la.norm(np.sum(abs(comp),axis=2))/la.norm(y))
previous=next(x for x in rd(ROOT/"audit/implementation_check.json") if x["band"]=="H" and x["at"]=="endpoint")
sort=np.argsort(ss)[::-1][:10]
f,z=m.modal(b);bt=m.truth();ft,zt=m.modal(bt)
result=dict(source=str(p),valid_estimate=False,gate="375-dimensional score < 1e-6",
    score=score,linear_score=float(max(ss[:m.nl])),
    highest_score=[dict(index=int(i),parameter=label(int(i),m),score=float(ss[i]),signed_moment=float(g[i])) for i in sort],
    frequency_hz=f,frequency_truth_hz=ft,frequency_error_percent=100*np.expm1(b[m.iw:m.iz]-bt[m.iw:m.iz]),
    damping_percent=z,damping_truth_percent=zt,damping_error_percent=100*np.expm1(b[m.iz:m.delay_slice.start]-bt[m.iz:m.delay_slice.start]),
    delay=old.delay_values(m.delay(b)),
    cancellation_factor=cancel,previous_QR_cancellation_factor=previous["cancellation_factor"],
    cancellation_reduction_factor=previous["cancellation_factor"]/cancel,
    max_abs_linear=float(max(abs(b[:m.nl]))),
    output_relative_residual=float(la.norm(pr.Y-y)/la.norm(pr.Y)),
    notes="Diagnostic of eval_040 only. Endpoint fails original score and is not a valid parameter estimate.")
js(ROOT/"audit/stable_endpoint_check.json",result)
print(result)
