"""Cross-validation: estimate on periods [1,5), validate on excited lines of periods [5,9) and vice versa."""
import numpy as np, os, sys, json
import run_real as RR
from realfit import *
def dft(p0, p1):
    P = RR.P; seg = slice(p0*P, p1*P); N = (p1-p0)*P; nper = p1-p0
    k = np.arange(1, N//2); f = k*fs/N
    Uk = np.fft.fft(RR.U[seg], axis=0)[k]; Yk = np.fft.fft(RR.Y[seg], axis=0)[k]
    lines = (f >= RR.FLO) & (f <= RR.FHI) & (k % nper == 0)
    return f[lines], Uk[lines], Yk[lines]
def fit_save(p0, p1, methods):
    o = RR.estimate.__wrapped__(p0, p1, methods) if hasattr(RR.estimate, '__wrapped__') else None
out = {}
for (ea, eb), (va, vb) in [((1, 5), (5, 9)), ((5, 9), (1, 5))]:
    res = RR.estimate(ea, eb, ['dir_all', 'dir_lines', 'ref_K'], return_beta=True)
    fv, Uv, Yv = dft(va, vb); c = Ctx(fv, RR.tau); W = np.ones((len(fv), 4))
    for m, b in res['beta'].items():
        b = np.array(b); bv = b.copy(); bv[ITR:] = 0
        bv = iv_solve(c, Uv, Yv, Uv, W, bv, free=np.arange(ITR, NB))[0]   # refit only transient coefficients
        Yh, _ = out_jac(c, bv, Uv, need_jac=False)
        e = float(np.linalg.norm(Yv-Yh)/np.linalg.norm(Yv))
        out.setdefault(m, []).append(e)
        print(f"est {ea}-{eb} -> val {va}-{vb}  {m:10s} validation rel. residual on excited lines: {100*e:.2f}%", flush=True)
json.dump(out, open(f"xval_z{os.environ.get('ZPOS','0')}.json", 'w'))
