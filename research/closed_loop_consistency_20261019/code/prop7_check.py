"""Prop 7 check: if C1 is fixed to zero (real modes, displacement output), are all 7 delays first-order identifiable?"""
import numpy as np, json
from delay_ident_check import build, resid_frac, rng
def drop_c1_cols(X, nmodes, rank1):
    # free residues: per mode 32 cols ordered (b0,b1) per (j,m) -> drop odd ones; rank1: per mode [k0,k1,psi/phi x8] -> drop k1
    keep = []; idx = 0
    for i in range(nmodes):
        if rank1 is None:
            for q in range(32): keep.append(idx+q) if q % 2 == 0 else None
            idx += 32
        else:
            for q in range(10): keep.append(idx+q) if q != 1 else None
            idx += 10
        keep += [idx, idx+1]; idx += 2          # pole columns
    keep += list(range(idx, X.shape[1]))          # D0/D2
    return X[:, keep]
out = {}
for label, f, modes, bg, r1flag in [
    ('H300 free C0, C1=0, no background', np.arange(220, 420, 0.5), [(2*np.pi*282, .0195), (2*np.pi*333, .0065), (2*np.pi*359, .0025)], 0.0, False),
    ('H300 rank-1, C1=0, no background', np.arange(220, 420, 0.5), [(2*np.pi*282, .0195), (2*np.pi*333, .0065), (2*np.pi*359, .0025)], 0.0, True),
    ('LOW 3-25 free C0, C1=0, background 10%', np.arange(3, 25, 1/6), [(2*np.pi*10.4, .049), (2*np.pi*11.5, .025), (2*np.pi*15.4, .035)], 0.1, False),
    ('LOW 3-25 rank-1, C1=0, background 10%', np.arange(3, 25, 1/6), [(2*np.pi*10.4, .049), (2*np.pi*11.5, .025), (2*np.pi*15.4, .035)], 0.1, True)]:
    tau = 300e-6 + 20e-6*rng.standard_normal(4)[:, None] + 20e-6*rng.standard_normal(4)[None]
    if r1flag:
        r1 = [(rng.standard_normal(4), rng.standard_normal(4), rng.standard_normal(), 0.0) for _ in modes]
        C0 = [r[2]*np.outer(r[0], r[1]) for r in r1]; C1 = [0*C for C in C0]
    else:
        r1 = None; C0 = [rng.standard_normal((4, 4)) for _ in modes]; C1 = [np.zeros((4, 4)) for _ in modes]
    D0 = bg*rng.standard_normal((4, 4)); D2 = bg*rng.standard_normal((4, 4))
    X, Dl = build(f, modes, C0, C1, D0, D2, tau, r1, False)
    Xr = drop_c1_cols(X, len(modes), r1)
    res = resid_frac(Xr, Dl); out[label] = res
    print(f"{label:45s}", {k: f"{v:.3f}" for k, v in res.items()}, flush=True)
json.dump(out, open('prop7_check_results.json', 'w'), indent=1)
