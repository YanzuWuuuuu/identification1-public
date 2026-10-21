import numpy as np, sys, time, os
from cl2 import *
from predict2 import setup

N = int(os.environ.get('NN', '30000')); NPER_ = int(os.environ.get('NPER', '10'))
SD = float(sys.argv[1]); SV = float(sys.argv[2]); NMC = int(sys.argv[3]); SEED = int(sys.argv[4]); OUT = sys.argv[5]
methods = sys.argv[6].split(',')
k, f, zk, band, full, bt, G0, Kz, wg = setup(N)
I4 = np.eye(4)

def ref_instr(zs, Dsel, Ksel):
    def Uh(beta):
        G = resp(zs, beta)
        return np.einsum('kij,kj->ki', np.linalg.inv(I4[None] + Ksel @ G), Dsel)
    return Uh

rng = np.random.default_rng(SEED)
res = {m: [] for m in methods}
for r in range(NMC):
    t0 = time.time()
    d, u, y = simulate_periodic(N//NPER_, NPER_, SD, SV, rng)
    D = np.fft.fft(d, axis=0)[k]; U = np.fft.fft(u, axis=0)[k]; Y = np.fft.fft(y, axis=0)[k]
    pert = rng.standard_normal(2*NM)
    Tfir = fir_tud(d, u) if any(m.startswith('fir') for m in methods) else None
    for m in methods:
        sel = full if '_full_' in m else band
        W = np.ones((sel.sum(), 4)) if m.endswith('W1') else wg[sel]
        zs = zk[sel]
        b0 = bt.copy(); b0[IW:IW+NM] += np.log(1 + 0.003*pert[:NM]); b0[IZ:IZ+NM] += np.log(1 + 0.15*pert[NM:])
        b0[:NLIN] = 0; b0[ITR:] = 0
        b0 = init_lin(zs, U[sel], Y[sel], W, b0)
        if m.startswith('dir'): Uh = U[sel]
        elif m.startswith('ref'): Uh = ref_instr(zs, D[sel], Kz[sel])
        else: Uh = np.einsum('kij,kj->ki', Tfir(zs), D[sel])
        try:
            b, it = iv_solve(zs, U[sel], Y[sel], Uh, W, b0, it=80, tol=1e-5)
        except Exception as e:
            b, it = b0*np.nan, -1
        res[m].append(np.r_[b[IW:ITR] - bt[IW:ITR], it])
    print(r, round(time.time()-t0, 1), flush=True)
    np.savez(OUT, **{m: np.array(v) for m, v in res.items()})
