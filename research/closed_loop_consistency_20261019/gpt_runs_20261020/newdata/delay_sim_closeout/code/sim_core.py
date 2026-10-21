"""Read-only reuse of the existing exact-ZOH modal and profiled IV code."""
import os
import sys
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
os.environ.update(NOFLEX="0", FC="60", AL="30", AH="100", SW=".03",
                  FSPR=".01", ZSPR=".4", FLP="120", FSC="3")
sys.dont_write_bytecode = True
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
NEW = ROOT.parent
RESEARCH = NEW.parents[1]
sys.path.insert(0, str(NEW / "reparam_tests/code"))
from profile_core import old, ProfileIV
sys.path.insert(0, str(RESEARCH / "code"))
import cl2 as cl
import numpy as np
import scipy.linalg as la
from scipy import signal
import time

js, rd = old.js, old.rd
FS, P0 = 5000., 30000
NI = np.array([0, 1, 0, 1])
NO = np.array([1, 0, 1, 0])
TRUTH_TAU = (NO[:, None] + NI[None, :]) / FS
BANDS = {"L": (3., 100.), "H": (100., 550.), "W": (3., 550.)}


class SimModel:
    """All eight original cl2 modes; free C0,C1,D0,D2 and modal transients."""
    def __init__(self, high=550.):
        self.full = old.DelayModel(8, np.full((4, 4), 200e-6), high)
        self.nm = 8
        self.nl = self.full.iw
        self.nq = 23
        self.nb = self.full.nb
        self.nd = 7
        self.nshape = 0
        self.iw, self.iz = self.nl, self.nl+8
        self.delay_slice = slice(self.nl+16, self.nb)
        # Non-overlapping modal boxes. Unobserved poles are not fixed at truth.
        lower = np.array([8, 11.1, 13, 240, 310, 420, 650, 737.])
        upper = np.array([11.1, 13, 20, 310, 400, 600, 737, 850.])
        self.lb = np.r_[np.log(2*np.pi*lower), np.full(8, np.log(1e-5)), np.full(7, -100.)]
        self.ub = np.r_[np.log(2*np.pi*upper), np.full(8, np.log(.3)), np.full(7, 100.)]
        self.eye = np.eye(self.nb)

    def expand(self, b, deriv=True):
        return b, self.eye if deriv else None

    def resp(self, c, b):
        return self.full.resp(c, b)

    def output(self, c, b, U, deriv=True):
        return self.full.output(c, b, U, deriv)

    def delay(self, b):
        return self.full.delay(b)

    def modal(self, b):
        return self.full.modal(b)

    def truth(self):
        b = np.zeros(self.nb)
        for o in range(4):
            for j in range(8):
                b[o*self.full.L+10*j:o*self.full.L+10*j+4] = cl.RES[j][o]
        b[self.iw:self.iz] = np.log(2*np.pi*cl.FM)
        b[self.iz:self.nl+16] = np.log(cl.ZM)
        d = TRUTH_TAU-self.full.tau
        b[self.nl+16] = d.mean()/old.UNIT
        b[self.nl+17:self.nl+20] = old.HELM.T@d.mean(1)/old.UNIT
        b[self.nl+20:] = old.HELM.T@d.mean(0)/old.UNIT
        return b


class ChunkProfileIV(ProfileIV):
    """Same equations; avoid a dense all-frequency 4 x 375 Jacobian."""
    def eq(self, b):
        uh = self.uh(b)
        A, g, score, ee = old.equations(self.m.full, self.c, self.U, self.Y, uh, b)
        if self.method == "M1":
            zz = np.diag(A).copy()
        else:
            zz = np.zeros(self.m.nb)
            f = self.m.full
            for start in range(0, len(uh), 1800):
                sl = slice(start, start+1800)
                cc = {k: v[sl] for k, v in self.c.items()}
                _, zl, zn = f.compact(cc, b, uh[sl])
                for o in range(4):
                    zz[o*f.L:(o+1)*f.L] += np.sum(abs(zl[:, o])**2, axis=0)
                zz[f.iw:] += np.sum(abs(zn)**2, axis=(0, 1))
        return A, g, score, ee, zz


def plant_with_delays():
    """Copied exact delay-state construction from delay_check/run_sim.py."""
    A, B, C = cl.ct_ss()
    Ad, Bd, Cd, _, _ = signal.cont2discrete((A, B, C, np.zeros((4, 4))), 1/FS)
    AK, BK, CK, DK = cl.controller()
    np_, nk = len(A), len(AK)
    n = np_+sum(NI)+sum(NO)+nk
    ii, oo, off = [], [], np_
    for nn in NI:
        ii.append(np.arange(off, off+nn)); off += nn
    for nn in NO:
        oo.append(np.arange(off, off+nn)); off += nn
    def update(x, e):
        raw = Cd@x[:np_]
        y = np.array([x[ix[-1]] if len(ix) else raw[j] for j, ix in enumerate(oo)])+e[4:8]
        u = e[:4]-CK@x[-nk:]-DK@y
        uw = u+e[8:]
        up = np.array([x[ix[-1]] if len(ix) else uw[j] for j, ix in enumerate(ii)])
        xn = np.zeros(n)
        xn[:np_] = Ad@x[:np_]+Bd@up
        xn[-nk:] = AK@x[-nk:]+BK@y
        for j, ix in enumerate(ii):
            if len(ix): xn[ix] = np.r_[uw[j], x[ix[:-1]]]
        for j, ix in enumerate(oo):
            if len(ix): xn[ix] = np.r_[raw[j], x[ix[:-1]]]
        return xn, np.r_[y, u]
    AA, BB, CC, DD = np.zeros((n, n)), np.zeros((n, 12)), np.zeros((8, n)), np.zeros((8, 12))
    for i in range(n):
        AA[:, i], CC[:, i] = update(np.eye(n)[i], np.zeros(12))
    for i in range(12):
        BB[:, i], DD[:, i] = update(np.zeros(n), np.eye(12)[i])
    return AA, BB, CC, DD


def ktf(f):
    AK, BK, CK, DK = cl.controller()
    z = np.exp(2j*np.pi*np.asarray(f)/FS)
    return CK[None]@np.linalg.solve(z[:, None, None]*np.eye(len(AK))-AK, np.broadcast_to(BK, (len(z),)+BK.shape))+DK


def generate(P, seed, noise=True):
    """Same noise and phase realization rule as delay_check, with variable P."""
    rng = np.random.default_rng(seed)
    kk = np.arange(1, P0//2)
    dp = np.zeros((P0//2+1, 4), complex)
    dp[kk] = .2*np.sqrt(P0*cl.dshape(np.exp(2j*np.pi*kk/P0)))[:, None]*np.exp(2j*np.pi*rng.random((len(kk), 4)))
    d = np.tile(np.fft.irfft(dp, n=P0, axis=0), (P+2, 1))
    nn = len(d)
    v = rng.normal(size=(nn, 4))*3e-7*noise
    w = rng.normal(size=(nn, 4))*.03*noise
    AA, BB, CC, DD = plant_with_delays()
    radius = float(max(abs(la.eigvals(AA))))
    if radius >= 1: raise ValueError(("unstable", radius))
    # scipy's exact state recursion; no fractional-delay approximation.
    _, out, _ = signal.dlsim((AA, BB, CC, DD, 1/FS), np.c_[d, v, w])
    N = P*P0
    k = np.arange(int(3*N/FS), int(550*N/FS)+1)
    a = dict(k=k, f=k*FS/N, N=N, P=P, seed=seed, radius=radius)
    for name, x in (("D", d), ("Y", out[:, :4]), ("U", out[:, 4:])):
        a[name] = np.fft.rfft(x[2*P0:], axis=0)[k]
    return a


def select(a, band, method):
    lo, hi = BANDS[band]
    ix = (a["f"] >= lo-1e-10) & (a["f"] <= hi+1e-10)
    if method == "M3": ix &= a["k"] % a["P"] == 0
    b = {s: a[s][ix] for s in ("f", "k", "D", "U", "Y")}
    b.update(line_indices=np.arange(len(b["f"])), K=ktf(b["f"]))
    return b


def pack(b, info, m):
    return dict(beta=b, solver=info, delay=old.delay_values(m.delay(b)),
                frequency_hz=m.modal(b)[0], damping_percent=m.modal(b)[1],
                error_percent=100*np.expm1(b[m.nl:m.nl+16]-m.truth()[m.nl:m.nl+16]),
                common_delay_error_us=(m.delay(b).mean()-TRUTH_TAU.mean())*1e6)


def pilot(P=20, band="L", method="M3", seed=1002, maxouter=5):
    start = time.perf_counter()
    a = generate(P, seed)
    m = SimModel(BANDS[band][1])
    pr = ChunkProfileIV(m, select(a, band, method), method)
    b0 = m.truth()
    # Diagnostic deliberately starts at the true poles/delays: favourable test.
    b, ok, _, inner = pr.profile(b0[m.nl:], b0[:m.nl])
    init = pack(b, dict(success=ok, inner_iterations=inner), m)
    b, info = pr.solve(b, maxouter=maxouter, hybrid=False)
    result = dict(P=P, band=band, method=method, seed=seed,
                  initialization="truth nonlinear parameters; all nuisance parameters re-estimated",
                  initial=init, endpoint=pack(b, info, m), seconds=time.perf_counter()-start)
    js(ROOT/"pilots"/f"{band}_{method}_P{P}_{seed}.json", result)
    print(band, method, info["status"], info["score"], result["seconds"], flush=True)
    return result


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="L"); ap.add_argument("--method", default="M3")
    ap.add_argument("--periods", type=int, default=20); ap.add_argument("--maxouter", type=int, default=5)
    args = ap.parse_args()
    pilot(args.periods, args.band, args.method, maxouter=args.maxouter)
