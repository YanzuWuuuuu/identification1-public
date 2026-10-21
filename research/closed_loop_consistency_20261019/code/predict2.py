"""Theorem-based prediction for cl2: first-order bias delta = M^-1 b of the direct FD-IV root,
and asymptotic std of direct / reference IV. Units: unnormalized DFT (Phi = N sigma^2)."""
import numpy as np, sys, os
from cl2 import *

NP = ITR   # non-transient parameters
def psqrt(A):
    w, V = np.linalg.eigh(A); return V*np.sqrt(np.clip(w, 0, None))[:, None, :]
NPER = int(os.environ.get('NPER', '0'))
def LINEMASK(zs):
    if NPER <= 1: return np.ones(len(zs))
    kk = np.rint(np.angle(zs)/(2*np.pi)*NN_).astype(int)
    return np.where(kk % NPER == 0, NPER, 0.0)
NN_ = int(os.environ.get('NN', '30000'))

def gl_tensor(zs, beta):
    Gl = np.zeros((len(zs), 4, 4, NP), complex)
    for mi in range(4):
        Ue = np.zeros((len(zs), 4)); Ue[:, mi] = 1
        Gl[:, :, mi, :] = out_jac(zs, beta, Ue)[1][:, :, :NP]
    return Gl

def predict(N, sd, sv, G, Kk, W, zs, beta, method, chunk=800):
    I4 = np.eye(4); Ktot = len(G); P_ = NP
    b = np.zeros(P_); M = np.zeros((P_, P_)); Cg = np.zeros((P_, P_))
    for s0 in range(0, Ktot, chunk):
        sl = slice(s0, min(Ktot, s0+chunk)); Gc = G[sl]; Kc = Kk[sl]; Wc = W[sl]; K = len(Gc)
        Gl = gl_tensor(zs[sl], beta)
        S = np.linalg.inv(I4[None] + Kc @ Gc); F = S @ Kc
        Pd = (N*sd**2*dshape(zs[sl])*LINEMASK(zs[sl]))[:, None, None]
        Pv = N*sv**2*I4[None] + N*SW**2*Gc @ Gc.conj().transpose(0, 2, 1)       # output-equivalent noise spectrum
        Phiu = Pd*S @ S.conj().transpose(0, 2, 1) + F @ Pv @ F.conj().transpose(0, 2, 1)
        GF = np.einsum('kijp,kjl->kilp', Gl, F)
        WPv = Wc[:, :, None]*Pv
        b += -np.real(np.einsum('kij,kijp->p', WPv, GF.conj()))
        Lx = psqrt(Phiu) if method == 'dir' else np.sqrt(Pd)*S
        Lm = Lx if method == 'dir' else np.sqrt(Pd)*S
        Q1 = np.einsum('ki,kijp,kjl->kilp', np.sqrt(Wc), Gl, Lx).reshape(K*16, -1)
        Q2 = Q1 if method == 'dir' else np.einsum('ki,kijp,kjl->kilp', np.sqrt(Wc), Gl, Lm).reshape(K*16, -1)
        M += np.real(Q1.conj().T @ Q2)
        Lv = psqrt(Pv)
        Qz = np.einsum('kji,kj,kjmp,kml->kilp', Lv.conj(), Wc, Gl, Lx).reshape(K*16, -1)
        Cg += 0.5*np.real(Qz.conj().T @ Qz)
    bb = b if method == 'dir' else 0*b
    sc = 1/np.sqrt(np.abs(np.diag(M))); Ms = M*sc[:, None]*sc[None, :]
    Mi = sc[:, None]*np.linalg.inv(Ms)*sc[None, :]
    P = Mi @ Cg @ Mi.T
    return Mi @ bb, np.sqrt(np.diag(P))

def setup(N, flo=None, fhi=None):
    flo = float(os.environ.get('FLO', '3')) if flo is None else flo; fhi = float(os.environ.get('FHI', '800')) if fhi is None else fhi
    k = np.arange(1, N//2); f = k*fs/N; zk = np.exp(2j*np.pi*k/N)
    band = (f >= flo) & (f <= fhi); full = np.ones_like(band)
    bt = truth_beta(); G0 = resp(zk, bt); Kz = Ktf(zk)
    wg = 1/np.sum(np.abs(G0)**2, axis=2); wg /= wg[band].mean(0, keepdims=True)
    return k, f, zk, band, full, bt, G0, Kz, wg

CASES = [('dir_full_W1', 'full', 'W1', 'dir'), ('dir_band_W1', 'band', 'W1', 'dir'),
         ('dir_band_Wg', 'band', 'Wg', 'dir'), ('ref_band_Wg', 'band', 'Wg', 'ref'), ('ref_band_W1', 'band', 'W1', 'ref')]

if __name__ == '__main__':
    N = int(os.environ.get('NN', '30000')); sd = float(sys.argv[1]); sv = float(sys.argv[2])
    k, f, zk, band, full, bt, G0, Kz, wg = setup(N)
    out = {}
    for name, sname, wname, meth in CASES:
        sel = full if sname == 'full' else band
        W = np.ones((sel.sum(), 4)) if wname == 'W1' else wg[sel]
        dl, sdv = predict(N, sd, sv, G0[sel], Kz[sel], W, zk[sel], bt, meth)
        out[name] = (dl, sdv)
        print(f"{name:12s} bias f%:", np.round(dl[IW:IW+NM]*100, 3), "\n             std  f%:", np.round(sdv[IW:IW+NM]*100, 3))
        print(f"{'':12s} bias z%:", np.round(dl[IZ:IZ+NM]*100, 2), "\n             std  z%:", np.round(sdv[IZ:IZ+NM]*100, 2), flush=True)
    np.save(f'pred_sd{sd}_sv{sv}.npy', {k_: (v[0], v[1]) for k_, v in out.items()}, allow_pickle=True)
