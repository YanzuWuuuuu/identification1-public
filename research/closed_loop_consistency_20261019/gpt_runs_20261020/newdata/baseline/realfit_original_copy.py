"""Real-data (old z=0) low-frequency support-mode identification: direct FD-IV vs reference IV.
Model in band [flo,fhi]:  G(s) = sum_{i=1..NMODE} [C0_i w_i^2 + C1_i w_i s]/(s^2+2 z_i w_i s + w_i^2) + D0 + D2 (s/wr)^2,
multiplied element-wise by fixed effective channel delays exp(-s tau_jm) (from earlier delay identification),
plus per-mode transient terms (t0 + t1 z^-1)/a_d,i(z^-1) (ZOH denominators of the modes).
"""
import numpy as np, scipy.linalg as sla
fs = 5000.0; Ts = 1/fs
NMODE = 3; WR = 2*np.pi*25.0
NLIN = 32*NMODE + 32; IW = NLIN; IZ = IW + NMODE; ITR = IZ + NMODE; NB = ITR + 8*NMODE

def ad_zoh(w, z_):
    A = np.array([[0, 1], [-w**2, -2*z_*w]]); Ad = sla.expm(A*Ts)
    return np.array([1, -np.trace(Ad), np.linalg.det(Ad)])

def unpack(b):
    lin = b[:32*NMODE].reshape(NMODE, 2, 4, 4); D0 = b[32*NMODE:32*NMODE+16].reshape(4, 4)
    D2 = b[32*NMODE+16:NLIN].reshape(4, 4)
    return lin, D0, D2, b[IW:IZ], b[IZ:ITR], b[ITR:].reshape(NMODE, 4, 2)

class Ctx:
    def __init__(s, f, tau):          # f in Hz (K,), tau (4x4) seconds
        s.f = f; s.s = 2j*np.pi*f; s.zi = np.exp(-2j*np.pi*f/fs)
        s.E = np.exp(-s.s[:, None, None]*tau[None])
    def sub(s, sl):
        c = Ctx.__new__(Ctx); c.f = s.f[sl]; c.s = s.s[sl]; c.zi = s.zi[sl]; c.E = s.E[sl]; return c

def blocks(c, lw, lz):
    out = []
    for i in range(NMODE):
        w, z_ = np.exp(lw[i]), np.exp(lz[i]); den = c.s**2 + 2*z_*w*c.s + w**2
        ad = ad_zoh(w, z_); adz = ad[0] + ad[1]*c.zi + ad[2]*c.zi**2
        out.append((w**2/den, w*c.s/den, adz))
    return out

def resp(c, b):
    lin, D0, D2, lw, lz, _ = unpack(b); bl = blocks(c, lw, lz)
    G = D0[None] + ((c.s/WR)**2)[:, None, None]*D2[None]
    for i in range(NMODE): G = G + bl[i][0][:, None, None]*lin[i, 0][None] + bl[i][1][:, None, None]*lin[i, 1][None]
    return G*c.E

def out_jac(c, b, U, need_jac=True, dl=1e-6):
    lin, D0, D2, lw, lz, tr = unpack(b); K = len(c.f); bl = blocks(c, lw, lz)
    EU = c.E*U[:, None, :]                       # K x out x in : E_jm U_m
    G = resp(c, b); Yh = np.einsum('kjm,km->kj', G, U)
    for i in range(NMODE): Yh += (tr[i][:, 0][None] + tr[i][:, 1][None]*c.zi[:, None])/bl[i][2][:, None]
    if not need_jac: return Yh, None
    J = np.zeros((K, 4, NB), complex)
    def put_lin(base, h):
        for o in range(4): J[:, o, base+4*o:base+4*o+4] = h[:, None]*EU[:, o, :]
    for i in range(NMODE):
        put_lin(32*i, bl[i][0]); put_lin(32*i+16, bl[i][1])
    put_lin(32*NMODE, np.ones(K)); put_lin(32*NMODE+16, (c.s/WR)**2)
    for i in range(NMODE):
        for which, idx in ((0, IW+i), (1, IZ+i)):
            lw2, lz2 = lw.copy(), lz.copy()
            if which == 0: lw2[i] += dl
            else: lz2[i] += dl
            b2 = blocks(c, lw2, lz2)[i]
            dG = ((b2[0]-bl[i][0])[:, None, None]*lin[i, 0][None] + (b2[1]-bl[i][1])[:, None, None]*lin[i, 1][None])*c.E/dl
            tn = tr[i][:, 0][None] + tr[i][:, 1][None]*c.zi[:, None]
            J[:, :, idx] = np.einsum('kjm,km->kj', dG, U) + (tn/b2[2][:, None] - tn/bl[i][2][:, None])/dl
    for i in range(NMODE):
        for o in range(4):
            J[:, o, ITR+8*i+2*o] = 1/bl[i][2]; J[:, o, ITR+8*i+2*o+1] = c.zi/bl[i][2]
    return Yh, J

LINIDX = np.r_[0:NLIN, ITR:NB]
def iv_solve(c, U, Y, Uh, W, b0, it=60, tol=1e-6, chunk=600, free=None, verbose=False):
    b = b0.copy(); K = len(c.f); fr = np.arange(NB) if free is None else free
    for itn in range(it):
        Uhb = Uh(b) if callable(Uh) else Uh
        A = np.zeros((len(fr), len(fr))); g = np.zeros(len(fr))
        for s0 in range(0, K, chunk):
            sl = slice(s0, min(K, s0+chunk)); cc = c.sub(sl)
            Yh, J = out_jac(cc, b, U[sl])
            Z = J if Uhb is U else out_jac(cc, b, Uhb[sl])[1]
            J = J[:, :, fr]; Z = Z[:, :, fr]
            ZW = (Z.conj()*W[sl][:, :, None]).reshape(-1, len(fr))
            A += np.real(ZW.T @ J.reshape(-1, len(fr))); g += np.real(ZW.T @ (Y[sl]-Yh).reshape(-1))
        sc = 1/np.sqrt(np.abs(np.diag(A)) + 1e-300)
        st = sc*np.linalg.solve(A*sc[:, None]*sc[None, :], g*sc)
        if free is not None:
            b[fr] += st; return b, 0
        mx = np.max(np.abs(st[IW:ITR]))
        if mx > 0.1: st *= 0.1/mx
        b = b + st
        if verbose: print(itn, np.round(np.exp(b[IW:IZ])/2/np.pi, 3), np.round(np.exp(b[IZ:ITR])*100, 2), flush=True)
        if np.max(np.abs(st[IW:ITR])) < tol: break
    return b, itn

def init_lin(c, U, Y, W, b):
    return iv_solve(c, U, Y, U, W, b, free=LINIDX)[0]

def modal(b):
    return np.exp(b[IW:IZ])/2/np.pi, np.exp(b[IZ:ITR])*100
