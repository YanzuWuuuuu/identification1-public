"""
v2: 4x4 platform with suspended rigid body (3 low-frequency support modes) + 5 flexible modes.
CT additive modal model evaluated through exact ZOH equivalents; DT controller closed loop.
Root of the frequency-domain IV equation  sum_k Re{ Z_k^H W_k (Y_k - Yhat_k(beta)) } = 0
(Z = model-output gradient evaluated with the instrument input Uhat) solved by IV-Gauss-Newton.
"""
import numpy as np, scipy.linalg as sla, scipy.signal as ss

fs = 5000.0; Ts = 1/fs
xa = np.array([0.2, -0.2, -0.2, 0.2]); ya = np.array([0.15, 0.15, -0.15, -0.15])
m, Jx, Jy = 5.0, 0.05, 0.08
Br = np.vstack([np.ones(4)/m, ya/Jx, -xa/Jy])
xs = np.array([0.18, -0.18, -0.18, 0.18]); ys = np.array([0.13, 0.13, -0.13, -0.13])
Cr = np.vstack([np.ones(4), ys, -xs]).T
fsup = np.array([10.6, 11.7, 16.1]); zsup = np.array([0.03, 0.04, 0.05])
fflex = np.array([286.573, 344.945, 490.39, 718.958, 755.501])
zflex = np.array([0.0119, 0.0132, 0.0186, 0.0322, 0.0346])
rng0 = np.random.default_rng(7)
import os
FLP = float(os.environ.get('FLP', '120')); FSC = float(os.environ.get('FSC', '3'))
bm = rng0.normal(size=(5, 4))*0.6/np.sqrt(m)*FSC; cm = rng0.normal(size=(5, 4))*0.6
NOFLEX = os.environ.get('NOFLEX', '0') == '1'
FM = fsup.copy() if NOFLEX else np.r_[fsup, fflex]; ZM = zsup.copy() if NOFLEX else np.r_[zsup, zflex]; NM = len(FM)
RES = [np.outer(Cr[:, j], Br[j])/(2*np.pi*fsup[j])**2 for j in range(3)] + \
      ([] if NOFLEX else [np.outer(cm[i], bm[i])/(2*np.pi*fflex[i])**2 for i in range(5)])
NLIN = 32*NM; IW = NLIN; IZ = NLIN+NM; ITR = NLIN+2*NM; NB = ITR + 8*NM

def ct_ss():
    n = 2*NM; A = np.zeros((n, n)); B = np.zeros((n, 4)); C = np.zeros((4, n))
    for j in range(NM):
        w = 2*np.pi*FM[j]; i = 2*j
        A[i, i+1] = 1; A[i+1, i] = -w**2; A[i+1, i+1] = -2*ZM[j]*w
        if j < 3: B[i+1] = Br[j]; C[:, i] = Cr[:, j]
        else: B[i+1] = bm[j-3]; C[:, i] = cm[j-3]
    return A, B, C

def controller(fc=None):
    fc = float(os.environ.get('FC', '25')) if fc is None else fc
    Ty = np.linalg.pinv(Cr); Tu = np.linalg.pinv(Br)
    wc = 2*np.pi*fc; wz, wp = wc/3, wc*3; kp = wc**2/3
    num = kp*np.array([1/wz, 1.0]); den = np.array([1/wp, 1.0])
    wl = 2*np.pi*FLP; n = np.polymul(num, [wl**2]); d = np.polymul(den, [1, 2*0.5*wl, wl**2])
    nd, dd, _ = ss.cont2discrete((n, d), Ts, method='bilinear'); nd = np.atleast_1d(np.squeeze(nd))
    Ad, Bd, Cd, Dd = ss.tf2ss(nd, dd)
    AK = sla.block_diag(*[Ad]*3); BK = sla.block_diag(*[Bd]*3); CK = sla.block_diag(*[Cd]*3); DK = np.kron(np.eye(3), Dd)
    return AK, BK @ Ty, Tu @ CK, Tu @ DK @ Ty

def closed_loop():
    A, B, C = ct_ss()
    Ad, Bd, Cd, _, _ = ss.cont2discrete((A, B, C, np.zeros((4, 4))), Ts, method='zoh')
    AK, BK, CK, DK = controller(); nk = AK.shape[0]
    Acl = np.block([[Ad - Bd @ DK @ Cd, -Bd @ CK], [BK @ Cd, AK]])
    Bcl = np.block([[Bd, -Bd @ DK, Bd], [np.zeros((nk, 4)), BK, np.zeros((nk, 4))]])   # [d; v; w(input disturbance)]
    Cy = np.block([[Cd, np.zeros((4, nk))]]); Dy = np.hstack([np.zeros((4, 4)), np.eye(4), np.zeros((4, 4))])
    Cu = np.block([[-DK @ Cd, -CK]]); Du = np.hstack([np.eye(4), -DK, np.zeros((4, 4))])
    return Acl, Bcl, np.vstack([Cy, Cu]), np.vstack([Dy, Du]), np.max(np.abs(np.linalg.eigvals(Acl)))

_CL = None
SW = float(os.environ.get('SW', '0')); AL = float(os.environ.get('AL', '1')); AH = float(os.environ.get('AH', '1')); FX = float(os.environ.get('FX', '40')); NO = int(os.environ.get('NO', '8'))
_bl = ss.butter(NO, FX, 'low', fs=fs, output='sos'); _bh = ss.butter(NO, FX, 'high', fs=fs, output='sos')
def dshape(zk):
    """|d spectrum| amplitude per unit sd: AL*H_low + AH*H_high (independent components) -> returns power factor"""
    w = np.angle(zk)
    hl = ss.sosfreqz(_bl, worN=w)[1]; hh = ss.sosfreqz(_bh, worN=w)[1]
    return AL**2*np.abs(hl)**2 + AH**2*np.abs(hh)**2
def simulate(N, sd, sv, rng, Ncut=6000):
    global _CL
    if _CL is None: _CL = closed_loop()
    Acl, Bcl, Cc, Dc, _ = _CL
    e1 = rng.standard_normal((N+Ncut, 4)); e2 = rng.standard_normal((N+Ncut, 4))
    d = sd*(AL*ss.sosfilt(_bl, e1, axis=0) + AH*ss.sosfilt(_bh, e2, axis=0))
    v = sv*rng.standard_normal((N+Ncut, 4)); w = SW*rng.standard_normal((N+Ncut, 4))
    _, out, _ = ss.dlsim(ss.dlti(Acl, Bcl, Cc, Dc, dt=Ts), np.hstack([d, v, w]))
    return d[Ncut:], out[Ncut:, 4:8], out[Ncut:, 0:4]

def zoh_2nd(w, z_, zk):
    A = np.array([[0, 1], [-w**2, -2*z_*w]]); B = np.array([[0], [1.0]])
    E = sla.expm(np.block([[A, B], [np.zeros((1, 3))]])*Ts); Ad, Bd = E[:2, :2], E[:2, 2:]
    a1 = -np.trace(Ad); a2 = np.linalg.det(Ad); det = zk**2 + a1*zk + a2
    x0 = ((zk - Ad[1, 1])*Bd[0, 0] + Ad[0, 1]*Bd[1, 0])/det
    x1 = (Ad[1, 0]*Bd[0, 0] + (zk - Ad[0, 0])*Bd[1, 0])/det
    return w**2*x0, w*x1, np.array([1, a1, a2])

def unpack(beta):
    lin = beta[:NLIN].reshape(NM, 2, 4, 4)
    return lin, beta[IW:IW+NM], beta[IZ:IZ+NM], beta[ITR:].reshape(NM, 4, 2)

def blocks(zk, lw, lz):
    h0 = []; h1 = []; ad = []
    for j in range(NM):
        a, b, c = zoh_2nd(np.exp(lw[j]), np.exp(lz[j]), zk); h0.append(a); h1.append(b); ad.append(c)
    return h0, h1, ad

def resp(zk, beta):
    lin, lw, lz, _ = unpack(beta); h0, h1, _ = blocks(zk, lw, lz)
    return sum(h0[j][:, None, None]*lin[j, 0][None] + h1[j][:, None, None]*lin[j, 1][None] for j in range(NM))

def out_jac(zk, beta, U, need_jac=True, dlog=1e-6, nonlin=True):
    lin, lw, lz, tr = unpack(beta); K = len(zk); zi = 1/zk
    h0, h1, ad = blocks(zk, lw, lz)
    adz = [ad[j][0] + ad[j][1]*zi + ad[j][2]*zi**2 for j in range(NM)]
    UL0 = [U @ lin[j, 0].T for j in range(NM)]; UL1 = [U @ lin[j, 1].T for j in range(NM)]
    Yh = np.zeros((K, 4), complex)
    for j in range(NM):
        Yh += h0[j][:, None]*UL0[j] + h1[j][:, None]*UL1[j]
        Yh += (tr[j][:, 0][None] + tr[j][:, 1][None]*zi[:, None])/adz[j][:, None]
    if not need_jac: return Yh, None
    J = np.zeros((K, 4, NB), complex)
    for j in range(NM):
        for q, h in ((0, h0[j]), (1, h1[j])):
            base = 32*j + 16*q
            for o in range(4): J[:, o, base+4*o:base+4*o+4] = h[:, None]*U
    if nonlin:
        for j in range(NM):
            for which, idx in ((0, IW+j), (1, IZ+j)):
                w2, z2 = lw[j] + (dlog if which == 0 else 0), lz[j] + (dlog if which == 1 else 0)
                a0, a1_, a2_ = zoh_2nd(np.exp(w2), np.exp(z2), zk)
                adz2 = a2_[0] + a2_[1]*zi + a2_[2]*zi**2
                tnum = tr[j][:, 0][None] + tr[j][:, 1][None]*zi[:, None]
                J[:, :, idx] = ((a0-h0[j])[:, None]*UL0[j] + (a1_-h1[j])[:, None]*UL1[j]
                                + tnum/adz2[:, None] - tnum/adz[j][:, None])/dlog
    for j in range(NM):
        for o in range(4):
            J[:, o, ITR+8*j+2*o] = 1/adz[j]; J[:, o, ITR+8*j+2*o+1] = zi/adz[j]
    return Yh, J

def iv_solve(zk, U, Y, Uh, W, beta0, it=40, tol=1e-9, chunk=1200, free=None):
    beta = beta0.copy(); K = len(zk)
    fr = np.arange(NB) if free is None else free
    for itn in range(it):
        Uhb = Uh(beta) if callable(Uh) else Uh
        A = np.zeros((len(fr), len(fr))); g = np.zeros(len(fr))
        for s0 in range(0, K, chunk):
            sl = slice(s0, min(K, s0+chunk))
            Yh, J = out_jac(zk[sl], beta, U[sl], nonlin=free is None)
            Z = J if Uhb is U else out_jac(zk[sl], beta, Uhb[sl], nonlin=free is None)[1]
            J = J[:, :, fr]; Z = Z[:, :, fr]
            ZW = (Z.conj()*W[sl][:, :, None]).reshape(-1, len(fr))
            A += np.real(ZW.T @ J.reshape(-1, len(fr))); g += np.real(ZW.T @ (Y[sl]-Yh).reshape(-1))
        sc = 1/np.sqrt(np.abs(np.diag(A)) + 1e-300)
        step = sc*np.linalg.solve(A*sc[:, None]*sc[None, :], g*sc)
        if free is not None:
            beta[fr] += step; return beta, 0
        mx = np.max(np.abs(step[IW:ITR]))
        if mx > 0.05: step *= 0.05/mx
        beta = beta + step
        if np.max(np.abs(step[IW:ITR])) < tol: break
    return beta, itn

LINIDX = np.r_[0:NLIN, ITR:NB]
def init_lin(zk, U, Y, W, beta):
    return iv_solve(zk, U, Y, U, W, beta, free=LINIDX)[0]

def truth_beta():
    b = np.zeros(NB)
    for j in range(NM): b[32*j:32*j+16] = RES[j].ravel()
    b[IW:IW+NM] = np.log(2*np.pi*FM); b[IZ:IZ+NM] = np.log(ZM)
    return b

def Ktf(zk):
    AK, BK, CK, DK = controller(); n = AK.shape[0]
    return np.array([CK @ np.linalg.solve(z*np.eye(n) - AK, BK) + DK for z in zk])

def fir_tud(d, u, L=600):
    N = d.shape[0]; Phi = np.zeros((N-L, 4*L))
    for l in range(L): Phi[:, 4*l:4*l+4] = d[L-l:N-l]
    Th = np.linalg.lstsq(Phi, u[L:], rcond=None)[0].reshape(L, 4, 4)
    return lambda zk: np.einsum('kl,lmo->kom', zk[:, None]**(-np.arange(L))[None], Th)

def simulate_periodic(P0, nper, sd, sv, rng, settle=2):
    """random-phase multisine d with per-line power P0*sd^2*dshape on every bin of one period; nper recorded periods"""
    global _CL
    if _CL is None: _CL = closed_loop()
    Acl, Bcl, Cc, Dc, _ = _CL
    kk = np.arange(1, P0//2); zz = np.exp(2j*np.pi*kk/P0)
    amp = sd*np.sqrt(P0*dshape(zz))
    Dl = np.zeros((P0, 4), complex)
    Dl[kk] = amp[:, None]*np.exp(2j*np.pi*rng.random((len(kk), 4)))
    Dl[P0-kk] = np.conj(Dl[kk])
    d1 = np.real(np.fft.ifft(Dl, axis=0))
    ntot = (nper+settle)*P0
    d = np.tile(d1, (nper+settle, 1)); v = sv*rng.standard_normal((ntot, 4)); w = SW*rng.standard_normal((ntot, 4))
    _, out, _ = ss.dlsim(ss.dlti(Acl, Bcl, Cc, Dc, dt=Ts), np.hstack([d, v, w]))
    c = settle*P0
    return d[c:], out[c:, 4:8], out[c:, 0:4]
