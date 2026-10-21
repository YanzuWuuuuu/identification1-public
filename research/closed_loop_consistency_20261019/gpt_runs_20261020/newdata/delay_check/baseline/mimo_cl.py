"""
Closed-loop frequency-domain additive SRIVC/IV consistency study (4x4 rigid + flexible).
Plant: CT, 3 rigid DOF (rank-3 double integrator) + 5 flexible modes, ZOH-sampled at fs.
Controller: DT decoupled PD-lead in logical coordinates (u = d - K y).
Model: CT additive, evaluated through its exact ZOH equivalent (no intersample error),
plus additive transient terms. Root of the IV equation found by IV-Gauss-Newton
(same roots as SRIVC fixed points).
"""
import numpy as np, scipy.linalg as sla, scipy.signal as ss, sys, time, json, os

fs = 5000.0; Ts = 1/fs
# ---------------- truth ----------------
xa = np.array([0.2, -0.2, -0.2, 0.2]); ya = np.array([0.15, 0.15, -0.15, -0.15])
m, Jx, Jy = 5.0, 0.05, 0.08
Br = np.vstack([np.ones(4)/m, ya/Jx, -xa/Jy])          # 3x4  (z, rx, ry) accelerations per unit force
xs = np.array([0.18, -0.18, -0.18, 0.18]); ys = np.array([0.13, 0.13, -0.13, -0.13])
Cr = np.vstack([np.ones(4), ys, -xs]).T                  # 4x3 sensor = z + y*rx - x*ry
M0 = Cr @ Br                                             # rank 3
fm = np.array([286.573, 344.945, 490.39, 718.958, 755.501])
zm = np.array([0.0119, 0.0132, 0.0186, 0.0322, 0.0346])
rng0 = np.random.default_rng(7)
bm = rng0.normal(size=(5, 4)) * 0.6 / np.sqrt(m)
cm = rng0.normal(size=(5, 4)) * 0.6
# flexible block i: cm_i bm_i^T * w^2/(s^2+2 z w s + w^2) / w^2  -> static gain cm bm^T / w^2
C0true = [np.outer(cm[i], bm[i]) for i in range(5)]      # numerator multiplying w^2/(den)  => gain/w^2 ... see h0
nb = 6  # blocks: rigid + 5

def ct_ss():
    # rigid: states q(3), qd(3)
    A = np.zeros((6 + 10, 6 + 10)); B = np.zeros((16, 4)); C = np.zeros((4, 16))
    A[0:3, 3:6] = np.eye(3); B[3:6, :] = Br; C[:, 0:3] = Cr
    for i in range(5):
        w = 2*np.pi*fm[i]; j = 6 + 2*i
        A[j, j+1] = 1; A[j+1, j] = -w**2; A[j+1, j+1] = -2*zm[i]*w
        B[j+1, :] = bm[i]; C[:, j] = cm[i]
    return A, B, C

# ---------------- controller ----------------
def controller():
    Ty = np.linalg.pinv(Cr)          # 3x4 sensor -> logical
    Tu = np.linalg.pinv(Br)          # 4x3 logical force -> actuators (in units giving unit accel)
    wc = 2*np.pi*25.0
    # per axis: plant ~ 1/s^2 (after Tu, Ty) ; lead: kp*(1+s/wz)/(1+s/wp)
    wz, wp = wc/3, wc*3; kp = wc**2/3
    num = kp*np.array([1/wz, 1.0]); den = np.array([1/wp, 1.0])
    # extra 2nd-order lowpass at 400 Hz
    wl = 2*np.pi*400; numl = [wl**2]; denl = [1, 2*0.7*wl, wl**2]
    n = np.polymul(num, numl); d = np.polymul(den, denl)
    nd, dd, _ = ss.cont2discrete((n, d), Ts, method='bilinear')
    nd = np.atleast_1d(np.squeeze(nd))
    Ad, Bd, Cd, Dd = ss.tf2ss(nd, dd)
    # 3 axes block diag
    k = Ad.shape[0]
    AK = sla.block_diag(*[Ad]*3); BK = sla.block_diag(*[Bd]*3); CK = sla.block_diag(*[Cd]*3); DK = np.kron(np.eye(3), Dd)
    # K: y(4) -> u(4) : u = Tu * C(z) * Ty y
    return AK, BK @ Ty, Tu @ CK, Tu @ DK @ Ty, (nd, dd), Ty, Tu

def closed_loop():
    A, B, C = ct_ss()
    Ad, Bd, Cd, _, _ = ss.cont2discrete((A, B, C, np.zeros((4, 4))), Ts, method='zoh')
    AK, BK, CK, DK, _, _, _ = controller()
    n, nk = Ad.shape[0], AK.shape[0]
    # u = d - (CK xk + DK y),  y = Cd x + v ; xk+ = AK xk + BK y
    # x+ = Ad x + Bd u
    Acl = np.block([[Ad - Bd @ DK @ Cd, -Bd @ CK], [BK @ Cd, AK]])
    Bd_ = np.block([[Bd, -Bd @ DK, Bd @ DK], [np.zeros((nk, 4)), BK, -BK]])   # inputs [d; v; r]
    Cy = np.block([[Cd, np.zeros((4, nk))]]); Dy = np.hstack([np.zeros((4, 4)), np.eye(4), np.zeros((4, 4))])
    Cu = np.block([[-DK @ Cd, -CK]]); Du = np.hstack([np.eye(4), -DK, DK])
    rho = np.max(np.abs(np.linalg.eigvals(Acl)))
    return Acl, Bd_, Cy, Dy, Cu, Du, rho

def simulate(N, sd, sv, rng, Ncut=4000, sr=0.0):
    Acl, Bcl, Cy, Dy, Cu, Du, _ = closed_loop()
    d = sd*rng.standard_normal((N+Ncut, 4)); v = sv*rng.standard_normal((N+Ncut, 4))
    r = sr*rng.standard_normal((N+Ncut, 4))
    w = np.hstack([d, v, r])
    sysd = ss.dlti(Acl, Bcl, np.vstack([Cy, Cu]), np.vstack([Dy, Du]), dt=Ts)
    _, out, _ = ss.dlsim(sysd, w)
    out = out[Ncut:]; d = d[Ncut:]; r = r[Ncut:]
    if sr > 0: return r, out[:, 4:8], out[:, 0:4]   # reference-injected: return r as the external signal
    return d, out[:, 4:8], out[:, 0:4]   # d, u, y

# ---------------- ZOH-equivalent model blocks ----------------
def zoh_2nd(w, z_, zk):
    """ZOH equivalents of h0 = w^2/(s^2+2 z w s+w^2), h1 = w s/(s^2+...). returns h0(zk), h1(zk), a_d coeffs"""
    A = np.array([[0, 1], [-w**2, -2*z_*w]]); B = np.array([[0], [1.0]])
    E = sla.expm(np.block([[A, B], [np.zeros((1, 3))]])*Ts)
    Ad, Bd = E[:2, :2], E[:2, 2:]
    # (zI-Ad)^-1 Bd
    a1 = -np.trace(Ad); a2 = np.linalg.det(Ad)
    det = zk**2 + a1*zk + a2
    adj00 = zk - Ad[1, 1]; adj01 = Ad[0, 1]; adj10 = Ad[1, 0]; adj11 = zk - Ad[0, 0]
    x0 = (adj00*Bd[0, 0] + adj01*Bd[1, 0])/det
    x1 = (adj10*Bd[0, 0] + adj11*Bd[1, 0])/det
    return w**2*x0, w*x1, np.array([1, a1, a2])

def zoh_rigid(zk):
    return Ts**2*(zk+1)/(2*(zk-1)**2), np.array([1, -2, 1.0])

class Model:
    """beta = [vec(M) (16), for each mode: vec(C0)(16), vec(C1)(16)] linear ; nonlinear: log w, log zeta per mode;
       transient: per block, per output: 2 real coeffs."""
    def __init__(s, zk):
        s.zk = zk; s.zi = 1/zk
        s.hr, s.ar = zoh_rigid(zk)
    def blocks(s, lw, lz):
        H = [s.hr]; Ad = [s.ar]
        for i in range(5):
            h0, h1, ad = zoh_2nd(np.exp(lw[i]), np.exp(lz[i]), s.zk)
            H += [h0, h1]; Ad.append(ad)
        return H, Ad
    def resp(s, lin, lw, lz):
        H, _ = s.blocks(lw, lz)
        G = H[0][:, None, None]*lin[0][None]
        for i in range(5):
            G = G + H[1+2*i][:, None, None]*lin[1+2*i][None] + H[2+2*i][:, None, None]*lin[2+2*i][None]
        return G

def unpack(beta):
    lin = [beta[16*j:16*(j+1)].reshape(4, 4) for j in range(11)]
    lw = beta[176:181]; lz = beta[181:186]; tr = beta[186:186+48].reshape(6, 4, 2)
    return lin, lw, lz, tr
NB = 186 + 48

def model_out_and_jac(mod, beta, U, need_jac=True, dlog=1e-6):
    """Yhat_k (K x 4) incl. transient, and Jacobian dYhat/dbeta (K x 4 x NB)"""
    lin, lw, lz, tr = unpack(beta)
    H, Ad = mod.blocks(lw, lz)
    K = U.shape[0]; zi = mod.zi
    Hs = [H[0]] + H[1:]   # 11 scalar responses aligned with lin
    Yh = np.zeros((K, 4), complex)
    for j in range(11): Yh += Hs[j][:, None]*(U @ lin[j].T)
    # transient: block b denominator ad (in z^-1): ad0 + ad1 z^-1 + ad2 z^-2
    adz = [Ad[b][0] + Ad[b][1]*zi + Ad[b][2]*zi**2 for b in range(6)]
    for b in range(6):
        Yh += (tr[b][:, 0][None, :] + tr[b][:, 1][None, :]*zi[:, None])/adz[b][:, None]
    if not need_jac: return Yh, None
    J = np.zeros((K, 4, NB), complex)
    for j in range(11):
        for o in range(4):
            # d Yh[:,o] / d lin[j][o,m] = H_j * U[:,m]
            J[:, o, 16*j+4*o:16*j+4*o+4] = Hs[j][:, None]*U
    # nonlinear params by finite differences (only the flexible blocks change)
    for i in range(5):
        for which, idx in ((0, 176+i), (1, 181+i)):
            lw2, lz2 = lw.copy(), lz.copy()
            if which == 0: lw2[i] += dlog
            else: lz2[i] += dlog
            h0, h1, ad = zoh_2nd(np.exp(lw2[i]), np.exp(lz2[i]), mod.zk)
            dh0 = (h0 - H[1+2*i])/dlog; dh1 = (h1 - H[2+2*i])/dlog
            adz2 = ad[0] + ad[1]*zi + ad[2]*zi**2
            tnum = tr[1+i][:, 0][None, :] + tr[1+i][:, 1][None, :]*zi[:, None]
            dtr = (tnum/adz2[:, None] - tnum/adz[1+i][:, None])/dlog
            J[:, :, idx] = dh0[:, None]*(U @ lin[1+2*i].T) + dh1[:, None]*(U @ lin[2+2*i].T) + dtr
    for b in range(6):
        for o in range(4):
            J[:, o, 186+8*b+2*o] = 1/adz[b]
            J[:, o, 186+8*b+2*o+1] = zi/adz[b]
    return Yh, J

def iv_solve(mod, U, Y, Uh, W, beta0, it=25, tol=1e-9, chunk=1500, verbose=False):
    """Solve sum_k Re{Z_k^H W_k (Y_k - Yhat_k(beta))} = 0, Z = dYhat/dbeta evaluated with Uh (instrument input).
       IV-Gauss-Newton: beta += (Re sum Z^H W J)^-1 Re sum Z^H W E. Uh may be callable(beta)->Uh."""
    beta = beta0.copy(); K = U.shape[0]
    for itn in range(it):
        Uhb = Uh(beta) if callable(Uh) else Uh
        A = np.zeros((NB, NB)); g = np.zeros(NB)
        for s0 in range(0, K, chunk):
            sl = slice(s0, min(K, s0+chunk))
            m2 = Model.__new__(Model); m2.zk = mod.zk[sl]; m2.zi = mod.zi[sl]; m2.hr = mod.hr[sl]; m2.ar = mod.ar
            Yh, J = model_out_and_jac(m2, beta, U[sl])
            if Uhb is U: Z = J
            else: _, Z = model_out_and_jac(m2, beta, Uhb[sl])
            E = Y[sl] - Yh
            Wk = W[sl]                                  # K x 4 diagonal weights
            ZW = Z.conj()*Wk[:, :, None]
            ZW2 = ZW.reshape(-1, NB); A += np.real(ZW2.T @ J.reshape(-1, NB))
            g += np.real(ZW2.T @ E.reshape(-1))
        sc = 1/np.sqrt(np.abs(np.diag(A)) + 1e-300)
        step = sc*np.linalg.solve(A*sc[:, None]*sc[None, :], g*sc)
        mx = np.max(np.abs(step[176:186]))
        if mx > 0.05: step *= 0.05/mx
        beta = beta + step
        rel = np.max(np.abs(step[176:186]))
        if verbose: print(itn, rel, flush=True)
        if rel < tol: break
    return beta, itn

def truth_beta():
    beta = np.zeros(NB)
    beta[0:16] = M0.ravel()
    for i in range(5):
        beta[16*(1+2*i):16*(2+2*i)] = (C0true[i]/(2*np.pi*fm[i])**2).ravel()
    beta[176:181] = np.log(2*np.pi*fm); beta[181:186] = np.log(zm)
    return beta

def init_beta(mod, U, Y, W, rng):
    b = truth_beta()
    b[176:181] += np.log(1 + 0.003*rng.standard_normal(5))
    b[181:186] += np.log(1 + 0.15*rng.standard_normal(5))
    b[:176] = 0; b[186:] = 0
    # linear LS for linear params given denominators (one IV step with Z=J is LS)
    beta, _ = iv_solve_linear(mod, U, Y, W, b)
    return beta

def iv_solve_linear(mod, U, Y, W, beta, chunk=1500):
    K = U.shape[0]; idx = np.r_[0:176, 186:NB]
    A = np.zeros((len(idx), len(idx))); g = np.zeros(len(idx))
    for s0 in range(0, K, chunk):
        sl = slice(s0, min(K, s0+chunk))
        m2 = Model.__new__(Model); m2.zk = mod.zk[sl]; m2.zi = mod.zi[sl]; m2.hr = mod.hr[sl]; m2.ar = mod.ar
        Yh, J = model_out_and_jac(m2, beta, U[sl]); J = J[:, :, idx]
        ZW = J.conj()*W[sl][:, :, None]
        ZW2 = ZW.reshape(-1, len(idx)); A += np.real(ZW2.T @ J.reshape(-1, len(idx))); g += np.real(ZW2.T @ (Y[sl]-Yh).reshape(-1))
    sc = 1/np.sqrt(np.abs(np.diag(A)) + 1e-300)
    beta = beta.copy(); beta[idx] += sc*np.linalg.solve(A*sc[:, None]*sc[None, :], g*sc)
    return beta, 0

# ---------------- helpers for instruments ----------------
def Ktf(zk):
    AK, BK, CK, DK, _, _, _ = controller()
    n = AK.shape[0]; out = np.zeros((len(zk), 4, 4), complex)
    for k, z in enumerate(zk):
        out[k] = CK @ np.linalg.solve(z*np.eye(n) - AK, BK) + DK
    return out

def fir_tud(d, u, L=400):
    """two-stage: FIR LS estimate of d->u (4x4), returns function of zk"""
    N = d.shape[0]
    Phi = np.zeros((N-L, 4*L))
    for l in range(L): Phi[:, 4*l:4*l+4] = d[L-l:N-l]
    Th, *_ = np.linalg.lstsq(Phi, u[L:], rcond=None)   # (4L x 4)
    Th = Th.reshape(L, 4, 4)                            # Th[l][m_in, o_out]
    def T(zk):
        zz = zk[:, None]**(-np.arange(L))[None, :]
        return np.einsum('kl,lmo->kom', zz, Th)         # K x out x in
    return T
