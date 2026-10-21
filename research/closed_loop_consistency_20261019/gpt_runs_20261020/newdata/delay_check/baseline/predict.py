"""First-order asymptotic bias of the direct FD-SRIVC root (Theorem: delta = M^-1 b) and
asymptotic covariance of direct/reference IV, for the 4x4 closed-loop example."""
import numpy as np, sys
from mimo_cl import *

def gl_tensor(zsel, beta):
    mod = Model(zsel); K = len(zsel)
    Gl = np.zeros((K, 4, 4, 186), complex)
    for mi in range(4):
        Ue = np.zeros((K, 4)); Ue[:, mi] = 1
        _, J = model_out_and_jac(mod, beta, Ue)
        Gl[:, :, mi, :] = J[:, :, :186]
    return Gl

def predict(N, sd, sv, sel, W, Gl, G0, Kz, method='dir', sr=0.0):
    """returns delta (186,), covariance (186x186) per unit (DFT unnormalized: Phi = N*sigma^2)"""
    K = sel.sum(); I4 = np.eye(4)
    G = G0[sel]; Kk = Kz[sel]
    S = np.linalg.inv(I4[None] + Kk @ G)               # d -> u
    F = S @ Kk                                          # v -> -u
    Phid = N*sd**2; Phiv = N*sv**2; Phir = N*sr**2
    Sx = S if sr == 0 else S @ Kk                    # external -> u
    Phx = Phid if sr == 0 else Phir
    Phiu = Phx*Sx @ Sx.conj().transpose(0, 2, 1) + Phiv*F @ F.conj().transpose(0, 2, 1)
    Gls = Gl                                            # K x 4 x 4 x P
    # b_l = - sum_k Re tr(F^H G_l^H W Phiv) = - Phiv * sum_k Re sum_i w_i conj((G_l F)_ii)
    GF = np.einsum('kijp,kjl->kilp', Gls, F)
    diag = np.einsum('kiip->kip', GF)
    b = -Phiv*np.real(np.einsum('ki,kip->p', W, diag.conj()))
    # Z = G_l Uhat. direct: Uhat=U ; ref: Uhat = S D
    Wh = np.sqrt(W)
    if method == 'dir':
        Lu = np.linalg.cholesky(Phiu + 1e-30*np.eye(4)[None])
        Q = np.einsum('ki,kijp,kjl->kilp', Wh, Gls, Lu).reshape(K, 16, -1)
        M = np.real(np.einsum('kap,kaq->pq', Q.conj(), Q))
        b_used = b
    else:
        # A = sum Re E[(G_l S D)^H W G_m U] = sum Re tr(G_l^H W G_m S Phid S^H)
        Ls = np.sqrt(Phx)*Sx
        Q = np.einsum('ki,kijp,kjl->kilp', Wh, Gls, Ls).reshape(K, 16, -1)
        M = np.real(np.einsum('kap,kaq->pq', Q.conj(), Q))
        b_used = np.zeros_like(b)
    sc = 1/np.sqrt(np.diag(M))
    Ms = M*sc[:, None]*sc[None, :]
    delta = sc*np.linalg.solve(Ms, b_used*sc)
    # covariance of IV root: M^-1 Cov(sum Re Z^H W V) M^-T ; Cov ~ 0.5*sum Re tr(Zcov..) -> use Q-type
    # Cov(g) where g_l = sum Re (Z_l^H W V): for complex circular V: Cov = 0.5 * sum Re[ Zt^H W Phiv W Zt ]
    if method == 'dir':
        Lz = np.linalg.cholesky(Phiu + 1e-30*np.eye(4)[None])
    else:
        Lz = np.sqrt(Phx)*Sx
    Qz = np.einsum('ki,kijp,kjl->kilp', W*np.sqrt(Phiv), Gls, Lz).reshape(K, 16, -1)
    Cg = 0.5*np.real(np.einsum('kap,kaq->pq', Qz.conj(), Qz))
    Mi = sc[:, None]*np.linalg.inv(Ms)*sc[None, :]
    P = Mi @ Cg @ Mi.T
    return delta, P

if __name__ == '__main__':
    N = 20000; sv = float(sys.argv[1]); sr = float(sys.argv[2]) if len(sys.argv)>2 else 0.0; sd = 0.0 if sr>0 else 1.0
    k = np.arange(1, N//2); f = k*fs/N; zk = np.exp(2j*np.pi*k/N)
    band = (f >= 10) & (f <= 800); full = np.ones_like(band)
    bt = truth_beta(); lin0, lw0, lz0, _ = unpack(bt)
    G0 = Model(zk).resp(lin0, lw0, lz0)
    wg = 1/np.sum(np.abs(G0)**2, axis=2); wg /= wg[band].mean(0, keepdims=True)
    Kz = Ktf(zk)
    for name, sel, W, meth in [('dir_full_W1', full, None, 'dir'), ('dir_band_W1', band, None, 'dir'),
                               ('dir_band_Wg', band, 'g', 'dir'), ('ref_band_Wg', band, 'g', 'ref')]:
        Ws = np.ones((sel.sum(), 4)) if W is None else wg[sel]
        Gl = gl_tensor(zk[sel], bt)
        dl, P = predict(N, sd, sv, sel, Ws, Gl, G0, Kz, meth, sr)
        sd_ = np.sqrt(np.diag(P))
        print(f"{name:12s} bias f%:", np.round(dl[176:181]*100, 4), " z%:", np.round(dl[181:186]*100, 2))
        print(f"{'':12s} std  f%:", np.round(sd_[176:181]*100, 4), " z%:", np.round(sd_[181:186]*100, 2),
              " rigid M bias/|M|:", np.round(np.linalg.norm(dl[:16])/np.linalg.norm(bt[:16])*100, 3), '%', flush=True)
