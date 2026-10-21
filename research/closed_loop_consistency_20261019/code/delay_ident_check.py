"""Numerical checks for the delay extension (Lemma 3, Prop 4, Prop 5).
Additive real-residue modal model G = E(tau) o [sum_i (C0_i w_i^2 + C1_i w_i s)/den_i + D0 + D2 (s/wr)^2],
tau_jm = tau_c + a_j + b_m (zero-sum a, b). For each delay direction we report the fraction of its
Jacobian column that is NOT explained by all other parameter columns:
   r = ||(I-P_other) j||^2 / ||j||^2   (r = 0: first-order unidentifiable, absorbed).
"""
import numpy as np, scipy.linalg as sla, json, sys
rng = np.random.default_rng(1)
fs = 5000.0; T = 1/fs

def helmert(n):
    H = np.zeros((n-1, n))
    for k in range(1, n):
        H[k-1, :k] = 1/k; H[k-1, k] = -1
        H[k-1] /= np.linalg.norm(H[k-1])
    return H  # rows: orthonormal zero-sum vectors

def zoh_basis(w, z, om):
    A = np.array([[0, 1], [-w**2, -2*z*w]]); B = np.array([0., 1.])
    Ad = sla.expm(A*T); Bd = np.linalg.solve(A, (Ad - np.eye(2)) @ B)
    zz = np.exp(1j*om*T); out0 = []; out1 = []
    for q in zz:
        x = np.linalg.solve(q*np.eye(2) - Ad, Bd)
        out0.append(w**2*x[0]); out1.append(w*x[1])
    return np.array(out0), np.array(out1)

def basis(w, z, om, zoh):
    s = 1j*om
    if zoh: return zoh_basis(w, z, om)
    den = s**2 + 2*z*w*s + w**2
    return w**2/den, w*s/den

def build(f, modes, C0, C1, D0, D2, tau, rank1=None, zoh=False, wr=2*np.pi*25):
    om = 2*np.pi*f; s = 1j*om; K = len(f); E = np.exp(-s[:, None, None]*tau[None])
    cols = []  # each column: K x 4 x 4 complex dG
    Gt = D0[None] + (s/wr)[:, None, None]**2*D2[None]
    for i, (w, z) in enumerate(modes):
        b0, b1 = basis(w, z, om, zoh)
        Gt = Gt + b0[:, None, None]*C0[i][None] + b1[:, None, None]*C1[i][None]
    G = E*Gt
    for i, (w, z) in enumerate(modes):
        b0, b1 = basis(w, z, om, zoh)
        if rank1 is None:
            for j in range(4):
                for m in range(4):
                    e = np.zeros((4, 4)); e[j, m] = 1
                    cols.append(E*b0[:, None, None]*e); cols.append(E*b1[:, None, None]*e)
        else:
            psi, phi, k0, k1 = rank1[i]
            R = np.outer(psi, phi)
            cols.append(E*b0[:, None, None]*R); cols.append(E*b1[:, None, None]*R)       # k0, k1
            for j in range(4):
                dR = np.outer(np.eye(4)[j], phi)
                cols.append(E*(k0*b0+k1*b1)[:, None, None]*dR)                          # psi_j
                dR = np.outer(psi, np.eye(4)[j])
                cols.append(E*(k0*b0+k1*b1)[:, None, None]*dR)                          # phi_j
        for which in (0, 1):                                                              # log w, log z
            h = 1e-6; ww, zz_ = (w*np.exp(h), z) if which == 0 else (w, z*np.exp(h))
            c0p, c1p = basis(ww, zz_, om, zoh)
            Ci0, Ci1 = (C0[i], C1[i]) if rank1 is None else (rank1[i][2]*np.outer(rank1[i][0], rank1[i][1]), rank1[i][3]*np.outer(rank1[i][0], rank1[i][1]))
            cols.append(E*(((c0p-b0)/h)[:, None, None]*Ci0 + ((c1p-b1)/h)[:, None, None]*Ci1))
    for j in range(4):
        for m in range(4):
            e = np.zeros((4, 4)); e[j, m] = 1
            cols.append(E*e[None]*np.ones(K)[:, None, None])                               # D0
            cols.append(E*((s/wr)**2)[:, None, None]*e)                                  # D2
    Hh = helmert(4)
    dcols = {'tau_c': -s[:, None, None]*G}
    for k in range(3):
        dcols[f'a{k+1}'] = -s[:, None, None]*(Hh[k][None, :, None]*G)
        dcols[f'b{k+1}'] = -s[:, None, None]*(G*Hh[k][None, None, :])
    st = lambda c: np.concatenate([c.reshape(K, -1).real.ravel(), c.reshape(K, -1).imag.ravel()])
    X = np.array([st(c) for c in cols]).T
    Dl = {k: st(v) for k, v in dcols.items()}
    return X, Dl

def resid_frac(X, Dl):
    out = {}
    names = list(Dl)
    for k in names:
        others = np.column_stack([X] + [Dl[o] for o in names if o != k])
        sc = np.linalg.norm(others, axis=0); sc[sc == 0] = 1
        coef, *_ = np.linalg.lstsq(others/sc, Dl[k], rcond=1e-12)
        r = Dl[k] - (others/sc) @ coef
        out[k] = float(np.dot(r, r)/np.dot(Dl[k], Dl[k]))
    return out

def identity_check():
    w, z = 2*np.pi*333, 0.007; C0, C1 = 0.7, -0.3
    f = np.linspace(200, 500, 2001); s = 2j*np.pi*f; den = s**2 + 2*z*w*s + w**2
    lhs = s*(C0*w**2 + C1*w*s)/den
    rhs = C1*w + w*(C0 - 2*z*C1)*(w*s/den) - C1*w*(w**2/den)
    return float(np.max(np.abs(lhs - rhs))/np.max(np.abs(lhs)))

def case(f, modes, bg, rank1, zoh, label):
    tau = 300e-6 + 20e-6*rng.standard_normal(4)[:, None] + 20e-6*rng.standard_normal(4)[None]
    if rank1:
        r1 = [(rng.standard_normal(4), rng.standard_normal(4), rng.standard_normal(), 0.3*rng.standard_normal()) for _ in modes]
        C0 = [r[2]*np.outer(r[0], r[1]) for r in r1]; C1 = [r[3]*np.outer(r[0], r[1]) for r in r1]
    else:
        r1 = None; C0 = [rng.standard_normal((4, 4)) for _ in modes]; C1 = [0.3*rng.standard_normal((4, 4)) for _ in modes]
    D0 = bg*rng.standard_normal((4, 4)); D2 = bg*rng.standard_normal((4, 4))
    X, Dl = build(f, modes, C0, C1, D0, D2, tau, r1, zoh)
    res = resid_frac(X, Dl)
    print(f"{label:58s}", {k: f"{v:.2e}" for k, v in res.items()}, flush=True)
    return res

if __name__ == '__main__':
    out = {'identity_rel_err': identity_check()}
    print('Lemma 3 identity relative error:', out['identity_rel_err'])
    fH = np.arange(220, 420, 0.5); mH = [(2*np.pi*282, 0.0195), (2*np.pi*333, 0.0065), (2*np.pi*359, 0.0025)]
    fL = np.arange(3, 25, 1/6); mL = [(2*np.pi*10.4, 0.049), (2*np.pi*11.5, 0.025), (2*np.pi*15.4, 0.035)]
    C = {}
    C['H_free_noBG_CT'] = case(fH, mH, 0.0, False, False, 'H300 free residues, D0=D2=0, CT')
    C['H_free_BG_CT'] = case(fH, mH, 0.1, False, False, 'H300 free residues, background 10%, CT')
    C['H_free_noBG_ZOH'] = case(fH, mH, 0.0, False, True, 'H300 free residues, D0=D2=0, ZOH')
    C['H_rank1_noBG_CT'] = case(fH, mH, 0.0, True, False, 'H300 rank-1 shared shapes, D0=D2=0, CT')
    C['H_rank1_BG_CT'] = case(fH, mH, 0.1, True, False, 'H300 rank-1 shared shapes, background 10%, CT')
    C['L_free_BG_CT'] = case(fL, mL, 0.1, False, False, 'LOW 3-25 Hz free residues, background 10%, CT')
    C['L_rank1_BG_CT'] = case(fL, mL, 0.1, True, False, 'LOW 3-25 Hz rank-1, background 10%, CT')
    out['cases'] = C
    json.dump(out, open('delay_ident_check_results.json', 'w'), indent=1)
