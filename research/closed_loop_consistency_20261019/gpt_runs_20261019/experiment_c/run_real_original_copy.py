import numpy as np, sys, os, json
from load import load
from realfit import *
S = load(); D, U, Y, V, Qa, P = S['D'], S['U'], S['Y'], S['V'], S['Q'], S['P']
tau = np.loadtxt('/mnt/user-data/uploads/identifiation1/experiment_data/srivc_no_lrm_20261015/delay_identification_20261016/effective_delay_matrix_us.csv', delimiter=',')*1e-6
FLO = float(os.environ.get('FLO', '3')); FHI = float(os.environ.get('FHI', '25'))
F0 = [float(x) for x in os.environ.get('F0', '10.6,11.7,16.1').split(',')]

# ---- controller: exact ARX v <- y (deterministic), from a long segment
def fit_controller(y, v, na=6, nb=6):
    L = max(na, nb)
    Phi = np.hstack([v[L-i:len(v)-i] for i in range(1, na+1)] + [y[L-j:len(y)-j] for j in range(0, nb+1)])
    th, *_ = np.linalg.lstsq(Phi, v[L:], rcond=None)
    r = v[L:] - Phi @ th
    A = [th[3*(i):3*(i+1)].T for i in range(na)]            # each 3x3 : v(t) = sum A_i v(t-i) + sum B_j y(t-j)
    B = [th[3*na+4*j:3*na+4*(j+1)].T for j in range(nb+1)]  # each 3x4
    def Kv(zi):
        out = np.zeros((len(zi), 3, 4), complex)
        for k, x in enumerate(zi):
            Az = np.eye(3) - sum(A[i]*x**(i+1) for i in range(na)); Bz = sum(B[j]*x**j for j in range(nb+1))
            out[k] = np.linalg.solve(Az, Bz)
        return out
    return Kv, np.abs(r).max()/np.abs(v).max()

def estimate(p0, p1, methods, verbose=False, return_beta=False):
    seg = slice(p0*P, p1*P); N = (p1-p0)*P; nper = p1-p0
    u, y, d, v = U[seg], Y[seg], D[seg], V[seg]
    k = np.arange(1, N//2); f = k*fs/N
    Uk = np.fft.fft(u, axis=0)[k]; Yk = np.fft.fft(y, axis=0)[k]; Dk = np.fft.fft(d, axis=0)[k]; Vk = np.fft.fft(v, axis=0)[k]
    band = (f >= FLO) & (f <= FHI); lines = band & (k % nper == 0)
    Kv, kres = fit_controller(y, v)
    out = {'kres': kres}
    res = {}; betas = {}
    for m in methods:
        sel = band if m == 'dir_all' else lines
        c = Ctx(f[sel], tau); Us, Ys, Ds = Uk[sel], Yk[sel], Dk[sel]
        W = np.ones((sel.sum(), 4))
        b0 = np.zeros(NB); b0[IW:IZ] = np.log(2*np.pi*np.array(F0)); b0[IZ:ITR] = np.log(0.03)
        b0 = init_lin(c, Us, Ys, W, b0)
        if m in ('dir_all', 'dir_lines'): Uh = Us
        elif m == 'ref_K':
            Kz = -Qa[None] @ Kv(c.zi)                      # u = d + Q v = d - K y  ->  K = -Q Kv
            Uh = lambda b: np.einsum('kij,kj->ki', np.linalg.inv(np.eye(4)[None] + Kz @ resp(c, b)), Ds)
            out['K_check'] = float(np.linalg.norm(Vk[sel] - np.einsum('kij,kj->ki', Kv(c.zi), Ys))/np.linalg.norm(Vk[sel]))
        elif m == 'ref_loo':
            idx = np.where(sel)[0]; Uh = np.zeros_like(Us); HW = 6
            for ii in range(len(idx)):
                js = [j for j in range(max(0, ii-HW), min(len(idx), ii+HW+1)) if j != ii]
                T = np.linalg.lstsq(Ds[js], Us[js], rcond=None)[0]
                Uh[ii] = Ds[ii] @ T
        b, it = iv_solve(c, Us, Ys, Uh, W, b0, verbose=verbose)
        fz = modal(b); betas[m] = b.tolist()
        res[m] = dict(f=fz[0].tolist(), zeta=fz[1].tolist(), it=int(it))
        # residual on lines (one-step fit quality)
        Yh, _ = out_jac(c, b, Us, need_jac=False)
        res[m]['rel_resid'] = float(np.linalg.norm(Ys-Yh)/np.linalg.norm(Ys))
        print(p0, p1, m, 'f', np.round(fz[0], 3), 'zeta%', np.round(fz[1], 2), 'it', it, 'resid', round(res[m]['rel_resid'], 4), flush=True)
    out['res'] = res
    if return_beta: out['beta'] = betas
    return out

if __name__ == '__main__':
    p0, p1 = int(sys.argv[1]), int(sys.argv[2]); methods = sys.argv[3].split(',')
    o = estimate(p0, p1, methods, verbose=os.environ.get('VB') == '1')
    print(json.dumps({'kres': o['kres'], 'K_check': o.get('K_check')}))
    json.dump(o, open(f'res_z{os.environ.get("ZPOS","0")}_{p0}_{p1}_{FLO}_{FHI}.json', 'w'))
