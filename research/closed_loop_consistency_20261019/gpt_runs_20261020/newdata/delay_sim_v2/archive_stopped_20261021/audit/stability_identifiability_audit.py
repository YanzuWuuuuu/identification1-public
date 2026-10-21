"""Read-only source audit; writes only this new audit directory.

No fits are performed.  The exact old delay-register state update is extracted
from the source AST without importing run_sim (whose module overrides NOFLEX).
"""
import os
import sys
for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "1"
os.environ.update(NOFLEX="0", FC="60", AL="30", AH="100", SW=".03",
                  FLP="120", FSC="3", NO="8", FX="40")
sys.dont_write_bytecode = True
from pathlib import Path
import ast
import json
import hashlib
from datetime import datetime, timezone
import numpy as np
import scipy.linalg as la
from scipy import signal, optimize

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "gpt_runs_20261020/newdata/delay_check/code"))
import cl2
import core

FS = 5000.
NI = np.array([0, 1, 0, 1])
NO = np.array([1, 0, 1, 0])
TAU = (NO[:, None] + NI[None, :]) / FS
BANDS = {"low_3_100": (3., 100.), "high_100_550": (100., 550.),
         "wide_3_550": (3., 550.)}


def encode(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(type(obj).__name__)


def delayed_state_space():
    path = ROOT / "gpt_runs_20261020/newdata/delay_check/code/run_sim.py"
    source = path.read_text(encoding="utf8")
    node = next(n for n in ast.parse(source).body
                if isinstance(n, ast.FunctionDef) and n.name == "plant_with_delays")
    scope = dict(cl2=cl2, signal=signal, np=np, FS=FS, NI=NI, NO=NO)
    exec(compile(ast.Module(body=[node], type_ignores=[]),
                 str(path) + ":plant_with_delays", "exec"), scope)
    return scope["plant_with_delays"]()


def loop(f):
    f = np.atleast_1d(f)
    z = np.exp(2j * np.pi * f / FS)
    G = cl2.resp(z, cl2.truth_beta()) * np.exp(
        -2j * np.pi * f[:, None, None] * TAU)
    K = cl2.Ktf(z)
    return G @ K, G, K


def metrics(L):
    return dict(sigma_max=np.linalg.svd(L, compute_uv=False)[:, 0],
                spectral_radius=np.abs(np.linalg.eigvals(L)).max(axis=1),
                max_abs_element=np.abs(L).max(axis=(1, 2)))


def loop_audit():
    A, _, _, _ = delayed_state_space()
    eig = np.linalg.eigvals(A)
    f = np.unique(np.r_[np.arange(3., 550.0001, .05), cl2.FM,
                        [3., 100., 550.]])
    f = f[(f >= 3) & (f <= 550)]
    L, _, _ = loop(f)
    vals = metrics(L)
    result = dict(closed_loop_radius=float(abs(eig).max()),
                  stable=bool(np.all(abs(eig) < 1)),
                  no_delay_closed_loop_radius=float(cl2.closed_loop()[-1]),
                  unstable_poles=[
                      dict(real=x.real, imag=x.imag, modulus=abs(x),
                           frequency_hz=abs(np.angle(x)) * FS / (2 * np.pi))
                      for x in eig if abs(x) >= 1], bands={})
    for band, (lo, hi) in BANDS.items():
        ix = np.flatnonzero((f >= lo) & (f <= hi))
        out = {}
        for name, val in vals.items():
            k = ix[np.argmax(val[ix])]
            bounds = (max(lo, f[k] - .2), min(hi, f[k] + .2))
            def objective(ff):
                ll, _, _ = loop(ff)
                return -float(metrics(ll)[name][0])
            opt = optimize.minimize_scalar(objective, bounds=bounds,
                method="bounded", options={"xatol": 1e-9})
            choices = [(float(val[k]), float(f[k])),
                       (-float(opt.fun), float(opt.x))]
            for endpoint in (lo, hi):
                choices.append((-objective(endpoint), endpoint))
            value, freq = max(choices)
            out[name] = dict(value=value, frequency_hz=freq,
                             below_point_one=bool(value < .1))
        result["bands"][band] = out
    fp = np.r_[100., cl2.FM[3:6], 550.]
    LL, GG, KK = loop(fp)
    vv = metrics(LL)
    result["modal_and_edge_points"] = [
        dict(frequency_hz=ff, sigma_GK=vv["sigma_max"][i],
             rho_GK=vv["spectral_radius"][i],
             max_element_GK=vv["max_abs_element"][i],
             sigma_G=la.svdvals(GG[i])[0], sigma_K=la.svdvals(KK[i])[0])
        for i, ff in enumerate(fp)]
    np.savez_compressed(OUT / "loop_gain_grid.npz", frequency_hz=f, **vals)
    return result


def information_audit():
    m = core.DelayModel(8, TAU, 550.)
    b = np.zeros(m.nb)
    for o in range(4):
        for j in range(8):
            b[o*m.L+10*j:o*m.L+10*j+4] = cl2.RES[j][o]
    b[m.iw:m.iz] = np.log(2 * np.pi * cl2.FM)
    b[m.iz:m.base_nb] = np.log(cl2.ZM)
    f = np.arange(18, 3301) / 6.
    z = np.exp(2j * np.pi * f / FS)
    ctx = m.context(f)
    G = m.resp(ctx, b)
    K = cl2.Ktf(z)
    rng = np.random.default_rng(20261021)
    D = .2 * np.sqrt(30000 * cl2.dshape(z))[:, None] * np.exp(
        2j * np.pi * rng.random((len(f), 4)))
    U = np.linalg.solve(np.eye(4)[None] + K @ G, D[:, :, None])[:, :, 0]
    _, J = m.output(ctx, b, U)
    result = dict(parameters=m.nb, linear_parameters=m.iw,
        nuisance_parameters=m.base_nb, delay_parameters=7,
        true_C1_D0_D2_zero=True, random_phase_seed=20261021,
        grid_hz=1/6,
        diagnostic="Noiseless reference U; column-normalized real output "
                   "Jacobian, including all modal transients. This is not "
                   "a fit or a noisy-direct-IV bias calculation.",
        normalization="Xn=X diag(1/||X_j||); Mn=Xn^T Xn is the two-sided "
                      "diagonally normalized deterministic Gram.",
        bands={})
    for name, (lo, hi) in BANDS.items():
        ix = (f >= lo) & (f <= hi)
        X = core.realvec(J[ix]).reshape(-1, m.nb)
        Xn = X / np.linalg.norm(X, axis=0)
        V = Xn[:, :m.base_nb]
        Z = Xn[:, m.base_nb:]
        coef, _, rank, _ = la.lstsq(V, Z, cond=1e-12,
                                   lapack_driver="gelsy")
        R = Z - V @ coef
        sr = la.svdvals(R)
        sx = la.svdvals(Xn)
        sv = la.svdvals(V)
        Mn = Xn.T @ Xn
        sm = la.svdvals(Mn)
        em = la.eigvalsh(Mn)
        thresholds = (1e-6, 1e-8, 1e-10, 1e-12)
        result["bands"][name] = dict(rows=len(X),
            other_rank_pivoted_qr_cond_1e12=rank,
            other_normalized_J_singular_values=sv,
            delay_unexplained_energy_fractions=np.sum(R**2, axis=0),
            delay_projected_singular_values=sr,
            full_normalized_J_singular_values=sx,
            full_normalized_J_condition=sx[0]/sx[-1],
            full_normalized_Gram_condition_from_J=(sx[0]/sx[-1])**2,
            full_normalized_Gram_singular_values=sm,
            full_normalized_Gram_condition_direct=sm[0]/sm[-1],
            full_normalized_Gram_eigenvalue_min=em[0],
            full_normalized_J_rank={
                str(t): int(np.sum(sx > t*sx[0])) for t in thresholds},
            full_two_sided_normalized_Gram_rank={
                str(t): int(np.sum(sm > t*sm[0])) for t in thresholds},
            caution="Forming a Gram squares conditioning and loses tiny "
                    "singular directions. Low-band numerical projection "
                    "residual near roundoff is not a proof of an exact "
                    "ZOH structural identity.")
    return result


def main():
    files = [ROOT / "code/cl2.py",
        ROOT / "gpt_runs_20261020/newdata/delay_check/code/run_sim.py",
        ROOT / "gpt_runs_20261020/newdata/delay_check/code/core.py",
        ROOT / "gpt_runs_20261020/newdata/delay_check/baseline/modal_model.py"]
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        source_hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in files},
        config=dict(NOFLEX=0, FC=60, FLP=120, FSC=3, AL=30, AH=100,
                    SW=.03, fs=FS, NI=NI, NO=NO, true_delay_seconds=TAU,
                    natural_frequencies_hz=cl2.FM),
        loop=loop_audit(), information=information_audit())
    (OUT / "stability_identifiability_audit.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, default=encode),
        encoding="utf8")
    print(json.dumps(dict(loop=result["loop"],
        information={name: {
            "delay_energy_fraction": rr["delay_unexplained_energy_fractions"],
            "J_condition": rr["full_normalized_J_condition"],
            "J_rank": rr["full_normalized_J_rank"],
            "two_sided_Gram_rank": rr["full_two_sided_normalized_Gram_rank"]}
            for name, rr in result["information"]["bands"].items()}),
        indent=2, default=encode))


if __name__ == "__main__":
    main()
