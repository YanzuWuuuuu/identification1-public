"""Four-mode + real constant D: deterministic Jacobian/information checks.

No identification, Monte Carlo, fifth-mode dynamics, or MATLAB-file mutation.
D_demo is an arbitrary algebra example, NOT an estimated stage parameter.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import numpy as np
import scipy
from scipy.linalg import qr, svdvals

ROOT = Path(__file__).resolve().parents[2]
Q = 16
TAU = 10e-6


def realify(x):
    return np.sqrt(2) * np.concatenate((x.real, x.imag), axis=0)


def geometry():
    source = ROOT / "matlab/3/model.m"
    text = source.read_text(encoding="utf-8")
    def matrix(name):
        value = re.search(r"\b" + name + r"\s*=\s*\[([^\]]+)\]", text).group(1)
        return np.array([[float(x) for x in row.replace("\n", " ").split(",")]
                         for row in value.split(";") if row.strip()])
    br, bf, cr, cf = [matrix(x) for x in ("B1", "B2", "C1", "C2")]
    mass = (cr @ np.diag([1/.0439, 1/.0384, 1/4.2809]) @ br).ravel(order="F")
    wn = 2*np.pi*np.array([286.43, 344.58, 488.27, 705.67])
    residues = np.stack([np.outer(cf[:, i], bf[i]).ravel(order="F")
                         for i in range(4)])
    return mass, residues/wn[:, None]**2, .02/wn, 1/wn**2, wn


def evaluate(f, mass, d, b0, b1, a1, a2, tau=TAU):
    s = 2j*np.pi*f
    a = 1+a1[:, None]*s+a2[:, None]*s**2
    flexible = ((b0[:, None, :]+s[None, :, None]*b1[:, None, :])
                /a[:, :, None]).sum(axis=0)
    phase = np.exp(-s*tau)
    p = mass[None, :]/s[:, None]**2 + d[None, :] + flexible
    return phase[:, None]*p, flexible, a


def system(f, mass, d, b0, b1, a1, a2, general=True, unknown_d=True,
           weighted=False):
    s = 2j*np.pi*f
    n = len(f)
    phase = np.exp(-s*TAU)
    mu, flexible, a = evaluate(f, mass, d, b0, b1, a1, a2)
    blocks, names = [], {}
    def add(name, value):
        start = sum(x.shape[2] for x in blocks)
        if value.ndim == 2:
            value = value[:, :, None]
        names[name] = slice(start, start+value.shape[2])
        blocks.append(value)
    def diagonal(v):
        return v[:, None, None]*np.eye(Q)[None, :, :]
    add("M", diagonal(phase/s**2))
    if unknown_d:
        add("D", diagonal(phase))
    for i in range(4):
        numerator = b0[i][None, :]+s[:, None]*b1[i][None, :]
        add(f"a1_{i}", -phase[:, None]*s[:, None]*numerator/a[i, :, None]**2)
        add(f"a2_{i}", -phase[:, None]*s[:, None]**2*numerator/a[i, :, None]**2)
        add(f"B0_{i}", diagonal(phase/a[i]))
        if general:
            add(f"B1_{i}", diagonal(phase*s/a[i]))
    j = np.concatenate(blocks, axis=2)
    # Optional non-diagonal, complex, frequency-dependent whitening stress test.
    l = np.broadcast_to(np.eye(Q, dtype=complex), (n, Q, Q)).copy()
    if weighted:
        l += .12j*np.diag(np.ones(Q-1), 1)[None, :, :]
        l *= (1+.25*f/f.max())[:, None, None]
    def whiten(value):
        if value.ndim == 3:
            return realify(np.einsum("nab,nbp->nap", l, value).reshape(n*Q, -1))
        return realify(np.einsum("nab,nb->na", l, value).ravel())
    x = whiten(j)
    scales = np.linalg.norm(x, axis=0)
    scaled = x/scales
    basis, r, _ = qr(scaled, mode="economic", pivoting=True)
    singular = svdvals(r)
    threshold = max(scaled.shape)*np.finfo(float).eps*singular[0]
    rank = int(np.count_nonzero(singular > threshold))
    # All test systems below are full-column-rank. Fail rather than hide rank loss.
    assert rank == x.shape[1], (rank, x.shape[1], singular[-1], threshold)
    def residual(v):
        return v-basis@(basis.T@v)
    t = whiten(-s[:, None]*mu)
    tf = whiten(-s[:, None]*phase[:, None]*flexible)
    tr = whiten(-phase[:, None]*mass[None, :]/s[:, None])
    td = whiten(-s[:, None]*phase[:, None]*d[None, :])
    out = {"f": f, "x": x, "j": j, "names": names, "residual": residual,
           "whiten": whiten, "t": t, "tf": tf, "tr": tr, "td": td,
           "rank": rank, "condition_scaled": float(singular[0]/singular[-1])}
    rt = residual(t)
    out["info"] = float(rt@rt)
    if general and unknown_d:
        direction = np.zeros(x.shape[1])
        direction[names["D"]] = -(b1/a2[:, None]).sum(axis=0)
        for i in range(4):
            direction[names[f"B0_{i}"]] = b1[i]/a2[i]
            direction[names[f"B1_{i}"]] = a1[i]*b1[i]/a2[i]-b0[i]
        out["identity_relative"] = float(np.linalg.norm(x@direction-tf)/np.linalg.norm(tf))
        out["flexible_residual_relative"] = float(np.linalg.norm(residual(tf))/np.linalg.norm(tf))
        out["reduced_residual_relative"] = float(
            np.linalg.norm(rt-residual(tr+td))/np.linalg.norm(t))
    return out


def finite_difference_checks(f, mass, d, b0, b1, a1, a2):
    mu, _, _ = evaluate(f, mass, d, b0, b1, a1, a2)
    h = 1e-8
    plus, _, _ = evaluate(f, mass, d, b0, b1, a1, a2, TAU+h)
    minus, _, _ = evaluate(f, mass, d, b0, b1, a1, a2, TAU-h)
    exact = -(2j*np.pi*f)[:, None]*mu
    delay_error = np.linalg.norm((plus-minus)/(2*h)-exact)/np.linalg.norm(exact)
    obj = system(f, mass, d, b0, b1, a1, a2, weighted=True)
    rng = np.random.default_rng(20261013)
    step = np.zeros(obj["x"].shape[1])
    parameters = {"M": mass, "D": d}
    parameters.update({f"a1_{i}": np.array([a1[i]]) for i in range(4)})
    parameters.update({f"a2_{i}": np.array([a2[i]]) for i in range(4)})
    parameters.update({f"B0_{i}": b0[i] for i in range(4)})
    parameters.update({f"B1_{i}": b1[i] for i in range(4)})
    for name, values in parameters.items():
        step[obj["names"][name]] = .1*np.linalg.norm(values)/np.sqrt(values.size)*rng.normal(size=values.size)
    def shifted(sign):
        values = {name: value+sign*1e-4*step[obj["names"][name]]
                  for name, value in parameters.items()}
        return evaluate(f, values["M"], values["D"],
                        np.stack([values[f"B0_{i}"] for i in range(4)]),
                        np.stack([values[f"B1_{i}"] for i in range(4)]),
                        np.array([values[f"a1_{i}"][0] for i in range(4)]),
                        np.array([values[f"a2_{i}"][0] for i in range(4)]))[0]
    numerical = obj["whiten"]((shifted(1)-shifted(-1))/(2e-4))
    analytic = obj["x"]@step
    return {"delay_jacobian_relative": float(delay_error),
            "nuisance_directional_jacobian_relative": float(np.linalg.norm(numerical-analytic)/np.linalg.norm(analytic))}


def pure_d_check(d):
    f = np.arange(5., 801.)
    s = 2j*np.pi*f
    phase = np.exp(-s*TAU)
    x = realify((phase[:, None, None]*np.eye(Q)[None]).reshape(len(f)*Q, Q))
    t = realify((-s[:, None]*phase[:, None]*d[None]).ravel())
    basis, _ = qr(x, mode="economic")
    r = t-basis@(basis.T@t)
    expected = 2*float(d@d)*float(np.sum((2*np.pi*f)**2))
    return {"projected_information": float(r@r), "closed_form_information": expected,
            "relative_error": float(abs(r@r-expected)/expected),
            "tangent_orthogonality_relative": float(np.linalg.norm(basis.T@t)/np.linalg.norm(t))}


def main():
    mass, b0, a1, a2, wn = geometry()
    f = np.arange(5., 801.)
    d_demo = (1e-8*np.eye(4)).ravel(order="F")
    b1_zero = np.zeros_like(b0)
    # Deliberately nonzero real first-order numerators: tests the stronger identity.
    b1_stress = .17*b0/wn[:, None]
    checks = finite_difference_checks(f[::3], mass, d_demo, b0, b1_stress, a1, a2)
    stress = system(f, mass, d_demo, b0, b1_stress, a1, a2, weighted=True)
    for key in ("identity_relative", "flexible_residual_relative", "reduced_residual_relative"):
        checks["nonzero_B1_weighted_"+key] = stress[key]
    zero_bg = system(f, np.zeros(Q), np.zeros(Q), b0, b1_stress, a1, a2)
    checks["M_D_zero_information_fraction"] = zero_bg["info"]/float(zero_bg["t"]@zero_bg["t"])
    # Same plant and same observations: adding D as nuisance must not add information.
    known_d = system(f, mass, d_demo, b0, b1_zero, a1, a2, unknown_d=False)
    unknown_d = system(f, mass, d_demo, b0, b1_zero, a1, a2)
    real_model = system(f, mass, d_demo, b0, b1_zero, a1, a2, general=False)
    inequalities = {"general_D_known": known_d["info"],
                    "general_D_unknown": unknown_d["info"],
                    "real_numerator_D_unknown": real_model["info"]}
    # For fixed remaining truth, X does not depend on true D: verify quadratic dependence.
    residual = unknown_d["residual"]
    u, v = residual(unknown_d["tr"]), residual(unknown_d["td"])
    coefficients = [float(u@u), float(2*u@v), float(v@v)]
    quadratic_rows = []
    for alpha in (-1., 0., 1., 2.):
        predicted = coefficients[0]+alpha*coefficients[1]+alpha**2*coefficients[2]
        actual = float(np.linalg.norm(residual(unknown_d["tf"]+unknown_d["tr"]+alpha*unknown_d["td"]))**2)
        quadratic_rows.append({"alpha": alpha, "information": actual,
                               "quadratic_relative_error": abs(actual-predicted)/actual})
    rows = []
    for d_label, d in (("zero", np.zeros(Q)), ("illustrative_1e-8_I4", d_demo)):
        for general in (False, True):
            local = []
            for low in (5, 20, 200):
                if d_label == "illustrative_1e-8_I4" and low == 5:
                    obj = unknown_d if general else real_model
                else:
                    obj = system(f[f >= low], mass, d, b0, b1_zero, a1, a2, general=general)
                row = {"D_truth": d_label, "general_numerator": general, "fmin_Hz": low,
                       "fmax_Hz": 800, "N_frequency": int(len(obj["f"])),
                       "information_for_sigma2_1": obj["info"], "nuisance_rank": obj["rank"],
                       "condition_after_column_scaling": obj["condition_scaled"]}
                if general:
                    for key in ("identity_relative", "flexible_residual_relative", "reduced_residual_relative"):
                        row[key] = obj[key]
                local.append(row)
            for row in local:
                row["information_fraction_of_5Hz"] = row["information_for_sigma2_1"]/local[0]["information_for_sigma2_1"]
            rows.extend(local)
    pure_d = pure_d_check(d_demo)
    assert checks["delay_jacobian_relative"] < 1e-8
    assert checks["nuisance_directional_jacobian_relative"] < 1e-6
    assert all(checks["nonzero_B1_weighted_"+key] < 1e-11 for key in
               ("identity_relative", "flexible_residual_relative", "reduced_residual_relative"))
    assert checks["M_D_zero_information_fraction"] < 1e-24
    assert unknown_d["info"] <= known_d["info"]*(1+1e-10)
    assert unknown_d["info"] <= real_model["info"]*(1+1e-10)
    assert max(x["quadratic_relative_error"] for x in quadratic_rows) < 1e-10
    assert pure_d["relative_error"] < 1e-12
    for row in rows:
        if row["general_numerator"]:
            assert row["flexible_residual_relative"] < 1e-11
    report = {"scope": "Deterministic four-mode Jacobian and Fisher geometry checks; not fitting or Monte Carlo.",
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "model": {"frequencies_Hz": (wn/(2*np.pi)).tolist(), "damping": .01,
                        "B0": "Ri / wi^2 using only the first four columns/rows of matlab/3/model.m",
                        "D_demo": "1e-8 I_4, arbitrary illustrative value, NOT measured/estimated or derived from mode 5",
                        "noise": "Proper complex Gaussian, covariance I for tables (sigma^2=1); all parameters real",
                        "weights_stress": "Non-diagonal complex per-frequency whitening, fixed and nonsingular",
                        "constraints": "Free M,D,B0 matrices; free B1 only in general model; four shared pole pairs",
                        "grid": "5:1:800 Hz, subsetted by each lower cutoff; positive frequencies only",
                        "tau_seconds": TAU},
              "algebra_checks": checks, "nuisance_space_comparisons": inequalities,
              "D_quadratic": {"coefficients_c0_c1_c2": coefficients, "checks": quadratic_rows},
              "pure_real_D_check": pure_d, "frequency_band_examples": rows,
              "all_assertions_passed": True,
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (ROOT/"matlab/3/model.m", Path(__file__))}}
    destination = Path(__file__).with_name("delay_information_with_D_checks.json")
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
