"""Independent algebra and information audit; reads, but never changes, MATLAB files."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from scipy.linalg import qr

ROOT = Path(__file__).resolve().parents[2]


def realify(x):
    return np.concatenate((x.real, x.imag), axis=0)


def geometry():
    text = (ROOT / "matlab/3/model.m").read_text(encoding="utf-8")
    def matrix(name):
        value = re.search(r"\b" + name + r"\s*=\s*\[([^\]]+)\]", text).group(1)
        return np.array([[float(x) for x in row.replace("\n", " ").split(",")]
                         for row in value.split(";") if row.strip()])
    br, bf, cr, cf = [matrix(x) for x in ("B1", "B2", "C1", "C2")]
    mass = cr @ np.diag([1/.0439, 1/.0384, 1/4.2809]) @ br
    residues = np.stack([np.outer(cf[:, i], bf[i]).ravel(order="F")
                         for i in range(5)], axis=1)
    return mass.ravel(order="F"), residues


def information(fmin, general, m, residues):
    f = np.arange(5., 801., 2.)
    f = f[f >= fmin]
    s = 2j*np.pi*f
    wn = 2*np.pi*np.array([286.43, 344.58, 488.27, 705.67, 738.83])
    a1, a2 = .02/wn, 1/wn**2
    q, n = len(m), len(f)
    eye = np.eye(q)
    blocks = [np.einsum("n,ab->nab", 1/s**2, eye).reshape(n*q, q)]
    gf = np.zeros((q, n), complex)
    coefficient = [np.zeros(q)]
    for i in range(5):
        b0 = residues[:, i]*a2[i]
        a = 1+a1[i]*s+a2[i]*s**2
        p = b0[:, None]/a
        gf += p
        block = [(-s*p/a).T.reshape(-1, 1), (-s**2*p/a).T.reshape(-1, 1),
                 np.einsum("n,ab->nab", 1/a, eye).reshape(n*q, q)]
        c = [np.zeros(2+q)]
        if general:
            block.append(np.einsum("n,ab->nab", s/a, eye).reshape(n*q, q))
            c.append(-b0)
        blocks.append(np.concatenate(block, axis=1))
        coefficient.append(np.concatenate(c))
    j = np.concatenate(blocks, axis=1)
    jr = realify(j)
    scale = np.linalg.norm(jr, axis=0)
    basis, _ = qr(jr/scale, mode="economic")
    g = m[:, None]/s**2+gf
    dt = realify((-s*g).T.reshape(-1))
    drb = realify((-m[:, None]/s).T.reshape(-1))
    dflex = realify((-s*gf).T.reshape(-1))
    residual = lambda x: x-basis@(basis.T@x)
    total, rigid, flex = [residual(x) for x in (dt, drb, dflex)]
    out = dict(fmin=fmin, general_numerator=general, N=n,
               profiled_norm_squared=float(total@total),
               removed_fraction=float(1-(total@total)/(dt@dt)),
               flexible_residual_relative=float(np.linalg.norm(flex)/np.linalg.norm(dflex)),
               full_vs_rigid_residual_relative=float(np.linalg.norm(total-rigid)/np.linalg.norm(dt)))
    if general:
        c = np.concatenate(coefficient)
        out["exact_flexible_tangent_identity_relative"] = float(
            np.linalg.norm(jr@c-dflex)/np.linalg.norm(dflex))
    return out


def block_identity():
    rng = np.random.default_rng(3159)
    s = 1j*np.linspace(.3, 5., 50)
    q, k = 2, 2
    theta = [np.array([.4, -.2]), np.array([.15, .7, 1., .3]),
             np.array([.2, .13, -.4, .9])]
    a = [1+t[0]*s+t[1]*s**2 for t in theta[1:]]
    p = [theta[i+1][2:, None]/a[i] for i in range(k)]
    g = theta[0][:, None]/s**2+sum(p)
    g += .03*(rng.normal(size=g.shape)+1j*rng.normal(size=g.shape))
    xrows, zrows, yrows, erows = [], [], [], []
    e = g-theta[0][:, None]/s**2-sum(p)
    block_sizes = [q]+[q+2]*k
    diagonal = np.zeros((sum(block_sizes), k+1))
    offset = 0
    for i, size in enumerate(block_sizes):
        diagonal[offset:offset+size, i] = theta[i]
        offset += size
    for n, sn in enumerate(s):
        # Non-diagonal Hermitian positive-definite metric via a whitening factor.
        factor = np.array([[1+.01*n, .2j], [.1, .8]])
        xs, zs = [np.eye(q)/sn**2], [np.eye(q)/sn**2]
        ys = [g[:, n]-sum(pi[:, n] for pi in p)]
        for i in range(k):
            gt = e[:, n]+p[i][:, n]
            xs.append(np.column_stack((-sn*gt/a[i][n], -sn**2*gt/a[i][n], np.eye(q)/a[i][n])))
            zs.append(np.column_stack((-sn*p[i][:, n]/a[i][n],
                                       -sn**2*p[i][:, n]/a[i][n], np.eye(q)/a[i][n])))
            ys.append(gt/a[i][n])
        xrows.append(factor@np.concatenate(xs, axis=1))
        zrows.append(factor@np.concatenate(zs, axis=1))
        yrows.append(factor@np.column_stack(ys))
        erows.append(factor@e[:, n])
    x, z, y = map(lambda rows: realify(np.concatenate(rows)), (xrows, zrows, yrows))
    er = realify(np.concatenate(erows))
    score = z.T@er
    delta = np.linalg.solve(z.T@x, score)
    solution = np.linalg.solve(z.T@x, z.T@y)
    reconstructed = diagonal+delta[:, None]
    return dict(joint_rhs_identity_relative=float(np.linalg.norm(y-x@diagonal-er[:, None])/np.linalg.norm(y)),
                joint_update_identity_relative=float(np.linalg.norm(solution-reconstructed)/np.linalg.norm(solution)),
                note="Identity checked off a solution, with nonzero residual and non-diagonal whitening.")


def rational_counterexample():
    s = 1j*np.linspace(.01, 12., 250)
    pair = lambda u, v: ((u+v)/(u*v)+(2/(u*v))*s)/(1+(u+v)/(u*v)*s+s*s/(u*v))
    p1, p2 = pair(1, 2)+pair(3, 4), pair(1, 3)+pair(2, 4)
    return dict(relative_transfer_discrepancy=float(np.linalg.norm(p1-p2)/np.linalg.norm(p1)),
                first_natural_frequencies=[2**.5, 12**.5],
                second_natural_frequencies=[3**.5, 8**.5],
                poles=[-1, -2, -3, -4],
                note="Both have positive denominator coefficients and distinct poles, but different overdamped pairings.")


def main():
    m, residues = geometry()
    rows = [information(lo, general, m, residues)
            for general in (False, True) for lo in (5, 10, 20, 50, 100, 200)]
    for general in (False, True):
        matching = [x for x in rows if x["general_numerator"] == general]
        base = matching[0]["profiled_norm_squared"]
        for x in matching:
            x["information_relative_to_5Hz"] = x["profiled_norm_squared"]/base
    report = dict(scope="Independent deterministic audit, not a MATLAB rerun or Monte Carlo.",
                  joint_fixed_point_algebra=block_identity(),
                  information=rows,
                  parameter_identifiability_counterexample=rational_counterexample(),
                  data_weight_counterexample={"weight":"1+Y^2; Y=mu+N(0,1)",
                     "pseudo_true_when_mu_1":5/3, "pseudo_true_when_mu_0":0,
                     "note":"Data dependence can cause bias, but dependence alone does not imply bias."},
                  source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in [ROOT/"数学推导_时延增广频域模态辨识.md",
                              ROOT/"matlab/research_experiments/theory_check.m",
                              ROOT/"matlab/research_experiments/theory_check2.m",
                              ROOT/"matlab/research_experiments/joint_srivc2.m"]})
    assert report["joint_fixed_point_algebra"]["joint_rhs_identity_relative"] < 1e-12
    assert report["joint_fixed_point_algebra"]["joint_update_identity_relative"] < 1e-10
    assert report["parameter_identifiability_counterexample"]["relative_transfer_discrepancy"] < 1e-12
    for row in rows:
        if row["general_numerator"]:
            assert row["exact_flexible_tangent_identity_relative"] < 1e-12
            assert row["flexible_residual_relative"] < 1e-11
            assert row["full_vs_rigid_residual_relative"] < 1e-11
    out = Path(__file__).with_name("delay_theory_review_checks.json")
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
