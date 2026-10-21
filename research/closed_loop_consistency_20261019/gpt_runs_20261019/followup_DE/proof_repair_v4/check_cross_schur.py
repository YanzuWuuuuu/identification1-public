"""Deterministic formula checks only; no Monte Carlo or engineering validation.

Run: python check_cross_schur.py
Writes cross_schur_checks.json next to this script.
"""
import json
from pathlib import Path
import numpy as np


def response_and_gradient(theta, z, ts=0.01):
    """Exact ZOH of G(s)=b*r/(s+r), theta=(log(b), log(r))."""
    b, rate = np.exp(theta)
    a = np.exp(-rate * ts)
    g = b * (1 - a) / (z - a)
    dg_dloga = b * (1 - z) / (z - a)**2 * (-rate * ts * a)
    return g, np.column_stack([g, dg_dloga])


def main():
    omega = np.linspace(0.02, 2.8, 12001)
    z = np.exp(1j * omega)
    theta0 = np.log([2.0, 10.0])
    true_k = 1.5
    assumed_k = 0.2
    weight = 1.0 + 0.35 * np.sin(omega)
    pd = 1.0 + 0.2 * np.cos(0.7 * omega)
    g0, dg0 = response_and_gradient(theta0, z)
    s0 = 1 / (1 + true_k * g0)

    def population_moment(theta, auxiliary_k):
        g, dg = response_and_gradient(theta, z)
        shat = 1 / (1 + auxiliary_k * g)
        instrument = dg * shat[:, None]
        mean_residual = (g0 - g) * s0
        return np.mean(np.real(instrument.conj() *
            (weight * pd * mean_residual)[:, None]), axis=0)

    rows = []
    for kh in [assumed_k, true_k, 0.0]:
        shat = 1 / (1 + kh * g0)
        zz = dg0 * shat[:, None]
        jd = dg0 * s0[:, None]
        cross = np.real(np.einsum('ki,k,kj->ij', zz.conj(), weight*pd, jd)) / len(z)
        gram = np.real(np.einsum('ki,k,kj->ij', zz.conj(), weight*pd, zz)) / len(z)
        jac = np.empty((2, 2))
        h = 1e-5
        for j in range(2):
            step = np.eye(2)[j] * h
            jac[:, j] = (population_moment(theta0+step, kh)-
                         population_moment(theta0-step, kh))/(2*h)
        cross_relerr = np.linalg.norm(jac+cross)/np.linalg.norm(cross)
        gram_relerr = np.linalg.norm(jac+gram)/np.linalg.norm(cross)
        assert cross_relerr < 1e-7
        if kh == true_k:
            assert np.allclose(cross, gram, rtol=1e-13, atol=1e-15)
        else:
            assert gram_relerr > 0.05
        rows.append(dict(true_K=true_k, auxiliary_K=kh,
            theta0=theta0.tolist(), cross_A=cross.tolist(), incorrect_Gram=gram.tolist(),
            finite_difference_Dh=jac.tolist(), relative_error_cross=cross_relerr,
            relative_error_Gram=gram_relerr, determinant_A=float(np.linalg.det(cross)),
            condition_A=float(np.linalg.cond(cross)),
            norm_h_at_truth=float(np.linalg.norm(population_moment(theta0, kh)))))

    # General nonsymmetric cross matrix: 3 structural and 2 nuisance parameters.
    # A=Z^T J/n, Q=Z^T diag(variance) Z/n, using a fixed random numerical fixture.
    rng = np.random.default_rng(20261019)
    zreg = rng.normal(size=(80, 5))
    jreg = zreg + .35*rng.normal(size=(80, 5))
    variance = np.exp(.4*rng.normal(size=80))
    aa = zreg.T@jreg/80
    qq = zreg.T@(variance[:, None]*zreg)/80
    aii = np.linalg.inv(aa)
    full = aii@qq@aii.T
    ntheta = 3
    att, ate = aa[:ntheta, :ntheta], aa[:ntheta, ntheta:]
    aet, aee = aa[ntheta:, :ntheta], aa[ntheta:, ntheta:]
    tmp = np.linalg.solve(aee.T, ate.T).T
    ae = att-tmp@aet
    ll = np.column_stack([np.eye(ntheta), -tmp])
    aei = np.linalg.inv(ae)
    schur = aei@(ll@qq@ll.T)@aei.T
    wrong = np.linalg.inv(att)@qq[:ntheta, :ntheta]@np.linalg.inv(att).T
    serr = np.linalg.norm(schur-full[:ntheta, :ntheta])/np.linalg.norm(full[:ntheta, :ntheta])
    werr = np.linalg.norm(wrong-full[:ntheta, :ntheta])/np.linalg.norm(full[:ntheta, :ntheta])
    assert serr < 1e-12
    assert werr > 1e-3
    result = dict(purpose="deterministic algebra check; not new MC or hardware validation",
        cross_checks=rows,
        schur_check=dict(seed=20261019, A=aa.tolist(), Q=qq.tolist(),
            full_sandwich_structural_block=full[:ntheta, :ntheta].tolist(),
            schur_structural_covariance=schur.tolist(),
            relative_error=serr,
            deleting_nuisance_columns_relative_error=werr),
        passed=True)
    path=Path(__file__).with_name('cross_schur_checks.json')
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(dict(cross_relative_errors=[r['relative_error_cross'] for r in rows],
        wrong_Gram_relative_errors=[r['relative_error_Gram'] for r in rows],
        schur_relative_error=serr, deleting_nuisance_error=werr), indent=2))


if __name__ == '__main__':
    main()
