"""Read-only geometry checks for the project's redundant actuator/sensor model.

This is a deterministic algebra check, not an identification performance study.
It reads the nominal B1/B2/C1/C2 matrices; no MATLAB data or code is modified.
"""
import json
import pathlib
import re

import numpy as np


def main():
    root = pathlib.Path(__file__).resolve().parents[2]
    source = (root / "matlab/0/Multiple_main_Control.m").read_text(
        encoding="utf-8", errors="replace"
    )

    def matrix(name):
        text = re.search(r"\b" + name + r"\s*=\s*\[([^\]]+)\]", source).group(1)
        return np.array([
            [float(v) for v in row.replace("\n", " ").split(",")]
            for row in text.split(";") if row.strip()
        ])

    br, bf, cr, cf = [matrix(n) for n in ("B1", "B2", "C1", "C2")]
    _, _, vh = np.linalg.svd(br, full_matrices=True)
    q_u = vh.T
    q_y, _, _ = np.linalg.svd(cr, full_matrices=True)
    n, m = q_u[:, -1], q_y[:, -1]

    # Unique cross entries: last column (4) and last row of first 3 columns.
    select = np.eye(16)[[12, 13, 14, 15, 3, 7, 11]]
    cross = select @ np.kron(q_u.T, q_y.T)
    rigid_map = np.kron(br.T, cr)
    assert np.linalg.matrix_rank(cross) == 7
    assert np.linalg.matrix_rank(rigid_map) == 9
    assert np.linalg.norm(cross @ rigid_map) < 1e-12

    modes = []
    for i in range(5):
        left, right = cf[:, i], bf[i, :]
        residue = np.outer(left, right)
        alpha = abs(right @ n) / np.linalg.norm(right)
        beta = abs(m @ left) / np.linalg.norm(left)
        pivot = m @ residue @ n
        recovered = np.outer(residue @ n, m @ residue) / pivot
        err = np.linalg.norm(recovered - residue) / np.linalg.norm(residue)
        union = np.linalg.norm(cross @ residue.ravel(order="F")) / np.linalg.norm(residue)
        assert abs(union**2 - (alpha**2 + beta**2 - alpha**2 * beta**2)) < 1e-13
        assert err < 1e-12
        # Normalize factor units before comparing local Jacobian singular values.
        lc, rr = left / np.linalg.norm(left), right / np.linalg.norm(right)
        jac = np.concatenate((np.kron(rr[:, None], np.eye(4)),
                              np.kron(np.eye(4), lc[:, None])), axis=1)
        singular = np.linalg.svd(cross @ jac, compute_uv=False)
        assert np.linalg.matrix_rank(cross @ jac) == 7
        modes.append({
            "mode": i + 1, "input_null_amplitude_fraction": float(alpha),
            "output_null_amplitude_fraction": float(beta),
            "double_null_amplitude_fraction": float(alpha * beta),
            "cross_amplitude_fraction": float(union),
            "cross_energy_fraction": float(union**2),
            "normalized_factor_jacobian_condition": float(singular[0] / singular[-1]),
            "exact_rank1_recovery_relative_error": float(err),
        })

    # Profile-nuisance GLS identity, for a random positive-definite covariance.
    rng = np.random.default_rng(20261012)
    noise_factor = rng.normal(size=(16, 16))
    covariance = noise_factor @ noise_factor.T + np.eye(16)
    weight = np.linalg.inv(covariance)
    omega_profile = weight - weight @ rigid_map @ np.linalg.solve(
        rigid_map.T @ weight @ rigid_map, rigid_map.T @ weight
    )
    omega_cross = cross.T @ np.linalg.solve(cross @ covariance @ cross.T, cross)
    weight_error = np.linalg.norm(omega_profile - omega_cross) / np.linalg.norm(omega_cross)
    assert weight_error < 1e-11

    print(json.dumps({
        "scope": "Nominal spatial algebra only; no noise Monte Carlo or fitted modal performance.",
        "rigid_input_rank": int(np.linalg.matrix_rank(br)),
        "rigid_output_rank": int(np.linalg.matrix_rank(cr)),
        "flexible_input_rank": int(np.linalg.matrix_rank(bf)),
        "input_null_vector": n.tolist(), "output_null_vector": m.tolist(),
        "rigid_annihilation_norm": float(np.linalg.norm(cross @ rigid_map)),
        "modes": modes, "profile_gls_relative_discrepancy": float(weight_error),
    }, indent=2))


if __name__ == "__main__":
    main()
