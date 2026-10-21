"""Deterministic counterexample to monotonic RIV cost decrease under underfitting.

This standalone audit neither imports nor modifies the MATLAB project.
Requires NumPy. Run: python docs/research/check_chapter3_riv_step.py
"""

import json
import numpy as np


def main():
    omega = np.linspace(0.05, 6.0, 200)
    s = 1j * omega
    observed = 1.0 / (s**2 + 0.1 * s + 1.0) + 0.3 / (
        s**2 + 0.2 * s + 4.0
    )
    initial = np.array([0.0025, 1.5625, 0.01])
    denominator = 1.0 + initial[0] * s + initial[1] * s**2
    fitted = initial[2] / denominator

    # Columns are a1, a2, b; all unknown parameters are REAL.
    regressor = np.column_stack(
        (-s * observed / denominator, -s**2 * observed / denominator,
         1.0 / denominator)
    )
    instruments = np.column_stack(
        (-s * fitted / denominator, -s**2 * fitted / denominator,
         1.0 / denominator)
    )
    filtered_data = observed / denominator
    moment = np.real(instruments.conj().T @ regressor)
    rhs = np.real(instruments.conj().T @ filtered_data)
    updated = np.linalg.solve(moment, rhs)
    updated_fitted = updated[2] / (1.0 + updated[0] * s + updated[1] * s**2)
    cost_before = float(np.sum(np.abs(observed - fitted)**2))
    cost_after = float(np.sum(np.abs(observed - updated_fitted)**2))
    poles_before = np.roots([initial[1], initial[0], 1.0])
    poles_after = np.roots([updated[1], updated[0], 1.0])

    # Two independent ways to implement the complex-to-real equations agree.
    x_real = np.vstack((regressor.real, regressor.imag))
    m_real = np.vstack((instruments.real, instruments.imag))
    y_real = np.concatenate((filtered_data.real, filtered_data.imag))
    updated_real = np.linalg.solve(m_real.T @ x_real, m_real.T @ y_real)
    assert np.allclose(updated, updated_real, rtol=1e-11, atol=1e-11)
    assert np.all(poles_before.real < 0) and np.all(poles_after.real < 0)
    assert cost_after > 2 * cost_before

    # The local directional derivative and full-step cost are different checks.
    # In this example the direction is locally descending but the full step overshoots.
    direction = updated - initial
    gradient = -2 * np.real(instruments.conj().T @ (observed - fitted))
    report = {
        'initial_a1_a2_b': initial.tolist(),
        'updated_a1_a2_b': updated.tolist(),
        'initial_poles_real_imag': [[float(p.real), float(p.imag)] for p in poles_before],
        'updated_poles_real_imag': [[float(p.real), float(p.imag)] for p in poles_after],
        'unweighted_complex_cost_before': cost_before,
        'unweighted_complex_cost_after': cost_after,
        'cost_ratio': cost_after / cost_before,
        'gradient_dot_riv_direction': float(gradient @ direction),
        'complex_real_solution_discrepancy': float(np.linalg.norm(updated-updated_real)),
        'scope': 'A deterministic underfit example; not a Monte Carlo study or a proof about all RIV iterations.',
    }
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
