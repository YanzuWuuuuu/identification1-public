"""Deterministic algebra checks only; not an identification-performance experiment."""
from pathlib import Path
import json
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / 'matlab/3/model.m').read_text(encoding='utf-8-sig')

def matrix(name):
    body = re.search(r'\b' + name + r'\s*=\s*\[([^\]]+)\]', source).group(1)
    return np.array([[float(x) for x in row.replace(',', ' ').split()]
                     for row in body.strip().split(';')])

C, B = matrix('C1'), matrix('B1')
Cflex, Bflex = matrix('C2')[:, :4], matrix('B2')[:4, :]
H0 = np.diag([1 / .0439, 1 / .0384, 1 / 4.2809])
M = C @ H0 @ B
ny = np.linalg.svd(C.T)[2][-1]
nu = np.linalg.svd(B)[2][-1]

def rank(a):
    return int(np.linalg.matrix_rank(a, tol=1e-10))

def lag_design(m):
    p, q = m.shape
    return np.column_stack(
        [(np.diag(np.eye(p)[i]) @ m).ravel() for i in range(p)] +
        [(m @ np.diag(np.eye(q)[j])).ravel() for j in range(q)])

def effective_rank(c, b):
    m = c @ H0 @ b
    nuisance = np.column_stack([(c @ np.eye(3)[i:i+1].T @
                                np.eye(3)[j:j+1] @ b).ravel()
                               for i in range(3) for j in range(3)])
    u, s, _ = np.linalg.svd(nuisance, full_matrices=False)
    u = u[:, s > 1e-10]
    lag = lag_design(m)
    return rank(lag - u @ (u.T @ lag))

Ay = np.column_stack([ny[i] * M[i, :] for i in range(4)])
Au = np.column_stack([nu[j] * M[:, j] for j in range(4)])
dy = np.array([0., 8., 17., 5.]) * 1e-6
du = np.array([0., 13., 6., 21.]) * 1e-6
H1 = np.array([[2., -1., .2], [.3, .8, -.1], [-.2, .1, .4]]) * 1e-4
L = C @ H1 @ B - np.diag(dy) @ M - M @ np.diag(du)
dyhat = np.r_[0., np.linalg.lstsq(Ay[:, 1:], -ny @ L, rcond=None)[0]]
duhat = np.r_[0., np.linalg.lstsq(Au[:, 1:], -L @ nu, rcond=None)[0]]

f = 37.
s = 2j * np.pi * f
qy = np.exp(-s * dy) / (1 + s * np.array([20., 31., 9., 25.]) * 1e-6)
qu = np.exp(-s * du) / (1 + s * np.array([17., 11., 23., 35.]) * 1e-6)
z = np.diag(qy) @ C @ (H0 + s * H1) @ B @ np.diag(qu)
U, sv, Vh = np.linalg.svd(z)
ell, rv = U[:, -1].conj(), Vh[-1, :].conj()
qry = ny * ell[0] / (ny[0] * ell)
qru = nu * rv[0] / (nu[0] * rv)

freqs = [1., 2., 5., 10., 20., 40.]
wn = 2 * np.pi * np.array([286.43, 344.58, 488.27, 705.67])
pollution = []
for f in freqs:
    s = 2j * np.pi * f
    qy = np.exp(-s * dy) / (1 + s * np.array([20., 31., 9., 25.]) * 1e-6)
    qu = np.exp(-s * du) / (1 + s * np.array([17., 11., 23., 35.]) * 1e-6)
    flex = Cflex @ np.diag(1 / (s*s + .02*wn*s + wn*wn)) @ Bflex
    z = np.diag(qy) @ (M/s**2 + flex + 1e-8*np.eye(4)) @ np.diag(qu)
    U, sv, Vh = np.linalg.svd(z)
    ell, rv = U[:, -1].conj(), Vh[-1, :].conj()
    ry = ny * ell[0] / (ny[0] * ell)
    ru = nu * rv[0] / (nu[0] * rv)
    pollution.append({'frequency_hz': f,
                      'relative_output_response_max_error': float(np.max(np.abs(ry-qy/qy[0]))),
                      'relative_input_response_max_error': float(np.max(np.abs(ru-qu/qu[0]))),
                      'sigma4_over_sigma3': float(sv[-1]/sv[-2])})

result = {
    'status': 'deterministic algebra only; artificial channel filters and delays; no noise, fitting or measured platform evidence',
    'geometry_source': 'matlab/3/model.m C1, B1',
    'left_null_vector': ny.tolist(), 'right_null_vector': nu.tolist(),
    'rank_rigid_matrix': rank(M),
    'sensor_relative_lag_rank': rank(Ay), 'actuator_relative_lag_rank': rank(Au),
    'strict_free_body_full_slope_rank_without_nuisance': rank(lag_design(M)),
    'effective_lag_ranks_with_arbitrary_3x3_rigid_slope_nuisance': {
        '4_sensors_4_actuators': effective_rank(C, B),
        '4_sensors_3_actuators': effective_rank(C, B[:, :3]),
        '3_sensors_4_actuators': effective_rank(C[:3, :], B),
        '3_sensors_3_actuators': effective_rank(C[:3, :], B[:, :3])},
    'slope_recovery_max_abs_error_seconds': {
        'sensor': float(np.max(np.abs(dyhat-dy))),
        'actuator': float(np.max(np.abs(duhat-du)))},
    'four_mode_plus_D_pollution': pollution,
    'caveat': 'D=1e-8*I and filter parameters are chosen checks, not estimates. An arbitrary rigid slope is a mathematical robustness class, not an established physical feature of this platform.'
}
# Recompute the earlier rigid-only target because the frequency loop changed qy/qu.
s = 2j*np.pi*37.
target_y = np.exp(-s*dy)/(1+s*np.array([20.,31.,9.,25.])*1e-6)
target_u = np.exp(-s*du)/(1+s*np.array([17.,11.,23.,35.])*1e-6)
result['rigid_only_complex_response_ratio_error'] = {
    'sensor': float(np.max(np.abs(qry-target_y/target_y[0]))),
    'actuator': float(np.max(np.abs(qru-target_u/target_u[0])))}
out = ROOT/'docs/research/redundant_calibration_algebra_checks_20261013.json'
out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=False, indent=2))
