"""Deterministic controller-design audit only: no MC, fitting or simulation."""
import os
import sys
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.update(NOFLEX='0', FC='60', FLP='120', FSC='3', AL='30', AH='100', SW='.03')
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import scipy.linalg as la
from scipy import signal, optimize

OUT = Path(__file__).resolve().parent
NEW = OUT.parents[1]
RESEARCH = NEW.parents[1]
sys.path.insert(0, str(RESEARCH / 'code'))
import cl2

FS = 5000.
NI = np.array([0, 1, 0, 1])
NO = np.array([1, 0, 1, 0])
TAU = (NO[:, None] + NI[None, :]) / FS
F = np.unique(np.r_[np.arange(3., 550.0001, .02), cl2.FM[:6], 3., 100., 550.])
POINTS = np.array([3., 10., 30., 60., 100., 286.573, 344.945, 490.39, 550.])


def plain(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    raise TypeError(type(x).__name__)


def plant(f, nm=6):
    f = np.atleast_1d(f)
    z = np.exp(2j*np.pi*f/FS)
    G = np.zeros((len(f), 4, 4), complex)
    for i in range(nm):
        h0 = cl2.zoh_2nd(2*np.pi*cl2.FM[i], cl2.ZM[i], z)[0]
        G += h0[:, None, None] * cl2.RES[i]
    return G*np.exp(-2j*np.pi*f[:, None, None]*TAU)


def kresp(f, fc=60.):
    A, B, C, D = cl2.controller(fc)
    f = np.atleast_1d(f)
    z = np.exp(2j*np.pi*f/FS)
    return C[None] @ np.linalg.solve(z[:, None, None]*np.eye(len(A))-A,
        np.broadcast_to(B, (len(z),)+B.shape)) + D


def rolloff_response(f, cutoff):
    f = np.atleast_1d(f)
    if cutoff is None: return np.ones(len(f), complex)
    sos = signal.butter(2, cutoff, fs=FS, output='sos')
    return signal.sosfreqz(sos, worN=f, fs=FS)[1]


def cascaded_controller(fc, cutoff):
    A, B, C, D = cl2.controller(fc)
    if cutoff is None: return A, B, C, D
    # Identical scalar digital low-pass on all four controller outputs.
    sos = signal.butter(2, cutoff, fs=FS, output='sos')
    an, bn = signal.sos2tf(sos)
    Af0, Bf0, Cf0, Df0 = signal.tf2ss(an, bn)
    Af = la.block_diag(*([Af0]*4))
    Bf = la.block_diag(*([Bf0]*4))
    Cf = la.block_diag(*([Cf0]*4))
    Df = np.eye(4)*Df0.item()
    AA = np.block([[A, np.zeros((len(A), len(Af)))], [Bf@C, Af]])
    BB = np.vstack([B, Bf@D])
    CC = np.hstack([Df@C, Cf])
    return AA, BB, CC, Df@D


def closed_loop(fc=60., cutoff=None, nm=6):
    A, B, C = cl2.ct_ss()
    A, B, C = A[:2*nm, :2*nm], B[:2*nm], C[:, :2*nm]
    Ad, Bd, Cd, _, _ = signal.cont2discrete((A, B, C, np.zeros((4, 4))), 1/FS)
    AK, BK, CK, DK = cascaded_controller(fc, cutoff)
    np_, nk = len(A), len(AK)
    n = np_ + sum(NI) + sum(NO) + nk
    ii, oo, off = [], [], np_
    for nn in NI:
        ii.append(np.arange(off, off+nn)); off += nn
    for nn in NO:
        oo.append(np.arange(off, off+nn)); off += nn
    def update(x):
        raw = Cd@x[:np_]
        y = np.array([x[ix[-1]] if len(ix) else raw[j] for j, ix in enumerate(oo)])
        u = -CK@x[-nk:] - DK@y
        up = np.array([x[ix[-1]] if len(ix) else u[j] for j, ix in enumerate(ii)])
        xn = np.zeros(n)
        xn[:np_] = Ad@x[:np_] + Bd@up
        xn[-nk:] = AK@x[-nk:] + BK@y
        for j, ix in enumerate(ii):
            if len(ix): xn[ix] = np.r_[u[j], x[ix[:-1]]]
        for j, ix in enumerate(oo):
            if len(ix): xn[ix] = np.r_[raw[j], x[ix[:-1]]]
        return xn
    Ac = np.column_stack([update(row) for row in np.eye(n)])
    poles = la.eigvals(Ac)
    worst = poles[np.argmax(abs(poles))]
    return dict(radius=float(abs(worst)), stable=bool(abs(worst)<1),
        dominant_pole=dict(real=float(worst.real), imag=float(worst.imag),
                          frequency_hz=float(abs(np.angle(worst))*FS/(2*np.pi))),
        unstable_pole_count=int(np.sum(abs(poles)>=1)), states=n)


def spectrum(f, fc=60., cutoff=None, nm=6):
    L = plant(f, nm) @ kresp(f, fc)
    h = rolloff_response(f, cutoff)
    return dict(sigma=np.linalg.svd(L, compute_uv=False)[:, 0]*abs(h),
                rho=np.max(abs(np.linalg.eigvals(L)), axis=1)*abs(h))


def refine_peak(f, vals, fun):
    k = np.argmax(vals)
    lb, ub = max(100., f[k]-.04), min(550., f[k]+.04)
    best = (float(vals[k]), float(f[k]))
    if ub > lb:
        opt = optimize.minimize_scalar(lambda x: -float(fun(x)), bounds=(lb, ub),
              method='bounded', options={'xatol': 1e-9})
        best = max(best, (-float(opt.fun), float(opt.x)))
    for end in (100., 550.): best = max(best, (float(fun(end)), end))
    return dict(value=best[0], frequency_hz=best[1])


def evaluate(fc, cutoff, nm=6):
    vals = spectrum(F, fc, cutoff, nm)
    hi = F>=100
    result = dict(controller_fc_hz=fc, added_rolloff_fc_hz=cutoff, modes=nm,
        closed_loop=closed_loop(fc, cutoff, nm))
    for key in ('sigma', 'rho'):
        result['high_max_'+key] = refine_peak(F[hi], vals[key][hi],
            lambda ff: spectrum([ff], fc, cutoff, nm)[key][0])
    h = rolloff_response(POINTS, cutoff)
    oldK = kresp(POINTS, 60.)
    newK = kresp(POINTS, fc)*h[:, None, None]
    ratio = np.sum(newK*oldK.conj(), axis=(1, 2))/np.sum(abs(oldK)**2, axis=(1, 2))
    points = spectrum(POINTS, fc, cutoff, nm)
    result['points'] = [dict(frequency_hz=float(f),
        rolloff_magnitude_db=float(20*np.log10(abs(hh))),
        rolloff_phase_deg=float(np.angle(hh)*180/np.pi),
        total_controller_change_from_fc60_db=float(20*np.log10(abs(ratio[i]))),
        total_controller_change_from_fc60_deg=float(np.angle(ratio[i])*180/np.pi),
        sigma_GK=float(points['sigma'][i]), rho_GK=float(points['rho'][i]))
        for i, (f, hh) in enumerate(zip(POINTS, h))]
    result['high_condition_pass'] = result['high_max_sigma']['value'] < .1
    # Only a quantitative description: no invented tolerance for "not obvious".
    low = (F>=3)&(F<=100)
    lh = rolloff_response(F[low], cutoff)
    result['rolloff_change_3_100'] = dict(
        largest_attenuation_db=float(-np.min(20*np.log10(abs(lh)))),
        largest_added_phase_lag_deg=float(-np.min(np.unwrap(np.angle(lh))*180/np.pi)))
    return result, vals


def main():
    deps = [RESEARCH/'code/cl2.py', NEW/'delay_sim_closeout/code/sim_core.py',
            NEW/'delay_sim_closeout/audit/stability_identifiability_audit.py']
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in deps}
    (OUT/'dependency_hashes_before.json').write_text(json.dumps(hashes, indent=2), encoding='utf8')
    results = []
    for cutoff in (None, 500., 300., 200., 150., 120., 100., 80., 60., 50., 45., 40., 39., 35., 30., 20.):
        row, vals = evaluate(60., cutoff)
        results.append(row)
        print('FC60', cutoff, row['high_max_sigma'], row['closed_loop'], flush=True)
    # This deterministic scalar bisection locates a design boundary, not a data fit.
    base = spectrum(F, 60., None)['sigma']
    hi = F>=100
    upper = optimize.brentq(lambda x: np.max(base[hi]*abs(rolloff_response(F[hi], x)))-.1,
                           5., 150., xtol=1e-10)
    row, _ = evaluate(60., float(upper))
    results.append(row)
    # Dense deterministic stability scan tests whether a sharper cutoff can help.
    scan = []
    for cutoff in np.arange(5., 501., 1.):
        stability = closed_loop(60., float(cutoff))
        mx = float(np.max(base[hi]*abs(rolloff_response(F[hi], cutoff))))
        scan.append(dict(cutoff_hz=float(cutoff), high_max_sigma_grid=mx, **stability))
    alternatives = []
    for fc in (15., 20., 25., 30.):
        for cutoff in (None, 150., 100., 60.):
            alt, _ = evaluate(fc, cutoff)
            alternatives.append(alt)
    # Validate the delayed-state construction against the existing eight-mode audit.
    check8, _ = evaluate(60., None, nm=8)
    # Independent checks: cascade transfer equals scalar product; augmented open
    # plant with integer delays equals the modal-frequency construction.
    A, B, C = cl2.ct_ss()
    A, B, C = A[:12, :12], B[:12], C[:, :12]
    Ad, Bd, Cd, _, _ = signal.cont2discrete((A, B, C, np.zeros((4, 4))), 1/FS)
    z = np.exp(2j*np.pi*POINTS/FS)
    gs = Cd[None]@np.linalg.solve(z[:, None, None]*np.eye(12)-Ad,
         np.broadcast_to(Bd, (len(z),)+Bd.shape))
    gs *= np.exp(-2j*np.pi*POINTS[:, None, None]*TAU)
    gp = plant(POINTS)
    cascade_errors = []
    for fc, cut in ((60., 150.), (60., upper), (15., None), (20., 100.)):
        ak, bk, ck, dk = cascaded_controller(fc, cut)
        ks = ck[None]@np.linalg.solve(z[:, None, None]*np.eye(len(ak))-ak,
              np.broadcast_to(bk, (len(z),)+bk.shape))+dk
        kp = kresp(POINTS, fc)*rolloff_response(POINTS, cut)[:, None, None]
        cascade_errors.append(dict(fc=fc, cutoff=cut,
            relative_error=float(np.linalg.norm(ks-kp)/np.linalg.norm(kp))))
    validation = dict(plant_state_vs_modal_relative_error=float(np.linalg.norm(gs-gp)/np.linalg.norm(gp)),
        cascade_state_vs_scalar_product=cascade_errors,
        eight_mode_original_radius=check8['closed_loop']['radius'],
        eight_mode_original_sigma=check8['high_max_sigma']['value'])
    np.savez_compressed(OUT/'six_mode_definition.npz', A=A, B=B, C=C,
        Ad=Ad, Bd=Bd, Cd=Cd, frequencies_hz=cl2.FM[:6], damping=cl2.ZM[:6],
        C0=np.asarray(cl2.RES[:6]), NI=NI, NO=NO, tau=TAU)
    original_sigma100 = float(spectrum([100.])['sigma'][0])
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        scope='No Monte Carlo, no parameter fitting, no time-domain simulation.',
        configuration=dict(fs=FS, original_fc_hz=60., original_FLP_hz=120.,
            mode_frequencies_hz=cl2.FM[:6], mode_damping=cl2.ZM[:6], NI=NI, NO=NO,
            delays_seconds=TAU, rolloff='second-order digital Butterworth, scipy.signal.butter, fs=5000',
            frequency_grid_hz=.02),
        source_hashes=hashes, original_and_rolloff=results,
        exact_grid_limit_cutoff_hz=upper,
        necessary_attenuation_at_100_db=20*np.log10(.1/original_sigma100),
        cutoff_scan=scan,
        feasible_stable_and_small_loop_cutoffs=[x for x in scan if x['stable'] and x['high_max_sigma_grid']<.1],
        low_band_change_tolerance='Not specified by the user; no pass criterion invented.',
        validation=validation,
        lower_controller_fc_alternatives=alternatives, existing_eight_mode_check=check8,
        source_hashes_after={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
    (OUT/'controller_check.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, default=plain), encoding='utf8')
    np.savez_compressed(OUT/'controller_base_grid.npz', frequency_hz=F,
                        sigma_GK=base, rho_GK=spectrum(F)['rho'])
    print('Small-loop cutoff boundary:', upper, 'Hz', flush=True)
    print('Feasible stable and small-loop cutoffs:', result['feasible_stable_and_small_loop_cutoffs'], flush=True)


if __name__ == '__main__': main()
