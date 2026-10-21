"""Read-only audit of archived experiment signals; recover candidate external input.

Run with system Python 3.11 (SciPy installed). Raw MAT files are never rewritten.
The exporter has invalid zlib Adler trailers; full raw-DEFLATE blocks are parsed
and their sizes, finite values, timestamps, and archived statistics are checked.
This recovery is an algebraic evidence audit, not a proof of noise exogeneity.
"""
from pathlib import Path
import hashlib
import io
import json
import struct
import zlib
import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import csd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
RAW = ROOT / "experiment_data"


def read_export(path):
    source = path.read_bytes()
    pos, arrays, blocks = 128, {}, []
    while pos + 8 <= len(source):
        tag, n = struct.unpack("<II", source[pos:pos + 8])
        assert tag == 15, (pos, tag)
        data = source[pos + 8:pos + 8 + n]
        assert len(data) == n
        decoder = zlib.decompressobj(-15)
        unpacked = decoder.decompress(data[2:-4]) + decoder.flush()
        assert decoder.eof and not decoder.unused_data
        inner_tag, inner_len = struct.unpack("<II", unpacked[:8])
        assert inner_tag == 14 and len(unpacked) == 8 + inner_len
        fields = loadmat(io.BytesIO(source[:128] + unpacked))
        real_fields = {k: v for k, v in fields.items() if not k.startswith("__")}
        assert len(real_fields) == 1
        name, value = next(iter(real_fields.items()))
        assert value.shape == (305001, 2) and np.isfinite(value).all()
        arrays[name] = value
        blocks.append({"name": name, "unpacked_bytes": len(unpacked),
                       "trailer_adler": data[-4:].hex(),
                       "computed_adler": zlib.adler32(unpacked).to_bytes(4, "big").hex()})
        pos += n + 8
    assert pos == len(source)
    return arrays, {"filename": str(path), "sha256": hashlib.sha256(source).hexdigest(),
                    "bytes": len(source), "blocks": blocks}


def matrices(a):
    time = a["u_Motor1"][:, 0]
    for value in a.values():
        assert np.array_equal(time, value[:, 0])
    return (time, np.column_stack([a[f"u_Motor{k}"][:, 1] for k in range(1, 5)]),
            np.column_stack([a[k][:, 1] for k in ("u_z", "u_Rx", "u_Ry")]),
            np.column_stack([a[f"pre_y_{k}__0_"][:, 1] for k in range(4)]))


a, meta = read_export(RAW / "auto_identify_motor1234_z0.mat")
time, U, V, Y = matrices(a)
t = time - time[0]
quiet = U[:, 3] == 0
train = quiet & (t < 30)
check = quiet & (t >= 30)
# Estimate only from first-half intervals where channel 4 has no input.
# No intercept: constant actuator loads are already included in V.
Q = np.linalg.lstsq(V[train], U[train], rcond=None)[0].T
D = U - V @ Q.T
assert np.max(abs(D[train])) < 1e-12
assert np.max(abs(D[check])) < 1e-12

archive = json.loads((RAW / "analysis_20261015/spectral_audit.json").read_text(encoding="utf-8"))
arch = next(v for v in archive if v["filename"] == "auto_identify_motor1234_z0.mat")
np.testing.assert_allclose(U.mean(0), arch["input_mean"], rtol=1e-11, atol=1e-12)
np.testing.assert_allclose(U.std(0, ddof=1), arch["input_std"], rtol=1e-11, atol=1e-12)
np.testing.assert_allclose(Y.std(0, ddof=1), arch["output_std"], rtol=1e-11, atol=1e-15)

Z = np.column_stack([a[f"y_0__{k}_"][:, 1] for k in range(3)])
P = np.linalg.lstsq(Y, Z, rcond=None)[0].T
vf = np.column_stack([U, np.ones(len(U))])
VfromU = vf @ np.linalg.lstsq(vf, V, rcond=None)[0]
periodic = []
for lag in (30000, 150000):
    i = np.arange(300000 - lag)
    periodic.append({"lag_samples": lag, "lag_s": lag / 5000,
                     "D_difference_rms_relative": (np.linalg.norm(D[i + lag] - D[i], axis=0) / np.linalg.norm(D[i], axis=0)).tolist(),
                     "D_max_abs_difference": np.max(abs(D[i + lag] - D[i]), axis=0).tolist(),
                     "V_difference_rms": np.sqrt(np.mean((V[i + lag] - V[i])**2, axis=0)).tolist()})

spectral = []
for name, sl in (("first_30_s", slice(5000, 155000)), ("second_29_s", slice(155000, 300000))):
    d = D[sl]
    S = np.zeros((5001, 4, 4), complex)
    for j in range(4):
        for k in range(4):
            f, S[:, j, k] = csd(d[:, k], d[:, j], fs=5000, nperseg=10000,
                                noverlap=5000, detrend="linear")
    scale = np.sqrt(np.real(np.diagonal(S, axis1=1, axis2=2)))
    Sn = S / scale[:, :, None] / scale[:, None, :]
    ev = np.linalg.eigvalsh(Sn)
    for lo, hi in ((3, 25), (25, 80), (80, 170), (260, 310), (25, 800)):
        ix = (f >= lo) & (f <= hi)
        cond = ev[ix, -1] / ev[ix, 0]
        spectral.append({"split": name, "band_hz": [lo, hi],
                         "normalized_eigenvalue_minimum": float(ev[ix, 0].min()),
                         "normalized_condition_median": float(np.median(cond)),
                         "normalized_condition_maximum": float(cond.max())})

cross = []
for name in ("auto_identify_motor1234_z-0.5.mat", "auto_identify_motor1234_z0.5.mat", "auto_identify_motor1234_z1.mat"):
    b, bmeta = read_export(RAW / name)
    tb, ub, vb, _ = matrices(b)
    db = ub - vb @ Q.T
    qb = ub[:, 3] == 0
    # Align one period using the exact channel-4 sequence, then check all four.
    x, y = D[:30000, 3], db[:30000, 3]
    cor = np.fft.ifft(np.fft.fft(x) * np.conj(np.fft.fft(y))).real
    shift = int(np.argmax(cor))
    aligned = np.roll(db[:30000], shift, axis=0)
    diff = abs(aligned - D[:30000])
    mismatch = np.max(diff, axis=1) > 1e-10
    cross.append({"filename": name, "sha256": bmeta["sha256"],
                  "quiet_max_abs_residual": np.max(abs(db[qb]), axis=0).tolist(),
                  "alignment_roll_samples": shift,
                  "aligned_period_max_abs_difference": np.max(diff, axis=0).tolist(),
                  "aligned_samples_different_above_1e_10": int(mismatch.sum()),
                  "aligned_mismatch_sample_indices_0_based": np.flatnonzero(mismatch).tolist(),
                  "other_samples_max_abs_difference": float(np.max(diff[~mismatch]))})

summary = {"source": meta, "sample_rate_hz": float(1 / np.median(np.diff(time))),
           "signal_names": list(a), "allocation_Q_motor_by_z_Rx_Ry": Q.tolist(),
           "logical_output_P_z_Rx_Ry_by_sensor": P.tolist(),
           "logical_output_max_error": np.max(abs(Y @ P.T - Z), axis=0).tolist(),
           "quiet_train_count": int(train.sum()), "quiet_check_count": int(check.sum()),
           "quiet_train_max_abs_residual": np.max(abs(D[train]), axis=0).tolist(),
           "quiet_check_max_abs_residual": np.max(abs(D[check]), axis=0).tolist(),
           "active_D_std": D[~quiet].std(0).tolist(),
           "active_D_correlation": np.corrcoef(D[~quiet].T).tolist(),
           "V_not_U_projection_relative_rms": ((VfromU - V).std(0) / V.std(0)).tolist(),
           "d_fields_nonzero_counts": {k: int(np.count_nonzero(a[k][:, 1])) for k in ("d_z", "d_Rx", "d_Ry")},
           "periodic_input_audit": periodic, "spectral_coverage_audit": spectral,
           "cross_position_fixed_allocation_check": cross,
           "caveat": "Strong algebraic evidence for archived pre-allocation feedback commands and additive repeatable excitation; recording block locations still need hardware provenance confirmation. Exact one-frequency repeated-period DFT remains rank 1; Welch local averaging rank must not be confused with independent spectral realizations."}
(OUT / "reference_recovery_audit.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
savemat(OUT / "recovered_reference_z0_candidate.mat", {
    "time": time[:, None], "time_relative": t[:, None], "fs": 5000.,
    "U": U, "Y": Y, "V_logic": V, "D_candidate": D, "Q_allocation": Q,
    "P_logical_output": P, "quiet_training_mask": train[:, None],
    "quiet_check_mask": check[:, None], "source_sha256": meta["sha256"],
    "status": "Recovered additive repeated input candidate; hardware block provenance unconfirmed.",
    "period_samples": 30000}, do_compression=True)
print(json.dumps({k: v for k, v in summary.items() if k not in ("source", "spectral_coverage_audit")}, indent=2))
