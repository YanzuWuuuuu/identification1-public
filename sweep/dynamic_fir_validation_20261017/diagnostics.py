"""Training-only ideal 3+1 cancellation diagnostic; never load Hcheck.

The frequency-by-frequency inverses below are not implemented filters.
Only diagnostics.py, diagnostic_limits.csv and diagnostic_limits.md belong
to this independent audit. Existing data/candidate files remain unchanged.
"""
from pathlib import Path
import csv
import numpy as np
from scipy.io import loadmat

OUT = Path(__file__).resolve().parent
BANDS = [(25., 80.), (260., 310.), (320., 380.), (450., 800.), (3., 800.)]


def physical(H):
    # Exact root-reviewed output-fast conversion, no conjugation.
    return H.reshape(len(H), 4, 4).transpose(0, 2, 1)


def energy(x, f):
    p = np.sum(np.abs(x.reshape(len(f), -1)) ** 2, axis=1)
    return float(np.trapezoid(p, f))


def project(H, Ty, Tu):
    return np.einsum("ab,kbc,cd->kad", Ty, physical(H), Tu)


def diagnostic(source, f, H, Ty, Tu, Py, Pu, ny, nu):
    P = project(H, Ty, Tu)
    A, b, c, d = P[:, :3, :3], P[:, :3, 3], P[:, 3, :3], P[:, 3, 3]
    # Solve actual inverses, not pseudo-inverses or regularized surrogates.
    Fu = np.linalg.solve(A, b[..., None])[..., 0]
    Fy = np.linalg.solve(A.transpose(0, 2, 1), c[..., None])[..., 0]
    S = d - np.sum(c * Fu, axis=1)
    singular = np.linalg.svd(A, compute_uv=False)
    cond = singular[:, 0] / singular[:, -1]
    pu = nu[None, :] - Fu @ Pu.T
    py = ny[None, :] - Fy @ Py
    gu, gy = np.linalg.norm(pu, axis=1), np.linalg.norm(py, axis=1)

    R = np.broadcast_to(np.eye(4, dtype=complex), P.shape).copy()
    L = R.copy()
    R[:, :3, 3] = -Fu
    L[:, 3, :3] = -Fy
    actual = L @ P @ R
    full_u = np.linalg.svd(Tu[None, :, :] @ R, compute_uv=False)[:, 0]
    full_y = np.linalg.svd(L @ Ty[None, :, :], compute_uv=False)[:, 0]
    full_u /= np.linalg.norm(Tu, 2)
    full_y /= np.linalg.norm(Ty, 2)
    rows = []
    for lo, hi in BANDS:
        ix = np.flatnonzero((f >= lo - 1e-9) & (f <= hi + 1e-9))
        ff = f[ix]
        iu = ix[np.argmax(gu[ix])]
        iy = ix[np.argmax(gy[ix])]
        ic = ix[np.argmax(cond[ix])]
        denergy = energy(d[ix], ff)
        assert denergy > 0
        er_u = np.sqrt(energy(actual[ix, :3, 3], ff) / energy(b[ix], ff))
        er_y = np.sqrt(energy(actual[ix, 3, :3], ff) / energy(c[ix], ff))
        schur_check = np.sqrt(energy(actual[ix, 3, 3] - S[ix], ff) / denergy)
        # Relative residual check ensures no transpose/phase/layout mistake.
        assert max(er_u, er_y, schur_check) < 1e-8
        rows.append(dict(
            source=source, fmin_Hz=lo, fmax_Hz=hi, samples=len(ix),
            ideal_actual_N_input_gain_max=float(gu[iu]),
            input_gain_max_frequency_Hz=float(f[iu]),
            ideal_actual_N_output_gain_max=float(gy[iy]),
            output_gain_max_frequency_Hz=float(f[iy]),
            ideal_full_input_operator_gain_ratio_max=float(np.max(full_u[ix])),
            ideal_full_output_operator_gain_ratio_max=float(np.max(full_y[ix])),
            Fu_norm_max=float(np.max(np.linalg.norm(Fu[ix], axis=1))),
            Fy_norm_max=float(np.max(np.linalg.norm(Fy[ix], axis=1))),
            Schur_relative_complex_change=np.sqrt(energy(S[ix]-d[ix], ff)/denergy),
            Schur_RMS_ratio_to_static_NN=np.sqrt(energy(S[ix], ff)/denergy),
            A_condition_min=float(np.min(cond[ix])),
            A_condition_median=float(np.median(cond[ix])),
            A_condition_max=float(cond[ic]),
            A_condition_max_frequency_Hz=float(f[ic]),
            A_smallest_singular_value_min=float(np.min(singular[ix, -1])),
            ideal_u_cross_relative_residual=float(er_u),
            ideal_y_cross_relative_residual=float(er_y),
            ideal_Schur_algebra_relative_residual=float(schur_check),
        ))
    return rows


def main():
    # Strict variable allow-list: check/validation data never enter memory.
    raw = loadmat(OUT / "input_z0.mat", squeeze_me=True,
                  variable_names=["f", "Htrain", "Hmodel"])
    cc = loadmat(OUT / "coordinates.mat", squeeze_me=True,
                 variable_names=["Py", "Pu", "ny", "nu"])
    f = np.asarray(raw["f"], dtype=float)
    Py, Pu = np.asarray(cc["Py"]), np.asarray(cc["Pu"])
    ny, nu = np.asarray(cc["ny"]).ravel(), np.asarray(cc["nu"]).ravel()
    Ty, Tu = np.vstack([Py, ny]), np.column_stack([Pu, nu])
    assert np.all(np.diff(f) > 0) and f.size == 1771
    assert np.isclose(np.linalg.norm(ny), 1) and np.isclose(np.linalg.norm(nu), 1)
    rows = []
    for source, name in [("nonparam_training", "Htrain"), ("SRIVC_training_model", "Hmodel")]:
        rows.extend(diagnostic(source, f, raw[name], Ty, Tu, Py, Pu, ny, nu))
    with (OUT / "diagnostic_limits.csv").open("w", encoding="utf-8-sig", newline="") as h:
        wr = csv.DictWriter(h, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)

    lines = [
        "# 训练数据上的理想逐频消元：目标响应与端口增益代价",
        "",
        "日期：2026-10-17。这里只加载 input_z0.mat 的 f、Htrain、Hmodel；"
        "coordinates.mat 只加载 Py、Pu、ny、nu。未读取检查数据，未操作硬件。",
        "",
        "所有计算使用同一强静态坐标 Ty=[Py;nyᵀ]、Tu=[Pu,nu]，"
        "不重拟合方向、不改变端口尺度。输入倍率修正已在 FRF 中，未重复缩放。",
        "",
        "令 P=Ty G Tu=[A,b;c,d]，逐频计算 Fu=A⁻¹b、Fy=cA⁻¹，"
        "则 [I,0;−Fy,1] P [I,−Fu;0,1]=diag(A,S)，S=d−cA⁻¹b。"
        "这里使用复数矩阵乘法，不把转置误写为共轭转置。",
        "",
        "**这是理想消交叉诊断，不是可实施滤波器，也不是可部署性能保证。**"
        "逐频逆尚未检查因果性、稳定性、非最小相位零点及频带外行为。"
        "加共同纯延迟会改变相位，不改变本表的端口增益或相对 NN 代价。",
        "",
        "## 25–80 Hz：理想完全消交叉需要的真实端口增益",
        "",
        "| 来源 | max‖nu−Pu Fu‖ | 所在 Hz | max‖nyᵀ−Fy Py‖ | 所在 Hz | A 条件数中位/最大 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        if r["fmin_Hz"] != 25:
            continue
        lines.append(
            f'| {r["source"]} | {r["ideal_actual_N_input_gain_max"]:.6g} | '
            f'{r["input_gain_max_frequency_Hz"]:.3g} | {r["ideal_actual_N_output_gain_max"]:.6g} | '
            f'{r["output_gain_max_frequency_Hz"]:.3g} | {r["A_condition_median"]:.6g} / '
            f'{r["A_condition_max"]:.6g} |'
        )
    lines += [
        "",
        "静态 N 端口范数为 1；本轮 FIR 约束上限为 1.05。"
        "上表的范数在物理电机/传感器方向计算，没有使用仅逻辑滤波器的系数范数替代。"
        "这是指令与单位独立白噪声下的增益代理，不是已测到的真实噪声谱。",
        "",
        "## 柔性目标频带：完全消交叉会改变多少 NN 响应",
        "",
        "复误差为 sqrt(∫|S−d|²df / ∫|d|²df)，RMS 比为 sqrt(∫|S|²df / ∫|d|²df)。"
        "它们是完整观测响应的变化，不是孤立模态的振型保留率。",
        "",
        "| 来源 | Hz | ‖S−d‖/‖d‖ | ‖S‖/‖d‖ | A 条件数中位/最大 |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        if r["fmin_Hz"] not in (260, 320, 450):
            continue
        lines.append(
            f'| {r["source"]} | {r["fmin_Hz"]:.0f}–{r["fmax_Hz"]:.0f} | '
            f'{r["Schur_relative_complex_change"]:.6g} | '
            f'{r["Schur_RMS_ratio_to_static_NN"]:.6g} | '
            f'{r["A_condition_median"]:.6g} / {r["A_condition_max"]:.6g} |'
        )
    lines += [
        "",
        "## 如何解读",
        "",
        "- 在当前固定 A、固定 NN 参考和双三角补偿结构下，"
        "若某频带的复变化超过 0.10，则该频带逐频完全消除 b、c 与"
        "“保持原 d 的复频响至 10% 内”不能同时成立；"
        "允许保留部分交叉、限定低频消除或重新定义目标，是不同的选择。",
        "- 该结论不表示仅要求 25–80 Hz 解耦的 FIR 一定无法保留高频 d。"
        "受约束设计可以在高频关闭修正，以牺牲高频消交叉换取目标保留。",
        "- 非参数与 SRIVC 的逆需求差异说明两类模型的弱块/"
        "相消组合可能不同；不能直接选看起来更好的模型数字来替代实际训练和检查验证。",
        "- A 条件数与最小奇异值都按当前单位端口坐标给出，"
        "它们是被逆块的数值诊断，不是整体系统可辨识性的充分判据。",
        "- CSV 同时给出全输入/输出操作器相对静态的谱范数增益，"
        "以及逐频代数消元的残差检查。",
        "",
        f"固定 Ty/Tu 条件数分别为 {np.linalg.cond(Ty):.8g} / {np.linalg.cond(Tu):.8g}。",
        "逐频矩阵乘法独立核对已通过：两向剩余交叉和 Schur 公式的归一化误差均小于 1e−8。",
        "",
        "复现：在本目录运行 python diagnostics.py。输出仅 diagnostic_limits.csv 和 diagnostic_limits.md。",
    ]
    (OUT / "diagnostic_limits.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for r in rows:
        print(r["source"], f'{r["fmin_Hz"]:g}–{r["fmax_Hz"]:g}',
              f'gain_u/y={r["ideal_actual_N_input_gain_max"]:.6g}/{r["ideal_actual_N_output_gain_max"]:.6g}',
              f'Sdiff={r["Schur_relative_complex_change"]:.6g}',
              f'Sratio={r["Schur_RMS_ratio_to_static_NN"]:.6g}',
              f'cond_max={r["A_condition_max"]:.6g}', flush=True)


if __name__ == "__main__":
    main()
