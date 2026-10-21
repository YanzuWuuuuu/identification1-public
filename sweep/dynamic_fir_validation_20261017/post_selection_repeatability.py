"""Post-selection diagnostics; Hcheck is inaccessible before frozen selection.

No fitting, candidate modification, or selection occurs in this program.
The selected JSON must exist and parse before any input MAT is loaded.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.io import loadmat

OUT = Path(__file__).resolve().parent
BANDS = [(25., 80.), (80., 170.), (260., 380.), (450., 800.), (25., 800.)]


def physical(H):
    return H.reshape(len(H), 4, 4).transpose(0, 2, 1)


def project(H, Ty, Tu):
    return np.einsum("ab,kbc,cd->kad", Ty, physical(H), Tu)


def blocks(H, Ty, Tu):
    P = project(H, Ty, Tu)
    return {"physical": physical(H), "A": P[:, :3, :3],
            "b_N_to_logical": P[:, :3, 3], "c_logical_to_N": P[:, 3, :3],
            "d_N_to_N": P[:, 3, 3]}


def energy(H, f):
    return float(np.trapezoid(np.sum(np.abs(H.reshape(len(f), -1)) ** 2, axis=1), f))


def main():
    selected = OUT / "selected_training_only.json"
    # This must precede loadmat of the file containing Hcheck.
    assert selected.is_file(), "Selection not frozen: Hcheck will not be loaded."
    frozen_bytes = selected.read_bytes()
    frozen = json.loads(frozen_bytes)
    assert isinstance(frozen, dict) and frozen, "Frozen selection is empty."
    selection_sha256 = hashlib.sha256(frozen_bytes).hexdigest()

    co = loadmat(OUT / "coordinates.mat", squeeze_me=True,
                 variable_names=["Py", "Pu", "ny", "nu"])
    raw = loadmat(OUT / "input_z0.mat", squeeze_me=True,
                  variable_names=["f", "Htrain", "Hcheck", "Hmodel"])
    f = np.asarray(raw["f"], dtype=float)
    Ty = np.vstack([co["Py"], co["ny"]])
    Tu = np.column_stack([co["Pu"], co["nu"]])
    data = {name: blocks(raw[key], Ty, Tu) for name, key in
            [("training", "Htrain"), ("check", "Hcheck"), ("model", "Hmodel")]}
    rows = []
    for lo, hi in BANDS:
        mask = (f >= lo - 1e-9) & (f <= hi + 1e-9)
        ff = f[mask]
        for kind, comparison in [("half_record_FRF_difference", "training"),
                                 ("SRIVC_vs_check_FRF", "model")]:
            for block, ref in data["check"].items():
                den = energy(ref[mask], ff)
                err = energy(data[comparison][block][mask] - ref[mask], ff)
                ratio = np.sqrt(err / den)
                rows.append(dict(diagnostic=kind, fmin_Hz=lo, fmax_Hz=hi,
                    block=block, ratio=float(ratio), percent=float(100 * ratio),
                    ratio_dB=float(20 * np.log10(max(ratio, 1e-300))),
                    reference_check_complex_RMS=float(np.sqrt(den / (ff[-1]-ff[0]))),
                    difference_complex_RMS=float(np.sqrt(err / (ff[-1]-ff[0])))))

    # Do not publish results if selection changed during this read-only audit.
    assert selected.read_bytes() == frozen_bytes, "Selection changed during diagnostic."
    with (OUT / "repeatability_diagnostic.csv").open("w", encoding="utf-8-sig", newline="") as h:
        wr = csv.DictWriter(h, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)

    lines = [
        "# 冻结选择后的分半 FRF 差异与 SRIVC 检查误差",
        "",
        "日期：2026-10-17。本程序先确认 selected_training_only.json 已存在且可解析，"
        "才读取检查 FRF。未修改候选、拟合、筛选或补偿器。本报告仅用于冻结后的误差诊断。",
        "",
        f"冻结选择 SHA256：{selection_sha256}",
        "",
        "坐标沿用同一训练得到的强静态 Ty/Tu。所有比值都用检查 FRF 的对应块能量作分母，"
        "对实际频率网格作梯形积分；物理4×4不做逻辑变换。",
        "",
        "分半差异 = sqrt(∫||Htrain−Hcheck||²df / ∫||Hcheck||²df)。"
        "模型检查误差 = sqrt(∫||Hmodel−Hcheck||²df / ∫||Hcheck||²df)。",
        "",
        "## 25–80 Hz 主要频段",
        "",
        "| 块 | 两半 FRF 差异 | SRIVC 对检查 FRF 差异 | 模型差异/分半差异 |",
        "|---|---:|---:|---:|",
    ]
    for block in data["check"]:
        gap = next(r for r in rows if r["fmin_Hz"] == 25 and r["fmax_Hz"] == 80
                   and r["block"] == block and r["diagnostic"] == "half_record_FRF_difference")
        model = next(r for r in rows if r["fmin_Hz"] == 25 and r["fmax_Hz"] == 80
                     and r["block"] == block and r["diagnostic"] == "SRIVC_vs_check_FRF")
        lines.append(f'| {block} | {gap["percent"]:.4f}% | {model["percent"]:.4f}% | '
                     f'{model["ratio"]/gap["ratio"]:.4g} |')
    lines += [
        "",
        "## 全部频段",
        "",
        "| Hz | 对比 | 物理4×4 | A | b：N→逻辑 | c：逻辑→N | d：N→N |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for lo, hi in BANDS:
        for kind, label in [("half_record_FRF_difference", "两半差异"),
                            ("SRIVC_vs_check_FRF", "模型与检查")]:
            rr = [r for r in rows if r["fmin_Hz"] == lo and r["fmax_Hz"] == hi
                  and r["diagnostic"] == kind]
            lines.append(f"| {lo:g}–{hi:g} | {label} | " +
                         " | ".join(f'{r["percent"]:.4f}%' for r in rr) + " |")
    lines += [
        "",
        "## 结论边界",
        "",
        "- 两半来自同一记录；后半此前已经被开发查看，因此这里不是独立重测或盲测。"
        "两半差异可以显示记录内的不稳定程度，不能直接当成置信区间、无噪真值误差或精确噪声下限。",
        "- 弱块由多个强物理响应相消而来。它的相对分半差异可能远大于物理整体，"
        "说明用于验证的小响应本身也比总体 FRF 难确定。",
        "- 如果弱块模型与检查的差异远大于分半差异，"
        "则不能仅凭“弱块噪声大”把全部模型失配解释掉；"
        "但该倍数不是统计显著性检验，也不能排除两半共享的系统性偏差。",
        "- 这里固定了数据入口、坐标和频带；"
        "没有用本诊断重新选择模型结构、FIR 阶数、延迟或正则权。",
        "",
        "复现：选择冻结后运行 python post_selection_repeatability.py。"
        "若冻结 JSON 不存在，程序在读取检查 FRF 之前终止。",
    ]
    (OUT / "repeatability_diagnostic.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for r in rows:
        if r["fmin_Hz"] == 25 and r["fmax_Hz"] == 80:
            print(r["diagnostic"], r["block"], f'{r["percent"]:.7f}%', flush=True)


if __name__ == "__main__":
    main()
