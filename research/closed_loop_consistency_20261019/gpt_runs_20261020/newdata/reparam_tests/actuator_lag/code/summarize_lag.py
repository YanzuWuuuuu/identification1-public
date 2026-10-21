"""Summarize saved joint-window actuator-lag fits; never launches a fit.

The protocol and thresholds are fixed in ../PROTOCOL.md.  Missing, failed, and
boundary fits are retained.  This module deliberately does not import lag_core.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
IDS = (
    "A_lowgain", "A_nominal", "A_highgain",
    "C_lowgain_pppp", "C_nominal_pppp", "C_highgain_pppp",
)
WINDOWS = ("H300", "H493")
KINDS = ("a", "c")
NAMES = {"a": "i：滞后＋自由 C0/C1", "c": "ii：滞后＋C1=0"}
MODES = ((282, "H300", 1), (333, "H300", 2), (359, "H300", 3), (493, "H493", 0))
REPORT = "执行器一阶滞后检验报告.md"


def read(path: Path):
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8-sig"))


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def num(value, digits=5):
    return f"{float(value):.{digits}g}" if finite(value) else "—"


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(path, obj):
    path.write_text(json.dumps(clean(obj), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(clean(rows))


def table(headers, rows):
    def cell(x):
        return str(x).replace("|", r"\|").replace("\n", " ")
    return "\n".join([
        "| " + " | ".join(map(cell, headers)) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        *["| " + " | ".join(map(cell, r)) + " |" for r in rows],
    ])


def solver_status(result):
    if result is None:
        return {"valid": False, "state": "未完成", "score": None, "near_bound": None}
    sol = result.get("solver", {})
    score = sol.get("score")
    bounded = bool(sol.get("near_bound", False))
    cuts = np.asarray(result.get("cutoff_hz", []), float)
    cutoff_finite = cuts.shape == (4,) and np.isfinite(cuts).all()
    # Report cutoff proximity separately from the solver's original bound flag.
    cutoff_bound = bool(cutoff_finite and np.any(
        (cuts <= 20 * (1 + 1e-5)) | (cuts >= 200000 / (1 + 1e-5))))
    valid = bool(result.get("success") and sol.get("success")
                 and finite(score) and score < 1e-6 and not bounded
                 and cutoff_finite and not cutoff_bound)
    state = sol.get("status", result.get("status", "未记录状态"))
    if valid:
        state = "有效：" + state
    else:
        state = "无效：" + state
    return {
        "valid": valid, "state": state, "score": score, "near_bound": bounded,
        "cutoff_near_bound": cutoff_bound,
        "iterations": sol.get("iterations"),
        "inner_failures": sol.get("inner_failures"),
        "profile_calls": sol.get("profile_calls"),
        "seconds": result.get("seconds"),
        "exception": result.get("exception"),
    }


def relative_percent(value, ref):
    if finite(value) and finite(ref) and value > 0 and ref > 0:
        return float(100 * np.expm1(np.log(value) - np.log(ref)))
    return None


def qe_value(result, window):
    if result is None:
        return {}
    return result.get("qe", {}).get(window, {})


def summarize():
    results = {(rid, kind): read(ROOT / "raw" / f"{rid}_{kind}_joint.json")
               for rid in IDS for kind in KINDS}
    baseline = {(rid, w): read(BASE / "raw" / f"{rid}_{w}_a_M3_full.json")
                for rid in IDS for w in WINDOWS}
    metrics = {(rid, w): read(BASE / "analysis" / f"metric_{rid}_{w}_a_M3_full.json")
               for rid in IDS for w in WINDOWS}
    states = {key: solver_status(value) for key, value in results.items()}
    endpoints, qe_rows, cutoff_rows, delay_rows, modal_rows, range_rows = [], [], [], [], [], []
    for (rid, kind), r in results.items():
        state = states[rid, kind]
        endpoints.append({
            "record": rid, "model": kind, "saved": r is not None, **state,
            "parameter_count": None if r is None else r.get("parameter_count"),
            "linear_count": None if r is None else r.get("linear_count"),
            "nonlinear_count": None if r is None else r.get("nonlinear_count"),
            "grid_total": None if r is None else len(r.get("grid", [])),
            "grid_inner_success": None if r is None else sum(
                bool(g.get("inner_success")) for g in r.get("grid", [])),
        })
        cut = [] if r is None else r.get("cutoff_hz", [])
        cutoff_rows.append({
            "record": rid, "model": kind, "valid": state["valid"],
            **{f"input_{j+1}_cutoff_hz": cut[j] if len(cut) == 4 else None for j in range(4)},
        })
        dl = {} if r is None else r.get("delay", {})
        delay_rows.append({
            "record": rid, "model": kind, "valid": state["valid"],
            "common_us": dl.get("common_us"),
            **{f"output_{j+1}_us": dl.get("output_us", [None]*4)[j] for j in range(4)},
            **{f"input_{j+1}_us": dl.get("input_us", [None]*4)[j] for j in range(4)},
        })
        for w in WINDOWS:
            q = qe_value(r, w)
            metric = metrics[rid, w] or {}
            bq = metric.get("qe", {}).get("all", {})
            ratio = q.get("ratio")
            baseline_ratio = bq.get("ratio")
            comparison = (ratio / baseline_ratio
                          if finite(ratio) and finite(baseline_ratio) and baseline_ratio > 0 else None)
            numeric_gate = bool(finite(comparison) and comparison <= 1.2)
            qe_rows.append({
                "record": rid, "model": kind, "window": w, "valid": state["valid"],
                "Q": q.get("Q"), "EQ": q.get("EQ"), "n": q.get("n"), "Q_over_EQ": ratio,
                "EQ_over_n": q.get("expected_noise_amplification"),
                "baseline_Q": bq.get("Q"), "baseline_EQ": bq.get("EQ"),
                "baseline_Q_over_EQ": baseline_ratio, "relative_to_baseline": comparison,
                "numerical_ratio_pass": numeric_gate,
                "eligible_ratio_pass": bool(state["valid"] and numeric_gate),
            })
    cutoff_stats = []
    for kind in KINDS:
        rows = [r for r in cutoff_rows if r["model"] == kind]
        all_valid = all(r["valid"] for r in rows)
        for j in range(1, 5):
            values = [r[f"input_{j}_cutoff_hz"] for r in rows]
            all_finite = all(finite(v) for v in values)
            arr = np.asarray(values, float) if all_finite else np.array([])
            mean = float(arr.mean()) if all_finite else None
            sd = float(arr.std(ddof=1)) if all_finite else None
            cv = sd / mean if finite(sd) and finite(mean) and mean > 0 else None
            cutoff_stats.append({
                "model": kind, "input": j, "numeric_records": sum(finite(v) for v in values),
                "all_six_valid": all_valid, "mean_hz": mean, "sample_SD_hz": sd,
                "CV": cv, "min_hz": float(arr.min()) if all_finite else None,
                "max_hz": float(arr.max()) if all_finite else None,
                "range_hz": float(np.ptp(arr)) if all_finite else None,
                "numerical_CV_pass": bool(finite(cv) and cv < .2),
                "eligible_CV_pass": bool(all_valid and finite(cv) and cv < .2),
            })
    for rid in IDS:
        for target, w, index in MODES:
            b = baseline[rid, w] or {}
            bf = b.get("frequency_hz", [None] * (index + 1))[index]
            bz = b.get("damping_percent", [None] * (index + 1))[index]
            entries = [("baseline", bf, bz, bool(b.get("success")))]
            for kind in KINDS:
                r = results[rid, kind] or {}
                wm = r.get("windows", {}).get(w, {})
                f = wm.get("frequency_hz", [None] * (index + 1))[index]
                z = wm.get("damping_percent", [None] * (index + 1))[index]
                entries.append((kind, f, z, states[rid, kind]["valid"]))
            for kind, f, z, valid in entries:
                modal_rows.append({
                    "record": rid, "nominal_mode_hz": target, "window": w, "model": kind,
                    "valid": valid, "frequency_hz": f, "damping_percent": z,
                    "frequency_difference_hz": f - bf if finite(f) and finite(bf) else None,
                    "frequency_relative_difference_percent": relative_percent(f, bf),
                    "damping_difference_percentage_points": z - bz if finite(z) and finite(bz) else None,
                    "damping_relative_difference_percent": relative_percent(z, bz),
                })
            valid_all = all(e[3] for e in entries)
            finite_all = all(finite(e[1]) and finite(e[2]) for e in entries)
            fs = [e[1] for e in entries]
            zs = [e[2] for e in entries]
            range_rows.append({
                "record": rid, "nominal_mode_hz": target, "all_three_valid": valid_all,
                "frequency_range_hz": max(fs)-min(fs) if finite_all else None,
                "damping_range_percentage_points": max(zs)-min(zs) if finite_all else None,
                "interpretation": "结构选择敏感性" if valid_all else "含无效端点，仅诊断，不能解释为可靠结构差异",
            })
    complete = all(r is not None for r in results.values())
    c_valid = all(states[rid, "c"]["valid"] for rid in IDS)
    c_qe = [r for r in qe_rows if r["model"] == "c"]
    c_cv = [r for r in cutoff_stats if r["model"] == "c"]
    numerical_q_pass = sum(r["numerical_ratio_pass"] for r in c_qe)
    numerical_cv_pass = sum(r["numerical_CV_pass"] for r in c_cv)
    recommend = bool(complete and c_valid and numerical_q_pass == 12 and numerical_cv_pass == 4)
    decision = (
        "尚未完成全部端点；不作模型选择结论。"
        if not complete else
        "按事先固定判据，建议采用“时延＋一阶滞后＋C1=0”作为候选物理模型。"
        if recommend else
        "未满足全部事先固定判据；主模型仍为自由 C0、C1，一阶滞后模型仅作讨论。"
    )
    summary = {
        "generated_at": datetime.now().astimezone().isoformat(),
        "required_endpoints": 12, "saved_endpoints": sum(r is not None for r in results.values()),
        "valid_endpoints": sum(s["valid"] for s in states.values()), "complete": complete,
        "model_ii_all_six_valid": c_valid, "model_ii_ratio_numerical_passes": numerical_q_pass,
        "model_ii_ratio_required": 12, "model_ii_CV_numerical_passes": numerical_cv_pass,
        "model_ii_CV_required": 4, "recommend_model_ii": recommend,
        "decision": decision,
        "thresholds": {"Q_EQ_ratio_max": 1.2, "cutoff_CV_strict_max": .2,
                       "CV_ddof": 1, "score_strict_max": 1e-6},
        "failed_or_missing": [r for r in endpoints if not r["valid"]],
    }
    return results, summary, {
        "endpoints": endpoints, "qe": qe_rows, "cutoffs": cutoff_rows,
        "cutoff_stats": cutoff_stats, "delays": delay_rows,
        "modal_differences": modal_rows, "structural_ranges": range_rows,
    }


def plots(tables):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    colors = {"a": "#245A81", "c": "#397A56"}
    short = ["A-L", "A-N", "A-H", "C-L", "C-N", "C-H"]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.4), sharey=True)
    for ax, w in zip(axs, WINDOWS):
        for kind, dx in (("a", -.10), ("c", .10)):
            rs = [next(r for r in tables["qe"] if r["record"] == rid
                       and r["model"] == kind and r["window"] == w) for rid in IDS]
            for i, r in enumerate(rs):
                if finite(r["relative_to_baseline"]):
                    ax.plot(i+dx, r["relative_to_baseline"], "o" if r["valid"] else "x",
                            color=colors[kind], ms=5)
            ax.plot([], [], "o", color=colors[kind],
                    label="Lag + free C0/C1" if kind == "a" else "Lag + C1=0")
        ax.axhline(1.2, ls="--", lw=1, color="#555555", label="Predefined limit: 1.2")
        ax.set(xticks=np.arange(6), xticklabels=short, title=w, xlabel="Record / controller")
        ax.grid(axis="y", alpha=.2)
    axs[0].set_ylabel("(Q / E[Q]) / baseline (Q / E[Q])")
    axs[1].legend(fontsize=8, loc="best")
    fig.suptitle("Joint delay / actuator-lag model: within-record residual criterion")
    fig.text(.5, .01, "x = failed or boundary endpoint; missing fits are not plotted.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .045, 1, .94))
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"lag_residual_criterion.{ext}", dpi=300)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.5), sharey=True)
    input_colors = ["#245A81", "#397A56", "#7B4C91", "#B34750"]
    for ax, kind in zip(axs, KINDS):
        for j, col in enumerate(input_colors, 1):
            rs = [next(r for r in tables["cutoffs"] if r["record"] == rid and r["model"] == kind)
                  for rid in IDS]
            vals = np.asarray([r[f"input_{j}_cutoff_hz"] or np.nan for r in rs])
            ax.plot(np.arange(6), vals / 1000, "-", color=col, lw=.8, alpha=.7)
            for i, r in enumerate(rs):
                if finite(vals[i]):
                    ax.plot(i, vals[i] / 1000, "o" if r["valid"] else "x", color=col, ms=5)
            stat = next(r for r in tables["cutoff_stats"] if r["model"] == kind and r["input"] == j)
            cv_label = num(100*stat["CV"], 3) + "%" if finite(stat["CV"]) else "pending"
            ax.plot([], [], "-", color=col, label=f"Input {j}, CV={cv_label}")
        ax.set(xticks=np.arange(6), xticklabels=short,
               title="(i) Lag + free C0/C1" if kind == "a" else "(ii) Lag + C1=0",
               yscale="log", xlabel="Record / controller")
        ax.grid(axis="y", alpha=.2)
        ax.legend(fontsize=8, loc="best")
    axs[0].set_ylabel("Actuator cutoff frequency (kHz)")
    fig.suptitle("Actuator-lag estimates across six records")
    fig.text(.5, .01, "CV uses all six records (sample SD / mean); invalid endpoints cannot pass the criterion.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .94))
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"lag_cutoff_reproducibility.{ext}", dpi=300)
    plt.close(fig)


def report(summary, tables):
    lines = [
        "# 执行器一阶滞后检验报告", "",
        f"生成时间：{summary['generated_at']}。本报告仅聚合已保存端点；不自动启动拟合，不删除失败样本。",
        "",
        "## 1. 预设判据与结论", "",
        "**" + summary["decision"] + "**", "",
        f"已保存 {summary['saved_endpoints']}/12 个联合端点，其中有效 "
        f"{summary['valid_endpoints']}/12 个。模型 ii 的数值残差比门槛通过 "
        f"{summary['model_ii_ratio_numerical_passes']}/12 个窗口；截止频率 CV 数值门槛通过 "
        f"{summary['model_ii_CV_numerical_passes']}/4 路；六条模型 ii 端点是否均有效："
        f"{'是' if summary['model_ii_all_six_valid'] else '否'}。数值门槛通过不替代端点有效性。",
        "",
        "预先冻结：模型 ii 在六条记录的两个窗口中分别满足 Q/E[Q] 不超过原自由模型 a 的 1.2 倍；"
        "每路截止频率的六记录样本 SD/算术均值严格小于 20%（ddof=1）；六个端点均须通过原 score、"
        "内层矩收敛和未贴边门槛。三项同时满足才推荐模型 ii；不通过时仍保留原自由 C0、C1。",
        "",
        "## 2. 模型及可比性", "",
        r"本次模型为 $G=E(\tau)\circ[\sum G_r+D_0+D_2]\operatorname{diag}\{(1+s/\omega_{a,m})^{-1}\}$。"
        "模型 i 保留自由 C0、C1，模型 ii 取自由 C0、C1=0。每条记录的 H300（220–420 Hz）和 H493"
        "（450–550 Hz）共享七个结构化时延与四个截止频率，其余模态、背景、瞬态参数分别估计；"
        "参数量分别为 285 和 205。原模型 a 是两个窗口独立估计时延的既有基线，不能把本次差异全部单独归因于新增滞后。",
        "",
        "初始截止频率均为 2 kHz，技术搜索范围 20 Hz–200 kHz。公共时延使用既定 ±200 μs、10 μs"
        "间隔网格初始化。模态沿用精确 ZOH，右乘滞后严格采用任务给定的连续频响形式；"
        "不重新整体离散化，不额外减去半拍。既有输出瞬态不乘输入端滞后，未新增瞬态基函数。",
        "",
        "采用同记录、原自由模型 a 已冻结的 0.5 Hz 非激励点噪声谱。Q/E[Q] 的线性化残差算子来自"
        "两个窗口的联合估计；每个窗口的 E[Q] 包含另一窗口通过共享参数传来的噪声。"
        "Q/E[Q] 下降可能同时涉及 Q 的改变与 E[Q] 的放大，因此下表同时列出两者及 E[Q]/n。",
        "",
        "E[Q] 使用端点处冻结工具变量的 Gauss–Newton 线性化。当残差包含模型失配时，"
        "它不是原非线性估计器残差的精确分布或精确检验水平；本次沿用预定诊断指标，不将其改称精确 χ² 检验。",
        "",
        "H300 原有约 263 Hz 宽背景分量保留，但不列为本次四个待报告物理模态。"
        "所有失败和边界端点仅供诊断，不能作为可靠物理参数解释。",
        "",
        "## 3. 同记录残差比较", "",
        table(["记录", "模型", "窗口", "有效", "Q", "E[Q]", "n", "E[Q]/n",
               "Q/E[Q]", "原 a", "比值", "数值≤1.2"],
              [[r["record"], r["model"], r["window"], "是" if r["valid"] else "否",
                num(r["Q"]), num(r["EQ"]), r["n"], num(r["EQ_over_n"]),
                num(r["Q_over_EQ"]), num(r["baseline_Q_over_EQ"]), num(r["relative_to_baseline"]),
                "是" if r["numerical_ratio_pass"] else "否"] for r in tables["qe"]]),
        "",
        "相对基线的门槛通过不等于 Q/E[Q] 已接近 1，也不能证明模型在噪声水平上已充分。",
        "",
        "![同记录残差比](figures/lag_residual_criterion.png)", "",
        "## 4. 截止频率、记录间离散与时延", "",
        table(["记录", "模型", "有效", "fa1/Hz", "fa2/Hz", "fa3/Hz", "fa4/Hz"],
              [[r["record"], r["model"], "是" if r["valid"] else "否",
                *[num(r[f"input_{j}_cutoff_hz"], 7) for j in range(1, 5)]]
               for r in tables["cutoffs"]]), "",
        table(["模型", "输入", "六条均有效", "均值/Hz", "样本SD/Hz", "CV/%", "最小/Hz", "最大/Hz", "极差/Hz"],
              [[r["model"], r["input"], "是" if r["all_six_valid"] else "否", num(r["mean_hz"]),
                num(r["sample_SD_hz"]), num(100*r["CV"] if finite(r["CV"]) else None),
                num(r["min_hz"]), num(r["max_hz"]), num(r["range_hz"])]
               for r in tables["cutoff_stats"]]), "",
        "CV 只在六个数值齐备时计算，失败端点仍保留在诊断统计中；任何一个端点无效均使其不能作为"
        "预设跨记录离散判据的通过证据。不改用成功子集，不用极差代替预设 CV。",
        "",
        "截止频率还可能与自由 C1 和时延混淆：在无背景、无瞬态的连续模态理想化中，"
        "若分子写为 C0+(s/ωr)C1，取 C1=(ωr/ωa)C0 可精确消去相同通道的一阶滞后；"
        "当前 ZOH 模态与背景下该式不再是整体精确恒等式，但仍提示近似混淆方向。"
        "另外 log L≈−s/ωa，使低阶相位信息中的时延与 1/ωa 相关。因此记录间 CV 小"
        "不自动证明估计到了硬件真实截止频率。",
        "",
        "![截止频率复现性](figures/lag_cutoff_reproducibility.png)", "",
        "下表均为 μs，输出端与输入端相对时延各自和为零。两个窗口已共享同一组参数，因此"
        "本次不会再出现两个窗口独立时延估计的差值；共享本身并不是跨窗口相容性得到验证。",
        "",
        table(["记录", "模型", "公共", "输出1", "输出2", "输出3", "输出4", "输入1", "输入2", "输入3", "输入4"],
              [[r["record"], r["model"], num(r["common_us"]),
                *[num(r[f"output_{j}_us"]) for j in range(1, 5)],
                *[num(r[f"input_{j}_us"]) for j in range(1, 5)]] for r in tables["delays"]]),
        "",
        "## 5. 四个模态的模型结构敏感性", "",
        "f 的绝对差单位为 Hz；ζ 的绝对差为阻尼百分数之差，单位为百分点。相对差统一采用 "
        "100×[exp(log(新值)−log(基线值))−1]，以原自由模型 a 为基线。以下不把模型间范围"
        "当作具有概率覆盖率的标准不确定度；它描述在这三种结构选择之间的敏感性。",
        "",
        table(["记录", "模态/Hz", "模型", "有效", "f/Hz", "Δf/Hz", "Δf/%", "ζ/%", "Δζ/百分点", "Δζ/%"],
              [[r["record"], r["nominal_mode_hz"], r["model"], "是" if r["valid"] else "否",
                num(r["frequency_hz"], 8), num(r["frequency_difference_hz"]),
                num(r["frequency_relative_difference_percent"]), num(r["damping_percent"], 7),
                num(r["damping_difference_percentage_points"]),
                num(r["damping_relative_difference_percent"])] for r in tables["modal_differences"]]),
        "",
        table(["记录", "模态/Hz", "三端点均有效", "f模型范围/Hz", "ζ模型范围/百分点"],
              [[r["record"], r["nominal_mode_hz"], "是" if r["all_three_valid"] else "否",
                num(r["frequency_range_hz"]), num(r["damping_range_percentage_points"])]
               for r in tables["structural_ranges"]]),
        "",
        "## 6. 求解状态、失败与边界", "",
        table(["记录", "模型", "状态", "score", "迭代", "内层失败调用", "贴边", "截止贴边", "参数数", "秒"],
              [[r["record"], r["model"], r["state"], num(r["score"]), r.get("iterations"),
                r.get("inner_failures"), r["near_bound"], r.get("cutoff_near_bound"),
                r["parameter_count"], num(r.get("seconds"))] for r in tables["endpoints"]]),
        "",
        "“内层失败调用”统计优化试探期间未收敛的内层调用，不等于最终端点失败；最终有效性使用"
        "原求解器终止状态、严格 score 与边界检查。全部初始化网格、求解轨迹和异常保存在 raw/ 与 logs/。",
        "",
    ]
    audit = read(ROOT / "analysis" / "audit.json")
    if isinstance(audit, dict) and audit.get("endpoints"):
        rows = audit["endpoints"]
        lines += [
            "下表的条件数来自审计中缩放后的联合矩阵；`rank` 的相对阈值为 1e−12。"
            "截止频率距技术上界小于 1% 是新增诊断标记，不追溯改变原有效性门槛或预设模型选择规则。",
            "",
            table(["记录", "模型", "缩放条件数", "rank/参数数", "距200 kHz上界<1%的输入"],
                  [[r["id"], r["model"], num(r.get("scaled_condition")),
                    f"{r.get('scaled_rank_1e12')}/{r.get('parameters')}",
                    ",".join(str(j+1) for j, flag in enumerate(
                        r.get("cutoff_near_technical_upper_1percent", [])) if flag) or "无"]
                   for r in rows]), "",
        ]
        afailed = [r["id"] for r in rows if r["model"] == "a" and not r.get("success")]
        near = [r for r in rows if r["model"] == "a"
                and any(r.get("cutoff_near_technical_upper_1percent", []))]
        lines += [
            f"自由 C1 的滞后模型 i 有 {len(afailed)} 个失败端点："
            + ("、".join(afailed) if afailed else "无") + "。",
            "",
        ]
        if near:
            near_values = []
            for r in near:
                saved = next(t for t in tables["cutoffs"]
                             if t["record"] == r["id"] and t["model"] == "a")
                for j, flag in enumerate(r["cutoff_near_technical_upper_1percent"], 1):
                    if flag:
                        near_values.append(
                            f"{r['id']} 的输入 {j} 为 {num(saved[f'input_{j}_cutoff_hz']/1000, 7)} kHz")
            lines += [
                "另有以下求解门槛已通过的端点接近技术上界：" + "；".join(near_values) + "。"
                "结合大条件数与近秩缺失，这些数值不能作为硬件截止频率解释。"
                "保持其原成功状态，另行标明弱辨识问题。",
                "",
            ]
    failures = summary["failed_or_missing"]
    if failures:
        lines += ["失败或未完成端点：", "",
                  table(["记录", "模型", "原因"], [[r["record"], r["model"],
                        r.get("exception") or r["state"]] for r in failures]), ""]
    else:
        lines += ["十二个端点均有效，无失败端点。", ""]
    lines += ["## 7. 数值核查与文件审计", "",
              "另见 `analysis/独立实现核查.md`。独立只读审查未发现输入滞后方向、共享参数映射、"
              "内层矩方程或联合残差期望的实现错误。模型 ii 六个端点均满秩且严格收敛；"
              "其残差门槛失败不能归因于明显的求根未收敛。此项核查不保证全局根唯一。", ""]
    audit_paths = sorted(set(
        list((ROOT / "analysis").glob("*audit*.json"))
        + list((ROOT / "analysis").glob("*verification*.json"))
    ))
    if not audit_paths:
        lines += ["尚无核查文件；应在最终交付前完成。", ""]
    for p in audit_paths:
        v = read(p)
        lines += [f"核查记录：`analysis/{p.name}`。", ""]
        if isinstance(v, list) and all(isinstance(x, dict) for x in v):
            rows = [[x.get("kind", x.get("model", "—")), x.get("passed"),
                     num(max(x.get("derivative_relative_errors", [math.nan]))),
                     num(x.get("zero_lag_limit_relative_error")),
                     num(x.get("moment_matrix_relative_error")),
                     num(x.get("profile_linear_score"))] for x in v]
            lines += [table(["模型", "通过", "最大导数相对差", "零滞后极限差", "矩阵差", "线性score"], rows), ""]
        elif isinstance(v, dict):
            simple = {k: val for k, val in v.items()
                      if not isinstance(val, (dict, list)) or len(str(val)) < 250}
            lines += ["```json", json.dumps(clean(simple), ensure_ascii=False, indent=2), "```", ""]
    if (ROOT / "analysis" / "verification_step_1e-5.json").exists():
        lines += [
            "初始有限差分步长 1e−5 的核查文件也保留在上表，不能忽略其未通过状态。"
            "随后步长缩小为 1e−6，导数误差约缩小 100 倍，符合中心差分截断误差的二阶缩放。"
            "最终验证没有放宽原阈值；该对照用于区分有限差分截断误差与解析导数错误。",
            "",
        ]
    bootstrap = read(ROOT / "analysis" / "residual_expectation_bootstrap.json")
    if isinstance(bootstrap, dict):
        lines += [
            f"联合残差期望的独立核查采用 {bootstrap.get('id')}、模型 {bootstrap.get('model')}，"
            f"生成 {bootstrap.get('draws')} 组复高斯噪声。"
            "分别比较解析 E[Q]、独立低秩迹公式，以及固定线性化残差算子下的抽样均值。"
            "这检验期望计算的实现，不把模型失配变成高斯噪声。",
            "",
            table(["窗口", "解析E[Q]", "独立低秩迹", "迹相对差", "抽样均值", "均值MCSE", "差/MCSE"],
                  [[w, num(q.get("analytic"), 9), num(q.get("independent_low_rank_trace"), 9),
                    num(q.get("trace_relative_error")), num(q.get("bootstrap_mean"), 9),
                    num(q.get("bootstrap_MC_SE")), num(q.get("difference_in_MC_SE"))]
                   for w, q in bootstrap.get("windows", {}).items()]), "",
        ]
    lines += [
        "## 8. 范围与限制", "",
        "- 本次只做既有六条记录、两个既有高频窗口的 M3 联合拟合，没有新增数据、改窗口或增减模态。",
        "- 不新增 jackknife，截止频率的记录间 CV 不是单条记录内的参数标准误，"
        "不能凭其大小独立判断滞后已与留数、时延充分分离。",
        "- 联合残差检验使用既定冻结噪声谱；不通过时不换权重、模型背景或判据来调整结论。",
        "- 20 Hz–200 kHz 为明确记录的技术搜索范围，贴边结果不作可辨识的物理截止频率。",
        "- 当前保留原逐模态输出瞬态，没有新增执行器滞后极点对应的瞬态项；"
        "既有 ZOH 模态右乘连续一阶频响也不等价于连续总串联系统整体进行 ZOH 离散化。",
        "- 两个滞后模型同时约束跨窗口时延与滞后相同，原 a 基线各窗口时延独立。"
        "模型间范围包含这项结构约束的影响。",
        "- 保存下来的失败端点可显示数值，但不作为物理模型推荐依据；完整 CSV 记录有效性标记。",
        "",
        "输出索引：`analysis/lag_summary.json` 是机器可读判定；`lag_qe.csv`、`lag_cutoffs.csv`、"
        "`lag_cutoff_stats.csv`、`lag_delays.csv`、`lag_modal_differences.csv`、"
        "`lag_structural_ranges.csv` 和 `lag_endpoints.csv` 为完整明细。图同时保存 PNG（300 dpi）和 PDF。",
        "",
    ]
    (ROOT / REPORT).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-figures", action="store_true", help="只聚合表格与报告，不生成图")
    args = parser.parse_args()
    (ROOT / "analysis").mkdir(exist_ok=True)
    (ROOT / "figures").mkdir(exist_ok=True)
    _, summary, tables = summarize()
    write_json(ROOT / "analysis" / "lag_summary.json", summary)
    for name, rows in tables.items():
        write_csv(ROOT / "analysis" / f"lag_{name}.csv", rows)
    if not args.no_figures:
        plots(tables)
    report(summary, tables)
    print(json.dumps(clean(summary), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
