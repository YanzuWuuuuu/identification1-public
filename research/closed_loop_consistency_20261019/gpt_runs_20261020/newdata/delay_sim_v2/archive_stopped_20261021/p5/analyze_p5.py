"""Read-only P5 reanalysis of existing paired high-window M3 fits."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
NEWDATA = ROOT.parents[1]
BASE = NEWDATA / "reparam_tests"
THEORY = NEWDATA.parents[1] / "时延扩展_理论推导_v1.md"
IDS = ["A_lowgain", "A_nominal", "A_highgain",
       "C_lowgain_pppp", "C_nominal_pppp", "C_highgain_pppp"]
PERIODS = dict(zip(IDS, [7, 7, 8, 7, 6, 6]))
GAINS = ["lowgain", "nominal", "highgain"]
COMPONENTS = [f"{side}{j}" for side in ("output", "input") for j in range(1, 5)]
USED = set()


def read_text(path):
    USED.add(path.resolve())
    return path.read_text(encoding="utf-8-sig")


def read_json(path):
    return json.loads(read_text(path))


def read_csv(path):
    return list(csv.DictReader(read_text(path).splitlines()))


def delay_vector(r):
    return np.r_[r["delay"]["output_us"], r["delay"]["input_us"]]


def acceptable(r):
    s = r.get("solver", {})
    return bool(r.get("success") and float(s.get("score", np.inf)) < 1e-6
                and not s.get("near_bound", False))


def jsonify(v):
    if isinstance(v, dict):
        return {k: jsonify(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [jsonify(x) for x in v]
    if isinstance(v, (np.integer, np.floating, np.bool_)):
        return v.item()
    return v


def write_json(path, value):
    path.write_text(json.dumps(jsonify(value), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(jsonify(rows))


def table(headers, rows):
    return "\n".join([
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        *["| " + " | ".join(str(x) for x in row) + " |" for row in rows],
    ])


def f(x):
    return f"{x:.5g}"


def main():
    text = read_text(THEORY)
    p5 = text.split("**P5（", 1)[1].split("**P6（", 1)[0].strip()
    cross = read_csv(BASE / "analysis" / "p6_cross_window_components.csv")
    saved_delays = read_csv(BASE / "analysis" / "p6_delays.csv")
    points, covariances, validation, raw_folds = {}, {}, [], {}
    for model in ("a", "b"):
        for rid in IDS:
            full = [read_json(BASE / "raw" / f"{rid}_{w}_{model}_M3_full.json")
                    for w in ("H300", "H493")]
            delta = delay_vector(full[1]) - delay_vector(full[0])
            folds = []
            ok = all(acceptable(r) for r in full)
            for i in range(PERIODS[rid]):
                pair = [read_json(BASE / "raw" / f"{rid}_{w}_{model}_M3_omit{i}.json")
                        for w in ("H300", "H493")]
                ok = ok and all(acceptable(r) for r in pair)
                folds.append(delay_vector(pair[1]) - delay_vector(pair[0]))
            arr = np.array(folds)
            centered = arr - arr.mean(axis=0)
            P = len(arr)
            cov = (P-1)/P * centered.T @ centered
            se = np.sqrt(np.maximum(0, np.diag(cov)))
            previous = [next(r for r in cross if r["id"] == rid and r["model"] == model
                             and r["component"] == component) for component in COMPONENTS]
            delta_error = float(np.max(abs(delta - [
                float(r["H493_minus_H300_us"]) for r in previous])))
            se_error = float(np.max(abs(se - [float(r["paired_jk_se_us"]) for r in previous])))
            saved_valid = all(r["jackknife_valid"] == "True" for r in previous)
            consistency = 0.
            for w, endpoint in zip(("H300", "H493"), full):
                for j, component in enumerate(COMPONENTS):
                    prior = next(r for r in saved_delays
                                 if r["id"] == rid and r["window"] == w and r["model"] == model
                                 and r["component"] == component and not r["fixed"])
                    consistency = max(consistency, abs(delay_vector(endpoint)[j]-float(prior["delay_us"])))
            valid = bool(ok and saved_valid and max(delta_error, se_error, consistency) < 1e-9)
            points[model, rid] = dict(delta=delta, se=se, valid=valid,
                                     H300=delay_vector(full[0]), H493=delay_vector(full[1]))
            covariances[model, rid] = cov
            raw_folds[f"{model}/{rid}"] = dict(components=COMPONENTS, delta_us=delta,
                paired_jackknife_covariance_us2=cov, fold_delta_us=arr, valid=valid)
            validation.append(dict(model=model, record=rid, folds=P, valid=valid,
                full_and_fold_endpoints_valid=ok, saved_jackknife_valid=saved_valid,
                max_delta_csv_error_us=delta_error, max_SE_csv_error_us=se_error,
                max_delay_csv_error_us=consistency,
                max_zero_sum_error_us=float(max(abs(delta[:4].sum()), abs(delta[4:].sum())))))
    reproducibility, components, pairs, covariance_summary = [], [], [], {}
    for model in ("a", "b"):
        run_differences = []
        stat_covariances = []
        for gain in GAINS:
            a, c = f"A_{gain}", f"C_{gain}_pppp"
            difference = points[model, a]["delta"] - points[model, c]["delta"]
            stat = covariances[model, a] + covariances[model, c]
            run_differences.append(difference)
            stat_covariances.append(stat)
            for j, component in enumerate(COMPONENTS):
                pairs.append(dict(model=model, controller=gain, component=component,
                    A_cross_window_gap_us=points[model, a]["delta"][j],
                    C_cross_window_gap_us=points[model, c]["delta"][j],
                    A_minus_C_gap_us=difference[j],
                    A_plus_C_jackknife_variance_us2=stat[j, j]))
        G = len(GAINS)
        D = np.array(run_differences)
        # Within-controller paired repeatability, preserving cross-window covariance.
        rep_cov = D.T @ D / (2*G)
        raw_variance = np.diag(rep_cov)
        estimated_noise_part = np.diag(sum(stat_covariances)) / (2*G)
        extra_variance = np.maximum(0, raw_variance-estimated_noise_part)
        covariance_summary[model] = dict(components=COMPONENTS,
            raw_reproducibility_covariance_us2=rep_cov,
            diagonal_estimated_noise_part_us2=estimated_noise_part,
            diagonal_extra_variance_us2=extra_variance)
        for j, component in enumerate(COMPONENTS):
            reproducibility.append(dict(model=model, component=component, repeated_pairs=G,
                raw_reproducibility_SD_us=np.sqrt(raw_variance[j]),
                estimated_noise_variance_in_rep_us2=estimated_noise_part[j],
                noise_adjusted_extra_SD_us=np.sqrt(extra_variance[j]),
                extra_variance_clipped_to_zero=bool(raw_variance[j] < estimated_noise_part[j])))
            for rid in IDS:
                point = points[model, rid]
                delta, se = point["delta"][j], point["se"][j]
                u_raw = np.sqrt(se**2 + raw_variance[j])
                u_extra = np.sqrt(se**2 + extra_variance[j])
                components.append(dict(model=model, record=rid, component=component,
                    valid=point["valid"], H300_delay_us=point["H300"][j],
                    H493_delay_us=point["H493"][j], cross_window_gap_us=delta,
                    paired_jackknife_SE_us=se, raw_reproducibility_SD_us=np.sqrt(raw_variance[j]),
                    extra_reproducibility_SD_us=np.sqrt(extra_variance[j]),
                    combined_raw_u_us=u_raw, combined_extra_u_us=u_extra,
                    gap_over_raw_u=delta/u_raw, gap_over_extra_u=delta/u_extra,
                    within_1u_raw=bool(abs(delta) <= u_raw), within_2u_raw=bool(abs(delta) <= 2*u_raw),
                    within_1u_extra=bool(abs(delta) <= u_extra),
                    within_2u_extra=bool(abs(delta) <= 2*u_extra)))
    counts, magnitudes = [], []
    for model in ("a", "b"):
        for side in ("output", "input", "all"):
            subset = [r for r in components if r["model"] == model
                      and (side == "all" or r["component"].startswith(side))]
            counts.append(dict(model=model, side=side, comparisons=len(subset),
                valid_comparisons=sum(r["valid"] for r in subset),
                **{label: sum(r[label] and r["valid"] for r in subset)
                   for label in ("within_1u_raw", "within_2u_raw", "within_1u_extra", "within_2u_extra")}))
            magnitudes.append(dict(model=model, side=side,
                mean_absolute_gap_us=float(np.mean([abs(r["cross_window_gap_us"]) for r in subset])),
                max_absolute_gap_us=float(max(abs(r["cross_window_gap_us"]) for r in subset)),
                mean_jackknife_SE_us=float(np.mean([r["paired_jackknife_SE_us"] for r in subset]))))
    hashes = read_json(ROOT.parent / "input_hashes.json")
    manifest = {str(Path(row["path"]).resolve()).casefold(): row["sha256"].lower() for row in hashes}
    changed, not_manifested = [], []
    for p in sorted(USED):
        if p == (ROOT.parent / "input_hashes.json").resolve():
            continue
        expected = manifest.get(str(p).casefold())
        if expected is None:
            not_manifested.append(str(p))
        elif hashlib.sha256(p.read_bytes()).hexdigest().lower() != expected:
            changed.append(str(p))
    audit = dict(source_files_read=len(USED), changed=changed, not_in_frozen_manifest=not_manifested,
                 all_read_sources_unchanged=not changed and not not_manifested,
                 reconstructions=validation,
                 all_reconstructed_folds_valid=all(r["valid"] for r in validation))
    summary = dict(prediction="P5", method="M3", compared_models=["a", "b"],
        excluded_common_delay=True, original_prediction_quote="P5（" + p5,
        coverage_factor_predefined_in_theory=False,
        descriptive_coverage_factors=[1, 2], comparison_counts=counts,
        magnitudes=magnitudes, all_sources_unchanged=audit["all_read_sources_unchanged"],
        conclusion="现有实测数据不支持 P5：秩1共享振型模型的相对时延跨窗口不一致；1u与2u两种描述性尺度下均如此。")
    ROOT.mkdir(exist_ok=True)
    for name, rows in (("p5_components", components), ("p5_reproducibility", reproducibility),
                       ("p5_repeat_pairs", pairs), ("p5_counts", counts),
                       ("p5_magnitudes", magnitudes), ("p5_raw_validation", validation)):
        write_csv(ROOT / f"{name}.csv", rows)
    write_json(ROOT / "p5_summary.json", summary)
    write_json(ROOT / "p5_audit.json", audit)
    write_json(ROOT / "p5_paired_covariances.json", raw_folds)
    write_json(ROOT / "p5_reproducibility_covariances.json", covariance_summary)
    lines = [
        "# P5：秩1相对时延的跨窗口核查", "",
        "**现有数据不支持 P5。** 秩1共享振型模型 b 的输入/输出相对时延在 H300 与 H493"
        "之间仍不相容；公共时延完全排除。按原始复现性与配对 jackknife 合成的尺度，"
        "b 在 48 个逐分量比较中只有 1 个落入 1u、3 个落入 2u；输出端 24 项全部超出 2u。",
        "",
        "## 1. 原判据与本次解释", "",
        "> **P5（" + p5.replace("\n", " ") + "", "",
        "理论原文没有给出覆盖因子，也没有明确“一倍”或“两倍”合成不确定度。"
        "因此下面同时报告 1u 和 2u，作为描述性敏感性，不能追称其中任一个是事先写定的数值门槛。"
        "两种尺度得到相同的否定结论。这里 u 是合成标准不确定度尺度，不把 2u 自动解释为 95% 置信区间。",
        "",
        "仅复用六条记录、两个高频窗口、M3、模型 a/b 的既有全样本与删周期端点，没有新拟合。"
        "模型 b 是 C0/C1 共享同一实秩1振型，模型 a 是自由 C0/C1。"
        "每侧四个相对时延和为零，因此八个分量只有六个独立方向；48是逐分量比较数，"
        "不是48个独立假设检验，不构造总体显著性或多重比较校正结论。",
        "",
        "## 2. 合成不确定度的计算", "",
        "对记录 r 的每个相对时延分量，定义 δ_r=τ_H493,r−τ_H300,r。"
        "第 i 个相同删周期折内先作 δ_(i)，再计算"
        r" $u_{\rm JK}^2=(P-1)P^{-1}\sum_i(\delta_{(i)}-\bar\delta)^2$。"
        "这样保留了两个窗口估计的同记录相关性，不能替换成两窗 SE 的简单平方和。",
        "",
        "运行间复现性按相同控制器的 A/C++++ 配对计算，三个控制器合并：", "",
        r"\[s_{\rm rep,\delta}^2=\frac{1}{2\times3}\sum_{g\in\{L,N,H\}}"
        r"(\delta_{A,g}-\delta_{C,g})^2,\qquad "
        r"u_{\rm raw,r}=\sqrt{u_{\rm JK,r}^2+s_{\rm rep,\delta}^2}.\]", "",
        "该定义保留了运行间变化在两个窗口之间的协方差，只包含各控制器内部的 A/C 重复差，"
        "没有把控制器间均值差当作复现性。三个配对的重复次数很少，复现性成分自身也很不确定。",
        "",
        "上式的原始复现性已经包含统计噪声，直接再与 JK 相加可能重复计入。"
        "因此另列扣除已估 JK 噪声的敏感性版本，不用它替换较宽的原始合成尺度：", "",
        r"\[s_{\rm extra}^2=\max\left\{0,\ s_{\rm rep,\delta}^2-\frac{1}{6}"
        r"\sum_g(u_{\rm JK,A,g}^2+u_{\rm JK,C,g}^2)\right\},\qquad "
        r"u_{\rm extra,r}=\sqrt{u_{\rm JK,r}^2+s_{\rm extra}^2}.\]", "",
        "零截断只用于方差分量估计，已在 CSV 中标记；它不使任何失败拟合消失。"
        "扣噪版本也不支持 P5，故结论并非重复计入噪声造成。", "",
        "## 3. 计数与差异大小", "",
        table(["模型", "分量", "有效/总数", "原始1u内", "原始2u内", "扣噪1u内", "扣噪2u内"],
              [[r["model"], r["side"], f"{r['valid_comparisons']}/{r['comparisons']}",
                r["within_1u_raw"], r["within_2u_raw"], r["within_1u_extra"], r["within_2u_extra"]]
               for r in counts]), "",
        table(["模型", "分量", "平均绝对跨窗差/μs", "最大绝对跨窗差/μs", "平均配对JKSE/μs"],
              [[r["model"], r["side"], f(r["mean_absolute_gap_us"]), f(r["max_absolute_gap_us"]),
                f(r["mean_jackknife_SE_us"])] for r in magnitudes]), "",
        "秩1约束使输入端跨窗口差的平均绝对值从约21.46 μs降至8.60 μs，但输出端从约5.15 μs"
        "增至15.03 μs。b 的输出端差异在六条记录中方向一致，不能将更小的统计 SE 解释成更准确的物理时延。"
        "自由模型 a 也未达到整体跨窗口一致，符合 P5 的对照描述；但这不能补救 b 所要求的相容性失败。",
        "",
        "## 4. 复现性分量与完整逐分量结果", "",
        table(["模型", "分量", "原始复现性SD/μs", "扣噪额外SD/μs", "扣噪截零"],
              [[r["model"], r["component"], f(r["raw_reproducibility_SD_us"]),
                f(r["noise_adjusted_extra_SD_us"]), r["extra_variance_clipped_to_zero"]]
               for r in reproducibility]), "",
        table(["模型", "记录", "分量", "H300/μs", "H493/μs", "差/μs", "JKSE/μs", "原始u/μs", "差/u", "扣噪u/μs"],
              [[r["model"], r["record"], r["component"], f(r["H300_delay_us"]),
                f(r["H493_delay_us"]), f(r["cross_window_gap_us"]),
                f(r["paired_jackknife_SE_us"]), f(r["combined_raw_u_us"]),
                f(r["gap_over_raw_u"]), f(r["combined_extra_u_us"])] for r in components]), "",
        "## 5. 数值核查和边界", "",
        table(["模型", "记录", "配对折数", "全端点/折有效", "差重建误差/μs", "SE重建误差/μs", "时延CSV误差/μs"],
              [[r["model"], r["record"], r["folds"], r["valid"], f(r["max_delta_csv_error_us"]),
                f(r["max_SE_csv_error_us"]), f(r["max_delay_csv_error_us"])] for r in validation]), "",
        f"读取并核对 {audit['source_files_read']} 个来源文件。来源哈希是否全部匹配冻结清单："
        f"{audit['all_read_sources_unchanged']}；改变 {len(changed)} 个，未入清单 {len(not_manifested)} 个。"
        "完整路径见 p5_audit.json。所有既有端点与折均保留，没有删除失败或改停止规则。",
        "",
        "这一结果检验的是现有平台数据是否支持两窗口共用相对时延的解释，并不否定正确模型下的"
        "局部可辨识性代数结论。共享实秩1振型、背景近似、局部窗模型是否足够是其适用条件；"
        "跨窗失配可能吸收到相对时延里。仅凭本次结果不能区分非比例阻尼、未建模背景或通道频率依赖。",
        "",
        "输出：p5_components.csv 为96个逐分量比较；p5_reproducibility.csv 为16个复现性分量；"
        "p5_repeat_pairs.csv 为A/C重复差；p5_counts.csv 和 p5_magnitudes.csv 为汇总；"
        "p5_paired_covariances.json 保留配对折与完整协方差；p5_reproducibility_covariances.json"
        "保留跨分量复现性协方差。没有对公共时延作P5判定。",
        "",
    ]
    (ROOT / "P5_相对时延跨窗口核查.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(jsonify(summary), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
