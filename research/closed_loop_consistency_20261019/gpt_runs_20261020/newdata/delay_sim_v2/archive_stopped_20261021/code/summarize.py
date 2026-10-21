"""Assemble existing delay-theory checks and saved MC endpoints, without fitting."""
from __future__ import annotations
import csv
from datetime import datetime
import json
import math
from pathlib import Path
import sys
import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
BANDS = {"L": "3–100 Hz", "H": "100–550 Hz", "W": "3–550 Hz"}
FM = [10.6, 11.7, 16.1, 286.573, 344.945, 490.39, 718.958, 755.501]
DELAY_NAMES = ["common", *[f"output{j}" for j in range(1, 5)], *[f"input{j}" for j in range(1, 5)]]
TRUTH_DELAY = np.array([200., 100., -100., 100., -100., -100., 100., -100., 100.])
GROUPS = [(20, b, "M1") for b in ("L", "H", "W")] + [
    (p, b, "M3") for p in (10, 20, 40) for b in ("L", "H")]


def read(path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except (json.JSONDecodeError, OSError) as exc:
        return {"read_error": repr(exc), "file": str(path)}


def finite(v):
    try:
        return math.isfinite(float(v))
    except (ValueError, TypeError):
        return False


def f(v, digits=5):
    return f"{float(v):.{digits}g}" if finite(v) else "—"


def clean(v):
    if isinstance(v, dict):
        return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [clean(x) for x in v]
    if isinstance(v, (np.integer, np.bool_)):
        return v.item()
    if isinstance(v, (float, np.floating)):
        return float(v) if math.isfinite(v) else None
    return v


def json_write(path, value):
    path.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def csv_write(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(clean(rows))


def table(head, rows):
    def c(x):
        return str(x).replace("|", r"\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(head) + " |",
                      "| " + " | ".join("---" for _ in head) + " |",
                      *["| " + " | ".join(c(x) for x in r) + " |" for r in rows]])


def state(record):
    if record is None:
        return dict(saved=False, success=False, status="未完成", score=None, near_bound=None,
                    iterations=None, seconds=None)
    sol = record.get("endpoint", {}).get("solver", {})
    ok = bool(record.get("success", sol.get("success", False)) and sol.get("success", False)
              and not sol.get("near_bound", False))
    return dict(saved=True, success=ok,
        status=record.get("read_error") or record.get("exception") or sol.get("status", record.get("status", "未知")),
        score=sol.get("score"), near_bound=sol.get("near_bound"), iterations=sol.get("iterations"),
        original_complete_score=sol.get("original_complete_score"),
        function_evaluations=sol.get("nfev"), seconds=record.get("seconds"))


def errors(record):
    if not record:
        return {}
    ep = record.get("endpoint", {})
    output = {}
    dl = ep.get("delay", {})
    if dl:
        arr = [dl.get("common_us"), *dl.get("output_us", []), *dl.get("input_us", [])]
        if len(arr) == 9:
            for name, val, true in zip(DELAY_NAMES, arr, TRUTH_DELAY):
                if finite(val):
                    output[f"delay_{name}_us"] = val-true
    if finite(ep.get("common_delay_error_us")):
        output["delay_common_us"] = ep["common_delay_error_us"]
    vals = ep.get("error_percent", [])
    if len(vals) == 16:
        for j in range(8):
            if finite(vals[j]):
                output[f"frequency_mode{j+1}_percent"] = vals[j]
            if finite(vals[j+8]):
                output[f"damping_mode{j+1}_percent"] = vals[j+8]
    return output


def stat(values):
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    n = len(x)
    return dict(n=n, bias=float(x.mean()) if n else None,
                MC_SE=float(x.std(ddof=1)/np.sqrt(n)) if n > 1 else None,
                SD=float(x.std(ddof=1)) if n > 1 else None,
                RMSE=float(np.sqrt(np.mean(x*x))) if n else None)


def mc_tables():
    endpoints, statistics, groups, p4 = [], [], [], []
    by_group = {}
    for P, band, method in GROUPS:
        rows = []
        for r in range(40):
            name = f"P{P:02d}_r{r:03d}_{band}_{method}.json"
            raw = read(ROOT/"mc"/"raw"/name)
            refined = read(ROOT/"mc"/"refinements"/name)
            trial = read(ROOT/"mc_qr"/"raw"/name)
            # The formal campaign uses exactly one frozen numerical route.
            # Never fill a missing official endpoint from a diagnostic route.
            chosen = read(ROOT/"mc_qr_final"/"raw"/name)
            source = "mc_qr_final/raw" if chosen is not None else "missing_official"
            final_state, original_state, root_state, trial_state = [
                state(x) for x in (chosen, raw, refined, trial)]
            row = dict(P=P, band=band, method=method, replicate=r, seed=1002+r,
                       source=source, refinement_saved=refined is not None,
                       original_saved=original_state["saved"],
                       original_success=original_state["success"],
                       original_status=original_state["status"], original_score=original_state["score"],
                       diagnostic_full_root_success=root_state["success"],
                       diagnostic_full_root_status=root_state["status"],
                       diagnostic_full_root_score=root_state["score"],
                       aborted_qr_trial_saved=trial_state["saved"],
                       aborted_qr_trial_success=trial_state["success"],
                       aborted_qr_trial_status=trial_state["status"],
                       **final_state, **errors(chosen))
            # Diagnostics remain available but never enter the formal statistics.
            row.update({f"normal_diagnostic_{k}": v for k, v in errors(raw).items()})
            row.update({f"full_root_diagnostic_{k}": v for k, v in errors(refined).items()})
            endpoints.append(row)
            rows.append(row)
        g = dict(P=P, band=band, method=method, planned=40,
                 raw_saved=sum(r["original_saved"] for r in rows),
                 raw_success=sum(r["original_success"] for r in rows),
                 refinements_saved=sum(r["refinement_saved"] for r in rows),
                 diagnostic_full_root_success=sum(r["diagnostic_full_root_success"] for r in rows),
                 aborted_qr_trial_saved=sum(r["aborted_qr_trial_saved"] for r in rows),
                 final_saved=sum(r["saved"] for r in rows),
                 final_success=sum(r["success"] for r in rows),
                 final_failed=sum(r["saved"] and not r["success"] for r in rows),
                 missing=sum(not r["saved"] for r in rows),
                 raw_failures_not_refined=sum(r["original_saved"] and not r["original_success"]
                                             and not r["refinement_saved"] for r in rows))
        groups.append(g)
        by_group[P, band, method] = g
        variables = [f"delay_{j}_us" for j in DELAY_NAMES] + [
            f"{what}_mode{j}_percent" for what in ("frequency", "damping") for j in range(1, 9)]
        for source_tag, source_rows in (("final", rows),):
            for variable in variables:
                values = [r[variable] for r in source_rows if r["success"] and finite(r.get(variable))]
                s = stat(values)
                statistics.append(dict(P=P, band=band, method=method, endpoint_version=source_tag,
                    variable=variable, **s, planned=40,
                    saved=sum(r["saved"] for r in source_rows),
                    failed=sum(r["saved"] and not r["success"] for r in source_rows),
                    selected_success_subset=True))
    for band in ("L", "H"):
        g = by_group[40, band, "M3"]
        for component in DELAY_NAMES:
            variable = f"delay_{component}_us"
            ss = {P: next(r for r in statistics if r["P"] == P and r["band"] == band
                          and r["method"] == "M3" and r["endpoint_version"] == "final"
                          and r["variable"] == variable) for P in (10, 20, 40)}
            current = ss[40]
            computed = current["n"] > 1 and finite(current["MC_SE"])
            conditional_pass = bool(computed and abs(current["bias"]) < 2*current["MC_SE"])
            whole_ensemble = g["final_saved"] == 40 and g["final_success"] == 40 and current["n"] == 40
            all_periods_saved = all(s["saved"] == 40 for s in ss.values())
            trend_computed = all(s["n"] == 40 for s in ss.values())
            trend_pass = (abs(ss[10]["bias"]) >= abs(ss[20]["bias"]) >= abs(ss[40]["bias"])) if trend_computed else None
            if not all_periods_saved:
                decision = "未完成，暂不能判定"
            elif not whole_ensemble:
                decision = "未通过（数值不可检验）"
            elif not conditional_pass:
                decision = "未通过"
            else:
                decision = "通过"
            p4.append(dict(band=band, component=component,
                success_P10=ss[10]["n"], success_P20=ss[20]["n"], success_P40=current["n"],
                bias_P10=ss[10]["bias"], bias_P20=ss[20]["bias"], bias_P40=current["bias"],
                MC_SE_P40=current["MC_SE"], conditional_success_subset_2MCSE_pass=conditional_pass,
                all_40_valid=whole_ensemble, all_three_period_groups_saved=all_periods_saved,
                saved_P10=ss[10]["saved"], saved_P20=ss[20]["saved"], saved_P40=current["saved"],
                trend_computed=trend_computed,
                absolute_bias_nonincreasing_P10_P20_P40=bool(trend_pass) if trend_computed else None,
                decision=decision))
    return endpoints, groups, statistics, p4


def main():
    out = ROOT/"analysis"
    out.mkdir(exist_ok=True)
    p5 = read(ROOT/"p5"/"p5_summary.json") or {}
    align = read(ROOT/"alignment"/"alignment_summary.json") or {}
    structural = read(ROOT/"audit"/"stability_identifiability_audit.json") or {}
    verification = read(ROOT/"audit"/"simulation_validation.json") or {}
    runtime = read(ROOT/"RUN_STATUS.json") or {}
    ep, groups, stats, p4 = mc_tables()
    summary = dict(generated_at=datetime.now().astimezone().isoformat(),
        planned=360, formal_source="mc_qr_final/raw only; no fallback",
        normal_diagnostic_saved=sum(g["raw_saved"] for g in groups),
        full_root_diagnostic_saved=sum(g["refinements_saved"] for g in groups),
        aborted_qr_trial_saved=sum(g["aborted_qr_trial_saved"] for g in groups),
        final_saved=sum(g["final_saved"] for g in groups),
        final_success=sum(g["final_success"] for g in groups),
        final_failed=sum(g["final_failed"] for g in groups),
        incomplete=sum(g["missing"] for g in groups),
        all_formal_runs_saved=all(g["final_saved"] == 40 for g in groups),
        batch_complete=bool(runtime.get("complete", False)),
        batch_running=runtime.get("running"),
        runtime_status_updated=runtime.get("updated"),
        diagnostic_note="Normal equations, separate full roots, and the aborted initial QR trial are diagnostic only.",
        groups=groups, P4_component_checks=p4,
        P5=p5.get("conclusion"), p5_complete=bool(p5))
    json_write(out/"report_summary.json", summary)
    csv_write(out/"mc_endpoint_inventory.csv", ep)
    csv_write(out/"mc_groups.csv", groups)
    csv_write(out/"mc_statistics.csv", stats)
    csv_write(out/"p4_component_checks.csv", p4)
    prediction_rows = []
    for path in sorted(out.glob("predictions_step_*.json")):
        r = read(path) or {}
        for band, d in r.get("diagnostics", {}).items():
            for x in d.get("cutoffs", []):
                prediction_rows.append(dict(file=path.name, step_hz=r.get("integration", {}).get("step_hz"),
                    band=band, rcond=x.get("rcond"), rank=x.get("rank"),
                    common_delay_prediction_us=x.get("common_delay_error_us"),
                    full_schur_delay_difference_us=x.get("formula_5_3", {}).get("full_vs_schur_delay_norm_us"),
                    factor_schur_available=d.get("factor_schur_check", {}).get("available"),
                    factor_common_prediction_us=d.get("factor_schur_check", {}).get("common_delay_error_us"),
                    max_linearized_phase_radians=d.get("factor_schur_check", {}).get(
                        "maximum_delay_phase_linearization_radians_550Hz")))
    csv_write(out/"prediction_diagnostics.csv", prediction_rows)
    loop = structural.get("loop", {})
    information = structural.get("information", {})
    lines = [
        "# 时延理论仿真验证报告", "",
        f"汇总时间：{summary['generated_at']}。统一 QR/SVD 正式 MC 已保存 {summary['final_saved']}/360 个端点，"
        f"当前有效 "
        f"{summary['final_success']} 个、失败 {summary['final_failed']} 个、未完成 {summary['incomplete']} 个。",
        "",
        "**当前批次状态必须区分未完成和已完成但无法检验。** 尚未保存全部所需正式端点时，"
        "写“未完成，暂不能判定”；全部运行完成后仍因求根失败或信息不足而无法检验，才写"
        "“未通过（数值不可检验）”。成功子集均值不能当作总体无偏证据；数值不可检验也不等于理论被证伪。",
        "",
        "## 1. 本轮范围、预设预测与计算约定", "",
        "保留完整八模态、自由 C0/C1、自由 D0/D2、原逐模态瞬态和七个结构化时延，共375参数；"
        "真对象 C1=D0=D2=0。频带 L=3–100 Hz、H=100–550 Hz、W=3–550 Hz。"
        "所有八个极点均估计，包括目标频带外的718.958/755.501 Hz模态；没有通过冻结弱方向来降低难度。",
        "",
        "原计划共360次：P=20时 M1 在 L/H/W 各40次；P=10/20/40时 M3 在 L/H 各40次。"
        "M1使用全部频点，M3使用激励谱线。每周期30000点、5kHz，先模拟两个稳定周期再保留P周期。"
        "数据种子1002+r，初值扰动种子100200+r；频率log扰动标准差0.01、阻尼log扰动标准差0.4。"
        "时延坐标从零偏移开始，公共初值200μs恰好是真公共时延，但相对时延未使用真值。"
        "本仿真已知局部分支，没有运行实测的41点公共时延网格；这限制了全局收敛结论。",
        "",
        "**正式主统计统一且仅使用 `mc_qr_final/raw`。** 该路径无同名端点就标记未完成，"
        "绝不回退到正规方程、额外完整求根或试版 QR 路线。"
        "`mc/raw` 与 `mc/refinements` 只作原始数值失败诊断；"
        "`mc_qr/raw` 是计算求和顺序校验修改前启动即中止的试版，只列状态，不纳入正式统计。",
        "",
        "正式版在线性剖面内层把QR/SVD的1e−13截断施加在原设计矩阵，而不是条件数已平方的正规矩阵上；"
        "内层收敛统一检查原始未压缩线性矩的score<1e−9。完整原375维矩的score<1e−6与贴边判据保持不变。"
        "旧内层2e−12在H频带过严，属于本次明确修复的数值求解容差，不能声称没有改动任何容差；"
        "它不是预先冻结的科学结论判据。正式外层仍包含剖面完整矩hybrid求根，"
        "而旧路径额外375维完整求根只是诊断，不是正式主估计。"
        "这一内层修复不代表外层Schur计算也已避免正规矩阵；第4节单列其新发现的数值瓶颈。",
        "",
        "P2/P3的偏差比0.67–1.5、带外回路增益<0.1与2MCSE门槛不改变。P4仍要求"
        "P=40时偏差小于2MCSE并检查P=10/20/40趋势。P5原文未指定覆盖因子，另列1u/2u；"
        "不能把新选的覆盖因子冒称冻结判据。P1的秩1仿真MC不在当前360次安排中，不能据本轮称其已验证。",
        "",
        "## 2. 新旧记录对齐核查", "",
    ]
    if align:
        lines += [
            "当前记录证据不支持用“新旧离线准备相差一拍”解释约196μs的有效公共时延差。"
            "新B组电机指令与控制器输出同拍且逐值相等；旧cl_data中U/Y/V与原始导出逐值相同；"
            "新旧控制器递推最早必要传感器输入均为一拍（200μs）。",
            "",
            "ARX保留一拍输入时，独立半段相对残差约1e−11；删去一拍、只保留两拍及以后时"
            "恶化至约5e−5至3e−4。旧准备流程未发现移拍。"
            "这只验证记录信号边界，不是ADC/DAC或物理施力总时延标定；"
            "仍不能把局部模型有效时延直接解释为硬件时延。",
            "",
            "完整证据与限制：[记录对齐核查](alignment/记录对齐核查.md)。", "",
        ]
    else:
        lines += ["记录对齐核查尚未完成，不填结论。", ""]
    lines += ["## 3. 回路条件与闭环稳定性", ""]
    if loop:
        lines += [
            f"含已知输入/输出整数拍时延的闭环谱半径为 {f(loop.get('closed_loop_radius'), 10)}，"
            f"稳定性检查为 {loop.get('stable')}。", "",
            table(["频带", "max σmax(GK)", "峰值频率/Hz", "max谱半径(GK)", "max单元素模", "全带<0.1"],
                  [[band, f(v.get("sigma_max", {}).get("value")),
                    f(v.get("sigma_max", {}).get("frequency_hz")),
                    f(v.get("spectral_radius", {}).get("value")),
                    f(v.get("max_abs_element", {}).get("value")),
                    v.get("sigma_max", {}).get("below_point_one")] for band, v in loop.get("bands", {}).items()]),
            "",
            "**H=100–550Hz并不满足带内回路增益小于0.1。** 最大奇异值约0.714，"
            "即使改用谱半径或最大元素模也超0.1。高频不能仅因在控制器标称带宽以外就当作开环极限。"
            "因此P2/P3的带外零偏差子预测在当前配置中缺少前提；不改控制器或截频段来使其看似满足。",
            "",
        ]
        if (ROOT/"figures/audit_high_band_loop_gain.png").exists():
            lines += [
                "![高频真实回路最大奇异值](figures/audit_high_band_loop_gain.png)", "",
                "图：已保存真模型回路审计的100–550 Hz全部频点，虚线为预设0.1上限；"
                "标出286.6、344.9、490.4 Hz模态。这里没有新增拟合或MC平均。"
                "矢量版：[PDF](figures/audit_high_band_loop_gain.pdf)。", "",
            ]
    else:
        lines += ["回路条件审计未完成。", ""]
    lines += ["## 4. 弱辨识与模型求值核查", ""]
    if information:
        lines += [
            "下表为无噪声参考输入下、逐列归一化实雅可比的诊断，保留全部375列。"
            "Gram条件数会将雅可比条件数平方；在极弱方向上，直接形成Gram所得的微小负特征值属于"
            "数值精度问题，不能将其解释成负信息量。低频投影余量接近舍入误差，也不是精确ZOH结构恒等式的证明。",
            "",
            table(["频带", "归一化J条件数", "由J推算Gram条件数", "J数值秩(1e−12)", "时延未解释能量最小", "最大"],
                  [[band, f(v.get("full_normalized_J_condition")),
                    f(v.get("full_normalized_Gram_condition_from_J")),
                    v.get("full_normalized_J_rank", {}).get("1e-12"),
                    f(min(v.get("delay_unexplained_energy_fractions", [np.nan]))),
                    f(max(v.get("delay_unexplained_energy_fractions", [np.nan])))]
                   for band, v in information.get("bands", {}).items()]), "",
        ]
        if (ROOT/"figures/audit_band_information.png").exists():
            lines += [
                "![三频带雅可比与时延投影诊断](figures/audit_band_information.png)", "",
                "图：左为逐列归一化雅可比的奇异谱，右为时延列投影到其余参数列之外后剩余的能量比例。"
                "三个a与三个b坐标是零和相对时延的独立坐标，不是直接选四路中的前三路。"
                "该图是数值投影与条件性诊断；低频舍入量级余量不证明精确ZOH模型严格不可辨识，"
                "也不据此宣称全局不可辨识。矢量版：[PDF](figures/audit_band_information.pdf)。", "",
            ]
        for audit_note in sorted((ROOT/"audit").glob("*.md")):
            if audit_note.name != "仿真求值与导数核查.md":
                lines += [f"补充数学边界说明：[{audit_note.stem}](audit/{audit_note.name})。", ""]
    if verification:
        noiseless = verification.get("noiseless", {})
        lines += [
            table(["核查项", "结果"], [
                ["八模态传函与独立状态空间相对差", f(verification.get("response_relative_error"))],
                ["闭环 d→y 独立状态空间相对差", f(verification.get("state_y_from_d_relative_error"))],
                ["闭环 d→u 独立状态空间相对差", f(verification.get("state_u_from_d_relative_error"))],
                ["解析导数核查", verification.get("jacobian_all_passed")],
                ["无噪声谱残差相对量", f(noiseless.get("total_relative_residual"))],
                ["总体求值核查", verification.get("all_passed")],
            ]), "",
            "这些是模型求值与导数核查，不是从任意初值恢复全部375真参数的证明。"
            "时延寄存器可引入不属于原模态瞬态空间的有限记录边界项；本核查先稳定两周期且无噪声周期重复，"
            "该误差在数值上可以忽略，但不能据此称任意初始状态或含噪边界上省略项严格为零。",
            "",
        ]
    schur_audit = read(ROOT/"audit/schur_rank_check.json")
    if isinstance(schur_audit, list):
        lines += [
            "### 外层求根的独立数值诊断", "",
            "当前失败不能全部归因于弱辨识，更不能作为不存在IV根的证据。"
            "独立诊断发现，原高维设计在内层QR中保留的方向，会在外层正规Schur计算中再次被截断；"
            "弱的窗外模态随后产生异常大步长，共同缩放会把带内更新一并压低。"
            "下表是少数审计点，不是全部MC样本的原因统计。",
            "",
            table(["频带", "审计点", "tall设计秩/352", "Schur线性秩/352",
                   "Schur非线性秩/23", "Schur条件数", "共同步长因子"],
                  [[v["band"], v["point"], v["tall_X_rank"], v["schur_linear_rank"],
                    v["schur_nonlinear_rank"], f(v["schur_condition"]),
                    f(v["step_factor"])] for v in schur_audit]),
            "",
            "H代表端点中，内层设计满秩352，而外层Schur线性块只保留332个方向；"
            "约755Hz的窗外模态阻尼产生约−7765的原始log步长，使共同步长因子缩至约5.15×10⁻⁶，"
            "带内频率的实际步长仅约10⁻⁸量级。W频带虽在这几个审计点未截掉线性方向，"
            "也出现窗外阻尼大步造成共同缩步。这是可具体定位的算法瓶颈，不是无根或不一致性的证明。",
            "",
            "原始诊断：[Schur秩与步长核查](audit/schur_rank_check.json)。"
            "任何独立稳定求解器短试验仅作诊断，不改当前正式批次，也不混入主统计；"
            "本轮保留当前正式求解来源，现有成功门槛和未完成语义保持不变。",
            "",
        ]
    pilot = read(ROOT/"stable_profile_pilot/pilot_summary.json")
    if isinstance(pilot, list):
        lines += [
            "### 稳定投影与信赖域短试验：仅作数值诊断", "",
            "独立短试验将外层正规Schur改为投影构造并使用信赖域，随后对H做完整中央差分矩收尾。"
            "模型、375维原矩验收门槛、窗外模态及贴边判据均保留。"
            "下列是同一或少数数据样本上的不同尝试，不是独立Monte Carlo样本；所有尝试都保留，"
            "也不把其中最小score的结果回填到正式批次。",
            "",
            table(["独立短试验", "完整375维score", "内层score", "通过原门槛",
                   "贴边", "求解/全部profile次数", "秒", "退出状态"],
                  [[v["name"], f(v["score"]), f(v["linear_score"]), v["success"],
                    v["near_bound"], f"{v['calls']}/{v['all_calls']}",
                    f(v["seconds"]), v["status"]] for v in pilot]),
            "",
            "H的score从约0.373降至投影40步后的6.8077×10⁻⁶，再经完整矩收尾降至4.4351×10⁻⁶；"
            "W的终值约0.0013743。它们均未达到原1×10⁻⁶门槛，全部继续按失败诊断保存。"
            "改善表明原停滞包含实际数值机制；仍不能据此称已求得有效根、已验证一致性，"
            "或认为全部模态都无法准确估计。",
            "",
        ]
        relevant = [v for v in pilot if v["name"] in
                    ("v2_P10_r000_H_M3", "central_P10_r000_H_M3")]
        lines += [
            "H的带内模态单独列出如下。误差是相对真值的百分比，阻尼使用exp−1换算，"
            "不是阻尼比百分点。带内小误差不能替代全参数验收，也不能作为跨样本精度统计。",
            "",
            table(["短试验", "模态真值/Hz", "频率误差/%", "阻尼误差/%"],
                  [[v["name"], FM[j], f(v["error_percent"][j]),
                    f(v["error_percent"][j+8])]
                   for v in relevant for j in (3, 4, 5)]),
            "",
            "独立v2端点审计的最大10个score方向均涉及窗外第7/8模态或时延；"
            "因此“完整375参数未过门槛”与“带内三个模态拟合较好”可以同时成立。"
            "该审计只针对v2端点，不能自动推广到所有失败样本。"
            "详情：[端点方向核查](audit/stable_endpoint_check.json)；"
            "[全部数值修复短试验与预算](stable_profile_pilot/数值修复pilot说明.md)。",
            "",
            "**本轮不切换正式来源，也不启动新的MC路线。** 当前统一批次仍以"
            "`mc_qr_final/raw`为唯一来源并继续保留全部失败。独立短试验不计入360项完成度。"
            "正式批次未完成时P4仍暂不能判定；即使最终全部数值失败，"
            "也只能给出“未通过（数值不可检验）”，不能作为一致性理论的反证。",
            "",
        ]
    lines += [
        "## 5. Monte Carlo 完成度与全部失败", "",
        table(["P", "频带", "方法", "正式保存", "正式有效", "正式失败", "未完成",
               "正规方程诊断保存/成功", "额外root诊断保存/成功", "中止QR试版"],
              [[g["P"], g["band"], g["method"], g["final_saved"], g["final_success"],
                  g["final_failed"], g["missing"], f"{g['raw_saved']}/{g['raw_success']}",
                  f"{g['refinements_saved']}/{g['diagnostic_full_root_success']}",
                  g["aborted_qr_trial_saved"]] for g in groups]), "",
        "正式统计只取统一 `mc_qr_final/raw`。所有失败的端点、score、状态、边界及诊断路线状态都保存在"
        "`analysis/mc_endpoint_inventory.csv`，不删除异常值。下面的均值、SD、RMSE只使用通过求解的子集，"
        "明确不是无条件总体性能；若失败与数据相关，会产生选择效应。",
        "",
        "## 6. 时延估计结果与 P4", "",
        table(["P", "频带", "方法", "有效统计n/40", "公共时延偏差/μs", "MCSE/μs", "SD/μs", "RMSE/μs"],
              [[s["P"], s["band"], s["method"], f"{s['n']}/40", f(s["bias"]),
                f(s["MC_SE"]), f(s["SD"]), f(s["RMSE"])] for s in stats
               if s["endpoint_version"] == "final" and s["variable"] == "delay_common_us"]), "",
        "公共时延真值为200μs，输出相对时延为[100,−100,100,−100]μs，"
        "输入相对时延为[−100,100,−100,100]μs。"
        "公共及八个零和相对分量共七独立方向；下表九行分量不是九个独立检验。",
        "",
        table(["频带", "时延分量", "n(P10/P20/P40)", "偏差P10", "偏差P20", "偏差P40", "2MCSE(P40)", "判定"],
              [[r["band"], r["component"],
                f"{r['success_P10']}/{r['success_P20']}/{r['success_P40']}",
                f(r["bias_P10"]), f(r["bias_P20"]), f(r["bias_P40"]),
                f(2*r["MC_SE_P40"] if finite(r["MC_SE_P40"]) else None), r["decision"]] for r in p4]),
        "",
        "上表判定保留P40的2MCSE数值门槛；P10/P20/P40偏差趋势另列作描述性检查，"
        "不事后增设MC均值严格单调的通过条件。对应频带的P10/P20/P40正式MC尚未全部保存时，"
        "主判定一律是“未完成，暂不能判定”。全部保存后仍缺乏可靠数值端点才记为"
        "“未通过（数值不可检验）”，不能称为一致性理论被证伪。"
        "成功子集的条件性2MCSE门槛另存CSV，不能掩盖失败。",
        "",
        "## 7. 模态频率、阻尼及偏差传播", "",
        "全部频率/阻尼误差按100×expm1(估计log参数−真log参数)计算，不与log单位混用。"
        "完整八模态频率与阻尼正式统计在 `analysis/mc_statistics.csv`；"
        "诊断路线的端点误差仅保留于inventory，不参与该统计。", "",
        table(["P", "频带", "方法", "模态/Hz", "有效n", "阻尼偏差/%", "MCSE", "SD", "RMSE"],
              [[s["P"], s["band"], s["method"], FM[int(s["variable"].split("mode")[1].split("_")[0])-1],
                s["n"], f(s["bias"]), f(s["MC_SE"]), f(s["SD"]), f(s["RMSE"])]
               for s in stats if s["endpoint_version"] == "final" and s["variable"].startswith("damping_mode")]),
        "",
        "高频目标为286.573/344.945/490.390Hz。P3需要把自身/其他耦合项和共享时延传播项分开；"
        "当前H窗口回路并不小，不能预先把自身项置零，再把所有偏差归因于共享时延。",
        "",
        "## 8. P2/P3总体偏差计算与可检验性", "",
    ]
    p23_saved = sum(g["final_saved"] for g in groups if g["method"] == "M1")
    if p23_saved < 120:
        lines += [
            f"**P2/P3所需M1正式MC目前保存 {p23_saved}/120，比较仍待完成。** "
            "此时不对MC偏差比或2MCSE作通过/未通过判定。"
            "第3节的高频小回路前提不满足是独立完成的条件核查，不依赖这批MC是否跑完。",
            "",
        ]
    else:
        lines += [
            "P2/P3所需120个M1正式端点均已保存。是否可作偏差比较还需检查求解成功率、"
            "一阶预测数值稳定性与展开前提；保存齐全本身不是通过证据。",
            "",
        ]
    if prediction_rows:
        lines += [
            table(["积分步长/Hz", "频带", "rcond", "秩/375", "公共τ预测/μs", "完整/Schur差/μs",
                   "因子Schur可算", "因子公共τ/μs", "最大线性相位/rad"],
                  [[f(r["step_hz"]), r["band"], f(r["rcond"]), r["rank"],
                    f(r["common_delay_prediction_us"]), f(r["full_schur_delay_difference_us"]),
                    r["factor_schur_available"], f(r["factor_common_prediction_us"]),
                    f(r["max_linearized_phase_radians"])] for r in prediction_rows]), "",
            "表中全部截断选择均保留。伪逆分别截断全矩阵与Schur子块未必对应同一解；"
            "若预测随截断/积分变化很大，或线性化相位改变量不小，不能选一个有利值套入0.67–1.5偏差比。"
            "计算方法与逐项诊断见 [总体偏差计算说明](analysis/predictions.md)。", "",
        ]
    else:
        lines += ["总体偏差文件尚未保存；预测计算待完成，不填理论预测数值或偏差比。", ""]
    lines += [
        "当前不得把“完成矩阵运算”等同于“得到可靠的一阶预测”。是否满足局部展开、"
        "稳定数值秩和足够成功MC样本必须与偏差比一起检查。高频零偏差的<0.1回路前提已确认不满足，"
        "这与MC统计门槛是否可算是两个不同问题。",
        "",
        "## 9. P5：既有实测相对时延", "",
    ]
    if p5:
        lines += [p5.get("conclusion", ""), "",
                  table(["模型", "分量", "比较数", "1u内", "2u内", "扣噪1u内", "扣噪2u内"],
                        [[r["model"], r["side"], r["comparisons"], r["within_1u_raw"],
                          r["within_2u_raw"], r["within_1u_extra"], r["within_2u_extra"]]
                         for r in p5.get("comparison_counts", [])]), "",
                  "P5完全排除公共时延。两个窗口的相同删周期折内先做差再计算JK；"
                  "运行间复现性由三个控制器各自A/C重复差合并，保留跨窗口相关性。"
                  "1u/2u均不支持秩1模型跨窗口一致，扣除重复计入统计噪声也不改变判断。"
                  "48是逐分量比较数，非48个独立检验。", "",
                  "完整推导与源文件复核：[P5核查报告](p5/P5_相对时延跨窗口核查.md)。", ""]
    else:
        lines += ["P5结果尚未保存，未作判断。", ""]
    lines += [
        "## 10. 当前可写与不可写", "",
        "- 可报告：新旧记录时序核查、稳定性和回路条件、信息矩阵病态性、"
        "P5现有数据不支持，以及实际保存的MC成功/失败数量。",
        "- 不可把未完成MC写成完成；不可把成功子集的低偏差写成总体无偏；"
        "不可把带外小回路条件不成立时的结果当作相应推论的反例。",
        "- 本轮自由留数375参数的弱辨识可能阻止直接检验预定时延偏差命题。"
        "保持全部原条件与失败，不通过缩模型、固定弱参数或改噪声来凑结论。",
        "- 参考IV的一致性是正确模型、充分激励、可辨识与找到真值附近根下的条件性结果，"
        "本轮求值核查和实测相容性检查都不能替代这些条件。",
        "",
        "## 11. 输出与更新", "",
        "`code/summarize.py` 不导入拟合器，仅聚合已有文件。每次统一正式MC完成更多样本后可重新运行，"
        "覆盖本目录汇总，不改原始数据与端点。机器可读状态为 `analysis/report_summary.json`；"
        "全360计划索引、正式主统计与P4分量表均在 `analysis/`；诊断路线只作单列，不混入主统计。",
        "",
        "`code/watch_run.py` 是当前有限MC批次的伴随进程，每新增8个正式端点自动重建报告；"
        "批次进程退出时再重建一次并核验原文件哈希，输出 `analysis/source_hash_verification.json`。"
        "`RUN_STATUS.json` 提供实时数量、是否仍运行与 complete 标志；本报告是生成时快照。",
        "",
        f"本次读取的运行状态：running={runtime.get('running', '未记录')}，"
        f"complete={runtime.get('complete', False)}，更新时间={runtime.get('updated', '未记录')}。"
        "complete=false 时不得宣称本批次全部完成。",
        "",
    ]
    during_hash = read(ROOT/"analysis/source_hash_during_batch.json")
    if during_hash:
        lines += [
            f"运行中原文件审计：已核对 {during_hash.get('checked')}/{during_hash.get('total')} 个文件，"
            f"改变 {len(during_hash.get('changed', []))}、缺失 {len(during_hash.get('missing', []))}、"
            f"读取错误 {len(during_hash.get('errors', []))}。"
            "这是运行中快照，不是批次结束后的最终审计；最终哈希核查仍由监看进程在退出时执行。"
            "详情：[运行中哈希核查](analysis/source_hash_during_batch.json)。", "",
        ]
    (ROOT/"时延理论仿真验证报告.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("groups", "P4_component_checks")},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    # watch_run.py decodes this subprocess stream as UTF-8 on Windows.
    sys.stdout.reconfigure(encoding="utf-8")
    main()
