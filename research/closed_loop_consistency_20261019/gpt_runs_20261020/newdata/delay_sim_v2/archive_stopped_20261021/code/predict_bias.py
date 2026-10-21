"""Population direct-IV bias diagnostics; no fits and no original-file writes.

Random-phase averaging is used for the excitation.  Transient columns are kept
as a separate zero-mean-independent population block.  This is not the
conditional moment matrix for one realised multisine record.
"""
from sim_core import *
import argparse
import warnings

SD = .2
SV = 3e-7
SW = .03
CUTS = (1e-10, 1e-12, 1e-14, 1e-16)
OUT = ROOT / "analysis"


def verify_kernel():
    m = SimModel(550.)
    bt = m.truth()
    f = np.array([12., 60., 286.573, 344.945, 490.39])
    c = m.full.context(f)
    G, Gj = m.full.gmodel(c, bt)
    K = ktf(f)
    S = np.linalg.solve(np.eye(4)[None]+K@G,
                        np.broadcast_to(np.eye(4), (len(f), 4, 4)))
    F = S@K
    Phi = SV**2*np.eye(4)[None]+SW**2*(G@G.conj().transpose(0, 2, 1))
    L = np.linalg.cholesky(Phi)
    U_noise = -F@L
    Jn = np.einsum("koip,kil->kolp", Gj, U_noise, optimize=True)
    _, Jzero = m.output(c, bt, np.zeros((len(f), 4)))
    dense = np.stack([
        m.output(c, bt, U_noise[:, :, i])[1]-Jzero
        for i in range(4)], axis=2)
    btrace = -np.einsum("koip,kio->kp", Gj, F@Phi, optimize=True).real
    bcolumns = np.einsum("kolp,kol->kp", Jn.conj(), L, optimize=True).real
    common_formula = -2*np.pi*f*np.imag(np.trace(G@F@Phi, axis1=1, axis2=2))*old.UNIT
    out = dict(frequencies_hz=f,
        Jacobian_noise_factor_relative_error=float(
            la.norm(dense-Jn)/la.norm(Jn)),
        bias_trace_vs_latent_noise_columns_relative_error=float(
            la.norm(btrace-bcolumns)/la.norm(btrace)),
        common_delay_formula_5_1_relative_error=float(
            la.norm(btrace[:, m.delay_slice.start]-common_formula) /
            la.norm(common_formula)))
    out["passed"] = all(out[k] < 1e-10 for k in (
        "Jacobian_noise_factor_relative_error",
        "bias_trace_vs_latent_noise_columns_relative_error",
        "common_delay_formula_5_1_relative_error"))
    js(OUT/"predictions_verification.json", out)
    assert out["passed"], out


def pseudo(A, cut):
    u, s, vh = la.svd(A, full_matrices=False)
    keep = s > cut * s[0] if s[0] else np.zeros(len(s), bool)
    inv = (vh[keep].T / s[keep]) @ u[:, keep].T
    return inv, int(sum(keep)), s


def transform_percent(x):
    with np.errstate(over="ignore", invalid="ignore"):
        return 100 * np.expm1(x)


def population(step, periods=20, chunk=160):
    begin = time.perf_counter()
    m = SimModel(550.)
    bt = m.truth()
    ti = np.array([o*m.full.L+10*j+q for o in range(4)
                   for j in range(8) for q in (8, 9)])
    di = np.setdiff1d(np.arange(m.nb), ti)
    n = int(round((550-3)/step))
    f = np.linspace(3., 550., n+1)
    assert np.max(abs(f / step - np.round(f / step))) < 1e-7
    terms = {band: {key: np.zeros((m.nb, m.nb)) for key in
                    ("M_signal", "M_feedback_noise", "M_transient")}
             for band in BANDS}
    vectors = {band: np.zeros(m.nb) for band in BANDS}
    N = int(periods * P0)
    for start in range(0, len(f), chunk):
        ff = f[start:start+chunk]
        c = m.full.context(ff)
        G, Gjfull = m.full.gmodel(c, bt)
        Gj = Gjfull[:, :, :, di]
        K = ktf(ff)
        S = np.linalg.solve(np.eye(4)[None]+K@G,
                            np.broadcast_to(np.eye(4), (len(ff), 4, 4)))
        F = S @ K
        Phi = SV**2*np.eye(4)[None]+SW**2*(G@G.conj().transpose(0, 2, 1))
        L = np.linalg.cholesky(Phi)
        shape = cl.dshape(np.exp(2j*np.pi*ff/FS))
        A_signal = S * (SD*np.sqrt(shape))[:, None, None]
        A_noise = F @ L
        Js = np.einsum("koip,kil->kolp", Gj, A_signal, optimize=True)
        Jn = np.einsum("koip,kil->kolp", Gj, A_noise, optimize=True)
        # b_l=-Re tr(G_l F Phi), with W=I and G including every delay.
        bv = -np.einsum("koip,kio->kp", Gj, F@Phi, optimize=True).real
        _, Jzero = m.output(c, bt, np.zeros((len(ff), 4)), True)
        T = Jzero[:, :, ti]
        for band, (lo, hi) in BANDS.items():
            ix = (ff >= lo-1e-10) & (ff <= hi+1e-10)
            if not ix.any():
                continue
            weights = np.full(sum(ix), step/FS)
            fi = ff[ix]
            weights[np.isclose(fi, lo, rtol=0, atol=1e-8)] *= .5
            weights[np.isclose(fi, hi, rtol=0, atol=1e-8)] *= .5
            for key, arr in (("M_signal", Js), ("M_feedback_noise", Jn)):
                X = (arr[ix]*np.sqrt(weights)[:, None, None, None]).reshape(-1, len(di))
                terms[band][key][np.ix_(di, di)] += (X.conj().T @ X).real
            # Unnormalised DFT moment matrices are divided by N^2.
            # Hence deterministic transient columns have integral weight 1/N.
            Xt = (T[ix]*np.sqrt(weights/N)[:, None, None]).reshape(-1, len(ti))
            terms[band]["M_transient"][np.ix_(ti, ti)] += (Xt.conj().T@Xt).real
            vectors[band][di] += weights @ bv[ix]
    all_result = {}
    for band, pieces in terms.items():
        # Match SimModel(high=100) coordinates without recomputing the FRF.
        scale = np.ones(m.nb)
        if BANDS[band][1] != 550.:
            for o in range(4):
                scale[o*m.full.L+84:o*m.full.L+88] = (550/BANDS[band][1])**2
        for key in pieces:
            pieces[key] *= scale[:, None]*scale[None, :]
        vectors[band] *= scale
        M = sum(pieces.values())
        all_result[band] = dict(M=M, b=vectors[band], **pieces)
    return m, all_result, dict(step_hz=step, frequency_points=len(f),
                               periods_for_transient_scaling=periods,
                               seconds=time.perf_counter()-begin)


def diagnose(m, data):
    M, b = data["M"], data["b"]
    scale = 1/np.sqrt(np.maximum(abs(np.diag(M)), 1e-300))
    A = M*scale[:, None]*scale[None, :]
    A = (A+A.T)/2
    g = b*scale
    eig = la.eigvalsh((A+A.T)/2)
    u, s, vh = la.svd(A)
    ta = np.arange(m.delay_slice.start, m.delay_slice.stop)
    nu = np.setdiff1d(np.arange(m.nb), ta)
    high = np.r_[m.iw+np.array([3, 4, 5]), m.iz+np.array([3, 4, 5])]
    q = np.r_[high, ta]
    other = np.setdiff1d(np.arange(m.nb), q)
    result = dict(normalized_singular_values=s, normalized_eigenvalue_min=eig[0],
                  normalized_condition=s[0]/s[-1],
                  parameter_count=m.nb, high_parameter_indices=high,
                  high_frequencies_hz=cl.FM[3:6], cutoffs=[])
    # A factor-level Schur computation avoids subtraction of M_nv^{-1} terms.
    # It is a check of the algebra, not a cure for a nearly singular model.
    try:
        L = la.cholesky((A+A.T)/2, lower=True)
        xc = la.cho_solve((L, True), g)
        perm = np.r_[nu, ta]
        LP = la.cholesky(A[np.ix_(perm, perm)], lower=True)
        nnu = len(nu)
        Lt = LP[nnu:, nnu:]
        rhs = g[ta]-LP[nnu:, :nnu]@la.solve_triangular(
            LP[:nnu, :nnu], g[nu], lower=True)
        xt = la.cho_solve((Lt, True), rhs)
        permh = np.r_[other, q]
        LH = la.cholesky(A[np.ix_(permh, permh)], lower=True)
        no = len(other)
        Lq = LH[no:, no:]
        bh = g[q]-LH[no:, :no]@la.solve_triangular(
            LH[:no, :no], g[other], lower=True)
        own = la.cho_solve((Lq[:6, :6], True), bh[:6])
        propagated = -la.solve_triangular(
            Lq[:6, :6].T, Lq[6:, :6].T@xc[ta], lower=False)
        physical = scale*xc
        result["factor_schur_check"] = dict(
            available=True,
            common_delay_error_us=physical[ta[0]]*old.UNIT*1e6,
            delay_coordinate_error_us=physical[ta]*old.UNIT*1e6,
            modal_error_percent=transform_percent(physical[m.iw:m.iz+8]),
            maximum_modal_log_change=float(np.max(abs(physical[m.iw:m.iz+8]))),
            normalized_equation_relative_residual=float(
                la.norm(A@xc-g)/max(la.norm(g), 1e-300)),
            formula_5_3_delay_difference_us=float(
                la.norm(scale[ta]*(xt-xc[ta]))*old.UNIT*1e6),
            formula_5_4_own_log=scale[high]*own,
            formula_5_4_delay_propagation_log=scale[high]*propagated,
            formula_5_4_full_solution_log=physical[high],
            formula_5_4_decomposition_error=float(
                la.norm(scale[high]*(own+propagated-xc[high]))),
            maximum_delay_phase_linearization_radians_550Hz=float(
                2*np.pi*550*np.max(abs(old.UNIT*np.einsum(
                    "oip,p->oi", old.DESIGN, physical[ta])))),
            caveat="Successful Cholesky is not evidence of reliable "
                   "information below floating-point resolution. Compare "
                   "cutoff and grid sensitivity; large phase/log changes "
                   "invalidate an interpretation as a small-bias prediction.")
    except la.LinAlgError as exc:
        result["factor_schur_check"] = dict(available=False, reason=str(exc))
    for cut in CUTS:
        keep = s > cut*s[0]
        x = vh[keep].T @ ((u[:, keep].T@g)/s[keep])
        delta = scale*x
        pnu, ranknu, _ = pseudo(A[np.ix_(nu, nu)], cut)
        St = A[np.ix_(ta, ta)]-A[np.ix_(ta, nu)]@pnu@A[np.ix_(nu, ta)]
        beff = g[ta]-A[np.ix_(ta, nu)]@pnu@g[nu]
        pt, rankt, st = pseudo(St, cut)
        dt = scale[ta]*(pt@beff)
        po, ranko, _ = pseudo(A[np.ix_(other, other)], cut)
        barM = A[np.ix_(q, q)]-A[np.ix_(q, other)]@po@A[np.ix_(other, q)]
        barb = g[q]-A[np.ix_(q, other)]@po@g[other]
        ph, rankh, sh = pseudo(barM[:6, :6], cut)
        own = scale[high]*(ph@barb[:6])
        propagated = -scale[high]*(ph@barM[:6, 6:]@x[ta])
        total = own+propagated
        result["cutoffs"].append(dict(
            rcond=cut, rank=int(sum(keep)),
            normalized_equation_relative_residual=float(np.linalg.norm(A@x-g)/max(np.linalg.norm(g), 1e-300)),
            delta=delta, modal_log_change=delta[m.iw:m.iz+8],
            modal_error_percent=transform_percent(delta[m.iw:m.iz+8]),
            common_delay_error_us=delta[ta[0]]*old.UNIT*1e6,
            delay_coordinate_error_us=delta[ta]*old.UNIT*1e6,
            formula_5_3=dict(nuisance_rank=ranknu, delay_schur_rank=rankt,
                delay_schur_singular_values=st,
                schur_delay_coordinate_error_us=dt*old.UNIT*1e6,
                full_vs_schur_delay_norm_us=float(np.linalg.norm(dt-delta[ta])*old.UNIT*1e6)),
            formula_5_4=dict(other_nuisance_rank=ranko, high_schur_rank=rankh,
                high_schur_singular_values=sh,
                own_and_other_coupling_log=own,
                shared_delay_propagation_log=propagated,
                total_log=total, full_solution_log=delta[high],
                decomposition_discrepancy_norm=float(np.linalg.norm(total-delta[high])),
                own_component_exp_minus_one_percent=transform_percent(own),
                propagated_component_exp_minus_one_percent=transform_percent(propagated),
                total_exp_minus_one_percent=transform_percent(total),
                note="Only the log-coordinate components add. Their separately "
                     "converted exp-minus-one percentages do not add.")))
    return result


def summarise(allruns):
    lines = [
        "# 时延偏差总体矩阵与截断敏感性",
        "",
        "这是原 375 参数模型的随机相位平均总体积分诊断，不是拟合结果，"
        "也不是某一条随机相位记录的条件预测。真值 C1=D0=D2=0，未人为添加背景。",
        "",
        r"采用 $\Phi_v=\sigma_v^2I+\sigma_w^2G_0G_0^H$、"
        r"$S=(I+KG_0)^{-1}$、$F=SK$；输入信号谱为"
        r"$\Phi_{u,d}=\sigma_d^2dshape\,SS^H$，反馈输入噪声谱为"
        r"$\Phi_{u,n}=F\Phi_vF^H$。总体 $M$ 同时包含两项，"
        r"$b_l=-\int\Re\operatorname{tr}(G_lF\Phi_v)\,df/f_s$。",
        "",
        "保留全部模态瞬态列。随机相位平均将它们与动态参数的交叉块置零，"
        "瞬态自身 Gram 按原始 DFT 的 N² 归一化得到 1/N 系数，b 的瞬态块为零。"
        "这忽略有限记录状态边界与噪声的相关项，不能当作严格有限样本恒等式。",
        "",
        "| 积分步长 Hz | 频带 | 截断 | 秩/375 | 公共时延预测 μs | 5.3与完整解差的范数 μs |",
        "|---:|:---:|---:|---:|---:|---:|"]
    for run in allruns:
        for band, d in run["diagnostics"].items():
            for row in d["cutoffs"]:
                lines.append(
                    f"| {run['integration']['step_hz']:.8g} | {band} | {row['rcond']:.0e} | "
                    f"{row['rank']}/375 | {row['common_delay_error_us']:.7g} | "
                    f"{row['formula_5_3']['full_vs_schur_delay_norm_us']:.7g} |")
    lines += ["", "## 高频参数分块",
        "",
        "高频目标为 286.573/344.945/490.390 Hz 的 log频率、log阻尼。"
        "先消去其余所有留数、背景、瞬态及其他模态极点，才分成",
        r"$\delta\theta_H=\bar M_{HH}^{-1}\bar b_H"
        r"-\bar M_{HH}^{-1}\bar M_{H\tau}\delta\tau$。",
        "JSON 中分别保留自身/其他耦合项、共享时延传播项、完整解以及分解误差。"
        "本平台 100–550 Hz 的回路增益并不满足 <0.1，因此不能预设第一项为零。",
        "",
        "伪逆对全矩阵与各 Schur 子块分别截断时不一定给出相同解。"
        "若 5.3 与完整解差大，或改变截断/积分步长即大幅改变预测，"
        "应报告预测不能稳定计算，不能择取一个截断值作为精确理论预测。",
        "",
        "另用 Cholesky 因子直接计算 Schur 块，避免显式逆与大数相减，"
        "核对理论分块与完整解。它是代数实现自检；病态矩阵偶然能完成"
        "Cholesky，并不证明所有弱方向可可靠辨识。",
        "",
        "| 步长 Hz | 频带 | 因子分块可算 | 公共τ μs | 5.3一致性 μs | 5.4分解差（log范数） |",
        "|---:|:---:|:---:|---:|---:|---:|"]
    for run in allruns:
        for band, dd in run["diagnostics"].items():
            cc = dd["factor_schur_check"]
            if cc["available"]:
                lines.append(
                    f"| {run['integration']['step_hz']:.8g} | {band} | 是 | "
                    f"{cc['common_delay_error_us']:.7g} | "
                    f"{cc['formula_5_3_delay_difference_us']:.3g} | "
                    f"{cc['formula_5_4_decomposition_error']:.3g} |")
            else:
                lines.append(
                    f"| {run['integration']['step_hz']:.8g} | {band} | 否 | — | — | — |")
    finest = allruns[-1]["diagnostics"]["W"]["factor_schur_check"]
    if finest["available"]:
        own = transform_percent(np.asarray(finest["formula_5_4_own_log"]))
        prop = transform_percent(np.asarray(finest["formula_5_4_delay_propagation_log"]))
        total = transform_percent(np.asarray(finest["formula_5_4_full_solution_log"]))
        lines += ["", "最细积分下的宽带形式分解如下；各项分别按 exp−1 换算，"
                  "因此百分比列不能直接相加。该表不构成局部展开有效性的证据。",
                  "",
                  "| 目标频率 Hz | 参数 | 自身/其他耦合项 % | 共享时延项 % | 完整解 % |",
                  "|---:|:---:|---:|---:|---:|"]
        for j in range(6):
            lines.append(
                f"| {cl.FM[3+j%3]:.3f} | {'频率' if j<3 else '阻尼'} | "
                f"{own[j]:.7g} | {prop[j]:.7g} | {total[j]:.7g} |")
        lines += ["",
            f"宽带解的最大通道时延相位改变量（550 Hz）为 "
            f"{finest['maximum_delay_phase_linearization_radians_550Hz']:.4g} rad；"
            f"全部八模态的最大 |δlog参数| 为 "
            f"{finest['maximum_modal_log_change']:.4g}。"
            "这些量不处于小扰动范围；特别是频带外两项模态可产生极大形式偏差。"]
    lines += [
        "",
        "宽带的形式完整解若给出毫秒级时延改变量，应同时检查 JSON 中"
        "最高频率相位改变量。该量远大于1时，局部一阶展开已不支持"
        "把这个数值当作真实估计偏差。不能用它检验MC偏差比是否处于"
        "0.67–1.5而省略展开前提。",
        "",
        "## 与单条实测 D 矩阵版本的区别",
        "",
        "给定一条长度 N 的多正弦 D_k 时，信号部分应使用实际 "
        "U_d,k=S_kD_k，构造动态参数列与瞬态列，累加 Re(J_dᴴJ_d)/N²；"
        "它包含跨参数及动态—瞬态交叉块。再叠加反馈噪声期望项和 b。"
        "本脚本使用独立随机相位平均的 E[D_kD_kᴴ]，"
        "因此不对应 seed=1002 的特定相位，但可以检查总体偏差预测的秩和病态性。",
        "",
        "## 文件",
        "",
        "- predictions_step_*.npz：每个频带的 M、b、M_signal、M_feedback_noise、M_transient。",
        "- predictions_step_*.json：全部截断、5.3/5.4分块结果。",
        "- predictions_convergence.json：两档积分的矩阵和向量收敛差。",
        "- 所有频率/阻尼预测百分比统一按 100×expm1(δlog参数) 输出。"]
    (OUT/"predictions.md").write_text("\n".join(lines)+"\n", encoding="utf8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", default="0.16666666666666666,0.08333333333333333")
    ap.add_argument("--periods", type=int, default=20)
    ap.add_argument("--reuse", action="store_true",
                    help="Reuse this script's saved population matrices.")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    verify_kernel()
    runs, matrices = [], []
    for step in map(float, args.steps.split(",")):
        tag = f"{step:.9f}".replace(".", "p")
        archive_path = OUT/f"predictions_step_{tag}.npz"
        if args.reuse and archive_path.exists():
            m = SimModel(550.)
            arr = np.load(archive_path)
            data = {band: {key: arr[f"{band}_{key}"] for key in
                          ("M", "b", "M_signal", "M_feedback_noise", "M_transient")}
                    for band in BANDS}
            meta = rd(OUT/f"predictions_step_{tag}.json")["integration"]
        else:
            m, data, meta = population(step, args.periods)
        archive = {f"{band}_{key}": value for band, x in data.items()
                   for key, value in x.items()}
        np.savez_compressed(OUT/f"predictions_step_{tag}.npz", **archive)
        diagnostics = {band: diagnose(m, x) for band, x in data.items()}
        run = dict(integration=meta, diagnostics=diagnostics)
        js(OUT/f"predictions_step_{tag}.json", run)
        runs.append(run); matrices.append(data)
        print("completed", meta, flush=True)
    convergence = []
    for i in range(1, len(matrices)):
        row = dict(coarse=runs[i-1]["integration"]["step_hz"],
                   fine=runs[i]["integration"]["step_hz"], bands={})
        for band in BANDS:
            a, b = matrices[i-1][band], matrices[i][band]
            d = 1/np.sqrt(np.maximum(np.diag(b["M"]), 1e-300))
            row["bands"][band] = dict(
                normalized_M_relative_difference=float(la.norm(
                    (a["M"]-b["M"])*d[:, None]*d[None, :]) /
                    la.norm(b["M"]*d[:, None]*d[None, :])),
                normalized_b_relative_difference=float(la.norm(
                    (a["b"]-b["b"])*d)/max(la.norm(b["b"]*d), 1e-300)))
        convergence.append(row)
    js(OUT/"predictions_convergence.json", convergence)
    summarise(runs)
    print("saved", OUT/"predictions.md", flush=True)


if __name__ == "__main__":
    main()
