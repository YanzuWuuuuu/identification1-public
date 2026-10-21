"""Stable tall-matrix linear IV profile; original full moments remain the check.

This file alone defines the new route.  No original model, data, solver or
saved result is modified.  Numerical rank loss is reported explicitly.
"""
from sim_core import *
import argparse
import traceback

COND = 1e-13
IMPLEMENTATION_VERSION = "tall_qr_v3_compact_residual"


def real_rows(X):
    return np.concatenate((X.real, X.imag), axis=0)


def tall_factor(X, cond=COND):
    """Column scaling then QR and SVD of R; never construct X.T @ X."""
    scale = 1/np.maximum(la.norm(X, axis=0), 1e-300)
    Q0, R0 = la.qr(X*scale, mode="economic", check_finite=False)
    u, s, vh = la.svd(R0, full_matrices=False, check_finite=False)
    keep = s > cond*s[0] if s[0] else np.zeros(len(s), bool)
    # X approximately equals Q @ diag(s_kept) @ vh_kept @ diag(1/scale).
    Q = Q0@u[:, keep]
    return dict(Q=Q, s=s[keep], vh=vh[keep], scale=scale,
                rank=int(sum(keep)), columns=X.shape[1],
                singular_values=s,
                scaled_condition=float(s[0]/max(s[-1], 1e-300)))


def coefficients(factor, gamma):
    return factor["scale"]*(factor["vh"].T@(gamma/factor["s"]))


def linear_moment_check(pr, b):
    """Original uncompressed complex linear IV moments and normalisation."""
    u = pr.uh(b)
    yh, XL, _ = pr.m.full.compact(pr.c, b, pr.U, False)
    if pr.method in ("M1", "M2"):
        ZL = XL
    else:
        _, ZL, _ = pr.m.full.compact(pr.c, b, u, False)
    L = pr.m.full.L
    ee, gs, norms = 0., [], []
    for o in range(4):
        sl = slice(o*L, (o+1)*L)
        e = pr.Y[:, o]-yh[:, o]
        gs.append((ZL[:, o].conj().T@e).real)
        norms.append(np.sum(abs(ZL[:, o])**2, axis=0))
        ee += float(np.sum(abs(e)**2))
    g = np.concatenate(gs)
    zz = np.concatenate(norms)
    score = float(np.max(abs(g)/np.sqrt(np.maximum(zz*ee, 1e-300))))
    return score, ee


class QRProfileIV(ChunkProfileIV):
    """Override only the linear nuisance solve at fixed nonlinear coordinates."""
    inner_tol = 1e-9

    def profile(self, q, lin, maxinner=60):
        started = time.perf_counter()
        b = np.r_[lin, q]
        m, L = self.m, self.m.full.L
        _, XL, _ = m.full.compact(self.c, b, self.U, False)
        factors = [tall_factor(real_rows(XL[:, o])) for o in range(4)]
        yr = [real_rows(self.Y[:, o]) for o in range(4)]
        meta = dict(route="scaled tall QR/SVD; original complex moment check",
                    implementation_version=IMPLEMENTATION_VERSION,
                    rcond_on_tall_matrix=COND, normal_equations_used=False,
                    inner_linear_tolerance=self.inner_tol,
                    original_reference_inner_tolerance=2e-12,
                    final_full_root_tolerance_unchanged=1e-6,
                    X=[{k: ff[k] for k in ("rank", "columns",
                                          "scaled_condition", "singular_values")}
                       for ff in factors],
                    iterations=[], numerical_rank_deficient=False)
        ok = False
        if self.method in ("M1", "M2"):
            for o, fx in enumerate(factors):
                sl = slice(o*L, (o+1)*L)
                b[sl] = coefficients(fx, fx["Q"].T@yr[o])
            score, ee = linear_moment_check(self, b)
            ok = bool(score < self.inner_tol and np.isfinite(b).all())
            meta["iterations"].append(dict(
                iteration=1, original_linear_score=score,
                residual_squared=ee))
            niter = 1
        else:
            history = []
            for it in range(maxinner):
                u = self.uh(b)
                _, ZL, _ = m.full.compact(self.c, b, u, False)
                zmeta = []
                for o, fx in enumerate(factors):
                    fz = tall_factor(real_rows(ZL[:, o]))
                    C = fz["Q"].T@fx["Q"]
                    rhs = fz["Q"].T@yr[o]
                    gamma, _, rank, singular = la.lstsq(
                        C, rhs, cond=COND, lapack_driver="gelsd",
                        check_finite=False)
                    sl = slice(o*L, (o+1)*L)
                    b[sl] = coefficients(fx, gamma)
                    zmeta.append(dict(
                        Z_rank=fz["rank"],
                        Z_condition=fz["scaled_condition"],
                        cross_rank=int(rank), cross_shape=list(C.shape),
                        cross_singular_values=singular,
                        cross_relative_residual=float(
                            la.norm(C@gamma-rhs)/max(la.norm(rhs), 1e-300))))
                # Recompute Uh at the new b before accepting a reference-IV root.
                score, ee = linear_moment_check(self, b)
                history.append(score)
                meta["iterations"].append(dict(
                    iteration=it+1, original_linear_score=score,
                    residual_squared=ee, Z_and_cross=zmeta))
                if score < self.inner_tol and np.isfinite(b).all():
                    ok = True
                    break
                if not np.isfinite(b).all():
                    break
                if len(history) > 8 and history[-1] > history[-5]*1.5:
                    break
            niter = len(history)
        meta["numerical_rank_deficient"] = any(
            item["rank"] < item["columns"] for item in meta["X"])
        for step in meta["iterations"]:
            for item in step.get("Z_and_cross", []):
                meta["numerical_rank_deficient"] |= (
                    item["Z_rank"] < L or
                    item["cross_rank"] < min(item["cross_shape"]))
        meta.update(inner_success=ok, seconds=time.perf_counter()-started,
                    final_original_linear_score=score,
                    retained_equations_note="Rank-deficient tall factors use "
                    "their numerically resolved column spaces. Such a solution "
                    "is not asserted unique; acceptance uses the original "
                    "uncompressed moment equation, not the reduced score.")
        self.last_qr_profile = meta
        self.calls += 1
        if not ok:
            self.inner_failures += 1
        return b, ok, None, niter


def initial(m, kind):
    b = m.truth()
    if kind == "mc":
        rng = np.random.default_rng(100200)
        b[m.iw:m.iz] += rng.normal(size=8)*.01
        b[m.iz:m.nl+16] += rng.normal(size=8)*.4
        b[m.delay_slice] = 0.
    return b


def benchmark(periods=(10, 20), init="mc", bands=("L", "H", "W"),
              methods=("M1", "M3"), maxinner=60):
    directory = ROOT/"analysis"
    directory.mkdir(exist_ok=True)
    summary = []
    for P in periods:
        path = ROOT/"mc/data"/f"P{P:02d}_r000.npz"
        with np.load(path) as f:
            data = {k: f[k] for k in f.files}
        for band in bands:
            m = SimModel(BANDS[band][1])
            b0 = initial(m, init)
            for method in methods:
                key = f"P{P:02d}_{band}_{method}_{init}"
                dest = directory/f"qr_profile_v3_{key}.json"
                if dest.exists():
                    summary.append(rd(dest))
                    continue
                result = dict(P=P, band=band, method=method,
                              fixed_nonlinear_initialization=init,
                              data_source=str(path), status="complete")
                aa = select(data, band, method)
                for name, klass in (("normal", ChunkProfileIV),
                                    ("qr", QRProfileIV)):
                    prior = directory/f"qr_profile_{key}.json"
                    if name == "normal" and prior.exists():
                        result[name] = rd(prior)[name]
                        result["old_strict_qr"] = rd(prior)["qr"]
                        if "beta" in result[name]:
                            pr = klass(m, aa, method)
                            result[name]["previous_matmul_linear_score"] = result[name]["original_linear_score"]
                            result[name]["original_linear_score"], _ = linear_moment_check(
                                pr, np.asarray(result[name]["beta"]))
                        continue
                    start = time.perf_counter()
                    try:
                        pr = klass(m, aa, method)
                        b, ok, _, it = pr.profile(
                            b0[m.nl:], b0[:m.nl], maxinner=maxinner)
                        elapsed = time.perf_counter()-start
                        linear_score, ee = linear_moment_check(pr, b)
                        _, _, full_score, _, _ = pr.eq(b)
                        result[name] = dict(
                            inner_success=ok, inner_iterations=it,
                            seconds=elapsed,
                            original_linear_score=linear_score,
                            original_full_moment_score=full_score,
                            residual_squared=ee, beta=b,
                            qr_diagnostics=getattr(pr, "last_qr_profile", None))
                    except Exception as exc:
                        result[name] = dict(
                            inner_success=False, seconds=time.perf_counter()-start,
                            exception=repr(exc), traceback=traceback.format_exc())
                js(dest, result)
                summary.append(result)
                print(key, {k: (v.get("inner_success"),
                      v.get("original_linear_score"), round(v.get("seconds", 0), 2))
                      for k, v in result.items() if k in ("normal", "qr")},
                      flush=True)
    js(directory/f"qr_profile_benchmark_v3_{init}.json", summary)
    lines = ["# 固定非线性参数的 QR 线性 IV 基准", "",
        "只更换线性求解表示，不更改模型、初值、数据或原始矩方程。"
        "完整非线性 score 不要求在固定非线性参数时为零，"
        "内层验收使用原始未压缩线性矩 score。",
        "新版统一内层线性score<1e-9；原M3内部门槛2e-12的结果保留"
        "在qr_profile_benchmark_mc.json及各旧端点中。最终完整375维"
        "score<1e-6和贴边规则不变。",
        "",
        "| P | 窗口 | 方法 | 原内层成功 | 原线性score | QR内层成功 | QR线性score | 原/QR秒 | QR秩亏 |",
        "|---:|:---:|:---:|:---:|---:|:---:|---:|:---:|:---:|"]
    for r in summary:
        a, b = r["normal"], r["qr"]
        lines.append(
            f"| {r['P']} | {r['band']} | {r['method']} | {a['inner_success']} | "
            f"{a.get('original_linear_score', float('nan')):.4g} | "
            f"{b['inner_success']} | {b.get('original_linear_score', float('nan')):.4g} | "
            f"{a['seconds']:.2f}/{b['seconds']:.2f} | "
            f"{b.get('qr_diagnostics', {}).get('numerical_rank_deficient', '未算')} |")
    lines += ["",
        "QR/SVD 的 1e-13 截断作用于缩放后的原始设计矩阵，"
        "不作用于已平方条件数的正规矩阵。秩亏时只求数值可分辨列空间，"
        "报告保留了每个输出的X/Z秩、交叉矩阵秩和未压缩矩方程误差；"
        "不能把通过score解释成物理参数唯一可辨识。"]
    (directory/f"qr_profile_benchmark_v3_{init}.md").write_text(
        "\n".join(lines)+"\n", encoding="utf8")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--periods", default="10,20")
    ap.add_argument("--initial", choices=("mc", "truth"), default="mc")
    ap.add_argument("--bands", default="L,H,W")
    ap.add_argument("--methods", default="M1,M3")
    ap.add_argument("--maxinner", type=int, default=60)
    args = ap.parse_args()
    benchmark(tuple(map(int, args.periods.split(","))), args.initial,
              tuple(args.bands.split(",")), tuple(args.methods.split(",")),
              args.maxinner)
