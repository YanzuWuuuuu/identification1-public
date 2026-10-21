"""Isolated numerical pilot: tall-space IV Schur + bounded trust-region.

The original 375-parameter model, IV moments, physical boxes and final gate
are unchanged.  These endpoints are NEVER merged into mc_qr_final.
"""
from qr_profile import *
from scipy.optimize import least_squares
import hashlib

VERSION = "projected_schur_trf_v2_jac_time_check"


class PilotLimit(RuntimeError):
    pass


def projected_schur(pr, b):
    """Eliminate linear nuisance directions without forming their Gram matrix.

    X_l=Q_x R_x, Z_l=Q_z R_z; C=Q_z.T Q_x.
    Pi=Q_x C^-1 Q_z.T.
    Jbar=(I-Pi)Jn, Zbar=(I-Pi.T)Zn, S=Zbar.T Jbar.
    All factorizations retain the same 1e-13 tall-matrix rank convention.
    """
    f = pr.m.full
    uh = pr.uh(b)
    yh, XL, JN = f.compact(pr.c, b, pr.U, True)
    if pr.method in ("M1", "M2"):
        ZL, ZN = XL, JN
    else:
        _, ZL, ZN = f.compact(pr.c, b, uh, True)
    S = np.zeros((pr.m.nq, pr.m.nq))
    direct_S = np.zeros_like(S)
    rhs = np.zeros(pr.m.nq)
    diagnostics = []
    Jbars, Zbars = [], []
    for o in range(4):
        X, Z = real_rows(XL[:, o]), real_rows(ZL[:, o])
        Jn, Zn = real_rows(JN[:, o]), real_rows(ZN[:, o])
        e = real_rows(pr.Y[:, o]-yh[:, o])
        fx = tall_factor(X)
        fz = fx if pr.method in ("M1", "M2") else tall_factor(Z)
        Qx, Qz = fx["Q"], fz["Q"]
        C = Qz.T@Qx
        c_rhs = Qz.T@Jn
        d_rhs = Qx.T@Zn
        B, _, cr, cs = la.lstsq(C, c_rhs, cond=COND,
                               lapack_driver="gelsd", check_finite=False)
        D, _, dr, _ = la.lstsq(C.T, d_rhs, cond=COND,
                              lapack_driver="gelsd", check_finite=False)
        Jb, Zb = Jn-Qx@B, Zn-Qz@D
        S += Zb.T@Jb
        direct_S += Zn.T@Jb
        rhs += Zb.T@e
        Jbars.append(Jb)
        Zbars.append(Zb)
        diagnostics.append(dict(
            output=o+1, X_rank=fx["rank"], Z_rank=fz["rank"],
            linear_columns=f.L, cross_rank=int(cr),
            cross_transpose_rank=int(dr), cross_singular_values=cs,
            X_condition=fx["scaled_condition"],
            Z_condition=fz["scaled_condition"],
            projected_J_orthogonality=float(
                la.norm(Qz.T@Jb)/max(la.norm(Jn), 1e-300)),
            projected_Z_orthogonality=float(
                la.norm(Qx.T@Zb)/max(la.norm(Zn), 1e-300)),
            full_linear_rank=bool(fx["rank"] == f.L
                                 and fz["rank"] == f.L
                                 and cr == f.L and dr == f.L)))
    Jbar, Zbar = np.vstack(Jbars), np.vstack(Zbars)
    sj = la.norm(Jbar, axis=0)
    sz = la.norm(Zbar, axis=0)
    scaled = S/np.maximum(sz[:, None]*sj[None, :], 1e-300)
    sv = la.svdvals(scaled)
    info = dict(
        formula="Pi=Qx solve(Qz.T Qx,Qz.T); S=Zbar.T Jbar",
        linear_Gram_constructed=False, outputs=diagnostics,
        full_linear_rank=all(x["full_linear_rank"] for x in diagnostics),
        Jbar_norm=sj, Zbar_norm=sz,
        S_scaled_singular_values=sv,
        S_scaled_rank_1e13=int(sum(sv > COND*sv[0])),
        S_scaled_condition=float(sv[0]/max(sv[-1], 1e-300)),
        double_vs_single_projection_relative=float(
            la.norm(S-direct_S)/max(la.norm(direct_S), 1e-300)),
        projected_rhs=rhs)
    return S, info


class StableProfileIV(QRProfileIV):
    """Trust-region root search with an explicitly approximate IV GN Jacobian."""

    def solve_stable(self, b0, output, max_nfev=12, max_profile_calls=24,
                     max_seconds=360., jacobian="projected", fd_step=1e-4):
        directory = Path(output)
        directory.mkdir(parents=True, exist_ok=True)
        start = time.perf_counter()
        nl, nq = self.m.nl, self.m.nq
        q0 = b0[nl:].copy()
        # Coordinate scaling, not a shared clipping multiplier.  The trust
        # region damps weak directions rather than scaling all by their max.
        qscale = np.r_[np.full(2*self.m.nm, .04),
                       np.full(self.m.nd, .4)]
        lb, ub = (self.m.lb-q0)/qscale, (self.m.ub-q0)/qscale
        records, jac_records, cache = [], [], {}
        total_profiles = 0
        last_accepted = None
        denom = None
        halted = None
        optinfo = None

        def budget():
            if total_profiles >= max_profile_calls:
                raise PilotLimit("total_profile_call_limit")
            if time.perf_counter()-start >= max_seconds:
                raise PilotLimit("wall_time_limit")

        def evaluate(x):
            nonlocal total_profiles, denom, last_accepted
            x = np.asarray(x, float)
            key = x.tobytes()
            if key in cache:
                return cache[key]
            budget()
            q = q0+qscale*x
            # Deterministic warm start from last accepted point; all trials,
            # including inner failures, are preserved.  No trial is silently
            # promoted to the final endpoint based on having a smaller score.
            seed = b0[:nl] if last_accepted is None else last_accepted["b"][:nl]
            total_profiles += 1
            b, inner_ok, _, inner_it = self.profile(q, seed)
            A, g, sc, ee, zz = self.eq(b)
            if denom is None:
                denom = np.sqrt(np.maximum(zz[nl:]*ee, 1e-300))
            residual = g[nl:]/denom
            linear_score = float(np.max(abs(g[:nl])/np.sqrt(
                np.maximum(zz[:nl]*ee, 1e-300))))
            record = dict(
                evaluation=total_profiles, x=x, q=q, beta=b,
                inner_ok=bool(inner_ok), inner_iterations=inner_it,
                full_375_score=sc, original_linear_score=linear_score,
                residual_squared=ee, fixed_scaled_moment_norm=la.norm(residual),
                seconds=time.perf_counter()-start,
                qr_numerical_rank_deficient=bool(
                    self.last_qr_profile["numerical_rank_deficient"]))
            js(directory/f"eval_{total_profiles:03d}.json", record)
            records.append(record)
            state = dict(b=b, x=x, fun=residual, inner_ok=inner_ok,
                         A=A, g=g, score=sc, ee=ee, zz=zz,
                         record=record)
            cache[key] = state
            js(directory/"progress.json", dict(
                version=VERSION, complete=False,
                evaluations=total_profiles, last_evaluation=record,
                jacobian=jacobian, actual_profile_call_limit=max_profile_calls,
                wall_time_limit_seconds=max_seconds))
            print("profile", total_profiles, "score", f"{sc:.4g}",
                  "inner", inner_ok, "seconds",
                  round(time.perf_counter()-start, 2), flush=True)
            return state

        def fun(x):
            r = evaluate(x)
            # Invalid inner profiles are not accepted as exact reductions.
            # An exception preserves the failed evaluation instead of inventing
            # a penalty function that could change the root.
            if not r["inner_ok"]:
                raise PilotLimit("inner_linear_moment_not_converged")
            return r["fun"]

        def jac(x):
            nonlocal last_accepted
            r = evaluate(x)
            last_accepted = r
            # LAPACK calls cannot be interrupted safely; stop at the next
            # boundary without starting another expensive projection.
            if time.perf_counter()-start >= max_seconds:
                raise PilotLimit("wall_time_limit")
            S, info = projected_schur(self, r["b"])
            info.update(jacobian_number=len(jac_records)+1,
                        evaluation=r["record"]["evaluation"],
                        full_score=r["score"],
                        frozen_instrument_GN=(jacobian == "projected"))
            js(directory/f"jac_{len(jac_records)+1:03d}.json", info)
            jac_records.append(info)
            if not info["full_linear_rank"]:
                raise PilotLimit("tall_linear_or_cross_rank_deficient")
            if time.perf_counter()-start >= max_seconds:
                raise PilotLimit("wall_time_limit")
            if jacobian == "projected":
                # This is the frozen-instrument GN derivative, not the full
                # derivative of Z(beta)^T e(beta).  The actual residual always
                # updates the instrument; convergence is checked with all
                # original moments.
                return -S*qscale[None, :]/denom[:, None]
            J = np.empty((nq, nq))
            for j in range(nq):
                h = fd_step*max(1., abs(x[j]))
                h = min(h, .49*(ub[j]-x[j]), .49*(x[j]-lb[j]))
                if h <= 1e-12:
                    raise PilotLimit("central_difference_near_boundary")
                xp, xm = x.copy(), x.copy()
                xp[j] += h
                xm[j] -= h
                J[:, j] = (fun(xp)-fun(xm))/(2*h)
            return J

        try:
            opt = least_squares(
                fun, np.zeros(nq), jac=jac, bounds=(lb, ub), method="trf",
                tr_solver="exact", x_scale=np.ones(nq),
                max_nfev=max_nfev, ftol=1e-12, xtol=1e-10, gtol=1e-12)
            end = evaluate(opt.x)
            optinfo = dict(status=int(opt.status), message=str(opt.message),
                           scipy_success=bool(opt.success), nfev=opt.nfev,
                           njev=opt.njev, optimality=float(opt.optimality),
                           cost=float(opt.cost))
            halted = "trf_finished"
        except PilotLimit as exc:
            halted = str(exc)
            end = last_accepted
            if end is None and records:
                end = cache[np.asarray(records[0]["x"]).tobytes()]
        except Exception as exc:
            halted = "exception"
            optinfo = dict(exception=repr(exc), traceback=traceback.format_exc())
            end = last_accepted
        if end is None:
            return b0, dict(success=False, status=halted,
                            optimizer=optinfo, seconds=time.perf_counter()-start,
                            profile_calls=total_profiles)
        b = end["b"]
        # Independent final call to the original, uncompressed 375 equations.
        _, g, score, ee, zz = self.eq(b)
        linear_score = float(np.max(abs(g[:nl])/np.sqrt(
            np.maximum(zz[:nl]*ee, 1e-300))))
        near = bool(np.min(np.minimum(
            b[nl:]-self.m.lb, self.m.ub-b[nl:])) < 1e-5)
        good = bool(np.isfinite(b).all() and score < 1e-6
                    and linear_score < self.inner_tol and not near)
        info = dict(
            success=good, status=halted, score=score,
            original_linear_score=linear_score, near_bound=near,
            seconds=time.perf_counter()-start, profile_calls=total_profiles,
            iterations=len(jac_records), optimizer=optinfo,
            version=VERSION, jacobian=jacobian,
            frozen_instrument_GN=(jacobian == "projected"),
            final_full_375_threshold=1e-6, linear_threshold=self.inner_tol,
            rcond_on_tall_matrix=COND,
            qscale=qscale, final_evaluation=end["record"]["evaluation"],
            history=[{k: r[k] for k in (
                "evaluation", "full_375_score", "original_linear_score",
                "inner_ok", "seconds")} for r in records],
            all_trials_preserved=True, formal_MC_inclusion=False)
        info["wall_limit_note"] = (
            "Cooperative check before/after each profile/projection; an in-flight "
            "LAPACK call may exceed the requested seconds, but no further "
            "expensive projection is started after expiry.")
        js(directory/"progress.json", dict(complete=True, solver=info))
        return b, info


def run_pilot(P, band, method, max_nfev=12, max_profile_calls=24,
              max_seconds=360., jacobian="projected", tag="v1"):
    dest = ROOT/"stable_profile_pilot"/f"{tag}_P{P:02d}_r000_{band}_{method}"
    dest.mkdir(parents=True, exist_ok=True)
    if (dest/"result.json").exists():
        raise FileExistsError(f"Pilot already exists: {dest}")
    data_path = ROOT/"mc/data"/f"P{P:02d}_r000.npz"
    with np.load(data_path) as fh:
        data = {k: fh[k] for k in fh.files}
    m = SimModel(BANDS[band][1])
    pr = StableProfileIV(m, select(data, band, method), method)
    b0 = initial(m, "mc")
    meta = dict(
        P=P, band=band, method=method, replicate=0, data_source=data_path,
        initial_beta=b0, version=VERSION, jacobian=jacobian,
        max_nfev=max_nfev, max_actual_profile_calls=max_profile_calls,
        wall_time_limit_seconds=max_seconds, formal_MC_inclusion=False,
        source_hashes={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), ROOT/"code/qr_profile.py",
                      ROOT/"code/sim_core.py")})
    js(dest/"protocol.json", meta)
    try:
        b, info = pr.solve_stable(
            b0, dest, max_nfev, max_profile_calls, max_seconds, jacobian)
        meta.update(status="complete", endpoint=pack(b, info, m))
    except Exception as exc:
        meta.update(status="exception", exception=repr(exc),
                    traceback=traceback.format_exc())
    js(dest/"result.json", meta)
    print("PILOT", P, band, method,
          meta.get("endpoint", {}).get("solver"), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--periods", type=int, required=True)
    ap.add_argument("--band", choices=tuple(BANDS), required=True)
    ap.add_argument("--method", choices=("M1", "M3"), required=True)
    ap.add_argument("--max-nfev", type=int, default=12)
    ap.add_argument("--max-profile-calls", type=int, default=24)
    ap.add_argument("--max-seconds", type=float, default=360.)
    ap.add_argument("--jacobian", choices=("projected", "central"),
                    default="projected")
    ap.add_argument("--tag", default="v1")
    ar = ap.parse_args()
    run_pilot(ar.periods, ar.band, ar.method, ar.max_nfev,
              ar.max_profile_calls, ar.max_seconds, ar.jacobian, ar.tag)
