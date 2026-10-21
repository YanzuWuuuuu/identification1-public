"""Uniform QR/SVD rerun, separate from all preserved normal-equation attempts."""
from sim_core import *
from run_mc import dataset
from qr_profile import QRProfileIV
import argparse
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed


def run(P, r):
    targets = [("L", "M3"), ("H", "M3")]
    if P == 20:
        targets += [("L", "M1"), ("H", "M1"), ("W", "M1")]
    paths = {t: ROOT/"mc_qr_final"/"raw"/f"P{P:02d}_r{r:03d}_{t[0]}_{t[1]}.json" for t in targets}
    if all(p.exists() for p in paths.values()):
        return
    a = dataset(P, r)
    for band, method in targets:
        dest = paths[band, method]
        if dest.exists():
            continue
        start = time.perf_counter()
        out = dict(P=P, replicate=r, seed=1002+r, band=band, method=method,
                   parameter_count=375, modal_count=8, status="exception",
                   solver_route="uniform tall QR/SVD nuisance projection; complete profiled moment hybrid",
                   initialization="same frequency/damping perturbation seed as delay_check; zero delay offsets",
                   common_grid=False, linear_moment_tolerance=1e-9,
                   final_complete_moment_tolerance=1e-6, tall_rcond=1e-13)
        try:
            m = SimModel(BANDS[band][1])
            pr = QRProfileIV(m, select(a, band, method), method)
            b0 = m.truth()
            rng = np.random.default_rng(100200+r)
            b0[m.iw:m.iz] += rng.normal(size=8)*.01
            b0[m.iz:m.nl+16] += rng.normal(size=8)*.4
            b0[m.delay_slice] = 0.
            out["initial_beta"] = b0
            js(dest.with_suffix(".started.json"), out)
            b, info = pr.solve(b0, maxouter=70, tol=2e-8, hybrid=True)
            # The acceptance check uses the original 375-dimensional equation,
            # not the QR coordinates or a selected subset of its moments.
            _, _, original_score, _, _ = pr.eq(b)
            info["original_complete_score"] = original_score
            info["success"] = bool(info["success"] and original_score < 1e-6)
            qmeta = getattr(pr, "last_qr_profile", {})
            out.update(endpoint=pack(b, info, m), status="complete",
                       success=info["success"], frequency_points=len(pr.U),
                       linear_profile_diagnostics=qmeta,
                       complete_moment_root_attempted=info.get("hybrid") is not None,
                       phase_branch_error_radians=float(2*np.pi*BANDS[band][1]*np.max(abs(m.delay(b)-TRUTH_TAU))),
                       complete_moment_root_note="profile hybrid differentiates updated instruments and re-estimated nuisances; unavailable if the inner original linear moments cannot be solved to 1e-9")
        except Exception as exc:
            out.update(exception=repr(exc), traceback=traceback.format_exc(), success=False)
        out["seconds"] = time.perf_counter()-start
        js(dest, out)
        print(P, r, band, method, out["success"],
              out.get("endpoint", {}).get("solver", {}).get("status"),
              out.get("endpoint", {}).get("solver", {}).get("score"),
              round(out["seconds"], 2), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--periods", default="10,20,40")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=40)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(run, P, r) for P in map(int, args.periods.split(","))
                for r in range(args.start, args.start+args.count)]
        for job in as_completed(jobs):
            job.result()
    js(ROOT/"mc_qr_final"/"batch_complete.json",
       dict(periods=args.periods, start=args.start, count=args.count,
            final_score=1e-6, linear_score=1e-9, tall_rcond=1e-13))
