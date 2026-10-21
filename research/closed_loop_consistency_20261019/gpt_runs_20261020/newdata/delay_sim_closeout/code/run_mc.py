"""Checkpointed requested MC. Failures and every solver endpoint are retained."""
from sim_core import *
import argparse
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed


def dataset(P, r):
    path = ROOT/"mc"/"data"/f"P{P:02d}_r{r:03d}.npz"
    if path.exists():
        with np.load(path) as f:
            return {k: f[k] for k in f.files}
    a = generate(P, 1002+r)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".tmp.npz")
    np.savez_compressed(temp, **a)
    temp.replace(path)
    return a


def run(P, r, methods=None):
    targets = [("L", "M3"), ("H", "M3")]
    if P == 20: targets += [("L", "M1"), ("H", "M1"), ("W", "M1")]
    if methods:
        targets = [x for x in targets if f"{x[0]}_{x[1]}" in methods]
    paths = {t: ROOT/"mc"/"raw"/f"P{P:02d}_r{r:03d}_{t[0]}_{t[1]}.json" for t in targets}
    if all(p.exists() for p in paths.values()):
        return
    a = dataset(P, r)
    for band, method in targets:
        dest = paths[band, method]
        if dest.exists(): continue
        start = time.perf_counter()
        out = dict(P=P, replicate=r, seed=1002+r, band=band, method=method,
                   parameter_count=375, modal_count=8, status="exception",
                   initialization="same frequency/damping perturbation seed as delay_check; zero delay offsets",
                   common_grid=False,
                   common_grid_note="known local simulation branch: initial 200us, true common 200us; no truth relative delays supplied")
        try:
            m = SimModel(BANDS[band][1])
            pr = ChunkProfileIV(m, select(a, band, method), method)
            b0 = m.truth()
            rng = np.random.default_rng(100200+r)
            b0[m.iw:m.iz] += rng.normal(size=8)*.01
            b0[m.iz:m.nl+16] += rng.normal(size=8)*.4
            b0[m.delay_slice] = 0.
            out["initial_beta"] = b0
            js(dest.with_suffix(".started.json"), out)
            b, info = pr.solve(b0, maxouter=70, tol=2e-8, hybrid=True)
            out.update(endpoint=pack(b, info, m), status="complete",
                       success=info["success"], frequency_points=len(pr.U),
                       complete_moment_root_attempted=info.get("hybrid") is not None,
                       complete_moment_root_note="profile hybrid differentiates updated instruments and re-estimated nuisances; unavailable if the inner linear moment equation has no converged root")
        except Exception as exc:
            out.update(exception=repr(exc), traceback=traceback.format_exc(), success=False)
        out["seconds"] = time.perf_counter()-start
        js(dest, out)
        print(P, r, band, method, out["success"],
              out.get("endpoint", {}).get("solver", {}).get("status"), round(out["seconds"], 2), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--periods", default="10,20,40")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=40)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--methods", default="")
    args = ap.parse_args()
    selected = args.methods.split(",") if args.methods else None
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(run, P, r, selected) for P in map(int, args.periods.split(","))
                for r in range(args.start, args.start+args.count)]
        for job in as_completed(jobs):
            job.result()
