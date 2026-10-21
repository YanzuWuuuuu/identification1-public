"""Complete 375-dimensional moment root when profiled inner solve fails.

Original attempts remain untouched. Finite-difference Newton/hybrid differentiates
the instruments as well as the output residual. No modal directions are removed.
"""
from sim_core import *
from run_mc import dataset
from root_solver import exact_root
import argparse
import traceback
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED


def refine(path, maxfev=900):
    source = rd(path)
    dest = ROOT/"mc"/"refinements"/path.name
    if dest.exists(): return rd(dest)
    start = time.perf_counter()
    result = dict(source_file=str(path), P=source["P"], replicate=source["replicate"],
                  seed=source["seed"], band=source["band"], method=source["method"],
                  original_success=source.get("success", False), status="exception")
    try:
        P, r, band, method = (source[k] for k in ("P", "replicate", "band", "method"))
        data = select(dataset(P, r), band, method)
        m = SimModel(BANDS[band][1])
        pr = ChunkProfileIV(m, data, method)
        b0 = np.asarray(source["endpoint"]["beta"])
        b, info = exact_root(m.full, pr.c, pr.U, pr.Y, pr.uh, b0, maxfev=maxfev)
        q = b[m.nl:]
        near = bool(np.any(q <= m.lb+1e-5) or np.any(q >= m.ub-1e-5))
        info.update(near_bound=near, success=bool(info["success"] and not near))
        result.update(endpoint=pack(b, info, m), success=info["success"], status="complete",
                      complete_moment_root_attempted=True,
                      frequency_points=len(pr.U))
    except Exception as exc:
        result.update(success=False, exception=repr(exc), traceback=traceback.format_exc())
    result["seconds"] = time.perf_counter()-start
    js(dest, result)
    print(path.name, result.get("success"), result.get("endpoint", {}).get("solver", {}).get("score"),
          round(result["seconds"], 2), flush=True)
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--file")
    ap.add_argument("--maxfev", type=int, default=900)
    ap.add_argument("--watch-all", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if args.file:
        refine(Path(args.file), args.maxfev)
    elif args.watch_all:
        expected = []
        for P in (10, 20, 40):
            for r in range(40):
                for band in ("L", "H"):
                    expected.append(ROOT/"mc/raw"/f"P{P:02d}_r{r:03d}_{band}_M3.json")
        submitted = set()
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            active = {}
            while len(submitted) < len(expected) or active:
                for path in expected:
                    if path in submitted or not path.exists(): continue
                    dest = ROOT/"mc/refinements"/path.name
                    source = rd(path)
                    if dest.exists() or source.get("success"):
                        submitted.add(path)
                        continue
                    if len(active) >= args.workers: break
                    submitted.add(path)
                    active[pool.submit(refine, path, args.maxfev)] = path
                if active:
                    done, _ = wait(active, timeout=5, return_when=FIRST_COMPLETED)
                    for task in done:
                        task.result()
                        del active[task]
                else:
                    time.sleep(5)
        js(ROOT/"mc/full_root_batch_complete.json",
           dict(total_reference_endpoints=len(expected), maxfev=args.maxfev))
    else:
        ap.error("supply --file or --watch-all")
