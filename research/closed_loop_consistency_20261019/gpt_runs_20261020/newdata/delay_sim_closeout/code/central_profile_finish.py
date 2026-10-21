"""Separate short complete-profile Jacobian finish; retains prior pilot."""
from stable_profile import *


def main():
    started = time.perf_counter()
    source = ROOT/"stable_profile_pilot/v2_P10_r000_H_M3/result.json"
    dest = ROOT/"stable_profile_pilot/central_P10_r000_H_M3"
    dest.mkdir(parents=True, exist_ok=True)
    if (dest/"result.json").exists():
        raise FileExistsError(dest)
    previous = rd(source)
    b0 = np.asarray(previous["endpoint"]["beta"])
    data_path = ROOT/"mc/data/P10_r000.npz"
    with np.load(data_path) as fh:
        data = {k: fh[k] for k in fh.files}
    model = SimModel(BANDS["H"][1])
    pr = StableProfileIV(model, select(data, "H", "M3"), "M3")
    qscale = np.r_[np.full(16, .04), np.full(7, .4)]
    _, _, _, ee, zz = pr.eq(b0)
    denom = np.sqrt(np.maximum(zz[model.nl:]*ee, 1e-300))
    directions = {}
    x = np.zeros(23)
    x[4] = 1.
    directions["frequency_345Hz"] = x
    x = np.zeros(23)
    x[16] = 1.
    directions["common_delay"] = x
    x = np.random.default_rng(20261021).normal(size=23)
    directions["mixed"] = x/la.norm(x)
    diagnostic = []
    calls = 0
    for name, direction in directions.items():
        slopes = []
        for h in (1e-3, 3e-3):
            values = []
            for sign in (-1, 1):
                q = b0[model.nl:]+sign*h*qscale*direction
                b, ok, _, it = pr.profile(q, b0[:model.nl])
                calls += 1
                _, g, sc, _, _ = pr.eq(b)
                record = dict(direction=name, h_x=h, sign=sign, beta=b,
                              inner_ok=ok, inner_iterations=it,
                              full_score=sc, moment=g[model.nl:]/denom)
                js(dest/f"fdcheck_{calls:02d}.json", record)
                values.append(record["moment"])
            slopes.append((values[1]-values[0])/(2*h))
        item = dict(
            direction=name, steps_x=[1e-3, 3e-3], derivatives=slopes,
            relative_difference=float(la.norm(slopes[0]-slopes[1])
                                      /max(la.norm(slopes[1]), 1e-300)),
            norm_1e3=float(la.norm(slopes[0])),
            norm_3e3=float(la.norm(slopes[1])))
        diagnostic.append(item)
        js(dest/"difference_step_check.json", diagnostic)
        print("FD_CHECK", name, item["relative_difference"], flush=True)
    protocol = dict(
        source_endpoint=str(source), source_score=previous["endpoint"]["solver"]["score"],
        starting_beta=b0, method="M3", band="H", P=10,
        jacobian="complete central difference of profiled IV moments",
        difference_step_x=1e-3, difference_diagnostic_calls=calls,
        max_solver_profile_calls=150, max_total_profile_calls=150+calls,
        total_wall_time_limit_seconds=480, formal_MC_inclusion=False,
        source_hashes={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), ROOT/"code/stable_profile.py",
                      ROOT/"code/qr_profile.py", ROOT/"code/sim_core.py")})
    js(dest/"protocol.json", protocol)
    try:
        b, info = pr.solve_stable(
            b0, dest, max_nfev=8, max_profile_calls=150,
            max_seconds=max(1., 480-(time.perf_counter()-started)),
            jacobian="central", fd_step=1e-3)
        info["difference_diagnostic_calls"] = calls
        info["total_profile_calls_including_diagnostics"] = calls+info["profile_calls"]
        info["total_seconds_including_diagnostics"] = time.perf_counter()-started
        result = dict(protocol=protocol, endpoint=pack(b, info, model),
                      status="complete", derivative_check=diagnostic)
    except Exception as exc:
        result = dict(protocol=protocol, status="exception",
                      exception=repr(exc), traceback=traceback.format_exc(),
                      derivative_check=diagnostic)
    js(dest/"result.json", result)
    print("CENTRAL_FINISH", result.get("endpoint", {}).get("solver"), flush=True)


if __name__ == "__main__":
    main()
