# Experiment B reproduction

This folder contains only new files. The original `../../code/cl2.py` is imported read-only. Its SHA256 and numerical equivalence checks are in `validation.json`.

Run from this directory with Python 3 and NumPy/SciPy:

```powershell
python run_b.py --validate
python run_b.py --reps 20 --workers 2
python summarize_b.py
python check_tolerance.py
```

The script fixes the experiment environment itself, before importing the original model:

```text
NOFLEX=0 FC=60 AL=30 AH=100 FHI=800 FLO=3 SW=0.03
FSPR=0.01 ZSPR=0.4 NN=300000 NPER=10
```

All BLAS thread counts are set to one. Two independent Monte Carlo records can run concurrently. Every method shares its record and pole initialization with the other two methods. Seeds are `52061019 + repetition_index`.

The primary stopping rule and step cap exactly follow the original IV–Gauss–Newton solver: at most 80 iterations; cap the largest pole-coordinate step at 0.05; stop below 1e-5. The full returned-parameter moment residual is separately evaluated and saved. `check_tolerance.py` performs a stricter 1e-8 re-solve of repetitions 0 and 19 without replacing the primary results.

The 336-parameter Jacobian is sparse by output. `run_b.py` avoids multiplying its explicit zeros by grouping 80 output-local coefficients and 16 common pole coordinates. It preserves the original model, forward finite difference, and normal equations. This speeds computation but is not an additional identification method.

## Raw records

`raw/repNNN_METHOD.npz` contains:

- `beta`, `beta0`, `truth_beta`: all 336 coordinates;
- `error_percent`: 8 frequency errors followed by 8 damping errors, all `100*expm1(delta_log_parameter)`;
- `history`: iteration, uncapped maximum pole step, capped maximum pole step, normalized moment diagnostic, residual sum of squares, iteration seconds;
- `perturbation`: shared initialization perturbation;
- `metadata_json`: seed, method, status, time, iteration count, score, settings and error traceback.

The companion JSON provides metadata without requiring NumPy.

An initial reference-method wrapper had a NumPy 2 batch `solve` RHS shape error before iteration zero. The fixed call supplies an explicit single-column RHS. These development failures are retained in the base files; the same seeded records are successfully rerun as `_retry1` files. `summarize_b.py` selects the retries and explicitly reports the interface errors separately from algorithm convergence failures. It never discards finite outliers. Rerunning the script skips completed records.

`summary.json` gives all-finite and converged-only summaries, min/max, full records and paired differences against reference IV. `report_b.md` contains complete support and flexible-mode tables.
