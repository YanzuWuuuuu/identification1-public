# Follow-up D/E reproduction

All code and outputs added by this campaign are in this directory. The parent simulation, measured records, prior results and theory document are read-only. Run from this folder with Python 3.11 + numpy/scipy. Scripts set one BLAS thread before importing numpy. No MATLAB dependency or model-structure change.

## Read first

- `后续实验D_E_报告.md`: integrated report; task0 first, all simulation statistics, scale and covariance checks, real data and half-window differences, E, fixed decision rules, deviations, all-finite failure-inclusive tables.
- `PROTOCOL.md`: decisions frozen before follow-up results.
- `decision_rules.json`: exact per-parameter rule evaluation, both converged-only and all-finite.
- `package_validation.json`: expected counts; `original_integrity_check.json`: original-file checksums.
- `audit_math.md`, `audit_e.md`: independent read-only implementation audits and limitations.

## Runtime

Used local Python `C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe`, 16 logical processors, approximately32 GB RAM, one BLAS thread per worker. Peak planned concurrency: task0 2 + D4 + E2 + real2 workers, with task0 finishing before D's bulk run. Every finished estimator writes JSON atomically; every newly generated simulation writes NPZ. Reruns skip existing trial results. To intentionally rerun a comparison without overwriting evidence, use another isolated output directory/code copy.

```powershell
$py = 'C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe'
Set-Location -LiteralPath 'C:/Users/Administrator/Desktop/my_article/identifiation1/research/closed_loop_consistency_20261019/gpt_runs_20261019/followup_DE'
```

## Task 0

```powershell
& $py task0/run_task0.py audit
& $py task0/run_task0.py long --workers 2 --count 40
& $py task0/run_task0.py population
& $py task0/run_task0.py population_analytic
& $py task0/check_initialization.py
& $py task0/build_task0_report.py
```

`task0/long_data/` and `long_results/` contain40long-P1 records, P0=240000, originalP1 seeds. `population/` retains both original finite-difference and analytic Frechet-derivative diagnostics; `initialization_results/` retains additional short-P1 checks including one failure at the tighter diagnostic tolerance. See task0's own report for exact CLI details and caveats.

## D simulation and shared functions

```powershell
& $py sim_d_run.py --gate
& $py sim_d_validate.py
& $py sim_d_run.py --periods 2,4,8,16 --count 40 --workers 4
& $py sim_d_run.py --summary
& $py sim_d_analyze.py
```

- `sim_d_run.py`: reuse140existingpairedrecords, generate20additionalP16records; fit10variants perrecord; calculateconditionalfullmatrixoraclecovariance.
- `weighted_core.py`: fullmatrixprecision/whitening, originalIV-GN stoppingrule, nearestnonexcitedbincovariance.
- `sml_core_copy.py`: copiedoldprofileSML, plus externallysuppliedcovarianceadapter. Includes model-dependentcovariancederivative.
- `sim_d_validate.py`: scale/whitening/localprojectiondiagnostics. Different G at neighboring frequencies is explicitly distinguished from fixed-central-G covariance projection.
- `sim_d_analyze.py`: complete summary, paired M5/M5n checks and fixed-rule evaluation.
- `sim_d/results/`:1600JSONendpoints, all40replicates of allmethods, includingfailures.
- `sim_d/diagnostics/`:160oraclecovarianceNPZ, localwindowcoverage, scaleandmatrixchecks.
- `sim_d/data/`:20newP16NPZonly. OlderAinputNPZfilesareloadedread-only; theirpaths andhashes are retained in `original_inventory.json`.
- `sim_d/run.log`, `sim_d/summary.json`, `sim_d/oracle_theory_summary.json`, `sim_d/RESULTS_D.md`.

Primary variant is `M3w_n16_o1` versus `M5n_n16`. Never select `n8` or `o2` after inspecting results to change the main conclusion. SML nfev is not comparable to IV outeriterations as a computational-efficiency claim.

## Real D

```powershell
& $py real_d/run_real_d.py --all --workers 2
& $py real_d/summarize_real_d.py
```

See `real_d/README.md` for command flags. `raw/` contains140unique method-window records:116actualfits and24explicitly unavailable covariance fits on one-period halves. The report expands repeated identical windows into180requested contexts. Fourpositions -.5,0,.5,1; P2,4,8; firstwindowstarts atsecondrecordedperiod. ExactZOHmodel, fixed100usdelaycorrection, D0/D2 andmodaltransients unchangedfrompriorC. No ground-truth accuracy ranking.

## E short-record robustness

```powershell
& $py validate_e.py
& $py run_e.py --periods 1,2 --count 40 --workers 2
```

- `run_e.py`: E1originalinitialization, E2directM2initialization, E3dampedIVmomentroot; M3P1/P2 andM3wP2, 360formalrecords.
- `iv_lm.py`: momentobjective andfullmomentJacobian, including changingreferenceinstrument andmodelHessiancontraction. Outerfinite-difference step1e-4.
- `validation_e.json`: complexmatrixweight parity, actualmomentderivative check, priorAparity.
- `short_e/`: allformalendpoints andlog. `pilot_quasi_jacobian/` retains one preliminary failedsimplified-Jacobian trial; it is not included in the40-replicate tables.
- The M3w E comparison applies eachstrategy to both preliminaryM3 andweightedIV; weights maythereforediffer, unlike a fixed-weight solver-only ablation.

## Rebuild and verify

```powershell
& $py build_report.py
& $py inventory.py --verify
```

Do not rerun `inventory.py` without `--verify` after modifications to originals: it creates the baselineinventory and must not be used to silently reset thatbaseline. All sourcehashdifferences, ifany, mustbe investigated. Statistics preserve nonfinite/failed records, use expm1percentageerrors, and publish both converged-only andall-finite subsets. Primary E3's unfavorable outcomes remain partofthedeliverable.
