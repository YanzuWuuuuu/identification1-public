# Run contract — 2026-10-19

Expand existing closed-loop frequency-domain additive IV experiments without modifying the baseline, theory, original data or prior results. Main question: how period count, EIV covariance estimation and indirect multi-experiment FRFs affect modal bias, dispersion and practical availability.

## Locked protocol
- A: cl2, NOFLEX=1, FC=60, AL=30, AH=100, SW=.03, 5000 Hz, 30000 samples/period, band 3–100 Hz. P=1,2,4,8: 40 paired trials each; P=16:20. M1–M4 same data and initial modal perturbations. Seed=202610190+1000*P+r. FSPR=.01 ZSPR=.4. M5 uses same per-period DFTs.
- A': four single-input experiments, P/4 recorded periods each, P=4,8,16; at least two settling periods each. M6 covariance weight and M6u uniform. Count settling overhead separately.
- B: full 8-mode model, band 3–800 Hz, 10 periods, 20 trials, M1/M2/M3.
- C: four specified real records, second period onward, P=1,2,4,8, 3–25 Hz. CT versus exact ZOH support blocks with requested 100 us delay adjustment, M1/M2/M3, SML where definable. P=8 split/block estimates.
- Metrics: physical relative errors 100*expm1(log-error), mean±SE, SD, RMSE, all failures and iteration counts. No truth claims for real records. Numerical convergence is distinguished from statistical consistency.
- Each trial/method persisted immediately; complete parameters retained. No cherry-picked reruns. All finite endpoints and converged-only summaries distinguished.

## Comparator validity checks
- Sample covariance output rank <= P-1, so raw 4-output SML not invertible at P=2 or 4; report unavailable and distinguish regularized diagnostics.
- Parameter-dependent covariance requires appropriate objective derivatives. Literature verification determines SML implementation.
- M6 has one period per axis at P=4: period covariance unavailable. M6u remains definable; distinguish any covariance regularization/pooling.
- P=1 all DFT bins are excited, so M1 and M2 must agree exactly for matched starts.
- Reference with K=0 is exogenous but useful only with sufficient instrument relevance. Finite-sample unbiasedness is not assumed.

## Implementation map
- Root: experiment_a/ runner using unchanged cl2; shared DFT checkpoints; aggregate report and reproducibility package.
- Comparator agent: comparators/ SML and indirect FRF algorithms, validation, runs and citations.
- B agent: experiment_b/ full-model runs.
- C agent: experiment_c/ copied-and-adapted real fitting code.

## Smoke
Requested 3-trial command completed. Mean damping relative errors: direct [51.427924,17.535260,5.347962]%, reference [2.277899,-3.960363,1.043120]%; iterations 8–17. Baseline scale accepted.

## Evidence decision
Classify each user-provided expectation as supported, unsupported, or inconclusive from all recorded trials. Report rank/optimizer failures rather than imposing expected trends. No theoretical rewrite in this run.

## Final evidence and route
All requested run counts and positions completed. A:720; B:60 final; C:120 cases including covariance-unavailable entries; SML:140 main; indirect:160 canonicalIVGN plus160 optimizer checks. Score and rank failures are disclosed. Original source/data hashes unchanged. No-noise and independent algebra checks passed.

Conclusion: referenceIV reduces long-record support-mode bias, but P1/P2 reliability is not established and SML is more efficient in this comparison at P8/P16. Indirect FRF under the given budget has rank/numerical/ratio difficulties. Flexible bias is much smaller but not uniformly zero. ZOH real-data correction is small at P8.

Next route: use the measured report to narrow manuscript claims; no additional algorithm redesign or theory modification authorized by this experiment run has been performed.
