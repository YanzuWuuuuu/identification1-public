# Closed-loop validation, 2026-10-14

Status: completed scoped validation and manuscript integration. Original stage-0
and earlier experiments were not changed. Five noisy cases and one noiseless check
are saved. All 56 structured starts are feasible and stable; 42 pass the strict
solver check, 14 do not. The separate raw-score gate fails for every start, so no
normality / coverage claim is made. Built-in LaTeX compilation was attempted after
editing but failed before compilation due to the Windows sandbox helper; source
checks do not replace compilation. See outputs_structured/source_verification.json.

Second amendment, supported by the origin audit: the independent 32-coefficient
family produced a seventh integrator (one pole at z=1) after known rank-three
feedback. For the representative first fit, n_y' H_-1 n_u=-5.2734e-8, whereas
the physical truth gives 5.0822e-21. This is a model-family issue, not a c2d error.
Preserve the physically known additive a_jk=a_j1+a_1k-a_11 structure (7 linear
coordinates); only b_jk remain independent (16). Total free parameters drop from
134 to 125. The larger-family identifiability proof restricts to this smooth
subfamily; repeated roots remain allowed. Outputs for the original 56 fits remain
untouched in outputs/. The same configurations/seeds/budgets are rerun into
outputs_structured/. The first 125-parameter pilot meets constrained KKT tolerance
after 307 iterations and its fitted closed-loop spectral radius is 0.9964.
This change was chosen from the analytic low-frequency structure, before inspecting
heldout performance. Global-score / interior normality is still a separate gate.

2026-10-14 implementation amendment: analytic constraint derivatives fixed the
initial finite-difference evaluation-budget overhead. On the identical noiseless
pilot and 400-iteration budget, generic SQP reached cost 6.61e-6 in 85.1 s;
constrained least-squares interior-point reached 6.25e-7 in 28.2 s. Neither passed
the strict score tolerance. Main runs now lock `lsqnonlin` interior-point with
the SAME residual, analytic derivatives, constraints, bounds and 400-iteration
budget. This is a numerical solver change, not a different statistical estimator.
Retain both failed-to-converge noiseless runs; never relabel them as convergence.
Run 16 starts for distinct N=32, nested N=8, and repeated-filter N=32; two additional
noise realizations use four starts and are diagnostic rather than equal-budget
solver-success comparisons. No large-sample rate or coverage claim is made.

Question: can full-output fitting against exogenous actuator excitation recover
five-mode, D=0 mechanics and relative channel dynamics under known digital feedback?
Null: the implementation fails the matched-model fit or its calibrated maps do not
improve independent validation over static calibration. Alternative: it recovers
task targets and improves the relevant rigid leakage without erasing flexible modes.

## Locked scope
- Read numeric matrices from matlab/0/Multiple_main_Control.m; all five modes, D=0.
- Known rank-three stabilizing controller; excitation after force allocation.
- Exact ZOH and fractional-delay frequency predictor (matrix exponentials, no Pade).
- Independently generated sample-domain closed-loop records, synchronous d/u/y.
- Pilot data determine mechanical initialization and chart; no truth in estimator.
- Expanded denominator coefficients, rank-three M, rank-one complex residues.
- Constrained least-squares with analytic derivatives, fixed weights and logged failures.
- Nominal, estimated static-only, estimated dynamic-plus-static, and oracle maps.
- Separate output-only rigid-motion leakage from double-null transfer leakage.

## Minimal code map
New files only in this directory: source loader and known controller, exact sampled
predictor with derivatives, data generator, data-only initializer/chart, constrained
fit, independent metrics, run entry point and JSON/MAT/log outputs.

## Gates and bounded runs
1. Verify predictor against independent c2d and derivative finite differences.
   Verify digital closed-loop poles, d/u/y equation and noise correlation.
2. Noise-free matched-model run: data-derived pilot and joint estimator.
3. Main noisy distinct-channel run: 16 starts, 400 iterations / 4000 evaluations each.
4. Repeated-noise and common-filter controls, with the same fixed configuration;
   smaller diagnostic runs must be explicitly marked if computation requires it.
5. Evaluate on withheld frequencies and fresh records; do not select by test error.

No significance claim or consistency-rate claim from a handful of realizations.
Report per-run failures, parameter/relative-response/subspace errors, residuals,
leakage, retained target response, noise and identifiability conditioning.
Pilot-only success cannot close the request. Failures are valid results; change a
route only with a recorded reason, never tune on hidden truth or omit failed seeds.

## Outputs / stop condition
Save config, seeds, environment, raw-data examples, all start outcomes, chosen fits,
metrics and summary. Update the existing TeX with measured conclusions and limits.
Finish when implementation gates and a substantive noisy comparison are available,
or when a reproducible model/solver failure provides a concrete negative conclusion.
This local directory is the isolated experiment surface; external artifact tools
referenced by the experiment skill are unavailable in this session.
