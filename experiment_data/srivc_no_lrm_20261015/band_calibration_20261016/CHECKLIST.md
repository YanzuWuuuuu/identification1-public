# Completion and reproducibility record

- [x] Primary band 25–80 Hz and diagnostic bounds/sensitivity runs recorded before execution.
- [x] Previous four normal-position first/second-half LRM arrays reused, anomalous -1 mm excluded.
- [x] Old raw data, fit outputs and model scripts unchanged.
- [x] Real null-vector eigenprofile derived separately for output and input sides.
- [x] Analytic delay gradient finite-difference errors below 8e-8 in all eight primary side fits.
- [x] Three prescribed starts; training-only selection; all main side fits positive exit, KKT <1e-10, no boundary hit.
- [x] Nominal-static training bound and zero-delay feasible cost checks pass.
- [x] Both cross blocks and full-port residuals saved; N-to-N is not treated as the primary suppression target.
- [x] Same directions and scaling for D1/D2 delay ablation.
- [x] Old B4 results independently reproduced: maximum relative RMS discrepancy 1.26e-15.
- [x] Causal alignment sum check: maximum error 1.53e-20 seconds.
- [x] Unitary compensation operator-norm change <=4.45e-16; no claim of measured-noise PSD equality.
- [x] z0 matrices transferred without refitting to three other positions.
- [x] Prescribed training-window and delay-bound sensitivities reported, no validation winner substituted.
- [x] First-plate-band total response included without calling it true modal retention.
- [x] Nine figures saved in PNG/vector-PDF, figure manifest saved; key figures visually inspected.
- [x] Existing manuscript edited in place; previous source backed up under docs/research/history.
- [x] Built-in compiler attempted; Windows process-launch access denied. Previously authorized local Tectonic compilation succeeded.
- [x] Local manuscript PDF has 27 pages; added experimental pages visually inspected; no undefined-reference/overfull/missing-character warnings found.

Commands in MATLAB from workspace root:

```matlab
addpath(fullfile(pwd,'experiment_data','srivc_no_lrm_20261015'));
campaign = run_band_calibration_diagnostic;
files = plot_band_calibration_diagnostic;
```

The run has no random draws. Optimization starts and all numerical settings are in the source and saved MAT configuration. MATLAB version is recorded in the campaign. The separate independent_numeric_checks.mat contains the three comparison checks above. z0_band_static_matrices.m is a labeled export of the new D0 diagnostic matrices; it does not replace the previous model-based matrices.

Evidence status: supported for lower measured cross responses in this band and this recorded dataset; inconclusive for true hardware delay recovery, a substantial dynamic increment, separated rigid/flexible truth, statistical significance, the paper's complete estimator and online bidirectional control independence.

Next action: use the band-target diagnostic as a comparator and initializer when adapting the joint mechanical/channel estimator to the supported plant. Do not claim the diagnostic itself completes the main method's hardware validation.
