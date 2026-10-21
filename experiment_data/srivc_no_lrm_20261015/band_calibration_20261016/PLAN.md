# Band-limited calibration diagnostic, 2026-10-16

Question: is the weak two-sided result caused by estimating the null directions from poorly fitted support resonances, or is there little removable channel mismatch in the measured band?

This is an auxiliary diagnostic baseline, not a replacement for the paper's structured modal estimator or its consistency theorem. No new data or hardware changes. Previous data, splits, fits and metrics remain read-only.

Primary estimation band is fixed before execution at 25–80 Hz, above the observed support resonances and below the unresolved 100–170 Hz features. Use the original first-half LRM for fitting, original second-half LRM for checking. This previously inspected half is not a new blind test. The original four normal positions are retained; anomalous z=-1 remains excluded.

For each side, minimize the integrated full-port residual over real unit null vectors. At fixed relative delays this is the smallest eigenpair of a real symmetric covariance. Delays are profiled with reference channel 1 fixed to zero, bounds +/-0.25 ms, and three prescribed starts (zero, previous modal estimate, its negative). Select by training objective only; record boundaries and convergence. A zero-delay static solution is an explicit comparator. The common delay is unidentifiable and not estimated.

Comparators: B0 nominal; B1/B4 previous structured estimate without/with delay; D0 direct static; D1 same directions as D2 without delay; D2 direct directions with causal delay alignment. Keep nominal logical directions and unit row/column norms fixed. Retain original three cross/null-block RMS metrics in fixed bands 3–25,25–80,100–170,260–310,25–800. Add full-port null residuals to prevent hiding energy in omitted coordinates. Show absolute values, same-band B0 ratios and D2/D1 incremental benefit. Never call all total residuals rigid leakage or all high-band residuals flexible preservation.

Transfer check: reuse the z=0 D0/D1/D2 matrices without refitting at the other positions. Training band sensitivities 25–60 and 35–80 Hz are diagnostics at z=0, not candidates for validation-based selection. An expanded +/-1 ms delay-bound run at z=0 checks boundary dependence only. These are not new main test sets.

Validation: check eigenvector orientation, real-unit constraints, analytic delay gradient against finite differences, direct static training cost <= nominal, optimized training cost <= zero-delay candidate, causal-delay sign identity, and equal B1/B4 or D1/D2 operator norms. Retain all candidates and raw metrics. Stop after this prescribed campaign; report inconclusive/negative effects honestly.

Outputs: MAT inputs/config/results, CSV metrics/directions/delays, PNG/vector-PDF figures, concise Chinese result report and bounded experimental subsection in the existing paper draft if supported. No claims of physical delay identification, measured noise reduction, closed-loop deployment or statistical significance.
