# Offline structured calibration pilot — 2026-10-16

Objective: test whether joint estimation of a common low-frequency subspace and relative channel delays improves held-out nullspace responses.

Accepted comparison: current frozen training/validation LRM arrays. Same first/second half split, same frozen inverse-magnitude quadratic weights and low-band coefficient 0.03. No validation-based branch selection. Start at z=0 then repeat the normal positions -0.5/0.5/1 mm. Exclude anomalous -1 mm.

Model: SRIVC initializes 3 real rank-one second-order support modes, whose spatial factors share a three-dimensional left/right space. This is a proportional damping support-model subclass. Upper 12 pole pairs remain frozen at their training estimates, with free real numerators and D profiled as nuisance terms. Seven additive delays are jointly optimized with support shapes, poles and damping. This is constrained output-error estimation after SRIVC, not pure SRIVC, not the unchanged double-integrator theorem, and not a full physical high-mode realization.

Sequence: complex-residue SVD precheck; verify analytic profile gradient; fit support structure at the old delays; jointly update from two prespecified delay starts; select by training cost among feasible positive-exit candidates; evaluate on held-out LRM directly.

Fair comparison: B0 nominal geometry, B1 estimated static direction, B4 the same B1 direction plus causal all-pass delay alignment. Input command and output combination norms are held at one. Nominal Rx/Ry/Z reference directions fixed. Report output-side, input-side, double-sided residuals. Bands 3–25,25–80,100–170,260–310,25–800 Hz; no selection using best-looking validation band.

Success evidence: low-frequency held-out residual reduction with stable spatial directions and retained first-plate-band response. Delays should be judged as effective model parameters, not certified hardware delays. Failure or weak effects remain reported. No real feedback/controller changes.

Limits: closed-loop total command records lack a usable independent reference; second-half split is not an independent repeated experiment. No independent sensor-noise truth. Low-frequency total response is not ground-truth rigid leakage. Low-model exact nulling is only an algebra check.

Outputs: calibration_z*.mat, evaluation files and plots, consolidated Chinese report. Prior data/model/paper preserved.
