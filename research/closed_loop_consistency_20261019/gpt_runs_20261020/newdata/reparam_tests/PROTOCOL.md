# P6 / P7 frozen analysis protocol

The decision rules are those in `时延扩展_理论推导_v1.md`, version 1.1, section 8. The source is hashed in `baseline/input_hashes.json`. No parent files may be changed.

* P6: assess common-delay agreement between H300 and H493 against model (a), report the approximately 13 microsecond reference as well as each matched record. Model (d) must not increase the 359 Hz and 493 Hz Q/E[Q] by more than 20%, and must not increase run-to-run modal dispersion. A >20% increase at 333 Hz or 493 Hz under (c)/(d) is evidence to retain C1.
* The H300 model contains four terms, including its existing broad background term. The same constraints apply to all four terms; failure of the constrained background is not proof of a physical complex mode.
* H300 Q/E[Q] is also resolved into fixed local intervals: 333 Hz uses 307.5–346 Hz; 359 Hz uses 346–420 Hz. These intervals are evaluation intervals only. The estimation window stays 220–420 Hz. H493 uses its entire original window.
* For comparing Q/E[Q], all four models use the same covariance estimated from the free-model M3 residual on nonexcited lines, with 0.5 Hz smoothing. Model-specific re-estimated covariance is secondary.
* Report both six-run standard deviations and pooled paired A/C repeatability. Neither is assumed to measure noise alone.
* P7: damping bias must be less than two Monte Carlo standard errors, separately for each P and mode. Report failures and variance inflation; real-data agreement uses ONE combined jackknife SE, as specified in P7.
* Jackknife correction is applied to the fitted log frequency/damping coordinates. Report errors as 100 expm1(log estimate − log truth). This does not assert exact finite-sample unbiasedness.
* Simulation seeds and settings follow experiment A. Each P has 40 replicates, including the 20 previously unrun P=16 seeds.
* Linear IV elimination is closed form only for frozen instruments. Reference instruments are updated inside the profile until the linear moment equations hold. Rank-one shape directions are additional nonlinear coordinates; they cannot be eliminated as unrestricted linear residues.
* Rank-one gauge: a unit output shape with positive pivot component; input shape absorbs amplitude. Shared C0/C1 additionally use (kappa0,kappa1)=(cos gamma,sin gamma), avoiding a residual scale ambiguity.
* Root equivalence is a prerequisite. The verification reports absolute differences and scale-aware differences for all parameters, strict moment scores, and prediction differences. Failure to verify must not be described as completion.
* Computational environment: Python 3.11, NumPy/SciPy, 16 logical CPUs, 32 GB RAM; at most four independent workers, one BLAS thread each. Every fitted endpoint is saved independently, including failures.
