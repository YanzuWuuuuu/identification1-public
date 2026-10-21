# Modal residue and recorded-coordinate screening

Source result: `C:\Users\Administrator\Desktop\my_article\identifiation1\experiment_data\srivc_no_lrm_20261015\spring_identification_20261016\spring_selected.mat` (modified 16-10月-2026 09:29:00, 67642899 bytes).

Audit created: 2026-10-16 09:30:57. Fit converged=1, score=6.38891e-07. No poles, delays, or source data were changed.

For q=s/(2*pi*fref), A=1+a1*q+a2*q^2, the positive-imaginary pole p has residue R=(b0+b1*q(p))/((a1+2*a2*q(p))/(2*pi*fref)). Rows/columns are outputs/inputs; flattened column index is y+4*(u-1). A constant direct term has no pole residue.

Two views are reported: `rational_delay_removed` and `observed_delay_included`. The latter is R_observed=R.*exp(-p*ell), evaluated at the complex pole. Delay additivity residual is 1.0842e-19 s. Removing estimated delay is model dependent, not independently verified hardware calibration.

## Spatial references

- Nominal sensor coordinates come from `matlab/anlaysis/test.m:6-9`: (-.095,.160), (.095,.160), (.095,-.160), (-.095,-.160) m. The nominal null vector is [1,-1,1,-1]/2.
- C1 and five C2 columns are copied from `matlab/0/Multiple_main_Control.m:78-85`; their nominal modal frequencies are 286.43, 344.58, 488.27, 705.67, 738.83 Hz. These are the old simulation references, not experimental ground truth.
- Channel order, sensor gains, signs, timing, and exact hardware geometry are not confirmed. All geometric fractions and MAC values are exploratory.

## Meaning of metrics

`sigma2_over_sigma1` and `first_singular_energy_fraction` describe algebraic residue rank concentration. They do not establish a plate/suspension origin. Projection columns are amplitude norm ratios, not squared energy fractions. MAC compares the leading complex left singular vector with each original C2 column and is phase invariant. Near-equal leading singular values make the individual leading vector ambiguous.

Four sensors give only one nominal non-rigid direction. Small null projection does not prove rigid motion: a flexible shape may have rigid-compatible values at these four points. Large projection may also reflect geometry, gain, phase, or model errors. No component receives an automatic physical label.

## Recorded signal mapping

Raw source: `C:\Users\Administrator\Desktop\my_article\identifiation1\experiment_data\auto_identify_motor1234_z0.mat`. Status: completed.

Only the four pre_y records and three y_0__k_ records were loaded; metadata lists all 17 file variables. The first second was excluded. A static affine map Z=Y*L'+offset' was fitted on the first half and predicted on the second half. An all-sample map and both split maps are also saved. No coordinate names or physical units are inferred from variable names.

All-sample map rank=3; standardized sensor design rank=4, condition=25.0435. A deficient/ill-conditioned sensor design limits uniqueness/stability of the recovered map. Second-half relative residuals are [5.70707e-15 7.77427e-15 2.70341e-14], with each denominator the corresponding centered output norm.

The recorded-map row/null spaces represent an existing signal transformation. They do not independently establish the true rigid-body geometry. Their directly consistent comparison is the `observed_delay_included` residue view, because the affine map was fitted to raw sensor records. The rational view is also tabulated only as a coordinate-sensitivity check.

The same acquisition and previously reused split are being diagnosed; this is not a new independent validation experiment. Current calculations do not identify which hardware part creates the 124/128/138 Hz features.
