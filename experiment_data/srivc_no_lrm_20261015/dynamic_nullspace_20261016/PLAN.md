# Algorithm revision: constrained stable dynamic nullspace

User request: improve the mathematics/algorithm, not the presentation of modest static improvements.

Hypothesis: a fixed real direction plus pure delay does not represent the frequency-dependent interaction. Identify the two cross blocks directly using stable triangular dynamic coordinate corrections while preserving the nominal logical block and the baseline null-coordinate coefficient.

For the fixed nominal transformed FRF Z=[A b;c d], use L=[I 0; fy 1], R=[I fu;0 1]. Cross blocks are c+fy A and b+A fu; the target block changes by c fu+fy b+fy A fu. The complete transform is stably invertible if fy,fu are stable/proper. This is a mathematical extension to general dynamic direction calibration, not pure physical delay identification or the original modal estimator unchanged.

Primary candidates fixed before execution: (S) constant triangular corrections; (R40) two stable proper bases phi0=p^2/(s+p)^2 and phi1=p*s/(s+p)^2, p=2*pi*40, six real coefficients per side. Both are fit to the same 25–80 Hz training FRFs with unaltered quadrature RMS. No validation selection. z=0 sensitivities use p=20,80 only; the primary remains R40.

Constraints: each physical sensor-combination row and actuator-command column has norm <=1.05 on the frequency guard grid. The nominal null direction is orthogonal to the frozen logical directions, hence the anchor component equals one and the norm cannot fall below one. Protect the complex N-to-N response in 260–310,320–380,450–800 Hz: at most 10% training RMS deviation in each band. Use ridge 1e-6 for numerical regularization. Solve each side's convex quadratic subproblem with analytic gradients, alternating at most eight sweeps under the bilinear joint target-response constraints. Do not claim joint global optimality.

Data: preserve original four normal-position first-half LRM training/second-half checks. Exclude abnormal z=-1. The previously inspected second halves are not blind tests and do not eliminate closed-loop FRF bias. Evaluate all original bands/cross blocks/NN with fixed same-band B0 denominators and no output renormalization. Report original B1/B4 and D0/D2 alongside S/R40.

Checks: analytic objective and constraint gradients, stable poles/proper filters, exact unchanged A block, anchor and minimum norm, sensor/command H-infinity gain, target-band complex deviation, solver feasibility, training cost monotonicity, holdout and z0 cross-position transfer. The measured NN response is not separated flexible truth; preserve its complex response, not just peak amplitude. Do not deploy to hardware or claim closed-loop independence without controller verification.

Outputs: new code/MAT/CSV, complete magnitude/phase and cross-block plots, reproducible mathematical report with all failures. Preserve earlier code/data/manuscript; update paper only if the new method and evidence warrant it after inspection.

## Evidence-driven extension after primary z0 run

R40 training cross squared-RMS costs were 0.142/0.519, with Hinf 1.006/1.001. The gain budget was not limiting. A nonimplementable frequency-by-frequency inverse fitted on training and frozen on the second half gave cross RMS ratios 0.056/0.332 at z0; corresponding half-record differences were 0.059/0.326. Thus output-side low-order approximation has room, whereas input-side record disagreement is already substantial. These pointwise inverse scores are diagnostics, not performance bounds or deployable filters.

Add one nested R4 candidate: retain both R40 bases, append phi2=phi0/(1+s/(2*pi*20)), phi3=phi2/(1+s/(2*pi*80)). This increases coefficients from 6 to 12 per side with stable common denominator degree four. Warm-start by padding the R40 coefficients with zeros; all gain/target/data settings unchanged. Record unconstrained least-squares cost and protection-budget usage to distinguish basis limitation from protection constraints. Do not use validation to choose a winning pole/regularization combination. Original R40 remains reported.

R4 analysis: output cost 0.1121 equals its unconstrained linear least-squares value; input cost 0.4632 vs unconstrained 0.4628. Gain remains below budget. Thus low-order basis limitation, not the protection constraint, explains most of the remaining training residual. Add a bounded nested Laguerre complexity check with L=4 then L=8 and the same p=40 Hz. Keep phi0,phi1; append (1/(1+q))*((q-1)/(q+1))^ell for ell=2,...,L-1. Use an explicitly stable shared cascade state-space realization, not unstable polynomial roots or a long FIR. Same data and budgets; retain all candidates, do not choose order using held-out results. L8 is the final higher-complexity candidate for this pass.

## Causal-alignment extension after L8 evaluation

L8 held-out z0 ratios are 0.2858/0.6069, still modest versus D0=0.3538/0.6842. Output gain reaches 1.05. Test one fixed additional common null-channel alignment delay T=2 ms per side, with the same L8 bases and budgets. This is a new implementation hypothesis, not an order sweep: Ny=a ny0^T+fy P and Nu=a nu0+Q fu, a=exp(-sT). Fit a*c+fy*A and a*b+A*fu. Protect NN relative to the explicitly transported reference a^2*d (4 ms round trip), and report both transported-reference and raw complex differences. The actual implementation is causal and stable; its coordinate inverse need not be causal. Do not claim real physical channel delay identification or free online improvement. The added latency must be checked with the real controller before use. Evaluate T=0 and fixed T=2 ms separately; no post-validation latency selection.

## Full-band audit correction

Unprotected L8 increased 3–25 Hz held-out cross responses by factors 1.49–3.94 and 100–170 Hz by 1.07–1.25. Thus it cannot be promoted as broadly useful even though 25–80 Hz improved. Add fixed per-side convex quadratic cross-response constraints in 3–25 and 80–170 Hz: band RMS must not exceed the same nominal B0 training response. Use the already fixed L8, p=40, zero extra latency, same gain/target budgets. Tag these candidates _g. These guards are band-energy constraints, not pointwise or closed-loop stability guarantees. This correction follows a discovered failure; the already examined second-half data remain an exploratory diagnostic, not a blind validation set. Retain and report the failed unguarded candidates.
