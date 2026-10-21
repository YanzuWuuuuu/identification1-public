# Bounded response to Claude's second review, 2026-10-14

Parent: closed_loop_validation/outputs_structured (five-mode D=0, 125 parameters).
Use analysis-campaign, analysis-lite route. No new main fit or physical model.
Available: MATLAB R2025a, CPU, about 18 GiB free RAM. Each Jacobian is roughly
6272 x 125 real doubles. Serial run; no GPU or workers needed. Target under
10 minutes compute, saved checkpoints. Artifact backend is unavailable; local
files and the existing manuscript provide lineage.

S1 (C1/C2, sec.4 information): inspect the main selected fit with unchanged
data/weights/model. Compare existing dimensionless coordinates with the invertible
shear eta=ell+a (retain a and b); also column-normalize. Inspect weakest right
singular vectors, cancellation delta ell + delta a, and invariance of calibrated
target covariance. Do not claim coordinate change creates information or that
it makes a fixed numerical score threshold meaningful. No repeated optimization.

S2 (C1/C2, sec.6 auxiliary): true-object local mean-information scan at controller
design frequencies 15,30,60,100,150,200 Hz, derivative rolloff fixed 120 Hz;
same 5 kHz ZOH model, 20:4:800 Hz frequencies, four directions, 0.1 N line peaks,
100 nm sample output noise, N=32. Check full closed-loop stability FIRST. Reject
unstable designs; do not reinterpret these as valid experiments. Report design
frequency separately from loop crossover. Use full output noise spectral
covariance in the stationary-frequency approximation; this differs from the
finite-record diagonal fitting weights. Include complex versus proportional
real modal tangents, delay profile information and relative-response target
uncertainty. Compare fixed injected d with ideal feedforward-shaped d that fixes
Su D (report extra force cost; no saturation/deployment guarantee).

S3 (C3, sec.6 main support): reuse saved independent validation records, fit,
feedback and compensation; report 260--310 Hz and existing 202--798 Hz bands,
with explicit rigid/flexible/noise decomposition and mode retention. No refit,
new data or test-based selection. Note low-band flexibility is quasi-static.

Stop: all three checks saved and relevant manuscript wording corrected. Keep
null/negative results. Do not undertake global convergence, coverage studies,
new controller design or wholesale manuscript shortening in this response.
