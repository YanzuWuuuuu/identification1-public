function out = rcp_metrics(truth, estimate, nominal, f)
%RCP_METRICS Evaluate calibration on the true object, not its fitted model.
%   Four channels, rank-three rigid residue and four flexible pole pairs.
%   f is a sorted positive evaluation grid in Hz. All returned values are
%   real; NaN marks undefined ratios (JSON encodes these as null).
%   Target flexibility excludes D. The leakage bound uses oracle errors;
%   it is a deterministic check, not an estimated confidence certificate.

f = f(:);
assert(numel(f) >= 2 && all(f > 0) && all(diff(f) > 0), ...
    'Evaluation frequencies must be positive and strictly increasing.');
[ny, nu, st] = null_pair(truth.M);
[ey, eu, se] = null_pair(estimate.M);
[my, mu] = null_pair(nominal.M);
ey = align_sign(ey, ny); eu = align_sign(eu, nu);
my = align_sign(my, ny); mu = align_sign(mu, nu);
delta = norm(truth.M - estimate.M, 2);

out.rigid.relative_matrix_error = delta / norm(truth.M, 2);
out.rigid.error_over_true_gap = delta / st(3);
out.rigid.sigma3_truth = st(3);
out.rigid.sigma3_estimate = se(3);
out.rigid.estimate_sigma4_over_sigma3 = se(4) / se(3);
out.rigid.output_null_angle_deg = acosd(min(1, abs(ny.' * ey)));
out.rigid.input_null_angle_deg = acosd(min(1, abs(nu.' * eu)));
out.rigid.nominal_output_null_angle_deg = acosd(min(1, abs(ny.' * my)));
out.rigid.nominal_input_null_angle_deg = acosd(min(1, abs(nu.' * mu)));

dy = estimate.ty(:) - estimate.ty(1);
du = estimate.tu(:) - estimate.tu(1);
dy0 = truth.ty(:) - truth.ty(1);
du0 = truth.tu(:) - truth.tu(1);
Ly = max([0; dy]) + 10e-6;
Lu = max([0; du]) + 10e-6;
Ly0 = max([0; dy0]) + 10e-6;
Lu0 = max([0; du0]) + 10e-6;
out.delay.output_relative_error_us = 1e6 * (dy - dy0);
out.delay.input_relative_error_us = 1e6 * (du - du0);
out.delay.output_max_abs_error_us = max(abs(out.delay.output_relative_error_us));
out.delay.input_max_abs_error_us = max(abs(out.delay.input_relative_error_us));
out.compensation.padding_output_us = 1e6 * Ly;
out.compensation.padding_input_us = 1e6 * Lu;

N = numel(f);
rigidAmp = zeros(N, 4); flexAmp = zeros(N, 4);
etaY = zeros(N, 1); etaU = zeros(N, 1);
bound = nan(N, 1); qAmp = zeros(N, 1);
gainY = zeros(N, 1); gainU = zeros(N, 1);
relY = zeros(N, 4); relU = zeros(N, 4);
consistency = zeros(N, 1);
frfError = zeros(N, 1); frfTruth = zeros(N, 1);
targetIdeal = zeros(N, 1);
nullY = [my, ey, ey, ny]; nullU = [mu, eu, eu, nu];
boundValid = delta < se(3);

for k = 1:N
    s = 2i * pi * f(k);
    [hy, hu] = true_filters(truth, s);
    [rhY, rhU, dhat] = estimated_relative(estimate, s);
    rtY = hy / hy(1); rtU = hu / hu(1);
    errY = rtY ./ rhY - 1; errU = rtU ./ rhU - 1;
    relY(k, :) = abs(errY).'; relU(k, :) = abs(errU).';
    etaY(k) = max(abs(errY)); etaU(k) = max(abs(errU));
    cy = exp(-s * Ly) ./ rhY; cu = exp(-s * Lu) ./ rhU;
    cy0 = exp(-s * Ly0) ./ rtY; cu0 = exp(-s * Lu0) ./ rtU;
    CY = [ones(4, 2), cy, cy0]; CU = [ones(4, 2), cu, cu0];
    gainY(k) = max(abs(cy)); gainU(k) = max(abs(cu));
    F = flexible_part(truth, s);
    Gr = diag(hy) * (truth.M / s^2) * diag(hu);
    Gf = diag(hy) * F * diag(hu);
    for m = 1:4
        left = (nullY(:, m) .* CY(:, m)).';
        right = CU(:, m) .* nullU(:, m);
        rigidAmp(k, m) = abs(left * Gr * right);
        flexAmp(k, m) = abs(left * Gf * right);
    end
    q = exp(-s * (Ly + Lu)) * hy(1) * hu(1);
    qAmp(k) = abs(q);
    targetIdeal(k) = abs(q * (ny.' * F * nu));
    if boundValid
        bound(k) = abs(q) / abs(s)^2 * ( ...
            delta^2 / (se(3) - delta) + ...
            delta * (etaY(k) + etaU(k)) + ...
            norm(truth.M, 2) * etaY(k) * etaU(k));
    end
    % Exact two-sided denominators obey d_jk*d_11=d_j1*d_1k.
    residual = dhat .* dhat(1, 1) - dhat(:, 1) * dhat(1, :);
    scale = max(abs(dhat .* dhat(1, 1)), ...
        abs(dhat(:, 1) * dhat(1, :)));
    consistency(k) = max(abs(residual(:)) ./ max(scale(:), realmin));
    Gtrue = Gr + Gf + diag(hy) * truth.D * diag(hu);
    Pest = estimate.M / s^2 + estimate.D + flexible_part(estimate, s);
    Gest = exp(-s * (estimate.ty(:) + estimate.tu(:).')) .* Pest ./ dhat;
    frfError(k) = norm(Gest - Gtrue, 'fro');
    frfTruth(k) = norm(Gtrue, 'fro');
end

out.relative.output_max_error_per_channel = max(relY, [], 1).';
out.relative.input_max_error_per_channel = max(relU, [], 1).';
out.relative.output_max_error = max(etaY);
out.relative.input_max_error = max(etaU);
out.compensation.output_max_gain = max(gainY);
out.compensation.input_max_gain = max(gainU);
out.denominator.max_relative_separability_residual = max(consistency);
out.fit.truth_relative_frf_l2_error = sqrt(trapz(f, frfError.^2) / ...
    trapz(f, frfTruth.^2));
out.fit.note = 'Truth FRF error; held-out noise-level residual is evaluated by the runner.';

labels = {'nominal', 'static_only', 'joint', 'oracle'};
Er = trapz(f, rigidAmp.^2, 1); Ef = trapz(f, flexAmp.^2, 1);
idealEnergy = trapz(f, targetIdeal.^2);
for m = 1:4
    v.rigid_energy = Er(m);
    v.target_flexible_energy = Ef(m);
    v.rigid_l2 = sqrt(Er(m));
    v.target_flexible_l2 = sqrt(Ef(m));
    v.leak_to_flex_energy_ratio = safe_ratio(Er(m), Ef(m));
    v.leak_to_flex_amplitude_ratio = sqrt(v.leak_to_flex_energy_ratio);
    v.leak_to_flex_db = safe_db(Er(m), Ef(m));
    v.rigid_reduction_vs_nominal_db = safe_db(Er(1), Er(m));
    v.flex_energy_over_ideal = safe_ratio(Ef(m), idealEnergy);
    out.baseline.(labels{m}) = v;
end

% At a true mechanical pole, all channel filters are analytic. Evaluate
% their actual residue factors, without refitting the calibrated signal.
for i = 1:4
    s = truth.lambda(i);
    [hy, hu] = true_filters(truth, s);
    [rhY, rhU] = estimated_relative(estimate, s);
    rtY = hy / hy(1); rtU = hu / hu(1);
    CY = [ones(4, 2), exp(-s * Ly) ./ rhY, exp(-s * Ly0) ./ rtY];
    CU = [ones(4, 2), exp(-s * Lu) ./ rhU, exp(-s * Lu0) ./ rtU];
    q = hy(1) * hu(1) * [1, 1, exp(-s * (Ly + Lu)), ...
        exp(-s * (Ly0 + Lu0))];
    rho = ny.' * truth.R(:, :, i) * nu;
    out.mode(i).frequency_hz = imag(s) / (2 * pi);
    out.mode(i).ideal_projected_residue_abs = abs(rho);
    out.mode(i).ideal_projection_fraction = abs(rho) / norm(truth.R(:, :, i), 'fro');
    for m = 1:4
        r = (nullY(:, m) .* CY(:, m) .* hy).' * ...
            truth.R(:, :, i) * (hu .* CU(:, m) .* nullU(:, m));
        normalized = r / q(m);
        u.actual_residue_abs = abs(r);
        u.common_factor_removed_abs = abs(normalized);
        u.retention_amplitude_ratio = safe_ratio(abs(normalized), abs(rho));
        u.relative_complex_residue_error = safe_ratio(abs(normalized - rho), abs(rho));
        out.mode(i).(labels{m}) = u;
    end
end

out.bound.valid_gap_condition = boundValid;
out.bound.oracle_errors_used = true;
out.bound.is_confidence_certificate = false;
if boundValid
    out.bound.rigid_l2_upper = sqrt(trapz(f, bound.^2));
    out.bound.max_absolute_excess = max([0; rigidAmp(:, 3) - bound]);
    numericalFloor = eps(max(frfTruth));
    out.bound.max_actual_over_bound_with_roundoff_floor = ...
        max(rigidAmp(:, 3) ./ max(bound, numericalFloor));
else
    out.bound.rigid_l2_upper = NaN;
    out.bound.max_absolute_excess = NaN;
    out.bound.max_actual_over_bound_with_roundoff_floor = NaN;
end
out.frequency_min_hz = f(1); out.frequency_max_hz = f(end);
out.frequency_count = N;
out.baseline_order = labels;
out.curves = struct('frequency_hz', f, 'rigid_amplitude', rigidAmp, ...
    'target_flexible_amplitude', flexAmp, 'rigid_bound_joint', bound, ...
    'output_relative_error', relY, 'input_relative_error', relU, ...
    'joint_output_max_gain', gainY, 'joint_input_max_gain', gainU, ...
    'common_factor_amplitude', qAmp, 'ideal_target_amplitude', targetIdeal, ...
    'denominator_consistency', consistency);
end

function [y, u, sigma] = null_pair(M)
[U, S, V] = svd(M);
y = U(:, end); u = V(:, end); sigma = diag(S);
end

function a = align_sign(a, reference)
if a.' * reference < 0, a = -a; end
end

function [hy, hu] = true_filters(t, s)
hy = exp(-s * t.ty(:)) ./ (1 + s * t.filter_y(:));
hu = exp(-s * t.tu(:)) ./ (1 + s * t.filter_u(:));
end

function [ry, ru, d] = estimated_relative(t, s)
d = 1 + t.a * s + t.b * s^2;
ry = exp(-s * (t.ty(:) - t.ty(1))) .* (d(1, 1) ./ d(:, 1));
ru = exp(-s * (t.tu(:) - t.tu(1))) .* (d(1, 1) ./ d(1, :).');
end

function F = flexible_part(t, s)
F = complex(zeros(4));
for k = 1:4
    F = F + t.R(:, :, k) / (s - t.lambda(k)) + ...
        conj(t.R(:, :, k)) / (s - conj(t.lambda(k)));
end
end

function r = safe_ratio(a, b)
if b > 0, r = a / b; else, r = NaN; end
end

function r = safe_db(a, b)
if a > 0 && b > 0, r = 10 * log10(a / b); else, r = NaN; end
end
