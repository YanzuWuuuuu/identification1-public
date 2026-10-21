function out = d_residual_metrics(core5, estimate5, full8or5, nominal, f)
%D_RESIDUAL_METRICS Assess five target modes with an omitted high-mode tail.
%   CORE5 contains the five target modes and true channel filters. FULL8OR5
%   shares those components and may contain additional flexible modes.
%   ESTIMATE5 contains five fitted modes and an optional fitted constant D.
%   No fitted D is counted as target flexibility. A true nonzero D is not
%   needed: omitted high-mode dynamics are evaluated explicitly in FULL8OR5.
%   All four calibration baselines act on the SAME true rigid/five-mode
%   components and filters. FRF errors use integration in df (Hz), not a
%   weighted or equally spaced point norm. f must be increasing and positive.

truth = core5; estimate = estimate5; fullTruth = full8or5;
f = f(:);
assert(numel(f) >= 2 && all(f > 0) && all(diff(f) > 0), ...
    'Evaluation frequencies must be positive and strictly increasing.');
assert(numel(truth.lambda) == 5 && numel(estimate.lambda) == 5, ...
    'Exactly five core and fitted mechanical modes are required.');
assert(numel(fullTruth.lambda) >= 5, 'The full object must contain all five target modes.');
sharedFields = {'M','ty','tu','filter_y','filter_u'};
for k = 1:numel(sharedFields)
    name = sharedFields{k};
    assert(isequaln(truth.(name), fullTruth.(name)), ...
        'Core and full objects must share their true rigid matrix and channel dynamics.');
end
assert(isequaln(truth.lambda(:), reshape(fullTruth.lambda(1:5),[],1)) && ...
    isequaln(truth.R, fullTruth.R(:,:,1:5)), ...
    'The first five full-object modes must be the specified target modes.');
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
referenceFrequency = abs(truth.lambda(:))/(2*pi);
referenceDamping = -real(truth.lambda(:))./abs(truth.lambda(:));
out.reference_frequencies_hz = referenceFrequency;
out.estimated_frequencies_hz = abs(estimate.lambda(:))/(2*pi);
out.estimated_damping = -real(estimate.lambda(:))./abs(estimate.lambda(:));
out.frequency_relative_error = (out.estimated_frequencies_hz-referenceFrequency)./referenceFrequency;
out.damping_relative_error = (out.estimated_damping-referenceDamping)./referenceDamping;

dy = estimate.ty(:) - estimate.ty(1); du = estimate.tu(:) - estimate.tu(1);
dy0 = truth.ty(:) - truth.ty(1); du0 = truth.tu(:) - truth.tu(1);
Ly = max([0; dy]) + 10e-6; Lu = max([0; du]) + 10e-6;
Ly0 = max([0; dy0]) + 10e-6; Lu0 = max([0; du0]) + 10e-6;
out.delay.output_relative_error_us = 1e6 * (dy - dy0);
out.delay.input_relative_error_us = 1e6 * (du - du0);
out.delay.output_max_abs_error_us = max(abs(out.delay.output_relative_error_us));
out.delay.input_max_abs_error_us = max(abs(out.delay.input_relative_error_us));
out.compensation.padding_output_us = 1e6 * Ly;
out.compensation.padding_input_us = 1e6 * Lu;

N = numel(f);
rigidAmp = zeros(N,4); flexAmp = zeros(N,4); tailAmp = zeros(N,4);
calibratedFive = complex(zeros(N,4));
etaY = zeros(N,1); etaU = zeros(N,1);
bound = nan(N,1); qAmp = zeros(N,1);
gainY = zeros(N,1); gainU = zeros(N,1);
relY = zeros(N,4); relU = zeros(N,4); consistency = zeros(N,1);
fullError = zeros(N,1); fullNorm = zeros(N,1);
coreError = zeros(N,1); coreNorm = zeros(N,1);
fiveError = zeros(N,1); fiveNorm = zeros(N,1);
targetIdeal = zeros(N,1); targetIdealComplex = complex(zeros(N,1));
nullY = [my, ey, ey, ny]; nullU = [mu, eu, eu, nu];
boundValid = delta < se(3);
for k = 1:N
    s = 2i*pi*f(k);
    [hy, hu] = true_filters(truth,s);
    [rhY, rhU, dhat] = estimated_relative(estimate,s);
    rtY = hy/hy(1); rtU = hu/hu(1);
    errY = rtY./rhY-1; errU = rtU./rhU-1;
    relY(k,:) = abs(errY).'; relU(k,:) = abs(errU).';
    etaY(k) = max(abs(errY)); etaU(k) = max(abs(errU));
    cy = exp(-s*Ly)./rhY; cu = exp(-s*Lu)./rhU;
    cy0 = exp(-s*Ly0)./rtY; cu0 = exp(-s*Lu0)./rtU;
    CY = [ones(4,2),cy,cy0]; CU = [ones(4,2),cu,cu0];
    gainY(k) = max(abs(cy)); gainU(k) = max(abs(cu));
    F = flexible_part(truth,s);
    Ffull = flexible_part(fullTruth,s);
    Gr = diag(hy)*(truth.M/s^2)*diag(hu);
    Gf = diag(hy)*F*diag(hu);
    Gtail = diag(hy)*(Ffull-F)*diag(hu);
    for m = 1:4
        left = (nullY(:,m).*CY(:,m)).'; right = CU(:,m).*nullU(:,m);
        rigidAmp(k,m) = abs(left*Gr*right);
        calibratedFive(k,m) = left*Gf*right;
        flexAmp(k,m) = abs(calibratedFive(k,m));
        tailAmp(k,m) = abs(left*Gtail*right);
    end
    q = exp(-s*(Ly+Lu))*hy(1)*hu(1); qAmp(k) = abs(q);
    targetIdealComplex(k) = q*(ny.'*F*nu);
    targetIdeal(k) = abs(targetIdealComplex(k));
    if boundValid
        bound(k) = abs(q)/abs(s)^2*(delta^2/(se(3)-delta) + ...
            delta*(etaY(k)+etaU(k)) + norm(truth.M,2)*etaY(k)*etaU(k));
    end
    residual = dhat.*dhat(1,1)-dhat(:,1)*dhat(1,:);
    scale = max(abs(dhat.*dhat(1,1)),abs(dhat(:,1)*dhat(1,:)));
    consistency(k) = max(abs(residual(:))./max(scale(:),realmin));
    Gcore = Gr+Gf+diag(hy)*truth.D*diag(hu);
    Gfull = Gr+diag(hy)*(Ffull+fullTruth.D)*diag(hu);
    channelEstimate = exp(-s*(estimate.ty(:)+estimate.tu(:).'))./dhat;
    Fhat = flexible_part(estimate,s);
    Gest = channelEstimate.*(estimate.M/s^2+estimate.D+Fhat);
    GfiveEstimate = channelEstimate.*Fhat;
    fullError(k) = norm(Gest-Gfull,'fro'); fullNorm(k) = norm(Gfull,'fro');
    coreError(k) = norm(Gest-Gcore,'fro'); coreNorm(k) = norm(Gcore,'fro');
    fiveError(k) = norm(GfiveEstimate-Gf,'fro'); fiveNorm(k) = norm(Gf,'fro');
end

out.relative.output_max_error_per_channel = max(relY,[],1).';
out.relative.input_max_error_per_channel = max(relU,[],1).';
out.relative.output_max_error = max(etaY); out.relative.input_max_error = max(etaU);
out.compensation.output_max_gain = max(gainY);
out.compensation.input_max_gain = max(gainU);
out.denominator.max_relative_separability_residual = max(consistency);
out.fit = fit_band(f,fullError,fullNorm,coreError,coreNorm,fiveError,fiveNorm, ...
    calibratedFive(:,3)-targetIdealComplex,targetIdeal, true(N,1));
out.fit.band_200_800_hz = fit_band(f,fullError,fullNorm,coreError,coreNorm,fiveError,fiveNorm, ...
    calibratedFive(:,3)-targetIdealComplex,targetIdeal, f>=200 & f<=800);
out.fit.truth_relative_frf_l2_error = out.fit.full_object_relative_frf_l2_error;
out.fit.note = ['Unweighted continuous-frequency L2 ratios approximated by trapz in Hz. ' ...
    'Full-object FRF includes every true mode. Core-five FRF includes rigid motion and only five true modes; ' ...
    'its estimate includes fitted D, so this diagnostic also contains the deliberate tail approximation. ' ...
    'Target-five flexible FRF excludes rigid motion and D on both sides. ' ...
    'The joint calibrated-five error is measured against ideal five-mode projection with the same common delay.'];

labels = {'nominal','static_only','joint','oracle'};
Er = trapz(f,rigidAmp.^2,1); Ef = trapz(f,flexAmp.^2,1); Et = trapz(f,tailAmp.^2,1);
idealEnergy = trapz(f,targetIdeal.^2);
for m = 1:4
    v.rigid_energy = Er(m); v.target_flexible_energy = Ef(m);
    v.rigid_l2 = sqrt(Er(m)); v.target_flexible_l2 = sqrt(Ef(m));
    v.leak_to_flex_energy_ratio = safe_ratio(Er(m),Ef(m));
    v.leak_to_flex_amplitude_ratio = sqrt(v.leak_to_flex_energy_ratio);
    v.leak_to_flex_db = safe_db(Er(m),Ef(m));
    v.rigid_reduction_vs_nominal_db = safe_db(Er(1),Er(m));
    v.flex_energy_over_ideal = safe_ratio(Ef(m),idealEnergy);
    v.high_mode_tail_l2 = sqrt(Et(m));
    v.high_mode_tail_to_target_five_amplitude_ratio = sqrt(safe_ratio(Et(m),Ef(m)));
    out.baseline.(labels{m}) = v;
end

% True pole residues of the five target modes; appended modes never count as targets.
for i = 1:numel(truth.lambda)
    s = truth.lambda(i); [hy,hu] = true_filters(truth,s);
    [rhY,rhU] = estimated_relative(estimate,s);
    rtY = hy/hy(1); rtU = hu/hu(1);
    CY = [ones(4,2),exp(-s*Ly)./rhY,exp(-s*Ly0)./rtY];
    CU = [ones(4,2),exp(-s*Lu)./rhU,exp(-s*Lu0)./rtU];
    q = hy(1)*hu(1)*[1,1,exp(-s*(Ly+Lu)),exp(-s*(Ly0+Lu0))];
    rho = ny.'*truth.R(:,:,i)*nu;
    out.mode(i).frequency_hz = imag(s)/(2*pi);
    out.mode(i).undamped_frequency_hz = abs(s)/(2*pi);
    out.mode(i).ideal_projected_residue_abs = abs(rho);
    out.mode(i).ideal_projection_fraction = abs(rho)/norm(truth.R(:,:,i),'fro');
    for m = 1:4
        r = (nullY(:,m).*CY(:,m).*hy).'*truth.R(:,:,i)*(hu.*CU(:,m).*nullU(:,m));
        normalized = r/q(m);
        u.actual_residue_abs = abs(r); u.common_factor_removed_abs = abs(normalized);
        u.retention_amplitude_ratio = safe_ratio(abs(normalized),abs(rho));
        u.relative_complex_residue_error = safe_ratio(abs(normalized-rho),abs(rho));
        out.mode(i).(labels{m}) = u;
    end
end

out.bound.valid_gap_condition = boundValid; out.bound.oracle_errors_used = true;
out.bound.is_confidence_certificate = false;
if boundValid
    out.bound.rigid_l2_upper = sqrt(trapz(f,bound.^2));
    out.bound.max_absolute_excess = max([0;rigidAmp(:,3)-bound]);
    out.bound.max_actual_over_bound_with_roundoff_floor = ...
        max(rigidAmp(:,3)./max(bound,eps(max(fullNorm))));
else
    out.bound.rigid_l2_upper = NaN; out.bound.max_absolute_excess = NaN;
    out.bound.max_actual_over_bound_with_roundoff_floor = NaN;
end
out.frequency_min_hz = f(1); out.frequency_max_hz = f(end); out.frequency_count = N;
out.target_mode_count = numel(truth.lambda); out.full_mode_count = numel(fullTruth.lambda);
out.target_mode_indices = 1:5; out.baseline_order = labels;
out.units = struct('frequency','Hz','frf','m/N (physical port convention)', ...
    'integrated_squared_frf','(m/N)^2 Hz','integrated_frf_l2','(m/N) sqrt(Hz)', ...
    'relative_errors','dimensionless','delay','microseconds');
out.scope = ['All true rigid leakage and retained target flexibility are evaluated on the actual ' ...
    'core five-mode components with the common full-object channel dynamics. The added high modes ' ...
    'are reported as tail response, excluded from target flexibility, and included in full-object FRF error. ' ...
    'These are transfer-response measures; integrated L2 values are not measured displacement RMS without an input spectrum.'];
out.curves = struct('frequency_hz',f,'rigid_amplitude',rigidAmp, ...
    'target_flexible_amplitude',flexAmp,'high_mode_tail_amplitude',tailAmp, ...
    'rigid_bound_joint',bound,'output_relative_error',relY,'input_relative_error',relU, ...
    'joint_output_max_gain',gainY,'joint_input_max_gain',gainU, ...
    'common_factor_amplitude',qAmp,'ideal_target_amplitude',targetIdeal, ...
    'joint_target_complex_error_amplitude',abs(calibratedFive(:,3)-targetIdealComplex), ...
    'full_object_frf_error_amplitude',fullError,'full_object_frf_amplitude',fullNorm, ...
    'core_five_frf_error_amplitude',coreError,'core_five_frf_amplitude',coreNorm, ...
    'target_five_flexible_frf_error_amplitude',fiveError, ...
    'target_five_flexible_frf_amplitude',fiveNorm,'denominator_consistency',consistency);
end

function result = fit_band(f,eFull,nFull,eCore,nCore,eFive,nFive,eJoint,nJoint,mask)
result.requested_band_available = sum(mask)>=2;
if sum(mask)>=2
    result.frequency_min_hz = min(f(mask)); result.frequency_max_hz = max(f(mask));
    result.frequency_count = sum(mask);
    result.full_object_relative_frf_l2_error = relative_l2(f(mask),eFull(mask),nFull(mask));
    result.core_five_relative_frf_l2_error = relative_l2(f(mask),eCore(mask),nCore(mask));
    result.target_five_flexible_relative_frf_l2_error = relative_l2(f(mask),eFive(mask),nFive(mask));
    result.joint_calibrated_five_relative_complex_l2_error = relative_l2(f(mask),eJoint(mask),nJoint(mask));
else
    result.frequency_min_hz = NaN; result.frequency_max_hz = NaN; result.frequency_count = sum(mask);
    result.full_object_relative_frf_l2_error = NaN; result.core_five_relative_frf_l2_error = NaN;
    result.target_five_flexible_relative_frf_l2_error = NaN;
    result.joint_calibrated_five_relative_complex_l2_error = NaN;
end
end

function r = relative_l2(f,e,n)
r = sqrt(safe_ratio(trapz(f,abs(e).^2),trapz(f,abs(n).^2)));
end

function [y,u,sigma] = null_pair(M)
[U,S,V] = svd(M); y = U(:,end); u = V(:,end); sigma = diag(S);
end

function a = align_sign(a,reference)
if a.'*reference<0, a = -a; end
end

function [hy,hu] = true_filters(t,s)
hy = exp(-s*t.ty(:))./(1+s*t.filter_y(:));
hu = exp(-s*t.tu(:))./(1+s*t.filter_u(:));
end

function [ry,ru,d] = estimated_relative(t,s)
d = 1+t.a*s+t.b*s^2;
ry = exp(-s*(t.ty(:)-t.ty(1))).*(d(1,1)./d(:,1));
ru = exp(-s*(t.tu(:)-t.tu(1))).*(d(1,1)./d(1,:).');
end

function F = flexible_part(t,s)
F = complex(zeros(4));
for k = 1:numel(t.lambda)
    F = F+t.R(:,:,k)/(s-t.lambda(k))+conj(t.R(:,:,k))/(s-conj(t.lambda(k)));
end
end

function r = safe_ratio(a,b)
if b>0, r=a/b; else, r=NaN; end
end

function r = safe_db(a,b)
if a>0 && b>0, r=10*log10(a/b); else, r=NaN; end
end
