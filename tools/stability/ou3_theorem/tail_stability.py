"""Analytic lemmas for the single OU-III stability theorem.

H18 is a finite bridge to the implemented accelerometer-bias release.
Asymptotic contraction is required on the recurring magnetically informed A21
tail, not on the held-bias bridge.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

def finite_bridge_bound(v0: float, steps: int, gain: float, supply: float) -> float:
    """Iterate V+ <= gain*V+supply for a finite H18 bridge."""
    if not all(math.isfinite(x) for x in (v0,gain,supply)): raise ValueError("finite arguments required")
    if v0 < 0 or steps < 0 or gain < 0 or supply < 0: raise ValueError("nonnegative bridge data required")
    if steps == 0: return v0
    if gain == 1.0: return v0 + steps*supply
    return gain**steps*v0 + supply*(gain**steps-1.0)/(gain-1.0)

def release_bound_from_service(remaining_applied_updates: int, service_window_s: float,
                               guard_s: float, reference_release_s: float) -> float:
    """Finite gate-clear bound after finite magnetic-reference refinement."""
    if remaining_applied_updates < 0: raise ValueError("remaining update count must be nonnegative")
    if not all(math.isfinite(x) for x in (service_window_s,guard_s,reference_release_s)): raise ValueError("finite release data required")
    if service_window_s <= 0 or guard_s < 0 or reference_release_s < 0: raise ValueError("positive window and nonnegative times required")
    return reference_release_s + remaining_applied_updates*service_window_s + guard_s

@dataclass(frozen=True)
class A21TailPremises:
    linear_storage_ratio: float
    nonlinear_storage_lipschitz: float
    supply_gain: float
    coercivity_lower: float
    prefix_gain: float
    def __post_init__(self) -> None:
        v=(self.linear_storage_ratio,self.nonlinear_storage_lipschitz,self.supply_gain,self.coercivity_lower,self.prefix_gain)
        if not all(math.isfinite(x) for x in v): raise ValueError("finite premises required")
        if self.linear_storage_ratio < 0 or self.nonlinear_storage_lipschitz < 0: raise ValueError("nonnegative contraction data required")
        if self.supply_gain < 0 or self.coercivity_lower <= 0 or self.prefix_gain < 0: raise ValueError("valid supply/coercivity/prefix data required")

def nonlinear_tail_ratio(p: A21TailPremises) -> float:
    """Small-gain lift: rho=(sqrt(rho0)+eta)^2."""
    return (math.sqrt(p.linear_storage_ratio)+p.nonlinear_storage_lipschitz)**2

def practical_radius_bound(p: A21TailPremises, disturbance_bound: float) -> float:
    if not math.isfinite(disturbance_bound) or disturbance_bound < 0: raise ValueError("finite nonnegative disturbance bound required")
    rho=nonlinear_tail_ratio(p)
    if not rho < 1.0: raise ValueError("A21 tail is not strictly contractive")
    return math.sqrt(p.supply_gain/(p.coercivity_lower*(1.0-rho)))*disturbance_bound

def proof_route_status(*, reference_refinement_finite: bool, h18_bridge_retained: bool,
                       a21_linear_uniform: bool, a21_nonlinear_bound: bool,
                       a21_prefix_retained: bool) -> dict:
    flags=(reference_refinement_finite,h18_bridge_retained,a21_linear_uniform,a21_nonlinear_bound,a21_prefix_retained)
    if any(type(x) is not bool for x in flags): raise ValueError("literal booleans required")
    return {"h18_role":"finite bridge only","a21_role":"recurring asymptotic tail",
            "release_finite":reference_refinement_finite,"h18_bridge_closed":h18_bridge_retained,
            "a21_tail_closed":all(flags[2:]),"route_closed":all(flags)}


@dataclass(frozen=True)
class RefinementPremises:
    """Sufficient conditions for the deployed unweighted MagAutoTuner."""
    min_samples: int
    min_window_s: float
    service_window_s: float
    true_field_norm_lower: float
    true_field_norm_upper: float
    measurement_residual_norm: float
    true_horizontal_lower: float
    tilt_error_upper_rad: float
    max_norm_ratio_from_mean: float
    min_horizontal_fraction: float
    def __post_init__(self) -> None:
        if self.min_samples < 1: raise ValueError("positive sample count required")
        vals=(self.min_window_s,self.service_window_s,
              self.true_field_norm_lower,self.true_field_norm_upper,
              self.measurement_residual_norm,self.true_horizontal_lower,
              self.tilt_error_upper_rad,self.max_norm_ratio_from_mean,
              self.min_horizontal_fraction)
        if not all(math.isfinite(x) for x in vals): raise ValueError("finite refinement premises required")
        if self.min_window_s < 0 or self.service_window_s <= 0:
            raise ValueError("positive service-window timing required")
        if self.true_field_norm_lower <= self.measurement_residual_norm or self.true_field_norm_upper < self.true_field_norm_lower:
            raise ValueError("field must dominate residual")
        if self.measurement_residual_norm < 0 or self.true_horizontal_lower <= 0 or self.tilt_error_upper_rad < 0:
            raise ValueError("valid field/residual/tilt bounds required")
        if self.max_norm_ratio_from_mean < 0 or not 0 < self.min_horizontal_fraction < 1:
            raise ValueError("valid tuner gates required")

def refinement_gate_margins(p: RefinementPremises) -> dict:
    """Derive literal MagAutoTuner gates from physical field/residual bounds.

    Rotation preserves the true field norm. With ||r||<=R, every corrected
    sample norm lies in [B_min-R,B_max+R]. For one physical field magnitude the
    sharper running-norm variation is 2R/(B_min-R), independent of attitude.
    A tilt-frame error eps changes a vector by at most 2 B_max sin(eps/2);
    adding the measurement residual gives a conservative horizontal-mean loss.
    """
    norm_ratio=2.0*p.measurement_residual_norm/(p.true_field_norm_lower-p.measurement_residual_norm)
    rotation_loss=2.0*p.true_field_norm_upper*math.sin(min(math.pi,p.tilt_error_upper_rad)/2.0)
    horizontal_lower=p.true_horizontal_lower-rotation_loss-p.measurement_residual_norm
    sample_norm_upper=p.true_field_norm_upper+p.measurement_residual_norm
    horizontal_fraction=horizontal_lower/sample_norm_upper
    return {"norm_ratio_upper":norm_ratio,
            "norm_ratio_margin":p.max_norm_ratio_from_mean-norm_ratio,
            "horizontal_mean_lower":horizontal_lower,
            "horizontal_fraction_lower":horizontal_fraction,
            "horizontal_fraction_margin":horizontal_fraction-p.min_horizontal_fraction}

def refinement_tilt_limit_rad(p: RefinementPremises) -> float:
    """Largest tilt-frame error allowed by the literal horizontal gate."""
    required=p.min_horizontal_fraction*(p.true_field_norm_upper+p.measurement_residual_norm)
    budget=p.true_horizontal_lower-p.measurement_residual_norm-required
    if budget <= 0.0: return 0.0
    x=min(1.0,budget/(2.0*p.true_field_norm_upper))
    return 2.0*math.asin(x)

def refinement_sample_gate_uniform(p: RefinementPremises) -> bool:
    m=refinement_gate_margins(p)
    return m["norm_ratio_margin"] >= 0.0 and m["horizontal_fraction_margin"] >= 0.0

def refinement_completion_bound(start_s: float, p: RefinementPremises) -> float:
    """Completion bound from recurring MAGNETIC SERVICE.

    Once the physical margins make every service event acceptable to the
    unweighted tuner, at least one such event occurs per T_M. The count reaches
    N in at most N*T_M. Because accepted-window time telescopes between accepted
    callbacks, the elapsed-window gate reaches W in at most W+T_M.
    """
    if not math.isfinite(start_s) or start_s < 0: raise ValueError("finite nonnegative start required")
    if not refinement_sample_gate_uniform(p):
        raise ValueError("physical bounds do not guarantee MagAutoTuner acceptance")
    wait=max(p.min_samples*p.service_window_s,p.min_window_s+p.service_window_s)
    return start_s+wait

def a21_entry_bound(capture_s: float, refinement_not_before_s: float,
                    refinement: RefinementPremises,
                    accel_bias_unlock_updates: int, unlock_guard_s: float) -> float:
    """Conservative finite time to active-bias A21 after captured service starts."""
    if not all(math.isfinite(x) for x in (capture_s,refinement_not_before_s,unlock_guard_s)):
        raise ValueError("finite entry times required")
    if min(capture_s,refinement_not_before_s,unlock_guard_s) < 0 or accel_bias_unlock_updates < 0:
        raise ValueError("nonnegative entry data required")
    refine_start=max(capture_s,refinement_not_before_s)
    refine_done=refinement_completion_bound(refine_start,refinement)
    # Worst case assumes none of the internal accepted-update count was earned
    # before the outer reference hold is released.
    return release_bound_from_service(accel_bias_unlock_updates,
                                      refinement.service_window_s,
                                      unlock_guard_s,refine_done)

def bias_projection_inactive_margin(projection_radius: float, true_bias_bound: float) -> float:
    """Error radius that keeps the A21 estimate strictly inside projection."""
    if not all(math.isfinite(x) for x in (projection_radius,true_bias_bound)):
        raise ValueError("finite bias radii required")
    if projection_radius <= true_bias_bound or true_bias_bound < 0:
        raise ValueError("projection radius must exceed physical bias bound")
    return projection_radius-true_bias_bound

def information_contraction_ratio(info_floor: float) -> float:
    """Covariance-metric contraction from a normalized full-word information floor.

    In root-whitened coordinates, prediction with Q>=0 is nonexpansive because
    P- = F P F' + Q. A linear Kalman correction is likewise nonexpansive in its
    updated covariance metric. If the complete transported information over a
    word is bounded below by mu*I in the normalized root coordinates, the
    information-form comparison gives rho0 <= 1/(1+mu).
    """
    if not math.isfinite(info_floor) or info_floor <= 0:
        raise ValueError("positive finite information floor required")
    return 1.0/(1.0+info_floor)

def nonlinear_margin_from_information(info_floor: float) -> float:
    """Largest remainder Lipschitz gain compatible with strict small gain."""
    return 1.0-math.sqrt(information_contraction_ratio(info_floor))

def translational_observability_determinant(dt_s: float, tau_s: float) -> float:
    """Determinant of a five-row scalar OU-III observability minor.

    State order is (v,p,S,a_w,b_a). Use four consecutive S observations
    H_S F^j, j=0..3, and one accelerometer row after attitude normalization,
    H_a=(0,0,0,1,1). For the exact OU-III discrete primitives the determinant
    simplifies to -dt^3*tau^3*(1-exp(-dt/tau))^3. It is nonzero for every
    finite dt,tau>0. Thus each translational/bias axis is structurally
    observable; compact positive dt/tau bounds turn this into a uniform
    nonsingularity margin once measurement weights are bounded.
    """
    if not all(math.isfinite(x) for x in (dt_s,tau_s)) or dt_s <= 0 or tau_s <= 0:
        raise ValueError("positive finite dt and tau required")
    phi=math.exp(-dt_s/tau_s)
    return -(dt_s**3)*(tau_s**3)*(1.0-phi)**3

def translational_observability_nonsingular(dt_s: float, tau_s: float) -> bool:
    return translational_observability_determinant(dt_s,tau_s) != 0.0

@dataclass(frozen=True)
class ShippingScheduleBounds:
    """Compact A21 scheduler/tuner envelope used by the information proof."""
    dt_min_s: float
    dt_max_s: float
    tau_min_s: float
    tau_max_s: float
    pseudo_ratio: float
    pseudo_min_s: float
    pseudo_max_s: float
    def __post_init__(self) -> None:
        v=(self.dt_min_s,self.dt_max_s,self.tau_min_s,self.tau_max_s,
           self.pseudo_ratio,self.pseudo_min_s,self.pseudo_max_s)
        if not all(math.isfinite(x) for x in v): raise ValueError("finite schedule bounds required")
        if not 0 < self.dt_min_s <= self.dt_max_s: raise ValueError("ordered sample bounds required")
        if not 0 < self.tau_min_s <= self.tau_max_s: raise ValueError("ordered tau bounds required")
        if self.pseudo_ratio <= 0 or not 0 < self.pseudo_min_s <= self.pseudo_max_s:
            raise ValueError("valid pseudo cadence required")

def pseudo_period_for_tau(tau_s: float, s: ShippingScheduleBounds) -> float:
    if not math.isfinite(tau_s) or not s.tau_min_s <= tau_s <= s.tau_max_s:
        raise ValueError("tau outside certified schedule")
    return min(max(s.pseudo_ratio*tau_s,s.pseudo_min_s),s.pseudo_max_s)

def pseudo_decay_exponent_per_gap(s: ShippingScheduleBounds) -> float:
    """Uniform integral(dt/tau) bound between successive S updates.

    The shipping cadence is retargeted progress-preservingly whenever tau is
    committed. On each constant-tau segment a due event occurs no later than
    one selected pseudo period plus one IMU sample of scheduler overshoot.
    The ratio (period(tau)+dt_max)/tau is monotone on each clamp branch, so its
    maximum is attained at a clamp endpoint or tau_min/tau_max.
    """
    candidates={s.tau_min_s,s.tau_max_s,
                s.pseudo_min_s/s.pseudo_ratio,
                s.pseudo_max_s/s.pseudo_ratio}
    vals=[]
    for tau in candidates:
        tau=min(max(tau,s.tau_min_s),s.tau_max_s)
        vals.append((pseudo_period_for_tau(tau,s)+s.dt_max_s)/tau)
    return max(vals)

def time_varying_translation_minor_floor(s: ShippingScheduleBounds) -> float:
    """Strict determinant floor for the four-S-row time-varying OU minor.

    For arbitrary measurable positive lambda(t)=1/tau(t), the initial-a_w
    column of S(t) is K(t) with K'''(t)=exp(-integral lambda)>0. Thus
    {1,t,t^2/2,K(t)} is an extended complete Chebyshev system. The generalized
    mean-value formula gives determinant = Vandermonde*K'''(xi)/12.
    Four successive S events are separated by at least one IMU sample, so the
    Vandermonde/12 is at least dt_min^6. The cadence/tau coupling above bounds
    the decay over the three event gaps.
    """
    decay=math.exp(-3.0*pseudo_decay_exponent_per_gap(s))
    return (s.dt_min_s**6)*decay

def gravity_tilt_sensitivity_floor(gravity_mps2: float, acceleration_bound_mps2: float) -> float:
    """Lower nonzero singular value of -[R(a_w-g)]x on the captured domain."""
    if not all(math.isfinite(x) for x in (gravity_mps2,acceleration_bound_mps2)):
        raise ValueError("finite gravity/acceleration bounds required")
    if gravity_mps2 <= acceleration_bound_mps2 or acceleration_bound_mps2 < 0:
        raise ValueError("gravity must dominate admitted inertial acceleration")
    return gravity_mps2-acceleration_bound_mps2

def coupled_information_floor(block_a: float, block_b: float, cross_norm: float) -> float:
    """Lower eigenvalue from two certified information blocks and cross norm."""
    if not all(math.isfinite(x) for x in (block_a,block_b,cross_norm)):
        raise ValueError("finite information bounds required")
    if block_a <= 0 or block_b <= 0 or cross_norm < 0:
        raise ValueError("positive block floors and nonnegative cross bound required")
    return 0.5*(block_a+block_b-math.sqrt((block_a-block_b)**2+4.0*cross_norm**2))

def active_tilt_bias_minor(gravity_sensitivity: float, spacing_s: float,
                           bias_tau_s: float) -> float:
    """Three-row tilt/gyro-bias/accelerometer-bias separation minor.

    After the translational columns are eliminated, one fixed tilt axis has
    rows [g,0,1], [g,-g*T,phi,], [g,-2g*T,phi^2] in the worst constant-frame
    case, phi=exp(-T/tau_b). Its determinant is
    -T*g^2*(1-phi)^2. Rotation of the measured gravity direction is handled by
    the compact full-word argument; this minor excludes the stationary worst
    case from carrying a hidden tilt/bias mode.
    """
    if not all(math.isfinite(x) for x in (gravity_sensitivity,spacing_s,bias_tau_s)):
        raise ValueError("finite tilt/bias data required")
    if gravity_sensitivity <= 0 or spacing_s <= 0 or bias_tau_s <= 0:
        raise ValueError("positive tilt/bias data required")
    phi=math.exp(-spacing_s/bias_tau_s)
    return -spacing_s*(gravity_sensitivity**2)*((1.0-phi)**2)

def compact_uniform_information_exists(*, translation_minor_floor: float,
                                       gravity_tilt_floor: float,
                                       magnetic_information_floor: float,
                                       active_bias_minor_abs_floor: float,
                                       strict_branch_guard_margin: float) -> bool:
    """Compactness lemma for the finite A21 service word.

    With dt bounded away from zero, a finite service horizon contains a bounded
    number of hybrid operations. On each strict-guard branch cell the literal
    shipping transition and measurement maps are continuous in the compact
    tuner/sensor/state parameters. Positive structural minors for translation,
    tilt/active bias and magnetic heading exclude a zero-output direction on
    every cell. The normalized information Gramian is therefore positive
    definite pointwise. A finite union of compact cells has a positive minimum
    eigenvalue mu_A21>0.
    """
    vals=(translation_minor_floor,gravity_tilt_floor,magnetic_information_floor,
          active_bias_minor_abs_floor,strict_branch_guard_margin)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite compactness margins required")
    return all(x > 0.0 for x in vals)

def local_nonlinear_radius_exists(info_floor: float, derivative_remainder_at_zero: float = 0.0) -> bool:
    """Existence of a finite-error neighborhood satisfying the A21 small gain.

    Smooth MEKF/reset/tuner maps have derivative remainder tending to its value
    at the origin. Floating-point roundoff is not put in this multiplicative
    remainder; it is retained as additive supply. Thus a local radius exists
    whenever the limiting nonlinear derivative remainder is below the strict
    linear margin.
    """
    if not all(math.isfinite(x) for x in (info_floor,derivative_remainder_at_zero)):
        raise ValueError("finite local data required")
    if info_floor <= 0 or derivative_remainder_at_zero < 0:
        raise ValueError("positive information and nonnegative remainder required")
    return derivative_remainder_at_zero < nonlinear_margin_from_information(info_floor)

def tilt_gyro_bias_quotient_minor(gravity_sensitivity: float, spacing_s: float) -> float:
    """Two-row neutral-quotient minor for one tilt/gyro-bias axis.

    Accelerometer bias is excluded from this quotient because its active A21 OU
    predictor is strictly stable and its physical mismatch belongs to supply.
    After translational columns are eliminated, two gravity observations of
    theta(t)=theta0-t*b_g give determinant -g^2*T.
    """
    if not all(math.isfinite(x) for x in (gravity_sensitivity,spacing_s)):
        raise ValueError("finite quotient data required")
    if gravity_sensitivity <= 0 or spacing_s <= 0:
        raise ValueError("positive quotient data required")
    return -(gravity_sensitivity**2)*spacing_s

def neutral_quotient_uniform_mu_exists(*, translation_minor_floor: float,
                                       gravity_tilt_floor: float,
                                       tilt_gyro_minor_abs_floor: float,
                                       magnetic_information_floor: float,
                                       strict_guard_margin: float) -> bool:
    """Uniform observability of the non-decaying A21 quotient.

    The quotient contains attitude, gyro bias and (v,p,S,a_w); active
    accelerometer bias is removed because phi_b<1. Translation is uniformly
    observable from recurring S rows, tilt/gyro-bias from recurring gravity
    rows after translation elimination, and heading/axial gyro-bias from
    MAGNETIC SERVICE. Strict guard margins make each finite hybrid branch cell
    compact and continuous. Pointwise full rank on the finite union therefore
    gives a uniform normalized quotient Gramian floor mu_N>0.
    """
    vals=(translation_minor_floor,gravity_tilt_floor,tilt_gyro_minor_abs_floor,
          magnetic_information_floor,strict_guard_margin)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite quotient margins required")
    return all(x > 0.0 for x in vals)

def active_bias_homogeneous_ratio(word_s: float, bias_tau_s: float) -> float:
    """Squared homogeneous A21 accelerometer-bias ratio over a word."""
    if not all(math.isfinite(x) for x in (word_s,bias_tau_s)) or word_s <= 0 or bias_tau_s <= 0:
        raise ValueError("positive finite bias word required")
    return math.exp(-2.0*word_s/bias_tau_s)

def detectable_tail_margin(quotient_information_floor: float,
                           word_s: float, bias_tau_s: float) -> float:
    """Strict norm margin from observable neutral quotient plus stable bias mode."""
    q_obs=math.sqrt(information_contraction_ratio(quotient_information_floor))
    q_b=math.sqrt(active_bias_homogeneous_ratio(word_s,bias_tau_s))
    return 1.0-max(q_obs,q_b)

def uniform_controllability_exists(*, dt_min_s: float, tau_min_s: float,
                                   tau_max_s: float, aw_sigma_floor: float,
                                   gyro_noise_floor: float,
                                   gyro_bias_rw_floor: float,
                                   accel_bias_drive_floor: float) -> bool:
    """Compact UCC existence for the shipping A21 prediction.

    Integrated-OU process covariance is the finite-horizon controllability
    Gramian of a controllable chain and is SPD for dt,tau,sigma>0. The exact
    attitude/gyro-bias discretization is driven by positive gyro white noise
    and gyro-bias RW density. Active accelerometer bias has positive OU driving
    density. Continuity on the compact timing/tau envelope gives a uniform
    positive controllability floor.
    """
    vals=(dt_min_s,tau_min_s,tau_max_s,aw_sigma_floor,gyro_noise_floor,
          gyro_bias_rw_floor,accel_bias_drive_floor)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite UCC bounds required")
    return (dt_min_s > 0 and 0 < tau_min_s <= tau_max_s and
            min(aw_sigma_floor,gyro_noise_floor,gyro_bias_rw_floor,
                accel_bias_drive_floor) > 0)

def linear_kalman_tail_exponentially_stable(*, uniform_detectability: bool,
                                            uniform_controllability: bool,
                                            covariance_hard_events_retained: bool) -> bool:
    """Fail-closed composition for the LTV A21 Kalman error transition.

    Uniform detectability/observability bounds the Riccati covariance above;
    uniform complete controllability bounds it below. The covariance metric is
    then uniformly coercive, prediction plus Joseph corrections are
    nonexpansive, and the positive process/measurement dissipation over each
    service word gives a source-uniform rho0<1. Discrete covariance sync/release
    events must separately preserve the same compact bounds.
    """
    for x in (uniform_detectability,uniform_controllability,covariance_hard_events_retained):
        if type(x) is not bool: raise ValueError("literal proof flags required")
    return uniform_detectability and uniform_controllability and covariance_hard_events_retained

def smooth_a21_remainder_vanishes_locally(*, projection_inactive: bool,
                                          tuner_exogenous_same_history: bool,
                                          magnetic_reference_exogenous_same_history: bool,
                                          finite_operation_domain: bool) -> bool:
    """Whether the multiplicative nonlinear remainder eta(r) tends to zero.

    On the strict inner domain, quaternion exp/normalization, measurement maps
    and reset Jacobians are smooth. The tuner and post-refinement magnetic
    reference are measurement-only/exogenous, hence identical in a same-history
    error comparison and belong to the LTV schedule rather than the nonlinear
    remainder. Floating-point roundoff is additive supply, not eta.
    """
    flags=(projection_inactive,tuner_exogenous_same_history,
           magnetic_reference_exogenous_same_history,finite_operation_domain)
    if any(type(x) is not bool for x in flags): raise ValueError("literal proof flags required")
    return all(flags)

def nonlinear_small_gain_exists_from_linear(linear_rho0: float,
                                            smooth_remainder_vanishes: bool) -> bool:
    """Existential local closure: eta(r)->0 and rho0<1 imply some r*>0."""
    if not math.isfinite(linear_rho0) or not 0 <= linear_rho0 < 1:
        raise ValueError("strict finite linear ratio required")
    if type(smooth_remainder_vanishes) is not bool:
        raise ValueError("literal smoothness flag required")
    return smooth_remainder_vanishes

def shipping_covariance_hard_events_retain_compactness(*,
        aw_stationary_std_floor: float, aw_stationary_std_ceiling: float,
        release_bias_variance_ceiling: float) -> bool:
    """Retention lemma for the two A21 covariance reconfiguration events.

    Default periodic a_w synchronization queues a PSD increment inside the next
    prediction, bounded by the stationary covariance target. Bias release only
    raises the three P_ba diagonal entries to a finite configured initial
    variance; hold previously zeroed its cross blocks. Neither operation
    creates an unbounded covariance or destroys positive definiteness.
    """
    vals=(aw_stationary_std_floor,aw_stationary_std_ceiling,release_bias_variance_ceiling)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite hard-event bounds required")
    return (0 < aw_stationary_std_floor <= aw_stationary_std_ceiling and
            release_bias_variance_ceiling > 0)

def attitude_gyro_quotient_uniformly_observable(*, gravity_floor: float,
                                                magnetic_floor: float,
                                                sample_spacing_floor: float,
                                                angular_rate_ceiling: float) -> bool:
    """Uniform attitude/gyro-bias observability on the neutral quotient.

    Known body rotation contributes skew/orthogonal transport and cannot change
    singular values. Gravity supplies a rank-two attitude sensitivity with
    nonzero singular values >=gravity_floor; two separated observations expose
    the corresponding gyro-bias components. The sole instantaneous gravity
    null direction is heading/axial gyro bias, exactly the two-coordinate
    subspace controlled by MAGNETIC SERVICE. A changing gravity direction can
    add information; the constant-direction case is the rank-minimal case.
    Compact bounded angular transport and positive sample spacing make the
    resulting quotient Gramian floor uniform.
    """
    vals=(gravity_floor,magnetic_floor,sample_spacing_floor,angular_rate_ceiling)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite attitude bounds required")
    return gravity_floor > 0 and magnetic_floor > 0 and sample_spacing_floor > 0 and angular_rate_ceiling >= 0

def explicit_small_gain_radius(linear_rho0: float, quadratic_remainder_slope: float,
                               guard_radius: float) -> float:
    """Constructive radius when eta(r)<=L2*r on the strict inner domain."""
    vals=(linear_rho0,quadratic_remainder_slope,guard_radius)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite constructive bounds required")
    if not 0 <= linear_rho0 < 1 or quadratic_remainder_slope <= 0 or guard_radius <= 0:
        raise ValueError("strict linear margin and positive remainder/domain bounds required")
    margin=1.0-math.sqrt(linear_rho0)
    return min(guard_radius,margin/quadratic_remainder_slope)

def additive_supply_practical_radius(linear_rho0: float, nonlinear_eta: float,
                                     coercivity_lower: float,
                                     additive_supply_energy: float) -> float:
    """Practical radius after multiplicative small gain and additive arithmetic/physical supply."""
    vals=(linear_rho0,nonlinear_eta,coercivity_lower,additive_supply_energy)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite practical bounds required")
    if linear_rho0 < 0 or nonlinear_eta < 0 or coercivity_lower <= 0 or additive_supply_energy < 0:
        raise ValueError("valid practical bounds required")
    rho=(math.sqrt(linear_rho0)+nonlinear_eta)**2
    if not rho < 1.0: raise ValueError("small-gain condition not satisfied")
    return math.sqrt(additive_supply_energy/(coercivity_lower*(1.0-rho)))

def marine_magnetic_vector_diversity_floor(window_s: float, gravity_mps2: float,
                                           horizontal_field_min: float,
                                           field_norm_max: float,
                                           velocity_bound_mps: float) -> float:
    """Uniform f x B diversity forced by bounded marine velocity.

    With f=a-g and a=dv/dt for one persistent physical history,
      integral_0^T f x B dt = (v(T)-v(0)) x B - T g x B.
    Hence the average cross-product magnitude is at least
      g*B_h,min - 2*V_max*B_max/T.
    A positive value proves that some instant in every T-window has at least
    that much gravity/magnetic vector diversity. This rules out sustained
    accelerometer/magnetometer collinearity without adding an excitation
    assumption: it follows from bounded marine velocity and nonzero horizontal
    geomagnetic field.
    """
    vals=(window_s,gravity_mps2,horizontal_field_min,field_norm_max,velocity_bound_mps)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite diversity bounds required")
    if min(window_s,gravity_mps2,horizontal_field_min,field_norm_max) <= 0 or velocity_bound_mps < 0:
        raise ValueError("positive physical diversity bounds required")
    return gravity_mps2*horizontal_field_min - 2.0*velocity_bound_mps*field_norm_max/window_s

def marine_magnetic_diversity_window_min(gravity_mps2: float, horizontal_field_min: float,
                                         field_norm_max: float, velocity_bound_mps: float) -> float:
    """Strict threshold above which every window has positive vector diversity."""
    vals=(gravity_mps2,horizontal_field_min,field_norm_max,velocity_bound_mps)
    if not all(math.isfinite(x) for x in vals) or min(gravity_mps2,horizontal_field_min,field_norm_max) <= 0 or velocity_bound_mps < 0:
        raise ValueError("valid physical diversity bounds required")
    return 2.0*velocity_bound_mps*field_norm_max/(gravity_mps2*horizontal_field_min)

def attitude_information_from_marine_diversity(*, diversity_floor: float,
                                               magnetic_service_floor: float,
                                               magnetic_service_window_s: float,
                                               angular_rate_ceiling: float) -> bool:
    """Full attitude information over a diversity window.

    A magnetic event from each service interval can be transported through the
    known attitude transition to the diversity instant. Orthogonal attitude
    transport preserves the magnetic sensitivity norm. Positive f x B
    diversity makes the two rank-two vector sensitivities jointly full rank.
    Recurrence over consecutive diversity windows then exposes gyro bias through
    its attitude injection. Bounded angular rate keeps all finite transports
    continuous on the compact branch cells.
    """
    vals=(diversity_floor,magnetic_service_floor,magnetic_service_window_s,angular_rate_ceiling)
    if not all(math.isfinite(x) for x in vals): raise ValueError("finite attitude information bounds required")
    return diversity_floor > 0 and magnetic_service_floor > 0 and magnetic_service_window_s > 0 and angular_rate_ceiling >= 0
