"""Primary weaker target: compact bias error and full-filter motion ISS.

This builds the changed theorem contract and its exact composition arithmetic.
It does NOT certify an active-mode gain, infer true bias from the clamp, or
transfer a held-mode 18-state certificate to the active 21-state estimator.
"""
from __future__ import annotations

from fractions import Fraction

from ou3_interval import Interval, matrix_mul, matrix_sub, matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_mems_bias_contract as BIAS
import ou3_brmm_response_union as UNION

TARGET = "BOUNDED_BIAS_FULL_FILTER_18_ERROR_PRACTICAL_ISS"
QUALIFICATION = "OU3_BOUNDED_BIAS_MOTION_TARGET_V1"


def bias_error_bound(true_bias_bound, projection_radius=Fraction(2, 5)) -> Fraction:
    """Triangle bound, conditional on a true-bias envelope; no qualification."""
    true, radius = Fraction(str(true_bias_bound)), Fraction(str(projection_radius))
    if true < 0 or radius <= 0:
        raise ValueError("nonnegative true bound and positive projection radius required")
    return true + radius


def compose_bounds(*, rho, Gamma, true_bias_bound, duration, gamma, kappa,
                   source_energy, noise_energy, endpoint_level, chart_level) -> dict:
    """Exact conditional gain/retention arithmetic, not a gain certificate.

    gamma/kappa are (bias, source, noise). A mathematical example passed to
    this function cannot populate the missing source-uniform gain fields.
    """
    rho, Gamma, duration, source_energy, noise_energy, level, chart = map(
        lambda x: Fraction(str(x)),
        (rho, Gamma, duration, source_energy, noise_energy, endpoint_level, chart_level))
    gamma, kappa = [tuple(Fraction(str(x)) for x in g) for g in (gamma, kappa)]
    if len(gamma) != 3 or len(kappa) != 3:
        raise ValueError("three bias/source/noise gains required")
    if not (0 < rho < 1 and Gamma >= 1 and duration > 0 and level >= 0 and chart > 0):
        raise ValueError("invalid contraction, prefix gain, duration or levels")
    if any(x < 0 for x in (*gamma, *kappa, source_energy, noise_energy)):
        raise ValueError("gains and hard energies must be nonnegative")
    bound = bias_error_bound(true_bias_bound)
    energies = (duration*bound*bound, source_energy, noise_energy)
    forcing = sum(g*x for g, x in zip(gamma, energies))
    prefix_forcing = sum(g*x for g, x in zip(kappa, energies))
    endpoint_floor = forcing/(1-rho)
    return {
        "bias_error_bound": bound,
        "bias_window_energy_bound": energies[0],
        "endpoint_forcing_bound": forcing,
        "prefix_forcing_bound": prefix_forcing,
        "endpoint_ultimate_bound": endpoint_floor,
        "all_prefix_ultimate_bound": Gamma*endpoint_floor+prefix_forcing,
        "endpoint_level_invariant": rho*level+forcing <= level,
        "strict_prefix_chart_retention": Gamma*level+prefix_forcing < chart,
        "actual_motion_gains_certified_by_this_arithmetic": False,
    }


def motion_master(C0, CN, M0, MN, factor, weighted_energy_forms, *, prefix=False):
    """Outward L=CN'MN CN-factor C0'M0 C0-sum gain_i Q_i.

    factor is rho<1 at the endpoint, or finite Gamma>=1 at a prefix.

    C0/CN map ONE complete shipping graph to the 18 performance coordinates;
    they do not implement an 18-state Riccati recursion. Q_i are same-graph
    energy forms, not independently generated bias/event histories. Graph
    provenance, PSD energy weights, metrics and sectors remain caller proof
    obligations. The existing full augmented S-procedure/LDLT consumes L.
    """
    def shape(matrix):
        if not matrix or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix):
            raise ValueError("nonempty rectangular matrices required")
        return len(matrix), len(matrix[0])

    rows, size = shape(C0)
    if rows != 18 or size < 21 or shape(CN) != (18, size) or shape(M0) != (18, 18) or shape(MN) != (18, 18):
        raise ValueError("both storage maps must retain all 18 motion coordinates")
    if not isinstance(factor, Interval) or not (0 < factor.lo <= factor.hi < float("inf")):
        raise ValueError("finite positive outward factor required")
    if (prefix and factor.lo < 1) or (not prefix and factor.hi >= 1):
        raise ValueError("endpoint requires contraction; prefix requires finite Gamma >= 1")
    terminal = matrix_mul(matrix_transpose(CN), matrix_mul(MN, CN))
    initial = matrix_mul(matrix_transpose(C0), matrix_mul(M0, C0))
    out = matrix_sub(terminal, [[factor*x for x in row] for row in initial])
    if not weighted_energy_forms:
        raise ValueError("explicit bias/source/noise energy forms required")
    for gain, form in weighted_energy_forms:
        if not isinstance(gain, Interval) or not (0 <= gain.lo <= gain.hi < float("inf")):
            raise ValueError("finite nonnegative outward gain required")
        if shape(form) != (size, size):
            raise ValueError("energy form detached from common graph dimension")
        out = matrix_sub(out, [[gain*x for x in row] for row in form])
    return matrix_symmetric_hull(out)


def build(*, p3_contract: dict) -> dict:
    bias = p3_contract["mems_bias_preconditions"]
    radius = bias["BIAS0"]["filter_projection_radius_mps2"]
    coverage = {obligation: UNION.branch_mode_coverage({"H18": False, "A21": False})
                for obligation in ("motion_endpoint", "bias_to_motion_gain", "every_event_prefix_gain",
                                   "every_event_retention", "projected_bias_domain_source_coverage")}
    return {
        "qualification": QUALIFICATION,
        "primary_target": TARGET,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "response_branches": list(UNION.BRANCHES),
        "source_domain_change": "replace 0.35 estimated-bias interior by closed deployed projection ball for motion target only",
        "projection_radius_mps2": radius,
        "legacy_interior_domain_and_full_state_flags_unchanged": True,
        "motion_storage_dimension": 18,
        "executed_active_filter_dimension": 21,
        "full_active_covariance_cross_terms_gains_resets_and_actual_RS_retained": True,
        "held_H18_certificate_transferred_to_active_filter": False,
        "bias_error_convergence_required": False,
        "full_21_state_nonlinear_contraction_required": False,
        "bias_coupling_is_external_sensor_noise": False,
        "sensor_noise_zero_implies_model_mismatch_zero": False,
        "zero_error_floor_required": False,
        "empirical_accuracy_is_a_certified_theorem_constant": False,
        "same_history_bias_prediction_corrections_and_projection_retained": True,
        "estimate_ball_invariance_lemma_closed_conditionally": True,
        "bias_error_bound_lemma_closed_conditionally": True,
        "bias_error_bound_formula": "B_true + R_b",
        "conditional_true_bias_bound_mps2": None,
        "clamp_or_OU_PSD_supplies_true_bias_bound": False,
        "BIAS2_optional_gain_sharpening_requires_valid_sector": True,
        "BIAS2_positive_separation_required_for_compactness": False,
        "mems_bias_preconditions": bias,
        "P3_CONDITIONAL_BRMM_PASS_consumed": p3_contract["P3_CONDITIONAL_BRMM_PASS"],
        "P3_delta_consumed": {mode: p3_contract["modes"][mode]["relative_Riccati_injection_margin_lower"]
                              for mode in ("H18", "A21")},
        "existing_P3_automatically_covers_projection_boundary": False,
        "hardware_qualification_required_to_run_conditional_mathematics": False,
        "same_history_finite_error_source_forcing_attachment_required": True,
        "prefix_contraction_required": False,
        "all_completed_shipping_event_prefixes_required": True,
        "preprojection_auxiliary_coordinate_may_leave_ball": True,
        "motion_master_assembler_available": True,
        "full_state_to_motion_storage_is_estimator_reduction": False,
        "response_mode_coverage": coverage,
        "certified_motion_gains": None,
        "certified_accuracy_floor": None,
        "conditional_composition_arithmetic_closed": True,
        "P4_MOTION_PASS": False,
        "P5_MOTION_MAY_START": False,
        "P4_MOTION_FAIL_REASONS": [
            "same-history source/error/forcing attachment and projected-bias-domain source/P3 coverage remain open",
            "source-uniform actual-filter 18-error endpoint and bias/source/noise gains remain unproved",
            "every-event finite prefix gains, chart retention, useful accuracy floor and hybrid attachment remain unproved",
        ],
    }


def validate(d: dict) -> list[str]:
    bias = d.get("mems_bias_preconditions", {})
    failures = [f"BIAS: {x}" for x in BIAS.validate(bias)]
    # Freeze the current evidentiary status: arbitrary true booleans are not
    # certificates. A future promoting builder must consume actual witnesses.
    expected = build(p3_contract={
        "mems_bias_preconditions": BIAS.build(), "P3_CONDITIONAL_BRMM_PASS": True,
        "modes": {mode: {"relative_Riccati_injection_margin_lower": 1e-18}
                  for mode in ("H18", "A21")}})
    for key, value in expected.items():
        if key == "mems_bias_preconditions":
            continue
        if d.get(key) != value or key not in d:
            failures.append(f"{key} changed or falsely certified")
    return failures
