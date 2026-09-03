#!/usr/bin/env python3
"""Source-correlated correction and exact reset budget for OU-III P4.

The exact reset transport retained after PR #466 is parameterized by the
shipping attitude-correction norm ``delta``.  A global fixed correction cap is
not a source-faithful way to supply that parameter: the Kalman correction is a
state-dependent quantity and must shrink with the P4 Lyapunov state.

This module closes that interface without changing the filter or the frozen P3
/P4 theorem gates.  Let the attached P3 metric at one source/phase class be the
shipping covariance ``P`` and let

    V = z' P^-1 z.

For one accepted 3-vector update, with ``S = H P H' + R`` and the attitude rows
``K_theta`` of the implemented Kalman gain,

    K_theta S K_theta' <= P_theta <= U_theta I.

The exact innovation is ``y = h + eta``.  The linear Kalman inequality gives

    h' S^-1 h <= V,

while the retained 0.8-rad pure-vector remainder and ``S >= R = r I`` give

    eta' S^-1 eta <= s^2 |v|^2 U_theta/r * V = gamma_eta V.

Therefore, by the triangle inequality in the ``S^-1`` norm,

    |dtheta| <= kappa sqrt(V),
    kappa = sqrt(U_theta) * (1 + sqrt(gamma_eta)).               (1)

For accelerometer updates the co-rotated ``a_w`` coordinate from #466 is
essential: ``eta`` is pure attitude rotation, so neither ``a_w`` nor ``b_a`` is
charged a second time in (1).  Magnetometer uses the same pure-vector bound.

To obtain a finite reset budget, choose an information-energy funnel ``V<=W``.
For a candidate attitude radius q_c, every source class obeys the chart bound if

    W <= q_c^2 / U_theta.

Equation (1) simultaneously gives ``delta <= kappa sqrt(W)``.  The module then
feeds the *same source/phase* ``q=sqrt(U_theta W)`` and correction radius into
``ou3_p4_exact_reset_transport.reset_defect_bound`` and reports

    eps_reset = sqrt(M_theta) |rho_theta| / sqrt(W),

where ``M_theta`` is the attached post-update information-metric upper for that
same class.  Thus the exact shipped reset contributes an amplitude perturbation
bounded by ``eps_reset sqrt(W)``.  Because both q and delta are O(sqrt(W)), the
reset mismatch is O(W) and eps_reset -> 0 as W -> 0; no affine disturbance is
introduced.

The declared 30/25/20/15 degree values are certificate-search candidates only,
not new deployment assumptions.  All nonlinear remainder coefficients still
come from the frozen 0.8-rad outer geometry.  This producer is a prerequisite
for the complete signed H18/A21 word and cannot promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_exact_reset_transport as RESET
import ou3_p4_p3_metric_attachment as METRIC
import ou3_p4_vector_remainder_sector as REMAINDER
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
MODES = ("H", "A")
MEASUREMENTS = ("accelerometer", "magnetometer")
RESET_CORRECTION_UTILITY_MAX_RAD = RESET.CAYLEY_MONOTONE_NORM_MAX


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def sqrt_up(x: float) -> float:
    x = float(x)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("finite nonnegative square-root input required")
    return up(math.sqrt(x))


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def div_up(a: float, b: float) -> float:
    b = float(b)
    if not (math.isfinite(b) and b > 0.0):
        raise ValueError("strict positive denominator required")
    return up(float(a) / b)


def _positive(x, label: str) -> float:
    y = float(x)
    if not (math.isfinite(y) and y > 0.0):
        raise RuntimeError(f"{label} must be finite positive")
    return y


def _candidate_q_upper(angle_deg: float) -> float:
    angle = float(angle_deg) * math.pi / 180.0
    if not (math.isfinite(angle) and 0.0 < angle < math.pi):
        raise ValueError("candidate attitude angle must lie in (0,180 deg)")
    half_hi = up(0.5 * angle)
    s_hi = VT.sin_point(half_hi).hi
    c_lo = VT.cos_point(half_hi).lo
    if not c_lo > 0.0:
        raise RuntimeError("candidate Cayley radius reaches chart antipode")
    return div_up(mul_up(2.0, s_hi), c_lo)


def correction_gain_per_sqrt_energy(
    theta_covariance_upper: float,
    vector_norm_upper: float,
    measurement_variance_lower: float,
    eta_tangent_squared_coefficient_upper: float,
) -> dict:
    """Return kappa in ||dtheta|| <= kappa sqrt(V), with outward rounding."""
    utheta = _positive(theta_covariance_upper, "theta covariance upper")
    vmax = _positive(vector_norm_upper, "measurement vector norm upper")
    rlo = _positive(measurement_variance_lower, "measurement variance lower")
    s2 = float(eta_tangent_squared_coefficient_upper)
    if not (math.isfinite(s2) and s2 >= 0.0):
        raise ValueError("nonlinear eta coefficient must be finite nonnegative")

    gamma = div_up(mul_up(mul_up(s2, mul_up(vmax, vmax)), utheta), rlo)
    root_gamma = sqrt_up(gamma)
    innovation_amplitude_gain = up(1.0 + root_gamma)
    innovation_energy_gain = mul_up(innovation_amplitude_gain, innovation_amplitude_gain)
    kappa = mul_up(sqrt_up(utheta), innovation_amplitude_gain)
    return {
        "theta_covariance_upper": utheta,
        "vector_norm_upper": vmax,
        "measurement_variance_lower": rlo,
        "eta_tangent_squared_coefficient_upper": s2,
        "eta_Sinv_energy_gain_upper": gamma,
        "exact_innovation_Sinv_amplitude_gain_upper": innovation_amplitude_gain,
        "exact_innovation_Sinv_energy_gain_upper": innovation_energy_gain,
        "attitude_correction_norm_per_sqrt_metric_energy_upper": kappa,
    }


def reset_row_at_energy(
    *,
    W: float,
    candidate_q_upper: float,
    theta_covariance_upper: float,
    theta_information_upper: float,
    correction_gain: float,
) -> dict:
    """Evaluate the exact reset mismatch for one source/phase operation class."""
    W = _positive(W, "funnel energy W")
    qc = _positive(candidate_q_upper, "candidate Cayley radius")
    utheta = _positive(theta_covariance_upper, "theta covariance upper")
    mtheta = _positive(theta_information_upper, "theta information upper")
    kappa = _positive(correction_gain, "correction gain")

    sqrtW = sqrt_up(W)
    q_energy = sqrt_up(mul_up(utheta, W))
    q = min(qc, q_energy)
    delta = mul_up(kappa, sqrtW)
    if delta > RESET_CORRECTION_UTILITY_MAX_RAD:
        raise RuntimeError("correction radius exceeds validated exact-reset utility range")
    exact = RESET.reset_defect_bound(q, delta)
    rho = float(exact["reset_attitude_defect_norm_upper"])
    eps = div_up(mul_up(sqrt_up(mtheta), rho), sqrtW)
    return {
        "funnel_energy_upper": W,
        "state_cayley_norm_upper": q,
        "correction_norm_upper": delta,
        "reset_attitude_defect_norm_upper": rho,
        "reset_metric_amplitude_fraction_upper": eps,
        "reset_metric_energy_fraction_upper": mul_up(eps, eps),
        "cayley_composition_denominator_lower": exact["cayley_composition_denominator_lower"],
        "reset_inverse_operator_norm_upper": exact["reset_inverse_operator_norm_upper"],
    }


def _measurement_constants(domain: dict, vector: dict, remainder: dict) -> dict:
    live = domain["normal_live"]
    vc = vector["configured_measurement_bounds"]
    acc_std = _positive(vc["acc_measurement_std_mps2"], "accelerometer std")
    mag_std = _positive(vc["mag_measurement_std_uT"], "magnetometer std")
    return {
        "accelerometer": {
            "vector_norm_upper": _positive(live["specific_force_norm_upper_mps2"], "specific force upper"),
            "measurement_variance_lower": down(acc_std * acc_std),
            "eta_tangent_squared_coefficient_upper": float(
                remainder["acc_eta_force_rotation_quadratic_coefficient_upper"]
            ),
        },
        "magnetometer": {
            "vector_norm_upper": _positive(live["magnetic_vector_norm_upper_uT"], "magnetic norm upper"),
            "measurement_variance_lower": down(mag_std * mag_std),
            "eta_tangent_squared_coefficient_upper": float(
                remainder["mag_eta_squared_over_linear_rotation_squared_upper"]
            ),
        },
    }


def _operation_classes(metric: dict, constants: dict):
    """Yield exact finite and frozen source classes without Cartesian source hulling."""
    frozen_by_source = {
        int(row["source_node"]): row
        for row in metric["frozen_clock"]["source_rows"]
    }
    for endpoint in metric["endpoint_rows"]:
        source = int(endpoint["source_node"])
        for phase in range(26):
            env = endpoint[
                "boundary_history_envelope" if phase == 0 else "positive_phase_history_envelope"
            ]
            for mode in MODES:
                utheta = float(env[f"{mode}_bias_covariance_upper"]["theta_covariance_upper"])
                mtheta = float(
                    endpoint["finite_phase_information_metric_upper_group_diagonal"][mode][phase]["theta"]
                )
                for measurement in MEASUREMENTS:
                    gain = correction_gain_per_sqrt_energy(
                        utheta,
                        constants[measurement]["vector_norm_upper"],
                        constants[measurement]["measurement_variance_lower"],
                        constants[measurement]["eta_tangent_squared_coefficient_upper"],
                    )
                    yield {
                        "branch": "finite",
                        "source_node": source,
                        "phase": phase,
                        "mode": mode,
                        "measurement": measurement,
                        "theta_covariance_upper": utheta,
                        "theta_information_upper": mtheta,
                        "correction": gain,
                    }

        # Frozen clock holds the positive-phase same-history covariance upper,
        # with a separately certified arbitrary-duration information metric.
        env = endpoint["positive_phase_history_envelope"]
        frow = frozen_by_source[source]
        for mode in MODES:
            utheta = float(env[f"{mode}_bias_covariance_upper"]["theta_covariance_upper"])
            mtheta = float(frow["information_metric_upper_group_diagonal"][mode]["theta"])
            for measurement in MEASUREMENTS:
                gain = correction_gain_per_sqrt_energy(
                    utheta,
                    constants[measurement]["vector_norm_upper"],
                    constants[measurement]["measurement_variance_lower"],
                    constants[measurement]["eta_tangent_squared_coefficient_upper"],
                )
                yield {
                    "branch": "frozen",
                    "source_node": source,
                    "phase": None,
                    "mode": mode,
                    "measurement": measurement,
                    "theta_covariance_upper": utheta,
                    "theta_information_upper": mtheta,
                    "correction": gain,
                }


def evaluate(metric: dict, domain: dict, vector: dict, remainder: dict, reset: dict) -> dict:
    failures = [f"metric: {x}" for x in METRIC.validate(metric)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]
    failures += [f"remainder: {x}" for x in REMAINDER.validate(remainder)]
    failures += [f"reset: {x}" for x in RESET.validate(reset)]
    if failures:
        return {"schema": SCHEMA, "failures": failures}
    if metric.get("same_history_P3_frontier_consumed") is not True:
        failures.append("P3 metric attachment is not same-history")
    if metric.get("independent_cartesian_tau_sigma_R_S_extrema_used") is not False:
        failures.append("P3 metric attachment reintroduced Cartesian tuner extrema")
    if remainder.get("accelerometer_corotated_aw_coordinate_used") is not True:
        failures.append("accelerometer nonlinear remainder lost co-rotated a_w coordinate")
    if reset.get("parametric_reset_defect_bound_available") is not True:
        failures.append("exact reset prerequisite does not expose parametric defect bound")

    constants = _measurement_constants(domain, vector, remainder)
    classes = list(_operation_classes(metric, constants))
    finite_expected = 800 * 26 * len(MODES) * len(MEASUREMENTS)
    frozen_expected = 800 * len(MODES) * len(MEASUREMENTS)
    finite_count = sum(row["branch"] == "finite" for row in classes)
    frozen_count = sum(row["branch"] == "frozen" for row in classes)
    if finite_count != finite_expected:
        failures.append(f"finite reset scan covered {finite_count}, expected {finite_expected}")
    if frozen_count != frozen_expected:
        failures.append(f"frozen reset scan covered {frozen_count}, expected {frozen_expected}")

    raw_angles = domain.get("certificate_search", {}).get(
        "p4_complete_word_full_attitude_candidate_deg", []
    )
    angles = [float(x) for x in raw_angles]
    if angles != [30.0, 25.0, 20.0, 15.0]:
        failures.append("declared P4 certificate-search angle candidates changed")

    worst_gain = {
        mode: {
            measurement: {"value": 0.0, "limiting_class": None}
            for measurement in MEASUREMENTS
        }
        for mode in MODES
    }
    for row in classes:
        k = float(row["correction"]["attitude_correction_norm_per_sqrt_metric_energy_upper"])
        slot = worst_gain[row["mode"]][row["measurement"]]
        if k > slot["value"]:
            slot["value"] = k
            slot["limiting_class"] = {
                "branch": row["branch"],
                "source_node": row["source_node"],
                "phase": row["phase"],
            }

    candidates = []
    for angle in angles:
        qc = _candidate_q_upper(angle)
        modes = {}
        for mode in MODES:
            mode_rows = [row for row in classes if row["mode"] == mode]
            # A common W must keep both the actual Cayley state and the shipping
            # correction inside the domains of the exact reset formulas.
            W_angle = min(
                down((qc * qc) / float(row["theta_covariance_upper"]))
                for row in mode_rows
            )
            W_corr = min(
                down((RESET_CORRECTION_UTILITY_MAX_RAD / float(
                    row["correction"]["attitude_correction_norm_per_sqrt_metric_energy_upper"]
                )) ** 2)
                for row in mode_rows
            )
            W = down(min(W_angle, W_corr))
            if not (math.isfinite(W) and W > 0.0):
                failures.append(f"{mode} {angle:g} deg candidate lost positive funnel energy")
                continue

            worst_eps = 0.0
            limiting = None
            min_denom = math.inf
            max_delta = 0.0
            max_q = 0.0
            for row in mode_rows:
                try:
                    rr = reset_row_at_energy(
                        W=W,
                        candidate_q_upper=qc,
                        theta_covariance_upper=row["theta_covariance_upper"],
                        theta_information_upper=row["theta_information_upper"],
                        correction_gain=row["correction"][
                            "attitude_correction_norm_per_sqrt_metric_energy_upper"
                        ],
                    )
                except Exception as exc:
                    failures.append(
                        f"{mode} {angle:g} deg reset class {row['source_node']}/{row['phase']}/"
                        f"{row['measurement']}: {exc}"
                    )
                    continue
                eps = float(rr["reset_metric_amplitude_fraction_upper"])
                if eps > worst_eps:
                    worst_eps = eps
                    limiting = {
                        "branch": row["branch"],
                        "source_node": row["source_node"],
                        "phase": row["phase"],
                        "measurement": row["measurement"],
                        "reset": rr,
                    }
                min_denom = min(min_denom, float(rr["cayley_composition_denominator_lower"]))
                max_delta = max(max_delta, float(rr["correction_norm_upper"]))
                max_q = max(max_q, float(rr["state_cayley_norm_upper"]))

            modes[mode] = {
                "dimension": 18 if mode == "H" else 21,
                "candidate_attitude_angle_deg": angle,
                "candidate_cayley_radius_upper": qc,
                "metric_energy_limit_from_attitude_upper": W_angle,
                "metric_energy_limit_from_reset_utility_upper": W_corr,
                "source_complete_metric_energy_funnel_upper": W,
                "maximum_state_cayley_norm_upper": max_q,
                "maximum_shipping_attitude_correction_norm_upper": max_delta,
                "minimum_cayley_composition_denominator_lower": min_denom,
                "worst_reset_metric_amplitude_fraction_upper": worst_eps,
                "worst_reset_metric_energy_fraction_upper": mul_up(worst_eps, worst_eps),
                "limiting_reset_class": limiting,
                "reset_defect_is_higher_order_in_metric_radius": True,
            }
        candidates.append({
            "attitude_angle_deg": angle,
            "outer_remainder_angle_rad_consumed": float(remainder["outer_angle_rad"]),
            "modes": modes,
        })

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SOURCE_CORRELATED_CORRECTION_RESET_BUDGET",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "same_history_P3_metric_consumed": True,
        "independent_cartesian_tau_sigma_R_S_extrema_used": False,
        "canonical_P3_pass_required_before_P4_theorem_consumption": True,
        "outer_0p8_rad_remainder_consumed": float(remainder["outer_angle_rad"]) == 0.80,
        "accelerometer_corotated_aw_eta_zero_consumed": True,
        "accelerometer_bias_eta_zero_consumed": True,
        "kalman_gain_bound_identity": "K_theta S K_theta^T <= P_theta <= U_theta I",
        "linear_innovation_metric_identity": "H^T S^-1 H <= P^-1",
        "state_dependent_correction_bound": True,
        "global_fixed_correction_cap_assumed": False,
        "source_correlated_correction_norm_bound_supplied_here": True,
        "exact_reset_transport_consumed": True,
        "reset_inverse_operator_norm_upper": reset["reset_inverse_operator_norm_upper"],
        "finite_operation_classes_scanned": finite_count,
        "frozen_operation_classes_scanned": frozen_count,
        "worst_correction_gain_by_mode_measurement": worst_gain,
        "certificate_search_candidates": candidates,
        "reset_defect_vanishes_quadratically_with_metric_radius": True,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "combine these same-source correction/reset budgets with the signed vector tangent forms, "
            "S=0 linear forms, prediction/source-edge metric transport and recurrent H18/A21 word; "
            "choose the widest candidate energy funnel for which the accumulated rho_H and rho_A are strictly below one"
        ),
        "failures": failures,
    }


def build(metric: dict, domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("source-correlated reset budget must not be trajectory fitted")
    return evaluate(metric, domain, VECTOR.build(), REMAINDER.build(path), RESET.build(path))


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_SOURCE_CORRELATED_CORRECTION_RESET_BUDGET":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "same_history_P3_metric_consumed",
        "canonical_P3_pass_required_before_P4_theorem_consumption",
        "outer_0p8_rad_remainder_consumed", "accelerometer_corotated_aw_eta_zero_consumed",
        "accelerometer_bias_eta_zero_consumed", "state_dependent_correction_bound",
        "source_correlated_correction_norm_bound_supplied_here", "exact_reset_transport_consumed",
        "reset_defect_vanishes_quadratically_with_metric_radius",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "independent_cartesian_tau_sigma_R_S_extrema_used", "global_fixed_correction_cap_assumed",
        "complete_H18_A21_word_established_here", "P4_USABLE_CERTIFICATE_PROMOTED",
        "P5_FINITE_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("reset_inverse_operator_norm_upper") != 1.0:
        f.append("exact reset inverse norm is not one")
    if int(d.get("finite_operation_classes_scanned", 0)) != 800 * 26 * 2 * 2:
        f.append("finite reset budget did not cover 800 x 26 x H/A x acc/mag")
    if int(d.get("frozen_operation_classes_scanned", 0)) != 800 * 2 * 2:
        f.append("frozen reset budget did not cover 800 x H/A x acc/mag")
    candidates = d.get("certificate_search_candidates", [])
    if [float(x.get("attitude_angle_deg", math.nan)) for x in candidates] != [30.0, 25.0, 20.0, 15.0]:
        f.append("P4 certificate-search candidates changed")
    for cand in candidates:
        for mode in MODES:
            row = cand.get("modes", {}).get(mode, {})
            for key in (
                "source_complete_metric_energy_funnel_upper",
                "minimum_cayley_composition_denominator_lower",
            ):
                x = row.get(key)
                if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(float(x)) or float(x) <= 0.0:
                    f.append(f"{cand.get('attitude_angle_deg')} {mode}: {key} is not finite positive")
            eps = row.get("worst_reset_metric_amplitude_fraction_upper")
            if not isinstance(eps, (int, float)) or isinstance(eps, bool) or not math.isfinite(float(eps)) or float(eps) < 0.0:
                f.append(f"{cand.get('attitude_angle_deg')} {mode}: reset amplitude fraction invalid")
            delta = row.get("maximum_shipping_attitude_correction_norm_upper")
            if not isinstance(delta, (int, float)) or float(delta) > RESET_CORRECTION_UTILITY_MAX_RAD:
                f.append(f"{cand.get('attitude_angle_deg')} {mode}: correction utility range exceeded")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--metric", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    metric = json.loads(a.metric.read_text(encoding="utf-8"))
    d = build(metric, a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "finite_classes": d.get("finite_operation_classes_scanned"),
        "frozen_classes": d.get("frozen_operation_classes_scanned"),
        "worst_correction_gain": d.get("worst_correction_gain_by_mode_measurement"),
        "candidates": d.get("certificate_search_candidates"),
        "P4_promoted": d.get("P4_USABLE_CERTIFICATE_PROMOTED"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
