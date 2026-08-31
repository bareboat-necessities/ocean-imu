#!/usr/bin/env python3
"""Operation-matched finite-angle sector certificate for OU-III P4/P5.

The legacy P4 proof converted every accepted correction into one adverse
full-word defect and then multiplied it by the number of operations.  PR #441
proved that accounting has a hard milliradian ceiling independent of further
constant sharpening.  This producer certifies the structural replacement on a
physical finite-angle domain.

For a fixed vector and Cayley attitude c=q u, let y_R be the exact rotational
residual and h the tangent residual.  The exact Cayley identities give

    y_R^T (y_R-h) = 0,
    ||y_R-h||^2 = q^2/4 ||y_R||^2,
    y_R^T h / ||h||^2 = 4/(4+q^2).

Thus the nonlinear rotational map is strongly monotone on every finite Cayley
ball and its sector constant is explicit.  The shipping magnetometer further
annihilates the radial residual exactly (K_m v=0), while the accelerometer's
nonlinear residual is exactly representable as an effective a_w input because
J_aw=R_wb is orthogonal/full-row-rank.  S=0 has no nonlinear residual.  Joseph
information transport and the covariance reset are exact congruence identities.

The certified sector is chosen to cover a 0.80 rad full-attitude error, larger
than both source-faithful gauged P1 handoff branches.  0.80 rad is not inferred
from replay data; it is a proof-design radius below the q<1 Cayley chart already
used by the P5 source machinery.  This file deliberately does NOT claim the
complete P4 word contraction: the remaining numerical backend must pair each
operation's sector residual with that operation's own information decrease and
carry directional state blocks along source-reachable paths.  What is proved
here is that the nonlinear geometry itself is no longer microscopic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_cayley_eta_geometry as ETA
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_exact_correction_transport as CORR
import ou3_p5_heading_handoff_contract as HEADING

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DESIGN_THETA_RAD = 0.80


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def cayley_from_theta_upper(theta: float) -> float:
    """Outward upper q=2 tan(theta/2), for 0<=theta<pi."""
    theta = float(theta)
    if not (math.isfinite(theta) and 0.0 <= theta < math.pi):
        raise ValueError("finite attitude angle in [0,pi) required")
    return up(2.0 * math.tan(up(0.5 * theta)))


def theta_from_cayley_upper(q: float) -> float:
    """Outward upper theta=2 atan(q/2)."""
    q = float(q)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    return up(2.0 * math.atan(up(0.5 * q)))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4 finite-angle sector must not be trajectory fitted")

    heading = HEADING.build(path)
    eta = ETA.build(path)
    veff = VEFF.build(path)
    corr = CORR.build(path)
    failures = [f"heading: {x}" for x in HEADING.validate(heading)]
    failures += [f"eta: {x}" for x in ETA.validate(eta)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"correction-transport: {x}" for x in CORR.validate(corr)]

    q = cayley_from_theta_upper(DESIGN_THETA_RAD)
    if not q < 1.0:
        failures.append("0.80 rad design sector does not stay inside q<1 chart")

    residual_factor_lower = down(4.0 / up(4.0 + up(q * q)))
    eta_to_residual_info_upper = up(up(q * q) / 4.0)
    tangent_defect_ratio_upper = up(q / down(math.sqrt(down(4.0 + down(q*q)))))
    exact_residual_to_tangent_norm_lower = down(2.0 / up(math.sqrt(up(4.0 + up(q*q)))))

    normal_q = float(heading["gauged_quality_handoff"]["full_attitude_cayley_norm_upper"])
    timeout_q = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    normal_theta = theta_from_cayley_upper(normal_q)
    timeout_theta = theta_from_cayley_upper(timeout_q)
    covers_normal = normal_q <= q
    covers_timeout = timeout_q <= q

    # These are intentionally useful margins, not epsilon-level positivity.
    # q<1 alone implies eta/residual information ratio<1/4 and vector strong
    # monotonicity>4/5.  Pin stronger numerical floors for this 0.80-rad sector
    # so a future accidental return to a microscopic design radius is visible.
    if not residual_factor_lower > 0.80:
        failures.append("finite-angle vector monotonicity fell below 0.80")
    if not eta_to_residual_info_upper < 0.25:
        failures.append("finite-angle eta information ratio reached 0.25")
    if not tangent_defect_ratio_upper < 0.45:
        failures.append("finite-angle tangent defect ratio reached 0.45")
    if not (covers_normal and covers_timeout):
        failures.append("0.80 rad sector does not cover both gauged P1 handoffs")

    mag = veff.get("magnetometer", {})
    acc = veff.get("accelerometer", {})
    if mag.get("kalman_gain_radial_action_exact_zero") is not True:
        failures.append("magnetometer radial finite-angle residual is not annihilated")
    if acc.get("standalone_eta_information_penalty_required_for_state_correction") is not False:
        failures.append("accelerometer eta was not reduced to effective a_w input")
    if corr.get("condition_number_multiplier_used_for_reset_transport") is not False:
        failures.append("reset transport still uses a condition-number multiplier")

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "design_full_attitude_angle_rad": DESIGN_THETA_RAD,
        "design_full_attitude_angle_deg": math.degrees(DESIGN_THETA_RAD),
        "design_cayley_norm_upper": q,
        "design_cayley_chart_q_lt_1": q < 1.0,
        "exact_vector_strong_monotonicity_factor_lower": residual_factor_lower,
        "exact_eta_to_rotational_residual_information_ratio_upper": eta_to_residual_info_upper,
        "exact_residual_to_tangent_norm_ratio_lower": exact_residual_to_tangent_norm_lower,
        "exact_effective_tangent_defect_ratio_upper": tangent_defect_ratio_upper,
        "operation_classes": {
            "S_zero": {
                "nonlinear_measurement_residual": "IDENTICALLY_ZERO",
                "standalone_eta_penalty": False,
            },
            "magnetometer": {
                "radial_residual_gain_action": "EXACTLY_ZERO",
                "effective_coordinate_nonexpansive": True,
                "standalone_eta_penalty": False,
            },
            "accelerometer": {
                "nonlinear_residual_representation": "EXACT_EFFECTIVE_AW_INPUT",
                "J_aw_orthogonal_full_row_rank": True,
                "standalone_eta_penalty": False,
            },
            "joseph_reset": {
                "exact_information_identity": True,
                "exact_reset_congruence": True,
                "reset_inverse_operator_norm_upper": 1.0,
                "condition_number_multiplier_used": False,
            },
        },
        "P1_overlap": {
            "normal_gauged_cayley_norm_upper": normal_q,
            "normal_gauged_angle_upper_rad": normal_theta,
            "normal_gauged_inside_sector": covers_normal,
            "timeout_gauged_cayley_norm_upper": timeout_q,
            "timeout_gauged_angle_upper_rad": timeout_theta,
            "timeout_gauged_inside_sector": covers_timeout,
            "ungauged_timeout_route": heading["ungauged_timeout_subbranch"]["required_route"],
        },
        "global_packet_count_times_lipschitz_defect_used": False,
        "whole_word_weakest_P3_delta_used_as_attitude_sector_margin": False,
        "P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR_CERTIFICATE": "PASS" if passed else "FAIL",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "next_obligation": (
            "pair each accepted operation with its own Joseph information decrease, retain attitude/gyro-bias/translation directional margins, "
            "and compose only source-reachable path cells; do not return to N-times-global-defect accounting"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in ("source_generated_not_trajectory_fit",):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed",
                "global_packet_count_times_lipschitz_defect_used",
                "whole_word_weakest_P3_delta_used_as_attitude_sector_margin",
                "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if d.get("P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR_CERTIFICATE") != "PASS":
        failures.append("operation-matched finite-angle sector did not pass")
    if not float(d.get("design_full_attitude_angle_rad", 0.0)) >= 0.80:
        failures.append("P4 finite-angle sector regressed below 0.80 rad")
    if not float(d.get("design_cayley_norm_upper", math.inf)) < 1.0:
        failures.append("P4 sector left q<1 chart")
    if not float(d.get("exact_vector_strong_monotonicity_factor_lower", 0.0)) > 0.80:
        failures.append("P4 vector monotonicity is not useful")
    if not float(d.get("exact_eta_to_rotational_residual_information_ratio_upper", math.inf)) < 0.25:
        failures.append("P4 eta ratio is not useful")
    overlap = d.get("P1_overlap", {})
    if overlap.get("normal_gauged_inside_sector") is not True:
        failures.append("normal P1 handoff outside finite-angle sector")
    if overlap.get("timeout_gauged_inside_sector") is not True:
        failures.append("gauged timeout P1 handoff outside finite-angle sector")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR_CERTIFICATE"],
        "theta_rad": out["design_full_attitude_angle_rad"],
        "q": out["design_cayley_norm_upper"],
        "monotonicity": out["exact_vector_strong_monotonicity_factor_lower"],
        "eta_ratio": out["exact_eta_to_rotational_residual_information_ratio_upper"],
        "P1_overlap": out["P1_overlap"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
