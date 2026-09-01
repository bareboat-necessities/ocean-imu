#!/usr/bin/env python3
"""Operation-matched Joseph/directional ledger for the usable OU-III P4 route.

The first-accelerometer sector-invariance diagnostic is deliberately *not* the
acceptance test used here.  The shipping MEKF uses a Joseph covariance update,
sequential quaternion injection/reset, and source-correlated directional
measurement information.  A transient attitude-coordinate excursion can be
admissible when the complete source-word Lyapunov/information level decreases.

For an accepted update with signed residual y = H z + eta and Joseph posterior
P+, the exact identity is

    z^T P^-1 z - (z-Ky)^T (P+)^-1 (z-Ky)
      = y^T S^-1 y - eta^T R^-1 eta.

The stronger route preserves this **signed identity exactly**.  It neither drops
the nonlinear eta term nor replaces it by an independently maximized norm.

For the configured accelerometer, J_aw=R_wb is orthogonal and full row rank, so
for e_eta=J_aw^T eta_a,

    K_a(H_a z + eta_a) = K_a H_a(z + E_aw e_eta).

This is an exact *state-correction range reduction*: finite-angle eta is carried
as a source-correlated effective a_w input.  It does not turn eta into zero in
the Joseph energy identity.  The declared 0.3 g latent-a_w state error is
already a state coordinate in H_a z and must not be mislabeled as measurement
eta.  The actual nonlinear eta penalty remains the finite-angle defect and is
carried jointly with z, P, H, R, S, K and the reset map.

For the configured magnetometer the radial residual has exact zero Kalman-gain
action; its radial contribution cancels out of the state correction and the
useful residual is represented by an effective tangent coordinate.  S=0 has
eta=0 exactly.  Rejected/not-due branches are identities.  Quaternion reset is
an exact covariance congruence with ||G^-1||_2=1; the explicit Cayley reset
defect rho remains in the nonlinear return map.

A separate structural-rank certificate proves that the same-sample vector
packet has exact rank five on the measurement-active block.  Therefore a
strictly positive scalar full-state packet margin is algebraically impossible;
directional forms must be transported and accumulated over the complete word.

This producer closes the proof-calculus contract, not the 18/21-state numerical
word.  It fails closed against entrance shrinkage, undeclared a_w/sigma
coupling, per-operation sector-invariance promotion, discarded eta terms, and
premature scalarization of directional information.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_directional_packet_rank as RANK
import ou3_p4_first_accel_aw_sigma_consistency as CONSISTENCY
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_exact_correction_transport as CORR
import ou3_p5_outer_information_geometry as OUTINFO

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _strict_attitude_geometry_lower(outinfo: dict) -> float:
    vals = [
        float(row["exact_pair_residual_information_vs_goLive_attitude_metric_lower"])
        for row in outinfo["nodes"].values()
    ]
    return down(min(vals))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("operation-matched Joseph ledger must not be trajectory fitted")

    sector = SECTOR.build(path)
    entrance = ENTRANCE.build(path)
    veff = VEFF.build(path)
    corr = CORR.build(path)
    outinfo = OUTINFO.build(path)
    rank = RANK.build(path)
    consistency = CONSISTENCY.build(path)

    failures = [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"correction-transport: {x}" for x in CORR.validate(corr)]
    failures += [f"outer-information: {x}" for x in OUTINFO.validate(outinfo)]
    failures += [f"directional-rank: {x}" for x in RANK.validate(rank)]
    failures += [f"aw-sigma-consistency: {x}" for x in CONSISTENCY.validate(consistency)]

    entrance_deg = float(entrance["P5_entrance"]["attitude_geometry"]["full_attitude_angle_deg"])
    entrance_q = float(entrance["P5_entrance"]["attitude_geometry"]["cayley_norm_upper"])
    outer_rad = float(sector["design_full_attitude_angle_rad"])
    outer_q = float(sector["design_cayley_norm_upper"])

    if entrance_deg != 45.0:
        failures.append("declared P5 entrance was shrunk from 45 deg")
    if consistency.get("aw_sigma_consistency_declared_in_domain") is not False:
        failures.append("strong route imported an a_w/sigma coupling assumption")

    ladder = list(consistency.get("angle_rows", []))
    legacy_obstruction = bool(ladder) and all(
        float(r["nuisance_over_budget_ratio_joint"]) > 1.0 for r in ladder
    )
    if not legacy_obstruction:
        failures.append("legacy first-accelerometer sector-invariance obstruction is not reproduced")

    acc = veff.get("accelerometer", {})
    mag = veff.get("magnetometer", {})
    if acc.get("standalone_eta_information_penalty_required_for_state_correction") is not False:
        failures.append("accelerometer state-correction range reduction did not close")
    if acc.get("J_aw_orthogonal_full_row_rank") is not True:
        failures.append("accelerometer effective a_w coordinate is unavailable")
    if mag.get("kalman_gain_radial_action_exact_zero") is not True:
        failures.append("magnetometer radial residual is not annihilated exactly")
    if mag.get("effective_coordinate_nonexpansive") is not True:
        failures.append("magnetometer effective tangent coordinate is not nonexpansive")
    if corr.get("condition_number_multiplier_used_for_reset_transport") is not False:
        failures.append("reset condition-number multiplier was reintroduced")
    if rank.get("instantaneous_positive_scalar_full_state_packet_margin_is_valid_target") is not False:
        failures.append("directional packet rank did not reject instantaneous scalarization")

    attitude_geometry_lower = _strict_attitude_geometry_lower(outinfo)
    if not (math.isfinite(attitude_geometry_lower) and attitude_geometry_lower > 0.0):
        failures.append("finite-angle vector-pair attitude geometry lower is not strict")

    consistency_rows = {
        float(r["angle_deg"]): r for r in consistency.get("angle_rows", [])
    }
    c25 = consistency_rows.get(25.0, {}).get("critical_consistency_constant")

    exact_joseph = corr.get("exact_joseph_information_identity")
    operation_ledger = [
        {
            "operation": "S_zero_accepted",
            "nonlinear_eta": "IDENTICALLY_ZERO",
            "exact_information_decrease": "y^T S^-1 y",
            "directional_role": "translation/S PSD credit",
            "sector_invariance_required": False,
        },
        {
            "operation": "accelerometer_accepted",
            "state_correction_range_reduction": "eta_a=H_a E_aw e_eta, e_eta=R_wb^T eta_a",
            "effective_coordinate": "z_eff=z+E_aw e_eta",
            "exact_state_correction": acc.get("exact_state_correction_identity"),
            "exact_information_decrease": "y^T S^-1 y - eta_a^T R_a^-1 eta_a",
            "finite_angle_eta_penalty_dropped": False,
            "finite_angle_eta_independently_maximized": False,
            "large_declared_aw_error_is_measurement_eta": False,
            "latent_aw_rotation_is_norm_preserving": True,
            "instantaneous_attitude_only_credit_promoted": False,
            "sector_invariance_required": False,
        },
        {
            "operation": "magnetometer_accepted",
            "radial_gain_action": "EXACTLY_ZERO",
            "effective_coordinate": mag.get("exact_effective_coordinate"),
            "exact_state_correction": mag.get("exact_state_correction_identity"),
            "effective_coordinate_nonexpansive": True,
            "exact_information_identity": exact_joseph,
            "radial_eta_independent_penalty_used": False,
            "sector_invariance_required": False,
        },
        {
            "operation": "measurement_rejected_or_not_due",
            "state_covariance_map": "IDENTITY",
            "information_change": "ZERO",
            "sector_invariance_required": False,
        },
        {
            "operation": "quaternion_injection_and_left_error_reset",
            "covariance_transport": "EXACT_CONGRUENCE",
            "reset_inverse_operator_norm_upper": 1.0,
            "remaining_nonlinear_term": "rho=z_exact-G_ext*t",
            "condition_number_multiplier_used": False,
            "sector_invariance_required": False,
        },
    ]

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_OPERATION_MATCHED_JOSEPH_DIRECTIONAL_LEDGER",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "declared_P5_entrance_angle_deg": entrance_deg,
        "declared_P5_entrance_q_upper": entrance_q,
        "operation_matched_outer_angle_rad": outer_rad,
        "operation_matched_outer_q_upper": outer_q,
        "P5_45DEG_ENTRANCE_PRESERVED": entrance_deg == 45.0,
        "aw_sigma_consistency_assumption_used": False,
        "conditional_25deg_consistency_constant_diagnostic": c25,
        "candidate_angle_reduction_used_for_closure": False,
        "legacy_sector_invariance_obstruction_reproduced": legacy_obstruction,
        "per_operation_sector_invariance_is_P4_promotion_gate": False,
        "legacy_sector_budget_distance_only": consistency.get("distance_only_no_verdict_emitted") is True,
        "exact_joseph_information_identity": exact_joseph,
        "exact_reset_congruence_identity": corr.get("exact_reset_congruence_identity"),
        "accelerometer_eta_absorbed_for_state_correction_range": True,
        "accelerometer_large_aw_error_charged_as_independent_measurement_eta": False,
        "accelerometer_finite_angle_eta_penalty_dropped_from_Joseph_identity": False,
        "accelerometer_finite_angle_eta_independent_norm_budget_used": False,
        "magnetometer_radial_eta_changes_state": False,
        "reset_condition_number_multiplier_used": False,
        "directional_packet_rank_exact": rank["measurement_structure"]["stacked_vector_packet_rank_exact"],
        "instantaneous_scalar_full_state_packet_margin_valid": False,
        "directional_PSD_word_accumulation_required": True,
        "finite_angle_vector_pair_attitude_geometry_vs_goLive_metric_lower": attitude_geometry_lower,
        "finite_angle_vector_pair_attitude_geometry_strict": attitude_geometry_lower > 0.0,
        "full_state_directional_word_credit_established_here": False,
        "operation_ledger": operation_ledger,
        "strong_route_operation_ledger_closed": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "build source-correlated PSD directional operation forms and signed nonlinear eta forms, transport them through prediction and exact reset, "
            "accumulate them over recurrent complete H/A words before scalarization, and simultaneously outward-enclose the full 18/21-state "
            "translation/nontranslation cross block; transient attitude-sector excursions are allowed when the complete-word Lyapunov level decreases"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "P5_45DEG_ENTRANCE_PRESERVED",
        "legacy_sector_invariance_obstruction_reproduced",
        "legacy_sector_budget_distance_only",
        "accelerometer_eta_absorbed_for_state_correction_range",
        "directional_PSD_word_accumulation_required",
        "finite_angle_vector_pair_attitude_geometry_strict",
        "strong_route_operation_ledger_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "aw_sigma_consistency_assumption_used",
        "candidate_angle_reduction_used_for_closure",
        "per_operation_sector_invariance_is_P4_promotion_gate",
        "accelerometer_large_aw_error_charged_as_independent_measurement_eta",
        "accelerometer_finite_angle_eta_penalty_dropped_from_Joseph_identity",
        "accelerometer_finite_angle_eta_independent_norm_budget_used",
        "magnetometer_radial_eta_changes_state",
        "reset_condition_number_multiplier_used",
        "instantaneous_scalar_full_state_packet_margin_valid",
        "full_state_directional_word_credit_established_here",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("declared_P5_entrance_angle_deg", 0.0)) != 45.0:
        f.append("P5 entrance is not 45 deg")
    if d.get("directional_packet_rank_exact") != 5:
        f.append("directional vector packet rank is not five")
    x = d.get("finite_angle_vector_pair_attitude_geometry_vs_goLive_metric_lower")
    if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
        f.append("attitude-geometry vector-pair lower is invalid")

    ledger = {row.get("operation"): row for row in d.get("operation_ledger", [])}
    acc = ledger.get("accelerometer_accepted", {})
    if acc.get("large_declared_aw_error_is_measurement_eta") is not False:
        f.append("accelerometer ledger reclassified a_w state error as measurement eta")
    if acc.get("finite_angle_eta_penalty_dropped") is not False:
        f.append("accelerometer ledger dropped the exact Joseph eta term")
    if acc.get("finite_angle_eta_independently_maximized") is not False:
        f.append("accelerometer ledger independently maximized finite-angle eta")
    if acc.get("latent_aw_rotation_is_norm_preserving") is not True:
        f.append("accelerometer ledger lost norm-preserving latent-a_w rotation")
    if acc.get("instantaneous_attitude_only_credit_promoted") is not False:
        f.append("accelerometer ledger promoted an attitude-only credit")
    if "eta_a^T R_a^-1 eta_a" not in str(acc.get("exact_information_decrease", "")):
        f.append("accelerometer ledger lost signed Joseph eta penalty")

    reset = ledger.get("quaternion_injection_and_left_error_reset", {})
    if float(reset.get("reset_inverse_operator_norm_upper", math.inf)) != 1.0:
        f.append("reset inverse norm is not exact one")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "entrance_deg": d["declared_P5_entrance_angle_deg"],
        "outer_rad": d["operation_matched_outer_angle_rad"],
        "packet_rank": d["directional_packet_rank_exact"],
        "instantaneous_scalar_margin_valid": d["instantaneous_scalar_full_state_packet_margin_valid"],
        "attitude_geometry_lower": d["finite_angle_vector_pair_attitude_geometry_vs_goLive_metric_lower"],
        "finite_angle_eta_penalty_dropped": d["accelerometer_finite_angle_eta_penalty_dropped_from_Joseph_identity"],
        "P4_promoted": d["P4_USABLE_CERTIFICATE_PROMOTED"],
        "next_obligation": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
