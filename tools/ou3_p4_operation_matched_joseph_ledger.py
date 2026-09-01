#!/usr/bin/env python3
"""Operation-matched Joseph/directional ledger for the usable OU-III P4 route.

The first-accelerometer sector-invariance diagnostic is deliberately *not* the
acceptance test used here.  That diagnostic asks whether every accepted update,
viewed in isolation and with every non-attitude contribution charged as an
adverse correction, remains inside the 0.80 rad outer attitude sector.  The
shipping MEKF does not have that structure and PR #448 measured that the test
cannot close by shrinking the candidate angle alone.

The stronger route keeps the actual measurement/covariance algebra together.
For an accepted update with signed residual y = H z + eta and Joseph posterior
P+, the exact identity is

    z^T P^-1 z - (z-Ky)^T (P+)^-1 (z-Ky)
      = y^T S^-1 y - eta^T R^-1 eta.

For the configured accelerometer, J_aw=R_wb is orthogonal and full row rank.
The finite-angle residual is therefore represented exactly in the a_w state
coordinate: y_a = H_a(z + E_aw e_eta).  The large declared 0.3 g a_w error is a
state coordinate, not an independent measurement-space eta penalty.  The
latent-a_w rotation is norm preserving; only the source-correlated effective
state map remains for the complete-word numerical enclosure.

For the configured magnetometer the radial finite-angle residual is annihilated
exactly by the gain and the useful residual is an effective tangent coordinate.
For S=0, eta is identically zero.  Rejected/not-due measurements are exact
identity branches.  Quaternion reset is an exact covariance congruence with
||G^-1||_2=1; the only nonlinear reset term is the explicit Cayley defect rho.

This producer does not claim the complete 18/21-state P4 word.  It certifies the
strong-route ledger and makes CI reject three shortcuts:

* shrinking the declared 45 deg P5 entrance;
* adding an undeclared ||delta a_w|| <= c sigma_applied assumption; or
* reviving per-operation sector invariance as a P4 promotion gate.

The remaining numerical obligation is intentionally concrete: propagate the
source-correlated effective-state map, Joseph decrease, exact reset defect,
prediction, and full translation/nontranslation cross block over the same
complete H/A source paths.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_first_accel_aw_sigma_consistency as CONSISTENCY
import ou3_p4_first_accel_sector_budget as SECTOR_BUDGET
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_exact_correction_transport as CORR
import ou3_p5_outer_information_geometry as OUTINFO

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _strict_packet_lower(outinfo: dict) -> float:
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
    old_budget = SECTOR_BUDGET.build(path)
    consistency = CONSISTENCY.build(path)

    failures = [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"correction-transport: {x}" for x in CORR.validate(corr)]
    failures += [f"outer-information: {x}" for x in OUTINFO.validate(outinfo)]
    failures += [f"sector-budget: {x}" for x in SECTOR_BUDGET.validate(old_budget)]
    failures += [f"aw-sigma-consistency: {x}" for x in CONSISTENCY.validate(consistency)]

    entrance_deg = float(entrance["P5_entrance"]["attitude_geometry"]["angle_deg"])
    entrance_q = float(entrance["P5_entrance"]["attitude_geometry"]["cayley_norm_upper"])
    outer_rad = float(sector["design_full_attitude_angle_rad"])
    outer_q = float(sector["design_cayley_norm_upper"])

    if entrance_deg != 45.0:
        failures.append("declared P5 entrance was shrunk from 45 deg")
    if consistency.get("aw_sigma_consistency_declared_in_domain") is not False:
        failures.append("strong route imported an a_w/sigma coupling assumption")
    if old_budget.get("shrinking_the_candidate_angle_alone_can_close_the_budget") is not False:
        failures.append("legacy sector-invariance diagnostic unexpectedly became closable by angle alone")

    ladder = list(old_budget.get("ladder_rows", []))
    if not ladder or not all(float(r["nuisance_over_budget_ratio"]) > 1.0 for r in ladder):
        failures.append("legacy first-accelerometer obstruction is not reproduced")

    acc = veff.get("accelerometer", {})
    mag = veff.get("magnetometer", {})
    if acc.get("standalone_eta_information_penalty_required_for_state_correction") is not False:
        failures.append("accelerometer finite-angle residual was reintroduced as standalone eta")
    if acc.get("J_aw_orthogonal_full_row_rank") is not True:
        failures.append("accelerometer effective a_w coordinate is unavailable")
    if mag.get("kalman_gain_radial_action_exact_zero") is not True:
        failures.append("magnetometer radial residual is not annihilated exactly")
    if mag.get("effective_coordinate_nonexpansive") is not True:
        failures.append("magnetometer effective tangent coordinate is not nonexpansive")
    if corr.get("condition_number_multiplier_used_for_reset_transport") is not False:
        failures.append("reset condition-number multiplier was reintroduced")

    packet_lower = _strict_packet_lower(outinfo)
    if not (math.isfinite(packet_lower) and packet_lower > 0.0):
        failures.append("finite-angle vector-pair directional information lower is not strict")

    consistency_rows = {
        float(r["angle_deg"]): r for r in consistency.get("angle_rows", [])
    }
    c25 = consistency_rows.get(25.0, {}).get("critical_consistency_constant")

    operation_ledger = [
        {
            "operation": "S_zero_accepted",
            "effective_measurement_eta": "IDENTICALLY_ZERO",
            "joseph_information_change": "- y^T S^-1 y",
            "directional_role": "translation/S information decrease",
            "sector_invariance_required": False,
        },
        {
            "operation": "accelerometer_accepted",
            "effective_measurement_eta": "ABSORBED_EXACTLY_INTO_AW_STATE_COORDINATE",
            "effective_coordinate": "z_eff=z+E_aw e_eta, e_eta=R_wb^T eta_a",
            "exact_state_correction": acc.get("exact_state_correction_identity"),
            "joseph_information_change_in_effective_coordinate": "- y^T S^-1 y",
            "large_declared_aw_error_is_measurement_eta": False,
            "latent_aw_rotation_is_norm_preserving": True,
            "sector_invariance_required": False,
        },
        {
            "operation": "magnetometer_accepted",
            "radial_gain_action": "EXACTLY_ZERO",
            "effective_coordinate": mag.get("exact_effective_coordinate"),
            "exact_state_correction": mag.get("exact_state_correction_identity"),
            "effective_coordinate_nonexpansive": True,
            "joseph_information_change_in_effective_coordinate": "- y_tangent^T S^-1 y_tangent",
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
        "legacy_sector_invariance_obstruction_reproduced": bool(ladder) and all(
            float(r["nuisance_over_budget_ratio"]) > 1.0 for r in ladder
        ),
        "per_operation_sector_invariance_is_P4_promotion_gate": False,
        "legacy_sector_budget_distance_only": old_budget.get("distance_only_no_verdict_emitted") is True,
        "exact_joseph_information_identity": corr.get("exact_joseph_information_identity"),
        "exact_reset_congruence_identity": corr.get("exact_reset_congruence_identity"),
        "accelerometer_eta_absorbed_as_effective_aw_state": True,
        "accelerometer_large_aw_error_charged_as_independent_measurement_eta": False,
        "magnetometer_radial_eta_changes_state": False,
        "reset_condition_number_multiplier_used": False,
        "finite_angle_vector_pair_directional_information_vs_goLive_metric_lower": packet_lower,
        "finite_angle_vector_pair_directional_information_strict": packet_lower > 0.0,
        "operation_ledger": operation_ledger,
        "strong_route_first_vector_packet_ledger_closed": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "outward-propagate the source-correlated effective-state map z->z_eff, each operation's Joseph information decrease, "
            "the exact Cayley/reset defect rho, prediction, and the full 18/21-state translation/nontranslation cross block over the same recurrent H/A paths; "
            "allow transient attitude-sector excursions when the complete source-word Lyapunov level decreases"
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
        "accelerometer_eta_absorbed_as_effective_aw_state",
        "finite_angle_vector_pair_directional_information_strict",
        "strong_route_first_vector_packet_ledger_closed",
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
        "magnetometer_radial_eta_changes_state",
        "reset_condition_number_multiplier_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("declared_P5_entrance_angle_deg", 0.0)) != 45.0:
        f.append("P5 entrance is not 45 deg")
    x = d.get("finite_angle_vector_pair_directional_information_vs_goLive_metric_lower")
    if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
        f.append("directional vector-pair information lower is invalid")
    ledger = {row.get("operation"): row for row in d.get("operation_ledger", [])}
    acc = ledger.get("accelerometer_accepted", {})
    if acc.get("large_declared_aw_error_is_measurement_eta") is not False:
        f.append("accelerometer ledger reclassified a_w as measurement eta")
    if acc.get("latent_aw_rotation_is_norm_preserving") is not True:
        f.append("accelerometer ledger lost norm-preserving latent-a_w rotation")
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
        "vector_pair_directional_info_lower": d["finite_angle_vector_pair_directional_information_vs_goLive_metric_lower"],
        "legacy_invariance_obstruction_reproduced": d["legacy_sector_invariance_obstruction_reproduced"],
        "sector_invariance_is_gate": d["per_operation_sector_invariance_is_P4_promotion_gate"],
        "P4_promoted": d["P4_USABLE_CERTIFICATE_PROMOTED"],
        "next_obligation": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
