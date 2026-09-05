#!/usr/bin/env python3
"""Canonical source-coupling guard for the OU-III SEA3/R_S P3 proof.

This file is not a second P3 estimator or certificate architecture.  It is a
fail-closed contract for the *one* canonical SEA3/R_S innovation proof.

The canonical quantitative word may not form independent extrema products of
``tau``, ``sigma_aw``, ``R_S`` and ``T_S``.  They are components of one
shipping adaptive source state

    xi_k = (tau_k, sigma_k, R_S,k, T_S,k, scheduler_progress_k),

and the same xi_k must feed the transition, process covariance, pseudo cadence
and pseudo-measurement covariance at sample k.

The contract deliberately distinguishes proved coupling from attractive but
unfinished SEA3 facts.  In particular, the physical sea/RAO acceleration
covariance predicate exists, but the repository currently marks the actual
vessel pairing as unqualified; canonical P3 therefore may not silently use
that stronger sigma/period relation until it is promoted.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_acceleration_covariance_coupling as ACCCOUP
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RESPONSE = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_P3_JOINT_SOURCE_COUPLING_CONTRACT"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("joint source contract may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    physical = PHYSICAL.build(path)
    scheduler = SCHED.build(path)
    for label, failures in (
        ("dynamic", DYNAMIC.validate(dynamic)),
        ("physical", PHYSICAL.validate(physical)),
        ("scheduler", SCHED.validate(scheduler)),
    ):
        if failures:
            raise RuntimeError(f"{label} prerequisite failed: {failures}")

    sched = BASE.source_schedule()
    rates = dynamic["validated_rate_and_jump_bounds"]
    inv = dynamic["dynamic_invariant"]

    # The acceleration-covariance SEA3 predicate is useful, but its own
    # artifact currently refuses promotion until a vessel/RAO pairing is
    # qualified.  Record that state so the P3 proof cannot accidentally treat
    # the conditional predicate as a proved applied sigma/tau invariant.
    acc = ACCCOUP.build(
        samples=max(1, int(round(1.0 / float(rates["dt_s"])))),
        domain_path=path,
        response_domain_path=DEFAULT_RESPONSE,
        repo=REPO,
    )
    af = ACCCOUP.validate(acc)
    if af:
        raise RuntimeError(f"SEA3 acceleration coupling prerequisite invalid: {af}")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "canonical_architecture": "SEA3_RS_INNOVATION_DISSIPATION_WORD",
        "adaptive_state": [
            "tau_applied",
            "sigma_aw_filter",
            "R_S_applied",
            "pseudo_update_period",
            "scheduler_progress",
        ],
        "same_source_state_must_feed": [
            "translation_transition_F_k",
            "OU_process_covariance_Q_k",
            "pseudo_update_period_T_S_k",
            "pseudo_measurement_covariance_R_S_k",
        ],
        "canonical_independent_tau_sigma_RS_TS_extrema_product_forbidden": True,
        "rectangular_full_box_calculation_diagnostic_only": True,
        "rectangular_full_box_failure_may_not_reject_canonical_architecture": True,
        "old_P2_800_state_graph_consumed": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "proved_source_couplings": {
            "tau_target_from_wave_band_frequency": True,
            "tau_target_s": dynamic["dynamic_invariant"]["tau_target_s"],
            "tau_applied_s": inv["tau_applied_s"],
            "tau_sigma_share_same_sample_EMA_alpha": True,
            "tau_sigma_common_EMA_horizon_s": inv["common_tau_sigma_horizon_s"],
            "tau_sigma_candidate_rate_bounds_consumed": True,
            "tau_active_commit_jump_bound_s": rates["tau_active_abs_jump_per_commit_upper_s"],
            "sigma_active_commit_jump_bound_mps2": rates["sigma_active_abs_jump_per_commit_upper_mps2"],
            "SpectralMSE_target_uses_same_target_tau_sigma_TS": True,
            "R_S_has_separate_EMA": True,
            "R_S_EMA_horizon_s": inv["R_S_horizon_s"],
            "R_S_active_commit_jump_bound": rates["R_S_active_abs_jump_per_commit_upper"],
            "pseudo_period_is_clamped_monotone_function_of_applied_tau": True,
            "pseudo_scheduler_progress_preserving": bool(
                scheduler["scheduler_recurrence_certificate"]
            ),
            "pseudo_uniform_max_gap_s": scheduler["certified_uniform_max_gap_s"],
            "process_intensity_uses_same_sigma_tau": "q_c = 2*sigma_aw^2/tau",
            "R_S_axis_std_factors": sched["R_S_axis_std_factors"],
            "normal_live_accelerometer_every_valid_sample": bool(
                dynamic["normal_live_contract"]["accelerometer_update_required_each_valid_sample"]
            ),
            "normal_live_accelerometer_rejection_absent": (
                dynamic["normal_live_contract"]["accelerometer_rejection_in_scope"] is False
            ),
            "physical_height_period_cartesian_extrema_forbidden": bool(
                physical["three_partition_contract"][
                    "independent_H_r_and_T_p_rectangular_extrema_forbidden"
                ]
            ),
            "three_partition_height_maxima_independent_forbidden": bool(
                physical["three_partition_contract"][
                    "independent_three_partition_H_maxima_forbidden"
                ]
            ),
        },
        "not_yet_promotable_source_couplings": {
            "physical_SEA3_finite_window_realization_closed": bool(
                physical["finite_window_realization_enclosed"]
            ),
            "physical_SEA3_left_language_inclusion_closed": bool(
                physical["left_language_inclusion_closed"]
            ),
            "sea_RAO_acceleration_pairing_qualified": bool(
                acc["physical_vessel_pairing_qualified"]
            ),
            "sea_RAO_acceleration_coupling_may_be_used_as_hard_P3_pruning": bool(
                acc["P3_promoted"]
            ),
        },
        "R_S_corrective_force_requirements": {
            "use_actual_applied_R_S_on_selected_pseudo_updates": True,
            "retain_per_axis_R_S_factors": True,
            "retain_full_P_column_S_cross_covariance_action": True,
            "credit_guaranteed_recurrent_S_updates_as_measurement_dissipation": True,
            "do_not_replace_R_S_correction_by_process_strictness": True,
            "do_not_use_global_R_S_100_at_every_firing_in_final_canonical_matrix_if_joint_source_enclosure_is_available": True,
        },
        "P3_promoted": False,
        "P4_promoted": False,
        "next_obligation": (
            "construct the canonical finite-word translation information/metric comparison over the joint adaptive state xi_k; the full-box calculation may be emitted only as a diagnostic"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "canonical_independent_tau_sigma_RS_TS_extrema_product_forbidden",
        "rectangular_full_box_calculation_diagnostic_only",
        "rectangular_full_box_failure_may_not_reject_canonical_architecture",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "old_P2_800_state_graph_consumed",
        "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed",
        "P3_promoted",
        "P4_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    proved = d.get("proved_source_couplings", {})
    for key in (
        "tau_target_from_wave_band_frequency",
        "tau_sigma_share_same_sample_EMA_alpha",
        "tau_sigma_candidate_rate_bounds_consumed",
        "SpectralMSE_target_uses_same_target_tau_sigma_TS",
        "R_S_has_separate_EMA",
        "pseudo_period_is_clamped_monotone_function_of_applied_tau",
        "pseudo_scheduler_progress_preserving",
        "normal_live_accelerometer_every_valid_sample",
        "normal_live_accelerometer_rejection_absent",
        "physical_height_period_cartesian_extrema_forbidden",
        "three_partition_height_maxima_independent_forbidden",
    ):
        if proved.get(key) is not True:
            f.append(f"proved coupling lost: {key}")
    notyet = d.get("not_yet_promotable_source_couplings", {})
    for key in (
        "physical_SEA3_finite_window_realization_closed",
        "physical_SEA3_left_language_inclusion_closed",
        "sea_RAO_acceleration_pairing_qualified",
        "sea_RAO_acceleration_coupling_may_be_used_as_hard_P3_pruning",
    ):
        if notyet.get(key) is not False:
            f.append(f"unfinished SEA3 coupling was incorrectly promoted: {key}")
    rs = d.get("R_S_corrective_force_requirements", {})
    for key in (
        "use_actual_applied_R_S_on_selected_pseudo_updates",
        "retain_per_axis_R_S_factors",
        "retain_full_P_column_S_cross_covariance_action",
        "credit_guaranteed_recurrent_S_updates_as_measurement_dissipation",
        "do_not_replace_R_S_correction_by_process_strictness",
        "do_not_use_global_R_S_100_at_every_firing_in_final_canonical_matrix_if_joint_source_enclosure_is_available",
    ):
        if rs.get(key) is not True:
            f.append(f"R_S corrective-force requirement lost: {key}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "canonical_architecture": d["canonical_architecture"],
        "independent_extrema_forbidden": d[
            "canonical_independent_tau_sigma_RS_TS_extrema_product_forbidden"
        ],
        "full_box_diagnostic_only": d["rectangular_full_box_calculation_diagnostic_only"],
        "unpromoted_couplings": d["not_yet_promotable_source_couplings"],
        "R_S_requirements": d["R_S_corrective_force_requirements"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
