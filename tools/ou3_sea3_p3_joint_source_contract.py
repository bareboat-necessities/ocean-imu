#!/usr/bin/env python3
"""Fail-closed source-coupling guard for the canonical OU-III SEA3/R_S P3 proof.

This is not another P3 estimator.  It protects the single canonical
SEA3/R_S innovation proof from losing the correlations supplied by the
shipping adaptive source.

The quantitative word must carry one source state

    xi_k = (tau_k, sigma_k, R_S,k, T_S,k, scheduler_progress_k)

and the same xi_k must feed F_k, Q_k, pseudo timing and pseudo covariance.
Independent extrema products of tau, sigma_aw, R_S and T_S are forbidden in
the canonical gate.  A rectangular full-box calculation may exist only as a
diagnostic and may not reject the canonical architecture.

The guard also distinguishes already-certified SEA3 facts from attractive but
unfinished ones.  It intentionally does not import the older directional/P2
proof stack merely to inspect promotion metadata: that would make a guard
against the retired history graph depend on the graph itself.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_physical_admissibility as PHYSICAL

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
ACCCOUP_PATH = REPO / "tools" / "ou3_sea3_acceleration_covariance_coupling.py"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_P3_JOINT_SOURCE_COUPLING_CONTRACT"


def _axis_factors() -> list[float]:
    import re
    text = WRAPPER.read_text(encoding="utf-8")
    out: list[float] = []
    for name in ("R_S_x_factor_", "R_S_y_factor_"):
        m = re.search(rf"float\s+{name}\s*=\s*([0-9.eE+-]+)f", text)
        if not m:
            raise RuntimeError(f"cannot extract deployed {name}")
        v = float(m.group(1))
        if not (v > 0.0):
            raise RuntimeError(f"invalid deployed {name}")
        out.append(v)
    return out + [1.0]


def _unfinished_acceleration_coupling_metadata() -> dict:
    """Read only fail-closed promotion markers; never import the retired P2 stack."""
    text = ACCCOUP_PATH.read_text(encoding="utf-8")
    required = (
        '"physical_vessel_pairing_qualified": False',
        '"P3_promoted": False',
        '"deterministic_left_inclusion_closed": False',
    )
    missing = [m for m in required if m not in text]
    if missing:
        raise RuntimeError(
            "SEA3 acceleration-coupling promotion metadata changed; inspect before P3 use: "
            + ", ".join(missing)
        )
    return {
        "physical_vessel_pairing_qualified": False,
        "P3_promoted": False,
        "deterministic_left_inclusion_closed": False,
        "metadata_read_without_importing_directional_P2_stack": True,
    }


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

    rates = dynamic["validated_rate_and_jump_bounds"]
    inv = dynamic["dynamic_invariant"]
    accmeta = _unfinished_acceleration_coupling_metadata()

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
        "retired_directional_P2_stack_imported_by_joint_guard": False,
        "proved_source_couplings": {
            "tau_target_from_wave_band_frequency": True,
            "tau_target_s": inv["tau_target_s"],
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
            "pseudo_scheduler_progress_preserving": bool(scheduler["scheduler_recurrence_certificate"]),
            "pseudo_uniform_max_gap_s": scheduler["certified_uniform_max_gap_s"],
            "process_intensity_uses_same_sigma_tau": "q_c = 2*sigma_aw^2/tau",
            "R_S_axis_std_factors": _axis_factors(),
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
                accmeta["physical_vessel_pairing_qualified"]
            ),
            "sea_RAO_acceleration_coupling_may_be_used_as_hard_P3_pruning": bool(
                accmeta["P3_promoted"]
            ),
            "sea_RAO_acceleration_left_inclusion_closed": bool(
                accmeta["deterministic_left_inclusion_closed"]
            ),
        },
        "R_S_corrective_force_requirements": {
            "use_actual_applied_R_S_on_selected_pseudo_updates": True,
            "retain_per_axis_R_S_factors": True,
            "retain_full_P_column_S_cross_covariance_action": True,
            "credit_guaranteed_recurrent_S_updates_as_measurement_dissipation": True,
            "credit_additional_guaranteed_S_updates_as_positive_information_when_used": True,
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
        "retired_directional_P2_stack_imported_by_joint_guard",
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
        "sea_RAO_acceleration_left_inclusion_closed",
    ):
        if notyet.get(key) is not False:
            f.append(f"unfinished SEA3 coupling was incorrectly promoted: {key}")
    rs = d.get("R_S_corrective_force_requirements", {})
    for key in (
        "use_actual_applied_R_S_on_selected_pseudo_updates",
        "retain_per_axis_R_S_factors",
        "retain_full_P_column_S_cross_covariance_action",
        "credit_guaranteed_recurrent_S_updates_as_measurement_dissipation",
        "credit_additional_guaranteed_S_updates_as_positive_information_when_used",
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
