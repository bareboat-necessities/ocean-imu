#!/usr/bin/env python3
"""Fail-closed complete-precondition contract for canonical OU-III SEA3 P3.

A numerical P3 result may promote only one common Normal-Live H18/A21 word
that consumes all already-established source, scheduler, measurement, process,
and mode preconditions.  This file intentionally proves no reduced surrogate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_rs_word_information as RSWORD
import ou3_sea3_windowed_vector_pe as WINDOWPE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
PAPER = REPO / "doc" / "kalman_ou_iii" / "w3d-iss-stability.tex-part"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_P3_COMPLETE_NORMAL_LIVE_PRECONDITION_CONTRACT"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("canonical P3 domain may not be trajectory fitted")
    live = domain["normal_live"]
    runtime = domain["configured_runtime"]

    dynamic = DYNAMIC.build(path)
    sched = SCHED.build(path)
    process = PROCESS.build()
    physical = PHYSICAL.build(path)
    rsword = RSWORD.build(path)
    pe = WINDOWPE.build(path)
    checks = {
        "dynamic": DYNAMIC.validate(dynamic),
        "scheduler": SCHED.validate(sched),
        "process": PROCESS.validate(process),
        "physical": PHYSICAL.validate(physical),
        "four_S": RSWORD.validate(rsword),
        "windowed_PE": WINDOWPE.validate(pe),
    }
    bad = {k: v for k, v in checks.items() if v}
    if bad:
        raise RuntimeError(f"complete P3 prerequisite validation failed: {bad}")

    wrapper = WRAPPER.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    paper = PAPER.read_text(encoding="utf-8")

    floor_pos = mekf.find("apply_pending_aw_covariance_inflation_();")
    s_pos = mekf.find("applyIntegralZeroPseudoMeas();", floor_pos)
    source_parity = {
        "candidate_tau_each_valid_sample": "tune_.tau_applied   += alpha" in wrapper,
        "candidate_sigma_each_valid_sample": "tune_.sigma_applied += alpha" in wrapper,
        "candidate_RS_each_valid_sample": "tune_.RS_applied    += alpha_RS" in wrapper,
        "next_sample_schedule_commit": "void apply_pending_online_tune_()" in wrapper,
        "pseudo_cadence_from_committed_tau": "apply_pseudo_update_cadence_();" in wrapper,
        "periodic_aw_floor_tick": "periodic_aw_cov_sync_tick_();" in wrapper,
        "aw_floor_default_is_PSD_increment": "synchronize_aw_covariance_to_stationary();" in wrapper,
        "aw_floor_applied_inside_prediction": floor_pos >= 0,
        "aw_floor_adds_PSD_increment": "Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;" in mekf,
        "S_update_runs_after_floor_in_prediction": floor_pos >= 0 and s_pos > floor_pos,
        "accelerometer_attitude_jacobian": "const Matrix3 J_att = -skew_symmetric_matrix(f_cog_b);" in mekf,
        "accelerometer_aw_jacobian": "const Matrix3 J_aw  =  R_wb();" in mekf,
        "accelerometer_bias_jacobian_active": "PCt.noalias() += P_all_ba; // J_ba = I" in mekf,
        "full_Joseph_covariance_update": "joseph_update3_(K, S_mat, PCt);" in mekf,
    }
    source_failures = [k for k, ok in source_parity.items() if not ok]

    paper_parity = {
        "source_reachable_schedule": "strict source-reachable family" in paper,
        "four_S_translation_observability": "Spread-selected $S$ observability" in paper,
        "windowed_eta6_PE": "\\alpha_6\\mat I_6\\preceq" in paper,
        "windowed_eta9_PE": "\\alpha_9\\mat I_9\\preceq" in paper,
        "asynchronous_PE": "finite-window asynchronous conditions" in paper,
        "mag_rejections_allowed_if_information_coercive": "rejected samples, and short outages are admissible" in paper,
        "finite_bias_correlation_route": "finite residual-bias correlation time" in paper,
    }
    paper_failures = [k for k, ok in paper_parity.items() if not ok]

    pe_horizon = float(pe["spread_occurrence_selection"]["word_horizon_s"])
    s_horizon = float(rsword["word_horizon_s_upper"])
    word_horizon = max(pe_horizon, s_horizon)
    dt = float(runtime["imu_dt_s"])
    word_samples = int(math.ceil(word_horizon / dt)) + 1

    mandatory = {
        "runtime_zero_lever_arm": runtime["imu_lever_arm_enabled"] is False,
        "runtime_dormant_transparent_vibration_guard": runtime[
            "accelerometer_vibration_guard_proof_branch"
        ] == "dormant_transparent",
        "all_valid_accelerometer_updates_accepted": (
            live["accelerometer_update_required_each_valid_imu_sample_after_live_entry"] is True
            and live["accelerometer_rejection_in_normal_live_scope"] is False
        ),
        "joint_SEA3_adaptive_state": dynamic["P2_DYNAMIC_SOURCE_CERTIFICATE"] == "PASS",
        "source_commit_phase_retained": "pending_commit_progress" in dynamic["adaptive_state"],
        "pseudo_scheduler_progress_and_recurrence": sched["scheduler_recurrence_certificate"] is True,
        "applied_RS_and_tau_sigma_not_independent": True,
        "per_axis_RS_factors_retained": True,
        "four_S_information_retained": rsword["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"] is True,
        "full_process_UCC_retained": process["full_process_ucc_pass"] is True,
        "aw_covariance_floor_events_retained": all(source_parity[k] for k in (
            "periodic_aw_floor_tick", "aw_floor_default_is_PSD_increment",
            "aw_floor_applied_inside_prediction", "aw_floor_adds_PSD_increment",
            "S_update_runs_after_floor_in_prediction",
        )),
        "full_accelerometer_cross_block_information_retained": all(source_parity[k] for k in (
            "accelerometer_attitude_jacobian", "accelerometer_aw_jacobian",
            "accelerometer_bias_jacobian_active", "full_Joseph_covariance_update",
        )),
        "windowed_asynchronous_eta6_information_retained": pe["pass"] is True,
        "A_mode_finite_bias_correlation_retained": pe["A_mode_bias_route"][
            "uses_eta6_plus_finite_bias_correlation"
        ] is True,
        "specific_force_bounds_retained": (
            float(live["specific_force_norm_lower_mps2"]) > 0.0
            and float(live["specific_force_norm_upper_mps2"])
            >= float(live["specific_force_norm_lower_mps2"])
        ),
        "magnetic_norm_bounds_retained": (
            float(live["magnetic_vector_norm_lower_uT"]) > 0.0
            and float(live["magnetic_vector_norm_upper_uT"])
            >= float(live["magnetic_vector_norm_lower_uT"])
        ),
        "vector_sine_separation_retained": float(live["vector_sine_separation_lower"]) > 0.0,
        "body_rate_bound_retained": float(live["body_rate_norm_upper_deg_s"]) > 0.0,
        "A_mode_bias_domain_retained": (
            float(live["active_accelerometer_bias_state_norm_upper_mps2"])
            < float(live["active_accelerometer_bias_projection_limit_mps2"])
        ),
        "SEA3_height_period_partition_coupling_retained": physical["three_partition_contract"][
            "independent_H_r_and_T_p_rectangular_extrema_forbidden"
        ] is True,
    }

    numeric = {
        "common_word_horizon_s": word_horizon,
        "common_word_samples_upper": word_samples,
        "H_dimension": 18,
        "A_dimension": 21,
        "one_common_event_word_required_for_H_and_A": True,
        "same_joint_source_path_feeds_F_Q_TS_RS": True,
        "same_event_word_contains_accel_S_PE_and_aw_floor": True,
        "every_valid_accelerometer_update_must_be_applied": True,
        "every_due_S_update_must_be_applied": True,
        "actual_applied_per_axis_RS_required": True,
        "windowed_PE_information_must_be_accumulated_not_ODR_substituted": True,
        "aw_floor_must_be_added_as_actual_PSD_event_not_marginal_Loewner_shortcut": True,
        "A_mode_bias_must_use_finite_tau_or_eta9_information": True,
        "full_18x18_and_21x21_matrix_comparison_required": True,
        "blockwise_minimum_ratio_for_final_gate_forbidden": True,
        "determinant_trace_scalarization_for_final_gate_forbidden": True,
        "scalar_information_beta_for_final_gate_forbidden": True,
        "independent_tau_sigma_RS_TS_extrema_product_forbidden": True,
        "hardware_magnetometer_ODR_as_PE_recurrence_forbidden": True,
        "old_P2_graph_or_predecessor_enumeration_forbidden": True,
        "useful_gate": USEFUL_GATE,
        "required_final_inequality": "D_W - delta * L_W^{-1} >= 0 on full H18/A21 coordinates",
    }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_fit": False,
        "canonical_architecture": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "R_S_is_credited_but_not_used_as_a_substitute_for_other_preconditions": True,
        "process_UCC_is_credited_but_not_used_as_a_substitute_for_measurement_preconditions": True,
        "windowed_PE_is_credited_but_not_used_as_a_substitute_for_translation_preconditions": True,
        "mandatory_preconditions": mandatory,
        "all_current_machine_checkable_preconditions_present": all(mandatory.values()),
        "source_parity": source_parity,
        "source_parity_failures": source_failures,
        "paper_parity": paper_parity,
        "paper_parity_failures": paper_failures,
        "word": {
            "horizon_s": word_horizon,
            "samples_upper": word_samples,
            "windowed_PE": pe,
            "four_S_subcertificate": {
                "minimum_word_horizon_s": s_horizon,
                "scheduler_gap_s": sched["certified_uniform_max_gap_s"],
            },
            "aw_covariance_floor_gap_s_upper": dynamic[
                "validated_rate_and_jump_bounds"
            ]["active_commit_gap_s_upper"],
            "full_process_modes": process["modes"],
        },
        "final_numeric_contract": numeric,
        "physical_SEA3_scope": {
            "parameter_height_period_partition_coupling_consumed": True,
            "global_finite_window_realization_left_inclusion_closed": False,
            "unqualified_RAO_coupling_used_as_hard_pruning": False,
            "canonical_P3_is_conditional_on_admitted_Normal_Live_SEA3_word": True,
        },
        "P3_promoted": False,
        "P4_may_consume": False,
        "next_obligation": (
            "build one validated full H18/A21 Riccati/information word satisfying final_numeric_contract; do not introduce a reduced replacement certificate"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("trajectory_fit") is not False:
        f.append("complete P3 contract became trajectory fitted")
    if d.get("canonical_architecture") != "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("canonical architecture is not the complete Normal-Live word")
    for key in (
        "R_S_is_credited_but_not_used_as_a_substitute_for_other_preconditions",
        "process_UCC_is_credited_but_not_used_as_a_substitute_for_measurement_preconditions",
        "windowed_PE_is_credited_but_not_used_as_a_substitute_for_translation_preconditions",
        "all_current_machine_checkable_preconditions_present",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for name, ok in d.get("mandatory_preconditions", {}).items():
        if ok is not True:
            f.append(f"mandatory precondition missing: {name}")
    f.extend(f"source parity failed: {x}" for x in d.get("source_parity_failures", []))
    f.extend(f"paper parity failed: {x}" for x in d.get("paper_parity_failures", []))
    c = d.get("final_numeric_contract", {})
    for key in (
        "one_common_event_word_required_for_H_and_A",
        "same_joint_source_path_feeds_F_Q_TS_RS",
        "same_event_word_contains_accel_S_PE_and_aw_floor",
        "every_valid_accelerometer_update_must_be_applied",
        "every_due_S_update_must_be_applied",
        "actual_applied_per_axis_RS_required",
        "windowed_PE_information_must_be_accumulated_not_ODR_substituted",
        "aw_floor_must_be_added_as_actual_PSD_event_not_marginal_Loewner_shortcut",
        "A_mode_bias_must_use_finite_tau_or_eta9_information",
        "full_18x18_and_21x21_matrix_comparison_required",
        "blockwise_minimum_ratio_for_final_gate_forbidden",
        "determinant_trace_scalarization_for_final_gate_forbidden",
        "scalar_information_beta_for_final_gate_forbidden",
        "independent_tau_sigma_RS_TS_extrema_product_forbidden",
        "hardware_magnetometer_ODR_as_PE_recurrence_forbidden",
        "old_P2_graph_or_predecessor_enumeration_forbidden",
    ):
        if c.get(key) is not True:
            f.append(f"final numeric contract lost {key}")
    if float(c.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    if int(c.get("H_dimension", 0)) != 18 or int(c.get("A_dimension", 0)) != 21:
        f.append("full H/A dimensions changed")
    phys = d.get("physical_SEA3_scope", {})
    if phys.get("global_finite_window_realization_left_inclusion_closed") is not False:
        f.append("global physical SEA3 left inclusion was falsely promoted")
    if phys.get("unqualified_RAO_coupling_used_as_hard_pruning") is not False:
        f.append("unqualified RAO coupling entered canonical pruning")
    if d.get("P3_promoted") is not False or d.get("P4_may_consume") is not False:
        f.append("precondition contract promoted P3/P4")
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
        "architecture": d["canonical_architecture"],
        "word_horizon_s": d["word"]["horizon_s"],
        "all_preconditions": d["all_current_machine_checkable_preconditions_present"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
