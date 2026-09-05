#!/usr/bin/env python3
"""Canonical admitted-SEA3 Normal-Live source family for OU-III P3.

This module defines *which perturbations the full 3 s H18/A21 assembler is
allowed to propagate*.  It is deliberately not a Riccati certificate.

The primitive source variables are the admitted SEA3/Normal-Live physical and
measurement-only front-end variables.  The Kalman tuning coordinates

    tau_applied, sigma_aw, R_S_applied, T_S

are **derived state**, never independently selectable perturbations.  The same
front-end/tuner/candidate/commit state generates the transition F_k, process
Q_k, pseudo cadence and pseudo covariance used by the literal Riccati word.

The theorem remains conditional on an admitted finite SEA3 word.  Closing the
separate global physical-spectrum -> deterministic finite-window left inclusion
is not required to use these restrictions inside conditional P3, and this file
does not falsely claim that left inclusion.

The dynamic rectangular invariant retained elsewhere is a safety enclosure
only.  It may detect implementation/domain violations, but it may neither
generate canonical P3 perturbations nor reject the canonical architecture.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_wave_period_frontend as FRONTEND
import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_ADMITTED_NORMAL_LIVE_SOURCE_FAMILY"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("canonical SEA3 source family may not be replay fitted")

    physical = PHYSICAL.build(path)
    dynamic = DYNAMIC.build(path)
    scheduler = SCHED.build(path)
    frontend = FRONTEND.build(REPO)
    source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
    prereq = {
        "physical": PHYSICAL.validate(physical),
        "dynamic": DYNAMIC.validate(dynamic),
        "scheduler": SCHED.validate(scheduler),
        "frontend": FRONTEND.validate(frontend),
    }
    bad = {k: v for k, v in prereq.items() if v}
    if bad:
        raise RuntimeError(f"SEA3 source prerequisite failed: {bad}")

    live = domain["normal_live"]
    runtime = domain["configured_runtime"]
    inv = dynamic["dynamic_invariant"]
    rates = dynamic["validated_rate_and_jump_bounds"]
    phys = physical["three_partition_contract"]

    # These are predicates/coordinates of the *conditional* admitted word.
    # They are intentionally stronger than an arbitrary tuner box and exactly
    # match the assumptions the theorem says are in scope.
    sea3_predicates = {
        "partition_count_max_3": physical["sea_modes_max"] == 3,
        "partition_peak_steepness_enforced": (
            phys["active_partition_constraint"]
            == "2*pi*H_r/(g*T_p,r^2) <= S_p,max(T_p,r)"
        ),
        "partition_energy_coupling_enforced": (
            phys["total_energy_coupling"] == "H_s^2 = sum_r H_r^2"
        ),
        "independent_partition_H_T_extrema_forbidden": bool(
            phys["independent_H_r_and_T_p_rectangular_extrema_forbidden"]
        ),
        "independent_partition_height_maxima_forbidden": bool(
            phys["independent_three_partition_H_maxima_forbidden"]
        ),
        "total_Hs_cap_retained": float(phys["total_Hs_upper_m"]) == 8.5,
        "non_gravitational_CoG_acceleration_ball_retained": (
            float(live["non_gravitational_cog_acceleration_norm_upper_mps2"]) == 4.0
        ),
        "impact_slam_branch_excluded": (
            live["impact_slam_acceleration_in_normal_live_P1_P5_scope"] is False
        ),
        "specific_force_is_derived_not_independent": (
            live["specific_force_bounds_derived_from_gravity_and_non_gravitational_acceleration"]
            is True
        ),
        "all_valid_accelerometer_packets_accepted": (
            live["accelerometer_update_required_each_valid_imu_sample_after_live_entry"]
            is True
            and live["accelerometer_rejection_in_normal_live_scope"] is False
        ),
        "magnetic_norm_bounds_retained": (
            float(live["magnetic_vector_norm_lower_uT"]) > 0.0
            and float(live["magnetic_vector_norm_upper_uT"])
            >= float(live["magnetic_vector_norm_lower_uT"])
        ),
        "vector_sine_separation_retained": (
            float(live["vector_sine_separation_lower"]) > 0.0
        ),
        "asynchronous_vector_PE_recurrence_retained": (
            float(live["vector_pe_recurrence_window_s"]) > 0.0
        ),
        "body_rate_ball_retained": float(live["body_rate_norm_upper_deg_s"]) > 0.0,
        "zero_lever_arm_branch_retained": runtime["imu_lever_arm_enabled"] is False,
        "dormant_vibration_guard_branch_retained": (
            runtime["accelerometer_vibration_guard_proof_branch"]
            == "dormant_transparent"
        ),
        "no_hard_attitude_rewrite_inside_word": (
            live["hard_attitude_rewrite_inside_word"] is False
        ),
    }

    primitive_coordinates = [
        "SEA3 partition tuple {(H_r,T_p,r,gamma_r,beta_r,s_r)}_{r<=3} subject to partition energy/steepness predicates",
        "measurement-only vertical-leveling state",
        "WavePeriodEstimator high-pass/leaky-integrator/moment/log-period state",
        "adaptive sigma-band state",
        "acceleration-moment EMA state",
        "candidate tuner state (tau_candidate,sigma_candidate,R_S_candidate)",
        "pending/committed tuner state and commit clock",
        "pseudo scheduler elapsed phase",
        "physical CoG non-gravitational acceleration vector with ||a_ng||<=4 m/s^2",
        "body-rate vector inside declared Normal-Live ball",
        "accepted magnetic vector/reference state satisfying norm, separation and asynchronous PE predicates",
        "H/A same-mode bias state and declared projection/correlation constraints",
    ]

    derived_coordinates = {
        "tuning_frequency": (
            "shipping WavePeriodEstimator/front-end state; fixed prior only before startup-usable latch, which is excluded by Normal-Live entry"
        ),
        "tau_target": "shipping tau_coeff*0.5/f_tune with source clamps",
        "sigma_target": "shipping debiased acceleration-moment/sigma-band state with source clamps",
        "T_S_target": "clamp(pseudo_update_tau_ratio * tau_target, T_S_min, T_S_max)",
        "R_S_target": (
            "shipping selected law; deployed SpectralMSE uses the same tau_target, sigma_target and realized T_S_target"
        ),
        "candidate_tau_sigma": "same-sample common EMA driven by those SEA3-derived targets",
        "candidate_R_S": "shipping separate R_S EMA driven by the same SEA3-derived target state",
        "committed_tau_sigma_R_S": "next-sample staged commit of one candidate snapshot",
        "committed_T_S": "shipping monotone cadence map of committed tau",
        "OU_process_intensity": "q_c = 2*sigma_aw^2/tau using the same committed source state",
        "S_measurement_std": "[0.72*R_S_applied, 0.72*R_S_applied, R_S_applied] for deployed SpectralMSE law",
        "S_measurement_covariance": "diag(S_measurement_std^2)",
    }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_architecture": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "source_generated_not_trajectory_fit": True,
        "theorem_is_conditional_on_admitted_SEA3_word": True,
        "global_physical_SEA3_left_inclusion_closed_here": False,
        "global_physical_SEA3_left_inclusion_required_before_conditional_P3": False,
        "SEA3_predicates": sea3_predicates,
        "all_declared_SEA3_predicates_bound": all(sea3_predicates.values()),
        "primitive_source_coordinates": primitive_coordinates,
        "derived_filter_coordinates": derived_coordinates,
        "tau_sigma_RS_TS_are_primitive_independent_perturbations": False,
        "rectangular_dynamic_invariant_role": "SAFETY_ENVELOPE_ONLY",
        "rectangular_dynamic_invariant_may_generate_canonical_P3_words": False,
        "rectangular_dynamic_invariant_failure_may_reject_canonical_architecture": False,
        "target_jump_may_span_full_box_in_canonical_P3": False,
        "front_end_recurrences_must_generate_target_motion": True,
        "same_source_state_feeds": [
            "F_k", "Q_k", "pseudo_scheduler_period_T_S_k",
            "every_due_S_measurement_R_k", "accelerometer_H_k", "magnetic_H_k",
        ],
        "R_S_regularizer": {
            "primary_translation_regularizer": True,
            "actual_applied_R_S_required_at_every_due_S_update": True,
            "per_axis_std_factors": [0.72, 0.72, 1.0],
            "SpectralMSE_has_no_extra_information_rate_rescale": True,
            "pseudo_scheduler_recurrence_certificate_consumed": bool(
                scheduler["scheduler_recurrence_certificate"]
            ),
            "uniform_pseudo_gap_s": scheduler["certified_uniform_max_gap_s"],
            "selected_four_S_updates_are_only_rank_witness_subset": True,
            "all_additional_due_S_updates_must_remain_in_literal_word": True,
            "full_P_column_S_cross_covariance_action_must_be_retained": True,
        },
        "frontend": {
            "source_parity_pass": all(frontend["source_parity"].values()),
            "startup_usable_period_required_before_Normal_Live": True,
            "front_end_state_may_not_be_frozen_to_replay_value": True,
            "screened_tuning_frequency_is_safety_output_range_not_physical_Tp_domain": True,
        },
        "dynamic_safety_envelope": {
            "invariant": inv,
            "rate_and_jump_bounds": rates,
            "may_be_used_to_check_derived_state_stays_safe": True,
            "may_be_used_as_independent_cartesian_source_family": False,
        },
        "physical_parameter_subcertificate": {
            "finite_window_realization_enclosed": physical["finite_window_realization_enclosed"],
            "left_language_inclusion_closed": physical["left_language_inclusion_closed"],
            "conditional_P3_may_still_use_declared_admitted_SEA3_word": True,
        },
        "old_P2_800_state_graph_consumed": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "P3_promoted": False,
        "next_obligation": (
            "replace the point-source driver of the literal 3 s H18/A21 executor with a validated cover of these admitted SEA3 primitive coordinates, propagate the shipping front-end/tuner state, and apply actual R_S at every scheduler firing"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_architecture") != "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("canonical architecture changed")
    for key in (
        "source_generated_not_trajectory_fit",
        "theorem_is_conditional_on_admitted_SEA3_word",
        "all_declared_SEA3_predicates_bound",
        "front_end_recurrences_must_generate_target_motion",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "global_physical_SEA3_left_inclusion_closed_here",
        "global_physical_SEA3_left_inclusion_required_before_conditional_P3",
        "tau_sigma_RS_TS_are_primitive_independent_perturbations",
        "rectangular_dynamic_invariant_may_generate_canonical_P3_words",
        "rectangular_dynamic_invariant_failure_may_reject_canonical_architecture",
        "target_jump_may_span_full_box_in_canonical_P3",
        "old_P2_800_state_graph_consumed",
        "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("rectangular_dynamic_invariant_role") != "SAFETY_ENVELOPE_ONLY":
        f.append("dynamic rectangle regained canonical source role")
    rs = d.get("R_S_regularizer", {})
    for key in (
        "primary_translation_regularizer",
        "actual_applied_R_S_required_at_every_due_S_update",
        "SpectralMSE_has_no_extra_information_rate_rescale",
        "pseudo_scheduler_recurrence_certificate_consumed",
        "selected_four_S_updates_are_only_rank_witness_subset",
        "all_additional_due_S_updates_must_remain_in_literal_word",
        "full_P_column_S_cross_covariance_action_must_be_retained",
    ):
        if rs.get(key) is not True:
            f.append(f"R_S regularizer contract lost {key}")
    if rs.get("per_axis_std_factors") != [0.72, 0.72, 1.0]:
        f.append("deployed R_S axis factors changed")
    if not d.get("primitive_source_coordinates"):
        f.append("SEA3 primitive source coordinates missing")
    if not d.get("derived_filter_coordinates"):
        f.append("derived tuner/source coordinates missing")
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
        "qualification": d["qualification"],
        "SEA3_predicates": d["SEA3_predicates"],
        "rectangle_role": d["rectangular_dynamic_invariant_role"],
        "R_S_regularizer": d["R_S_regularizer"],
        "P3_promoted": d["P3_promoted"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
