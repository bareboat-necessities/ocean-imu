#!/usr/bin/env python3
"""Canonical OU-III P3 architecture: SEA3 + recurrent R_S innovation dissipation.

The dominant stabilizing mechanism is the shipping measurement correction, not
one-sample process-noise injection.  For every exact linear Kalman measurement
update with prior error e-, innovation r=H e-, innovation covariance
S=H P- H' + R, posterior covariance P+, and homogeneous posterior error e+,

    V- - V+ = r' S^-1 r,          V=e' P^-1 e.

Prediction is non-expansive because P-=F P+ F'+Q with Q>=0.  Hence over a
recurrent word the accumulated sequential-innovation energy is a strictness
candidate.  The innovations are the block-Cholesky whitening of the complete
batch record, so the correction information is

    D_W = O_W' Sigma_Y^-1 O_W.

Translation now uses the strongest cleanly decoupled mechanism available in the
shipping filter: four separated S=0 pseudo updates.  The dedicated
``ou3_sea3_rs_word_information`` certificate proves that these four S rows are
full rank on [v,p,S,a_w] even when tau varies arbitrarily inside the word.  The
proof uses c'''(t)=exp(-int 1/tau)>0 for the a_w->S response and a third divided
difference, not a frozen source word.  Thus translation strictness can be built
from R_S correction alone; accelerometer information is no longer required to
repair a fourth translation direction.

This is important for the full H/A proof because accelerometer residuals couple
attitude, a_w and active b_a.  Keeping the translation block on S=0 observations
leaves accelerometer/magnetometer PE available for the genuinely coupled
attitude/bias block instead of pretending the accelerometer observes a_w in
isolation.

The deployed SpectralMSE *target* couples tau, sigma and realized T_S, but the
applied R_S state has its own EMA.  Canonical P3 therefore never substitutes the
instantaneous target law for active R_S.  The four-S certificate currently uses
the safe applied R_S ceiling and the exact source-independent scheduler
recurrence.  A future SEA3 lag/reachability theorem may tighten that bound, but
P3 does not depend on such an unproved tightening.

Process UCC is retained only to provide a finite-memory covariance lower L_W in
the same observation/divided-difference coordinates.  The final useful gate is
one full H18/A21 matrix comparison at delta>=1e-18.  No one-step strictness,
per-sample SPD lower, commit-word propagation, determinant/trace eigenvalue
scalarization, scalar information beta, source-history graph, or predecessor
path is canonical.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_riccati_tube_factored as TUBE
import ou3_sea3_rs_word_information as RSWORD
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_RS_INNOVATION_DISSIPATION_P3_ARCHITECTURE"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P3 architecture may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    physical = PHYSICAL.build(path)
    vector = VECTOR.build()
    process = PROCESS.build()
    trans = TRANS.build(TRANS.DEFAULT_HEADER)
    rsword = RSWORD.build(path)
    prereq = (
        [f"dynamic: {x}" for x in DYNAMIC.validate(dynamic)]
        + [f"physical: {x}" for x in PHYSICAL.validate(physical)]
        + [f"vector: {x}" for x in VECTOR.validate(vector)]
        + [f"process: {x}" for x in PROCESS.validate(process)]
        + [f"translation: {x}" for x in TRANS.validate(trans)]
        + [f"R_S word: {x}" for x in RSWORD.validate(rsword)]
    )
    if prereq:
        raise RuntimeError(f"P3 prerequisites failed: {prereq}")

    if tube_path is None:
        tube = TUBE.build(path)
    else:
        tube = json.loads(Path(tube_path).read_text(encoding="utf-8"))
    tf = TUBE.validate(tube)
    if tf:
        raise RuntimeError(f"endpoint covariance prerequisite failed: {tf}")

    inv = dynamic["dynamic_invariant"]
    live = domain["normal_live"]
    timing = tube["covariance_memory"]

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "SEA3_RS_INNOVATION_DISSIPATION_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_consumed": True,
        "SEA3_physical_parameter_coupling_consumed": True,
        "SEA3_full_finite_window_realization_assumed_conditionally_not_falsely_promoted": True,
        "R_S_is_primary_translation_correction_mechanism": True,
        "pseudo_update_recurrence_is_primary_word_structure": True,
        "accelerometer_each_valid_live_sample_consumed": bool(
            live["accelerometer_update_required_each_valid_imu_sample_after_live_entry"]
        ),
        "accelerometer_rejection_branch_consumed": False,
        "vector_PE_recurrence_consumed": True,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "old_P2_800_state_graph_consumed": False,
        "per_sample_Riccati_lower_propagation_used": False,
        "per_sample_SPD_lower_required": False,
        "selected_process_mode_strictness_used": False,
        "determinant_trace_scalarization_used": False,
        "scalar_information_beta_used": False,
        "gain_history_enumeration_used": False,
        "exact_measurement_dissipation_identity": {
            "identity": "V_minus - V_plus = r^T S_innov^-1 r",
            "S_pseudo_residual": "r_S0 = -delta_S",
            "S_pseudo_innovation": "S_innov = P_SS + R_S",
            "whole_state_action": "K_S uses P(:,S), so the S=0 correction acts through all learned cross-covariances",
        },
        "batch_innovation_identity": {
            "correction_information": "D_W = O_W^T Sigma_Y^-1 O_W",
            "word_decrease": "V_0 - V_W >= e_0^T D_W e_0",
            "gain_history_not_required": True,
        },
        "SEA3_coupled_schedule_contract": {
            "tau_applied_s": inv["tau_applied_s"],
            "sigma_aw_filter_mps2": inv["sigma_aw_filter_mps2"],
            "R_S_applied_base": inv["R_S_applied"],
            "pseudo_update_period_s": inv["pseudo_update_period_s"],
            "tau_and_active_pseudo_cadence_source_coupled": True,
            "SpectralMSE_target_tau_sigma_TS_coupled": True,
            "applied_RS_has_separate_EMA": True,
            "instantaneous_target_formula_substituted_for_applied_RS": False,
            "independent_applied_state_extrema_used_only_when_required_for_safe_bound": True,
            "physical_height_period_rectangular_extrema_forbidden": physical[
                "three_partition_contract"
            ]["independent_H_r_and_T_p_rectangular_extrema_forbidden"],
            "future_RS_tightening_requires_lag_reachability_theorem": True,
        },
        "translation_correction_word": {
            "mechanism": "FOUR_SEPARATED_S_ZERO_INNOVATIONS",
            "selected_windows_s": rsword["four_S_windows"],
            "word_horizon_s_upper": rsword["word_horizon_s_upper"],
            "dimensionless_state": rsword["dimensionless_state"],
            "time_varying_tau_allowed_inside_word": True,
            "four_S_observation_operator_full_rank": rsword[
                "four_S_translation_observation_operator_full_rank"
            ],
            "aw_third_divided_difference_lower": rsword[
                "aw_scaled_third_divided_difference_lower"
            ],
            "rank_witness_det_lower_not_gate": rsword[
                "scaled_observation_determinant_abs_lower_rank_witness_only"
            ],
            "selected_S_record_noise": rsword["selected_S_record_noise"],
            "accelerometer_needed_to_close_translation": False,
            "full_matrix_observation_coordinate_comparison_required": True,
        },
        "attitude_bias_correction_word": {
            "vector_PE_information_lower": vector["gyro_bias_two_packet"]["alpha_6_information_lower"],
            "accelerometer_and_magnetometer_information_reserved_for_coupled_attitude_bias_block": True,
            "process_UCC_used_only_for_metric_lower_not_as_primary_strictness": True,
        },
        "metric_scaling": {
            "process_UCC_covariance_lower_required": True,
            "translation_UCC_available": bool(trans["process_ucc"]["pass"]),
            "full_process_UCC_available": True,
            "preferred_translation_coordinates": "four-S observation/divided-difference coordinates",
            "target_inequality": "D_W >= delta * L_W^-1",
            "consequence": "P_0 >= L_W => D_W >= delta P_0^-1 => V_W <= (1-delta)V_0",
        },
        "endpoint_covariance_upper_retained_for_boundedness_and_P4": True,
        "covariance_memory": timing,
        "useful_gate": USEFUL_GATE,
        "P3_ARCHITECTURE_READY": True,
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED": bool(
            rsword["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"]
        ),
        "P3_RS_BATCH_NOISE_UPPER_CLOSED": bool(rsword["P3_RS_BATCH_NOISE_UPPER_CLOSED"]),
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED": False,
        "P3_UCC_METRIC_LOWER_CLOSED": False,
        "P3_FULL_MATRIX_COMPARISON_CLOSED": False,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "next_obligation": (
            "build the finite-memory translation covariance lower directly in the four-S observation/divided-difference coordinates, then compose the independent vector-PE attitude/gyro-bias block and A-mode bias block into one H18/A21 full-matrix gate"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "SEA3_RS_INNOVATION_DISSIPATION_WORD":
        f.append("wrong P3 architecture")
    for key in (
        "source_generated_not_trajectory_fit",
        "SEA3_dynamic_source_consumed",
        "SEA3_physical_parameter_coupling_consumed",
        "SEA3_full_finite_window_realization_assumed_conditionally_not_falsely_promoted",
        "R_S_is_primary_translation_correction_mechanism",
        "pseudo_update_recurrence_is_primary_word_structure",
        "accelerometer_each_valid_live_sample_consumed",
        "vector_PE_recurrence_consumed",
        "endpoint_covariance_upper_retained_for_boundedness_and_P4",
        "P3_ARCHITECTURE_READY",
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED",
        "P3_RS_BATCH_NOISE_UPPER_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "accelerometer_rejection_branch_consumed", "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed", "old_P2_800_state_graph_consumed",
        "per_sample_Riccati_lower_propagation_used", "per_sample_SPD_lower_required",
        "selected_process_mode_strictness_used", "determinant_trace_scalarization_used",
        "scalar_information_beta_used", "gain_history_enumeration_used",
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED", "P3_UCC_METRIC_LOWER_CLOSED",
        "P3_FULL_MATRIX_COMPARISON_CLOSED", "P3_CANONICAL_PASS", "P4_MAY_CONSUME_P3",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    c = d.get("SEA3_coupled_schedule_contract", {})
    if c.get("tau_and_active_pseudo_cadence_source_coupled") is not True:
        f.append("tau/pseudo cadence coupling lost")
    if c.get("SpectralMSE_target_tau_sigma_TS_coupled") is not True:
        f.append("SpectralMSE target coupling lost")
    if c.get("applied_RS_has_separate_EMA") is not True:
        f.append("applied R_S lag distinction lost")
    if c.get("instantaneous_target_formula_substituted_for_applied_RS") is not False:
        f.append("instantaneous target was substituted for applied R_S")
    t = d.get("translation_correction_word", {})
    if t.get("mechanism") != "FOUR_SEPARATED_S_ZERO_INNOVATIONS":
        f.append("translation is not using the four-S R_S word")
    if t.get("four_S_observation_operator_full_rank") is not True:
        f.append("four-S translation rank did not close")
    if t.get("accelerometer_needed_to_close_translation") is not False:
        f.append("accelerometer was incorrectly made necessary for translation rank")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--tube", type=Path, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.tube)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P3_architecture"],
        "R_S_primary": d["R_S_is_primary_translation_correction_mechanism"],
        "four_S_translation_rank": d["translation_correction_word"]["four_S_observation_operator_full_rank"],
        "S_record_information_scalar_lower": d["translation_correction_word"]["selected_S_record_noise"]["Sigma_S_inverse_scalar_lower"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "next_obligation": d["next_obligation"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
