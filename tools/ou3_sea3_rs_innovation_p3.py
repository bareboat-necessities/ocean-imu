#!/usr/bin/env python3
"""R_S translation-information component of the complete OU-III SEA3 P3 word.

This file is intentionally *not* a P3 architecture.  Canonical P3 is the
complete `SEA3_FULL_NORMAL_LIVE_RICCATI_WORD`.  This component supplies the
recurrent S=0 translation correction block that the full H18/A21 word must
consume together with accelerometer updates, asynchronous vector PE, process
UCC, covariance-floor events, A-mode bias dynamics, and the joint adaptive
source path.

For each exact Kalman correction,

    V^- - V^+ = r^T S_innov^-1 r,   V=e^T P^-1 e.

Four separated S=0 updates give full translation rank for arbitrary legal
time-varying tau.  `ou3_sea3_rs_word_information` now also closes a genuine
4x4 information lower in Newton divided-difference coordinates using the
inverse divided-difference-table norm, rather than determinant/Frobenius
singular-value conversion.

The remaining translation obligation is the finite-memory covariance lower in
the same path-dependent coordinates.  Nothing here promotes P3 or P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_riccati_tube_factored as TUBE
import ou3_sea3_rs_word_information as RSWORD

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_SEA3_RS_TRANSLATION_INFORMATION_COMPONENT"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("R_S translation component may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    physical = PHYSICAL.build(path)
    rsword = RSWORD.build(path)
    prereq = (
        [f"dynamic: {x}" for x in DYNAMIC.validate(dynamic)]
        + [f"physical: {x}" for x in PHYSICAL.validate(physical)]
        + [f"R_S word: {x}" for x in RSWORD.validate(rsword)]
    )
    if prereq:
        raise RuntimeError(f"R_S translation prerequisites failed: {prereq}")

    tube = TUBE.build(path) if tube_path is None else json.loads(
        Path(tube_path).read_text(encoding="utf-8")
    )
    tf = TUBE.validate(tube)
    if tf:
        raise RuntimeError(f"endpoint covariance prerequisite failed: {tf}")

    inv = dynamic["dynamic_invariant"]
    live = domain["normal_live"]
    ni = rsword["newton_coordinate_information"]

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "component_of_canonical_P3": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_consumed": True,
        "SEA3_physical_parameter_coupling_consumed": True,
        "R_S_is_primary_translation_correction_mechanism": True,
        "accelerometer_each_valid_live_sample_acknowledged_but_not_replaced_here": bool(
            live["accelerometer_update_required_each_valid_imu_sample_after_live_entry"]
        ),
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
            "whole_state_action": "K_S uses P(:,S); the full P(:,S) correction is retained by the complete word",
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
            "independent_tau_sigma_RS_TS_extrema_product_is_not_final_source_word": True,
            "physical_height_period_rectangular_extrema_forbidden": physical[
                "three_partition_contract"
            ]["independent_H_r_and_T_p_rectangular_extrema_forbidden"],
        },
        "translation_correction_word": {
            "mechanism": "FOUR_SEPARATED_S_ZERO_INNOVATIONS",
            "selected_windows_s": rsword["four_S_windows"],
            "word_horizon_s_upper": rsword["word_horizon_s_upper"],
            "dimensionless_state": rsword["dimensionless_state"],
            "newton_divided_difference_state": rsword["newton_divided_difference_state"],
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
            "newton_coordinate_information": ni,
            "D_S_newton_lambda_min_lower": ni["D_S_newton_lambda_min_lower"],
            "accelerometer_needed_to_close_translation_rank": False,
            "complete_word_must_still_retain_accelerometer_cross_block_information": True,
        },
        "endpoint_covariance_upper_retained_for_boundedness_and_P4": True,
        "covariance_memory": tube["covariance_memory"],
        "useful_gate": USEFUL_GATE,
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED": bool(
            rsword["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"]
        ),
        "P3_RS_BATCH_NOISE_UPPER_CLOSED": bool(
            rsword["P3_RS_BATCH_NOISE_UPPER_CLOSED"]
        ),
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED": bool(
            rsword["P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED"]
        ),
        "P3_UCC_METRIC_LOWER_CLOSED": False,
        "P3_FULL_MATRIX_COMPARISON_CLOSED": False,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "next_obligation": (
            "construct the finite-memory translation covariance lower in the same Newton coordinates while the complete 3 s H18/A21 producer carries all remaining Normal-Live preconditions"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("component_of_canonical_P3") != "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("R_S file is not scoped as a component of the full P3 word")
    for key in (
        "source_generated_not_trajectory_fit",
        "SEA3_dynamic_source_consumed",
        "SEA3_physical_parameter_coupling_consumed",
        "R_S_is_primary_translation_correction_mechanism",
        "accelerometer_each_valid_live_sample_acknowledged_but_not_replaced_here",
        "endpoint_covariance_upper_retained_for_boundedness_and_P4",
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED",
        "P3_RS_BATCH_NOISE_UPPER_CLOSED",
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed", "per_sample_Riccati_lower_propagation_used",
        "per_sample_SPD_lower_required", "selected_process_mode_strictness_used",
        "determinant_trace_scalarization_used", "scalar_information_beta_used",
        "gain_history_enumeration_used", "P3_UCC_METRIC_LOWER_CLOSED",
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
    if t.get("accelerometer_needed_to_close_translation_rank") is not False:
        f.append("accelerometer was incorrectly made necessary for translation rank")
    ni = t.get("newton_coordinate_information", {})
    if ni.get("full_4x4_matrix_inequality_closed") is not True:
        f.append("four-S Newton information matrix did not close")
    if float(t.get("D_S_newton_lambda_min_lower", 0.0)) <= 0.0:
        f.append("four-S Newton information lower is not strict")
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
        "component_of": d["component_of_canonical_P3"],
        "R_S_primary": d["R_S_is_primary_translation_correction_mechanism"],
        "four_S_translation_rank": d["translation_correction_word"]["four_S_observation_operator_full_rank"],
        "D_S_newton_lambda_lower": d["translation_correction_word"]["D_S_newton_lambda_min_lower"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "next_obligation": d["next_obligation"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
