#!/usr/bin/env python3
"""Canonical OU-III P3 architecture: SEA3 + recurrent R_S innovation dissipation.

The dominant stabilizing mechanism is the shipping measurement correction, not
one-sample process-noise injection.  For every exact linear Kalman measurement
update with prior error e-, innovation r=H e-, innovation covariance
S=H P- H' + R, posterior covariance P+, and homogeneous posterior error e+,

    V- - V+ = r' S^-1 r,
    V=e' P^-1 e.

Prediction is non-expansive because

    P- = F P+ F' + Q,  Q>=0.

Hence over a recurrent word

    V_0 - V_W >= sum_i r_i' S_i^-1 r_i.

The sequential innovations are the block-Cholesky whitening of the complete
batch measurement record.  Therefore the summed correction energy equals

    e_0' D_W e_0,
    D_W = O_W' Sigma_Y^-1 O_W,

where O_W is the finite-word observation operator from the initial error and
Sigma_Y is the joint covariance of process/measurement nuisance in that word.
This formulation eliminates gain-history enumeration completely.

R_S is central.  The S=0 update has H_S selecting the integral-displacement
state and innovation covariance P_SS+R_S.  The shipping Joseph update uses the
whole cross-covariance P(:,S), so one S correction acts on the entire correlated
[v,p,S,a_w] chain and can also nudge attitude through cross-covariance.  SEA3
supplies the recurrence and the source coupling needed to quantify that effect.

The proof must NOT replace the adaptive source by independent Cartesian extrema.
The shipping schedule couples

  * tau to the measured wave period,
  * T_S to tau through the bounded tau-scaled pseudo cadence,
  * r_S to the selected law evaluated at the same tau/sigma/T_S operating point,
  * the Cubic law to C_R*sqrt(R_a)*tau^3 with cadence information-rate matching,
  * other Riccati/MSE laws to their explicit same-operating-point formulas.

Thus pathological combinations such as maximum r_S with unrelated minimum
process/sea scale are inadmissible unless the actual SEA3+tuner dynamics permit
them.  This is the main place where SEA3 strengthens the stability certificate.

Quantitative P3 will certify two full matrices over one recurrent word:

  L_W  <= P_0                       (UCC covariance lower),
  D_W  <=/=> correction information (lower bound on O' Sigma_Y^-1 O).

The useful contraction gate is the full-matrix inequality

    D_W >= delta * L_W^-1,     delta >= 1e-18,

which implies D_W >= delta P_0^-1 and therefore

    V_W <= (1-delta) V_0.

UCC is retained only to scale the moving covariance metric.  Strictness comes
from recurrent measurement correction, principally S=0 plus vector/accelerometer
information.  No determinant/trace scalarization, per-sample SPD lower, selected
process-mode posterior attenuation, source-history graph, or predecessor path
is canonical.
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
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
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
    prereq = (
        [f"dynamic: {x}" for x in DYNAMIC.validate(dynamic)]
        + [f"physical: {x}" for x in PHYSICAL.validate(physical)]
        + [f"vector: {x}" for x in VECTOR.validate(vector)]
        + [f"process: {x}" for x in PROCESS.validate(process)]
        + [f"translation: {x}" for x in TRANS.validate(trans)]
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
            "same_operating_point_required_for_tau_sigma_RS_TS": True,
            "independent_extrema_product_forbidden": True,
            "physical_height_period_rectangular_extrema_forbidden": physical[
                "three_partition_contract"
            ]["independent_H_r_and_T_p_rectangular_extrema_forbidden"],
            "R_S_law_role": (
                "consume the shipping selected R_S law and cadence normalization at the same SEA3/tuner operating point; "
                "do not use max R_S independently of tau/sigma/T_S"
            ),
        },
        "translation_correction_word": {
            "endpoint_lag_windows": timing.get("S_observation_window_layout", "[0,g],[2g,3g],[4g,5g]"),
            "pseudo_gap_s_upper": timing["pseudo_gap_s_upper"],
            "three_separated_S_observations_reconstruct_v_p_S": True,
            "stable_aw_plus_accelerometer_closes_fourth_translation_direction": True,
            "full_matrix_weighted_observability_required": True,
        },
        "attitude_bias_correction_word": {
            "vector_PE_information_lower": vector["gyro_bias_two_packet"]["alpha_6_information_lower"],
            "accelerometer_and_magnetometer_information_used_as_correction": True,
            "process_UCC_used_only_for_metric_lower_not_as_primary_strictness": True,
        },
        "metric_scaling": {
            "process_UCC_covariance_lower_required": True,
            "translation_UCC_available": bool(trans["process_ucc"]["pass"]),
            "full_process_UCC_available": True,
            "target_inequality": "D_W >= delta * L_W^-1",
            "consequence": "P_0 >= L_W => D_W >= delta P_0^-1 => V_W <= (1-delta)V_0",
        },
        "endpoint_covariance_upper_retained_for_boundedness_and_P4": True,
        "covariance_memory": timing,
        "useful_gate": USEFUL_GATE,
        "P3_ARCHITECTURE_READY": True,
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED": False,
        "P3_UCC_METRIC_LOWER_CLOSED": False,
        "P3_FULL_MATRIX_COMPARISON_CLOSED": False,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "next_obligation": (
            "derive the SEA3-coupled R_S/T_S weighted batch innovation matrix D_W and the UCC covariance lower L_W, then run one full-matrix D_W >= delta L_W^-1 gate; no alternative P3 architecture is permitted"
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
    if c.get("same_operating_point_required_for_tau_sigma_RS_TS") is not True:
        f.append("tau/sigma/R_S/T_S coupling lost")
    if c.get("independent_extrema_product_forbidden") is not True:
        f.append("independent SEA3/tuner extrema reintroduced")
    t = d.get("translation_correction_word", {})
    if t.get("three_separated_S_observations_reconstruct_v_p_S") is not True:
        f.append("three-S translation correction structure lost")
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
        "SEA3_coupling": d["SEA3_coupled_schedule_contract"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "next_obligation": d["next_obligation"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
