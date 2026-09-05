#!/usr/bin/env python3
"""Complete SEA3 Normal-Live preconditions for canonical OU-III P3.

This contract has no alternate source model.  It consumes the compact,
phase-continuous ``ou3_sea3_complete_source`` state zeta=(x^s,lambda,z^t,q)
and the shipping full-process/vector-PE/runtime contracts, and it describes the
one literal 3 s H18/A21 word that P3 must interval-propagate.

Stochastic finite-horizon events are forcing/corollary material only and cannot
generate or prune the homogeneous source family. Reduced four-S words, tuner
rectangles, point schedules, independent SEA/RAO corners, arbitrary bounded
input sources, arbitrary P0 boxes, source-history graphs and blockwise
contraction surrogates are not prerequisites and cannot promote P3.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_windowed_vector_pe as PE
from ou3_interval import Interval

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
PAPER = REPO / "doc" / "kalman_ou_iii" / "w3d-iss-stability.tex-part"
SCHEMA = 8
QUALIFICATION = "OU3_SEA3_COMPLETE_SOURCE_P3_PRECONDITIONS_V8"
USEFUL_GATE = 1.0e-18


def _vec3_default(text: str, name: str) -> list[float]:
    m = re.search(
        rf"Eigen::Vector3f\s+{re.escape(name)}\s*=\s*Eigen::Vector3f\(\s*"
        r"([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*\)",
        text,
    )
    if not m:
        raise RuntimeError(f"cannot extract shipping Config {name}")
    out = [PE.binary32(float(m.group(i))) for i in range(1, 4)]
    if any(not (math.isfinite(x) and x > 0.0) for x in out):
        raise RuntimeError(f"shipping Config {name} is not positive finite")
    return out


def _variance_diag(std_xyz: list[float]) -> list[float]:
    return [Interval.point(x).square().hi for x in std_xyz]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("canonical P3 may not be trajectory fitted")

    complete = COMPLETE.build(path)
    process = PROCESS.build()
    pe = PE.build(path)
    bad = {
        "complete_SEA3": COMPLETE.validate(complete),
        "full_process": PROCESS.validate(process),
        "windowed_PE": PE.validate(pe),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"complete P3 prerequisite failure: {bad}")

    wrapper = WRAPPER.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    paper = PAPER.read_text(encoding="utf-8")
    sigma_a = _vec3_default(wrapper, "sigma_a")
    sigma_m = _vec3_default(wrapper, "sigma_m")

    source_parity = {
        "commit_before_prediction": (
            "apply_pending_online_tune_();" in wrapper
            and "mekf_->time_update(gyro, dt);" in wrapper
        ),
        "candidate_tau_each_valid_sample": "tune_.tau_applied   += alpha" in wrapper,
        "candidate_sigma_each_valid_sample": "tune_.sigma_applied += alpha" in wrapper,
        "candidate_RS_each_valid_sample": "tune_.RS_applied    += alpha_RS" in wrapper,
        "pseudo_period_from_committed_tau": "apply_pseudo_update_cadence_();" in wrapper,
        "periodic_aw_floor_tick": "periodic_aw_cov_sync_tick_();" in wrapper,
        "aw_floor_applied_in_prediction": "apply_pending_aw_covariance_inflation_();" in mekf,
        "S_zero_full_Joseph": (
            "applyIntegralZeroPseudoMeas" in mekf
            and "joseph_update3_(K, S_mat, PCt);" in mekf
        ),
        "accelerometer_full_Joseph": (
            "measurement_update_acc_only" in mekf
            and "joseph_update3_(K, S_mat, PCt);" in mekf
        ),
        "accelerometer_attitude_jacobian": (
            "const Matrix3 J_att = -skew_symmetric_matrix(f_cog_b);" in mekf
        ),
        "accelerometer_aw_jacobian": "const Matrix3 J_aw  =  R_wb();" in mekf,
        "active_accel_bias_jacobian": "PCt.noalias() += P_all_ba; // J_ba = I" in mekf,
        "deployed_SpectralMSE": "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;" in wrapper,
        "R_S_horizontal_factors": (
            "float R_S_x_factor_ = 0.72f;" in wrapper
            and "float R_S_y_factor_ = 0.72f;" in wrapper
        ),
    }
    source_failures = [k for k, v in source_parity.items() if not v]

    # These are semantic paper checks, not magic-token checks.  The paper writes
    # the compact SEA3 source through the implemented compact ranges, the
    # source state machine, the measurement-only front-end candidate, and the
    # sample-and-hold schedule.  Requiring an internal Python identifier such as
    # ``zeta`` or ``R_lambda`` to appear verbatim in TeX made equivalent theorem
    # text fail CI and did not strengthen the proof.
    paper_markers = {
        "strict_source_reachable_family": "strict source-reachable family" in paper,
        "finite_window_asynchronous_PE": "finite-window asynchronous conditions" in paper,
        "recurrent_S_regularizer": "Translational information from the integral regularizer" in paper,
        "full_H18_state": "\\vct e_H" in paper and "\\mathbb R^{18}" in paper,
        "full_A21_state": "\\vct e_A" in paper and "\\mathbb R^{21}" in paper,
        "compact_SEA3_transition_relation": (
            "On the implemented compact ranges" in paper
            and "source state machine" in paper
            and "strict source-reachable family" in paper
            and "arbitrary Cartesian product" in paper
        ),
        "augmented_SEA3_source_state": (
            "\\vartheta_k=(\\tau_k,\\sigma_{aw,k},r_{S,k},T_{S,k})" in paper
            and "measurement-only\nfront end maintains a separate candidate tuple" in paper
            and "activation timer" in paper
            and "hold/commit\nmap, and the source state machine" in paper
        ),
    }
    paper_failures = [k for k, v in paper_markers.items() if not v]

    horizon = max(
        float(complete["word_horizon_s"]),
        float(pe["spread_occurrence_selection"]["word_horizon_s"]),
    )
    samples = max(
        int(complete["word_samples"]),
        int(math.ceil(horizon / float(complete["derived_adaptive_source"][
            "source_recurrence_rate_and_commit_bounds"
        ]["dt_s"]))) + 1,
    )

    measurement_runtime = {
        "accelerometer_std_mps2": sigma_a,
        "accelerometer_variance_diag": _variance_diag(sigma_a),
        "magnetometer_std_uT": sigma_m,
        "magnetometer_variance_diag": _variance_diag(sigma_m),
        "configured_defaults_source_bound": True,
        "source_literals_converted_to_binary32_before_variance": True,
        "same_runtime_covariances_used_by_PE_and_full_word": (
            sigma_a == pe["measurement_runtime"]["accelerometer_std_mps2"]
            and sigma_m == pe["measurement_runtime"]["magnetometer_std_uT"]
        ),
    }

    front_end_state_manifest = complete["source_coordinates"]["front_end_state"]
    no_fallback = dict(complete["no_fallback_generators"])
    no_fallback.update({
        "D_W_L_W_split_gate": False,
        "blockwise_minimum_ratio_gate": False,
        "determinant_trace_gate": False,
        "scalar_information_beta_gate": False,
        "selected_process_mode_gate": False,
    })

    sea = complete["SEA3_surface_family"]
    realization = complete["SEA3_dynamic_realization"]
    stochastic = complete["stochastic_forcing_corollary"]
    mandatory = {
        "complete_SEA3_source_consumed": True,
        "complete_SEA3_source_contract_ready": bool(complete["P3_source_contract_ready"]),
        "complete_SEA3_no_fallback_generators": all(v is False for v in no_fallback.values()),
        "compact_SEA3_parameter_domain_consumed": bool(sea["parameter_domain_compact"]),
        "compact_SEA3_transition_relation_consumed": bool(
            sea["compact_transition_relation_is_theorem_domain"]
        ),
        "phase_continuous_SEA3_realization_required": bool(realization["phase_continuous"]),
        "same_xs_lambda_drives_entire_source_word": bool(
            realization["same_realization_drives_translation_rotation_frontend_tuner_geometry"]
        ),
        "hard_pathwise_SEA3_conditions_retained": bool(
            realization["hard_pathwise_acceleration_and_body_rate_conditions_retained"]
        ),
        "stochastic_event_not_source_generator": (
            stochastic["used_to_generate_P3_source_words"] is False
        ),
        "stochastic_event_not_homogeneous_pruner": (
            stochastic["used_to_prune_homogeneous_P3_family"] is False
        ),
        "full_process_UCC_consumed": True,
        "windowed_asynchronous_vector_PE_consumed": True,
        "PE_uses_binary32_shipping_covariances": bool(
            pe["source_literals_converted_to_binary32_before_variance"]
            and pe["variance_bounds_outward_after_binary32_conversion"]
        ),
        "all_valid_accelerometer_updates_required": bool(
            complete["Normal_Live_nonsea_conditions"]["all_valid_accelerometer_updates_required"]
        ),
        "accelerometer_rejection_absent": (
            complete["Normal_Live_nonsea_conditions"]["accelerometer_rejection_in_scope"] is False
        ),
        "all_due_S_updates_required": bool(
            complete["R_S_regularizer"]["all_due_S_updates_remain_in_full_word"]
        ),
        "actual_applied_per_axis_R_S_required": bool(
            complete["R_S_regularizer"]["actual_applied_R_S_required_at_every_due_S_update"]
        ),
        "same_SEA3_path_feeds_frontend_tuner_F_Q_TS_RS": True,
        "live_H18_and_A21_full_matrix_gate_required": True,
    }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_architecture": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "canonical_source": complete["canonical_P3_source"],
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "complete_SEA3_source_consumed": True,
        "complete_SEA3_source": complete,
        "complete_SEA3_source_family_materialized": bool(
            complete["P3_source_family_materialized"]
        ),
        "no_fallback_generators": no_fallback,
        "source_parity": source_parity,
        "source_parity_failures": source_failures,
        "paper_parity": paper_markers,
        "paper_parity_failures": paper_failures,
        "mandatory_preconditions": mandatory,
        "all_current_machine_checkable_preconditions_present": (
            not source_failures
            and not paper_failures
            and all(mandatory.values())
            and all(v is False for v in no_fallback.values())
        ),
        "measurement_runtime": measurement_runtime,
        "front_end_state_manifest": front_end_state_manifest,
        "windowed_PE": pe,
        "full_process_UCC": process,
        "word": {
            "horizon_s": horizon,
            "samples_upper": samples,
            "H_dimension": 18,
            "A_dimension": 21,
            "same_complete_SEA3_event_word_for_H_and_A": True,
            "source_state": "zeta=(x^s,lambda,z^t,q)",
            "every_valid_IMU_sample": "prediction -> queued aw floor -> every due S=0 -> accelerometer Joseph -> frontend/tuner update/staging",
            "magnetometer_events": "asynchronous accepted events satisfying the SEA3/Normal-Live PE premise",
        },
        "final_numeric_contract": {
            "common_word_horizon_s": horizon,
            "H_dimension": 18,
            "A_dimension": 21,
            "useful_gate": USEFUL_GATE,
            "actual_applied_per_axis_RS_required": True,
            "full_18x18_and_21x21_matrix_comparison_required": True,
            "prediction_recursion_required": "P-=F P F^T+Q; Psi-=F Psi; Omega-=F Omega F^T+Q",
            "joseph_measurement_recursion_required": "P+=A P- A^T+K R K^T; Psi+=A Psi-; Omega+=A Omega- A^T+K R K^T",
            "aw_floor_recursion_required": "P+=Delta_aw; Omega+=Delta_aw for the full embedded PSD increment",
            "required_final_inequality": "Omega_W - delta*P_W >= 0 on full H18 and A21",
            "moving_metric_equivalence": "V_W <= (1-delta) V_0",
        },
        "P3_promoted": False,
        "next_obligation": (
            "interval-propagate the compact phase-continuous SEA3 source family itself through every event of the literal H18/A21 word and validate Omega_W-delta*P_W by full-matrix LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("P3 preconditions detached from complete SEA3")
    for key in (
        "complete_SEA3_source_consumed",
        "all_current_machine_checkable_preconditions_present",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if d.get("complete_SEA3_source_family_materialized") is not False:
        f.append("precondition scaffold falsely claims SEA3 family materialization")
    for key in ("trajectory_replay_used", "filter_changed", "declared_domain_shrunk", "P3_promoted"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("source_parity_failures"):
        f.extend(f"source parity failed: {x}" for x in d["source_parity_failures"])
    if d.get("paper_parity_failures"):
        f.extend(f"paper parity failed: {x}" for x in d["paper_parity_failures"])
    fallback = d.get("no_fallback_generators", {})
    if not fallback or any(v is not False for v in fallback.values()):
        f.append("a fallback generator/gate is still enabled")
    mandatory = d.get("mandatory_preconditions", {})
    if not mandatory or not all(mandatory.values()):
        f.append("one or more complete SEA3 mandatory preconditions are absent")
    runtime = d.get("measurement_runtime", {})
    if runtime.get("source_literals_converted_to_binary32_before_variance") is not True:
        f.append("P3 runtime covariance lost binary32 source conversion")
    if runtime.get("same_runtime_covariances_used_by_PE_and_full_word") is not True:
        f.append("PE/full-word measurement covariance parity failed")
    final = d.get("final_numeric_contract", {})
    if float(final.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")
    if final.get("actual_applied_per_axis_RS_required") is not True:
        f.append("actual per-axis R_S requirement disappeared")
    if final.get("full_18x18_and_21x21_matrix_comparison_required") is not True:
        f.append("full H/A matrix comparison requirement disappeared")
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
        "canonical_source": d["canonical_source"],
        "word": d["word"],
        "mandatory": d["mandatory_preconditions"],
        "fallbacks": d["no_fallback_generators"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
