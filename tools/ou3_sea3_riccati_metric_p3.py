#!/usr/bin/env python3
"""Canonical OU-III P3 gate: complete SEA3 -> full H18/A21 Riccati word.

There is one promotable route only.

1. ``ou3_sea3_complete_source`` defines the compact, phase-continuous admitted
   SEA3 Normal-Live source word zeta=(x^s,lambda,z^t,q).
2. ``ou3_sea3_p3_full_preconditions`` binds the shipping runtime and theorem
   conditions to that same source word.
3. ``ou3_sea3_full_normal_live_word`` provides the literal full-state Joseph /
   prediction / covariance-floor operations.
4. The SEA3 source-family executor must interval-propagate every complete SEA3
   word for 3 s and establish, for both modes,

       Omega_W - delta P_W >= 0,    delta >= 1e-18.

A stochastic finite-horizon concentration event is not a source generator and
cannot prune the homogeneous P3 family; stochastic sensor/model perturbations
remain forcing while configured Racc/Rmag remain in the Riccati word.

No four-S reduced word, point executor, tuner rectangle, independent sea/RAO
corner, arbitrary bounded-input source, D_W/L_W split, blockwise ratio, scalar
information beta, determinant/trace surrogate, source-history graph, arbitrary
P0 box or selected process strictness is accepted by this gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_full_word_riccati_backend as BACKEND
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_p3_full_preconditions as FULL

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 13
QUALIFICATION = "OU3_SEA3_COMPLETE_SOURCE_FULL_WORD_P3_GATE_V13"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    full = FULL.build(path)
    live = LIVE.build(path)
    word = WORD.build(path)
    failures = {
        "complete_SEA3": COMPLETE.validate(complete),
        "preconditions": FULL.validate(full),
        "live_seed": LIVE.validate(live),
        "literal_word": WORD.validate(word),
        "backend": BACKEND.validate_backend(),
    }
    failures = {k: v for k, v in failures.items() if v}
    if failures:
        raise RuntimeError(f"canonical P3 prerequisites failed: {failures}")

    fallback = dict(complete["no_fallback_generators"])
    fallback.update(full["no_fallback_generators"])
    if any(v is not False for v in fallback.values()):
        raise RuntimeError("canonical P3 gate found an enabled fallback route")

    sea = complete["SEA3_surface_family"]
    realization = complete["SEA3_dynamic_realization"]
    stochastic = complete["stochastic_forcing_corollary"]

    word_family_materialized = bool(word["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
    h_executed = bool(word["FULL_H18_WORD_EXECUTED"])
    a_executed = bool(word["FULL_A21_WORD_EXECUTED"])
    ldlt_closed = bool(word["FULL_H18_A21_LDLT_CLOSED"])
    p3_pass = bool(word["P3_CANONICAL_PASS"])
    if p3_pass and not (word_family_materialized and h_executed and a_executed and ldlt_closed):
        raise RuntimeError("P3 claimed PASS without complete SEA3 full-word closure")

    modes = {
        "H18": {
            "dimension": 18,
            "full_word_executed": h_executed,
            "Omega_minus_delta_P_ldlt_closed": ldlt_closed and h_executed,
            "certified_delta_lower": 0.0 if not p3_pass else USEFUL_GATE,
        },
        "A21": {
            "dimension": 21,
            "full_word_executed": a_executed,
            "Omega_minus_delta_P_ldlt_closed": ldlt_closed and a_executed,
            "certified_delta_lower": 0.0 if not p3_pass else USEFUL_GATE,
        },
    }

    fail_reasons = [] if p3_pass else [
        "the compact phase-continuous SEA3 3 s source family has not yet been interval-propagated through every literal H18/A21 shipping event",
        "the same SEA3 zeta=(x^s,lambda,z^t,q) word must generate physical response, exact front-end state, tuner targets, EMA/commit path, T_S, Q, actual applied per-axis R_S, accelerometer/vector geometry, magnetic PE events and covariance-floor events",
        "the full 18x18 and 21x21 matrices Omega_W-delta*P_W have not yet both passed validated LDLT at delta=1e-18",
    ]

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "COMPLETE_SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "canonical_source": complete["canonical_P3_source"],
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "useful_gate": USEFUL_GATE,
        "complete_SEA3_source_consumed": True,
        "complete_SEA3_source_validation_pass": True,
        "complete_SEA3_response_couplings_consumed": True,
        "complete_SEA3_compact_parameter_domain_consumed": bool(
            sea["parameter_domain_compact"]
        ),
        "complete_SEA3_compact_transition_relation_consumed": bool(
            sea["compact_transition_relation_is_theorem_domain"]
        ),
        "complete_SEA3_phase_continuous_realization_required": bool(
            realization["phase_continuous"]
        ),
        "same_xs_lambda_drives_entire_word": bool(
            realization["same_realization_drives_translation_rotation_frontend_tuner_geometry"]
        ),
        "hard_pathwise_SEA3_conditions_retained": bool(
            realization["hard_pathwise_acceleration_and_body_rate_conditions_retained"]
        ),
        "stochastic_forcing_does_not_generate_source_words": (
            stochastic["used_to_generate_P3_source_words"] is False
        ),
        "stochastic_forcing_does_not_prune_homogeneous_family": (
            stochastic["used_to_prune_homogeneous_P3_family"] is False
        ),
        "configured_measurement_covariance_retained_under_stochastic_forcing": bool(
            stochastic["configured_Racc_Rmag_remain_in_every_covariance_update"]
        ),
        "complete_SEA3_frontend_state_consumed": True,
        "complete_SEA3_adaptive_state_consumed": True,
        "actual_applied_per_axis_RS_consumed": True,
        "all_due_S_updates_required": True,
        "all_valid_accelerometer_updates_required": True,
        "asynchronous_vector_PE_required": True,
        "all_full_process_Q_required": True,
        "all_aw_covariance_floor_events_required": True,
        "same_complete_SEA3_word_used_for_H18_and_A21": True,
        "live_entry_covariance_seed_consumed": True,
        "live_entry_covariance_seed_source_generated": bool(
            live["live_entry_seed_is_source_generated_not_arbitrary_PSD"]
        ),
        "joint_P_Psi_Omega_backend_consumed": True,
        "joint_backend_validation_pass": True,
        "literal_full_word_assembler_consumed": True,
        "literal_full_word_assembler_validation_pass": True,
        "literal_shipping_event_order_pass": bool(word["shipping_event_order_parity_pass"]),
        "no_fallback_generators": fallback,
        "no_fallback_route_enabled": all(v is False for v in fallback.values()),
        "independent_tau_sigma_RS_TS_extrema_product_used": False,
        "independent_sea_x_RAO_product_used": False,
        "point_source_word_used": False,
        "selected_four_S_word_used": False,
        "gaussian_good_event_source_used": False,
        "spectral_moment_only_source_used": False,
        "arbitrary_bounded_input_source_used": False,
        "D_W_L_W_split_used": False,
        "blockwise_minimum_ratio_used": False,
        "scalar_information_beta_used": False,
        "determinant_trace_scalarization_used": False,
        "source_history_graph_used": False,
        "predecessor_path_enumeration_used": False,
        "arbitrary_P0_rectangle_used": False,
        "selected_process_mode_strictness_used": False,
        "common_word_horizon_s": float(full["final_numeric_contract"]["common_word_horizon_s"]),
        "required_final_inequality": full["final_numeric_contract"]["required_final_inequality"],
        "moving_metric_equivalence": full["final_numeric_contract"]["moving_metric_equivalence"],
        "complete_preconditions": full,
        "literal_full_word": word,
        "modes": modes,
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED": word_family_materialized,
        "P3_FULL_WORD_ENCLOSED": word_family_materialized and h_executed and a_executed,
        "P3_FULL_MATRIX_COMPARISON_CLOSED": ldlt_closed,
        "P3_CANONICAL_PASS": p3_pass,
        "P4_MAY_CONSUME_P3": p3_pass,
        "P3_CANONICAL_FAIL_REASONS": fail_reasons,
        "next_obligation": (
            "execute the compact phase-continuous SEA3 forward family itself and close both full-matrix LDLTs"
            if not p3_pass else "P3 closed; P4 may consume the certified complete-SEA3 word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "COMPLETE_SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("canonical P3 architecture changed")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical P3 source is not complete SEA3")
    for key in (
        "source_generated_not_trajectory_fit",
        "complete_SEA3_source_consumed",
        "complete_SEA3_source_validation_pass",
        "complete_SEA3_response_couplings_consumed",
        "complete_SEA3_compact_parameter_domain_consumed",
        "complete_SEA3_compact_transition_relation_consumed",
        "complete_SEA3_phase_continuous_realization_required",
        "same_xs_lambda_drives_entire_word",
        "hard_pathwise_SEA3_conditions_retained",
        "stochastic_forcing_does_not_generate_source_words",
        "stochastic_forcing_does_not_prune_homogeneous_family",
        "configured_measurement_covariance_retained_under_stochastic_forcing",
        "complete_SEA3_frontend_state_consumed",
        "complete_SEA3_adaptive_state_consumed",
        "actual_applied_per_axis_RS_consumed",
        "all_due_S_updates_required",
        "all_valid_accelerometer_updates_required",
        "asynchronous_vector_PE_required",
        "all_full_process_Q_required",
        "all_aw_covariance_floor_events_required",
        "same_complete_SEA3_word_used_for_H18_and_A21",
        "live_entry_covariance_seed_consumed",
        "live_entry_covariance_seed_source_generated",
        "joint_P_Psi_Omega_backend_consumed",
        "joint_backend_validation_pass",
        "literal_full_word_assembler_consumed",
        "literal_full_word_assembler_validation_pass",
        "literal_shipping_event_order_pass",
        "no_fallback_route_enabled",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "independent_tau_sigma_RS_TS_extrema_product_used",
        "independent_sea_x_RAO_product_used", "point_source_word_used",
        "selected_four_S_word_used", "gaussian_good_event_source_used",
        "spectral_moment_only_source_used", "arbitrary_bounded_input_source_used",
        "D_W_L_W_split_used", "blockwise_minimum_ratio_used",
        "scalar_information_beta_used", "determinant_trace_scalarization_used",
        "source_history_graph_used", "predecessor_path_enumeration_used",
        "arbitrary_P0_rectangle_used", "selected_process_mode_strictness_used",
    ):
        if d.get(key) is not False:
            f.append(f"forbidden fallback route enabled: {key}")
    fallback = d.get("no_fallback_generators", {})
    if not fallback or any(v is not False for v in fallback.values()):
        f.append("fallback generator map is not entirely disabled")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")
    p = bool(d.get("P3_CANONICAL_PASS"))
    if p:
        if d.get("SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED") is not True:
            f.append("P3 PASS without complete SEA3 source-family materialization")
        if d.get("P3_FULL_WORD_ENCLOSED") is not True:
            f.append("P3 PASS without full word enclosure")
        if d.get("P3_FULL_MATRIX_COMPARISON_CLOSED") is not True:
            f.append("P3 PASS without full-matrix closure")
    else:
        if not d.get("P3_CANONICAL_FAIL_REASONS"):
            f.append("fail-closed P3 gate has no reason")
    if d.get("P4_MAY_CONSUME_P3") is not p:
        f.append("P4 promotion does not exactly follow P3")
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
        "fallbacks": d["no_fallback_generators"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "P3_CANONICAL_FAIL_REASONS": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
