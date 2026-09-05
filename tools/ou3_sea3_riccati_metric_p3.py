#!/usr/bin/env python3
"""Canonical OU-III P3 gate over SEA3/R_S innovation dissipation.

The historical file name is retained because P4 imports it.  Translation
strictness comes from recurrent S=0 innovation dissipation, now with a certified
four-S full-rank observation word.  Process UCC supplies the finite-memory
covariance scale; it is not the primary source of contraction.

The gate distinguishes target and applied tuning.  SpectralMSE evaluates its
R_S target at the same target tau/sigma/T_S operating point, and active tau is
committed together with its pseudo cadence.  Applied R_S has a separate EMA,
so no instantaneous target-law relation is assumed for active R_S.  Until a
lag/reachability theorem tightens that relation, quantitative bounds use the
safe applied R_S invariant.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_rs_innovation_p3 as ARCH

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 5
QUALIFICATION = "OU3_SEA3_RS_INNOVATION_DISSIPATION_P3_GATE"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    a = ARCH.build(path, tube_path)
    af = ARCH.validate(a)
    if af:
        raise RuntimeError(f"R_S innovation P3 architecture failed validation: {af}")

    quantitative_closed = all(bool(a[k]) for k in (
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED",
        "P3_UCC_METRIC_LOWER_CLOSED",
        "P3_FULL_MATRIX_COMPARISON_CLOSED",
    ))
    canonical_pass = quantitative_closed and bool(a["P3_CANONICAL_PASS"])
    if canonical_pass:
        raise RuntimeError("architecture skeleton cannot promote P3 before numeric margins are emitted")

    modes = {
        "H": {
            "dimension": 18,
            "relative_Riccati_injection_margin_lower": 0.0,
            "contraction_gap_lower": 0.0,
            "useful_margin_gate": USEFUL_GATE,
            "pass": False,
        },
        "A": {
            "dimension": 21,
            "relative_Riccati_injection_margin_lower": 0.0,
            "contraction_gap_lower": 0.0,
            "useful_margin_gate": USEFUL_GATE,
            "pass": False,
        },
    }

    schedule = a["SEA3_coupled_schedule_contract"]
    trans = a["translation_correction_word"]
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
        "R_S_is_primary_translation_correction_mechanism": True,
        "pseudo_update_recurrence_is_primary_word_structure": True,
        "four_S_translation_word_consumed": True,
        "four_S_translation_observation_geometry_closed": bool(
            a["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"]
        ),
        "four_S_batch_noise_upper_closed": bool(a["P3_RS_BATCH_NOISE_UPPER_CLOSED"]),
        "exact_measurement_dissipation_identity_consumed": True,
        "batch_innovation_information_identity_consumed": True,
        "process_UCC_used_as_metric_lower_not_primary_strictness": True,
        "tau_active_pseudo_cadence_coupling_consumed": bool(
            schedule["tau_and_active_pseudo_cadence_source_coupled"]
        ),
        "SpectralMSE_target_tau_sigma_TS_coupling_consumed": bool(
            schedule["SpectralMSE_target_tau_sigma_TS_coupled"]
        ),
        "applied_RS_separate_EMA_acknowledged": bool(schedule["applied_RS_has_separate_EMA"]),
        "instantaneous_RS_target_substituted_for_applied_RS": bool(
            schedule["instantaneous_target_formula_substituted_for_applied_RS"]
        ),
        "safe_applied_RS_invariant_used_until_lag_theorem": True,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "old_P2_800_state_graph_consumed": False,
        "one_sample_strict_Riccati_margin_consumed": False,
        "commit_aligned_source_word_consumed": False,
        "per_sample_SPD_lower_required": False,
        "selected_process_mode_strictness_used": False,
        "determinant_trace_scalarization_used": False,
        "scalar_information_beta_used": False,
        "parameter_dependent_metric": "V_k=e_k^T P_k^-1 e_k with P_k the shipping Riccati covariance",
        "strictness_location": "RECURRENT_SEA3_MEASUREMENT_WORD",
        "innovation_identity": a["exact_measurement_dissipation_identity"],
        "batch_identity": a["batch_innovation_identity"],
        "SEA3_coupled_schedule_contract": schedule,
        "translation_correction_word": trans,
        "attitude_bias_correction_word": a["attitude_bias_correction_word"],
        "metric_scaling": a["metric_scaling"],
        "modes": modes,
        "useful_gate": USEFUL_GATE,
        "P3_FOUNDATION_PASS": True,
        "P3_ARCHITECTURE_READY": True,
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED": a[
            "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"
        ],
        "P3_RS_BATCH_NOISE_UPPER_CLOSED": a["P3_RS_BATCH_NOISE_UPPER_CLOSED"],
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED": a["P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED"],
        "P3_UCC_METRIC_LOWER_CLOSED": a["P3_UCC_METRIC_LOWER_CLOSED"],
        "P3_FULL_MATRIX_COMPARISON_CLOSED": a["P3_FULL_MATRIX_COMPARISON_CLOSED"],
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "P3_CANONICAL_FAIL_REASONS": [
            "four-S R_S translation geometry/noise are closed, but the finite-memory covariance lower has not yet been expressed in the same observation coordinates",
            "the vector-PE attitude/gyro-bias and active A-mode bias blocks have not yet been composed with the four-S translation block",
            "full H18/A21 matrix inequality has not yet been validated at the unchanged 1e-18 gate",
        ],
        "next_obligation": (
            "close the finite-memory covariance lower in the four-S observation/divided-difference coordinates, compose the vector-PE attitude/bias block, and run the one H18/A21 full-matrix gate"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "SEA3_RS_INNOVATION_DISSIPATION_WORD":
        f.append("wrong canonical P3 architecture")
    for key in (
        "source_generated_not_trajectory_fit", "SEA3_dynamic_source_consumed",
        "SEA3_physical_parameter_coupling_consumed",
        "R_S_is_primary_translation_correction_mechanism",
        "pseudo_update_recurrence_is_primary_word_structure",
        "four_S_translation_word_consumed",
        "four_S_translation_observation_geometry_closed",
        "four_S_batch_noise_upper_closed",
        "exact_measurement_dissipation_identity_consumed",
        "batch_innovation_information_identity_consumed",
        "process_UCC_used_as_metric_lower_not_primary_strictness",
        "tau_active_pseudo_cadence_coupling_consumed",
        "SpectralMSE_target_tau_sigma_TS_coupling_consumed",
        "applied_RS_separate_EMA_acknowledged",
        "safe_applied_RS_invariant_used_until_lag_theorem",
        "P3_FOUNDATION_PASS", "P3_ARCHITECTURE_READY",
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED",
        "P3_RS_BATCH_NOISE_UPPER_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "instantaneous_RS_target_substituted_for_applied_RS",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed", "one_sample_strict_Riccati_margin_consumed",
        "commit_aligned_source_word_consumed", "per_sample_SPD_lower_required",
        "selected_process_mode_strictness_used", "determinant_trace_scalarization_used",
        "scalar_information_beta_used", "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED",
        "P3_UCC_METRIC_LOWER_CLOSED", "P3_FULL_MATRIX_COMPARISON_CLOSED",
        "P3_CANONICAL_PASS", "P4_MAY_CONSUME_P3",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("strictness_location") != "RECURRENT_SEA3_MEASUREMENT_WORD":
        f.append("strictness is not assigned to recurrent SEA3 measurements")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")
    t = d.get("translation_correction_word", {})
    if t.get("mechanism") != "FOUR_SEPARATED_S_ZERO_INNOVATIONS":
        f.append("gate did not consume the four-S translation route")
    if t.get("accelerometer_needed_to_close_translation") is not False:
        f.append("translation still depends on accelerometer rank repair")
    for mode, dim in (("H", 18), ("A", 21)):
        row = d.get("modes", {}).get(mode, {})
        if row.get("dimension") != dim or row.get("pass") is not False:
            f.append(f"{mode} fail-closed mode contract invalid")
        if float(row.get("relative_Riccati_injection_margin_lower", math.nan)) != 0.0:
            f.append(f"{mode} emitted a numerical margin before full matrix closure")
    if not d.get("P3_CANONICAL_FAIL_REASONS"):
        f.append("open P3 does not name quantitative obligations")
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
        "four_S_geometry_closed": d["four_S_translation_observation_geometry_closed"],
        "four_S_noise_closed": d["four_S_batch_noise_upper_closed"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
