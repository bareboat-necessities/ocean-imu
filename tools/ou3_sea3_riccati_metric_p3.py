#!/usr/bin/env python3
"""Canonical OU-III P3 gate over SEA3/R_S innovation dissipation.

The historical file name is retained because P4 imports it.  Strictness now
comes from the exact finite-word innovation-energy identity, with recurrent
S=0 pseudo updates and SEA3-coupled R_S/T_S as the primary translation
correction mechanism.  Process UCC supplies the lower covariance scale needed
to compare correction information with the moving P^-1 metric; it is not the
primary source of contraction.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_rs_innovation_p3 as ARCH

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 4
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
        "exact_measurement_dissipation_identity_consumed": True,
        "batch_innovation_information_identity_consumed": True,
        "process_UCC_used_as_metric_lower_not_primary_strictness": True,
        "same_operating_point_tau_sigma_RS_TS_required": True,
        "independent_source_extrema_product_forbidden": True,
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
        "SEA3_coupled_schedule_contract": a["SEA3_coupled_schedule_contract"],
        "translation_correction_word": a["translation_correction_word"],
        "attitude_bias_correction_word": a["attitude_bias_correction_word"],
        "metric_scaling": a["metric_scaling"],
        "modes": modes,
        "useful_gate": USEFUL_GATE,
        "P3_FOUNDATION_PASS": True,
        "P3_ARCHITECTURE_READY": True,
        "P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED": a["P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED"],
        "P3_UCC_METRIC_LOWER_CLOSED": a["P3_UCC_METRIC_LOWER_CLOSED"],
        "P3_FULL_MATRIX_COMPARISON_CLOSED": a["P3_FULL_MATRIX_COMPARISON_CLOSED"],
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "P3_CANONICAL_FAIL_REASONS": [
            "SEA3-coupled R_S/T_S weighted batch innovation matrix D_W is not yet emitted",
            "source-uniform UCC covariance lower L_W is not yet emitted in the same coordinates",
            "full H18/A21 matrix inequality D_W >= delta L_W^-1 has not yet been validated",
        ],
        "next_obligation": (
            "close D_W and L_W for this fixed R_S/SEA3 architecture and run the unchanged 1e-18 full-matrix gate; do not introduce another P3 architecture"
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
        "exact_measurement_dissipation_identity_consumed",
        "batch_innovation_information_identity_consumed",
        "process_UCC_used_as_metric_lower_not_primary_strictness",
        "same_operating_point_tau_sigma_RS_TS_required",
        "independent_source_extrema_product_forbidden",
        "P3_FOUNDATION_PASS", "P3_ARCHITECTURE_READY",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
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
    for mode, dim in (("H", 18), ("A", 21)):
        row = d.get("modes", {}).get(mode, {})
        if row.get("dimension") != dim or row.get("pass") is not False:
            f.append(f"{mode} fail-closed mode contract invalid")
        if float(row.get("relative_Riccati_injection_margin_lower", math.nan)) != 0.0:
            f.append(f"{mode} emitted a numerical margin before D_W/L_W closure")
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
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
