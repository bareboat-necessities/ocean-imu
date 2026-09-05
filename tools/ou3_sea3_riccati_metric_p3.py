#!/usr/bin/env python3
"""Canonical OU-III P3 gate over the lifted finite-word architecture.

The file name is retained because P4 imports it, but the canonical proof no
longer obtains strictness from a one-sample Riccati process floor.  P3 now uses
``ou3_sea3_lifted_word_p3``:

* every shipping sample is non-expansive in the moving covariance metric;
* strictness is a recurrent finite-word statement;
* the word injection is lower-bounded by a finite selected process-noise basis
  in lifted Gaussian information space;
* the final quantitative test is a full H18/A21 matrix comparison against the
  endpoint covariance upper.

Until B_W and J_W are rigorously enclosed and the generalized-eigenvalue test is
emitted, this gate fails closed.  It must never fall back to the retired
one-sample tube margin, determinant/trace scalarization, source-history graph,
or another recursive covariance-lower experiment.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_lifted_word_p3 as LIFTED

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_SEA3_LIFTED_FINITE_WORD_P3_GATE"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    lifted = LIFTED.build(path, tube_path)
    failures = LIFTED.validate(lifted)
    if failures:
        raise RuntimeError(f"lifted P3 architecture failed validation: {failures}")

    modes = {}
    fail_reasons = []
    quantitative_closed = bool(lifted["P3_QUANTITATIVE_WORD_MATRIX_CLOSED"])
    for mode in ("H", "A"):
        src = lifted["modes"][mode]
        delta = src["relative_Riccati_injection_margin_lower"]
        mode_closed = all(bool(src[k]) for k in (
            "lifted_endpoint_map_B_closed",
            "lifted_measurement_information_J_upper_closed",
            "selected_mode_posterior_lower_closed",
            "Omega_selected_full_matrix_lower_closed",
            "generalized_eigenvalue_comparison_closed",
        ))
        if delta is None:
            delta_value = 0.0
            mode_pass = False
        else:
            delta_value = float(delta)
            if not (math.isfinite(delta_value) and delta_value > 0.0):
                raise RuntimeError(f"{mode} emitted invalid lifted margin {delta!r}")
            mode_pass = mode_closed and delta_value >= USEFUL_GATE
        modes[mode] = {
            **src,
            "relative_Riccati_injection_margin_lower": delta_value,
            "contraction_gap_lower": delta_value,
            "useful_margin_gate": USEFUL_GATE,
            "pass": mode_pass,
        }
        if not mode_closed:
            fail_reasons.append(
                f"{mode} lifted word matrices B_W/J_W and Omega_sel generalized comparison are not yet closed"
            )
        elif not mode_pass:
            fail_reasons.append(
                f"{mode} lifted finite-word margin {delta_value:.17g} is below useful gate {USEFUL_GATE:.17g}"
            )

    canonical_pass = quantitative_closed and all(modes[m]["pass"] for m in ("H", "A"))
    if canonical_pass != bool(lifted["P3_CANONICAL_PASS"]):
        raise RuntimeError("lifted architecture aggregate verdict disagrees with canonical gate")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "LIFTED_FINITE_WORD_SELECTED_PROCESS_MODES",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_consumed": True,
        "lifted_word_architecture_consumed": True,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "old_P2_800_state_graph_consumed": False,
        "old_P2_V1_history_frontier_consumed": False,
        "old_terminal_source_phase_metric_attachment_consumed": False,
        "one_sample_strict_Riccati_margin_consumed": False,
        "commit_aligned_source_word_consumed": False,
        "per_sample_SPD_lower_required": False,
        "determinant_trace_scalarization_used": False,
        "scalar_information_beta_used": False,
        "parameter_dependent_metric": "V_k=e_k^T P_k^-1 e_k with P_k the shipping Riccati covariance",
        "samplewise_nonexpansion_closed": True,
        "strictness_location": "RECURRENT_FINITE_WORD_ONLY",
        "lifted_word_identity": lifted["lifted_word_identity"],
        "source_uniform_matrix_enclosure_contract": lifted["source_uniform_matrix_enclosure_contract"],
        "retained_strict_subcertificates": lifted["retained_strict_subcertificates"],
        "word_horizon_s": lifted["word_horizon_s"],
        "modes": modes,
        "useful_gate": USEFUL_GATE,
        "P3_FOUNDATION_PASS": True,
        "P3_ARCHITECTURE_READY": True,
        "P3_QUANTITATIVE_WORD_MATRIX_CLOSED": quantitative_closed,
        "P3_CANONICAL_PASS": canonical_pass,
        "P4_MAY_CONSUME_P3": canonical_pass,
        "P3_CANONICAL_FAIL_REASONS": list(dict.fromkeys(fail_reasons)),
        "next_obligation": (
            "P3 is quantitatively closed; proceed to P4"
            if canonical_pass
            else "enclose only B_W and J_W for the fixed lifted word architecture, form Omega_sel, and run the full-matrix H18/A21 gate; do not introduce another P3 architecture"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "LIFTED_FINITE_WORD_SELECTED_PROCESS_MODES":
        f.append("wrong canonical P3 architecture")
    for key in (
        "source_generated_not_trajectory_fit", "SEA3_dynamic_source_consumed",
        "lifted_word_architecture_consumed", "samplewise_nonexpansion_closed",
        "P3_FOUNDATION_PASS", "P3_ARCHITECTURE_READY",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed", "old_P2_V1_history_frontier_consumed",
        "old_terminal_source_phase_metric_attachment_consumed",
        "one_sample_strict_Riccati_margin_consumed", "commit_aligned_source_word_consumed",
        "per_sample_SPD_lower_required", "determinant_trace_scalarization_used",
        "scalar_information_beta_used",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("strictness_location") != "RECURRENT_FINITE_WORD_ONLY":
        f.append("strictness moved away from recurrent finite word")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")

    expected = bool(d.get("P3_QUANTITATIVE_WORD_MATRIX_CLOSED"))
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        delta = row.get("relative_Riccati_injection_margin_lower")
        if not isinstance(delta, (int, float)) or not math.isfinite(float(delta)) or float(delta) < 0.0:
            f.append(f"{mode} lifted margin is invalid")
            expected = False
            continue
        mode_closed = all(bool(row.get(k)) for k in (
            "lifted_endpoint_map_B_closed",
            "lifted_measurement_information_J_upper_closed",
            "selected_mode_posterior_lower_closed",
            "Omega_selected_full_matrix_lower_closed",
            "generalized_eigenvalue_comparison_closed",
        ))
        mode_expected = mode_closed and float(delta) >= USEFUL_GATE
        if row.get("pass") is not mode_expected:
            f.append(f"{mode} verdict disagrees with lifted matrix closure/margin")
        expected = expected and mode_expected

    if d.get("P3_CANONICAL_PASS") is not expected:
        f.append("canonical P3 verdict disagrees with lifted H/A result")
    if d.get("P4_MAY_CONSUME_P3") is not expected:
        f.append("P4 handoff disagrees with canonical P3 verdict")
    reasons = d.get("P3_CANONICAL_FAIL_REASONS", [])
    if expected and reasons:
        f.append("passing P3 still reports failure reasons")
    if not expected and not reasons:
        f.append("failing P3 does not report the open lifted-matrix obligation")
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
        "P3_ARCHITECTURE_READY": d["P3_ARCHITECTURE_READY"],
        "P3_QUANTITATIVE_WORD_MATRIX_CLOSED": d["P3_QUANTITATIVE_WORD_MATRIX_CLOSED"],
        "H_delta": d["modes"]["H"]["relative_Riccati_injection_margin_lower"],
        "A_delta": d["modes"]["A"]["relative_Riccati_injection_margin_lower"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
