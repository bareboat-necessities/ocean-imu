#!/usr/bin/env python3
"""Canonical P3 architecture: moving Riccati-covariance Lyapunov metric.

The old canonical route enumerated long histories of a broad 800-state tuner
abstraction and then attached a worst-case covariance envelope to each terminal
source/phase class.  That proves a stronger switching theorem than the SEA3
deployment theorem and manufactures the wall we are trying to remove.

For the exact linearized KF error recursion over any recurrent word,

    e_+ = Phi e,
    P_+ = Phi P Phi' + Omega,

where P is the *shipping Riccati covariance* and Omega is the covariance injected
through process noise and the actual Joseph measurement operations.  With

    V(e,P) = e' P^{-1} e

and a certified comparison Omega >= delta P_+, 0 < delta <= 1,

    Phi P Phi' <= (1-delta) P_+

implies directly

    V_+ <= (1-delta) V.

P may vary with every SEA3-driven tuner update: it is the metric state, not a
fixed matrix that must survive arbitrary parameter jumps.  Therefore there is
no separate dM/dt penalty and no source-word history enumeration.  Uniform UCO,
UCC and the rate-bounded SEA3 adaptive-state tube are the ingredients needed to
certify a positive source-uniform delta.

This producer makes that architecture canonical and fail-closed.  Until the
validated Riccati tube P <= P_bar and word injection Omega >= Omega_bar are
closed numerically for H18 and A21, P3 remains OPEN; no legacy P2 history result
may substitute for it.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR
import ou3_full_process_ucc as PROCESS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_MOVING_RICCATI_METRIC_P3"
USEFUL_GATE = 1.0e-18


def _positive(x, label: str) -> float:
    y = float(x)
    if not (math.isfinite(y) and y > 0.0):
        raise RuntimeError(f"{label} must be finite positive")
    return y


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    trans = TRANS.build(TRANS.DEFAULT_HEADER)
    tf = TRANS.validate(trans)
    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    process = PROCESS.build()
    pf = PROCESS.validate(process)
    failures = [f"dynamic: {x}" for x in df]
    failures += [f"translation: {x}" for x in tf]
    failures += [f"vector: {x}" for x in vf]
    failures += [f"process: {x}" for x in pf]
    if failures:
        raise RuntimeError(f"moving-Riccati prerequisites failed: {failures}")

    qH = _positive(process["modes"]["H"]["prediction_Q_lambda_min_lower"], "H process UCC")
    qA = _positive(process["modes"]["A"]["prediction_Q_lambda_min_lower"], "A process UCC")
    alpha6 = _positive(vector["gyro_bias_two_packet"]["alpha_6_information_lower"], "vector UCO")
    s_info = _positive(trans["integrator_detectability"]["information_gramian_lambda_min_lower"], "translation detectability")
    aw_alpha = _positive(trans["integrator_detectability"]["stable_aw_alpha_upper"], "a_w stability")
    pe_window = _positive(dynamic["normal_live_contract"]["vector_PE_recurrence_window_s"], "PE window")

    # The architecture proof is exact, but the quantitative delta still needs a
    # validated full-state Riccati tube for the time-varying H/A systems.  Keep
    # the stage open rather than recycling a terminal-history covariance hull.
    modes = {}
    for mode, dim, q in (("H", 18, qH), ("A", 21, qA)):
        modes[mode] = {
            "dimension": dim,
            "prediction_process_UCC_lower": q,
            "vector_information_UCO_lower": alpha6,
            "translation_information_detectability_lower": s_info,
            "stable_aw_alpha_upper": aw_alpha,
            "riccati_covariance_upper_bound_closed": False,
            "word_injection_comparison_closed": False,
            "relative_Riccati_injection_margin_lower": None,
            "contraction_factor_upper": None,
            "useful_margin_gate": USEFUL_GATE,
            "pass": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "MOVING_SHIPPING_RICCATI_COVARIANCE_METRIC",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_consumed": True,
        "old_P2_800_state_graph_consumed": False,
        "old_P2_V1_history_frontier_consumed": False,
        "old_terminal_source_phase_metric_attachment_consumed": False,
        "parameter_dependent_metric": "V_k = e_k^T P_k^{-1} e_k, P_k is the shipping Riccati covariance",
        "metric_derivative_or_jump_penalty_required": False,
        "metric_change_handled_by_exact_Riccati_recursion": True,
        "linear_word_identity": {
            "error": "e_plus = Phi e",
            "covariance": "P_plus = Phi P Phi^T + Omega",
            "sufficient_comparison": "Omega >= delta P_plus",
            "consequence": "Phi P Phi^T <= (1-delta) P_plus",
            "Lyapunov_contraction": "V_plus <= (1-delta) V",
        },
        "recurrent_word_contract": {
            "vector_PE_window_s": pe_window,
            "accelerometer_each_valid_live_sample": dynamic["normal_live_contract"][
                "accelerometer_update_required_each_valid_sample"
            ],
            "adaptive_state_rate_bounded": True,
            "adaptive_state": dynamic["adaptive_state"],
        },
        "prerequisite_constants": {
            "vector_alpha6_information_lower": alpha6,
            "translation_detectability_information_lower": s_info,
            "stable_aw_alpha_upper": aw_alpha,
            "H_prediction_Q_lambda_min_lower": qH,
            "A_prediction_Q_lambda_min_lower": qA,
        },
        "modes": modes,
        "useful_gate": USEFUL_GATE,
        "P3_FOUNDATION_PASS": True,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "P3_CANONICAL_FAIL_REASONS": [
            "validated time-varying H18/A21 Riccati covariance tube P <= P_bar not yet emitted",
            "source-uniform word comparison Omega >= delta P_plus not yet emitted",
        ],
        "next_obligation": (
            "bound the shipping time-varying Riccati covariance directly over the compact SEA3 dynamic source tube, "
            "then certify delta_H and delta_A from Omega >= delta P_plus; do not enumerate source histories"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "MOVING_SHIPPING_RICCATI_COVARIANCE_METRIC":
        f.append("wrong canonical P3 architecture")
    for key in (
        "source_generated_not_trajectory_fit", "SEA3_dynamic_source_consumed",
        "metric_change_handled_by_exact_Riccati_recursion", "P3_FOUNDATION_PASS",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "old_P2_800_state_graph_consumed", "old_P2_V1_history_frontier_consumed",
        "old_terminal_source_phase_metric_attachment_consumed",
        "metric_derivative_or_jump_penalty_required", "P3_CANONICAL_PASS",
        "P4_MAY_CONSUME_P3",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")
    if not d.get("P3_CANONICAL_FAIL_REASONS"):
        f.append("open P3 route does not name remaining obligations")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        if row.get("pass") is not False:
            f.append(f"{mode} falsely promoted")
        for key in (
            "prediction_process_UCC_lower", "vector_information_UCO_lower",
            "translation_information_detectability_lower", "stable_aw_alpha_upper",
        ):
            x = row.get(key)
            if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
                f.append(f"{mode}.{key} is not finite positive")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P3_architecture"],
        "foundation_pass": d["P3_FOUNDATION_PASS"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
