#!/usr/bin/env python3
"""Canonical P3 gate: quantitative SEA3 moving-Riccati certificate.

The canonical linear theorem uses the actual shipping Riccati covariance as the
parameter-dependent Lyapunov metric

    V_k(e) = e^T P_k^{-1} e.

For one recurrent linearized word,

    e_+ = Phi e,
    P_+ = Phi P Phi^T + Omega.

The quantitative SEA3 Riccati-tube producer certifies a source-uniform upper
bound on P_+ and a source-uniform lower comparison

    Omega >= delta P_+.

Therefore

    Phi P Phi^T <= (1-delta) P_+

and hence

    V_+ <= (1-delta) V.

This file is the *gate*, not a second proof route.  It consumes the quantitative
result from the canonical SEA3 Riccati tube (with the algebraically factored
small-x numerical backend) and promotes P3 if and only if both H18 and A21
margins meet the unchanged 1e-18 useful-margin gate.  No 800-state P2 graph,
P2-V1 history frontier, terminal source/phase metric attachment, or source-
history enumeration may substitute for that verdict.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_riccati_tube_factored as TUBE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_MOVING_RICCATI_METRIC_P3"
USEFUL_GATE = 1.0e-18


def _positive(x, label: str) -> float:
    y = float(x)
    if not (math.isfinite(y) and y > 0.0):
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return y


def _load_tube(domain_path: Path, tube_path: Path | None) -> dict:
    if tube_path is None:
        tube = TUBE.build(domain_path)
    else:
        tube = json.loads(Path(tube_path).read_text(encoding="utf-8"))
    failures = TUBE.validate(tube)
    if failures:
        raise RuntimeError(f"quantitative moving-Riccati tube failed validation: {failures}")
    if float(tube.get("useful_gate", math.nan)) != USEFUL_GATE:
        raise RuntimeError("quantitative tube changed the canonical 1e-18 useful-margin gate")
    return tube


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    tube_path: Path | None = None,
) -> dict:
    path = Path(domain_path).resolve()

    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    if df:
        raise RuntimeError(f"SEA3 dynamic source prerequisite failed: {df}")

    tube = _load_tube(path, tube_path)

    modes = {}
    for mode, dim in (("H", 18), ("A", 21)):
        trow = tube["modes"][mode]
        delta = _positive(
            trow["relative_Riccati_injection_margin_lower"],
            f"{mode} relative Riccati injection margin",
        )
        margin_pass = delta >= USEFUL_GATE
        if bool(trow.get("useful_margin_pass")) != margin_pass:
            raise RuntimeError(
                f"{mode} tube useful-margin verdict is inconsistent with delta={delta!r}"
            )

        # Do not rely on a binary64 representation of 1-delta for strictness:
        # at the 1e-18 gate, 1-delta rounds to 1.0.  The certified contraction
        # is represented by its positive gap delta and the exact formula.
        modes[mode] = {
            "dimension": dim,
            "riccati_covariance_upper_bound_closed": True,
            "word_injection_comparison_closed": True,
            "Pbar_lambda_max_trace_upper": _positive(
                trow["Pbar_lambda_max_trace_upper"], f"{mode} Pbar trace upper"
            ),
            "relative_Riccati_injection_margin_lower": delta,
            "contraction_gap_lower": delta,
            "contraction_factor_upper_formula": "1-relative_Riccati_injection_margin_lower",
            "useful_margin_gate": USEFUL_GATE,
            "pass": margin_pass,
            "worst_current_source_cell": trow["worst_current_source_cell"],
        }

    canonical_pass = all(modes[m]["pass"] for m in ("H", "A"))
    tube_pass = bool(tube["RICCATI_TUBE_PASS"])
    if tube_pass != canonical_pass:
        raise RuntimeError(
            "quantitative tube aggregate verdict disagrees with H/A canonical margin checks"
        )

    fail_reasons = []
    if not canonical_pass:
        for mode in ("H", "A"):
            if not modes[mode]["pass"]:
                fail_reasons.append(
                    f"{mode} relative Riccati injection margin "
                    f"{modes[mode]['relative_Riccati_injection_margin_lower']:.17g} "
                    f"is below useful gate {USEFUL_GATE:.17g}"
                )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "MOVING_SHIPPING_RICCATI_COVARIANCE_METRIC",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_consumed": True,
        "quantitative_Riccati_tube_consumed": True,
        "quantitative_Riccati_tube_pass": tube_pass,
        "current_source_interval_cover_only": bool(tube["current_source_interval_cover_only"]),
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
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
            "vector_PE_window_s": dynamic["normal_live_contract"][
                "vector_PE_recurrence_window_s"
            ],
            "accelerometer_each_valid_live_sample": dynamic["normal_live_contract"][
                "accelerometer_update_required_each_valid_sample"
            ],
            "adaptive_state_rate_bounded": True,
            "adaptive_state": dynamic["adaptive_state"],
        },
        "cell_cover": tube["cell_cover"],
        "modes": modes,
        "useful_gate": USEFUL_GATE,
        "P3_FOUNDATION_PASS": True,
        "P3_CANONICAL_PASS": canonical_pass,
        "P4_MAY_CONSUME_P3": canonical_pass,
        "P3_CANONICAL_FAIL_REASONS": fail_reasons,
        "next_obligation": (
            "P3 is quantitatively closed; proceed to nonlinear moving-Riccati P4 on the full 0.8-rad sector"
            if canonical_pass
            else "tighten only the SEA3 dynamic-state/finite-memory Riccati bounds that limit the reported H/A margin; do not reintroduce source histories"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "MOVING_SHIPPING_RICCATI_COVARIANCE_METRIC":
        f.append("wrong canonical P3 architecture")

    for key in (
        "source_generated_not_trajectory_fit",
        "SEA3_dynamic_source_consumed",
        "quantitative_Riccati_tube_consumed",
        "current_source_interval_cover_only",
        "metric_change_handled_by_exact_Riccati_recursion",
        "P3_FOUNDATION_PASS",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")

    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_shrunk",
        "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed",
        "old_P2_V1_history_frontier_consumed",
        "old_terminal_source_phase_metric_attachment_consumed",
        "metric_derivative_or_jump_penalty_required",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")
    if int(d.get("cell_cover", {}).get("history_depth", -1)) != 0:
        f.append("canonical P3 reintroduced source-history depth")

    expected_pass = True
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        delta = row.get("relative_Riccati_injection_margin_lower")
        if not isinstance(delta, (int, float)) or not math.isfinite(float(delta)) or float(delta) <= 0.0:
            f.append(f"{mode} Riccati margin is not finite positive")
            expected_pass = False
            continue
        mode_expected = float(delta) >= USEFUL_GATE
        expected_pass = expected_pass and mode_expected
        if row.get("riccati_covariance_upper_bound_closed") is not True:
            f.append(f"{mode} covariance upper bound not closed")
        if row.get("word_injection_comparison_closed") is not True:
            f.append(f"{mode} word injection comparison not closed")
        if row.get("pass") is not mode_expected:
            f.append(f"{mode} P3 verdict disagrees with quantitative margin")
        pbar = row.get("Pbar_lambda_max_trace_upper")
        if not isinstance(pbar, (int, float)) or not math.isfinite(float(pbar)) or float(pbar) <= 0.0:
            f.append(f"{mode} Pbar trace upper is not finite positive")

    if d.get("quantitative_Riccati_tube_pass") is not expected_pass:
        f.append("tube aggregate verdict disagrees with H/A margins")
    if d.get("P3_CANONICAL_PASS") is not expected_pass:
        f.append("canonical P3 verdict is not the quantitative H/A verdict")
    if d.get("P4_MAY_CONSUME_P3") is not expected_pass:
        f.append("P4 handoff does not follow canonical P3 verdict")

    reasons = d.get("P3_CANONICAL_FAIL_REASONS", [])
    if expected_pass and reasons:
        f.append("passing P3 still reports failure reasons")
    if not expected_pass and not reasons:
        f.append("failing P3 does not name limiting margin")

    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument(
        "--tube",
        type=Path,
        default=None,
        help="consume an already-built canonical SEA3 Riccati-tube JSON artifact",
    )
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.tube)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P3_architecture"],
        "tube_pass": d["quantitative_Riccati_tube_pass"],
        "H_delta": d["modes"]["H"]["relative_Riccati_injection_margin_lower"],
        "A_delta": d["modes"]["A"]["relative_Riccati_injection_margin_lower"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "P4_MAY_CONSUME_P3": d["P4_MAY_CONSUME_P3"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
