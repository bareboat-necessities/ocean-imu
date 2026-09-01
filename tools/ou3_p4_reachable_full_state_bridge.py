#!/usr/bin/env python3
"""Bridge validated translation/full-word results into the reachable full-state P4 test.

This producer combines source-only objects on the same recurrent worst-cell route:

* the outward-validated one-second complete-word translation margin;
* the direct nontranslation generalized margin diagnostic;
* the source-dynamic reachability graph; and, optionally,
* an outward-validated complete-word translation/nontranslation cross-block
  certificate covering every required recurrent source path/cell.

For a normalized full-state endpoint comparison split into translation and
nontranslation blocks,

    G = [[A, C], [C^T, B]],

with A >= a I and B >= b I, a sufficient Schur/Young condition for G >> 0 is

    ||C||_2 < sqrt(a b).

The bridge emits that budget and SVD-free sufficient variants.  When a
cross-block certificate is supplied, it is accepted only if it is source-only,
outward validated, covers all required reachable/recurrent paths, has matching
18/21-state dimensions, and proves either

    ||C||_1 ||C||_inf < a b

or

    ||C||_F < sqrt(a b).

Without such an input the bridge remains a target-only object.  Even after the
full-state *linear* endpoint closes, the nonlinear finite-angle return-map lift
remains a separate P4 obligation, so this producer never promotes the final
usable P4 certificate by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


TRANSLATION_DIM = 12
NONTRANSLATION_DIM = {"H": 6, "A": 9}


def _positive_finite(x) -> float | None:
    try:
        q = float(x)
    except (TypeError, ValueError):
        return None
    return q if math.isfinite(q) and q > 0.0 else None


def _cross_mode_status(mode: str, cross: dict | None, product_budget: float,
                       frobenius_budget: float, reachable: int,
                       recurrent: int) -> dict:
    pending = {
        "provided": False,
        "outward_validated": False,
        "all_required_paths_checked": False,
        "one_inf_product_upper": None,
        "frobenius_norm_upper": None,
        "one_inf_test_pass": False,
        "frobenius_test_pass": False,
        "accepted": False,
        "reasons": ["no complete-word cross-block certificate supplied"],
    }
    if cross is None:
        return pending

    reasons: list[str] = []
    if cross.get("source_only") is not True or cross.get("trajectory_replay_used") is not False:
        reasons.append("cross-block certificate is not source-only")
    if cross.get("outward_validated") is not True:
        reasons.append("cross-block certificate is not outward validated")
    if cross.get("all_required_paths_checked") is not True:
        reasons.append("cross-block certificate does not cover all required paths")
    if int(cross.get("reachable_state_count", -1)) != reachable:
        reasons.append("cross-block reachable-state count does not match source graph")
    if int(cross.get("recurrent_state_count", -1)) != recurrent:
        reasons.append("cross-block recurrent-state count does not match source graph")

    m = cross.get("modes", {}).get(mode, {})
    if int(m.get("translation_block_dimension", -1)) != TRANSLATION_DIM:
        reasons.append(f"{mode}: wrong translation cross-block dimension")
    if int(m.get("nontranslation_block_dimension", -1)) != NONTRANSLATION_DIM[mode]:
        reasons.append(f"{mode}: wrong nontranslation cross-block dimension")
    if m.get("all_required_cells_checked") is not True:
        reasons.append(f"{mode}: not all required source cells checked")

    one = _positive_finite(m.get("normalized_cross_block_one_norm_upper"))
    inf = _positive_finite(m.get("normalized_cross_block_inf_norm_upper"))
    frob = _positive_finite(m.get("normalized_cross_block_frobenius_norm_upper"))
    product = None if one is None or inf is None else one * inf
    one_inf_pass = product is not None and math.isfinite(product) and product < product_budget
    frob_pass = frob is not None and frob < frobenius_budget
    if not (one_inf_pass or frob_pass):
        reasons.append(f"{mode}: validated cross-block norm does not fit the Schur budget")

    return {
        "provided": True,
        "outward_validated": cross.get("outward_validated") is True,
        "all_required_paths_checked": cross.get("all_required_paths_checked") is True,
        "one_norm_upper": one,
        "inf_norm_upper": inf,
        "one_inf_product_upper": product,
        "frobenius_norm_upper": frob,
        "one_inf_test_pass": one_inf_pass,
        "frobenius_test_pass": frob_pass,
        "accepted": not reasons,
        "reasons": reasons,
    }


def build(translation: dict, bottleneck: dict, path: dict,
          cross_block: dict | None = None) -> dict:
    failures: list[str] = []
    if translation.get("P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS") != "PASS":
        failures.append("validated complete-word translation input is not PASS")
    if path.get("path_graph_ready") is not True:
        failures.append("source-path graph is not ready")
    if int(path.get("recurrent_states", 0)) <= 0:
        failures.append("source-path graph has no recurrent states")

    reachable = int(path.get("partition", {}).get("states", 0))
    recurrent = int(path.get("recurrent_states", 0))
    modes: dict[str, dict] = {}
    all_linear = True
    for mode in ("H", "A"):
        tr = translation.get("modes", {}).get(mode, {})
        bn = bottleneck.get("modes", {}).get(mode, {})
        a = _positive_finite(tr.get("complete_word_translation_margin_lower"))
        b = _positive_finite(bn.get("existing_direct_nontranslation_margin_lower"))
        if a is None or b is None:
            failures.append(f"{mode}: missing positive block margins")
            all_linear = False
            continue
        budget = math.sqrt(a * b)
        product_budget = a * b
        n_non = NONTRANSLATION_DIM[mode]
        entry_budget = budget / math.sqrt(TRANSLATION_DIM * n_non)
        cross_status = _cross_mode_status(
            mode, cross_block, product_budget, budget, reachable, recurrent
        )
        linear_ok = cross_status["accepted"]
        all_linear = all_linear and linear_ok
        modes[mode] = {
            "validated_translation_margin_lower": a,
            "direct_nontranslation_margin_lower": b,
            "normalized_cross_block_spectral_norm_budget_upper_open": budget,
            "normalized_cross_block_one_inf_norm_product_budget_upper_open": product_budget,
            "normalized_cross_block_frobenius_norm_budget_upper_open": budget,
            "uniform_normalized_cross_block_entry_abs_budget_upper_open": entry_budget,
            "translation_block_dimension": TRANSLATION_DIM,
            "nontranslation_block_dimension": n_non,
            "full_state_dimension": TRANSLATION_DIM + n_non,
            "full_state_sufficient_condition": "validated ||C||_2 < sqrt(delta_translation*delta_nontranslation)",
            "svd_free_sufficient_condition": "validated ||C||_1*||C||_inf < delta_translation*delta_nontranslation",
            "uniform_entry_sufficient_condition": "every normalized |C_ij| < sqrt(delta_translation*delta_nontranslation)/sqrt(n_translation*n_nontranslation)",
            "cross_block": cross_status,
            "cross_block_bound_validated": linear_ok,
            "full_state_linear_certificate_established": linear_ok,
        }

    linear_status = "PASS" if all_linear and len(modes) == 2 else "PENDING"
    next_obligation = (
        "lift the validated reachable full-state linear endpoint through the exact finite-angle return map"
        if linear_status == "PASS" else
        "propagate the complete 18/21-state Phi/Omega endpoint on the reachable source graph and "
        "outward-validate the normalized translation/nontranslation cross block on every required cell/path"
    )
    return {
        "qualification": "OU3_P4_REACHABLE_FULL_STATE_CROSS_BLOCK_BRIDGE",
        "source_only": True,
        "trajectory_replay_used": False,
        "reachable_state_count": reachable,
        "recurrent_state_count": recurrent,
        "source_graph_strongly_connected_components": int(path.get("strongly_connected_components", 0)),
        "old_worst_corner_recurrent": bool(path.get("old_worst_corner_has_internal_recurrent_cycle", False)),
        "cross_block_certificate_supplied": cross_block is not None,
        "modes": modes,
        "P4_FULL_STATE_LINEAR_STATUS": linear_status,
        "P4_USABLE_CERTIFICATE_STATUS": "NOT_ESTABLISHED",
        "next_obligation": next_obligation,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("source_only") is not True or d.get("trajectory_replay_used") is not False:
        f.append("bridge is not source-only")
    if d.get("P4_USABLE_CERTIFICATE_STATUS") != "NOT_ESTABLISHED":
        f.append("bridge prematurely promoted P4")
    if int(d.get("recurrent_state_count", 0)) <= 0:
        f.append("bridge lost recurrent source graph")
    linear_status = d.get("P4_FULL_STATE_LINEAR_STATUS")
    if linear_status not in ("PASS", "PENDING"):
        f.append("invalid full-state linear status")
    accepted = []
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        q = m.get("normalized_cross_block_spectral_norm_budget_upper_open")
        p = m.get("normalized_cross_block_one_inf_norm_product_budget_upper_open")
        e = m.get("uniform_normalized_cross_block_entry_abs_budget_upper_open")
        if not isinstance(q, (int, float)) or not math.isfinite(float(q)) or float(q) <= 0.0:
            f.append(f"{mode}: invalid cross-block budget")
        if not isinstance(p, (int, float)) or not math.isfinite(float(p)) or float(p) <= 0.0:
            f.append(f"{mode}: invalid one/inf product budget")
        if not isinstance(e, (int, float)) or not math.isfinite(float(e)) or float(e) <= 0.0:
            f.append(f"{mode}: invalid uniform entry budget")
        if isinstance(q, (int, float)) and isinstance(p, (int, float)) and not math.isclose(float(p), float(q) * float(q), rel_tol=1e-14):
            f.append(f"{mode}: one/inf product budget is not square of spectral budget")
        if int(m.get("translation_block_dimension", 0)) != TRANSLATION_DIM:
            f.append(f"{mode}: wrong translation dimension")
        if int(m.get("nontranslation_block_dimension", 0)) != NONTRANSLATION_DIM[mode]:
            f.append(f"{mode}: wrong nontranslation dimension")
        if int(m.get("full_state_dimension", 0)) != TRANSLATION_DIM + NONTRANSLATION_DIM[mode]:
            f.append(f"{mode}: wrong full-state dimension")
        c = m.get("cross_block", {})
        ok = m.get("cross_block_bound_validated") is True
        accepted.append(ok)
        if ok != (c.get("accepted") is True):
            f.append(f"{mode}: inconsistent cross-block acceptance")
        if (m.get("full_state_linear_certificate_established") is True) != ok:
            f.append(f"{mode}: inconsistent full-state linear status")
    if linear_status == "PASS" and not all(accepted):
        f.append("full-state linear PASS without both validated mode cross blocks")
    if linear_status == "PENDING" and all(accepted):
        f.append("full-state linear status stayed pending after both modes passed")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--translation", type=Path, required=True)
    ap.add_argument("--bottleneck", type=Path, required=True)
    ap.add_argument("--path", type=Path, required=True)
    ap.add_argument("--cross-block", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(
        json.loads(a.translation.read_text(encoding="utf-8")),
        json.loads(a.bottleneck.read_text(encoding="utf-8")),
        json.loads(a.path.read_text(encoding="utf-8")),
        json.loads(a.cross_block.read_text(encoding="utf-8")) if a.cross_block else None,
    )
    f = validate(d)
    d["validation_failures"] = f
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_USABLE_CERTIFICATE_STATUS"],
        "linear_status": d["P4_FULL_STATE_LINEAR_STATUS"],
        "reachable_states": d["reachable_state_count"],
        "recurrent_states": d["recurrent_state_count"],
        "old_worst_corner_recurrent": d["old_worst_corner_recurrent"],
        "modes": {
            mode: {
                "translation": d.get("modes", {}).get(mode, {}).get("validated_translation_margin_lower"),
                "nontranslation": d.get("modes", {}).get(mode, {}).get("direct_nontranslation_margin_lower"),
                "cross_block_budget": d.get("modes", {}).get(mode, {}).get("normalized_cross_block_spectral_norm_budget_upper_open"),
                "one_inf_product_budget": d.get("modes", {}).get(mode, {}).get("normalized_cross_block_one_inf_norm_product_budget_upper_open"),
                "uniform_entry_budget": d.get("modes", {}).get(mode, {}).get("uniform_normalized_cross_block_entry_abs_budget_upper_open"),
                "cross_block_accepted": d.get("modes", {}).get(mode, {}).get("cross_block_bound_validated"),
            }
            for mode in ("H", "A")
        },
        "failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
