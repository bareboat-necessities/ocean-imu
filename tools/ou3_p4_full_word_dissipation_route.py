#!/usr/bin/env python3
"""P4 replacement route: quantify the whole word, not one covariance seed.

The current P3/P4 scalar frontier is numerically dominated by a deliberate
conservatism in P3: a strict generalized margin is established from one
source-reachable covariance injection and then merely *preserved* through the
remaining affine-PSD covariance operations.  Every later process Q and Joseph
K R K^T term is known to improve the exact covariance-decomposition margin,
but the present numerical delta does not quantify that improvement.

For a complete fixed-mode word,

    P_N = Phi_N P_0 Phi_N^T + Omega_N,

with

    Omega_{k+1} = A_k Omega_k A_k^T + B_k,   B_k >= 0.

The useful linear contraction is therefore determined by the *complete*
Omega_N (equivalently by the complete deterministic information dissipation),
not by the first positive seed retained by the old proof.

This producer is an executable theorem-route gate.  It binds the replacement
route to the shipping word algebra and refuses to call the existing tiny P4
level usable.  Numerical promotion is intentionally blocked until a validated
source-cell backend propagates Phi/Omega over the complete recurrent-PE word
and certifies the resulting path-dependent matrix inequality directly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_explicit_information_word_certificate as P3
import ou3_p3_word_algebra as ALG
import ou3_p4_frontier_combined_certificate as OLD
import ou3_p4_exact_word_map as WORD

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    p = Path(domain_path).resolve()
    p3 = P3.build(p)
    alg = ALG.build()
    word = WORD.build(p)
    old = OLD.build(p)
    failures = [f"P3: {x}" for x in P3.validate(p3)]
    failures += [f"ALG: {x}" for x in ALG.validate(alg)]
    failures += [f"WORD: {x}" for x in WORD.validate(word)]
    failures += [f"OLD: {x}" for x in OLD.validate(old)]

    modes = {}
    for mode in ("H", "A"):
        row = p3["modes"][mode]
        arg = row["matrix_comparison"]["word_endpoint_information_argument"]
        oldm = old["modes"][mode]
        seed_only = (
            arg.get("endpoint_noise_lower_source") ==
            "validated matrix lower contribution from a source-reachable prediction/correction stage"
        )
        if not seed_only:
            failures.append(f"{mode}: P3 endpoint-noise provenance changed; re-audit full-word route")
        modes[mode] = {
            "current_P3_delta_lower": row["word_endpoint_relative_Riccati_injection_margin_lower"],
            "current_P4_W": oldm["certified_level_W"],
            "current_P4_prefix_norm": oldm["prefix_canonical_error_norm_upper"],
            "current_numeric_margin_is_single_seed_then_preserved": seed_only,
            "complete_word_covariance_identity": alg["covariance_decomposition_invariant"]["identity"],
            "strict_margin_preservation_identity": alg["strict_margin_preservation"]["identity"],
            "later_additive_PSD_terms_exist": True,
            "later_additive_PSD_terms_quantified_in_current_delta": False,
            "replacement_linear_test": (
                "validate Omega_word - delta_cell*P_word >> 0 on each reachable source/path cell, "
                "or equivalently validate Phi_word^T M_end Phi_word <= rho_cell M_start with rho_cell<1"
            ),
            "replacement_nonlinear_test": (
                "validated exact finite-angle return-map enclosure on the same reachable source/path cells; "
                "no scalar B*W endpoint reduction"
            ),
            "usable_certificate_status": "NOT_ESTABLISHED",
        }

    return {
        "qualification": "OU3_P4_FULL_WORD_DISSIPATION_REPLACEMENT_ROUTE",
        "source_only": True,
        "trajectory_replay_used_for_promotion": False,
        "old_scalar_frontier_is_regression_baseline_only": True,
        "word_horizon_s": word["source_word_horizon_s"],
        "word_samples_upper": word["word_samples_upper"],
        "shipping_operation_order_bound": True,
        "full_word_additive_covariance_terms_to_accumulate": [
            "every prediction Q",
            "every accepted S-zero Joseph K R K^T",
            "every accepted accelerometer Joseph K R K^T",
            "every accepted magnetometer Joseph K R K^T",
            "every PSD a_w covariance synchronization increment",
        ],
        "path_metric_required": True,
        "source_reachability_required": True,
        "scalar_uniform_min_delta_for_all_cells_forbidden_as_final_route": True,
        "scalar_BW_small_gain_forbidden_as_final_route": True,
        "next_numeric_backend": (
            "validated source-cell complete-word Phi/Omega propagation with outward-rounded matrices; "
            "use numerical LMIs only to design candidate path metrics, then independently validate every cell"
        ),
        "modes": modes,
        "P4_USABLE_CERTIFICATE_STATUS": "NOT_ESTABLISHED",
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("source_only") is not True:
        f.append("replacement route is not source-only")
    if d.get("old_scalar_frontier_is_regression_baseline_only") is not True:
        f.append("old tiny scalar frontier was not demoted to baseline")
    if d.get("path_metric_required") is not True or d.get("source_reachability_required") is not True:
        f.append("replacement route lost path/source correlation")
    if d.get("scalar_uniform_min_delta_for_all_cells_forbidden_as_final_route") is not True:
        f.append("replacement route permits global worst-cell scalarization")
    if d.get("scalar_BW_small_gain_forbidden_as_final_route") is not True:
        f.append("replacement route permits scalar B*W final gate")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("current_numeric_margin_is_single_seed_then_preserved") is not True:
            f.append(f"{mode}: current seed-only conservatism not established")
        if m.get("later_additive_PSD_terms_quantified_in_current_delta") is not False:
            f.append(f"{mode}: route incorrectly says current delta accumulates the whole word")
        if m.get("usable_certificate_status") != "NOT_ESTABLISHED":
            f.append(f"{mode}: route prematurely promoted a usable certificate")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_failures"] = f
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True))
    print(json.dumps({
        "status": d["P4_USABLE_CERTIFICATE_STATUS"],
        "route": d["qualification"],
        "modes": {
            mode: {
                "old_delta": d["modes"][mode]["current_P3_delta_lower"],
                "old_W": d["modes"][mode]["current_P4_W"],
                "seed_only": d["modes"][mode]["current_numeric_margin_is_single_seed_then_preserved"],
            }
            for mode in ("H", "A")
        },
        "failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
