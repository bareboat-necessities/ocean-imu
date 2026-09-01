#!/usr/bin/env python3
"""Combined P4 frontier certificate with a self-consistent prefix bootstrap.

This producer keeps the successful theorem-preserving P4 refinements already
validated on the predecessor branches:

* operation-specific nonlinear defect constants;
* directional measurement operators and class-local Kalman-gain bounds;
* exact S=0 / effective-vector correction structure from P5;
* proof-radius continuation; and
* the direct strict-return condition

      sqrt(1-delta) + B sqrt(W) < 1.

The remaining predecessor bound still inherited the legacy bootstrap
W_prefix <= 4 W0.  That factor was convenient for closing the old proof, but it
need not be frozen at four once the strict endpoint gap is solved explicitly.

Write b0 for the word-defect coefficient before the bootstrap factor.  If
W_prefix <= gamma W0, then

      ||r_word||_M <= gamma b0 W0.

At the strict endpoint boundary choose x=sqrt(W0) so that

      gamma b0 x < g,   g = 1-sqrt(1-delta).

Every prefix then satisfies

      sqrt(W_prefix)/x <= 1 + gamma b0 x < 1+g.

Thus the bootstrap closes self-consistently with

      gamma >= (1+g)^2,

rather than the fixed gamma=4.  The implementation evaluates this with
outward-safe arithmetic and independently rechecks the design-radius, Cayley,
quaternion, projection, and strict-endpoint constraints.  It fails closed if
it does not strictly improve the predecessor direct certificate.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as L
import ou3_p4_thirdgen_combined_certificate as T
import ou3_p4_exact_correction_structure_certificate as E
import ou3_p4_direct_word_contraction_certificate as DIRECT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def _bootstrap_gamma_upper(gap_lower: float) -> float:
    one_plus = L.add_up(1.0, float(gap_lower))
    return L.mul_up(one_plus, one_plus)


def _candidate(mode: str, base: dict, domain: dict, q: float) -> dict:
    # E._candidate contains all successful exact-correction/directional/gain
    # defect constants at this proof radius.  Its B still contains the legacy
    # PREFIX_BOOTSTRAP_W_FACTOR=4, so divide that factor out first.
    erow = E._candidate(mode, base, domain, q)
    B4 = float(erow["B"])
    b0 = L.div_up(B4, float(L.PREFIX_BOOTSTRAP_W_FACTOR))
    gap = DIRECT.strict_gap(float(base["P3_word_endpoint_delta_lower"]))
    gamma = _bootstrap_gamma_upper(gap)
    if not gamma < float(L.PREFIX_BOOTSTRAP_W_FACTOR):
        raise RuntimeError("self-consistent bootstrap did not improve fixed factor four")
    B = L.mul_up(gamma, b0)

    sm = L.GROUP.sqrt_point(float(base["metric_lambda_min_lower"])).lo
    caps = {
        "strict_endpoint": L.mul_down(0.999999999999, L.div_down(gap, B)),
        "self_consistent_bootstrap": L.div_down(gap, B),
        "design_radius": T._positive_root(B, L.mul_down(q, sm)),
        "cayley_chart": T._positive_root(
            B, L.mul_down(float(base["cayley_norm_limit"]), sm)
        ),
    }
    proj = base.get("active_bias_projection")
    if mode == "A" and proj:
        caps["bias_projection"] = T._positive_root(
            B,
            L.mul_down(float(proj["interior_margin_lower_mps2"]), sm),
        )

    sw = min(caps.values())
    W = L.mul_down(sw, sw)
    defect_ratio = L.mul_up(B, sw)
    prefix_factor = L.add_up(1.0, defect_ratio)
    qprefix = L.div_up(L.mul_up(prefix_factor, sw), sm)

    if not defect_ratio < gap:
        raise RuntimeError("strict endpoint gap does not close")
    if not L.mul_up(prefix_factor, prefix_factor) <= gamma:
        raise RuntimeError("self-consistent prefix bootstrap does not close")
    if not qprefix <= q:
        raise RuntimeError("prefix leaves selected design radius")
    if not qprefix < float(base["cayley_norm_limit"]):
        raise RuntimeError("prefix leaves Cayley chart")

    return {
        "q": q,
        "B_fixed4": B4,
        "b0_without_bootstrap": b0,
        "bootstrap_gamma_upper": gamma,
        "B": B,
        "sqrtW": sw,
        "W": W,
        "strict_gap": gap,
        "defect_ratio_upper": defect_ratio,
        "prefix_factor_upper": prefix_factor,
        "qprefix": qprefix,
        "active_cap": min(caps, key=caps.get),
        "caps": caps,
    }


def _refine(mode: str, base: dict, domain: dict) -> dict:
    q0 = float(base["thirdgen_selected_design_norm"])
    rows = []
    for q in T._grid(q0):
        try:
            rows.append(_candidate(mode, base, domain, q))
        except Exception:
            pass
    if not rows:
        raise RuntimeError("no self-consistent frontier cell certified")
    best = max(rows, key=lambda x: x["W"])
    before = float(base["certified_level_W"])
    if not best["W"] > before:
        raise RuntimeError("self-consistent bootstrap did not strictly widen direct P4")

    m = copy.deepcopy(base)
    m.update(
        {
            "frontier_self_consistent_bootstrap": True,
            "frontier_candidates": rows,
            "frontier_selected_design_norm": best["q"],
            "frontier_selected_active_cap": best["active_cap"],
            "frontier_bootstrap_gamma_upper": best["bootstrap_gamma_upper"],
            "frontier_B_before_fixed4": best["B_fixed4"],
            "frontier_B_upper": best["B"],
            "frontier_B_reduction_factor_lower": L.div_down(
                best["B_fixed4"], best["B"]
            ),
            "certified_level_W_before_frontier": before,
            "certified_level_W": best["W"],
            "certified_level_sqrt_W": best["sqrtW"],
            "prefix_canonical_error_norm_upper": best["qprefix"],
            "frontier_W_widening_factor_lower": L.div_down(best["W"], before),
            "frontier_claim": (
                "strict exact-word contraction W_end<W0 with a self-consistent "
                "prefix bootstrap; fixed delta/2 decrease is not claimed at the enlarged level"
            ),
        }
    )
    return m


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    p = Path(domain_path).resolve()
    domain = json.loads(p.read_text())
    base = DIRECT.build(p)
    failures = list(DIRECT.validate(base))
    modes = {}
    if not failures:
        for mode in ("H", "A"):
            try:
                modes[mode] = _refine(mode, base["modes"][mode], domain)
            except Exception as exc:
                failures.append(f"{mode}: {exc}")
    out = copy.deepcopy(base)
    out["modes"] = modes
    out["frontier_source_only"] = True
    out["P4_FRONTIER_COMBINED_CERTIFICATE"] = (
        "PASS" if not failures and len(modes) == 2 else "FAIL"
    )
    out["failures"] = failures
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if not m.get("frontier_self_consistent_bootstrap"):
            failures.append(f"{mode}: missing self-consistent bootstrap")
            continue
        if not float(m["certified_level_W"]) > float(
            m["certified_level_W_before_frontier"]
        ):
            failures.append(f"{mode}: frontier did not strictly widen direct P4")
        if not float(m["frontier_bootstrap_gamma_upper"]) < 4.0:
            failures.append(f"{mode}: bootstrap factor did not improve legacy factor four")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_failures"] = failures
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True))
    print(
        json.dumps(
            {
                "status": d["P4_FRONTIER_COMBINED_CERTIFICATE"],
                "modes": {
                    mode: {
                        "W_before": d.get("modes", {}).get(mode, {}).get(
                            "certified_level_W_before_frontier"
                        ),
                        "W_after": d.get("modes", {}).get(mode, {}).get(
                            "certified_level_W"
                        ),
                        "factor": d.get("modes", {}).get(mode, {}).get(
                            "frontier_W_widening_factor_lower"
                        ),
                        "B_factor": d.get("modes", {}).get(mode, {}).get(
                            "frontier_B_reduction_factor_lower"
                        ),
                        "gamma": d.get("modes", {}).get(mode, {}).get(
                            "frontier_bootstrap_gamma_upper"
                        ),
                        "cap": d.get("modes", {}).get(mode, {}).get(
                            "frontier_selected_active_cap"
                        ),
                    }
                    for mode in ("H", "A")
                },
                "failures": failures,
            },
            indent=2,
        )
    )
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
