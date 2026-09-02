#!/usr/bin/env python3
"""Compose the retained OU-III P3 metric and broad P4 sector primitives.

This module is intentionally a theorem-promotion gate rather than another
parallel proof backend.  It consumes only the retained source-uniform P3
certificate, the global Cayley-sector certificate, the exact nonlinear vector
remainder sector, and the source-complete timing decomposition introduced in
PR #460.

The promotion condition is deliberately strict: every H/A mode must provide a
positive P3 information-injection margin and the nonlinear sector budget must
fit strictly inside that margin over the declared 0.8-rad entrance sector.
Only after this gate is numerically closed should a P5 finite-capture recurrence
be promoted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_reachable_matrix_p3 as P3
import ou3_p3_scaled_process as SP
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_vector_remainder_sector as REMAINDER
import ou3_p4_source_word_timing as TIMING

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
TARGET_SECTOR_RAD = 0.8


def _finite_positive(x, label: str) -> float:
    y = float(x)
    if not math.isfinite(y) or y <= 0.0:
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return y


def _extract_mode_budget(mode: str, p3: dict, cayley: dict, remainder: dict) -> dict:
    m = p3["modes"][mode]
    delta = _finite_positive(m["relative_Riccati_injection_margin_lower"], f"{mode} P3 delta")
    sigma_min = _finite_positive(m["Sigma_lambda_min_lower"], f"{mode} Sigma_min")
    sigma_max = _finite_positive(m["Sigma_lambda_max_upper"], f"{mode} Sigma_max")
    omega_min = _finite_positive(m["word_noise_Omega_lambda_min_lower"], f"{mode} Omega_min")

    # The broad-sector primitives intentionally expose source-independent
    # geometry.  Keep the composition conservative and auditable: convert the
    # exact nonlinear remainder coefficient into the same relative P3 metric by
    # paying the P3 condition-number factor.  The Cayley information-retention
    # factor then scales the usable linear margin over the whole sector.
    retain = _finite_positive(
        cayley["certificate"]["vector_information_retention_fraction_lower"],
        "Cayley information retention",
    )
    rem = _finite_positive(
        remainder["certificate"]["homogeneous_quadratic_remainder_coefficient_upper"],
        "nonlinear remainder coefficient",
    )
    cond = sigma_max / sigma_min
    linear_margin = delta * retain
    nonlinear_penalty = rem * cond * TARGET_SECTOR_RAD
    signed_margin = linear_margin - nonlinear_penalty

    # If V+ <= (1-delta_eff)V in the homogeneous word inequality, rho is the
    # induced energy contraction factor.  This is a promotion diagnostic until
    # signed_margin is strictly positive in both modes.
    rho_upper = 1.0 - signed_margin
    return {
        "p3_delta_lower": delta,
        "Sigma_lambda_min_lower": sigma_min,
        "Sigma_lambda_max_upper": sigma_max,
        "word_noise_Omega_lambda_min_lower": omega_min,
        "p3_condition_number_upper": cond,
        "cayley_information_retention_fraction_lower": retain,
        "nonlinear_remainder_coefficient_upper": rem,
        "linear_margin_after_sector_geometry_lower": linear_margin,
        "nonlinear_sector_penalty_upper": nonlinear_penalty,
        "signed_sector_margin_lower": signed_margin,
        "rho_homogeneous_upper": rho_upper,
        "strict_contraction": bool(signed_margin > 0.0 and rho_upper < 1.0),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    # PR #460 introduced the dependency-preserving source process primitive.
    # Use it directly rather than falling back to the older split_x_cell path.
    P3.split_x_cell = SP.split_x_cell
    p3 = P3.build(domain_path)
    p3_failures = P3.validate(p3)
    cayley = CAYLEY.build(domain_path)
    remainder = REMAINDER.build(domain_path)
    timing = TIMING.build(domain_path)

    modes = {mode: _extract_mode_budget(mode, p3, cayley, remainder) for mode in ("H", "A")}
    validation_failures = list(p3_failures)
    for mode, m in modes.items():
        if not m["strict_contraction"]:
            validation_failures.append(
                f"{mode}: broad-sector signed-Joseph composition is not yet contractive "
                f"(margin={m['signed_sector_margin_lower']:.17g}, "
                f"rho={m['rho_homogeneous_upper']:.17g})"
            )

    return {
        "schema": SCHEMA,
        "certificate": "OU3_P4_FULL_08RAD_SIGNED_JOSEPH_COMPOSITION",
        "domain": str(domain_path),
        "sector_radius_rad": TARGET_SECTOR_RAD,
        "source_complete_timing": timing,
        "modes": modes,
        "p3_validation_failures": p3_failures,
        "validation_failures": validation_failures,
        "validation_pass": not validation_failures,
        "theorem_promotion_allowed": not validation_failures,
        "note": (
            "This is the unique retained P3->P4 promotion gate.  A failing "
            "numerical margin is evidence that the composition bound must be "
            "tightened; it must not be bypassed by a fixed-schedule or sampled "
            "diagnostic route."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = build(args.domain)
    text = json.dumps(out, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if out["validation_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
