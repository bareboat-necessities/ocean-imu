#!/usr/bin/env python3
"""Canonical source-uniform OU-III P3 certificate entry point.

PR #460 repaired the one-axis integrated-OU process primitive, but the workflow
only selected that implementation by monkey-patching ``split_x_cell`` inside an
inline CI script.  That made the advertised P3 certificate depend on how it was
invoked rather than on a stable proof producer.

This module is the retained P3 entry point.  It deliberately reuses the existing
source-reachable matrix construction while binding its process-cell splitter to
the dependency-preserving implementation in ``ou3_p3_scaled_process`` before any
certificate is built.  The resulting artifact therefore has one reproducible
meaning from the command line, tests, or CI.

No nonlinear P4 claim is made here.  The output exposes the H/A delta, Sigma and
Omega quantities that the signed-Joseph P4 composition must consume.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p3_scaled_process as SCALED
import ou3_source_reachable_matrix_p3 as MATRIX

SCHEMA = 1
DEFAULT_DOMAIN = MATRIX.DEFAULT_DOMAIN


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()

    # The old local splitter in MATRIX is intentionally not a proof dependency.
    # Bind the repaired, dependency-preserving process primitive before MATRIX
    # constructs any x=h/tau source cell.
    original = MATRIX.split_x_cell
    MATRIX.split_x_cell = SCALED.split_x_cell
    try:
        matrix = MATRIX.build(path)
    finally:
        MATRIX.split_x_cell = original

    failures = MATRIX.validate(matrix)
    modes = {}
    for mode in ("H", "A"):
        src = matrix.get("modes", {}).get(mode, {})
        modes[mode] = {
            "dimension": src.get("dimension"),
            "relative_Riccati_injection_margin_lower": src.get(
                "relative_Riccati_injection_margin_lower"
            ),
            "Sigma_lambda_min_lower": src.get("Sigma_lambda_min_lower"),
            "Sigma_lambda_max_upper": src.get("Sigma_lambda_max_upper"),
            "word_noise_Omega_lambda_min_lower": src.get(
                "word_noise_Omega_lambda_min_lower"
            ),
            "prefix_information_gain_upper": src.get("prefix_information_gain_upper"),
            "comparison_scale_diagonal_squared": src.get("matrix_comparison", {}).get(
                "comparison_scale_diagonal_squared"
            ),
            "Sigma_diagonal_upper": src.get("matrix_comparison", {}).get(
                "Sigma_diagonal_upper"
            ),
            "pass": src.get("pass") is True,
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_SOURCE_UNIFORM_DEPENDENCY_PRESERVING_MATRIX_CERTIFICATE",
        "source_generated_not_trajectory_fit": matrix.get(
            "source_generated_not_trajectory_fit"
        ) is True,
        "dependency_preserving_scaled_process_backend": True,
        "scaled_process_qualification": "FACTORED_SOURCE_BRANCH_WITH_CORRELATED_EXACT_EXPONENTIAL_SERIES",
        "matrix_certificate": matrix,
        "modes": modes,
        "P3_LINEAR_CERTIFICATE_ESTABLISHED": not failures,
        "P4_NONLINEAR_WORD_ESTABLISHED_HERE": False,
        "validation_failures": failures,
        "pass": not failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("validation_failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("P3 source-generation flag missing")
    if d.get("dependency_preserving_scaled_process_backend") is not True:
        failures.append("repaired scaled-process backend is not bound")
    if d.get("P3_LINEAR_CERTIFICATE_ESTABLISHED") is not True:
        failures.append("P3 linear certificate not established")
    if d.get("P4_NONLINEAR_WORD_ESTABLISHED_HERE") is not False:
        failures.append("P3 producer must not promote P4")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        for key in (
            "relative_Riccati_injection_margin_lower",
            "Sigma_lambda_min_lower",
            "Sigma_lambda_max_upper",
            "word_noise_Omega_lambda_min_lower",
        ):
            value = m.get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)) or float(value) <= 0.0:
                failures.append(f"{mode}.{key} is not finite positive")
        if m.get("prefix_information_gain_upper") != 1.0:
            failures.append(f"{mode}.prefix_information_gain_upper changed")
        if not m.get("comparison_scale_diagonal_squared"):
            failures.append(f"{mode} comparison scaling missing")
        if not m.get("Sigma_diagonal_upper"):
            failures.append(f"{mode} directional Sigma upper bound missing")
        if m.get("pass") is not True:
            failures.append(f"{mode} P3 mode failed")
    return list(dict.fromkeys(failures))


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
        "qualification": d["qualification"],
        "P3": d["P3_LINEAR_CERTIFICATE_ESTABLISHED"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
