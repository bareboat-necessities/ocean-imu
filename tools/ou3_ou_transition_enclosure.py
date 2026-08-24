#!/usr/bin/env python3
"""Validated scalar OU transition enclosure for the configured OU-III runtime.

This stage closes the cancellation-sensitive scalar transition formulas used by
``KalmanOUCoreMath.h`` over the implementation-derived tau range. It consumes
only source constants and the exact-rational/outward-rounded transcendental
backend. It deliberately does not claim the matrix Riccati word or nonlinear
SO(3) theorem yet.

The current shipping wrapper is designed around ``FREQ_SMOOTHER_DT`` (200 Hz)
but ``updateTime(dt, ...)`` itself accepts arbitrary positive finite ``dt``.
Therefore this certificate is explicitly a configured-runtime theorem layer,
not an unconditional theorem for every caller-supplied positive ``dt``. A
later deployment gate must either bind the runtime scheduler to this interval or
add an implementation guard before calling the full source domain closed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
from ou3_validated_transcendentals import ou_discrete_coefficients

SCHEMA = 1
QUALIFICATION = "VALIDATED_CONFIGURED_RUNTIME_OU_SCALAR_TRANSITION"


def _I(bounds) -> Interval:
    return Interval(float(bounds[0]), float(bounds[1]))


def build(header: Path = SOURCE.DEFAULT_HEADER.resolve()) -> dict:
    header = header.resolve()
    source = SOURCE.build(header)
    text = header.read_text(encoding="utf-8")
    dt = SOURCE.parse_const(text, "FREQ_SMOOTHER_DT")
    h = Interval.outward_bounds(dt, dt)
    tau = _I(source["validated_parameter_box"]["continuous_parameters"]["tau_aw_s"])
    coeff = ou_discrete_coefficients(h, tau)
    x = coeff["x"]
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "implementation_header": source["implementation_header"],
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "arithmetic_backend": (
            "EXACT_BINARY32_SOURCE_PLUS_EXACT_RATIONAL_TAYLOR_PLUS_BINARY64_NEXTAFTER"
        ),
        "runtime_timing_contract": {
            "kind": "CONFIGURED_WRAPPER_SAMPLE_INTERVAL",
            "source_constant": "FREQ_SMOOTHER_DT",
            "imu_dt_s": h.as_list(),
            "arbitrary_positive_api_dt_covered": False,
            "closure_requirement": (
                "deployment must bind scheduler dt to this interval or enforce a source guard"
            ),
        },
        "tau_aw_s": tau.as_list(),
        "x_h_over_tau": x.as_list(),
        "x_supported_by_validated_series": bool(x.lo >= 0.0 and x.hi <= 1.0),
        "coefficients": {name: value.as_list() for name, value in coeff.items() if name != "x"},
        "continuous_ou_scalar_transition_enclosed": True,
        "continuous_matrix_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate these scalar transition intervals through the implemented attitude/linear "
            "matrix, Joseph/Riccati corrections and source-reachable measurement schedule"
        ),
    }


def validate(payload: dict, header: Path = SOURCE.DEFAULT_HEADER.resolve()) -> list[str]:
    expected = build(header.resolve())
    failures = []
    for key in (
        "schema",
        "qualification",
        "implementation_header",
        "source_generated_not_trajectory_fit",
        "validated_arithmetic",
        "outward_rounded",
        "runtime_timing_contract",
        "tau_aw_s",
        "x_h_over_tau",
        "x_supported_by_validated_series",
        "coefficients",
        "continuous_ou_scalar_transition_enclosed",
        "continuous_matrix_word_enclosed",
        "nonlinear_word_enclosed",
        "theorem_promotion",
    ):
        if payload.get(key) != expected.get(key):
            failures.append(f"field {key!r} does not match source-derived validated enclosure")
    if payload.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("scalar transition stage must not promote the theorem")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=SOURCE.DEFAULT_HEADER)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    payload = build(args.header.resolve())
    failures = validate(payload, args.header.resolve())
    payload["validation_pass"] = not failures
    payload["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
