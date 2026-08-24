#!/usr/bin/env python3
"""Produce the first validated outward-rounded OU-III source enclosure layer.

This certificate is intentionally narrower than the deployment theorem.  It
binds the implementation-derived source-domain contract to explicit
outward-rounded binary64 intervals for every continuous parameter and timing
constant.  The later Riccati/Taylor backend must consume these intervals rather
than reconstructing ordinary floating-point extrema.

A PASS here means only that the source parameter box is source-derived and
outward-rounded.  It does *not* claim that the continuous H/A word, nonlinear
remainder, hybrid funnel, or stochastic sensitivities have been enclosed yet.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE

SCHEMA = 1
QUALIFICATION = "SOURCE_DERIVED_OUTWARD_ROUNDED_PARAMETER_BOX"


def _box(bounds: list[float]) -> list[float]:
    if len(bounds) != 2:
        raise ValueError(f"expected two endpoints, got {bounds!r}")
    return Interval.outward_bounds(float(bounds[0]), float(bounds[1])).as_list()


def build(header: Path = SOURCE.DEFAULT_HEADER.resolve()) -> dict:
    source = SOURCE.build(header.resolve())
    continuous = {
        name: _box(bounds)
        for name, bounds in source["continuous_parameters"].items()
    }
    timing = {
        name: Interval.outward_bounds(float(value), float(value)).as_list()
        for name, value in source["timing_constants_s"].items()
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_claim": source["claim"],
        "implementation_header": source["implementation_header"],
        "source_generated_not_trajectory_fit": bool(
            source["source_generated_not_trajectory_fit"]
        ),
        "source_complete_parameter_domain": bool(
            source["source_complete_parameter_domain"]
        ),
        "validated_arithmetic": True,
        "outward_rounded": True,
        "arithmetic_backend": "IEEE754_BINARY64_BASIC_OPS_NEXTAFTER_OUTWARD",
        "continuous_parameters": continuous,
        "timing_constants_s": timing,
        "discrete_source_branches": source["discrete_source_branches"],
        "hybrid_obligations": source["hybrid_obligations"],
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate this box with validated transcendental and matrix interval/Taylor "
            "arithmetic through the H/A Riccati words and nonlinear SO(3) remainder"
        ),
    }


def validate(payload: dict, header: Path = SOURCE.DEFAULT_HEADER.resolve()) -> list[str]:
    expected_source = SOURCE.build(header.resolve())
    failures: list[str] = []
    if payload.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if payload.get("qualification") != QUALIFICATION:
        failures.append("qualification mismatch")
    if payload.get("implementation_header") != expected_source["implementation_header"]:
        failures.append("implementation header mismatch")
    if payload.get("validated_arithmetic") is not True:
        failures.append("validated arithmetic flag is not true")
    if payload.get("outward_rounded") is not True:
        failures.append("outward-rounded flag is not true")
    if payload.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("parameter-box stage must not promote the theorem")

    got = payload.get("continuous_parameters", {})
    for name, source_bounds in expected_source["continuous_parameters"].items():
        bounds = got.get(name)
        if not isinstance(bounds, list) or len(bounds) != 2:
            failures.append(f"missing interval for {name}")
            continue
        I = Interval(float(bounds[0]), float(bounds[1]))
        if not I.contains(float(source_bounds[0])) or not I.contains(float(source_bounds[1])):
            failures.append(f"interval for {name} excludes source-domain endpoint")

    got_timing = payload.get("timing_constants_s", {})
    for name, value in expected_source["timing_constants_s"].items():
        bounds = got_timing.get(name)
        if not isinstance(bounds, list) or len(bounds) != 2:
            failures.append(f"missing timing interval for {name}")
            continue
        if not Interval(float(bounds[0]), float(bounds[1])).contains(float(value)):
            failures.append(f"timing interval for {name} excludes source value")
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
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
