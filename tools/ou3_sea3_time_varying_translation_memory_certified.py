#!/usr/bin/env python3
"""Source-bound facade for the complete-SEA3 time-varying translation memory.

The numerical backend is parameterized by the number of complete committed-tune
intervals retained from the 3 s word.  The dynamic-source certificate's outward
sample clock gives a 23-sample maximum commit gap.  After discarding one
boundary interval at each end of 600 samples, at least

    ceil((600 - 2*23)/23) = 25

complete intervals remain.  Bind the backend to that proved integer here.
Nothing else in its arithmetic/source semantics is changed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_time_varying_translation_memory as BACKEND

SCHEMA = BACKEND.SCHEMA
QUALIFICATION = "OU3_COMPLETE_SEA3_TIME_VARYING_TRANSLATION_PROCESS_MEMORY_CERTIFIED_GEOMETRY"
DEFAULT_DOMAIN = BACKEND.DEFAULT_DOMAIN
CERTIFIED_MACRO_INTERVALS = 25


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    old = BACKEND.MACRO_INTERVALS
    try:
        BACKEND.MACRO_INTERVALS = CERTIFIED_MACRO_INTERVALS
        BACKEND._build_cached.cache_clear()
        d = BACKEND.build(domain_path)
    finally:
        BACKEND.MACRO_INTERVALS = old
        BACKEND._build_cached.cache_clear()
    out = dict(d)
    out["backend_qualification"] = d["qualification"]
    out["qualification"] = QUALIFICATION
    out["certified_macro_intervals_bound_by_facade"] = CERTIFIED_MACRO_INTERVALS
    out["commit_interval_count_formula"] = "ceil((600-2*23)/23)=25"
    return out


def validate(d: dict) -> list[str]:
    probe = dict(d)
    probe["qualification"] = d.get("backend_qualification")
    failures = BACKEND.validate(probe)
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("facade schema/qualification mismatch")
    if d.get("certified_macro_intervals_bound_by_facade") != CERTIFIED_MACRO_INTERVALS:
        failures.append("certified macro-interval count changed")
    g = d.get("commit_geometry", {})
    if int(g.get("max_commit_interval_samples_certified", 0)) != 23:
        failures.append("dynamic source max commit gap is not 23 samples")
    if int(g.get("complete_constant_tune_intervals_lower", 0)) < CERTIFIED_MACRO_INTERVALS:
        failures.append("source geometry does not guarantee 25 complete intervals")
    if int(g.get("intervals_retained", 0)) != CERTIFIED_MACRO_INTERVALS:
        failures.append("backend did not retain exactly 25 intervals")
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
        "closed": d["full_4x4_time_varying_translation_memory_closed"],
        "commit_geometry": d["commit_geometry"],
        "induction": d["candidate_induction"],
        "tail": d["terminal_suffix"],
        "lower": d["word_endpoint_translation_process_measurement_noise_covariance_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
