#!/usr/bin/env python3
"""Source-bound facade for the complete-SEA3 time-varying translation memory.

The numerical backend is parameterized by the number of complete committed-tune
intervals retained from the 3 s word.  The dynamic-source certificate's outward
sample clock gives a 23-sample maximum commit gap.  After discarding one
boundary interval at each end of 600 samples, at least

    ceil((600 - 2*23)/23) = 25

complete intervals remain.  Bind the backend to that proved integer here.
Nothing else in its arithmetic/source semantics is changed.

The backend historically exposed a convenience matrix made from the lower
endpoint of each final interval entry.  That table is diagnostic only: choosing
entrywise lower endpoints is not, in general, a Loewner-lower operation when
off-diagonal entries are signed.  This facade therefore also reconstructs and
exports the *entire certified interval candidate matrix*.  Downstream matrix
proofs must consume that interval matrix (or prove a particular point selection
from it), never the diagnostic entrywise-lower table.
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


def _reconstruct_final_interval_candidate(path: Path):
    dynamic = BACKEND.DYNAMIC.build(path)
    four = BACKEND.FOUR.build(path)
    pe = BACKEND.PE.build(path)
    rates = dynamic["validated_rate_and_jump_bounds"]
    inv = dynamic["dynamic_invariant"]
    h = float(rates["dt_s"])
    min_commit_samples = int(__import__("math").floor(0.1 / h + 1.0e-9))
    sigma_floor = float(inv["sigma_aw_filter_mps2"][0])
    tau_hi = float(inv["tau_applied_s"][1])
    ra = float(pe["measurement_runtime"]["accelerometer_variance_upper"])
    rs_lo = float(inv["R_S_applied"][0])
    axis_factor_min = min(map(float, four["R_S_axis_std_factors"]))
    rs_std_min = BACKEND.TUBE.down(rs_lo * axis_factor_min)
    info_S = BACKEND.TUBE.up(1.0 / BACKEND.TUBE.down(rs_std_min * rs_std_min))
    info_aw = BACKEND.TUBE.up(3.0 / ra)
    refs = BACKEND._reference_sequence(
        h=h,
        tau_reference=tau_hi,
        sigma_floor=sigma_floor,
        macro_samples=min_commit_samples,
        info_S=info_S,
        info_aw=info_aw,
    )
    return BACKEND._candidate(refs, CERTIFIED_MACRO_INTERVALS)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    old = BACKEND.MACRO_INTERVALS
    try:
        BACKEND.MACRO_INTERVALS = CERTIFIED_MACRO_INTERVALS
        BACKEND._build_cached.cache_clear()
        d = BACKEND.build(path)
        L = _reconstruct_final_interval_candidate(path)
    finally:
        BACKEND.MACRO_INTERVALS = old
        BACKEND._build_cached.cache_clear()
    out = dict(d)
    out["backend_qualification"] = d["qualification"]
    out["qualification"] = QUALIFICATION
    out["certified_macro_intervals_bound_by_facade"] = CERTIFIED_MACRO_INTERVALS
    out["commit_interval_count_formula"] = "ceil((600-2*23)/23)=25"
    out["word_endpoint_translation_process_measurement_noise_interval_lower"] = [
        [x.as_list() for x in row] for row in L
    ]
    out["entrywise_lower_endpoint_table_is_Loewner_certificate"] = False
    out["downstream_must_consume_full_interval_candidate"] = True
    return out


def validate(d: dict) -> list[str]:
    probe = dict(d)
    probe["qualification"] = d.get("backend_qualification")
    old = BACKEND.MACRO_INTERVALS
    try:
        BACKEND.MACRO_INTERVALS = CERTIFIED_MACRO_INTERVALS
        failures = BACKEND.validate(probe)
    finally:
        BACKEND.MACRO_INTERVALS = old
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
    M = d.get("word_endpoint_translation_process_measurement_noise_interval_lower")
    if not isinstance(M, list) or len(M) != 4 or any(
        not isinstance(row, list) or len(row) != 4 for row in (M or [])
    ):
        failures.append("full interval translation-memory candidate is not 4x4")
    else:
        for row in M:
            for x in row:
                if not isinstance(x, list) or len(x) != 2 or float(x[0]) > float(x[1]):
                    failures.append("invalid interval entry in translation-memory candidate")
                    break
    if d.get("entrywise_lower_endpoint_table_is_Loewner_certificate") is not False:
        failures.append("entrywise lower diagnostic was promoted to Loewner certificate")
    if d.get("downstream_must_consume_full_interval_candidate") is not True:
        failures.append("full interval candidate is not required downstream")
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
        "interval_lower": d["word_endpoint_translation_process_measurement_noise_interval_lower"],
        "entrywise_lower_diagnostic": d["word_endpoint_translation_process_measurement_noise_covariance_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())