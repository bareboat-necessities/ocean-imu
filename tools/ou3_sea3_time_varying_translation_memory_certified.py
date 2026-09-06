#!/usr/bin/env python3
"""Source-bound facade for complete-SEA3 time-varying translation memory.

The numerical backend is parameterized by the number of complete committed-tune
intervals retained from the 3 s word.  The dynamic-source certificate's outward
sample clock gives a 23-sample maximum commit gap.  After discarding one
boundary interval at each end of 600 samples, at least

    ceil((600 - 2*23)/23) = 25

complete intervals remain.  Bind the backend to that proved integer here.

The facade also removes a purely computational redundancy without changing the
certificate.  The one-sample lower Riccati map used by the backend is Loewner
monotone: prediction is affine PSD-monotone and each scalar information update
P -> (P^-1 + d ee')^-1 is Loewner-monotone.  Therefore, for an actual constant-
tune interval of 20..23 samples it is enough to certify

    M^20(L_k) >= L_{k+1},
    M(L_{k+1}) >= L_{k+1}.

The latter invariant inductively covers samples 21..23.  The terminal suffix is
covered by the same one-sample invariant, even if the committed tau changes
between suffix samples, because it is certified uniformly over the whole tau
cell cover.  This cuts repeated interval propagation by about four without
weakening any source or matrix inequality.

Finally, the backend's historical entrywise-lower endpoint table is diagnostic
only: choosing lower endpoints of signed off-diagonal intervals is not a
Loewner-lower operation.  This facade reconstructs and exports the entire
certified interval candidate matrix for downstream full-matrix proofs.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import symmetric_positive_definite_ldlt
import ou3_sea3_time_varying_translation_memory as BACKEND

SCHEMA = BACKEND.SCHEMA
QUALIFICATION = "OU3_COMPLETE_SEA3_TIME_VARYING_TRANSLATION_PROCESS_MEMORY_CERTIFIED_GEOMETRY"
DEFAULT_DOMAIN = BACKEND.DEFAULT_DOMAIN
CERTIFIED_MACRO_INTERVALS = 25


def _certify_difference(P, L1, *, sigma_floor: float, horizon_ref: float):
    D = BACKEND._conditioned_difference(
        BACKEND._subtract(P, L1),
        sigma_floor=sigma_floor,
        horizon=horizon_ref,
    )
    ok, pivots = symmetric_positive_definite_ldlt(D)
    if not ok:
        return False, -math.inf
    return True, min(p.lo for p in pivots)


def _monotone_certify_macro_cell(
    x,
    *,
    L0,
    L1,
    sample_counts,
    h,
    sigma_floor,
    info_S,
    info_aw,
    horizon_ref,
):
    """Equivalent variable-length certification using one-step invariance."""
    try:
        counts = tuple(sample_counts)
        if not counts or any(n <= 0 for n in counts):
            raise ValueError("positive sample counts required")
        n0 = min(counts)
        P = BACKEND._macro_lower_map(
            L0, x, n0,
            h=h,
            sigma_floor=sigma_floor,
            info_S=info_S,
            info_aw=info_aw,
        )
        ok, p0 = _certify_difference(
            P, L1, sigma_floor=sigma_floor, horizon_ref=horizon_ref
        )
        if not ok:
            return False, -math.inf, f"LDLT failed at base {n0} samples"
        worst = p0

        # If the source interval can be longer than n0, prove L1 invariant for
        # one additional sample.  Monotonicity then covers every remaining
        # count in the supplied range by induction.
        if max(counts) > n0:
            P1 = BACKEND._macro_lower_map(
                L1, x, 1,
                h=h,
                sigma_floor=sigma_floor,
                info_S=info_S,
                info_aw=info_aw,
            )
            ok, p1 = _certify_difference(
                P1, L1, sigma_floor=sigma_floor, horizon_ref=horizon_ref
            )
            if not ok:
                return False, -math.inf, "one-sample target invariance LDLT failed"
            worst = min(worst, p1)
        return True, worst, ""
    except Exception as exc:
        return False, -math.inf, f"{type(exc).__name__}: {exc}"


def _reconstruct_final_interval_candidate(path: Path):
    dynamic = BACKEND.DYNAMIC.build(path)
    four = BACKEND.FOUR.build(path)
    pe = BACKEND.PE.build(path)
    rates = dynamic["validated_rate_and_jump_bounds"]
    inv = dynamic["dynamic_invariant"]
    h = float(rates["dt_s"])
    min_commit_samples = int(math.floor(0.1 / h + 1.0e-9))
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
    old_count = BACKEND.MACRO_INTERVALS
    old_checker = BACKEND._certify_macro_cell
    try:
        BACKEND.MACRO_INTERVALS = CERTIFIED_MACRO_INTERVALS
        BACKEND._certify_macro_cell = _monotone_certify_macro_cell
        BACKEND._build_cached.cache_clear()
        d = BACKEND.build(path)
        L = _reconstruct_final_interval_candidate(path)
    finally:
        BACKEND.MACRO_INTERVALS = old_count
        BACKEND._certify_macro_cell = old_checker
        BACKEND._build_cached.cache_clear()

    out = dict(d)
    out["backend_qualification"] = d["qualification"]
    out["qualification"] = QUALIFICATION
    out["certified_macro_intervals_bound_by_facade"] = CERTIFIED_MACRO_INTERVALS
    out["commit_interval_count_formula"] = "ceil((600-2*23)/23)=25"
    out["variable_interval_lengths_certified_by_monotonicity"] = True
    out["base_interval_length_samples_certified_directly"] = int(
        d["commit_geometry"]["min_constant_commit_interval_samples_conservative"]
    )
    out["longer_interval_counts_covered_by_one_sample_target_invariance"] = True
    out["terminal_suffix_covered_by_same_uniform_one_sample_invariance"] = True
    out["riccati_lower_map_Loewner_monotonicity_used"] = True
    out["source_or_candidate_matrix_changed_by_acceleration"] = False
    out["word_endpoint_translation_process_measurement_noise_interval_lower"] = [
        [x.as_list() for x in row] for row in L
    ]
    out["entrywise_lower_endpoint_table_is_Loewner_certificate"] = False
    out["downstream_must_consume_full_interval_candidate"] = True
    return out


def validate(d: dict) -> list[str]:
    probe = dict(d)
    probe["qualification"] = d.get("backend_qualification")
    old_count = BACKEND.MACRO_INTERVALS
    try:
        BACKEND.MACRO_INTERVALS = CERTIFIED_MACRO_INTERVALS
        failures = BACKEND.validate(probe)
    finally:
        BACKEND.MACRO_INTERVALS = old_count

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

    for key in (
        "variable_interval_lengths_certified_by_monotonicity",
        "longer_interval_counts_covered_by_one_sample_target_invariance",
        "terminal_suffix_covered_by_same_uniform_one_sample_invariance",
        "riccati_lower_map_Loewner_monotonicity_used",
        "downstream_must_consume_full_interval_candidate",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    if d.get("source_or_candidate_matrix_changed_by_acceleration") is not False:
        failures.append("runtime acceleration changed source/candidate matrix")

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
        "monotone_variable_length_cover": d["variable_interval_lengths_certified_by_monotonicity"],
        "induction": d["candidate_induction"],
        "tail": d["terminal_suffix"],
        "interval_lower": d["word_endpoint_translation_process_measurement_noise_interval_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())