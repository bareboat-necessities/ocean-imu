#!/usr/bin/env python3
"""Same-history translation covariance upper from the frozen P2 V1 interface.

The retained ``translation_upper`` theorem is monotone in four source-history
quantities: the largest pseudo-update cadence, largest applied sigma^2, largest
OU process intensity q_c, and largest S-measurement variance.  The old
source-uniform use formed these from independent global source extrema.  This
module instead obtains all four from ONE legal
``OU3_P2_CORRELATED_STAGE_TRANSFER_V1`` history and evaluates the same
observability/covariance argument with those path statistics.

For a source-varying history, replacing the time-varying coefficients by maxima
actually attained somewhere on that same history is conservative:

* with progress-preserving pseudo-period retargeting, every S-observation gap is bounded by the path maximum cadence plus one sample;
* q_c(t) <= max_path q_c, so process covariance is bounded by the constant
  intensity maximum over the same interval;
* sigma(t)^2 <= max_path sigma^2 for nuisance/initial-wave terms; and
* every S measurement covariance is <= the path maximum R_S variance.

The summarized source history MUST cover the complete covariance-observability
word produced by these same statistics.  A short history is rejected rather
than silently extrapolated.  This preserves history correlation without using
replay values or independently mixing global source extrema.

This module validates the single-history sufficient-statistic theorem and emits
representative legal full-word histories.  It does not yet enumerate all source
histories and cannot promote P3/P4/P5 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p2_correlation_path_memory as CORR
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def _path_segments(start_pair: tuple[int, int], transitions: list[tuple[int, int]], rt) -> list[dict]:
    c, s = map(int, start_pair)
    if not CORR.legal_pair(c, s, rt):
        raise ValueError("start pair is not reachable in P2 V1")
    out: list[dict] = []
    for gap, t in transitions:
        tr = CORR.transition(c, s, int(gap), int(t), rt)
        out.append(dict(tr["following_segment"]))
        c, s = s, int(t)
    if not out:
        raise ValueError("at least one applied V1 segment is required")
    return out


def summarize_segments(segments: list[dict], sched: dict) -> dict:
    if not segments:
        raise ValueError("nonempty source-history segments required")
    cadence_hi = 0.0
    cadence_lo = math.inf
    sigma2_hi = 0.0
    qc_hi = 0.0
    rS_variance_hi = 0.0
    total_duration_lo = 0.0
    total_duration_hi = 0.0
    axis_hi = max(sched.get("R_S_axis_std_factors", BASE.source_rs_axis_std_factors()))

    source_nodes: list[int] = []
    for seg in segments:
        if seg.get("tau_sigma_R_S_from_same_physical_cell") is not True:
            raise RuntimeError("P2 segment lost same-cell tau/sigma/R_S correlation")
        source_nodes.append(int(seg["applied_source_node"]))
        period_lo, period_hi = map(float, seg["pseudo_update_period_s"])
        cadence_lo = min(cadence_lo, period_lo)
        cadence_hi = max(cadence_hi, period_hi)
        sigma2_hi = max(sigma2_hi, float(seg["sigma_squared"][1]))
        qc_hi = max(qc_hi, float(seg["q_c_m2ps5"][1]))
        rs_hi = float(seg["R_S_filter_std"][1])
        rS_variance_hi = max(rS_variance_hi, BASE.up((rs_hi * axis_hi) ** 2))
        total_duration_lo = BASE.down(total_duration_lo + float(seg["duration_s"][0]))
        total_duration_hi = BASE.up(total_duration_hi + float(seg["duration_s"][1]))

    if not all(math.isfinite(x) and x > 0.0 for x in (cadence_lo, cadence_hi, sigma2_hi, qc_hi, rS_variance_hi)):
        raise RuntimeError("same-history source summary lost positive finite bounds")
    return {
        "source_nodes": source_nodes,
        "segments": len(segments),
        "history_duration_s": [total_duration_lo, total_duration_hi],
        "pseudo_update_cadence_s": [BASE.down(cadence_lo), BASE.up(cadence_hi)],
        "sigma_squared_upper": BASE.up(sigma2_hi),
        "q_c_upper": BASE.up(qc_hi),
        "S_measurement_variance_upper": BASE.up(rS_variance_hi),
        "all_statistics_from_one_legal_P2_history": True,
        "independent_global_source_extrema_used": False,
    }


def translation_upper_from_summary(summary: dict, Tpe: float, sched: dict,
                                   *, require_history_cover: bool = True) -> tuple[list[float], dict]:
    if summary.get("all_statistics_from_one_legal_P2_history") is not True:
        raise ValueError("same-history source summary required")
    if summary.get("independent_global_source_extrema_used") is not False:
        raise ValueError("independent global source extrema are forbidden")
    if sched.get("pseudo_period_retarget_progress_preserving") is not True:
        raise RuntimeError(
            "same-history finite S-observation gap requires progress-preserving pseudo-period retargeting"
        )

    h = float(sched["dt_s"])
    cadence = list(map(float, summary["pseudo_update_cadence_s"]))
    sigma2 = float(summary["sigma_squared_upper"])
    qc = float(summary["q_c_upper"])
    rmax = float(summary["S_measurement_variance_upper"])
    for name, x in (("h",h),("Tpe",Tpe),("sigma2",sigma2),("qc",qc),("rmax",rmax)):
        if not (math.isfinite(x) and x > 0.0):
            raise ValueError(f"{name} must be positive finite")

    gap = BASE.up(cadence[1] + h)
    spacing = BASE.up(max(Tpe, 2.0 * gap))
    Tobs = BASE.up(2.0 * spacing + gap)
    Tword = BASE.up(Tobs + Tpe)
    history_lo = float(summary["history_duration_s"][0])
    history_covers_word = history_lo >= Tword
    if require_history_cover and not history_covers_word:
        raise ValueError(
            f"same-history source summary covers only {history_lo:.17g}s but covariance word requires {Tword:.17g}s"
        )

    Binv = BASE.integrator_inverse(gap, spacing)
    s_nuis = BASE.up(sigma2 * (Tobs ** 3 / 6.0) ** 2)
    s_proc = BASE.up(qc * Tobs ** 7 / 252.0)
    rstack = BASE.up(3.0 * (rmax + s_nuis + s_proc))
    R = [[BASE.I(rstack if i == j else 0.0) for j in range(3)] for i in range(3)]
    Cspv = BASE.matrix_symmetric_hull(BASE.matrix_mul(BASE.matrix_mul(Binv, R), BASE.matrix_transpose(Binv)))
    order = (2, 1, 0)
    Cvps = [[Cspv[order[i]][order[j]] for j in range(3)] for i in range(3)]
    t = Interval.outward_bounds(0.0, Tword)
    F = [
        [BASE.I(1), BASE.I(0), BASE.I(0)],
        [t, BASE.I(1), BASE.I(0)],
        [BASE.I(0.5) * t.square(), t, BASE.I(1)],
    ]
    Cend = BASE.matrix_symmetric_hull(BASE.matrix_mul(BASE.matrix_mul(F, Cvps), BASE.matrix_transpose(F)))
    u = BASE.diagonal_dominator(Cend)
    variances = [
        BASE.up(sigma2 * Tword * Tword + qc * Tword ** 3 / 3.0),
        BASE.up(sigma2 * Tword ** 4 / 4.0 + qc * Tword ** 5 / 20.0),
        BASE.up(sigma2 * Tword ** 6 / 36.0 + qc * Tword ** 7 / 252.0),
        BASE.up(sigma2),
    ]
    roots = [math.sqrt(v) for v in variances]
    total = BASE.up(sum(roots))
    noise = [BASE.up(r * total) for r in roots]

    gap_lo = BASE.down(cadence[0] + h)
    spacing_lo = BASE.down(max(Tpe, 2.0 * gap_lo))
    Tword_lo = BASE.down(BASE.down(2.0 * spacing_lo + gap_lo) + Tpe)

    return [BASE.up(u[i] + noise[i]) for i in range(3)] + [noise[3]], {
        "cadence_s": cadence,
        "gap_s_upper": gap,
        "word_horizon_s_upper": Tword,
        "word_horizon_s_lower": Tword_lo,
        "summarized_history_duration_lower_s": history_lo,
        "summarized_history_covers_covariance_word": history_covers_word,
        "sigma_squared_upper_from_same_history": sigma2,
        "q_c_upper_from_same_history": qc,
        "S_measurement_variance_upper_from_same_history": rmax,
    }


def _representative_history(rt, seed: int, length: int, gap: int) -> tuple[tuple[int,int], list[tuple[int,int]]]:
    c = int(seed)
    outs = sorted(int(x) for x in rt["union_successors"][c])
    if not outs:
        raise RuntimeError("seed source has no reachable staged node")
    s = outs[len(outs)//2]
    start = (c, s)
    transitions: list[tuple[int,int]] = []
    for _ in range(int(length)):
        succ = CORR.successors(s, int(gap), rt)
        if not succ:
            raise RuntimeError("representative history hit dead gap-labelled source state")
        t = succ[len(succ)//2]
        transitions.append((int(gap), int(t)))
        c, s = s, t
    return start, transitions


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("same-history covariance upper must not be trajectory fitted")
    rt = CORR.runtime(path)
    corr = CORR.build(path)
    cf = CORR.validate(corr)
    if cf:
        raise RuntimeError(f"P2 V1 correlation certificate failed: {cf}")
    sched = BASE.source_schedule()
    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")

    # Each representative history has a lower physical duration above the
    # largest deployed ~3.17 s covariance word before it is accepted.
    cases = ((0,50,13),(137,31,21),(729,25,26),(799,31,21))
    rows = []
    for seed, length, gap in cases:
        start, trans = _representative_history(rt, seed, length, gap)
        segments = _path_segments(start, trans, rt)
        summary = summarize_segments(segments, sched)
        upper, timing = translation_upper_from_summary(summary, Tpe, sched, require_history_cover=True)
        rows.append({
            "start_pair": list(start),
            "transition_count": len(trans),
            "gap_samples": gap,
            "history_summary": summary,
            "Sigma_translation_diagonal_upper": upper,
            "timing": timing,
        })

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_P2_V1_SAME_HISTORY_TRANSLATION_COVARIANCE_UPPER",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "same_history_sufficient_statistics_used": True,
        "independent_cartesian_tau_sigma_R_S_extrema_used": False,
        "retained_translation_observability_theorem_reused": True,
        "progress_preserving_scheduler_required_for_gap_bound": True,
        "monotone_path_maxima_only": True,
        "full_covariance_word_history_required": True,
        "all_representative_histories_cover_their_covariance_word": all(
            row["timing"]["summarized_history_covers_covariance_word"] for row in rows
        ),
        "full_source_history_family_enumerated": False,
        "representative_rows": rows,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "next_obligation": (
            "propagate a Pareto/invariant enclosure of these same-history sufficient statistics over all P2 V1 histories covering the canonical covariance word, and pair each history class with the recent full-matrix selected-process lower before taking a worst-case ratio"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_P2_V1_SAME_HISTORY_TRANSLATION_COVARIANCE_UPPER":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed",
        "same_history_sufficient_statistics_used",
        "retained_translation_observability_theorem_reused",
        "progress_preserving_scheduler_required_for_gap_bound",
        "monotone_path_maxima_only", "full_covariance_word_history_required",
        "all_representative_histories_cover_their_covariance_word",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "independent_cartesian_tau_sigma_R_S_extrema_used",
        "full_source_history_family_enumerated", "P3_PROMOTED", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("same-history upper lost P2 V1 binding")
    rows = d.get("representative_rows", [])
    if not rows:
        f.append("no representative legal histories emitted")
    for row in rows:
        s = row.get("history_summary", {})
        if s.get("all_statistics_from_one_legal_P2_history") is not True:
            f.append("representative row lost same-history statistics")
        if s.get("independent_global_source_extrema_used") is not False:
            f.append("representative row used independent global source extrema")
        if row.get("timing", {}).get("summarized_history_covers_covariance_word") is not True:
            f.append("representative source history does not cover covariance word")
        u = row.get("Sigma_translation_diagonal_upper", [])
        if len(u) != 4 or any(not (isinstance(x,(int,float)) and math.isfinite(float(x)) and float(x)>0.0) for x in u):
            f.append("representative row has invalid translation covariance upper")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P2_correlation_interface_version": d["P2_correlation_interface_version"],
        "representative_histories": len(d["representative_rows"]),
        "all_cover_word": d["all_representative_histories_cover_their_covariance_word"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
