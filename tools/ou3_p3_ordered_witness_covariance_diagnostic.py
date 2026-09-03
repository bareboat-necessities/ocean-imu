#!/usr/bin/env python3
"""Point diagnostic for the legal P3 four-max witness with shipping S cadence.

This is deliberately non-promoting.  It keeps the exact P2 source order,
extends the 101-sample global-label witness to the 635-sample word using exact
gap-labelled edges, carries the shipping pseudo timer through source changes,
and applies only the S=0 covariance update (no accelerometer acceptance is
assumed).  Source cells use one real upper corner.  A second run uses the
nonphysical tuple formed from the four independent global maxima.  The result
measures whether time order/cadence is worth pursuing in a certified upper; it
is not itself a uniform covariance certificate.
"""
from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_covariance_upper as CUPPER
import ou3_p3_four_max_global_label_witness as WIT
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_p3_frozen_full_matrix_translation as FROZEN
import ou3_p3_pseudo_scheduler_starvation_witness as TIMER
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
PHASES = (0.0, 0.25, 0.5, 0.75, 0.999999999999)


def _mid(x):
    return float(x.lo + 0.5 * (x.hi - x.lo))


def _pmat(A):
    return [[_mid(x) for x in row] for row in A]


def _zero():
    return [[0.0] * 4 for _ in range(4)]


def _diag(v):
    A = _zero()
    for i, x in enumerate(v):
        A[i][i] = float(x)
    return A


def _T(A):
    return [list(x) for x in zip(*A)]


def _mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def _add(A, B):
    return [[A[i][j] + B[i][j] for j in range(4)] for i in range(4)]


def _sym(A):
    B = [[float(x) for x in row] for row in A]
    for i in range(4):
        for j in range(i + 1, 4):
            x = 0.5 * (B[i][j] + B[j][i])
            B[i][j] = B[j][i] = x
    return B


def _meas(P, k, R):
    den = P[k][k] + R
    if not (math.isfinite(den) and den > 0.0):
        raise RuntimeError("measurement denominator lost positivity")
    c = [P[i][k] for i in range(4)]
    return _sym([[P[i][j] - c[i] * c[j] / den for j in range(4)] for i in range(4)])


def _set_period(elapsed, period):
    """Use the exact binary32 setter transcription shared with the scheduler witness."""
    return TIMER._set_period(elapsed, period)


def _due(dt, period, elapsed):
    """Use the exact binary32 `periodic_update_due<float>` transcription."""
    return TIMER._due(dt, period, elapsed)


def _corner(node, axis_hi):
    tau = float(node["tau_s"][1])
    sigma = float(node["sigma_filter_committed_mps2"][1])
    return {
        "tau": tau,
        "sigma": sigma,
        "qc": 2.0 * sigma * sigma / tau,
        "Rstd": float(node["R_S_filter_std"][1]) * axis_hi,
        "period": TIMER._f32(node["pseudo_update_period_s"][1]),
        "node": int(node["index"]),
    }


@functools.lru_cache(maxsize=256)
def _kernel(tau, sigma, h):
    x = Interval.point(h / tau)
    F = _pmat(FROZEN._transition(x))
    Q0 = _pmat(FROZEN._scaled_Q(x))
    s2 = sigma * sigma
    return F, _sym([[s2 * Q0[i][j] for j in range(4)] for i in range(4)])


def _step(P, p, h, elapsed):
    F, Q = _kernel(p["tau"], p["sigma"], h)
    P = _sym(_add(_mm(_mm(F, P), _T(F)), Q))
    fire, elapsed = _due(h, p["period"], elapsed)
    if fire:
        P = _meas(P, 2, p["Rstd"] ** 2 / h ** 6)
    return P, elapsed, fire


def _physical_diag(P, h):
    d = (h, h * h, h ** 3, 1.0)
    out = [d[i] * d[i] * P[i][i] for i in range(4)]
    if any(not (math.isfinite(x) and x > 0.0) for x in out):
        raise RuntimeError("endpoint covariance diagonal is not positive finite")
    return out


def _z_seed(physical, h):
    d = (h, h * h, h ** 3, 1.0)
    return _diag([physical[i] / (d[i] * d[i]) for i in range(4)])


def _choose_edge(s, rt):
    self_edges, all_edges = [], []
    for gi, gap in enumerate(rt["gaps"]):
        for t0 in rt["labelled_successors"][s][gi]:
            e = (int(gap), int(t0))
            all_edges.append(e)
            if int(t0) == s:
                self_edges.append(e)
    if self_edges:
        return min(self_edges)
    if not all_edges:
        raise RuntimeError("dead P2 source in witness extension")
    return min(all_edges)


def extend_witness_to_target(path, rt, target):
    segs, n = [], 0
    for r in path:
        s, t, g = int(r["source"]), int(r["successor"]), int(r["gap_samples"])
        if t not in CORR.successors(s, g, rt):
            raise RuntimeError("four-max witness edge is no longer legal")
        n += g
        if n != int(r["cumulative_samples"]):
            raise RuntimeError("witness cumulative cost drifted")
        segs.append((s, g, g, t, True))
    s = int(path[-1]["successor"])
    while n < target:
        g, t = _choose_edge(s, rt)
        used = min(g, target - n)
        n += used
        segs.append((s, used, g, t, used == g))
        if used == g:
            s = t
    return segs


def _summary(rt, sched, target):
    st = HIST._stat_tables(rt, sched)
    g = WIT.global_rank_tuple(st["node_ranks"])
    tb = st["tables"]
    return {
        "source_nodes": [], "segments": None,
        "history_duration_s": [target["history_duration_lower_s"], target["terminal_history_duration_upper_s"]],
        "pseudo_update_cadence_s": [st["cadence_lower_global_safe"], tb[0][g[0]]],
        "sigma_squared_upper": tb[1][g[1]], "q_c_upper": tb[2][g[2]],
        "S_measurement_variance_upper": tb[3][g[3]],
        "all_statistics_from_one_legal_P2_history": True,
        "independent_global_source_extrema_used": False,
        "history_label_generated_by_exact_gap_successors": True,
        "dominance_pruning_only_removed_no_more_adverse_same_state_labels": True,
    }, list(map(int, g)), st


def _synthetic(s):
    sig2, qc = float(s["sigma_squared_upper"]), float(s["q_c_upper"])
    return {
        "tau": 2.0 * sig2 / qc, "sigma": math.sqrt(sig2), "qc": qc,
        "Rstd": math.sqrt(float(s["S_measurement_variance_upper"])),
        "period": TIMER._f32(s["pseudo_update_cadence_s"][1]), "node": None,
    }


def _run_ordered(P0, segs, rt, axis_hi, h, phase):
    p0 = _corner(rt["nodes"][segs[0][0]], axis_hi)
    elapsed = _set_period(phase * p0["period"], p0["period"])
    P, fires, prev, trace = [row[:] for row in P0], 0, None, []
    for s, used, support, t, complete in segs:
        p = _corner(rt["nodes"][s], axis_hi)
        if prev is None or s != prev:
            elapsed = _set_period(elapsed, p["period"])
        local = 0
        for _ in range(used):
            P, elapsed, fire = _step(P, p, h, elapsed)
            fires += int(fire); local += int(fire)
        trace.append({"source": s, "samples": used, "supporting_gap": support,
                      "successor": t, "complete": complete, "fires": local,
                      "tau": p["tau"], "sigma": p["sigma"], "qc": p["qc"],
                      "Rstd": p["Rstd"], "period": p["period"]})
        prev = s
    d = _physical_diag(P, h)
    return {"phase": phase, "fires": fires, "diag": d, "std": [math.sqrt(x) for x in d], "trace": trace}


def _run_fixed(P0, p, N, h, phase):
    elapsed = _set_period(phase * p["period"], p["period"])
    P, fires = [row[:] for row in P0], 0
    for _ in range(N):
        P, elapsed, fire = _step(P, p, h, elapsed); fires += int(fire)
    d = _physical_diag(P, h)
    return {"phase": phase, "fires": fires, "diag": d, "std": [math.sqrt(x) for x in d]}


def _env(rows):
    return [max(r["diag"][i] for r in rows) for i in range(4)]


def build(domain_path=DEFAULT_DOMAIN):
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text())
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("diagnostic must not be trajectory fitted")
    rt = CORR.runtime(path); sched = BASE.source_schedule(); h = TIMER._f32(rt["clock"]["dt_binary32_s"])
    target = HIST._global_word_target(domain, sched, h); N = int(target["target_samples"])
    summary, rank, stats = _summary(rt, sched, target)
    witness = WIT.shortest_global_label_witness(rt, stats["node_ranks"], N)
    if not witness["reachable_within_target"]:
        raise RuntimeError("global four-max label is no longer reachable")
    segs = extend_witness_to_target(witness["witness_path"], rt, N)
    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    old, timing = CUPPER.translation_upper_from_summary(summary, Tpe, sched, require_history_cover=True)
    old = list(map(float, old)); P0 = _z_seed(old, h); axis_hi = max(BASE.source_rs_axis_std_factors())
    ordered = [_run_ordered(P0, segs, rt, axis_hi, h, q) for q in PHASES]
    synp = _synthetic(summary); synthetic = [_run_fixed(P0, synp, N, h, q) for q in PHASES]
    oe, se = _env(ordered), _env(synthetic)
    og, sg = [old[i] / oe[i] for i in range(4)], [old[i] / se[i] for i in range(4)]
    return {
        "schema": SCHEMA, "qualification": "OU3_P3_ORDERED_FOUR_MAX_WITNESS_POINT_DIAGNOSTIC",
        "diagnostic_only": True, "trajectory_replay_used": False, "filter_changed": False,
        "declared_domain_changed": False, "canonical_gate_changed": False,
        "P2_correlation_interface_consumed": True, "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "exact_witness_source_order_retained": True, "exact_gap_labelled_legal_extension_used": True,
        "pseudo_period_change_uses_fmod_semantics": True, "periodic_update_due_shipping_semantics_transcribed": True,
        "pseudo_scheduler_numeric_type": "binary32/float",
        "accelerometer_measurement_updates_credited": False, "source_cells_use_one_real_upper_corner": True,
        "interval_certificate": False, "uniform_covariance_upper_certificate": False,
        "target_samples": N, "dt_s": h, "four_max_global_rank": rank,
        "four_max_witness_minimum_samples": witness["minimum_cost_samples"],
        "whole_word_segment_count": len(segs), "whole_word_terminal_partial_samples": segs[-1][1] if not segs[-1][4] else 0,
        "old_four_max_upper_diagonal": old, "old_four_max_upper_std": [math.sqrt(x) for x in old], "old_four_max_timing": timing,
        "initial_seed_note": "diag(old four-max upper) in fixed z metric; point comparison only, not a Loewner upper claim",
        "phase_fractions": list(PHASES), "ordered_runs": ordered, "ordered_phase_envelope_diagonal": oe,
        "ordered_phase_envelope_std": [math.sqrt(x) for x in oe], "old_upper_to_ordered_variance_gain": og,
        "synthetic_four_max_parameters": synp, "synthetic_runs": synthetic, "synthetic_phase_envelope_diagonal": se,
        "synthetic_phase_envelope_std": [math.sqrt(x) for x in se], "old_upper_to_synthetic_variance_gain": sg,
        "classification": "ORDERED_POINT_DIAGNOSTIC_COMPLETE_DO_NOT_PROMOTE",
        "matched_margin_computed": False, "P3_PROMOTED": False, "P4_PROMOTED": False, "P5_PROMOTED": False,
        "next_obligation": "select a certified time-ordered covariance-upper quotient carrying source order and pseudo-scheduler phase; do not promote from this point diagnostic",
        "failures": [],
    }


def validate(d):
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA or d.get("qualification") != "OU3_P3_ORDERED_FOUR_MAX_WITNESS_POINT_DIAGNOSTIC": f.append("schema/qualification mismatch")
    for k in ("diagnostic_only", "P2_correlation_interface_consumed", "exact_witness_source_order_retained",
              "exact_gap_labelled_legal_extension_used", "pseudo_period_change_uses_fmod_semantics",
              "periodic_update_due_shipping_semantics_transcribed", "source_cells_use_one_real_upper_corner"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("trajectory_replay_used", "filter_changed", "declared_domain_changed", "canonical_gate_changed",
              "accelerometer_measurement_updates_credited", "interval_certificate", "uniform_covariance_upper_certificate",
              "matched_margin_computed", "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    if d.get("pseudo_scheduler_numeric_type") != "binary32/float": f.append("ordered diagnostic is not using shipping float scheduler arithmetic")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION: f.append("lost P2 V1 binding")
    if int(d.get("target_samples", 0)) != 635 or int(d.get("four_max_witness_minimum_samples", 9999)) > 635: f.append("word/witness contract changed")
    for k in ("old_four_max_upper_diagonal", "ordered_phase_envelope_diagonal", "synthetic_phase_envelope_diagonal"):
        v = d.get(k, [])
        if len(v) != 4 or any(not (math.isfinite(float(x)) and float(x) > 0.0) for x in v): f.append(f"invalid {k}")
    return list(dict.fromkeys(f))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN); ap.add_argument("--output", type=Path, required=True); a = ap.parse_args()
    d = build(a.domain); vf = validate(d); d["validation_pass"] = not vf; d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(json.dumps(d, indent=2, sort_keys=True))
    print(json.dumps({"witness_samples": d["four_max_witness_minimum_samples"], "old_std": d["old_four_max_upper_std"],
                      "ordered_std": d["ordered_phase_envelope_std"], "ordered_gain": d["old_upper_to_ordered_variance_gain"],
                      "synthetic_std": d["synthetic_phase_envelope_std"], "synthetic_gain": d["old_upper_to_synthetic_variance_gain"],
                      "synthetic_parameters": d["synthetic_four_max_parameters"], "validation_failures": vf, "next": d["next_obligation"]}, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
