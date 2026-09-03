#!/usr/bin/env python3
"""Exact Pareto sufficient quotient for full-word P2 V1 source histories.

Canonical P3 needs source-history correlation over the complete covariance
observability word, but carrying every staged/committed path explicitly is
unnecessary.  The retained same-history covariance upper is monotone in four
path maxima:

    pseudo cadence, sigma^2, q_c=2 sigma^2/tau, and S variance.

This module propagates exactly those four finite-valued maxima on the frozen
``OU3_P2_CORRELATED_STAGE_TRANSFER_V1`` graph.  At an intermediate state with
the same elapsed sample count and the same current staged source node, a label
that is componentwise no larger than another label can be discarded: every
future V1 continuation is identical from that source node and componentwise
max preserves the dominance relation.  This is therefore a certified
sufficient quotient, not a flat ancestor hull and not a Cartesian tuner box.

The horizon is chosen from the global retained covariance-word upper and padded
by one deployed sample.  Every emitted terminal label consequently represents
at least one legal V1 history whose physical duration covers its own (no
larger) same-history covariance word.

No trajectory/replay values, filter changes, domain shrink, or theorem-gate
relaxation are used here.  This module only certifies the source-history
quotient; it cannot promote P3/P4/P5 by itself.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from ou3_interval import Interval
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_covariance_upper as CUPPER
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
EMPTY_LABEL = (-1, -1, -1, -1)
MAX_FRONTIER_PER_STATE = 20000


def dominates(a: tuple[int, ...], b: tuple[int, ...]) -> bool:
    """Return True when a is componentwise at least as adverse as b."""
    if len(a) != len(b):
        raise ValueError("Pareto labels must have the same dimension")
    return all(int(x) >= int(y) for x, y in zip(a, b))


def pareto_insert(front: set[tuple[int, ...]], label: tuple[int, ...]) -> bool:
    """Insert one adverse-max label, removing labels it dominates."""
    x = tuple(map(int, label))
    if any(dominates(y, x) for y in front):
        return False
    dead = [y for y in front if dominates(x, y)]
    for y in dead:
        front.remove(y)
    front.add(x)
    if len(front) > MAX_FRONTIER_PER_STATE:
        raise RuntimeError(
            f"P2 V1 history Pareto frontier exceeded {MAX_FRONTIER_PER_STATE} labels; "
            "refuse silent truncation"
        )
    return True


def update_label(label: tuple[int, int, int, int],
                 node_rank: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple(max(int(a), int(b)) for a, b in zip(label, node_rank))


def _ranked(values: list[float]) -> tuple[list[float], dict[float, int]]:
    uniq = sorted(set(map(float, values)))
    if not uniq or any(not (math.isfinite(x) and x > 0.0) for x in uniq):
        raise RuntimeError("P2 V1 history statistic lost positive finite values")
    return uniq, {x: i for i, x in enumerate(uniq)}


def _stat_tables(rt: dict, sched: dict) -> dict:
    axis_hi = max(sched.get("R_S_axis_std_factors", BASE.source_rs_axis_std_factors()))
    raw = []
    for node in rt["nodes"]:
        tau_lo = float(node["tau_s"][0])
        sig_hi = float(node["sigma_filter_committed_mps2"][1])
        rs_hi = float(node["R_S_filter_std"][1])
        cadence_hi = float(node["pseudo_update_period_s"][1])
        sigma2_hi = BASE.up(sig_hi * sig_hi)
        qc_hi = BASE.up(2.0 * sig_hi * sig_hi / tau_lo)
        rs_var_hi = BASE.up((rs_hi * axis_hi) ** 2)
        raw.append((cadence_hi, sigma2_hi, qc_hi, rs_var_hi))

    cols = list(zip(*raw))
    tables = []
    maps = []
    for col in cols:
        table, rank = _ranked(list(col))
        tables.append(table)
        maps.append(rank)
    ranks = [
        tuple(maps[j][raw[i][j]] for j in range(4))
        for i in range(len(raw))
    ]
    cadence_lo_global = min(float(n["pseudo_update_period_s"][0]) for n in rt["nodes"])
    return {
        "raw": raw,
        "tables": tables,
        "node_ranks": ranks,
        "cadence_lower_global_safe": BASE.down(cadence_lo_global),
    }


def _global_word_target(domain: dict, sched: dict, h: float) -> dict:
    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    tau = Interval(*map(float, sched["tau_applied_invariant_s"]))
    sigma = Interval(*map(float, sched["sigma_aw_applied_safety"]))
    rs = Interval(*map(float, sched["R_S_applied_invariant"]))
    _, timing = BASE.translation_upper(tau, sigma, rs, Tpe, sched)
    Tword = float(timing["word_horizon_s_upper"])
    target = int(math.ceil(Tword / h - 1.0e-14)) + 1
    if target <= 0:
        raise RuntimeError("invalid full-word history target")
    h_lo = math.nextafter(float(h), -math.inf)
    history_lo = BASE.down(target * h_lo)
    if history_lo <= Tword:
        target += 1
        history_lo = BASE.down(target * h_lo)
    if history_lo <= Tword:
        raise RuntimeError("one-sample history padding did not cover global covariance word")
    h_hi = math.nextafter(float(h), math.inf)
    history_hi = BASE.up((target + 25) * h_hi)
    return {
        "Tpe_s": Tpe,
        "global_covariance_word_upper_s": Tword,
        "target_samples": target,
        "history_duration_lower_s": history_lo,
        "terminal_history_duration_upper_s": history_hi,
    }


def _start_nodes(rt: dict) -> tuple[int, ...]:
    indegree = [0] * len(rt["nodes"])
    for out in rt["union_successors"]:
        for s in out:
            indegree[int(s)] += 1
    starts = tuple(i for i, d in enumerate(indegree) if d > 0)
    if not starts:
        raise RuntimeError("P2 V1 has no legal staged start nodes")
    return starts


def _frontier_digest(terminals: dict[tuple[int, int, int], set[tuple[int, ...]]]) -> str:
    h = hashlib.sha256()
    for key in sorted(terminals):
        h.update(("%d,%d,%d:" % key).encode())
        for label in sorted(terminals[key]):
            h.update(("%d,%d,%d,%d;" % label).encode())
    return h.hexdigest()


@functools.lru_cache(maxsize=4)
def frontier_runtime(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P2 V1 history quotient must not be trajectory fitted")
    rt = CORR.runtime(path)
    corr = CORR.build(path)
    cf = CORR.validate(corr)
    if cf:
        raise RuntimeError(f"P2 V1 correlation interface failed: {cf}")

    sched = BASE.source_schedule()
    h = float(rt["clock"]["dt_binary32_s"])
    target = _global_word_target(domain, sched, h)
    stats = _stat_tables(rt, sched)
    ranks = stats["node_ranks"]
    N = int(target["target_samples"])

    layers: list[dict[int, set[tuple[int, int, int, int]]]] = [
        {} for _ in range(N)
    ]
    for s in _start_nodes(rt):
        layers[0][s] = {EMPTY_LABEL}

    terminals: dict[tuple[int, int, int], set[tuple[int, int, int, int]]] = defaultdict(set)
    generated = 0
    dominated = 0
    max_state_front = 1

    for n in range(N):
        if not layers[n]:
            continue
        for s, front in list(layers[n].items()):
            nr = ranks[s]
            for gi, gap in enumerate(rt["gaps"]):
                succ = rt["labelled_successors"][s][gi]
                if not succ:
                    raise RuntimeError("P2 V1 history quotient encountered dead labelled edge")
                n2 = n + int(gap)
                updated = [update_label(label, nr) for label in front]
                for t0 in succ:
                    t = int(t0)
                    if n2 >= N:
                        dst = terminals[(int(s), int(gap), t)]
                    else:
                        dst = layers[n2].setdefault(t, set())
                    for label in updated:
                        generated += 1
                        before = len(dst)
                        kept = pareto_insert(dst, label)
                        if not kept or len(dst) <= before:
                            dominated += 1
                    max_state_front = max(max_state_front, len(dst))

    if not terminals:
        raise RuntimeError("P2 V1 full-word history quotient produced no terminal histories")

    # Every terminal label is an actual path label: we only seed legal staged
    # nodes and only extend with exact gap-labelled successors.  Dominance
    # pruning never invents a label.  Once the full-word target is covered,
    # predecessor identity is irrelevant to the covariance upper; merge again
    # by current staged endpoint using the same continuation-safe dominance.
    labels_total = sum(len(v) for v in terminals.values())
    endpoint_frontiers: dict[int, set[tuple[int, int, int, int]]] = defaultdict(set)
    for (_, _, end), labels in terminals.items():
        dst = endpoint_frontiers[int(end)]
        for label in labels:
            pareto_insert(dst, label)
    endpoint_nodes = sorted(endpoint_frontiers)
    return {
        "path": path,
        "domain": domain,
        "rt": rt,
        "sched": sched,
        "target": target,
        "stats": stats,
        "terminals": dict(terminals),
        "terminal_keys": len(terminals),
        "terminal_labels": labels_total,
        "endpoint_frontiers": dict(endpoint_frontiers),
        "endpoint_nodes": endpoint_nodes,
        "generated_candidates": generated,
        "dominated_or_replaced_candidates": dominated,
        "max_frontier_per_state": max_state_front,
        "frontier_digest_sha256": _frontier_digest(terminals),
    }


def label_summary(label: tuple[int, int, int, int], fr: dict,
                  *, include_node: int | None = None) -> dict:
    x = tuple(map(int, label))
    if include_node is not None:
        x = update_label(x, fr["stats"]["node_ranks"][int(include_node)])
    if any(v < 0 for v in x):
        raise ValueError("full-word history label is incomplete")
    tables = fr["stats"]["tables"]
    return {
        "source_nodes": [],
        "segments": None,
        "history_duration_s": [
            float(fr["target"]["history_duration_lower_s"]),
            float(fr["target"]["terminal_history_duration_upper_s"]),
        ],
        "pseudo_update_cadence_s": [
            float(fr["stats"]["cadence_lower_global_safe"]),
            float(tables[0][x[0]]),
        ],
        "sigma_squared_upper": float(tables[1][x[1]]),
        "q_c_upper": float(tables[2][x[2]]),
        "S_measurement_variance_upper": float(tables[3][x[3]]),
        "all_statistics_from_one_legal_P2_history": True,
        "independent_global_source_extrema_used": False,
        "history_label_generated_by_exact_gap_successors": True,
        "dominance_pruning_only_removed_no_more_adverse_same_state_labels": True,
    }


def endpoint_labels(fr: dict, endpoint_node: int) -> set[tuple[int, int, int, int]]:
    t = int(endpoint_node)
    labels = fr["endpoint_frontiers"].get(t)
    if labels is None:
        return set()
    return set(labels)


def endpoint_phase_upper(endpoint_node: int, phase_samples: int, fr: dict) -> dict:
    """Uniform covariance upper after mapping each retained label separately.

    The coordinatewise envelope is formed *after* each Pareto label has been
    converted to a same-history covariance upper.  It is therefore an outer
    envelope of legal-history covariance results, not a Cartesian tau/sigma/R_S
    construction.
    """
    t = int(endpoint_node)
    r = int(phase_samples)
    if not 0 <= r <= max(fr["rt"]["gaps"]) - 1:
        raise ValueError("endpoint phase outside 0..25 samples")
    labels = endpoint_labels(fr, t)
    if not labels:
        raise RuntimeError(f"endpoint source node {t} has no full-word histories")
    Tpe = float(fr["target"]["Tpe_s"])
    uppers = []
    mapped: set[tuple[int, int, int, int]] = set()
    for label in labels:
        q = update_label(label, fr["stats"]["node_ranks"][t]) if r > 0 else label
        pareto_insert(mapped, q)
    for label in mapped:
        summary = label_summary(label, fr)
        upper, timing = CUPPER.translation_upper_from_summary(
            summary, Tpe, fr["sched"], require_history_cover=True
        )
        if timing.get("summarized_history_covers_covariance_word") is not True:
            raise RuntimeError("retained P2 V1 history label does not cover covariance word")
        uppers.append(tuple(map(float, upper)))
    env = [max(row[i] for row in uppers) for i in range(4)]
    return {
        "endpoint_source_node": t,
        "phase_samples": r,
        "input_history_labels": len(labels),
        "phase_pareto_labels": len(mapped),
        "Sigma_translation_diagonal_upper_envelope": env,
        "same_history_upper_evaluated_before_endpoint_envelope": True,
        "raw_tuner_cartesian_extrema_used": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    fr = frontier_runtime(Path(domain_path).resolve())
    sample_nodes = [fr["endpoint_nodes"][0], fr["endpoint_nodes"][len(fr["endpoint_nodes"]) // 2], fr["endpoint_nodes"][-1]]
    reps = [
        endpoint_phase_upper(t, r, fr)
        for t in sample_nodes
        for r in (0, 12, 25)
    ]
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_P2_V1_FULL_WORD_PARETO_HISTORY_QUOTIENT",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "exact_gap_labelled_successors_used": True,
        "flat_800_node_ancestor_hull_used": False,
        "raw_tuner_cartesian_extrema_used": False,
        "certified_sufficient_quotient_used": True,
        "quotient_state": "(elapsed_samples,current_staged_source,Pareto path maxima)",
        "dominance_rule": "same elapsed samples/source: discard label iff another is componentwise >= in all four adverse maxima",
        "dominance_continuation_safe": True,
        "full_source_history_family_enumerated_at_sufficient_statistic_level": True,
        "global_word_target": fr["target"],
        "terminal_history_classes": fr["terminal_keys"],
        "terminal_pareto_labels": fr["terminal_labels"],
        "endpoint_source_nodes": len(fr["endpoint_nodes"]),
        "max_frontier_per_state": fr["max_frontier_per_state"],
        "generated_candidates": fr["generated_candidates"],
        "dominated_or_replaced_candidates": fr["dominated_or_replaced_candidates"],
        "frontier_digest_sha256": fr["frontier_digest_sha256"],
        "same_history_upper_evaluated_before_endpoint_envelope": True,
        "representative_endpoint_phase_rows": reps,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "next_obligation": (
            "pair each endpoint/phase covariance envelope with the common one-complete-segment "
            "full-matrix selected-process lower; canonical P3 may move only if all phases and H/A exceed 1e-18"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_P2_V1_FULL_WORD_PARETO_HISTORY_QUOTIENT":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed", "exact_gap_labelled_successors_used",
        "certified_sufficient_quotient_used", "dominance_continuation_safe",
        "full_source_history_family_enumerated_at_sufficient_statistic_level",
        "same_history_upper_evaluated_before_endpoint_envelope",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "flat_800_node_ancestor_hull_used", "raw_tuner_cartesian_extrema_used",
        "P3_PROMOTED", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("history quotient lost frozen P2 V1 binding")
    target = d.get("global_word_target", {})
    if not (
        int(target.get("target_samples", 0)) > 0
        and float(target.get("history_duration_lower_s", 0.0))
        > float(target.get("global_covariance_word_upper_s", math.inf))
    ):
        f.append("full-word history target does not strictly cover covariance word")
    if int(d.get("terminal_history_classes", 0)) <= 0 or int(d.get("terminal_pareto_labels", 0)) <= 0:
        f.append("empty full-word Pareto history quotient")
    if int(d.get("endpoint_source_nodes", 0)) != 800:
        f.append("not every P2 source node appears as a full-word endpoint")
    digest = str(d.get("frontier_digest_sha256", ""))
    if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        f.append("invalid deterministic frontier digest")
    for row in d.get("representative_endpoint_phase_rows", []):
        if row.get("same_history_upper_evaluated_before_endpoint_envelope") is not True:
            f.append("representative endpoint envelope lost same-history ordering")
        if row.get("raw_tuner_cartesian_extrema_used") is not False:
            f.append("representative endpoint envelope used Cartesian tuner extrema")
        u = row.get("Sigma_translation_diagonal_upper_envelope", [])
        if len(u) != 4 or any(not (math.isfinite(float(x)) and float(x) > 0.0) for x in u):
            f.append("representative endpoint envelope invalid")
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
        "target": d["global_word_target"],
        "terminal_history_classes": d["terminal_history_classes"],
        "terminal_pareto_labels": d["terminal_pareto_labels"],
        "max_frontier_per_state": d["max_frontier_per_state"],
        "digest": d["frontier_digest_sha256"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
