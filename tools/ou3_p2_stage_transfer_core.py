#!/usr/bin/env python3
"""Frozen P2 -> P3 correlation/path-memory interface for deployed OU-III.

P2's 800 physical tuner cells are a useful source partition, but the 800-node
endpoint quotient is not a sufficient interface to P3: over a long word its
ancestor relation can become complete, after which independently taking extrema
of tau, sigma_aw and R_S manufactures parameter combinations that no tuner
history can realize.

This module freezes the stronger interface that downstream proof stages must
consume.  It does NOT replace the existing P2 source/timing mathematics.  It
packages the already-certified clock-phase staged/committed automaton as a
factorized transfer system with source-correlated segment kernels.

A stage-boundary state is the ordered pair (c,s):

* c is the physical tuner cell applied on the staging/boundary sample;
* s is the staged cell committed before the following valid sample.

A legal edge is

    (c,s) --g--> (s,t),  g in {13,...,26},

where t is in the exact P2 successor set for (s,g).  The following g applied
samples all use the *same physical source cell s*.  Therefore every quantity in
one segment kernel is evaluated from one common (tau,sigma,R_S) cell; consumers
must not select tau from one node, sigma from another and R_S from a third.

The exported sufficient statistics are intentionally source-only and useful to
multiple later proof stages:

* integral bounds for lambda=1/tau over the segment;
* q_c=2 sigma^2/tau and its time mass;
* sigma^2, R_S^-2, and normalized S-measurement strength
  (sigma*h^3/R_S)^2;
* the complete tuple cell itself and exact gap/successor relation.

P3 may project this interface for a scalar theorem quantity only when it proves
that the projection can only add paths (for example the existing tau-decay
budget).  Correlated covariance/information propagation must retain the pair
state, or use a separately certified sufficient quotient.  The old 800-node
ancestor hull is explicitly forbidden as a correlated P3 interface.

No replay values, trajectory fitting, theorem-domain shrink, filter change, or
powf/sqrtf target tightening is introduced here.
"""
from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path

import ou3_p2_clock_phase_tuner_graph as GRAPH
import ou3_p4_source_node_cells as NODES

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
INTERFACE_VERSION = "OU3_P2_CORRELATED_STAGE_TRANSFER_V1"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _positive_box(values, name: str) -> tuple[float, float]:
    lo, hi = map(float, values)
    if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 < lo <= hi):
        raise RuntimeError(f"{name} lost positive finite interval")
    return lo, hi


def _stats(node: dict, gap: int, dt: float) -> dict:
    tau_lo, tau_hi = _positive_box(node["tau_s"], "tau")
    sig_lo, sig_hi = _positive_box(node["sigma_filter_committed_mps2"], "sigma")
    rs_lo, rs_hi = _positive_box(node["R_S_filter_std"], "R_S")
    duration_lo = down(int(gap) * dt)
    duration_hi = up(int(gap) * dt)

    lam_lo = down(1.0 / tau_hi)
    lam_hi = up(1.0 / tau_lo)
    qc_lo = down(2.0 * sig_lo * sig_lo / tau_hi)
    qc_hi = up(2.0 * sig_hi * sig_hi / tau_lo)
    sigma2_lo = down(sig_lo * sig_lo)
    sigma2_hi = up(sig_hi * sig_hi)
    inv_rs2_lo = down(1.0 / (rs_hi * rs_hi))
    inv_rs2_hi = up(1.0 / (rs_lo * rs_lo))

    h3 = dt * dt * dt
    s_info_lo = down((sig_lo * h3 / rs_hi) ** 2)
    s_info_hi = up((sig_hi * h3 / rs_lo) ** 2)

    return {
        "applied_source_node": int(node["index"]),
        "samples": int(gap),
        "duration_s": [duration_lo, duration_hi],
        "tau_s": [tau_lo, tau_hi],
        "sigma_filter_committed_mps2": [sig_lo, sig_hi],
        "R_S_filter_std": [rs_lo, rs_hi],
        "pseudo_update_period_s": list(map(float, node["pseudo_update_period_s"])),
        "lambda_inv_tau_per_s": [lam_lo, lam_hi],
        "decay_exponent_integral": [down(duration_lo * lam_lo), up(duration_hi * lam_hi)],
        "q_c_m2ps5": [qc_lo, qc_hi],
        "q_c_time_mass": [down(duration_lo * qc_lo), up(duration_hi * qc_hi)],
        "sigma_squared": [sigma2_lo, sigma2_hi],
        "inverse_R_S_variance": [inv_rs2_lo, inv_rs2_hi],
        "normalized_S_information_per_packet": [s_info_lo, s_info_hi],
        "tau_sigma_R_S_from_same_physical_cell": True,
        "independent_coordinate_extrema_used": False,
    }


@functools.lru_cache(maxsize=4)
def runtime(domain_path: Path = DEFAULT_DOMAIN):
    """Return the non-JSON factorized correlation interface for proof consumers."""
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P2 correlation interface must not be trajectory fitted")

    states, gaps, labelled, union, clock = GRAPH._build_labelled_successors(path)
    nodes_payload = NODES.build()
    nf = NODES.validate(nodes_payload)
    if nf:
        raise RuntimeError(f"P2 source-node materialization failed: {nf}")
    nodes = nodes_payload["nodes"]
    if len(states) != len(nodes) or len(states) != NODES.EXPECTED_STATES:
        raise RuntimeError("P2 correlation interface state cardinality mismatch")
    if gaps != list(range(13, 27)):
        raise RuntimeError("P2 correlation interface clock alphabet changed")

    # Verify that the graph's nested-loop tuple ordering and the materialized
    # node ordering are exactly the same partition before binding them.
    for i, ((tau, sigma_raw, rs), node) in enumerate(zip(states, nodes)):
        if int(node["index"]) != i:
            raise RuntimeError("P2 node index ordering changed")
        nt = tuple(map(float, node["tau_s"]))
        nr = tuple(map(float, node["R_S_filter_std"]))
        if not (float(tau[0]) == nt[0] and float(tau[1]) == nt[1]):
            raise RuntimeError(f"P2 tau partition mismatch at node {i}")
        if not (float(rs[0]) == nr[0] and float(rs[1]) == nr[1]):
            raise RuntimeError(f"P2 R_S partition mismatch at node {i}")
        # The graph state stores raw tuner sigma; node materialization keeps both
        # raw and filter-side clamped sigma.  Match the raw state here.
        ns = tuple(map(float, node["sigma_tuner_raw_mps2"]))
        if not (float(sigma_raw[0]) == ns[0] and float(sigma_raw[1]) == ns[1]):
            raise RuntimeError(f"P2 raw sigma partition mismatch at node {i}")

    dt = float(clock["dt_binary32_s"])
    kernels = [[_stats(nodes[s], gap, dt) for gap in gaps] for s in range(len(nodes))]
    return {
        "path": path,
        "nodes": nodes,
        "gaps": gaps,
        "labelled_successors": labelled,
        "union_successors": union,
        "clock": clock,
        "segment_kernels": kernels,
        "dt": dt,
    }


def legal_pair(c: int, s: int, rt=None) -> bool:
    r = runtime() if rt is None else rt
    c, s = int(c), int(s)
    return 0 <= c < len(r["nodes"]) and s in r["union_successors"][c]


def successors(staged: int, gap: int, rt=None) -> tuple[int, ...]:
    r = runtime() if rt is None else rt
    s = int(staged)
    try:
        gi = r["gaps"].index(int(gap))
    except ValueError as exc:
        raise ValueError("gap outside certified 13..26 alphabet") from exc
    if not 0 <= s < len(r["nodes"]):
        raise IndexError("staged source node outside P2 partition")
    return tuple(sorted(int(x) for x in r["labelled_successors"][s][gi]))


def segment_kernel(staged: int, gap: int, rt=None) -> dict:
    r = runtime() if rt is None else rt
    s = int(staged)
    try:
        gi = r["gaps"].index(int(gap))
    except ValueError as exc:
        raise ValueError("gap outside certified 13..26 alphabet") from exc
    if not 0 <= s < len(r["nodes"]):
        raise IndexError("staged source node outside P2 partition")
    return dict(r["segment_kernels"][s][gi])


def transition(c: int, s: int, gap: int, t: int, rt=None) -> dict:
    """Return one certified pair-state transition and its correlated segment."""
    r = runtime() if rt is None else rt
    c, s, t, gap = int(c), int(s), int(t), int(gap)
    if not legal_pair(c, s, r):
        raise ValueError("start (committed,staged) pair is not P2 reachable")
    if t not in successors(s, gap, r):
        raise ValueError("destination staged node is not a legal gap-labelled successor")
    return {
        "start_pair": [c, s],
        "gap_samples": gap,
        "end_pair": [s, t],
        "boundary_sample_applied_node": c,
        "following_segment_applied_node": s,
        "following_segment": segment_kernel(s, gap, r),
        "pair_shift_exact": True,
        "same_staged_node_becomes_next_committed": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    r = runtime(Path(domain_path).resolve())
    nodes = r["nodes"]
    gaps = r["gaps"]
    labelled = r["labelled_successors"]
    union = r["union_successors"]

    pair_count = sum(len(x) for x in union)
    labelled_edges = sum(
        len(labelled[s][gi])
        for s in range(len(nodes))
        for gi in range(len(gaps))
    )
    second_order = 0
    indegree = [0] * len(nodes)
    for out in union:
        for s in out:
            indegree[int(s)] += 1
    for s in range(len(nodes)):
        second_order += indegree[s] * sum(len(labelled[s][gi]) for gi in range(len(gaps)))

    segment_rows = []
    for s in range(len(nodes)):
        row = {"staged_source_node": s, "by_gap": {}}
        for gi, gap in enumerate(gaps):
            row["by_gap"][str(gap)] = {
                "kernel": r["segment_kernels"][s][gi],
                "next_staged_nodes": sorted(int(x) for x in labelled[s][gi]),
            }
        segment_rows.append(row)

    failures = []
    if len(nodes) != NODES.EXPECTED_STATES:
        failures.append("physical P2 source partition changed")
    if pair_count <= 0:
        failures.append("no reachable P2 stage-boundary pairs")
    if labelled_edges <= 0 or second_order <= 0:
        failures.append("P2 correlation transition system is empty")
    if any(not labelled[s][gi] for s in range(len(nodes)) for gi in range(len(gaps))):
        failures.append("P2 correlation interface has a dead finite-stage kernel")

    return {
        "schema": SCHEMA,
        "interface_version": INTERFACE_VERSION,
        "qualification": "OU3_P2_SOURCE_CORRELATION_PATH_MEMORY_INTERFACE",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "powf_sqrtf_target_tightening_used": False,
        "physical_source_states": len(nodes),
        "stage_boundary_state": "ordered_pair(committed_on_boundary_sample,staged_for_following_segment)",
        "stage_boundary_pair_states": pair_count,
        "clock_gap_alphabet_samples": gaps,
        "gap_labelled_first_order_edges": labelled_edges,
        "gap_labelled_pair_state_edges_factorized_count": second_order,
        "transition_rule": "(c,s)--g-->(s,t) iff t in next_staged_nodes[s,g]",
        "boundary_sample_uses_committed_c": True,
        "following_g_samples_use_staged_s": True,
        "same_staged_s_becomes_next_committed": True,
        "tau_sigma_R_S_joint_cell_retained_per_segment": True,
        "EMA_stage_commit_history_retained": True,
        "clock_phase_retained": True,
        "arbitrary_cartesian_tuner_switching_used": False,
        "old_800_node_ancestor_hull_allowed_for_correlated_P3": False,
        "marginal_projection_requires_monotone_adds_paths_proof": True,
        "P3_consumer_must_retain_pair_state_or_certified_sufficient_quotient": True,
        "segment_sufficient_statistics": [
            "tau_s", "sigma_filter_committed_mps2", "R_S_filter_std",
            "lambda_inv_tau_per_s", "decay_exponent_integral", "q_c_m2ps5",
            "q_c_time_mass", "sigma_squared", "inverse_R_S_variance",
            "normalized_S_information_per_packet",
        ],
        "segment_transfer_table": segment_rows,
        "consumer_contract": {
            "correlated_quantities_must_come_from_same_segment_node": True,
            "legal_word_must_follow_pair_shift_transition_rule": True,
            "independent_tau_sigma_R_S_extremization_before_propagation": "FORBIDDEN",
            "global_800_ancestor_hull_as_P3_covariance_information_input": "FORBIDDEN",
            "scalar_projection_exception": "allowed only with a proof that projection can only add source paths",
            "freeze_policy": "changes to this interface version or semantics require explicit P2 certificate/test update before P3/P4/P5 consumers may change",
        },
        "P2_TIMING_SOURCE_MATHEMATICS_RETAINED": True,
        "P2_CORRELATION_INTERFACE_CERTIFICATE": "PASS" if not failures else "FAIL",
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "next_obligation": (
            "make P3 covariance/information propagation consume this pair-state transfer interface; only certified scalar projections may discard coordinates, and no correlated P3 bound may return to the flat 800-node ancestor hull"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("interface_version") != INTERFACE_VERSION:
        f.append("P2 correlation interface version changed")
    if d.get("qualification") != "OU3_P2_SOURCE_CORRELATION_PATH_MEMORY_INTERFACE":
        f.append("wrong qualification")
    for key in (
        "source_only", "boundary_sample_uses_committed_c",
        "following_g_samples_use_staged_s", "same_staged_s_becomes_next_committed",
        "tau_sigma_R_S_joint_cell_retained_per_segment", "EMA_stage_commit_history_retained",
        "clock_phase_retained", "marginal_projection_requires_monotone_adds_paths_proof",
        "P3_consumer_must_retain_pair_state_or_certified_sufficient_quotient",
        "P2_TIMING_SOURCE_MATHEMATICS_RETAINED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "powf_sqrtf_target_tightening_used", "arbitrary_cartesian_tuner_switching_used",
        "old_800_node_ancestor_hull_allowed_for_correlated_P3",
        "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("physical_source_states", 0)) != NODES.EXPECTED_STATES:
        f.append("P2 physical source state count changed")
    if d.get("clock_gap_alphabet_samples") != list(range(13, 27)):
        f.append("P2 correlation clock alphabet changed")
    if int(d.get("stage_boundary_pair_states", 0)) <= 0:
        f.append("P2 correlation pair-state set is empty")
    if int(d.get("gap_labelled_first_order_edges", 0)) <= 0:
        f.append("P2 correlation first-order relation is empty")
    if int(d.get("gap_labelled_pair_state_edges_factorized_count", 0)) <= 0:
        f.append("P2 correlation pair-state relation is empty")
    c = d.get("consumer_contract", {})
    if c.get("correlated_quantities_must_come_from_same_segment_node") is not True:
        f.append("P2 consumer contract lost same-node correlation")
    if c.get("legal_word_must_follow_pair_shift_transition_rule") is not True:
        f.append("P2 consumer contract lost pair-shift rule")
    if c.get("independent_tau_sigma_R_S_extremization_before_propagation") != "FORBIDDEN":
        f.append("P2 consumer contract permits independent source extrema")
    if c.get("global_800_ancestor_hull_as_P3_covariance_information_input") != "FORBIDDEN":
        f.append("P2 consumer contract permits flat ancestor hull")
    if d.get("P2_CORRELATION_INTERFACE_CERTIFICATE") != "PASS":
        f.append("P2 correlation interface did not pass")
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
        "interface_version": d["interface_version"],
        "status": d["P2_CORRELATION_INTERFACE_CERTIFICATE"],
        "physical_source_states": d["physical_source_states"],
        "stage_boundary_pair_states": d["stage_boundary_pair_states"],
        "gap_labelled_first_order_edges": d["gap_labelled_first_order_edges"],
        "gap_labelled_pair_state_edges_factorized_count": d["gap_labelled_pair_state_edges_factorized_count"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
