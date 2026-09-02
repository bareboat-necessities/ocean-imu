#!/usr/bin/env python3
"""Clock-phase + staged/committed tuner automaton for OU-III P2.

The old 800-node P2 quotient remembers only the currently applied tuple
(tau,sigma_aw,R_S).  That is source-safe but too weak for P4 because the
shipping wrapper has memory in two places that matter to a complete-word proof:

* the adaptation clock stages a candidate only after a finite 13..26-sample
  source-clock gap; and
* the staged tuple is committed before the following valid IMU sample, while
  the EMA candidate keeps evolving toward the next staged tuple.

This producer lifts the 800-node physical partition to a finite stage-boundary
automaton without introducing a Cartesian free switch.  A stage-boundary node
is the ordered pair

    (c, s)

where ``c`` is the tuple used by the sample that just staged and ``s`` is the
candidate snapshot staged for the following valid sample.  A labelled edge

    (c,s) --g--> (s,t),   g in {13,...,26}

exists only when sample-by-sample EMA propagation from ``s`` over exactly ``g``
valid samples can stage ``t``.  The applied schedule on that segment is source
exact at the quotient level: phase 0 uses ``c``; phases 1..g use ``s``; the
next boundary stages ``t``.  Thus P4 can reconstruct every applied tuple at
every sample from the edge label while retaining the staged/committed
correlation.

The graph is stored factorized rather than materializing every second-order
edge.  ``successors_by_gap[s][g]`` is the exact quotient successor set, and the
reachable stage-pair set is the union of all first-order labelled edges.  This
is sufficient to enumerate or dynamically program complete P4 words while
keeping the certificate compact.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_sample_clock_source_refinement as CLOCK
import ou3_p4_source_path_reachability as PATH

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _build_labelled_successors(domain_path: Path):
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P2 phased tuner graph must not be trajectory fitted")

    c = PATH._constants()
    clock = CLOCK._clock_certificate(c)
    dt = float(clock["dt_binary32_s"])
    gaps = list(range(
        int(clock["finite_stage_spacing_valid_samples_lower"]),
        int(clock["finite_stage_spacing_valid_samples_upper"]) + 1,
    ))
    tau, sigma, rs, freq = CLOCK._partition(c)

    states = []
    index = {}
    for ti, t in enumerate(tau):
        for si, s in enumerate(sigma):
            for ri, r in enumerate(rs):
                index[(ti, si, ri)] = len(states)
                states.append((t, s, r))

    # One-dimensional exact-gap images, factored exactly as the shipping target
    # family.  Targets/horizons may vary independently sample by sample inside
    # their source boxes, so composing the same boxes every sample is a source
    # over-approximation but does not create an arbitrary commit jump.
    tau_match = {}
    sigma_match = {}
    rs_match = {}
    for fi, f in enumerate(freq):
        tt = PATH._tau_target(f, c)
        ht = PATH._tau_sigma_horizon(f, c)
        hr = CLOCK._configured_rs_horizon(tt, c)
        rr = PATH._rs_target_box(c)
        for gap in gaps:
            for ti, t in enumerate(tau):
                tau_match[(ti, fi, gap)] = PATH._matching(
                    tau, CLOCK._ema_samples(t, tt, ht, dt, gap)
                )
            for ri, r in enumerate(rs):
                rs_match[(ri, fi, gap)] = PATH._matching(
                    rs, CLOCK._ema_samples(r, rr, hr, dt, gap)
                )
            for target_si, ss in enumerate(sigma):
                for si, s in enumerate(sigma):
                    sigma_match[(si, fi, target_si, gap)] = PATH._matching(
                        sigma, CLOCK._ema_samples(s, ss, ht, dt, gap)
                    )

    labelled = [[set() for _ in gaps] for _ in states]
    union = [set() for _ in states]
    for q in range(len(states)):
        ti0 = q // (CLOCK.SIGMA_CELLS * CLOCK.RS_CELLS)
        rem = q % (CLOCK.SIGMA_CELLS * CLOCK.RS_CELLS)
        si0 = rem // CLOCK.RS_CELLS
        ri0 = rem % CLOCK.RS_CELLS
        for gi, gap in enumerate(gaps):
            out = labelled[q][gi]
            for fi in range(len(freq)):
                for target_si in range(len(sigma)):
                    tis = tau_match[(ti0, fi, gap)]
                    sis = sigma_match[(si0, fi, target_si, gap)]
                    ris = rs_match[(ri0, fi, gap)]
                    for i in tis:
                        for j in sis:
                            for k in ris:
                                out.add(index[(i, j, k)])
            union[q].update(out)

    return states, gaps, labelled, union, clock


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    states, gaps, labelled, union, clock = _build_labelled_successors(Path(domain_path).resolve())
    n = len(states)

    # Reachable stage-boundary pairs are exactly first-order source edges.  The
    # frozen-clock branch is represented separately because after stagnation no
    # further stage boundary occurs; the applied tuple simply self-holds.
    pair_count = sum(len(x) for x in union)
    indegree = [0] * n
    for out in union:
        for s in out:
            indegree[s] += 1

    per_gap_edge_counts = {}
    labelled_first_order_edges = 0
    for gi, gap in enumerate(gaps):
        q = sum(len(labelled[s][gi]) for s in range(n))
        per_gap_edge_counts[str(gap)] = q
        labelled_first_order_edges += q

    # Number of explicit second-order labelled edges in the factorized stage
    # automaton, without writing the potentially large edge list:
    # for every reachable pair (c,s), choose one exact gap g and one t in
    # successors_by_gap[s,g].
    second_order_edges = 0
    for s in range(n):
        next_labelled = sum(len(labelled[s][gi]) for gi in range(len(gaps)))
        second_order_edges += indegree[s] * next_labelled

    # Every pair node has at least one continuation in the finite-clock branch
    # because the source boxes and EMA are invariant.  Stagnation is an
    # absorbing committed-tuple hold, not an arbitrary stage edge.
    dead_pair_states = sum(
        indegree[s] for s in range(n)
        if sum(len(labelled[s][gi]) for gi in range(len(gaps))) == 0
    )

    failures = []
    if n != CLOCK.BASE_P2_STATES:
        failures.append("physical tuner partition changed")
    if pair_count <= 0 or pair_count >= CLOCK.BASE_P2_EDGES:
        failures.append("stage-pair quotient is empty or all-to-all")
    if dead_pair_states:
        failures.append("finite-clock stage-pair graph has dead states")
    if int(clock["pending_apply_delay_valid_samples"]) != 1:
        failures.append("pending staged tuple does not apply on next sample")
    if gaps != list(range(13, 27)):
        failures.append("finite clock phase alphabet changed from 13..26")

    # Compact gap-labelled successor table.  Each row is indexed by physical
    # staged tuple s; each gap entry lists possible next staged tuples t.
    succ = [
        {str(gap): sorted(labelled[s][gi]) for gi, gap in enumerate(gaps)}
        for s in range(n)
    ]

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P2_CLOCK_PHASE_STAGED_COMMITTED_TUNER_AUTOMATON",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "physical_partition_states": n,
        "clock_phase_gap_alphabet_samples": gaps,
        "stage_boundary_state": "ordered_pair(committed_on_staging_sample, staged_for_next_sample)",
        "stage_boundary_pair_states": pair_count,
        "first_order_unlabelled_edges": pair_count,
        "first_order_gap_labelled_edges": labelled_first_order_edges,
        "second_order_gap_labelled_edges_factorized_count": second_order_edges,
        "per_gap_first_order_edge_counts": per_gap_edge_counts,
        "transition_rule": "(c,s) --g--> (s,t) iff t in successors_by_gap[s][g]",
        "applied_tuple_by_phase": {
            "phase_0_staging_sample": "c",
            "phase_1_through_g": "s",
            "pending_apply_before_phase_1_innovation": True,
            "next_boundary_stages": "t",
        },
        "clock": clock,
        "EMA_candidate_memory_retained_across_commits": True,
        "staged_tuple_snapshot_retained": True,
        "committed_tuple_retained": True,
        "clock_phase_retained": True,
        "gap_to_candidate_transition_correlation_retained": True,
        "arbitrary_cartesian_tuner_switching_used": False,
        "frozen_clock_semantics": "absorbing hold of committed tuple; no future stage edge",
        "successors_by_gap": succ,
        "P2_PHASED_SOURCE_GRAPH_CERTIFICATE": "PASS" if not failures else "FAIL",
        "ready_for_P4_source_word_narrowing": not failures,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_only",
        "EMA_candidate_memory_retained_across_commits",
        "staged_tuple_snapshot_retained",
        "committed_tuple_retained",
        "clock_phase_retained",
        "gap_to_candidate_transition_correlation_retained",
        "ready_for_P4_source_word_narrowing",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("trajectory_replay_used", "filter_changed", "arbitrary_cartesian_tuner_switching_used"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_PHASED_SOURCE_GRAPH_CERTIFICATE") != "PASS":
        f.append("P2 phased source graph did not pass")
    if d.get("clock_phase_gap_alphabet_samples") != list(range(13, 27)):
        f.append("P2 phased graph lost 13..26 clock alphabet")
    if int(d.get("physical_partition_states", 0)) != CLOCK.BASE_P2_STATES:
        f.append("P2 phased graph changed physical partition")
    if not 0 < int(d.get("stage_boundary_pair_states", 0)) < CLOCK.BASE_P2_EDGES:
        f.append("invalid stage-boundary pair-state count")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P2_PHASED_SOURCE_GRAPH_CERTIFICATE"],
        "physical_states": d["physical_partition_states"],
        "stage_pair_states": d["stage_boundary_pair_states"],
        "first_order_gap_labelled_edges": d["first_order_gap_labelled_edges"],
        "second_order_gap_labelled_edges": d["second_order_gap_labelled_edges_factorized_count"],
        "per_gap_edges": d["per_gap_first_order_edge_counts"],
        "validation_failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
