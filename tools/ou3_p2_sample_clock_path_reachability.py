#!/usr/bin/env python3
"""Source-faithful sample-clock P2 reachability for the deployed OU-III tuner.

The older :mod:`ou3_p4_source_path_reachability` graph was written for a
0.1-second online-tune commit model.  The shipping wrapper no longer has that
schedule: at the beginning of *every* valid IMU sample it commits coefficients
staged by the previous sample, and later in the same sample ``update_tuner``
stages the next coefficients.  ``ADAPT_EVERY_SECS`` now belongs only to
posterior ``a_w`` covariance maintenance.

Consequently one tuner-state edge is one configured IMU interval, not an
arbitrarily late interval >= 0.1 s.  This producer reuses the existing audited
800-cell partition, target boxes, horizon bounds, binary32 literals, validated
EMA exponential, raw-tuner sigma state and filter-side sigma floor, but fixes
that timing semantics by evaluating every EMA image with

    min_elapsed = max_elapsed = dt.

The target boxes deliberately remain conservative (in particular R_S retains
its full deployed clamp until powf/sqrtf implementation error is enclosed).
No replay data or trajectory pruning enters.  The result is a replacement P2
source language for P3/P4; it does not itself promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_source_path_reachability as LEGACY

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _shipping_markers() -> list[str]:
    text = WRAPPER.read_text(encoding="utf-8")
    required = (
        "apply_pending_online_tune_();",
        "update_tuner(dt, a_vert_measurement, tuner_frequency_hz_());",
        "periodic_aw_cov_sync_tick_();",
        "constexpr float ADAPT_EVERY_SECS           = 0.1f;",
    )
    return [m for m in required if m not in text]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    dom = json.loads(domain_path.read_text(encoding="utf-8"))
    if dom.get("trajectory_fit") is not False:
        raise RuntimeError("sample-clock P2 domain must not be trajectory fitted")

    missing = _shipping_markers()
    if missing:
        raise RuntimeError(f"shipping online-tune order changed: {missing}")

    c = LEGACY._constants()
    dt = float(c["dt"])
    if not (math.isfinite(dt) and dt > 0.0):
        raise RuntimeError("configured IMU dt is not finite positive")

    tau_lo = max(c["min_tau"], LEGACY.down(c["tau_coeff"] * 0.5 / c["max_freq"]))
    tau = LEGACY._cells(tau_lo, c["max_tau"], 10)
    sigma_raw = LEGACY._cells(LEGACY.RAW_SIGMA_GRAPH_LOWER, c["max_sigma"], 8)
    rs = LEGACY._cells(c["min_RS"], c["max_RS"], 10)
    freq = LEGACY._cells(c["min_freq"], c["max_freq"], 8)

    states = []
    idx = {}
    for i, t in enumerate(tau):
        for j, s in enumerate(sigma_raw):
            for k, r in enumerate(rs):
                idx[(i, j, k)] = len(states)
                states.append((t, s, r))

    targets = []
    for f in freq:
        tt = LEGACY._tau_target(f, c)
        ht = LEGACY._tau_sigma_horizon(f, c)
        hr = LEGACY._rs_horizon(tt, c)
        rr = LEGACY._rs_target_box(c)
        for ss in sigma_raw:
            targets.append((tt, ss, rr, ht, hr))

    graph = [set() for _ in states]
    for q, (t, s, r) in enumerate(states):
        out = graph[q]
        for tt, ss, rr, ht, hr in targets:
            ti = LEGACY._matching(
                tau, LEGACY._ema_image(t, tt, ht, dt, max_elapsed=dt)
            )
            si = LEGACY._matching(
                sigma_raw, LEGACY._ema_image(s, ss, ht, dt, max_elapsed=dt)
            )
            ri = LEGACY._matching(
                rs, LEGACY._ema_image(r, rr, hr, dt, max_elapsed=dt)
            )
            for i in ti:
                for j in si:
                    for k in ri:
                        out.add(idx[(i, j, k)])

    gl = [sorted(x) for x in graph]
    comps = LEGACY._scc(gl)
    recurrent = set()
    for cc in comps:
        if len(cc) > 1 or (cc and cc[0] in graph[cc[0]]):
            recurrent.update(cc)

    # Keep the historical weak P3 corner visible so P3 can measure whether the
    # corrected source clock permits it to persist, rather than deleting it by
    # assumption.
    bad = []
    for q, (t, s_raw, r) in enumerate(states):
        x = (LEGACY.down(dt / t[1]), LEGACY.up(dt / t[0]))
        s_filter = LEGACY._filter_sigma_box(s_raw)
        if (LEGACY._overlap(s_filter, (0.05, 0.13025855423486765))
                and LEGACY._overlap(r, (149.21548743644342, 400.0))
                and LEGACY._overlap(x, (0.00041666665735344083, 0.0004837652693428343))):
            bad.append(q)
    bad_cycle = LEGACY._induced_cycle(gl, bad)
    bad_steps = LEGACY._longest_bad_residence(gl, bad)

    edge_count = sum(map(len, gl))
    full_cartesian = len(states) * len(states)
    failures = []
    if not edge_count:
        failures.append("sample-clock graph has no edges")
    if edge_count >= full_cartesian:
        failures.append("sample-clock graph collapsed to the full Cartesian transition relation")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P2_SAMPLE_CLOCK_SOURCE_PATH_REACHABILITY",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "shipping_online_tune_commit_at_sample_start": True,
        "shipping_online_tune_stage_after_current_sample": True,
        "aw_covariance_maintenance_cadence_not_tuner_commit_cadence": True,
        "configured_sample_dt_s": dt,
        "edge_elapsed_interval_s": [dt, dt],
        "arbitrary_late_online_tune_commit_modeled": False,
        "legacy_0p1s_online_tune_commit_model_used": False,
        "validated_exponential_used_for_ema": True,
        "source_float_literals_rounded_as_binary32": True,
        "raw_tuner_sigma_subfloor_states_included": True,
        "raw_tuner_sigma_partition_lower": LEGACY.RAW_SIGMA_GRAPH_LOWER,
        "filter_sigma_floor_mps2": LEGACY.FILTER_SIGMA_FLOOR,
        "filter_sigma_floor_separate_from_tuner_state": True,
        "RS_target_full_deployed_clamp_overapprox": True,
        "RS_target_powf_tightening_used": False,
        "RS_discrepancy_slew_horizon_covered": True,
        "partition": {
            "tau": len(tau),
            "sigma_tuner_raw": len(sigma_raw),
            "R_S": len(rs),
            "states": len(states),
            "target_boxes": len(targets),
        },
        "transition_edges": edge_count,
        "full_cartesian_transition_edges": full_cartesian,
        "transition_density": edge_count / full_cartesian,
        "strongly_connected_components": len(comps),
        "recurrent_states": len(recurrent),
        "old_worst_corner_state_count": len(bad),
        "old_worst_corner_states_in_any_recurrent_SCC": sum(q in recurrent for q in bad),
        "old_worst_corner_has_internal_recurrent_cycle": bad_cycle,
        "old_worst_corner_max_consecutive_sample_steps_upper": bad_steps,
        "old_worst_corner_max_residence_s_upper": (
            None if bad_steps is None else LEGACY.up(bad_steps * dt)
        ),
        # JSON-safe adjacency is intentionally emitted.  P3/P4 consumers need
        # the exact certified g->h language rather than rebuilding it with a
        # potentially different timing convention.
        "adjacency": gl,
        "path_graph_ready": not failures,
        "P2_SAMPLE_CLOCK_SOURCE_PATH_CERTIFICATE": "PASS" if not failures else "FAIL",
        "usable_P4_promoted": False,
        "next_obligation": (
            "attach the node-conditioned whole-word P3 information margin to this exact sample-clock graph; "
            "then accumulate directional Joseph credit and nonlinear sector defects only along certified paths"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_only",
        "shipping_online_tune_commit_at_sample_start",
        "shipping_online_tune_stage_after_current_sample",
        "aw_covariance_maintenance_cadence_not_tuner_commit_cadence",
        "validated_exponential_used_for_ema",
        "source_float_literals_rounded_as_binary32",
        "raw_tuner_sigma_subfloor_states_included",
        "filter_sigma_floor_separate_from_tuner_state",
        "RS_target_full_deployed_clamp_overapprox",
        "RS_discrepancy_slew_horizon_covered",
        "path_graph_ready",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "arbitrary_late_online_tune_commit_modeled",
        "legacy_0p1s_online_tune_commit_model_used",
        "RS_target_powf_tightening_used",
        "usable_P4_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    dt = float(d.get("configured_sample_dt_s", math.nan))
    elapsed = d.get("edge_elapsed_interval_s", [])
    if len(elapsed) != 2 or not all(float(x) == dt for x in elapsed):
        f.append("edge elapsed interval is not exactly one configured sample")
    states = int(d.get("partition", {}).get("states", 0))
    if states != 800:
        f.append("P2 source partition is not 800 states")
    edges = int(d.get("transition_edges", 0))
    full = int(d.get("full_cartesian_transition_edges", 0))
    if not (0 < edges < full == states * states):
        f.append("sample-clock transition relation is empty or Cartesian-complete")
    adjacency = d.get("adjacency", [])
    if len(adjacency) != states:
        f.append("adjacency does not contain every source node")
    elif any(any(not (0 <= int(w) < states) for w in ws) for ws in adjacency):
        f.append("adjacency contains an invalid node index")
    if d.get("P2_SAMPLE_CLOCK_SOURCE_PATH_CERTIFICATE") != "PASS":
        f.append("sample-clock P2 certificate did not pass")
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
        "status": d["P2_SAMPLE_CLOCK_SOURCE_PATH_CERTIFICATE"],
        "states": d["partition"]["states"],
        "edges": d["transition_edges"],
        "full_cartesian_edges": d["full_cartesian_transition_edges"],
        "transition_density": d["transition_density"],
        "scc": d["strongly_connected_components"],
        "recurrent_states": d["recurrent_states"],
        "old_worst_corner_cycle": d["old_worst_corner_has_internal_recurrent_cycle"],
        "old_worst_corner_max_residence_s": d["old_worst_corner_max_residence_s_upper"],
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
