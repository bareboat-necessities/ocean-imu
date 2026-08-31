#!/usr/bin/env python3
"""Source-dynamic reachability backend for the OU-III P2/P4 path certificate.

This graph is a conservative source-language over-approximation of the deployed
online tuner.  It exists to prevent later P3/P4/P5 certificates from selecting
independent worst-case ``tau``, ``sigma_aw`` and ``R_S`` values on every sample
when the shipping code can only move those states through its EMA dynamics.

Two details are intentionally fail-closed here:

* ``tune_.sigma_applied`` is a *raw tuner state*.  Once the variance estimate is
  ready it can fall below the 0.05 m/s^2 floor applied when the MEKF stationary
  OU standard deviation is committed.  The graph therefore contains sub-floor
  tuner states and records the separate filter-side 0.05 floor.
* the SpectralMSE ``R_S`` target contains implementation ``powf``/``sqrtf``
  operations.  Until their implementation error is independently enclosed, the
  path graph uses the full deployed ``[MIN_R_S, MAX_R_S]`` target clamp.  This
  loses some R_S correlation but cannot omit a real source transition.

EMA images use validated exponential enclosures.  A graph edge represents any
inter-commit delay >= the deployed commit threshold; arbitrary late commits are
covered by allowing the total decay factor down to zero.  The R_S horizon also
covers the shipping discrepancy/slew acceleration by using the full admissible
horizon interval rather than assuming the nominal ``adapt_RS_mult*tau`` value.

The result is a P2 path-language certificate and an input to P4/P5.  It does not
by itself promote a nonlinear funnel.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import struct
from pathlib import Path

from ou3_interval import Interval
import ou3_source_reachable_matrix_p3 as P3
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
LIMITS = REPO / "src" / "tuner" / "SeaStateAdaptationLimits.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2

# Conservative raw-tuner graph floor.  Shipping var_wave is floored at 1e-6
# before sqrt and the deployed sigma coefficient is positive, so 1e-6 is well
# below every physical target while still permitting geometric partitioning.
RAW_SIGMA_GRAPH_LOWER = 1.0e-6
# The shipping MEKF commit floor is a different object from the raw tuner state.
FILTER_SIGMA_FLOOR = 0.05


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _f32(x: float) -> float:
    return struct.unpack("!f", struct.pack("!f", float(x)))[0]


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _literal_member(text: str, name: str) -> float:
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f", text)
    if not m:
        raise RuntimeError(f"cannot extract deployed literal member {name}")
    return _f32(float(m.group(1)))


def _const(text: str, name: str) -> float:
    return _f32(float(SOURCE.parse_const(text, name)))


def _constants() -> dict:
    text = WRAPPER.read_text(encoding="utf-8")
    lim = LIMITS.read_text(encoding="utf-8")
    if "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;" not in text:
        raise RuntimeError("path backend requires deployed SpectralMSE law")
    return {
        "dt": _const(text, "FREQ_SMOOTHER_DT"),
        "commit": _const(text, "ADAPT_EVERY_SECS"),
        "tau_coeff": _literal_member(text, "tau_coeff_"),
        "sigma_coeff": _literal_member(text, "sigma_coeff_"),
        "adapt_tau_sea_periods": _const(text, "ADAPT_TAU_SEA_PERIODS"),
        "adapt_RS_mult": _const(text, "ADAPT_RS_MULT"),
        "min_tau": _const(text, "MIN_TAU_S"),
        "max_tau": _const(text, "MAX_TAU_S"),
        "max_sigma": _const(text, "MAX_SIGMA_A"),
        "min_RS": _const(text, "MIN_R_S"),
        "max_RS": _const(text, "MAX_R_S"),
        "min_freq": _const(text, "MIN_TUNE_FREQ_HZ"),
        "max_freq": _const(text, "MAX_TUNE_FREQ_HZ"),
        "horizon_min": _const(lim, "kDynamicEmaHorizonMinSec"),
        "horizon_max": _const(lim, "kDynamicEmaHorizonMaxSec"),
        "time_scale_min": _const(lim, "kDynamicEmaTimeScaleMinSec"),
        "time_scale_max": _const(lim, "kDynamicEmaTimeScaleMaxSec"),
    }


def _iv(pair) -> Interval:
    return Interval(float(pair[0]), float(pair[1]))


def _clamp_interval(x: Interval, lo: float, hi: float) -> Interval:
    if not (lo <= hi):
        raise RuntimeError("invalid clamp")
    return Interval(
        down(max(lo, min(hi, x.lo))),
        up(max(lo, min(hi, x.hi))),
    )


def _overlap(a, b) -> bool:
    return not (a[1] < b[0] or b[1] < a[0])


def _cells(lo: float, hi: float, n: int):
    if not (0.0 < lo < hi):
        raise RuntimeError(f"invalid geometric partition [{lo}, {hi}]")
    e = P3.geom_edges(float(lo), float(hi), int(n))
    return [(down(e[i]), up(e[i + 1])) for i in range(len(e) - 1)]


def _matching(cells, image):
    return [i for i, cell in enumerate(cells) if _overlap(cell, image)]


def _tau_target(freq, c):
    raw = I(c["tau_coeff"]) * I(0.5) / _iv(freq)
    return _clamp_interval(raw, c["min_tau"], c["max_tau"]).as_list()


def _tau_sigma_horizon(freq, c):
    # Shipping uses sea_time = 0.5/f_tune, independently clamps that dynamic
    # time scale, then multiplies by adapt_tau_sea_periods and clamps the final
    # horizon.  Do not infer it back from a possibly clamped tau_target.
    sea = I(0.5) / _iv(freq)
    sea = _clamp_interval(sea, c["time_scale_min"], c["time_scale_max"])
    h = I(c["adapt_tau_sea_periods"]) * sea
    return _clamp_interval(h, c["horizon_min"], c["horizon_max"]).as_list()


def _rs_horizon(tau_target, c):
    # adaptiveSmoothingHorizonSec starts from mult*clamped(tau) and may only
    # shorten it through the discrepancy gate before the final horizon clamp.
    # Hence [horizon_min, nominal_upper] covers every shipping slew branch.
    safe_tau = _clamp_interval(
        _iv(tau_target), c["time_scale_min"], c["time_scale_max"]
    )
    nominal = I(c["adapt_RS_mult"]) * safe_tau
    upper = min(c["horizon_max"], nominal.hi)
    return (down(c["horizon_min"]), up(max(c["horizon_min"], upper)))


def _rs_target_box(c):
    # Deliberately broad until the implementation powf/sqrtf error is itself
    # source-qualified.  The shipping clamp proves this box regardless of the
    # transcendental implementation.
    return (down(c["min_RS"]), up(c["max_RS"]))


def _ema_image(x, target, horizon, min_elapsed):
    """Conservative image after any elapsed time >= ``min_elapsed``.

    For a positive first-order EMA, the total old-state weight is in
    [0, exp(-min_elapsed/h_max)].  Allowing zero covers an arbitrarily late
    commit.  Variable targets inside one target box remain inside its convex
    hull, so the interval affine image is source-complete.
    """
    h = _iv(horizon)
    if not h.lo > 0.0:
        raise RuntimeError("invalid EMA horizon")
    exponent = -(I(min_elapsed) / I(h.hi))
    if exponent.lo < -VT.MAX_ABS_ARGUMENT or exponent.hi > 0.0:
        raise RuntimeError("EMA exponent left validated exponential range")
    a_hi = VT.exp_interval(exponent).hi
    a = Interval(0.0, up(a_hi))
    one_minus_a = I(1.0) - a
    image = a * _iv(x) + one_minus_a * _iv(target)
    return image.as_list()


def _filter_sigma_box(raw_sigma):
    s = _iv(raw_sigma)
    return (
        down(max(FILTER_SIGMA_FLOOR, s.lo)),
        up(max(FILTER_SIGMA_FLOOR, s.hi)),
    )


def _scc(graph):
    n = len(graph)
    seen = [False] * n
    order = []

    def dfs(v):
        seen[v] = True
        for w in graph[v]:
            if not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)
    rg = [[] for _ in range(n)]
    for v, ws in enumerate(graph):
        for w in ws:
            rg[w].append(v)
    seen = [False] * n
    comps = []

    def rdfs(v, acc):
        seen[v] = True
        acc.append(v)
        for w in rg[v]:
            if not seen[w]:
                rdfs(w, acc)

    for v in reversed(order):
        if not seen[v]:
            acc = []
            rdfs(v, acc)
            comps.append(acc)
    return comps


def _induced_cycle(graph, nodes):
    ns = set(nodes)
    sub = {v: [w for w in graph[v] if w in ns] for v in ns}
    visiting = set()
    done = set()

    def dfs(v):
        if v in visiting:
            return True
        if v in done:
            return False
        visiting.add(v)
        for w in sub[v]:
            if dfs(w):
                return True
        visiting.remove(v)
        done.add(v)
        return False

    return any(dfs(v) for v in list(ns) if v not in done)


def _longest_bad_residence(graph, bad):
    ns = set(bad)
    sub = {v: [w for w in graph[v] if w in ns] for v in ns}
    if _induced_cycle(graph, bad):
        return None
    memo = {}

    def longest(v):
        if v not in memo:
            memo[v] = 1 + max((longest(w) for w in sub[v]), default=0)
        return memo[v]

    return max((longest(v) for v in ns), default=0)


def build(domain_path=DEFAULT_DOMAIN):
    dom = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    if dom.get("trajectory_fit") is not False:
        raise RuntimeError("path domain must not be trajectory fitted")
    c = _constants()

    tau_lo = max(c["min_tau"], down(c["tau_coeff"] * 0.5 / c["max_freq"]))
    tau = _cells(tau_lo, c["max_tau"], 10)
    sigma_raw = _cells(RAW_SIGMA_GRAPH_LOWER, c["max_sigma"], 8)
    rs = _cells(c["min_RS"], c["max_RS"], 10)
    freq = _cells(c["min_freq"], c["max_freq"], 8)

    states = []
    idx = {}
    for i, t in enumerate(tau):
        for j, s in enumerate(sigma_raw):
            for k, r in enumerate(rs):
                idx[(i, j, k)] = len(states)
                states.append((t, s, r))

    targets = []
    for f in freq:
        tt = _tau_target(f, c)
        ht = _tau_sigma_horizon(f, c)
        hr = _rs_horizon(tt, c)
        rr = _rs_target_box(c)
        for ss in sigma_raw:
            targets.append((tt, ss, rr, ht, hr))

    graph = [set() for _ in states]
    min_elapsed = c["commit"]
    for q, (t, s, r) in enumerate(states):
        out = graph[q]
        for tt, ss, rr, ht, hr in targets:
            ti = _matching(tau, _ema_image(t, tt, ht, min_elapsed))
            si = _matching(sigma_raw, _ema_image(s, ss, ht, min_elapsed))
            ri = _matching(rs, _ema_image(r, rr, hr, min_elapsed))
            for i in ti:
                for j in si:
                    for k in ri:
                        out.add(idx[(i, j, k)])

    gl = [sorted(x) for x in graph]
    comps = _scc(gl)
    recurrent = set()
    for cc in comps:
        if len(cc) > 1 or (cc and cc[0] in graph[cc[0]]):
            recurrent.update(cc)

    # Historic weak P3 cell remains explicit.  Compare against the *filter*
    # sigma box, not the raw tuner sigma state, because P3 sees the 0.05 floor.
    bad = []
    for q, (t, s_raw, r) in enumerate(states):
        x = (down(c["dt"] / t[1]), up(c["dt"] / t[0]))
        s_filter = _filter_sigma_box(s_raw)
        if (_overlap(s_filter, (0.05, 0.13025855423486765))
                and _overlap(r, (149.21548743644342, 400.0))
                and _overlap(x, (0.00041666665735344083, 0.0004837652693428343))):
            bad.append(q)
    bad_cycle = _induced_cycle(gl, bad)
    steps = _longest_bad_residence(gl, bad)

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P2_SOURCE_DYNAMIC_PATH_REACHABILITY",
        "source_only": True,
        "trajectory_replay_used": False,
        "deployed_default_law": "SpectralMSE",
        "source_float_literals_rounded_as_binary32": True,
        "validated_exponential_used_for_ema": True,
        "arbitrary_late_commit_overapproximated": True,
        "inter_commit_elapsed_upper_assumed_s": None,
        "raw_tuner_sigma_subfloor_states_included": True,
        "raw_tuner_sigma_partition_lower": RAW_SIGMA_GRAPH_LOWER,
        "filter_sigma_floor_mps2": FILTER_SIGMA_FLOOR,
        "filter_sigma_floor_separate_from_tuner_state": True,
        "RS_target_full_deployed_clamp_overapprox": True,
        "RS_target_powf_tightening_used": False,
        "RS_discrepancy_slew_horizon_covered": True,
        "commit_period_s": min_elapsed,
        "partition": {
            "tau": len(tau),
            "sigma_tuner_raw": len(sigma_raw),
            "sigma": len(sigma_raw),
            "R_S": len(rs),
            "states": len(states),
            "target_boxes": len(targets),
        },
        "transition_edges": sum(map(len, gl)),
        "strongly_connected_components": len(comps),
        "recurrent_states": len(recurrent),
        "old_worst_corner_state_count": len(bad),
        "old_worst_corner_states_in_any_recurrent_SCC": sum(q in recurrent for q in bad),
        "old_worst_corner_has_internal_recurrent_cycle": bad_cycle,
        "old_worst_corner_max_consecutive_commit_steps_upper": steps,
        "old_worst_corner_max_residence_s_upper": None if steps is None else up(steps * min_elapsed),
        "path_graph_ready": True,
        "P2_SOURCE_PATH_CERTIFICATE": "PASS",
        "usable_P4_promoted": False,
        "next_obligation": (
            "propagate complete-word Phi/Omega and the exact nonlinear return map on this source-complete graph; "
            "R_S path tightening may be added only after powf/sqrtf implementation error is independently enclosed"
        ),
        "failures": [],
    }


def validate(d):
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_only") is not True or d.get("trajectory_replay_used") is not False:
        failures.append("path graph is not source-only")
    if d.get("source_float_literals_rounded_as_binary32") is not True:
        failures.append("source float literals are not modeled as binary32")
    if d.get("validated_exponential_used_for_ema") is not True:
        failures.append("EMA graph does not use validated exponential bounds")
    if d.get("arbitrary_late_commit_overapproximated") is not True:
        failures.append("graph assumes an unsupported upper inter-commit delay")
    if d.get("inter_commit_elapsed_upper_assumed_s") is not None:
        failures.append("graph introduced an unsupported upper inter-commit delay")
    if d.get("raw_tuner_sigma_subfloor_states_included") is not True:
        failures.append("sub-floor tuner sigma states omitted")
    if not float(d.get("raw_tuner_sigma_partition_lower", math.inf)) < FILTER_SIGMA_FLOOR:
        failures.append("raw tuner sigma partition does not extend below filter floor")
    if d.get("filter_sigma_floor_separate_from_tuner_state") is not True:
        failures.append("tuner/filter sigma states were conflated")
    if d.get("RS_target_full_deployed_clamp_overapprox") is not True:
        failures.append("R_S target is not conservatively source-complete")
    if d.get("RS_target_powf_tightening_used") is not False:
        failures.append("unqualified powf tightening entered path graph")
    if d.get("RS_discrepancy_slew_horizon_covered") is not True:
        failures.append("R_S discrepancy/slew horizon branch omitted")
    if d.get("path_graph_ready") is not True:
        failures.append("path graph not ready")
    if d.get("P2_SOURCE_PATH_CERTIFICATE") != "PASS":
        failures.append("P2 source path certificate did not pass")
    if int(d.get("partition", {}).get("states", 0)) <= 0 or int(d.get("transition_edges", 0)) <= 0:
        failures.append("empty path graph")
    if d.get("usable_P4_promoted") is not False:
        failures.append("reachability prematurely promoted P4")
    if int(d.get("old_worst_corner_state_count", 0)) <= 0:
        failures.append("old worst corner not represented")
    return list(dict.fromkeys(failures))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        k: d[k] for k in (
            "partition", "transition_edges", "strongly_connected_components",
            "recurrent_states", "old_worst_corner_state_count",
            "old_worst_corner_states_in_any_recurrent_SCC",
            "old_worst_corner_has_internal_recurrent_cycle",
            "old_worst_corner_max_residence_s_upper",
            "P2_SOURCE_PATH_CERTIFICATE", "failures",
        )
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
