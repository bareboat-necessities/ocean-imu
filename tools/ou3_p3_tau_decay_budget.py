#!/usr/bin/env python3
"""Clock-phase source bound for the LTV OU decay integral used by P3.

A source-uniform LTV controllability bound needs exp(-int lambda dt), with
lambda=1/tau.  Replacing that integral by H/tau_min is valid but unnecessarily
forgets the deployed tuner slew.  The complete 800-node source graph is not
needed for this scalar: lambda depends only on the ten-cell tau projection.

This producer rebuilds exactly that projection from the retained 13..26-sample
staging model.  For each endpoint tau cell and requested window length N IMU
samples it computes an upper bound on

    integral_window dt / tau(t).

The DP uses complete stage segments.  An arbitrary physical window may start and
end inside segments, so it is covered by extending to the preceding and next
stage boundary.  The resulting full-segment cover contains at most two extra
maximum-gap segments.  Maximizing over all compatible projected paths therefore
upper-bounds the actual integral.  The infinite-time floating-clock stagnation
branch is included separately as an exact constant-tau hold candidate.

Projecting away sigma and R_S only adds tau paths; it cannot remove a real one.
This is source evidence, not a P3 promotion by itself.
"""
from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path

import ou3_p4_sample_clock_source_refinement as SAMPLE
import ou3_p4_source_path_reachability as PATH

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


@functools.lru_cache(maxsize=4)
def _tau_projection(domain_path: Path = DEFAULT_DOMAIN):
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("tau decay budget must not be trajectory fitted")
    c = PATH._constants()
    clock = SAMPLE._clock_certificate(c)
    dt = float(clock["dt_binary32_s"])
    gaps = list(range(
        int(clock["finite_stage_spacing_valid_samples_lower"]),
        int(clock["finite_stage_spacing_valid_samples_upper"]) + 1,
    ))
    tau, _sigma, _rs, freq = SAMPLE._partition(c)

    labelled = [[set() for _ in gaps] for _ in tau]
    for f in freq:
        target = PATH._tau_target(f, c)
        horizon = PATH._tau_sigma_horizon(f, c)
        for gi, gap in enumerate(gaps):
            for ti, cell in enumerate(tau):
                image = SAMPLE._ema_samples(cell, target, horizon, dt, gap)
                labelled[ti][gi].update(PATH._matching(tau, image))

    if len(tau) != 10 or gaps != list(range(13, 27)):
        raise RuntimeError("retained tau source projection changed")
    if any(not labelled[s][gi] for s in range(len(tau)) for gi in range(len(gaps))):
        raise RuntimeError("tau projection has a dead finite-clock transition")
    return c, clock, tau, gaps, labelled


def _lambda_upper(cell) -> float:
    lo = float(cell[0])
    if not (math.isfinite(lo) and lo > 0.0):
        raise RuntimeError("tau cell lost positive lower endpoint")
    return math.nextafter(1.0 / lo, math.inf)


@functools.lru_cache(maxsize=512)
def decay_budget(endpoint_tau_index: int, window_samples: int,
                 domain_path: Path = DEFAULT_DOMAIN) -> dict:
    """Maximum source-compatible decay exponent for an arbitrary endpoint phase."""
    domain_path = Path(domain_path).resolve()
    _c, clock, tau, gaps, labelled = _tau_projection(domain_path)
    e = int(endpoint_tau_index)
    N = int(window_samples)
    if not 0 <= e < len(tau):
        raise IndexError("endpoint tau index outside source partition")
    if N <= 0:
        raise ValueError("positive window sample count required")

    dt = float(clock["dt_binary32_s"])
    gmax = max(gaps)
    cap = N + 2 * gmax
    lam = [_lambda_upper(cell) for cell in tau]

    # dp[n][s] is the largest exponent of any complete-segment path containing
    # exactly n applied samples and ending at the stage boundary whose next
    # applied tuple is s.  A proof window may begin from any source state.
    neg = -math.inf
    dp = [[neg] * len(tau) for _ in range(cap + 1)]
    for s in range(len(tau)):
        dp[0][s] = 0.0

    for n in range(cap + 1):
        row = dp[n]
        for s, base in enumerate(row):
            if not math.isfinite(base):
                continue
            for gi, gap in enumerate(gaps):
                n2 = n + gap
                if n2 > cap:
                    continue
                weight = math.nextafter(gap * dt * lam[s], math.inf)
                value = math.nextafter(base + weight, math.inf)
                for t in labelled[s][gi]:
                    if value > dp[n2][t]:
                        dp[n2][t] = value

    # To cover an arbitrary endpoint phase in state e, extend right through one
    # complete e segment.  The left extension is already represented by allowing
    # the complete cover to be up to 2*gmax samples longer than the true window.
    finite = neg
    finite_samples = None
    finite_final_gap = None
    outgoing_gaps = [gap for gi, gap in enumerate(gaps) if labelled[e][gi]]
    for n in range(cap + 1):
        base = dp[n][e]
        if not math.isfinite(base):
            continue
        for gap in outgoing_gaps:
            total = n + gap
            if not (N <= total <= cap):
                continue
            value = math.nextafter(base + gap * dt * lam[e], math.inf)
            if value > finite:
                finite = value
                finite_samples = total
                finite_final_gap = gap

    # Infinite-time binary64 clock stagnation: no future staging, so the
    # committed tuple can hold for the entire requested physical window.
    frozen = math.nextafter(N * dt * lam[e], math.inf)
    budget = max(frozen, finite)
    if not math.isfinite(budget) or budget <= 0.0:
        raise RuntimeError("failed to produce positive tau decay budget")

    return {
        "endpoint_tau_index": e,
        "endpoint_tau_s": [float(tau[e][0]), float(tau[e][1])],
        "window_samples": N,
        "window_s_nominal": N * dt,
        "clock_dt_binary32_s": dt,
        "gap_alphabet_samples": gaps,
        "full_segment_cover_samples_upper": cap,
        "finite_clock_decay_exponent_upper": finite,
        "finite_clock_cover_samples": finite_samples,
        "finite_clock_final_gap_samples": finite_final_gap,
        "frozen_clock_decay_exponent_upper": frozen,
        "decay_exponent_upper": budget,
        "exp_minus_decay_lower_role": "evaluate with validated exponential backend in consumer",
    }


def build(domain_path: Path = DEFAULT_DOMAIN, window_samples=(50, 100, 200, 300)) -> dict:
    domain_path = Path(domain_path).resolve()
    _c, clock, tau, gaps, labelled = _tau_projection(domain_path)
    rows = {}
    for N in map(int, window_samples):
        rows[str(N)] = [decay_budget(i, N, domain_path) for i in range(len(tau))]
    failures = []
    if gaps != list(range(13, 27)):
        failures.append("gap alphabet changed")
    if len(tau) != 10:
        failures.append("tau partition is not ten cells")
    if any(not labelled[s][gi] for s in range(10) for gi in range(len(gaps))):
        failures.append("tau projection has dead transition")
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_CLOCK_PHASE_TAU_DECAY_BUDGET",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "physical_tau_cells": len(tau),
        "clock_phase_gap_alphabet_samples": gaps,
        "clock": clock,
        "arbitrary_window_phase_covered_by_full_segment_extension": True,
        "frozen_clock_hold_branch_covered": True,
        "sigma_RS_projection_only_adds_paths": True,
        "windows": rows,
        "P3_PROMOTED": False,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_only", "arbitrary_window_phase_covered_by_full_segment_extension",
        "frozen_clock_hold_branch_covered", "sigma_RS_projection_only_adds_paths",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("trajectory_replay_used", "filter_changed", "P3_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("physical_tau_cells") != 10:
        f.append("tau partition changed")
    if d.get("clock_phase_gap_alphabet_samples") != list(range(13, 27)):
        f.append("clock gap alphabet changed")
    for rows in d.get("windows", {}).values():
        if len(rows) != 10:
            f.append("window missing tau rows")
            continue
        for row in rows:
            x = row.get("decay_exponent_upper")
            if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
                f.append("invalid decay exponent")
                break
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--window-samples", type=int, nargs="*", default=[50, 100, 200, 300])
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.window_samples)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    summary = {}
    for N, rows in d["windows"].items():
        summary[N] = {
            "tau0": rows[0]["decay_exponent_upper"],
            "tau9": rows[9]["decay_exponent_upper"],
            "max": max(r["decay_exponent_upper"] for r in rows),
        }
    print(json.dumps({"validation_failures": vf, "summary": summary}, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
