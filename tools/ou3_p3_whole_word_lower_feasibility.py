#!/usr/bin/env python3
"""Feasibility diagnostic for a whole-word translation covariance *lower*.

PR #476 established that the four-max "same-history" covariance upper is not a
sufficient source-correlation quotient: a legal P2 V1 path attains the full
global adverse label in 101 samples, far inside the 635-sample covariance word.
Its recorded repair has two structural parts, and it claims only the second:

1. a whole-word covariance **lower** valid for every admissible PSD initial
   covariance and every legal source history;
2. a time-ordered, duration-aware covariance **upper**.

This module is evidence for part 1 and deliberately does not build part 2.

What it measures
----------------
Canonical ``delta`` is same-history matched: ``phase_row`` pairs the floor image
of source node ``t`` with ``endpoint_phase_upper(t, ...)``, that same node's own
upper, and then minimises over nodes.  The analogue measured here is, per source
configuration, the horizon-matched ratio

    min_i  P_word[i,i](P0 = 0)  /  P_word[i,i](P0 = infinity on v,p,S)

with both sides propagated over the *same* whole word under the *same* source
parameters, so the ratio isolates how much of the word's covariance is forced by
the word itself rather than inherited from the initial covariance.  A ratio far
above the ``1e-18`` gate means a horizon-matched lower is not self-defeating; it
does **not** certify one.

Two source models are provided:

``dwell``
    Hold one physical source node for the whole word.  For a *lower* the natural
    adversary is a word that never leaves the slowest-cadence, weakest-S cell,
    because mixing in any faster-cadence cell only adds S information.

``ordered``
    Step a supplied legal ordered source word, carrying the shipping pseudo
    timer across source commits.  ``set_pseudo_update_period_s`` applies
    ``elapsed = fmod(elapsed, new_period)``, so fixed-cell firing counts cannot
    be pasted onto a changing-source word.

Non-promoting.  These are point diagnostics at one corner of each source cell,
in double precision; they are not interval-certified, do not quantify over legal
histories, and cannot set any theorem promotion flag.  The deployed filter,
operating domain, P2 V1 source language and the canonical ``1e-18`` threshold
are untouched.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import ou3_p4_source_node_cells as NODES

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1

DT_S = 0.004999999888241291
WORD_SAMPLES = 635
PSEUDO_RATIO = 0.015 / 1.1
PSEUDO_MIN_S = DT_S
PSEUDO_MAX_S = 0.25
CHANNELS = ("v", "p", "S", "a_w")
GATE = 1.0e-18


def cadence_s(tau: float) -> float:
    """Deployed ``pseudo_update_period_for_`` including both safety clamps."""
    return min(max(PSEUDO_RATIO * float(tau), PSEUDO_MIN_S), PSEUDO_MAX_S)


def periodic_update_due(dt: float, period: float, elapsed: float):
    """Exact port of ``ou_detail::periodic_update_due``."""
    total = elapsed + dt
    tol = 16.0 * sys.float_info.epsilon * max(1.0, period)
    if total + tol < period:
        return False, total
    e = math.fmod(total, period) if total >= period else 0.0
    if not (e >= 0.0) or not math.isfinite(e) or e >= period:
        e = 0.0
    return True, e


def commit_period(elapsed: float, period: float) -> float:
    """Exact port of ``set_pseudo_update_period_s``'s timer rebase."""
    e = math.fmod(float(elapsed), float(period))
    if not (math.isfinite(e) and 0.0 <= e < period):
        raise RuntimeError("invalid fmod pseudo-timer state")
    return e


def _mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def _sym(A):
    B = [[float(x) for x in row] for row in A]
    for i in range(4):
        for j in range(i + 1, 4):
            v = 0.5 * (B[i][j] + B[j][i])
            B[i][j] = B[j][i] = v
    return B


def _transition_and_noise(tau: float, sigma: float, h: float):
    """Exact discrete map of dv=a, dp=v, dS=p, da=-a/tau dt + sqrt(2 sig^2/tau) dW.

    Series-summed rather than matrix-exponentiated so the module has no
    numerical dependency beyond the standard library.
    """
    a = h / float(tau)
    e = math.exp(-a)
    # Integrals of the OU kernel against 1, t, t^2/2 over one sample.
    t = float(tau)
    i0 = t * (1.0 - e)
    i1 = t * (h - i0)
    i2 = t * (0.5 * h * h - i1)
    Phi = [[1.0, 0.0, 0.0, i0],
           [h, 1.0, 0.0, i1],
           [0.5 * h * h, h, 1.0, i2],
           [0.0, 0.0, 0.0, e]]
    q = 2.0 * float(sigma) * float(sigma) / t
    # Van Loan by direct quadrature: Q = int_0^h Phi(s) G q G' Phi(s)' ds with
    # G = e_4.  Simpson on a fixed fine grid is exact enough for a diagnostic.
    n = 64
    Q = [[0.0] * 4 for _ in range(4)]
    for k in range(n + 1):
        s = h * k / n
        w = (1.0 if k in (0, n) else (4.0 if k % 2 else 2.0)) * (h / n) / 3.0
        es = math.exp(-s / t)
        c0 = t * (1.0 - es)
        c1 = t * (s - c0)
        c2 = t * (0.5 * s * s - c1)
        col = (c0, c1, c2, es)
        for i in range(4):
            for j in range(4):
                Q[i][j] += w * q * col[i] * col[j]
    return Phi, _sym(Q)


def _measure_S(P, R: float):
    den = P[2][2] + float(R)
    if not (math.isfinite(den) and den > 0.0):
        raise RuntimeError("S measurement denominator lost positivity")
    c = [P[i][2] for i in range(4)]
    return _sym([[P[i][j] - c[i] * c[j] / den for j in range(4)] for i in range(4)])


def run_word(segments, P0_diag, *, word_samples: int = WORD_SAMPLES):
    """Propagate ``P0_diag`` over ``segments`` = [(tau, sigma, R_S_std, gap)]."""
    P = [[0.0] * 4 for _ in range(4)]
    for i, x in enumerate(P0_diag):
        P[i][i] = float(x)
    elapsed = 0.0
    fires = 0
    used = 0
    kernels: dict = {}
    for (tau, sigma, rs_std, gap) in segments:
        period = cadence_s(tau)
        elapsed = commit_period(elapsed, period)
        key = (round(float(tau), 12), round(float(sigma), 12))
        if key not in kernels:
            kernels[key] = _transition_and_noise(tau, sigma, DT_S)
        Phi, Q = kernels[key]
        R = float(rs_std) ** 2
        PhiT = [list(r) for r in zip(*Phi)]
        for _ in range(int(gap)):
            if used >= word_samples:
                break
            M = _mm(_mm(Phi, P), PhiT)
            P = _sym([[M[i][j] + Q[i][j] for j in range(4)] for i in range(4)])
            used += 1
            due, elapsed = periodic_update_due(DT_S, period, elapsed)
            if due:
                fires += 1
                P = _measure_S(P, R)
        if used >= word_samples:
            break
    return P, fires, used


P0_SCALE_START = 1.0e3
P0_SCALE_STEP = 1.0e3
P0_SCALE_MAX = 1.0e15
P0_TOLERANCE = 1.0e-4


def _ratio(segments, sigma_for_aw: float):
    """Horizon-matched ratio with an adaptive initial-covariance probe.

    The upper side wants ``P0 = infinity`` on v, p and S, which double precision
    cannot represent.  Two walls bracket the usable range and they move from cell
    to cell:

    * below the lower wall the word has not yet forgotten ``P0``, so the endpoint
      still grows with it -- the slowest-cadence cells need more than 1e9;
    * above the upper wall the covariance-form update loses the endpoint to
      cancellation, and the tell is that the endpoint starts to *decrease* as
      ``P0`` grows, which the monotone Riccati map forbids.

    The gap between the walls spans about eight orders across the 800 source
    cells, so no fixed probe pair serves all of them.  This escalates ``P0``
    relative to the word's own zero-start covariance until two consecutive probes
    agree, and stops early if monotonicity breaks.  The accepted scale and the
    achieved spread are reported per row so a caller can see which wall a
    non-converged row hit.
    """
    lo, fires, used = run_word(segments, [0.0, 0.0, 0.0, 0.0])
    aw = float(sigma_for_aw) ** 2
    base = [abs(lo[i][i]) for i in range(4)]
    if any(not (math.isfinite(b) and b > 0.0) for b in base):
        raise RuntimeError("zero-start whole-word lower lost positivity")

    def probe(scale):
        seed = [scale * base[0], scale * base[1], scale * base[2], aw]
        return run_word(segments, seed)[0]

    scale = P0_SCALE_START
    prev = probe(scale)
    best, spread, accepted, converged = prev, float("inf"), scale, False
    while scale < P0_SCALE_MAX:
        scale *= P0_SCALE_STEP
        cur = probe(scale)
        sp = 0.0
        decreased = False
        for i in range(4):
            a, b = abs(prev[i][i]), abs(cur[i][i])
            if not (math.isfinite(b) and b > 0.0):
                decreased = True
                break
            sp = max(sp, abs(a - b) / max(a, b))
            if b < a * (1.0 - P0_TOLERANCE):
                decreased = True
        if decreased:
            # Upper wall: keep the last value produced below it.
            break
        best, spread, accepted = cur, sp, scale
        if sp <= P0_TOLERANCE:
            converged = True
            break
        prev = cur

    r = [base[i] / abs(best[i][i]) for i in range(4)]
    k = r.index(min(r))
    return {"channel_ratios": r, "ratio": r[k], "binding_channel": CHANNELS[k],
            "S_firings": fires, "samples": used,
            "P0_probe_relative_spread": spread,
            "P0_accepted_scale": accepted,
            "P0_independent": bool(converged)}


def _adverse_corner(node) -> tuple:
    """Slowest cadence, largest sigma, weakest S measurement in the cell."""
    return (float(node["tau_s"][1]),
            float(node["sigma_filter_committed_mps2"][1]),
            float(node["R_S_filter_std"][1]))


def build(domain_path: Path = DEFAULT_DOMAIN, *, stride: int = 1) -> dict:
    nodes = NODES.build()["nodes"]
    rows = []
    for node in nodes[::max(1, int(stride))]:
        tau, sigma, rs = _adverse_corner(node)
        res = _ratio([(tau, sigma, rs, WORD_SAMPLES)], sigma)
        rows.append({"source_node": int(node["index"]),
                     "tau_index": int(node["tau_index"]),
                     "sigma_raw_index": int(node["sigma_raw_index"]),
                     "R_S_index": int(node["R_S_index"]),
                     "tau_s": tau, "sigma_aw": sigma, "R_S_filter_std": rs,
                     **res})
    rows.sort(key=lambda x: x["ratio"])
    worst = rows[0] if rows else None
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_WHOLE_WORD_LOWER_FEASIBILITY_POINT_DIAGNOSTIC",
        "non_promoting": True,
        "certifies_theorem_stage": False,
        "interval_certified": False,
        "quantifies_over_legal_histories": False,
        "source_model": "single_cell_dwell_adverse_corner",
        "word_samples": WORD_SAMPLES,
        "dt_s": DT_S,
        "canonical_gate": GATE,
        "nodes_evaluated": len(rows),
        "stride": int(stride),
        "worst": worst,
        "worst_ratio": None if worst is None else worst["ratio"],
        "worst_clears_canonical_gate": bool(worst is not None and worst["ratio"] > GATE),
        "rows": rows,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA:
        f.append("schema changed")
    for key in ("non_promoting", "interval_certified", "quantifies_over_legal_histories",
                "certifies_theorem_stage"):
        if key not in d:
            f.append(f"missing honesty flag {key}")
    if d.get("non_promoting") is not True:
        f.append("diagnostic must remain non-promoting")
    if d.get("certifies_theorem_stage") is not False:
        f.append("diagnostic must not claim a theorem stage")
    if d.get("interval_certified") is not False:
        f.append("point diagnostic must not claim interval certification")
    if float(d.get("canonical_gate", 0.0)) != GATE:
        f.append("canonical usefulness gate must remain exactly 1e-18")
    if int(d.get("word_samples", 0)) != WORD_SAMPLES:
        f.append("covariance word sample count changed")
    rows = d.get("rows") or []
    if not rows:
        f.append("no source rows evaluated")
    for r in rows:
        if not (math.isfinite(r["ratio"]) and r["ratio"] > 0.0):
            f.append(f"node {r['source_node']} produced a non-positive ratio")
            break
        if r["binding_channel"] not in CHANNELS:
            f.append(f"node {r['source_node']} reported an unknown binding channel")
            break
        if int(r["samples"]) != WORD_SAMPLES:
            f.append(f"node {r['source_node']} did not cover the whole word")
            break
        if not r.get("P0_independent"):
            f.append(f"node {r['source_node']} did not forget its initial covariance "
                     f"(probe spread {r.get('P0_probe_relative_spread'):.3g})")
            break
    return list(dict.fromkeys(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--stride", type=int, default=1,
                    help="evaluate every Nth source node (CI budget control)")
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, stride=a.stride)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({k: v for k, v in d.items() if k != "rows"}, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
