#!/usr/bin/env python3
"""Optimize the paper's spread-selected four-S translation UCO bound.

For a declared source-word horizon T_W, every integer q satisfying
3 q Delta_max <= T_W gives four guaranteed usable S updates separated by at
least q Delta_min.  The determinant grows as q^6, while OU decay and the
Frobenius bound worsen with the selected subwindow.  This producer evaluates
all admissible integer q with outward-rounded/validated primitives and chooses
the largest certified information lower bound.  The choice is a proof-design
choice over deterministic source bounds, not trajectory fitting.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_validated_transcendentals as VT

DEFAULT_HEADER = SOURCE.DEFAULT_HEADER
SCHEMA = 1


def _point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _pow_nonnegative(x: Interval, n: int) -> Interval:
    if x.lo < 0.0 or n < 0:
        raise ValueError("nonnegative integer power required")
    y = Interval.point(1.0)
    for _ in range(n):
        y = y * x
    return y


def _exp_negative_wide(x: float) -> float:
    if not math.isfinite(x) or x < 0.0:
        raise ValueError("finite nonnegative exponential argument required")
    scale = 1
    while x / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    z = Interval.outward_bounds(x / scale, x / scale)
    y = VT.exp_interval(-z)
    s = scale
    while s > 1:
        y = y.square()
        s //= 2
    return y.lo


def _candidate(q: int, delta_min: float, delta_max: float, tau_min: float,
               rs_upper: float) -> dict:
    selected_window = math.nextafter(3.0 * q * delta_max, math.inf)
    spacing = math.nextafter(q * delta_min, -math.inf)
    decay = _exp_negative_wide(selected_window / tau_min)
    det = math.nextafter((spacing ** 6 / 12.0) * decay, -math.inf)

    T = Interval(0.0, math.nextafter(selected_window, math.inf))
    row_norm2 = (
        _point(1.0) + T.square()
        + _pow_nonnegative(T, 4) / _point(4.0)
        + _pow_nonnegative(T, 6) / _point(36.0)
    )
    frob = math.nextafter(2.0 * math.sqrt(row_norm2.hi), math.inf)
    sigma_min = math.nextafter(det / (frob ** 3), -math.inf)
    info = math.nextafter((sigma_min ** 2) / (rs_upper ** 2), -math.inf)
    return {
        "q": q,
        "selected_window_s_upper": selected_window,
        "selected_spacing_s_lower": spacing,
        "exp_minus_window_over_tau_min_lower": decay,
        "observation_det_lower": det,
        "observation_frobenius_upper": frob,
        "observation_sigma_min_lower": sigma_min,
        "information_gramian_lambda_min_lower": info,
    }


def build(word_horizon_s: float, header: Path = DEFAULT_HEADER.resolve()) -> dict:
    Tword = float(word_horizon_s)
    if not math.isfinite(Tword) or Tword <= 0.0:
        raise ValueError("word horizon must be finite positive")
    base = TRANS.build(Path(header).resolve())
    failures = TRANS.validate(base)
    if failures:
        raise RuntimeError(f"translation prerequisite failed: {failures}")
    s = base["S_observation_uco"]
    delta_min = float(s["pseudo_gap_min_s"])
    delta_max = float(s["pseudo_gap_max_s"])
    tau_min = float(base["process_ucc"]["tau_s"][0])
    rs_upper = float(s["R_S_filter_std_upper"])
    qmax = int(math.floor(Tword / (3.0 * delta_max)))
    if qmax < 1:
        raise RuntimeError("word horizon does not contain four source-guaranteed S firings")

    candidates = [_candidate(q, delta_min, delta_max, tau_min, rs_upper)
                  for q in range(1, qmax + 1)]
    usable = [c for c in candidates if c["information_gramian_lambda_min_lower"] > 0.0]
    if not usable:
        raise RuntimeError("no spread candidate retained a positive validated information bound")
    best = max(usable, key=lambda c: c["information_gramian_lambda_min_lower"])
    adjacent = candidates[0]
    ratio = math.nextafter(
        best["information_gramian_lambda_min_lower"] /
        adjacent["information_gramian_lambda_min_lower"], -math.inf
    )
    return {
        "schema": SCHEMA,
        "qualification": "VALIDATED_OPTIMAL_SPREAD_FOUR_S_COMPLETE_TRANSLATION_UCO",
        "source_generated_not_trajectory_fit": True,
        "outward_rounded": True,
        "word_horizon_s": Tword,
        "admissible_q_max": qmax,
        "candidate_count": len(candidates),
        "selection_rule": "maximize certified information_gramian_lambda_min_lower over every integer q with 3 q Delta_max <= T_word",
        "best": best,
        "adjacent_q1": adjacent,
        "information_widening_factor_vs_adjacent_lower": ratio,
        "three_S_detectability_used_for_this_UCO": False,
        "pass": best["information_gramian_lambda_min_lower"] > 0.0,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True or d.get("outward_rounded") is not True:
        failures.append("spread UCO is not source-generated outward-rounded")
    if d.get("three_S_detectability_used_for_this_UCO") is not False:
        failures.append("spread four-S UCO incorrectly uses three-S detectability")
    best = d.get("best", {})
    if not isinstance(best.get("q"), int) or best.get("q", 0) < 1:
        failures.append("best spread index invalid")
    for key in ("selected_spacing_s_lower", "observation_det_lower",
                "observation_sigma_min_lower", "information_gramian_lambda_min_lower"):
        x = best.get(key)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
            failures.append(f"best.{key} is not finite positive")
    ratio = d.get("information_widening_factor_vs_adjacent_lower")
    if not isinstance(ratio, (int, float)) or not math.isfinite(float(ratio)) or float(ratio) < 1.0:
        failures.append("optimized spread is worse than adjacent q=1")
    if d.get("pass") is not True:
        failures.append("spread UCO did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--word-horizon-s", type=float, required=True)
    ap.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.word_horizon_s, args.header.resolve())
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"best": d["best"], "widening": d["information_widening_factor_vs_adjacent_lower"], "failures": failures}, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
