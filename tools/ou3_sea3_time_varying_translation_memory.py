#!/usr/bin/env python3
"""Time-varying complete-SEA3 translation process-memory lower bound.

This is a process/noise lemma inside the canonical complete SEA3 Normal-Live
word.  It does not freeze the tuner over a multi-second word and it does not
construct a predecessor/source-history graph.

Shipping commits the active OU schedule only at the staged adaptation clock.
Between two active commits, tau/sigma/R_S are constant.  The dynamic-source
certificate gives an active-commit gap of at most 22 IMU samples.  The source
condition ``time-last > ADAPT_EVERY_SECS`` with ADAPT_EVERY_SECS=0.1 s and the
5 ms IMU clock gives a conservative lower of 20 samples even if binary32 time
rounding makes the strict boundary fire one sample early.  In a 600-sample
(3 s) word, discarding at most one <=22-sample boundary interval at each end
leaves at least ceil((600-44)/22)=26 complete constant-tune intervals.

For one translation axis x=[v,p,S,a_w], keep a deterministic covariance lower
L_k after k such intervals.  On each possible current tau cell and for every
allowed interval length n in {20,21,22}, execute the source-uniform lower
Riccati recursion sample by sample:

    L- = F(tau) L F(tau)' + Q(tau,sigma_floor)
    L+ = (L-^-1 + D_upper)^-1.

The optimal-posterior map is a lower bound for the implemented Joseph update
for any gain.  D_upper deliberately assimilates an S=0 measurement on *every*
IMU sample and the factor-three accelerometer a_w information on every sample;
therefore it is at least as informative as the actual translation part of the
complete A21 word.  Scalar diagonal updates are used instead of a matrix
inverse, preserving interval conditioning.

A fixed candidate sequence is generated from the worst-reference tau=12 s,
20-sample macrostep and then weakened by alpha_k=0.8*0.98^(k-1).  Candidates
are not evidence: every induction step is recertified over the full tau cover
with outward-rounded interval LDLT after the same rational RL^-1 congruence used
by the shipping integrated-OU process certificate.

The final <=22-sample suffix after the last selected complete interval is also
certified not to reduce L_26.  Thus the emitted 4x4 matrix is a word-endpoint
lower bound for the actual process/measurement-noise covariance for arbitrary
legal time-varying committed tau/sigma/R_S inside the same complete SEA3 word.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import functools
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_add,
    matrix_mul,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_four_s_translation_information_tight as FOUR
import ou3_sea3_riccati_tube as TUBE
import ou3_sea3_riccati_tube_factored as FACTORED
import ou3_sea3_windowed_vector_pe as PE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = COMPLETE.DEFAULT_DOMAIN
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_TIME_VARYING_TRANSLATION_PROCESS_MEMORY"
WORD_HORIZON_S = 3.0
MACRO_INTERVALS = 26
ALPHA_FIRST = 0.8
ALPHA_DECAY = 0.98
BASE_X_CELLS = 48
MAX_X_SPLIT_DEPTH = 6
SERIES_ORDER = 12


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _zero(n: int) -> IntervalMatrix:
    z = I(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _diag(values: Sequence[float]) -> IntervalMatrix:
    M = _zero(len(values))
    for i, v in enumerate(values):
        M[i][i] = I(float(v))
    return M


def _frac_I(q: Fraction) -> Interval:
    return I(float(q))


# Same exact rational conditioning used by the established process primitive.
F = Fraction
_L_INV = (
    (F(1), F(0), F(0), F(0)),
    (-F(3, 8), F(1), F(0), F(0)),
    (F(1, 15), -F(4, 9), F(1), F(0)),
    (-F(15, 2), F(30), -F(105, 2), F(1)),
)
_R_DIAG = (F(1), F(10), F(100), F(2))
_C_RATIONAL = tuple(
    tuple(_R_DIAG[i] * _L_INV[i][j] for j in range(4))
    for i in range(4)
)
_C = [[_frac_I(v) for v in row] for row in _C_RATIONAL]
_CT = matrix_transpose(_C)


def _normalized_integral_series(x: Interval, offset: int) -> Interval:
    """Validated f_offset(x)=sum_n (-1)^n x^n/(n+offset)!.

    offset=1 gives (1-exp(-x))/x, offset=2 gives
    (x-1+exp(-x))/x^2, and offset=3 gives
    (x^2/2-x+1-exp(-x))/x^3.
    """
    if offset not in (1, 2, 3) or x.lo <= 0.0:
        raise ValueError("positive x and offset 1..3 required")
    y = I(((-1.0) ** SERIES_ORDER) / math.factorial(SERIES_ORDER + offset))
    for n in range(SERIES_ORDER - 1, -1, -1):
        c = ((-1.0) ** n) / math.factorial(n + offset)
        y = I(c) + x * y
    # Alternating terms decrease throughout the deployed x<0.01 sample domain.
    rem = TUBE.up(
        (x.hi ** (SERIES_ORDER + 1))
        / math.factorial(SERIES_ORDER + 1 + offset)
    )
    return y + Interval.outward_bounds(-rem, rem)


def _sample_transition(x: Interval, h: float) -> IntervalMatrix:
    if not (h > 0.0 and x.lo > 0.0):
        raise ValueError("positive sample time/x required")
    e = VT.exp_interval(-x)
    f1 = _normalized_integral_series(x, 1)
    f2 = _normalized_integral_series(x, 2)
    f3 = _normalized_integral_series(x, 3)
    hI = I(h)
    h2 = hI * hI
    h3 = h2 * hI
    z, o = I(0.0), I(1.0)
    return [
        [o, z, z, hI * f1],
        [hI, o, z, h2 * f2],
        [I(0.5) * h2, hI, o, h3 * f3],
        [z, z, z, e],
    ]


def _sample_process(x: Interval, h: float, sigma_floor: float) -> IntervalMatrix:
    """Exact shipping integrated-OU Q at sigma_floor on one tau cell."""
    qscaled = FACTORED.step_scaled_q(x)
    s = float(sigma_floor)
    scales = [s * h, s * h * h, s * h * h * h, s]
    return matrix_symmetric_hull([
        [I(scales[i]) * qscaled[i][j] * I(scales[j]) for j in range(4)]
        for i in range(4)
    ])


def _scalar_information_update(P: IntervalMatrix, index: int, information: float) -> IntervalMatrix:
    """Exact covariance update P-Pe e'P/(R+e'Pe), R=1/information."""
    d = float(information)
    if not d > 0.0:
        return P
    R = I(TUBE.up(1.0 / TUBE.down(d)))
    den = R + P[index][index]
    if not den.lo > 0.0:
        raise RuntimeError("scalar information update denominator lost positivity")
    n = len(P)
    out = _zero(n)
    for i in range(n):
        for j in range(n):
            out[i][j] = P[i][j] - (P[i][index] * P[index][j]) / den
    return matrix_symmetric_hull(out)


def _sample_lower_map(
    P: IntervalMatrix,
    x: Interval,
    *,
    h: float,
    sigma_floor: float,
    info_S: float,
    info_aw: float,
) -> IntervalMatrix:
    Fm = _sample_transition(x, h)
    Qm = _sample_process(x, h, sigma_floor)
    pred = matrix_symmetric_hull(
        matrix_add(matrix_mul(matrix_mul(Fm, P), matrix_transpose(Fm)), Qm)
    )
    # Deliberately assume S is due on every sample, then accelerometer a_w info.
    post = _scalar_information_update(pred, 2, info_S)
    post = _scalar_information_update(post, 3, info_aw)
    return post


def _macro_lower_map(
    P0: IntervalMatrix,
    x: Interval,
    samples: int,
    *,
    h: float,
    sigma_floor: float,
    info_S: float,
    info_aw: float,
) -> IntervalMatrix:
    P = [[v for v in row] for row in P0]
    for _ in range(samples):
        P = _sample_lower_map(
            P, x, h=h, sigma_floor=sigma_floor, info_S=info_S, info_aw=info_aw
        )
    return P


def _float_mul(A, B):
    return [
        [sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
        for i in range(len(A))
    ]


def _float_transpose(A):
    return [list(x) for x in zip(*A)]


def _float_add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A))] for i in range(len(A))]


def _float_sym(A):
    return [[0.5 * (A[i][j] + A[j][i]) for j in range(len(A))] for i in range(len(A))]


def _float_scalar_update(P, index: int, information: float):
    if not information > 0.0:
        return P
    R = 1.0 / information
    den = R + P[index][index]
    out = [[0.0] * len(P) for _ in range(len(P))]
    for i in range(len(P)):
        for j in range(len(P)):
            out[i][j] = P[i][j] - P[i][index] * P[index][j] / den
    return _float_sym(out)


def _mid(A: IntervalMatrix):
    return [[0.5 * (a.lo + a.hi) for a in row] for row in A]


def _reference_sequence(
    *,
    h: float,
    tau_reference: float,
    sigma_floor: float,
    macro_samples: int,
    info_S: float,
    info_aw: float,
):
    x = Interval.outward_bounds(h / tau_reference, h / tau_reference)
    Fm = _mid(_sample_transition(x, h))
    Qm = _mid(_sample_process(x, h, sigma_floor))
    P = [[0.0] * 4 for _ in range(4)]
    refs = [P]
    for _macro in range(MACRO_INTERVALS):
        for _ in range(macro_samples):
            P = _float_add(_float_mul(_float_mul(Fm, P), _float_transpose(Fm)), Qm)
            P = _float_scalar_update(P, 2, info_S)
            P = _float_scalar_update(P, 3, info_aw)
        P = _float_sym(P)
        refs.append(P)
    return refs


def _alpha(k: int) -> float:
    if k <= 0:
        return 0.0
    return math.nextafter(ALPHA_FIRST * (ALPHA_DECAY ** (k - 1)), -math.inf)


def _candidate(refs, k: int) -> IntervalMatrix:
    a = _alpha(k)
    return [[I(a * refs[k][i][j]) for j in range(4)] for i in range(4)]


def _conditioned_difference(A: IntervalMatrix, *, sigma_floor: float, horizon: float) -> IntervalMatrix:
    s = float(sigma_floor)
    T = float(horizon)
    if not (s > 0.0 and T > 0.0):
        raise ValueError("positive conditioning scale required")
    inv = [1.0 / (s * T), 1.0 / (s * T * T), 1.0 / (s * T * T * T), 1.0 / s]
    G = matrix_mul(_C, _diag(inv))
    return matrix_symmetric_hull(matrix_mul(matrix_mul(G, A), matrix_transpose(G)))


def _subtract(A: IntervalMatrix, B: IntervalMatrix) -> IntervalMatrix:
    return matrix_symmetric_hull([
        [A[i][j] - B[i][j] for j in range(len(A))]
        for i in range(len(A))
    ])


def _split_x(x: Interval) -> tuple[Interval, Interval]:
    mid = math.sqrt(x.lo * x.hi)
    return Interval.outward_bounds(x.lo, mid), Interval.outward_bounds(mid, x.hi)


def _certify_macro_cell(
    x: Interval,
    *,
    L0: IntervalMatrix,
    L1: IntervalMatrix,
    sample_counts: range,
    h: float,
    sigma_floor: float,
    info_S: float,
    info_aw: float,
    horizon_ref: float,
) -> tuple[bool, float, str]:
    worst = math.inf
    try:
        for n in sample_counts:
            P = _macro_lower_map(
                L0, x, n, h=h, sigma_floor=sigma_floor, info_S=info_S, info_aw=info_aw
            )
            D = _conditioned_difference(
                _subtract(P, L1), sigma_floor=sigma_floor, horizon=horizon_ref
            )
            ok, pivots = symmetric_positive_definite_ldlt(D)
            if not ok:
                return False, -math.inf, f"LDLT failed at {n} samples"
            worst = min(worst, min(p.lo for p in pivots))
    except Exception as exc:
        return False, -math.inf, f"{type(exc).__name__}: {exc}"
    return True, worst, ""


def _certify_with_splitting(
    x: Interval,
    *,
    depth: int,
    kwargs: dict,
    stats: dict,
) -> bool:
    ok, pivot, err = _certify_macro_cell(x, **kwargs)
    if ok:
        stats["leaves"] += 1
        stats["worst_pivot"] = min(stats["worst_pivot"], pivot)
        stats["max_depth"] = max(stats["max_depth"], depth)
        return True
    if depth >= MAX_X_SPLIT_DEPTH:
        stats["failures"].append({"x": x.as_list(), "depth": depth, "error": err})
        return False
    stats["splits"] += 1
    left, right = _split_x(x)
    a = _certify_with_splitting(left, depth=depth + 1, kwargs=kwargs, stats=stats)
    b = _certify_with_splitting(right, depth=depth + 1, kwargs=kwargs, stats=stats)
    return a and b


@functools.lru_cache(maxsize=4)
def _build_cached(path_text: str) -> dict:
    path = Path(path_text)
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    four = FOUR.build(path)
    pe = PE.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "four_S": FOUR.validate(four),
        "PE": PE.validate(pe),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"time-varying translation prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("translation memory detached from complete SEA3")
    if float(complete["word_horizon_s"]) != WORD_HORIZON_S:
        raise RuntimeError("complete SEA3 word horizon changed")

    rates = dynamic["validated_rate_and_jump_bounds"]
    inv = dynamic["dynamic_invariant"]
    h = float(rates["dt_s"])
    max_commit_samples = int(rates["active_commit_gap_samples_upper"])
    word_samples = int(round(WORD_HORIZON_S / h))
    # 20 samples is conservative for the strict >0.1 s source clock even if
    # binary32 accumulated time moves the equality boundary by one ulp.
    min_commit_samples = int(math.floor(0.1 / h + 1.0e-9))
    interior_samples = word_samples - 2 * max_commit_samples
    full_intervals_lower = int(math.ceil(interior_samples / max_commit_samples))
    if full_intervals_lower < MACRO_INTERVALS:
        raise RuntimeError("3 s word no longer guarantees 26 complete commit intervals")

    wrapper = WRAPPER.read_text(encoding="utf-8")
    parity = {
        "staging_condition_strict_0p1s_clock": (
            "time_ - last_adapt_time_sec_ > adapt_every_secs_" in wrapper
        ),
        "pending_commit_applied_next_sample": (
            "online_tune_apply_pending_ = true" in wrapper
            and "void apply_pending_online_tune_()" in wrapper
        ),
        "active_tau_and_sigma_commit_together": (
            "mekf_->set_aw_time_constant(tune_.tau_applied);" in wrapper
            and "const float sZ = std::max(sigma_floor, tune_.sigma_applied);" in wrapper
        ),
    }
    if not all(parity.values()):
        raise RuntimeError(f"shipping commit parity failed: {parity}")

    sigma_floor = float(inv["sigma_aw_filter_mps2"][0])
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    if not (sigma_floor > 0.0 and tau_lo > 0.0 and tau_hi >= tau_lo):
        raise RuntimeError("invalid dynamic process invariant")

    meas = pe["measurement_runtime"]
    ra = float(meas["accelerometer_variance_upper"])
    rs_lo = float(inv["R_S_applied"][0])
    axis_factor_min = min(map(float, four["R_S_axis_std_factors"]))
    rs_std_min = TUBE.down(rs_lo * axis_factor_min)
    if not all(x > 0.0 for x in (ra, rs_std_min)):
        raise RuntimeError("measurement-information upper lost positivity")
    info_S = TUBE.up(1.0 / TUBE.down(rs_std_min * rs_std_min))
    # Factor three is the full-A accelerometer theta/a_w/b_a cross-term payment.
    info_aw = TUBE.up(3.0 / ra)

    refs = _reference_sequence(
        h=h,
        tau_reference=tau_hi,
        sigma_floor=sigma_floor,
        macro_samples=min_commit_samples,
        info_S=info_S,
        info_aw=info_aw,
    )

    xlo, xhi = h / tau_hi, h / tau_lo
    edges = TUBE.geom_edges(xlo, xhi, BASE_X_CELLS)
    base_cells = TUBE.interval_cells(edges)
    stats = {
        "leaves": 0,
        "splits": 0,
        "max_depth": 0,
        "worst_pivot": math.inf,
        "failures": [],
    }
    all_ok = True
    for k in range(MACRO_INTERVALS):
        L0 = _candidate(refs, k)
        L1 = _candidate(refs, k + 1)
        kwargs = {
            "L0": L0,
            "L1": L1,
            "sample_counts": range(min_commit_samples, max_commit_samples + 1),
            "h": h,
            "sigma_floor": sigma_floor,
            "info_S": info_S,
            "info_aw": info_aw,
            "horizon_ref": (k + 1) * min_commit_samples * h,
        }
        for cell in base_cells:
            if not _certify_with_splitting(cell, depth=0, kwargs=kwargs, stats=stats):
                all_ok = False

    final_lower = _candidate(refs, MACRO_INTERVALS)

    # The selected 26 intervals may end at the final complete commit before the
    # word endpoint.  The remaining suffix is at most max_commit_samples.  Show
    # that any 1..max gap of worst-case measured samples cannot reduce L_26;
    # zero suffix is equality and needs no strict LDLT.
    tail_stats = {
        "leaves": 0,
        "splits": 0,
        "max_depth": 0,
        "worst_pivot": math.inf,
        "failures": [],
    }
    tail_ok = True
    tail_kwargs = {
        "L0": final_lower,
        "L1": final_lower,
        "sample_counts": range(1, max_commit_samples + 1),
        "h": h,
        "sigma_floor": sigma_floor,
        "info_S": info_S,
        "info_aw": info_aw,
        "horizon_ref": MACRO_INTERVALS * min_commit_samples * h,
    }
    for cell in base_cells:
        if not _certify_with_splitting(cell, depth=0, kwargs=tail_kwargs, stats=tail_stats):
            tail_ok = False

    final_float = [[v.lo for v in row] for row in final_lower]
    closed = all_ok and tail_ok and not stats["failures"] and not tail_stats["failures"]
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "complete_word_horizon_s": WORD_HORIZON_S,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "constant_tau_over_memory_assumed": False,
        "time_varying_committed_tau_sigma_RS_allowed": True,
        "shipping_commit_parity": parity,
        "commit_geometry": {
            "word_samples": word_samples,
            "min_constant_commit_interval_samples_conservative": min_commit_samples,
            "max_commit_interval_samples_certified": max_commit_samples,
            "boundary_intervals_discarded": 2,
            "interior_samples_lower": interior_samples,
            "complete_constant_tune_intervals_lower": full_intervals_lower,
            "intervals_retained": MACRO_INTERVALS,
        },
        "measurement_lower_recursion": {
            "optimal_posterior_lower_for_any_Joseph_gain": True,
            "S_zero_assumed_due_every_sample": True,
            "S_information_per_sample_upper": info_S,
            "accelerometer_aw_information_per_sample_upper": info_aw,
            "accelerometer_cross_block_factor_three_paid": True,
            "measurements_moved_to_word_endpoint": False,
            "measurements_executed_sample_by_sample": True,
        },
        "process": {
            "sigma_floor_mps2": sigma_floor,
            "tau_applied_s": [tau_lo, tau_hi],
            "sample_x_h_over_tau": [xlo, xhi],
            "exact_shipping_integrated_OU_Q_consumed": True,
        },
        "candidate_induction": {
            "alpha_first": ALPHA_FIRST,
            "alpha_decay_per_commit_interval": ALPHA_DECAY,
            "final_alpha": _alpha(MACRO_INTERVALS),
            "candidate_is_only_verified_target_not_source_assumption": True,
            "base_x_cells": BASE_X_CELLS,
            "certified_leaves": stats["leaves"],
            "adaptive_splits": stats["splits"],
            "max_split_depth": stats["max_depth"],
            "worst_conditioned_LDLT_pivot_lower": (
                stats["worst_pivot"] if math.isfinite(stats["worst_pivot"]) else None
            ),
            "failures": stats["failures"],
        },
        "terminal_suffix": {
            "max_samples": max_commit_samples,
            "zero_sample_suffix_is_identity": True,
            "positive_suffix_preserves_final_lower": tail_ok,
            "certified_leaves": tail_stats["leaves"],
            "adaptive_splits": tail_stats["splits"],
            "max_split_depth": tail_stats["max_depth"],
            "worst_conditioned_LDLT_pivot_lower": (
                tail_stats["worst_pivot"] if math.isfinite(tail_stats["worst_pivot"]) else None
            ),
            "failures": tail_stats["failures"],
        },
        "word_endpoint_translation_process_measurement_noise_covariance_lower": final_float,
        "full_4x4_time_varying_translation_memory_closed": closed,
        "P3_promoted": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    return _build_cached(str(Path(domain_path).resolve()))


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "time_varying_committed_tau_sigma_RS_allowed",
        "full_4x4_time_varying_translation_memory_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_family_replaced", "trajectory_replay_used", "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed", "constant_tau_over_memory_assumed", "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"forbidden/fail-closed flag {key} changed")
    geom = d.get("commit_geometry", {})
    if int(geom.get("complete_constant_tune_intervals_lower", 0)) < MACRO_INTERVALS:
        f.append("fewer than 26 complete constant-tune intervals certified")
    meas = d.get("measurement_lower_recursion", {})
    for key in (
        "optimal_posterior_lower_for_any_Joseph_gain",
        "S_zero_assumed_due_every_sample",
        "accelerometer_cross_block_factor_three_paid",
        "measurements_executed_sample_by_sample",
    ):
        if meas.get(key) is not True:
            f.append(f"measurement lower recursion lost {key}")
    if meas.get("measurements_moved_to_word_endpoint") is not False:
        f.append("measurements were illegally moved to the word endpoint")
    M = d.get("word_endpoint_translation_process_measurement_noise_covariance_lower")
    if not isinstance(M, list) or len(M) != 4 or any(not isinstance(r, list) or len(r) != 4 for r in (M or [])):
        f.append("final translation lower is not 4x4")
    if d.get("candidate_induction", {}).get("failures"):
        f.append("time-varying candidate induction has failed cells")
    if d.get("terminal_suffix", {}).get("failures"):
        f.append("terminal suffix has failed cells")
    if d.get("terminal_suffix", {}).get("positive_suffix_preserves_final_lower") is not True:
        f.append("terminal suffix does not preserve translation lower")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "closed": d["full_4x4_time_varying_translation_memory_closed"],
        "commit_geometry": d["commit_geometry"],
        "induction": d["candidate_induction"],
        "tail": d["terminal_suffix"],
        "lower": d["word_endpoint_translation_process_measurement_noise_covariance_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
