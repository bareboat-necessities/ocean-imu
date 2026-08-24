#!/usr/bin/env python3
"""Validated one-step scalar/axis OU transition enclosure for OU-III.

This is the next proof layer after ``ou3_source_interval_box.py``. It consumes
the source-derived outward-rounded tau box and encloses the scalar OU quantities
and the exact per-axis ``(v,p,S,a_w)`` transition used by
``KalmanOUCoreMath.h``. The shipping small-x polynomial branch and the expm1
branch are covered separately, including the threshold.

No replay, fitted extrema, Monte Carlo result, ordinary libm transcendental,
ordinary eigensolver, or SVD is used in the enclosure. Domain subdivision only
reduces interval dependency; every cell is widened before proof arithmetic.

``SeaStateFusionFilter_OU_III::updateTime`` currently accepts any positive finite
caller ``dt``. This stage therefore certifies the nominal 200 Hz configured
schedule only. The source-domain contract records that scope explicitly.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_interval as IA
from ou3_interval import Interval, hull
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_HEADER = SOURCE.DEFAULT_HEADER
CORE_MATH = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
SCHEMA = 2
BRANCH_X = 1.0e-2
CELLS_PER_BRANCH = 64


def _I(value: float) -> Interval:
    return Interval.outward_bounds(float(value), float(value))


def _small_pa_kernel(x: Interval) -> Interval:
    x2 = x.square()
    x3 = x2 * x
    x4 = x3 * x
    return x2 * _I(0.5) - x3 * _I(1.0 / 6.0) + x4 * _I(1.0 / 24.0)


def _small_Sa_kernel(x: Interval) -> Interval:
    x2 = x.square()
    x3 = x2 * x
    x4 = x3 * x
    x5 = x4 * x
    return x3 * _I(1.0 / 6.0) - x4 * _I(1.0 / 24.0) + x5 * _I(1.0 / 120.0)


def _geometric_edges(lo: float, hi: float, count: int) -> list[float]:
    if not (0.0 < lo < hi) or count < 1:
        raise ValueError(f"bad geometric partition [{lo}, {hi}] / {count}")
    ratio = (hi / lo) ** (1.0 / count)
    edges = [lo]
    value = lo
    for _ in range(1, count):
        value *= ratio
        if not (edges[-1] < value < hi):
            value = math.nextafter(edges[-1], math.inf)
        edges.append(min(value, math.nextafter(hi, -math.inf)))
    edges.append(hi)
    for a, b in zip(edges, edges[1:]):
        if not a < b:
            raise RuntimeError("partition failed to make strictly increasing cuts")
    return edges


def _partition(lo: float, hi: float, count: int) -> list[Interval]:
    if math.isclose(lo, hi, rel_tol=0.0, abs_tol=0.0):
        return [Interval.outward_bounds(lo, hi)]
    edges = _geometric_edges(lo, hi, count)
    return [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def _serialize_matrix(A: IA.IntervalMatrix) -> list[list[list[float]]]:
    return [[x.as_list() for x in row] for row in A]


def _axis_transition(h: Interval, tau: Interval, alpha: Interval,
                     em1: Interval, phi_pa: Interval,
                     phi_Sa: Interval) -> IA.IntervalMatrix:
    z = Interval.point(0.0)
    one = Interval.point(1.0)
    phi_va = -(tau * em1)
    half_h2 = h.square() * _I(0.5)
    return [
        [one, z, z, phi_va],
        [h, one, z, phi_pa],
        [half_h2, h, one, phi_Sa],
        [z, z, z, alpha],
    ]


def _cell(x: Interval, h: Interval, branch: str) -> dict:
    if x.lo <= 0.0:
        raise ValueError("h/tau cell must stay strictly positive")
    if x.hi > VT.MAX_ABS_ARGUMENT:
        raise ValueError(f"h/tau={x.as_list()} exceeds audited range {VT.MAX_ABS_ARGUMENT}")

    tau = h / x
    alpha = VT.exp_interval(-x)
    em1 = VT.expm1_interval(-x)
    if branch == "small_x_polynomial":
        pa_kernel = _small_pa_kernel(x)
        Sa_kernel = _small_Sa_kernel(x)
    elif branch == "expm1":
        pa_kernel = VT.ou_phi_pa_kernel_interval(x)
        Sa_kernel = VT.ou_phi_Sa_kernel_interval(x)
    elif branch == "threshold_hull":
        pa_kernel = hull(_small_pa_kernel(x), VT.ou_phi_pa_kernel_interval(x))
        Sa_kernel = hull(_small_Sa_kernel(x), VT.ou_phi_Sa_kernel_interval(x))
    else:
        raise ValueError(f"unknown source branch {branch!r}")

    tau2 = tau.square()
    tau3 = tau2 * tau
    phi_pa = tau2 * pa_kernel
    phi_Sa = tau3 * Sa_kernel
    F = _axis_transition(h, tau, alpha, em1, phi_pa, phi_Sa)
    return {
        "branch": branch,
        "x_h_over_tau": x.as_list(),
        "tau_aw_s": tau.as_list(),
        "alpha": alpha.as_list(),
        "em1": em1.as_list(),
        "phi_va_s": (-(tau * em1)).as_list(),
        "phi_pa_s2": phi_pa.as_list(),
        "phi_Sa_s3": phi_Sa.as_list(),
        "transition_axis_v_p_S_aw": _serialize_matrix(F),
        "transition_axis_spectral_norm_upper": IA.matrix_spectral_norm_upper(F),
    }


def _as_interval(bounds: list[float]) -> Interval:
    if not isinstance(bounds, list) or len(bounds) != 2:
        raise ValueError(f"expected two interval endpoints, got {bounds!r}")
    return Interval(float(bounds[0]), float(bounds[1]))


def _hull_key(cells: list[dict], key: str) -> list[float]:
    return hull(*[_as_interval(c[key]) for c in cells]).as_list()


def build(header: Path = DEFAULT_HEADER.resolve(), cells_per_branch: int = CELLS_PER_BRANCH) -> dict:
    header = header.resolve()
    source = SOURCE.build(header)
    validated_box = source["validated_parameter_box"]
    tau = _as_interval(validated_box["continuous_parameters"]["tau_aw_s"])
    runtime = source["configured_runtime_assumption"]
    nominal_dt = float(runtime["imu_dt_s"])
    h = _as_interval(runtime["imu_dt_outward_interval_s"])
    x_all = h / tau
    if x_all.lo <= 0.0 or x_all.hi > VT.MAX_ABS_ARGUMENT:
        raise RuntimeError(f"source h/tau box outside audited range: {x_all.as_list()}")

    cells: list[dict] = []
    if x_all.lo < BRANCH_X:
        small_hi = min(x_all.hi, BRANCH_X)
        if x_all.lo < small_hi:
            small_cells = _partition(x_all.lo, small_hi, cells_per_branch)
            for i, x in enumerate(small_cells):
                branch = "threshold_hull" if i == len(small_cells) - 1 and x.hi >= BRANCH_X else "small_x_polynomial"
                cells.append(_cell(x, h, branch))
    if x_all.hi >= BRANCH_X:
        large_lo = max(x_all.lo, BRANCH_X)
        if large_lo < x_all.hi:
            cells.extend(_cell(x, h, "expm1") for x in _partition(large_lo, x_all.hi, cells_per_branch))
        elif large_lo == x_all.hi:
            cells.append(_cell(Interval.outward_bounds(large_lo, x_all.hi), h, "expm1"))
    if not cells:
        raise RuntimeError("source box produced no scalar OU cells")

    return {
        "schema": SCHEMA,
        "qualification": "VALIDATED_CONFIGURED_ONE_STEP_OU_AXIS_TRANSITION_ENCLOSURE",
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "transcendental_backend": "EXACT_RATIONAL_TAYLOR_REMAINDER_PLUS_BINARY64_NEXTAFTER_OUTWARD",
        "matrix_backend": "OUTWARD_INTERVAL_BASIC_OPS_PLUS_ABSOLUTE_SUM_NORM_BOUND",
        "implementation_header": str(header.relative_to(REPO)),
        "implementation_core_math": str(CORE_MATH.relative_to(REPO)),
        "source_parameter_box_qualification": validated_box["qualification"],
        "configured_runtime_assumption": runtime,
        "tau_aw_source_box_s": tau.as_list(),
        "nominal_imu_dt_source_value_s": nominal_dt,
        "nominal_imu_dt_box_s": h.as_list(),
        "x_h_over_tau_source_box": x_all.as_list(),
        "implementation_branch_threshold_x": BRANCH_X,
        "cell_count": len(cells),
        "cells_per_regular_branch": cells_per_branch,
        "cells": cells,
        "global_bounds": {
            "alpha": _hull_key(cells, "alpha"),
            "em1": _hull_key(cells, "em1"),
            "phi_va_s": _hull_key(cells, "phi_va_s"),
            "phi_pa_s2": _hull_key(cells, "phi_pa_s2"),
            "phi_Sa_s3": _hull_key(cells, "phi_Sa_s3"),
            "transition_axis_spectral_norm_upper": max(
                float(c["transition_axis_spectral_norm_upper"]) for c in cells
            ),
        },
        "one_step_scalar_ou_enclosed": True,
        "one_step_axis_transition_enclosed": True,
        "deployment_timing_complete": bool(runtime["sample_period_contract"] == "FIXED_SOURCE_NOMINAL"),
        "deployment_timing_scope": runtime["theorem_scope_note"],
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate the validated one-step prediction/correction/reset matrices and Riccati "
            "covariance over complete H/A words under the explicit PE/gating operating envelope"
        ),
    }


def validate(payload: dict) -> list[str]:
    failures: list[str] = []
    for key, expected in (
        ("schema", SCHEMA),
        ("validated_arithmetic", True),
        ("outward_rounded", True),
        ("one_step_scalar_ou_enclosed", True),
        ("one_step_axis_transition_enclosed", True),
        ("deployment_timing_complete", True),
        ("continuous_word_enclosed", False),
        ("nonlinear_word_enclosed", False),
        ("theorem_promotion", "NOT_ESTABLISHED"),
    ):
        if payload.get(key) != expected:
            failures.append(f"{key} mismatch")
    x = _as_interval(payload.get("x_h_over_tau_source_box", [math.nan, math.nan]))
    if not (x.lo > 0.0 and x.hi <= VT.MAX_ABS_ARGUMENT):
        failures.append("h/tau source box is outside audited range")
    cells = payload.get("cells")
    if not isinstance(cells, list) or not cells:
        failures.append("no proof cells")
        return failures
    xs = sorted((_as_interval(c["x_h_over_tau"]) for c in cells), key=lambda I: (I.lo, I.hi))
    if xs[0].lo > x.lo or xs[-1].hi < x.hi:
        failures.append("proof cells do not cover source h/tau endpoints")
    for a, b in zip(xs, xs[1:]):
        if a.hi < b.lo:
            failures.append(f"gap in h/tau proof cells: {a.hi} < {b.lo}")
            break
    g = payload.get("global_bounds", {})
    for key in ("alpha", "em1", "phi_va_s", "phi_pa_s2", "phi_Sa_s3"):
        try:
            I = _as_interval(g[key])
        except (KeyError, TypeError, ValueError):
            failures.append(f"missing/invalid global interval {key}")
            continue
        if any(not I.contains_interval(_as_interval(c[key])) for c in cells):
            failures.append(f"global {key} does not contain a proof cell")
    if not failures:
        alpha = _as_interval(g["alpha"])
        em1 = _as_interval(g["em1"])
        va = _as_interval(g["phi_va_s"])
        pa = _as_interval(g["phi_pa_s2"])
        Sa = _as_interval(g["phi_Sa_s3"])
        if not (0.0 < alpha.lo <= alpha.hi <= 1.0): failures.append("alpha sign/range invalid")
        if not (-1.0 < em1.lo <= em1.hi < 0.0): failures.append("em1 sign/range invalid")
        if va.lo < 0.0 or pa.lo < 0.0 or Sa.lo < 0.0:
            failures.append("OU integral transition coefficient enclosure is negative")
        norm = float(g.get("transition_axis_spectral_norm_upper", math.nan))
        if not math.isfinite(norm) or norm < 1.0:
            failures.append("axis transition spectral-norm upper bound is invalid")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    ap.add_argument("--cells-per-branch", type=int, default=CELLS_PER_BRANCH)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    payload = build(args.header.resolve(), args.cells_per_branch)
    failures = validate(payload)
    payload["validation_pass"] = not failures
    payload["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": payload["qualification"],
        "validation_pass": payload["validation_pass"],
        "cell_count": payload["cell_count"],
        "x_h_over_tau_source_box": payload["x_h_over_tau_source_box"],
        "global_bounds": payload["global_bounds"],
        "deployment_timing_complete": payload["deployment_timing_complete"],
        "theorem_promotion": payload["theorem_promotion"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
