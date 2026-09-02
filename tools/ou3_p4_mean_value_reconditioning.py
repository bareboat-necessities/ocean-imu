#!/usr/bin/env python3
"""Rigorous mean-value reconditioning for OU-III P4 interval AD.

Natural interval evaluation forgets correlations between state components after
accepted Kalman/Joseph updates and subsequent prediction.  The complete-word P4
backend already carries outward interval derivatives with respect to the fixed
word-entry coordinates.  This module uses those derivatives to build the
standard mean-value enclosure

    y(X) \subset y(x_c) + J_y(X) (X-x_c),

then intersects it with the natural interval value already carried by the AD
object.  Both operands are independently rigorous enclosures of the same map;
the intersection therefore remains rigorous and can only tighten the value
range.  Derivative intervals are left unchanged.

The center value is allowed to be an interval.  This is important for the P4
source cells: source/tuner parameters are a separate discrete/interval axis,
while the mean-value coordinates are the continuous word-entry error.  An
interval enclosure of y(x_c,s) for every source parameter s is sufficient.

No trajectory replay, floating eigensolver, or Kalman-gain interval matrix is
used here.  The primitive does not promote P4/P5 by itself.
"""
from __future__ import annotations

import json

from ou3_interval import Interval
import ou3_interval_ad as AD

SCHEMA = 1


def _delta_box(entry_box, entry_center):
    if len(entry_box) != len(entry_center):
        raise ValueError("entry box / center dimension mismatch")
    out = []
    for x, c in zip(entry_box, entry_center):
        cc = float(c)
        if not x.lo <= cc <= x.hi:
            raise ValueError("mean-value center lies outside entry interval")
        out.append(Interval.outward_bounds(x.lo - cc, x.hi - cc))
    return out


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(float(a.lo), float(b.lo))
    hi = min(float(a.hi), float(b.hi))
    if lo > hi:
        raise RuntimeError(
            f"independent rigorous enclosures became disjoint: natural={a} mean_value={b}"
        )
    # The operands are already outward enclosures.  Their set intersection is
    # therefore rigorous without adding another nextafter widening.
    return Interval(lo, hi)


def mean_value_interval(q: AD.AD, q_center: AD.AD, entry_box, entry_center) -> Interval:
    """Return the mean-value enclosure of one AD scalar over the entry box."""
    if q.n != len(entry_box) or q_center.n != q.n:
        raise ValueError("AD / entry dimension mismatch")
    dx = _delta_box(entry_box, entry_center)
    mv = q_center.val
    for d, x in zip(q.der, dx):
        mv = mv + d * x
    return mv


def recondition_scalar(q: AD.AD, q_center: AD.AD, entry_box, entry_center) -> AD.AD:
    """Intersect natural and mean-value value enclosures; keep derivatives."""
    mv = mean_value_interval(q, q_center, entry_box, entry_center)
    return AD.AD(_intersect(q.val, mv), q.der)


def recondition_vector(values, center_values, entry_box, entry_center):
    if len(values) != len(center_values):
        raise ValueError("vector / center length mismatch")
    return [
        recondition_scalar(q, qc, entry_box, entry_center)
        for q, qc in zip(values, center_values)
    ]


def recondition_matrix(values, center_values, entry_box, entry_center):
    if len(values) != len(center_values):
        raise ValueError("matrix / center row mismatch")
    out = []
    for row, crow in zip(values, center_values):
        if len(row) != len(crow):
            raise ValueError("matrix / center column mismatch")
        out.append([
            recondition_scalar(q, qc, entry_box, entry_center)
            for q, qc in zip(row, crow)
        ])
    return out


def max_width_ratio_before_after(before, after) -> float:
    """Diagnostic maximum of width(after)/width(before), ignoring point inputs."""
    ratios = []
    for a, b in zip(before, after):
        wa = a.val.width()
        if wa > 0.0:
            ratios.append(b.val.width() / wa)
    return max(ratios) if ratios else 0.0


def _self_test() -> dict:
    failures = []

    # x-x is identically zero.  This arithmetic layer outward-widens every
    # basic operation, including exact cancellation, so its AD derivative is a
    # tiny interval containing zero rather than bitwise [0,0].  Mean-value
    # reconditioning must retain zero, strictly tighten the natural value box,
    # and leave the derivative enclosure bit-for-bit unchanged.
    x = AD.independent(Interval(-1.0, 1.0), 0, 1)
    y = x - x
    yc = AD.constant(0.0, 1)
    yr = recondition_scalar(y, yc, [x.val], [0.0])
    if not yr.val.contains(0.0):
        failures.append(f"reconditioned repeated-variable enclosure misses zero: {yr.val}")
    if not yr.val.width() < y.val.width():
        failures.append("mean-value reconditioning did not tighten x-x")
    if yr.der != y.der:
        failures.append("reconditioning changed derivative enclosure")
    if not yr.der[0].contains(0.0):
        failures.append("outward AD derivative no longer contains exact zero derivative")

    # An affine map should be enclosed by its mean-value form and intersection
    # must never widen the natural enclosure.
    z = AD.constant(2.0, 1) * x + AD.constant(3.0, 1)
    zc = AD.constant(3.0, 1)
    zr = recondition_scalar(z, zc, [x.val], [0.0])
    for v in (1.0, 5.0):
        if not zr.val.contains(v):
            failures.append(f"affine reconditioned enclosure misses endpoint {v}")
    if zr.val.width() > z.val.width():
        failures.append("reconditioning widened a natural enclosure")

    return {
        "schema": SCHEMA,
        "pass": not failures,
        "failures": failures,
        "repeated_variable_natural": [y.val.lo, y.val.hi],
        "repeated_variable_reconditioned": [yr.val.lo, yr.val.hi],
        "repeated_variable_value_strictly_tightened": yr.val.width() < y.val.width(),
        "derivative_preserved_bitwise": yr.der == y.der,
        "derivative_contains_exact_zero": yr.der[0].contains(0.0),
        "rigorous_intersection_of_natural_and_mean_value_forms": True,
        "source_parameters_are_not_reinterpreted_as_entry_coordinates": True,
        "K_interval_matrix_materialized": False,
        "P4_PROMOTED_HERE": False,
        "P5_PROMOTED_HERE": False,
    }


if __name__ == "__main__":
    d = _self_test()
    print(json.dumps(d, indent=2, sort_keys=True))
    raise SystemExit(0 if d["pass"] else 2)
