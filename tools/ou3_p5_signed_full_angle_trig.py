#!/usr/bin/env python3
"""Validated signed full-angle interval trigonometry for OU-III P5 proof gauges.

This helper is deliberately separate from the deployed-quaternion half-angle
radial primitive.  For a signed interval X=[lo,hi], write x=m+h and combine the
validated exact-rational point sin/cos enclosures at m with global elementary
bounds on sin(h), cos(h).  No range reduction or libm call is used for proof
arithmetic.  Midpoints outside the audited point range fail wide to [-1,1].
"""
from __future__ import annotations

import math

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14

FULL = V14.FULL
AUDITED_POINT_ABS_MAX = 4.5


def signed_full_angle_trig_interval(x: Interval) -> tuple[Interval, Interval, bool]:
    """Return rigorous (sin(X), cos(X), broad_fallback_used) for signed X."""
    lo = float(x.lo)
    hi = float(x.hi)
    if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
        raise ValueError("signed proof-gauge angle interval must be finite and ordered")

    mid = 0.5 * (lo + hi)
    if not math.isfinite(mid) or abs(mid) > AUDITED_POINT_ABS_MAX:
        whole = Interval(-1.0, 1.0)
        return whole, whole, True

    delta = FULL.up(max(abs(hi - mid), abs(mid - lo)))
    sm = V14.CAYLEY2._sin_point(mid)
    cm = V14.CAYLEY2._cos_point(mid)

    srad = min(1.0, delta)
    sh = Interval(FULL.down(-srad), FULL.up(srad))
    d2 = FULL.up(delta * delta)
    ch_lo = max(-1.0, FULL.down(1.0 - FULL.up(0.5 * d2)))
    ch = Interval(ch_lo, 1.0)

    sinx = sm * ch + cm * sh
    cosx = cm * ch - sm * sh
    sinx = Interval(max(-1.0, sinx.lo), min(1.0, sinx.hi))
    cosx = Interval(max(-1.0, cosx.lo), min(1.0, cosx.hi))
    return sinx, cosx, False
