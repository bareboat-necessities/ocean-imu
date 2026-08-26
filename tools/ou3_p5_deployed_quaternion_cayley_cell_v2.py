#!/usr/bin/env python3
"""Radially subdivided exact deployed-quaternion/Cayley primitive for P5.

The original outer primitive is rigorous through 6 rad because on its complete
axis-angle branch ||d||/2 < 3 < pi, so cos and sinc are monotone.  The sample-1
source witness now reaches corrections above 6 rad.  The shipping filter itself
has no 6-rad clamp: it evaluates

    q_d ~ [ cos(||d||/2), sin(||d||/2)/||d|| d ]

and normalizes the quaternion.  Past 2*pi the half-angle passes pi, so simply
raising the old constant would invalidate its monotonicity argument.

This module extends the proof primitive to 9 rad by requiring a source-certified
radial subcell [d_lo,d_hi] whenever d_hi>6.  On one such cell the half-angle is
centered at m with |h|<=delta and uses validated point sin/cos at m together with

    |sin h| <= |h|,
    1-h^2/2 <= cos h <= 1.

The point values are enclosed by exact-rational Taylor polynomials with a
Lagrange remainder.  Thus no monotonicity of sin, cos, or sinc is used across a
winding.  The homogeneous quaternion is composed directly with the pre-error
Cayley quaternion [2,c]; only a zero crossing of the *resulting* homogeneous
scalar is rejected.

Cells above 6 rad with no positive radial lower bound fail closed and must be
subdivided.  Existing <=6-rad calls delegate to the already validated V1
primitive unchanged.  This file changes no filter and promotes no P5 word.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_p5_deployed_quaternion_cayley_cell as V1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 2
MAX_CORRECTION_NORM = 9.0
MAX_HALF_ANGLE = 0.5 * MAX_CORRECTION_NORM
TRIG_ORDER = 42


def I(x: float) -> Interval:
    return Interval.point(float(x))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _down_fraction(q: Fraction) -> float:
    f = float(q)
    if Fraction.from_float(f) > q:
        f = math.nextafter(f, -math.inf)
    return f


def _up_fraction(q: Fraction) -> float:
    f = float(q)
    if Fraction.from_float(f) < q:
        f = math.nextafter(f, math.inf)
    return f


def _sin_point(x: float) -> Interval:
    x = float(x)
    if not (math.isfinite(x) and abs(x) <= MAX_HALF_ANGLE):
        raise ValueError("extended quaternion half-angle outside validated range")
    q = Fraction.from_float(x)
    total = Fraction(0, 1)
    for k in range(1, TRIG_ORDER + 1, 2):
        n = (k - 1) // 2
        total += (-1 if n & 1 else 1) * q**k / Fraction(math.factorial(k), 1)
    rem = abs(q) ** (TRIG_ORDER + 1) / Fraction(math.factorial(TRIG_ORDER + 1), 1)
    return Interval(_down_fraction(total - rem), _up_fraction(total + rem))


def _cos_point(x: float) -> Interval:
    x = float(x)
    if not (math.isfinite(x) and abs(x) <= MAX_HALF_ANGLE):
        raise ValueError("extended quaternion half-angle outside validated range")
    q = Fraction.from_float(x)
    total = Fraction(0, 1)
    for k in range(0, TRIG_ORDER + 1, 2):
        n = k // 2
        total += (-1 if n & 1 else 1) * q**k / Fraction(math.factorial(k), 1)
    rem = abs(q) ** (TRIG_ORDER + 1) / Fraction(math.factorial(TRIG_ORDER + 1), 1)
    return Interval(_down_fraction(total - rem), _up_fraction(total + rem))


def _trig_interval(lo: float, hi: float) -> tuple[Interval, Interval]:
    """Validated sin/cos enclosure on one nonnegative half-angle radial cell."""
    lo, hi = float(lo), float(hi)
    if not (0.0 <= lo <= hi <= MAX_HALF_ANGLE):
        raise ValueError("invalid extended half-angle interval")
    mid = 0.5 * (lo + hi)
    # Outward radius covers binary midpoint rounding as well as both endpoints.
    delta = up(max(abs(hi - mid), abs(mid - lo)))
    sm = _sin_point(mid)
    cm = _cos_point(mid)
    sh = Interval(down(-delta), up(delta))  # |sin h| <= |h|
    ch_lo = down(1.0 - up(0.5 * up(delta * delta)))
    ch = Interval(max(-1.0, ch_lo), 1.0)   # cos h >= 1-h^2/2
    sinx = sm * ch + cm * sh
    cosx = cm * ch - sm * sh
    # Universal range intersection; the raw enclosure already contains truth.
    sinx = Interval(max(-1.0, sinx.lo), min(1.0, sinx.hi))
    cosx = Interval(max(-1.0, cosx.lo), min(1.0, cosx.hi))
    return sinx, cosx


def _check_vec(v: Sequence[Interval], label: str) -> list[Interval]:
    out = list(v)
    if len(out) != 3 or not all(isinstance(x, Interval) for x in out):
        raise ValueError(f"{label} must be a three-vector of intervals")
    return out


def _norm_lower(v: Sequence[Interval]) -> float:
    v = _check_vec(v, "v")
    s = 0.0
    for x in v:
        if x.lo <= 0.0 <= x.hi:
            a = 0.0
        else:
            a = min(abs(x.lo), abs(x.hi))
        term = down(a * a) if a > 0.0 else 0.0
        s = down(s + term) if term > 0.0 else s
    return down(math.sqrt(max(0.0, s))) if s > 0.0 else 0.0


def correction_homogeneous_quaternion(
    d: Sequence[Interval], *, d_norm_lower: float | None = None,
    d_norm_upper: float | None = None,
):
    d = _check_vec(d, "d")
    box_hi = V1._norm_upper(d)
    hi = box_hi if d_norm_upper is None else float(d_norm_upper)
    if not (math.isfinite(hi) and 0.0 <= hi <= MAX_CORRECTION_NORM):
        raise ValueError(f"deployed correction norm upper outside validated range [0,{MAX_CORRECTION_NORM:g}]")
    if hi > box_hi + 64.0 * math.ulp(max(1.0, box_hi)):
        raise ValueError("supplied correction norm upper exceeds component-box bound")
    if hi <= V1.MAX_CORRECTION_NORM:
        w, v, dn, _ = V1.correction_homogeneous_quaternion(d, d_norm_upper=hi)
        return w, v, 0.0, dn, box_hi, "V1_MONOTONE_THROUGH_6_RAD"

    box_lo = _norm_lower(d)
    lo = box_lo if d_norm_lower is None else float(d_norm_lower)
    if not (math.isfinite(lo) and 0.0 < lo <= hi):
        raise ValueError("correction cells above 6 rad require a positive radial lower bound")
    # A supplied lower bound must not exceed the component-box guaranteed norm.
    if lo > box_lo + 64.0 * math.ulp(max(1.0, box_lo)):
        raise ValueError("supplied correction norm lower is not certified by the component box")

    half = Interval(down(0.5 * lo), up(0.5 * hi))
    sin_half, cos_half = _trig_interval(half.lo, half.hi)
    # k=sin(theta/2)/theta = 0.5*sin(x)/x, x=theta/2>0.
    k = I(0.5) * sin_half / half
    v = [k * x for x in d]
    return cos_half, v, lo, hi, box_hi, "RADIAL_SUBCELL_NONMONOTONE_TRIG"


def compose_cell(
    c: Sequence[Interval], d: Sequence[Interval], *,
    c_norm_upper: float | None = None,
    d_norm_lower: float | None = None,
    d_norm_upper: float | None = None,
) -> dict:
    c = _check_vec(c, "c")
    d = _check_vec(d, "d")
    cbox = V1._norm_upper(c)
    cn = cbox if c_norm_upper is None else float(c_norm_upper)
    if not (math.isfinite(cn) and cn >= 0.0 and cn <= cbox + 64.0 * math.ulp(max(1.0, cbox))):
        raise ValueError("supplied Cayley norm upper is invalid for component cell")

    wd, vd, dlo, dhi, dbox, backend = correction_homogeneous_quaternion(
        d, d_norm_lower=d_norm_lower, d_norm_upper=d_norm_upper
    )
    vdotc = V1._dot(vd, c)
    W = I(2.0) * wd - vdotc
    if W.lo <= 0.0 <= W.hi:
        raise RuntimeError("resulting deployed quaternion product can reach Cayley antipode")
    vxc = V1._cross(vd, c)
    V = [wd*c[i] + I(2.0)*vd[i] + vxc[i] for i in range(3)]
    cp = [I(2.0)*x/W for x in V]
    return {
        "c_plus": cp,
        "product_scalar": W,
        "signed_vd_dot_c": vdotc,
        "correction_norm_lower": dlo,
        "correction_norm_upper": dhi,
        "correction_component_box_norm_upper": dbox,
        "c_norm_upper_input": cn,
        "c_plus_norm_upper": V1._norm_upper(cp),
        "quaternion_enclosure_backend": backend,
        "radial_subdivision_required_above_6_rad": True,
        "correction_cayley_coordinate_formed": False,
        "source_quaternion_normalization_cancels_exactly": True,
    }


def _serial(r: dict) -> dict:
    return {
        "product_scalar": r["product_scalar"].as_list(),
        "correction_norm_lower": r["correction_norm_lower"],
        "correction_norm_upper": r["correction_norm_upper"],
        "c_plus_norm_upper": r["c_plus_norm_upper"],
        "quaternion_enclosure_backend": r["quaternion_enclosure_backend"],
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).resolve().read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("extended deployed quaternion primitive must not be trajectory fitted")
    base = V1.build(domain_path)
    failures = [f"V1: {x}" for x in V1.validate(base)]

    # Exercises a correction beyond 2*pi and beyond the old 6-rad range.
    extended = compose_cell(
        [I(0.0), I(0.0), I(0.0)],
        [Interval.outward_bounds(8.0, 8.2), I(0.0), I(0.0)],
        d_norm_lower=down(8.0), d_norm_upper=up(8.2),
    )
    if not extended["correction_norm_lower"] > 2.0 * math.pi:
        failures.append("extended audit did not exceed 2*pi")
    if not extended["correction_norm_upper"] > V1.MAX_CORRECTION_NORM:
        failures.append("extended audit did not exceed old 6-rad proof range")
    if extended["product_scalar"].lo <= 0.0 <= extended["product_scalar"].hi:
        failures.append("extended audit reaches resulting Cayley antipode")
    if extended["quaternion_enclosure_backend"] != "RADIAL_SUBCELL_NONMONOTONE_TRIG":
        failures.append("extended audit did not use radial nonmonotone backend")

    # Fail-closed audit: a symmetric >6-rad component box has no positive
    # radial lower bound and must be subdivided before proof use.
    broad_refused = False
    try:
        compose_cell(
            [I(0.0), I(0.0), I(0.0)],
            [Interval.outward_bounds(-8.2, 8.2), I(0.0), I(0.0)],
            d_norm_upper=up(8.2),
        )
    except ValueError:
        broad_refused = True
    if not broad_refused:
        failures.append("unsubdivided winding correction box was not refused")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_RADIAL_SUBCELL_DEPLOYED_QUATERNION_TO_RESULT_CAYLEY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "maximum_validated_correction_norm_rad": MAX_CORRECTION_NORM,
        "old_six_rad_monotonicity_not_extended_past_its_domain": True,
        "radial_subdivision_required_above_6_rad": True,
        "nonmonotone_half_angle_trig_enclosed_without_monotonicity": True,
        "only_resulting_error_antipode_is_gate": True,
        "unsubdivided_winding_box_refused": broad_refused,
        "complete_word_promoted_here": False,
        "audit_extended_over_2pi": _serial(extended),
        "P5_DEPLOYED_QUATERNION_CAYLEY_CELL_V2_PRIMITIVE": "PASS" if not failures else "FAIL",
        "next_obligation": "feed source-correlated radial correction subcells into V2 at sample 1, then add reset/process/tangent-force perturbations",
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "old_six_rad_monotonicity_not_extended_past_its_domain",
        "radial_subdivision_required_above_6_rad",
        "nonmonotone_half_angle_trig_enclosed_without_monotonicity",
        "only_resulting_error_antipode_is_gate",
        "unsubdivided_winding_box_refused",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "complete_word_promoted_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("maximum_validated_correction_norm_rad", 0.0)) < 9.0:
        f.append("extended correction range below 9 rad")
    a = d.get("audit_extended_over_2pi", {})
    if not float(a.get("correction_norm_lower", 0.0)) > 2.0 * math.pi:
        f.append("extended audit lower is not beyond 2*pi")
    W = a.get("product_scalar", [0.0, 0.0])
    if float(W[0]) <= 0.0 <= float(W[1]):
        f.append("extended audit product scalar crosses zero")
    if not f and d.get("P5_DEPLOYED_QUATERNION_CAYLEY_CELL_V2_PRIMITIVE") != "PASS":
        f.append("V2 primitive did not pass")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_DEPLOYED_QUATERNION_CAYLEY_CELL_V2_PRIMITIVE"],
        "max_rad": out["maximum_validated_correction_norm_rad"],
        "extended": out["audit_extended_over_2pi"],
        "broad_refused": out["unsubdivided_winding_box_refused"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
