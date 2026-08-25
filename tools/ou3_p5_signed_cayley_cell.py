#!/usr/bin/env python3
"""Signed outward Cayley correction cells for the OU-III P5 enclosure.

Norm-only correction bounds replace the exact denominator
``1-a^T c/4`` by ``1-|a||c|/4``.  That is safe for chart exclusion, but it loses
the most important source correlation of a corrective vector update: the
signed correction/Cayley inner product.  The complete P5 word backend therefore
needs a componentwise interval primitive.

For a source correction error-state vector ``d`` the deployed quaternion helper
produces a correction Cayley vector ``a=s(|d|)d``.  The normalization in the C++
helper cancels from ``2 v/w``.  Below ``1e-2`` the source polynomial gives

    s(x) = 2(1/2-x^2/48+x^4/3840)/(1-x^2/8+x^4/384),

and above it ``s(x)=2 tan(x/2)/x``.  Both branches are enclosed outward here.
For finite component boxes ``C,D`` the exact Cayley composition is then

    c+ = (a+c+0.5 a x c)/(1-0.25 a^T c).

The denominator is evaluated with the signed interval dot product, not an
independent norm product.  This is the primitive required by source-correlated
subdivision of the later vector/S prefixes.  It neither assumes a favorable
correction direction nor promotes P5 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, down, hull, up
import ou3_p4_group_algebra as GROUP
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _check_vec(v: Sequence[Interval], label: str) -> list[Interval]:
    out = list(v)
    if len(out) != 3 or not all(isinstance(x, Interval) for x in out):
        raise ValueError(f"{label} must be a 3-vector of intervals")
    return out


def vector_norm_upper(v: Sequence[Interval]) -> float:
    v = _check_vec(v, "v")
    s = 0.0
    for x in v:
        a = x.abs_upper()
        s = up(s + up(a * a))
    return up(math.sqrt(s))


def dot(a: Sequence[Interval], b: Sequence[Interval]) -> Interval:
    a, b = _check_vec(a, "a"), _check_vec(b, "b")
    y = I(0.0)
    for x, z in zip(a, b):
        y = y + x * z
    return y


def cross(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    a, b = _check_vec(a, "a"), _check_vec(b, "b")
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]


def _series_coeff_interval(x_hi: float) -> Interval:
    """Enclose 2k/w on the source polynomial branch, 0<=x<=x_hi<=1e-2."""
    if not (0.0 <= x_hi <= GROUP.SERIES_BRANCH_NORM):
        raise ValueError("series coefficient range invalid")
    x2 = Interval.outward_bounds(0.0, up(x_hi*x_hi))
    x4 = x2.square()
    w = I(1.0) - I(1.0/8.0)*x2 + I(1.0/384.0)*x4
    k = I(0.5) - I(1.0/48.0)*x2 + I(1.0/3840.0)*x4
    if not w.lo > 0.0 or not k.lo > 0.0:
        raise RuntimeError("source polynomial correction coefficient lost positivity")
    return I(2.0)*k/w


def _axis_coeff_point_bounds(x: float) -> Interval:
    if not (GROUP.SERIES_BRANCH_NORM <= x <= GROUP.CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError("axis coefficient point outside validated range")
    half = up(0.5*x)
    s = VT.sin_point(half)
    c = VT.cos_point(half)
    if not c.lo > 0.0:
        raise RuntimeError("correction quaternion scalar interval crosses zero")
    # 2 sin(x/2)/(x cos(x/2)); all factors are positive here.
    lo = down((2.0*s.lo)/(up(x*c.hi)))
    hi = up((2.0*s.hi)/(down(x*c.lo)))
    return Interval(lo, hi)


def correction_cayley_scale_interval(d_norm_upper: float) -> Interval:
    """Enclose positive a/d scale for every deployed correction with |d|<=upper."""
    d_norm_upper = float(d_norm_upper)
    if not (math.isfinite(d_norm_upper) and 0.0 <= d_norm_upper <= GROUP.CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError("correction norm upper outside validated range [0,3]")
    pieces = [_series_coeff_interval(min(d_norm_upper, GROUP.SERIES_BRANCH_NORM))]
    if d_norm_upper >= GROUP.SERIES_BRANCH_NORM:
        # 2 tan(x/2)/x is increasing on (0,pi).  The endpoint enclosures thus
        # cover the entire axis-angle source branch.
        lo = _axis_coeff_point_bounds(GROUP.SERIES_BRANCH_NORM).lo
        hi = _axis_coeff_point_bounds(max(d_norm_upper, GROUP.SERIES_BRANCH_NORM)).hi
        pieces.append(Interval(lo, hi))
    s = hull(*pieces)
    if not s.lo > 0.0:
        raise RuntimeError("deployed correction Cayley scale is not positive")
    return s


def correction_cayley_vector(d: Sequence[Interval]) -> tuple[list[Interval], Interval, float]:
    d = _check_vec(d, "d")
    dn = vector_norm_upper(d)
    scale = correction_cayley_scale_interval(dn)
    return [scale*x for x in d], scale, dn


def compose_cell(c: Sequence[Interval], d: Sequence[Interval]) -> dict:
    """Outward exact deployed Cayley composition of one source correction cell."""
    c, d = _check_vec(c, "c"), _check_vec(d, "d")
    a, scale, dn = correction_cayley_vector(d)
    ad = dot(a, c)
    denom = I(1.0) - I(0.25)*ad
    if denom.lo <= 0.0:
        raise RuntimeError("signed correction cell can reach Cayley antipodal denominator")
    axc = cross(a, c)
    num = [a[i] + c[i] + I(0.5)*axc[i] for i in range(3)]
    cp = [x/denom for x in num]
    return {
        "c_plus": cp,
        "a": a,
        "correction_scale": scale,
        "correction_norm_upper": dn,
        "a_dot_c": ad,
        "denominator": denom,
        "c_plus_norm_upper": vector_norm_upper(cp),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).resolve().read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("signed Cayley cell primitive must not be trajectory fitted")

    # Audited cells exercise both source quaternion branches and, critically,
    # demonstrate that a signed corrective dot product is retained.  These are
    # arithmetic self-tests, not source samples used in a theorem claim.
    small = compose_cell(
        [Interval.outward_bounds(0.19, 0.21), Interval.outward_bounds(-0.02, 0.02), Interval.outward_bounds(-0.02, 0.02)],
        [Interval.outward_bounds(-0.006, -0.004), Interval.outward_bounds(-1e-4, 1e-4), Interval.outward_bounds(-1e-4, 1e-4)],
    )
    large = compose_cell(
        [Interval.outward_bounds(0.55, 0.65), Interval.outward_bounds(-0.05, 0.05), Interval.outward_bounds(-0.05, 0.05)],
        [Interval.outward_bounds(-0.45, -0.35), Interval.outward_bounds(-0.02, 0.02), Interval.outward_bounds(-0.02, 0.02)],
    )
    failures = []
    for name, row in (("small", small), ("large", large)):
        if not row["denominator"].lo > 0.0:
            failures.append(f"{name}: denominator not positive")
        if not row["correction_scale"].lo > 0.0:
            failures.append(f"{name}: correction Cayley scale not positive")
    if not small["correction_norm_upper"] < GROUP.SERIES_BRANCH_NORM:
        failures.append("small audit cell did not exercise polynomial branch")
    if not large["correction_norm_upper"] > GROUP.SERIES_BRANCH_NORM:
        failures.append("large audit cell did not exercise axis-angle branch")
    if not large["a_dot_c"].hi < 0.0:
        failures.append("signed corrective audit cell lost negative a dot c")
    if not large["denominator"].lo > 1.0:
        failures.append("signed corrective denominator did not retain >1 widening")

    def serial(row: dict) -> dict:
        return {
            "c_plus": [x.as_list() for x in row["c_plus"]],
            "a": [x.as_list() for x in row["a"]],
            "correction_scale": row["correction_scale"].as_list(),
            "correction_norm_upper": row["correction_norm_upper"],
            "a_dot_c": row["a_dot_c"].as_list(),
            "denominator": row["denominator"].as_list(),
            "c_plus_norm_upper": row["c_plus_norm_upper"],
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SIGNED_OUTWARD_CAYLEY_CORRECTION_CELL_PRIMITIVE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "signed_a_dot_c_retained": True,
        "independent_abs_a_abs_c_denominator_used": False,
        "source_polynomial_and_axis_angle_branches_enclosed": True,
        "complete_word_promoted_here": False,
        "audit_cells": {"small_series": serial(small), "large_corrective": serial(large)},
        "P5_SIGNED_CAYLEY_CELL_PRIMITIVE": "PASS" if not failures else "FAIL",
        "next_obligation": "feed source-correlated interval K*r correction cells into this primitive at every later S/vector prefix",
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in ("source_generated_not_trajectory_fit", "signed_a_dot_c_retained", "source_polynomial_and_axis_angle_branches_enclosed"):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed", "independent_abs_a_abs_c_denominator_used", "complete_word_promoted_here"):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    large = d.get("audit_cells", {}).get("large_corrective", {})
    adot = large.get("a_dot_c", [math.inf, math.inf])
    den = large.get("denominator", [-math.inf, math.inf])
    if not (float(adot[1]) < 0.0 and float(den[0]) > 1.0):
        failures.append("signed corrective audit relation not retained")
    if not failures and d.get("P5_SIGNED_CAYLEY_CELL_PRIMITIVE") != "PASS":
        failures.append("signed Cayley primitive did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SIGNED_CAYLEY_CELL_PRIMITIVE"],
        "large_corrective": out["audit_cells"]["large_corrective"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
