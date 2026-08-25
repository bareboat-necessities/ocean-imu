#!/usr/bin/env python3
"""Exact deployed-quaternion/Cayley composition for large P5 corrections.

The local P4 helper represents the correction itself by a Cayley vector.  That
is ideal near the identity but introduces an artificial singularity when an
otherwise finite deployed correction passes through ``|d|=pi``.  The filter has
no singularity there: it applies a normalized quaternion correction.

For P5 outer capture this module therefore composes the *shipping quaternion*
directly and converts only the resulting physical error rotation back to the
Cayley chart.  A normalization-free homogeneous representation makes the
identity particularly simple.  The pre-correction Cayley coordinate ``c`` is
represented by the quaternion ``[2,c]``.  The deployed correction is represented
by the source's unnormalized quaternion ``[w_d,v_d]`` (normalization cancels in
Cayley coordinates).  Their left product is

    W = 2 w_d - v_d^T c,
    V = w_d c + 2 v_d + v_d x c,
    c+ = 2 V / W.

Thus a correction at ``|d|=pi`` is allowed whenever the *resulting* product has
``W != 0``.  The only promotion failure is the genuine Cayley antipode of the
resulting error rotation, not a coordinate singularity of the correction.

Both shipping correction branches are retained:

* ``|d| < 1e-2``: the source polynomial quaternion coefficients;
* otherwise: axis-angle ``[cos(|d|/2), 0.5 sinc(|d|/2) d]``.

The trusted transcendental backend is used on ``|d|/2 <= 3``.  Hence this outer
primitive validates correction norms through 6 rad, strictly wider than the
3-rad local helper and enough to test the current P5 source family.  It does not
change the filter, use replay, linearize the correction, or replace the signed
product term by ``|a||c|``.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, hull
import ou3_p4_group_algebra as GROUP
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
MAX_CORRECTION_NORM = 6.0


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _check_vec(v: Sequence[Interval], label: str) -> list[Interval]:
    out = list(v)
    if len(out) != 3 or not all(isinstance(x, Interval) for x in out):
        raise ValueError(f"{label} must be a 3-vector of intervals")
    return out


def _norm_upper(v: Sequence[Interval]) -> float:
    """Validated Euclidean-norm upper without manufacturing a negative zero sum."""
    v = _check_vec(v, "v")
    s = 0.0
    for x in v:
        a = x.abs_upper()
        term = a*a
        if term > 0.0:
            term = math.nextafter(term, math.inf)
        s += term
        if s > 0.0:
            s = math.nextafter(s, math.inf)
    if s == 0.0:
        return 0.0
    return math.nextafter(math.sqrt(s), math.inf)


def _dot(a: Sequence[Interval], b: Sequence[Interval]) -> Interval:
    a = _check_vec(a, "a")
    b = _check_vec(b, "b")
    y = I(0.0)
    for x, z in zip(a, b):
        y = y + x*z
    return y


def _cross(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    a = _check_vec(a, "a")
    b = _check_vec(b, "b")
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]


def _intersect(a: Interval, b: Interval) -> Interval | None:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    return None if lo > hi else Interval(lo, hi)


def _series_homogeneous(d: Sequence[Interval], norm_upper: float):
    """Homogeneous source polynomial quaternion on the strict-small branch."""
    hi = min(float(norm_upper), GROUP.SERIES_BRANCH_NORM)
    cap = Interval(-hi, hi)
    ds = []
    for x in d:
        z = _intersect(x, cap)
        if z is None:
            return None
        ds.append(z)
    t2 = Interval(0.0, math.nextafter(hi*hi, math.inf))
    t4 = t2.square()
    w = I(1.0) - I(1.0/8.0)*t2 + I(1.0/384.0)*t4
    k = I(0.5) - I(1.0/48.0)*t2 + I(1.0/3840.0)*t4
    if not w.lo > 0.0 or not k.lo > 0.0:
        raise RuntimeError("source polynomial quaternion lost positive coefficients")
    return w, [k*x for x in ds]


def _axis_homogeneous(d: Sequence[Interval], norm_upper: float):
    """Homogeneous axis-angle source quaternion over the non-small branch."""
    hi = float(norm_upper)
    if hi < GROUP.SERIES_BRANCH_NORM:
        return None
    if not (math.isfinite(hi) and hi <= MAX_CORRECTION_NORM):
        raise ValueError(
            f"deployed correction norm upper outside validated range [0,{MAX_CORRECTION_NORM:g}]"
        )
    # The axis branch starts at 1e-2.  For hi<=6, half-angle <=3<pi, so cos
    # decreases and sinc decreases throughout this complete branch interval.
    half_lo = 0.5*GROUP.SERIES_BRANCH_NORM
    half_hi = 0.5*hi
    clo = VT.cos_point(half_hi)
    chi = VT.cos_point(half_lo)
    w = Interval(clo.lo, chi.hi)
    sinc = VT.sinc_interval(Interval(half_lo, half_hi))
    k = I(0.5)*sinc
    return w, [k*x for x in d]


def _hull_quaternion(parts):
    parts = [x for x in parts if x is not None]
    if not parts:
        raise RuntimeError("no deployed quaternion branch intersects correction cell")
    w = hull(*(x[0] for x in parts))
    v = [hull(*(x[1][i] for x in parts)) for i in range(3)]
    return w, v


def correction_homogeneous_quaternion(
    d: Sequence[Interval], *, d_norm_upper: float | None = None
):
    d = _check_vec(d, "d")
    box_norm = _norm_upper(d)
    dn = box_norm if d_norm_upper is None else float(d_norm_upper)
    if not (math.isfinite(dn) and 0.0 <= dn <= box_norm + 64.0*math.ulp(max(1.0, box_norm))):
        raise ValueError("supplied correction norm upper is invalid for component cell")
    if dn > MAX_CORRECTION_NORM:
        raise ValueError(
            f"deployed correction norm upper outside validated range [0,{MAX_CORRECTION_NORM:g}]"
        )
    parts = [_series_homogeneous(d, dn)]
    if dn >= GROUP.SERIES_BRANCH_NORM:
        parts.append(_axis_homogeneous(d, dn))
    w, v = _hull_quaternion(parts)
    return w, v, dn, box_norm


def compose_cell(
    c: Sequence[Interval],
    d: Sequence[Interval],
    *,
    c_norm_upper: float | None = None,
    d_norm_upper: float | None = None,
) -> dict:
    """Compose the exact deployed correction with one Cayley error cell.

    ``c_norm_upper`` is optional metadata used to preserve a separately proved
    Euclidean radius; the component intervals remain the actual signed operands.
    ``d_norm_upper`` may tighten a component-box norm only when it is itself a
    certified upper bound for the same cell.
    """
    c = _check_vec(c, "c")
    d = _check_vec(d, "d")
    cbox = _norm_upper(c)
    cn = cbox if c_norm_upper is None else float(c_norm_upper)
    if not (math.isfinite(cn) and cn >= 0.0 and cn <= cbox + 64.0*math.ulp(max(1.0, cbox))):
        raise ValueError("supplied Cayley norm upper is invalid for component cell")

    wd, vd, dn, dbox = correction_homogeneous_quaternion(d, d_norm_upper=d_norm_upper)
    vdotc = _dot(vd, c)
    W = I(2.0)*wd - vdotc
    if W.lo <= 0.0 <= W.hi:
        raise RuntimeError("resulting deployed quaternion product can reach Cayley antipode")
    vxc = _cross(vd, c)
    V = [wd*c[i] + I(2.0)*vd[i] + vxc[i] for i in range(3)]
    cp = [I(2.0)*x/W for x in V]
    cpbox = _norm_upper(cp)
    return {
        "c_plus": cp,
        "product_scalar": W,
        # Compatibility name consumed by the existing prefix reporter.  Unlike
        # the local correction-Cayley denominator this is the homogeneous
        # *resulting quaternion scalar*.  Only exclusion of zero is required.
        "denominator": W,
        "correction_homogeneous_scalar": wd,
        "correction_homogeneous_vector": vd,
        "signed_vd_dot_c": vdotc,
        "correction_norm_upper": dn,
        "correction_component_box_norm_upper": dbox,
        "c_norm_upper_input": cn,
        "c_component_box_norm_upper_input": cbox,
        "c_plus_norm_upper": cpbox,
        "correction_cayley_coordinate_formed": False,
        "source_quaternion_normalization_cancels_exactly": True,
    }


def _serial(row: dict) -> dict:
    return {
        "c_plus": [x.as_list() for x in row["c_plus"]],
        "product_scalar": row["product_scalar"].as_list(),
        "signed_vd_dot_c": row["signed_vd_dot_c"].as_list(),
        "correction_norm_upper": row["correction_norm_upper"],
        "c_plus_norm_upper": row["c_plus_norm_upper"],
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).resolve().read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("deployed quaternion P5 cell must not be trajectory fitted")

    small = compose_cell(
        [Interval.outward_bounds(0.19, 0.21), Interval.outward_bounds(-0.01, 0.01), Interval.outward_bounds(-0.01, 0.01)],
        [Interval.outward_bounds(-0.006, -0.004), Interval.outward_bounds(-1e-4, 1e-4), Interval.outward_bounds(-1e-4, 1e-4)],
    )
    # This correction is strictly larger than the retired 3-rad correction-Cayley
    # helper range and crosses pi in the source axis-angle branch.  The chosen
    # corrective pre-error keeps the *product* strictly off the antipode.
    crossing_pi = compose_cell(
        [Interval.outward_bounds(-0.62, -0.58), Interval.outward_bounds(-0.005, 0.005), Interval.outward_bounds(-0.005, 0.005)],
        [Interval.outward_bounds(3.18, 3.22), Interval.outward_bounds(-0.002, 0.002), Interval.outward_bounds(-0.002, 0.002)],
    )
    failures = []
    if not small["product_scalar"].lo > 0.0:
        failures.append("small deployed quaternion audit lost positive product scalar")
    if not crossing_pi["correction_norm_upper"] > 3.0:
        failures.append("large audit did not exceed retired 3-rad helper range")
    if crossing_pi["product_scalar"].lo <= 0.0 <= crossing_pi["product_scalar"].hi:
        failures.append("large correction audit product reaches antipode")
    if crossing_pi["correction_cayley_coordinate_formed"] is not False:
        failures.append("large correction audit formed correction Cayley coordinate")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_DEPLOYED_QUATERNION_TO_RESULT_CAYLEY_CELL",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "maximum_validated_correction_norm_rad": MAX_CORRECTION_NORM,
        "shipping_polynomial_and_axis_angle_branches_retained": True,
        "correction_cayley_singularity_at_pi_used": False,
        "only_resulting_error_antipode_is_gate": True,
        "signed_quaternion_product_term_retained": True,
        "source_quaternion_normalization_cancels_in_resulting_cayley": True,
        "complete_word_promoted_here": False,
        "audit_cells": {"small": _serial(small), "crossing_pi": _serial(crossing_pi)},
        "P5_DEPLOYED_QUATERNION_CAYLEY_CELL_PRIMITIVE": "PASS" if not failures else "FAIL",
        "next_obligation": "replace the retired correction-Cayley compose in the full 18x18 prefix backend and re-run the first source-complete correction cell",
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "shipping_polynomial_and_axis_angle_branches_retained",
        "only_resulting_error_antipode_is_gate",
        "signed_quaternion_product_term_retained",
        "source_quaternion_normalization_cancels_in_resulting_cayley",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "correction_cayley_singularity_at_pi_used",
        "complete_word_promoted_here",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    large = d.get("audit_cells", {}).get("crossing_pi", {})
    if not float(large.get("correction_norm_upper", 0.0)) > 3.0:
        failures.append("large audit no longer exceeds 3 rad")
    W = large.get("product_scalar", [0.0, 0.0])
    if float(W[0]) <= 0.0 <= float(W[1]):
        failures.append("large audit product scalar crosses zero")
    if not failures and d.get("P5_DEPLOYED_QUATERNION_CAYLEY_CELL_PRIMITIVE") != "PASS":
        failures.append("deployed quaternion/Cayley primitive did not pass")
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
        "status": out["P5_DEPLOYED_QUATERNION_CAYLEY_CELL_PRIMITIVE"],
        "crossing_pi": out["audit_cells"]["crossing_pi"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
