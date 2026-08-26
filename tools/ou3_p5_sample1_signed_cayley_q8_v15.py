#!/usr/bin/env python3
"""V15: hybrid SO(3)-geodesic + signed-product closure for sample-1 q<8.

The fine V14D run closes the signed/radial prerequisite but leaves many q<8
cells open because it bounds the correction/current quaternion scalar product
componentwise.  In particular, the first witness has current Cayley norm
q<=0.6594 and correction radial upper <=1.7249 rad.  Those two rotations cannot
approach the Cayley antipode: their SO(3) geodesic angles add to <pi.

For an axis-angle correction cell (the deployed branch above the 1e-2 series
threshold), let phi_d be the maximum *principal* SO(3) rotation angle represented
by its radial interval.  For current Cayley radius q,

    phi_c <= 2 atan(q/2).

Bi-invariance of the SO(3) geodesic metric gives

    phi_plus <= phi_c + phi_d.

Whenever that sum is <pi, this independently proves both a positive product
scalar separation and

    q_plus <= 2 tan((phi_c + phi_d)/2).

V15 installs that bound only as an additional closure route inside the unchanged
V14D build.  The original signed v_d,x*c_x + yz-Cauchy product bound remains in
force and is used whenever it is stronger or the geodesic sum is not useful.
The normalized polynomial correction branch below 1e-2 rad is deliberately not
changed by this metric shortcut.

No estimator, source domain, deployed correction limit, q<8 target, source
branch, or theorem promotion state is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D

DEFAULT_DOMAIN = V14D.DEFAULT_DOMAIN
SCHEMA = 1500
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET
SERIES = V14.SERIES
TWO_PI = 2.0 * math.pi


def _principal_axis_angle_upper(radial_lower: float, radial_upper: float) -> float | None:
    """Upper principal SO(3) angle for one axis-angle radial interval <=9 rad."""
    lo = float(radial_lower)
    hi = float(radial_upper)
    if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 <= lo <= hi <= 9.0):
        raise ValueError("invalid correction radial interval")
    # The tiny source-polynomial branch is not literally the axis-angle map.
    if lo < SERIES:
        return None
    pi = math.pi
    if hi <= pi:
        ans = hi
    elif lo <= pi <= hi:
        ans = pi
    elif hi <= TWO_PI:
        # On [pi,2pi], the principal angle decreases from pi to zero.
        ans = TWO_PI - lo
    elif lo < TWO_PI < hi:
        # V-shaped around 2pi; 9 rad is still below 3pi.
        ans = max(TWO_PI - lo, hi - TWO_PI)
    else:
        # On [2pi,9] (<3pi), the principal angle increases from zero.
        ans = hi - TWO_PI
    return FULL.up(max(0.0, min(pi, ans)))


def _geodesic_q_and_scalar_lower(q_current_upper: float,
                                 radial_lower: float,
                                 radial_upper: float) -> tuple[float, float] | None:
    """Return rigorous (q_plus upper, |W| lower) from SO(3) triangle geometry."""
    q = float(q_current_upper)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("invalid current Cayley radius")
    phi_d = _principal_axis_angle_upper(radial_lower, radial_upper)
    if phi_d is None:
        return None
    phi_c = FULL.up(2.0 * math.atan(FULL.up(0.5 * q)))
    phi = FULL.up(phi_c + phi_d)
    if not phi < math.pi:
        return None
    half = FULL.up(0.5 * phi)
    q_plus = FULL.up(2.0 * math.tan(half))
    # V14's unnormalized current quaternion is [2,c], so its norm is >=2.
    # For a product principal angle <=phi, |W|/||[2,c]|| >= cos(phi/2).
    # Thus the source-independent lower bound below is rigorous.
    w_lower = FULL.down(2.0 * math.cos(half))
    if not w_lower > 0.0:
        return None
    return q_plus, w_lower


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    context = {
        "radial": None,
        "calls": 0,
        "geodesic_available": 0,
        "geodesic_improved": 0,
        "geodesic_newly_closed": 0,
        "first_geodesic_newly_closed": None,
    }
    original_quat = V14D.radial_sinc_normalized_shipping_quaternion
    original_qplus = V14._qplus_from_product_scalar

    def tracked_quat(dbox, *, radial_lower: float, radial_upper: float):
        context["radial"] = (float(radial_lower), float(radial_upper))
        return original_quat(
            dbox, radial_lower=radial_lower, radial_upper=radial_upper)

    def hybrid_qplus(q1: float, W):
        context["calls"] += 1
        w_signed, q_signed = original_qplus(q1, W)
        radial = context.get("radial")
        if radial is None:
            return w_signed, q_signed
        geo = _geodesic_q_and_scalar_lower(q1, radial[0], radial[1])
        if geo is None:
            return w_signed, q_signed
        context["geodesic_available"] += 1
        q_geo, w_geo = geo
        q_best = min(q_signed, q_geo)
        w_best = max(w_signed, w_geo)
        if q_geo < q_signed:
            context["geodesic_improved"] += 1
        signed_closed = math.isfinite(q_signed) and q_signed < Q_TARGET and w_signed > 0.0
        hybrid_closed = math.isfinite(q_best) and q_best < Q_TARGET and w_best > 0.0
        if hybrid_closed and not signed_closed:
            context["geodesic_newly_closed"] += 1
            if context["first_geodesic_newly_closed"] is None:
                context["first_geodesic_newly_closed"] = {
                    "q_current_upper": float(q1),
                    "correction_radial_lower_rad": radial[0],
                    "correction_radial_upper_rad": radial[1],
                    "signed_product_q_upper": q_signed,
                    "geodesic_q_upper": q_geo,
                    "signed_product_abs_W_lower": w_signed,
                    "geodesic_abs_W_lower": w_geo,
                }
        return w_best, q_best

    V14D.radial_sinc_normalized_shipping_quaternion = tracked_quat
    V14._qplus_from_product_scalar = hybrid_qplus
    try:
        core = V14D.build(
            Path(domain_path).resolve(),
            source_pieces=source_pieces,
            source_cell_index=source_cell_index,
            p_pieces=p_pieces,
            tangent_pieces=tangent_pieces,
            axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
        )
    finally:
        V14._qplus_from_product_scalar = original_qplus
        V14D.radial_sinc_normalized_shipping_quaternion = original_quat

    inherited = V14D.validate(core)
    status = core.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_GEODESIC_SIGNED_CAYLEY_Q8_V15",
        "V14D_radial_sinc_signed_product_parent_retained": True,
        "SO3_biinvariant_geodesic_triangle_bound_added": True,
        "current_cayley_to_principal_angle_exact_map_used": True,
        "axis_angle_radial_to_principal_SO3_angle_used": True,
        "series_branch_geodesic_shortcut_used": False,
        "signed_product_bound_retained_when_stronger": True,
        "geodesic_product_scalar_lower_is_independent_valid_bound": True,
        "geodesic_qplus_calls": int(context["calls"]),
        "geodesic_bound_available_cells": int(context["geodesic_available"]),
        "geodesic_bound_improved_cells": int(context["geodesic_improved"]),
        "geodesic_bound_newly_closed_cells": int(context["geodesic_newly_closed"]),
        "first_geodesic_newly_closed_cell": context["first_geodesic_newly_closed"],
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_GEODESIC_SIGNED_CAYLEY_Q8_V15": (
            "PASS" if status == "PASS" and not inherited else "NOT_ESTABLISHED"
        ),
        "next_obligation": (
            "LIFT_CLOSED_SAMPLE1_PREFIX_TO_ALL_SOURCE_PHASE_CELLS_AND_CONTINUE_SAMPLE2_PREFIX"
            if status == "PASS" and not inherited else
            "REFINE_REMAINING_HIGH_PRINCIPAL_ANGLE_CELLS_WITH_SOURCE_CORRELATED_CURRENT_DIRECTION"
        ),
        "failures": list(dict.fromkeys(inherited)),
    })
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_GEODESIC_SIGNED_CAYLEY_Q8_V15":
        failures.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V14D_radial_sinc_signed_product_parent_retained",
        "SO3_biinvariant_geodesic_triangle_bound_added",
        "current_cayley_to_principal_angle_exact_map_used",
        "axis_angle_radial_to_principal_SO3_angle_used",
        "signed_product_bound_retained_when_stronger",
        "geodesic_product_scalar_lower_is_independent_valid_bound",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "series_branch_geodesic_shortcut_used",
        "deployed_correction_limit_increased", "q8_word_promoted_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction limit changed")
    calls = int(d.get("geodesic_qplus_calls", -1))
    available = int(d.get("geodesic_bound_available_cells", -1))
    improved = int(d.get("geodesic_bound_improved_cells", -1))
    newly = int(d.get("geodesic_bound_newly_closed_cells", -1))
    if not (calls >= available >= improved >= newly >= 0):
        failures.append("invalid geodesic accounting")
    st = d.get("P5_SAMPLE1_GEODESIC_SIGNED_CAYLEY_Q8_V15")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            failures.append("V15 PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < Q_TARGET:
            failures.append("V15 PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is True:
            failures.append("V15 nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not failures:
            failures.append("V15 numerical nonclosure lacks q8 witness")
    else:
        failures.append("invalid V15 status")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=6)
    ap.add_argument("--parallel-pieces", type=int, default=6)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain,
        source_pieces=args.source_pieces,
        source_cell_index=args.source_cell_index,
        p_pieces=args.p_pieces,
        tangent_pieces=args.tangent_pieces,
        axial_pieces=args.axial_pieces,
        residual_x_pieces=args.residual_x_pieces,
        parallel_pieces=args.parallel_pieces,
    )
    vf = validate(out)
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_GEODESIC_SIGNED_CAYLEY_Q8_V15"],
        "cells": out.get("evaluated_signed_cayley_cells"),
        "unclosed": out.get("unclosed_q8_cells"),
        "antipode_interval_cells": out.get("product_scalar_antipode_cells"),
        "geodesic_available": out["geodesic_bound_available_cells"],
        "geodesic_improved": out["geodesic_bound_improved_cells"],
        "geodesic_newly_closed": out["geodesic_bound_newly_closed_cells"],
        "first_geodesic_newly_closed": out["first_geodesic_newly_closed_cell"],
        "first_unclosed": out.get("first_unclosed_q8_cell"),
        "worst": out.get("worst_q8_cell"),
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
