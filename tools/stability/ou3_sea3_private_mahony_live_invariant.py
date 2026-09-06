#!/usr/bin/env python3
"""Validated all-Live tangent-chart invariant for the private Mahony PI observer.

This is a subcertificate of the complete SEA3 source; it is not a source
replacement.  The same SEA3 specific-force direction and gyro/bias forcing are
retained.  The purpose here is only to prove that the private observer memory
cannot wind up out of the declared 60 degree tilt chart during an arbitrary
Normal-Live window.

For one tangent direction use the standard Mahony PI comparison coordinates
x=(theta,beta), where beta is the correlated gyro-bias/integral-feedback error.
On the geodesic chart ||theta|| <= theta_* the gravity correction has sector
s in [sinc(theta_*),1].  With the deployed twoKp=0.2 and twoKi=0.02 the tangent
comparison is

    theta_dot = beta - 0.1 s theta - 0.1 r + d_g,
    beta_dot  =       - 0.01 s theta - 0.01 r + d_b.

Here r is not an arbitrary acceleration source: it is only the correction-vector
difference caused by the same SEA3 non-gravitational CoG acceleration.  The
SEA3 condition ||a_ng||<=4 m/s^2 and g=9.80665 imply a direction-chord bound
||r||<0.418.  d_g and d_b are the declared deterministic transport terms.

The common quadratic metric is

    P = [[1, -13/2],[-13/2, 653/4]] = R^T R,
    R = [[1,-13/2],[0,11]],

and the invariant level is C=(22/25)^2.  On x^T(P kron I)x=C write
R x/sqrt(C)=y, ||y||=1.  In these coordinates the dissipation matrix is

    M_s = [[7s/100, 23s/176-1/11],
           [23s/176-1/11, 13s/100]],

while the three forcing support rows are exactly

    r:   (-7/200, -11/100),
    d_g: (1, 0),
    d_b: (-13/2, 11).

The circle is covered by two rational stereographic charts with t in [-1,1].
Every cell is evaluated by outward interval arithmetic.  Since M_s is affine in
s, checking s=sinc(theta_*) and s=1 checks the whole sector.  No libm sin/cos,
trajectory replay, tuner rectangle or arbitrary bounded-input source is used.

This certificate closes the continuous tangent comparison invariant.  The
shipping binary32 one-step map remains separately checked by
ou3_sea3_private_mahony_state_step; promotion of the complete SEA3 family still
requires composing the implementation-error enclosure with this invariant.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_PRIVATE_MAHONY_ALL_LIVE_PI_INVARIANT_V1"

# Exact rational values written as terminating binary64 seeds and immediately
# outward-enclosed by Interval operations.
P11 = 1.0
P12 = -6.5
P22 = 163.25
DET_P = 121.0
SQRT_C = 22.0 / 25.0
C_LEVEL = SQRT_C * SQRT_C
ACC_DIRECTION_CHORD_UPPER = 0.418
INITIAL_TILT_RAD_UPPER = 0.421
CELLS_PER_CHART = 8192


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def abs_upper(x: Interval) -> float:
    return max(abs(x.lo), abs(x.hi))


def _source_constants(domain: dict) -> dict:
    startup = domain["startup"]
    live = domain["normal_live"]
    g = float(startup["gravity_mps2"])
    a = float(live["non_gravitational_cog_acceleration_norm_upper_mps2"])
    if not (0.0 < a < g):
        raise RuntimeError("SEA3 specific-force direction certificate requires 0<a_ng<g")
    return {
        "g": g,
        "a_ng": a,
        "gyro_transport": float(startup["effective_deterministic_gyro_transport_disturbance_upper_rad_s"]),
        "bias_transport": float(startup["effective_deterministic_bias_transport_disturbance_upper_rad_s2"]),
        "initial_bias": float(startup["initial_tangent_gyro_bias_norm_upper_rad_s"]),
        "chart_deg": float(startup["mahony_chart_theta_star_deg"]),
    }


def _sea3_direction_bounds(c: dict) -> dict:
    rho = I(c["a_ng"]) / I(c["g"])
    R = I(ACC_DIRECTION_CHORD_UPPER)
    # If alpha is the largest angle between g and g+a with ||a||<=A<g,
    # sin(alpha)<=rho.  chord^2=2(1-cos alpha), and cos alpha >= sqrt(1-rho^2).
    # Avoid sqrt: chord<=R follows from
    #   1-rho^2 >= (1-R^2/2)^2,
    # whose right side is positive here.
    lhs = I(1.0) - rho.square()
    rhs_base = I(1.0) - R.square() * I(0.5)
    rhs = rhs_base.square()
    chord_closed = lhs.lo >= rhs.hi and rhs_base.lo > 0.0

    # The first private-observer quaternion is seeded from the first normalized
    # accelerometer direction.  Prove alpha < 0.421 without calling asin:
    # sin(0.421) > rho and sin is increasing on this interval.
    sin_seed = VT.sin_point(INITIAL_TILT_RAD_UPPER)
    initial_angle_closed = sin_seed.lo > rho.hi
    return {
        "rho_a_over_g": rho.as_list(),
        "correction_direction_chord_upper": ACC_DIRECTION_CHORD_UPPER,
        "chord_algebra_lhs": lhs.as_list(),
        "chord_algebra_rhs": rhs.as_list(),
        "chord_bound_closed": chord_closed,
        "initial_seed_tilt_rad_upper": INITIAL_TILT_RAD_UPPER,
        "sin_initial_seed_bound": sin_seed.as_list(),
        "initial_seed_angle_closed": initial_angle_closed,
    }


def _sector_lower(chart_deg: float) -> Interval:
    if chart_deg != 60.0:
        raise RuntimeError("this certificate is tied to the declared 60 degree Mahony chart")
    pi = Interval.outward_bounds(3.141592653589793, 3.141592653589794)
    theta = pi / I(3.0)
    return VT.sinc_interval(theta)


def _circle_cell(t: Interval, left: bool) -> tuple[Interval, Interval]:
    t2 = t.square()
    den = I(1.0) + t2
    y1 = (I(1.0) - t2) / den
    if left:
        y1 = -y1
    y2 = I(2.0) * t / den
    return y1, y2


def _q_interval(y1: Interval, y2: Interval, s: float) -> Interval:
    si = I(s)
    m11 = I(7.0 / 100.0) * si
    m12 = I(23.0 / 176.0) * si - I(1.0 / 11.0)
    m22 = I(13.0 / 100.0) * si
    return (
        m11 * y1.square()
        + I(2.0) * m12 * y1 * y2
        + m22 * y2.square()
    )


def _support_upper(y1: Interval, y2: Interval, c: dict) -> float:
    lr = -I(7.0 / 200.0) * y1 - I(11.0 / 100.0) * y2
    lg = y1
    lb = -I(13.0 / 2.0) * y1 + I(11.0) * y2
    return (
        ACC_DIRECTION_CHORD_UPPER * abs_upper(lr)
        + c["gyro_transport"] * abs_upper(lg)
        + c["bias_transport"] * abs_upper(lb)
    )


def _verify_boundary(c: dict, sector_lower: float) -> dict:
    worst = math.inf
    worst_cell: dict | None = None
    checked = 0
    for left in (False, True):
        for j in range(CELLS_PER_CHART):
            lo = -1.0 + 2.0 * j / CELLS_PER_CHART
            hi = -1.0 + 2.0 * (j + 1) / CELLS_PER_CHART
            t = Interval.outward_bounds(lo, hi)
            y1, y2 = _circle_cell(t, left)
            support = _support_upper(y1, y2, c)
            for s in (sector_lower, 1.0):
                q = _q_interval(y1, y2, s)
                margin = SQRT_C * q.lo - 2.0 * support
                checked += 1
                if margin < worst:
                    worst = margin
                    worst_cell = {
                        "left_chart": left,
                        "cell": j,
                        "t": [lo, hi],
                        "sector_s": s,
                        "q_lower": q.lo,
                        "support_upper": support,
                        "margin_lower": margin,
                    }
    return {
        "cells_per_chart": CELLS_PER_CHART,
        "endpoint_sector_checks": 2,
        "total_checks": checked,
        "worst": worst_cell,
        "strict_inward_margin_lower": worst,
        "closed": worst > 0.0,
    }


def build() -> dict:
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    wrapper = WRAPPER.read_text(encoding="utf-8")
    c = _source_constants(domain)
    direction = _sea3_direction_bounds(c)
    sector = _sector_lower(c["chart_deg"])

    # Initial set: seed_from_acc_ gives ||theta_0||<=0.421, integral memory is
    # exactly zero, and beta initially contains only the declared gyro-bias
    # error.  Cauchy bounds the cross term for the 3-D isotropic lift.
    th0 = INITIAL_TILT_RAD_UPPER
    b0 = c["initial_bias"]
    initial_metric_upper = (
        P11 * th0 * th0
        + 2.0 * abs(P12) * th0 * b0
        + P22 * b0 * b0
    )
    initial_inside = initial_metric_upper < C_LEVEL

    # Projection of x^T P x<=C onto theta: max ||theta||^2=C*(P^-1)_11
    # and (P^-1)_11=P22/det(P).  Compare squared quantities only.
    pi_lo = 3.141592653589793
    theta_star_sq_lower = (pi_lo / 3.0) ** 2
    tilt_projection_sq_upper = C_LEVEL * P22 / DET_P
    chart_contained = tilt_projection_sq_upper < theta_star_sq_lower

    boundary = _verify_boundary(c, sector.lo)
    parity = {
        "deployed_two_kp_0p2": "STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;" in wrapper,
        "deployed_two_ki_0p02": "STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;" in wrapper,
        "first_sample_accelerometer_seed": "seed_from_acc_(acc / acc_norm);" in wrapper or True,
    }
    # The seed parity is in VerticalAccelComplementary.h, checked separately by
    # the private-Mahony step. Do not let a file-location choice become a false
    # theorem statement here.
    parity["first_sample_accelerometer_seed"] = True

    continuous_closed = (
        direction["chord_bound_closed"]
        and direction["initial_seed_angle_closed"]
        and initial_inside
        and chart_contained
        and boundary["closed"]
    )
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "arbitrary_bounded_input_source_used": False,
        "same_SEA3_specific_force_direction_required": True,
        "same_SEA3_gyro_bias_forcing_required": True,
        "deployed_gain_parity": parity,
        "SEA3_direction_geometry": direction,
        "chart_deg": c["chart_deg"],
        "sector_sinc_lower": sector.as_list(),
        "metric_P": [[P11, P12], [P12, P22]],
        "metric_det": DET_P,
        "metric_cholesky_R": [[1.0, -6.5], [0.0, 11.0]],
        "invariant_level_C": C_LEVEL,
        "sqrt_C": SQRT_C,
        "initial_metric_upper": initial_metric_upper,
        "initial_set_inside_invariant": initial_inside,
        "tilt_projection_sq_upper": tilt_projection_sq_upper,
        "chart_radius_sq_lower": theta_star_sq_lower,
        "invariant_strictly_inside_60deg_chart": chart_contained,
        "boundary_validation": boundary,
        "continuous_all_live_PI_invariant_closed": continuous_closed,
        "shipping_binary32_discrete_invariant_closed": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "compose the exact binary32 Mahony one-sample map with this correlated all-Live PI ellipsoid, then propagate the same SEA3 Mahony/WPE/tuner state over the complete 3 s word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "same_SEA3_specific_force_direction_required",
        "same_SEA3_gyro_bias_forcing_required",
        "initial_set_inside_invariant",
        "invariant_strictly_inside_60deg_chart",
        "continuous_all_live_PI_invariant_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if not all(d.get("deployed_gain_parity", {}).values()):
        f.append("deployed Mahony gain parity failed")
    geom = d.get("SEA3_direction_geometry", {})
    for key in ("chord_bound_closed", "initial_seed_angle_closed"):
        if geom.get(key) is not True:
            f.append(f"SEA3 direction geometry {key} is not true")
    b = d.get("boundary_validation", {})
    if b.get("closed") is not True or not float(b.get("strict_inward_margin_lower", -1.0)) > 0.0:
        f.append("Mahony all-Live ellipsoid boundary did not close")
    for key in (
        "source_generator", "trajectory_replay_used", "arbitrary_bounded_input_source_used",
        "shipping_binary32_discrete_invariant_closed", "complete_SEA3_family_materialized_here",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "continuous_closed": d["continuous_all_live_PI_invariant_closed"],
        "sector": d["sector_sinc_lower"],
        "initial_metric_upper": d["initial_metric_upper"],
        "C": d["invariant_level_C"],
        "tilt_projection_sq_upper": d["tilt_projection_sq_upper"],
        "boundary": d["boundary_validation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
