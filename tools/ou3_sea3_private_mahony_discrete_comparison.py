#!/usr/bin/env python3
"""Discrete 5 ms closure of the SEA3 private-Mahony PI comparison ellipsoid.

The all-Live invariant producer proves a strict continuous boundary margin for
exactly the same SEA3 correction-direction/gyro/bias forcing.  This module
checks that the deployed 5 ms explicit step cannot consume that margin.

For x^T P x=C and x+=h(A_s x+f),

  V(x+) - V(x)
   = h[-x^T Q_s x + 2 x^T P f]
     + h^2 ||A_s x+f||_P^2.

The first bracket is bounded by the validated rational-circle sweep from
``ou3_sea3_private_mahony_live_invariant``.  The h^2 term is bounded in the
same Cholesky coordinates y=R x/sqrt(C), so no independent theta/beta box is
introduced.  This remains the ideal tangent comparison step; exact shipping
binary32/quaternion-map error is a separate, fail-closed composition step.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_sea3_private_mahony_live_invariant as LIVE

REPO = Path(__file__).resolve().parents[1]
DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_PRIVATE_MAHONY_DISCRETE_PI_COMPARISON_V1"

# Rational-safe uppers for the transformed operator/input norms.  The
# Frobenius interval is evaluated over the whole sector at once, so entry-wise
# interval dependency makes its enclosure slightly wider than the correlated
# endpoint value (whose true maximum is below 0.150).  0.152 covers that
# outward interval hull and still leaves more than an order of magnitude of
# the 5 ms inward metric margin.  These are theorem calculations, not fitted
# trajectory values.
N_FROBENIUS_UPPER = 0.152
RB_R_NORM_UPPER = 0.116
RB_G_NORM_UPPER = 1.0
RB_B_NORM_UPPER = 12.78


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _norm_bounds(sector_lo: float) -> dict:
    s = Interval.outward_bounds(sector_lo, 1.0)
    # N_s = R A_s R^{-1} exactly for R=[[1,-6.5],[0,11]].
    n11 = -I(7.0 / 200.0) * s
    n12 = I(1.0 / 11.0) - I(91.0 / 4400.0) * s
    n21 = -I(11.0 / 100.0) * s
    n22 = -I(13.0 / 200.0) * s
    fro2 = n11.square() + n12.square() + n21.square() + n22.square()
    n_closed = fro2.hi < N_FROBENIUS_UPPER * N_FROBENIUS_UPPER

    # R B_r=(-0.035,-0.11), R B_g=(1,0), R B_b=(-6.5,11).
    rb_r2 = I(7.0 / 200.0).square() + I(11.0 / 100.0).square()
    rb_b2 = I(13.0 / 2.0).square() + I(11.0).square()
    inputs_closed = (
        rb_r2.hi < RB_R_NORM_UPPER * RB_R_NORM_UPPER
        and 1.0 <= RB_G_NORM_UPPER
        and rb_b2.hi < RB_B_NORM_UPPER * RB_B_NORM_UPPER
    )
    return {
        "N_frobenius_sq_upper": fro2.hi,
        "N_frobenius_declared_upper": N_FROBENIUS_UPPER,
        "N_bound_closed": n_closed,
        "RB_r_norm_sq_upper": rb_r2.hi,
        "RB_r_norm_declared_upper": RB_R_NORM_UPPER,
        "RB_g_norm_declared_upper": RB_G_NORM_UPPER,
        "RB_b_norm_sq_upper": rb_b2.hi,
        "RB_b_norm_declared_upper": RB_B_NORM_UPPER,
        "input_norm_bounds_closed": inputs_closed,
    }


def build() -> dict:
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    dt = float(domain["configured_runtime"]["imu_dt_s"])
    if dt != 0.005:
        raise RuntimeError("discrete comparison is tied to shipping 5 ms theorem runtime")

    live = LIVE.build()
    lf = LIVE.validate(live)
    if lf:
        raise RuntimeError(f"continuous all-Live invariant prerequisite failed: {lf}")
    sector_lo = float(live["sector_sinc_lower"][0])
    norm = _norm_bounds(sector_lo)

    geom = live["SEA3_direction_geometry"]
    r = float(geom["correction_direction_chord_upper"])
    startup = domain["startup"]
    dg = float(startup["effective_deterministic_gyro_transport_disturbance_upper_rad_s"])
    db = float(startup["effective_deterministic_bias_transport_disturbance_upper_rad_s2"])

    f_p_norm_upper = (
        r * RB_R_NORM_UPPER
        + dg * RB_G_NORM_UPPER
        + db * RB_B_NORM_UPPER
    )
    # ||A x||_P <= ||R A R^-1||_2 sqrt(C), and ||.||_2<=||.||_F.
    velocity_p_norm_upper = (
        N_FROBENIUS_UPPER * LIVE.SQRT_C + f_p_norm_upper
    )
    h2_term_upper = dt * dt * velocity_p_norm_upper * velocity_p_norm_upper

    continuous_margin = float(
        live["boundary_validation"]["strict_inward_margin_lower"]
    )
    # The live module's normalized boundary inequality is
    # sqrt(C) q - 2 support >= m. Multiplying by sqrt(C) gives the
    # first-order metric derivative margin xQx-2xPf >= sqrt(C)*m.
    first_order_metric_margin = LIVE.SQRT_C * continuous_margin
    first_order_step_decrease = dt * first_order_metric_margin
    discrete_metric_decrease_lower = first_order_step_decrease - h2_term_upper
    closed = (
        norm["N_bound_closed"]
        and norm["input_norm_bounds_closed"]
        and discrete_metric_decrease_lower > 0.0
    )
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "same_SEA3_forcing_as_continuous_invariant": True,
        "sample_period_s": dt,
        "continuous_invariant_consumed": True,
        "norm_bounds": norm,
        "forcing_P_norm_upper": f_p_norm_upper,
        "boundary_velocity_P_norm_upper": velocity_p_norm_upper,
        "h2_metric_increase_upper": h2_term_upper,
        "continuous_normalized_boundary_margin_lower": continuous_margin,
        "first_order_step_decrease_lower": first_order_step_decrease,
        "discrete_metric_decrease_lower": discrete_metric_decrease_lower,
        "ideal_5ms_discrete_PI_invariant_closed": closed,
        "shipping_binary32_quaternion_map_error_composed": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "bound the difference between the actual binary32 Mahony quaternion/integral step and this same tangent PI comparison inside the invariant, then carry that correlated state into WPE/tuner over the 3 s SEA3 word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "same_SEA3_forcing_as_continuous_invariant",
        "continuous_invariant_consumed",
        "ideal_5ms_discrete_PI_invariant_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    n = d.get("norm_bounds", {})
    if n.get("N_bound_closed") is not True or n.get("input_norm_bounds_closed") is not True:
        f.append("discrete comparison norm bound did not close")
    if not float(d.get("discrete_metric_decrease_lower", -1.0)) > 0.0:
        f.append("5 ms h^2 term consumed continuous invariant margin")
    for key in (
        "source_generator",
        "shipping_binary32_quaternion_map_error_composed",
        "complete_SEA3_family_materialized_here",
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
        "norm_bounds": d["norm_bounds"],
        "h2_upper": d["h2_metric_increase_upper"],
        "first_order_step_decrease": d["first_order_step_decrease_lower"],
        "discrete_metric_decrease": d["discrete_metric_decrease_lower"],
        "closed": d["ideal_5ms_discrete_PI_invariant_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
