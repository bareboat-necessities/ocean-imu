#!/usr/bin/env python3
"""Analytic physical tilt bound at the first valid startup-proxy sample.

For true rest-specific-force direction -d and non-gravitational acceleration a,
shipping initializes the proxy from

    f = -g d + a,  ||a|| <= A < g.

The largest possible angle between f and -d is asin(A/g).  This is attained
when a has the tangent/axial combination that makes the ray from the origin
tangent to the radius-A ball centered at -g d.  It is a hard geometric bound,
not a replay statistic.

Sensor and finite-precision direction errors are deliberately separate charges;
this module closes the physical BRMM part only and must not be used to claim the
whole 150 s Mahony capture tube.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_STARTUP_PROXY_INITIAL_PHYSICAL_TILT_BOUND_V1"


def build() -> dict:
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    s = domain["startup"]
    g = float(s["gravity_mps2"])
    A = float(s["initial_non_gravitational_specific_force_norm_upper_mps2"])
    if not (0.0 <= A < g):
        raise RuntimeError("initial acceleration bound must satisfy 0<=A<g")
    theta = math.asin(A / g)
    chart = math.radians(float(s["mahony_chart_theta_star_deg"]))
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_fit": False,
        "gravity_mps2": g,
        "non_gravitational_acceleration_norm_upper_mps2": A,
        "physical_initial_tilt_error_upper_rad": theta,
        "physical_initial_tilt_error_upper_deg": math.degrees(theta),
        "proof": "max angle from -g*d to any point in ball B(-g*d,A) is asin(A/g)",
        "mahony_chart_upper_rad": chart,
        "physical_initialization_strictly_inside_chart": theta < chart,
        "physical_chart_margin_rad_before_sensor_fp_charges": chart - theta,
        "sensor_direction_error_charge_closed_here": False,
        "binary32_direction_error_charge_closed_here": False,
        "whole_startup_trajectory_capture_closed_here": False,
        "P5_PASS": False,
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("trajectory_fit") is not False:
        f.append("trajectory fit used")
    if d.get("physical_initialization_strictly_inside_chart") is not True:
        f.append("physical initialization does not fit Mahony chart")
    if not (0.0 < float(d.get("physical_chart_margin_rad_before_sensor_fp_charges", -1.0))):
        f.append("no positive physical chart margin")
    if d.get("whole_startup_trajectory_capture_closed_here") is not False or d.get("P5_PASS") is not False:
        f.append("initialization lemma over-promoted startup capture")
    return f


if __name__ == "__main__":
    d=build(); f=validate(d)
    print(json.dumps({**d,"validation_pass":not f,"validation_failures":f},indent=2))
    raise SystemExit(1 if f else 0)
