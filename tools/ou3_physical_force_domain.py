#!/usr/bin/env python3
"""Validate the physical accelerometer-force domain used by OU-III P1-P5.

The theorem's primary marine-motion assumption is a bound on the CoG
non-gravitational acceleration, not an arbitrary accelerometer magnitude cap.
With lever arm disabled, the predicted accelerometer vector before bias/noise is

    f = R_wb (a_w - g_w),

so rotation preserves its norm and, for ||a_w|| <= a_max,

    g-a_max <= ||f|| <= g+a_max.

For the declared normal-Live marine envelope a_max=4 m/s^2 (0.408 g), yielding
5.80665..13.80665 m/s^2 (0.592..1.408 g).  Impact/slam dynamics and the active
vibration-guard branch are explicitly outside this P1-P5 source branch.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    d = json.loads(path.read_text(encoding="utf-8"))
    live = d["normal_live"]
    startup = d["startup"]
    runtime = d["configured_runtime"]
    g = float(startup["gravity_mps2"])
    amax = float(live["non_gravitational_cog_acceleration_norm_upper_mps2"])
    derived_lo = g - amax
    derived_hi = g + amax
    failures = []
    if not (g > 0.0 and 0.0 <= amax < g):
        failures.append("marine CoG acceleration bound must satisfy 0 <= a_max < g")
    if runtime.get("imu_lever_arm_enabled") is not False:
        failures.append("derived CoG force envelope requires disabled lever arm")
    if live.get("specific_force_bounds_derived_from_gravity_and_non_gravitational_acceleration") is not True:
        failures.append("specific-force bounds are not declared derived")
    if live.get("impact_slam_acceleration_in_normal_live_P1_P5_scope") is not False:
        failures.append("impact/slam dynamics were silently admitted to normal-Live P1-P5")
    if not math.isclose(float(live["specific_force_norm_lower_mps2"]), derived_lo, rel_tol=0.0, abs_tol=1e-12):
        failures.append("specific-force lower bound is not g-a_max")
    if not math.isclose(float(live["specific_force_norm_upper_mps2"]), derived_hi, rel_tol=0.0, abs_tol=1e-12):
        failures.append("specific-force upper bound is not g+a_max")
    if not math.isclose(float(live["specific_force_norm_lower_fraction_g"]), derived_lo/g, rel_tol=0.0, abs_tol=1e-15):
        failures.append("specific-force lower g-fraction drifted")
    if not math.isclose(float(live["specific_force_norm_upper_fraction_g"]), derived_hi/g, rel_tol=0.0, abs_tol=1e-15):
        failures.append("specific-force upper g-fraction drifted")
    if not math.isclose(float(live["non_gravitational_cog_acceleration_fraction_g"]), amax/g, rel_tol=0.0, abs_tol=1e-15):
        failures.append("non-gravitational acceleration g-fraction drifted")
    if not math.isclose(float(startup["initial_non_gravitational_specific_force_norm_upper_mps2"]), amax, rel_tol=0.0, abs_tol=1e-12):
        failures.append("startup and normal-Live non-gravitational acceleration envelopes differ")
    qamax = float(live["gravity_quotient"]["non_gravitational_specific_force_norm_upper_mps2"])
    if not math.isclose(qamax, amax, rel_tol=0.0, abs_tol=1e-12):
        failures.append("gravity-quotient acceleration envelope differs from full-heading normal-Live envelope")
    if float(live["specific_force_norm_upper_mps2"]) >= 2.0*g:
        failures.append("normal-Live total specific-force envelope unexpectedly reaches >=2 g")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_MARINE_COG_FORCE_DOMAIN",
        "source_generated_not_replay_fit": True,
        "primary_physical_assumption": "||a_non-grav,CoG|| <= a_max",
        "gravity_mps2": g,
        "non_gravitational_cog_acceleration_norm_upper_mps2": amax,
        "non_gravitational_cog_acceleration_fraction_g": amax/g,
        "derived_specific_force_norm_lower_mps2": derived_lo,
        "derived_specific_force_norm_upper_mps2": derived_hi,
        "derived_specific_force_norm_lower_fraction_g": derived_lo/g,
        "derived_specific_force_norm_upper_fraction_g": derived_hi/g,
        "old_independent_30_mps2_cap_used": False,
        "impact_slam_in_scope": False,
        "lever_arm_enabled": bool(runtime.get("imu_lever_arm_enabled")),
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "failures": failures,
    }


def validate(x: dict) -> list[str]:
    f = list(x.get("failures", []))
    if x.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if x.get("qualification") != "OU3_MARINE_COG_FORCE_DOMAIN":
        f.append("wrong qualification")
    if x.get("source_generated_not_replay_fit") is not True:
        f.append("source_generated_not_replay_fit is not true")
    for key in ("old_independent_30_mps2_cap_used", "impact_slam_in_scope", "lever_arm_enabled", "P4_USABLE_CERTIFICATE_PROMOTED"):
        if x.get(key) is not False:
            f.append(f"{key} is not false")
    if not (0.0 < float(x.get("derived_specific_force_norm_lower_mps2", 0.0)) < float(x.get("gravity_mps2", 0.0))):
        f.append("derived lower force bound is not below one g and positive")
    if not (float(x.get("gravity_mps2", math.inf)) < float(x.get("derived_specific_force_norm_upper_mps2", 0.0)) < 2.0*float(x.get("gravity_mps2", 0.0))):
        f.append("derived upper force bound is not between one and two g")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(d, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
