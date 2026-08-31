#!/usr/bin/env python3
"""Non-microscopic P5 startup capture into the operation-matched outer sector.

P5 used to ask whether the P1 handoff was already inside the tiny perturbative
P4 Lyapunov level.  PR #441 proved that the old uniform-defect route can never
answer yes: its absolute attitude ceiling is milliradians.  The replacement P4
geometry is a finite-angle sector, so the correct first P5 question is whether
startup enters that sector source-faithfully.

For gauged branches the full SO(3) Cayley radii come from the heading handoff
contract.  For an ungauged timeout, full yaw is a gauge and no fictitious
full-heading radius is assigned.  P1 supplies a lower bound on the true gravity
cosine.  Since cosine decreases on [0,pi], inclusion in the 0.80 rad tilt sector
is proved only when that P1 lower bound exceeds the *upper* validated enclosure
of cos(0.80).  Using the lower enclosure would be fail-open at the boundary.

This certificate establishes capture into the outer finite-angle sector with
N_outer=0.  It does not claim capture into the legacy microscopic P4 inner
Lyapunov seed, and it does not claim that complete-word sector dissipation is
already numerically composed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p5_heading_handoff_contract as HEADING
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P5 outer capture must not be trajectory fitted")

    p1 = P1.build(path)
    heading = HEADING.build(path)
    sector = SECTOR.build(path)
    failures = [f"P1: {x}" for x in P1.validate(p1)]
    failures += [f"heading: {x}" for x in HEADING.validate(heading)]
    failures += [f"sector: {x}" for x in SECTOR.validate(sector)]

    theta_sector = float(sector["design_full_attitude_angle_rad"])
    q_sector = float(sector["design_cayley_norm_upper"])
    sector_cos_lower = float(sector["design_full_attitude_cosine_lower"])
    sector_cos_upper = float(sector["design_full_attitude_cosine_upper"])
    normal = sector["P1_overlap"]

    timeout_tilt_cos = float(p1["timeout_handoff"]["combined_true_gravity_cosine_lower"])
    # P1 proves cos(theta_tilt) >= timeout_tilt_cos.  The sector producer proves
    # cos(theta_sector) <= sector_cos_upper.  Therefore timeout_tilt_cos >=
    # sector_cos_upper implies cos(theta_tilt) >= cos(theta_sector), hence
    # theta_tilt <= theta_sector on [0,pi].
    ungauged_tilt_inside = timeout_tilt_cos >= sector_cos_upper

    physical = dict(p1["go_live"]["physical_coordinate_bounds"])
    for key, value in physical.items():
        if not (
            isinstance(value, (int, float))
            and math.isfinite(float(value))
            and float(value) >= 0.0
        ):
            failures.append(f"invalid P1 physical handoff bound {key}")

    if normal.get("normal_gauged_inside_sector") is not True:
        failures.append("normal gauged handoff misses outer sector")
    if normal.get("timeout_gauged_inside_sector") is not True:
        failures.append("gauged timeout handoff misses outer sector")
    if not ungauged_tilt_inside:
        failures.append("ungauged timeout tilt misses outer gravity sector")

    branches = {
        "normal_gauged": {
            "attitude_representation": "FULL_SO3_CAYLEY",
            "cayley_norm_upper": float(normal["normal_gauged_cayley_norm_upper"]),
            "inside_outer_sector": bool(normal["normal_gauged_inside_sector"]),
        },
        "timeout_gauged": {
            "attitude_representation": "FULL_SO3_CAYLEY",
            "cayley_norm_upper": float(normal["timeout_gauged_cayley_norm_upper"]),
            "inside_outer_sector": bool(normal["timeout_gauged_inside_sector"]),
        },
        "timeout_ungauged": {
            "attitude_representation": "GRAVITY_DIRECTION_QUOTIENT",
            "full_heading_radius_assigned": False,
            "tilt_cosine_lower": timeout_tilt_cos,
            "outer_sector_cosine_lower": sector_cos_lower,
            "outer_sector_cosine_upper": sector_cos_upper,
            "boundary_cosine_direction_used": "UPPER_ENCLOSURE",
            "inside_outer_gravity_sector": ungauged_tilt_inside,
            "yaw_role": "GAUGE_UNTIL_MAGNETIC_REGAUGE",
            "gravity_parallel_gyro_bias_role": "BOUNDED_NEUTRAL_INPUT_UNTIL_EXCITATION_OR_MAGNETIC_REGAUGE",
            "required_route": heading["ungauged_timeout_subbranch"]["required_route"],
        },
    }

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_STARTUP_TO_OPERATION_MATCHED_OUTER_SECTOR_CAPTURE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "outer_sector_angle_rad": theta_sector,
        "outer_sector_angle_deg": math.degrees(theta_sector),
        "outer_sector_cosine_lower": sector_cos_lower,
        "outer_sector_cosine_upper": sector_cos_upper,
        "outer_sector_cayley_norm_upper": q_sector,
        "validated_sector_boundary_consumed": True,
        "upper_cosine_enclosure_used_for_ungauged_inclusion": True,
        "physical_handoff_product_box": physical,
        "branches": branches,
        "all_source_handoff_branches_enter_outer_sector": passed,
        "N_outer_words": 0 if passed else None,
        "legacy_microscopic_inner_seed_used_as_outer_capture_target": False,
        "legacy_uniform_transport_route_used": False,
        "P5_OUTER_SECTOR_CAPTURE_CERTIFICATE": "PASS" if passed else "FAIL",
        "P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "finish operation-matched complete-word dissipation on the outer sector; once it overlaps the existing inner P4 seed, "
            "derive a finite word count to the inner stochastic localization level without using the retired uniform-defect ceiling"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("P5 outer capture is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("P5 outer capture uses replay")
    if d.get("filter_changed") is not False:
        failures.append("P5 outer capture changes the filter")
    if d.get("validated_sector_boundary_consumed") is not True:
        failures.append("P5 did not consume the validated sector boundary")
    if d.get("upper_cosine_enclosure_used_for_ungauged_inclusion") is not True:
        failures.append("ungauged timeout did not use upper sector cosine enclosure")
    clo = d.get("outer_sector_cosine_lower")
    chi = d.get("outer_sector_cosine_upper")
    if not (
        isinstance(clo, (int, float))
        and isinstance(chi, (int, float))
        and math.isfinite(float(clo))
        and math.isfinite(float(chi))
        and float(clo) <= float(chi)
    ):
        failures.append("invalid outer-sector cosine enclosure")
    if d.get("P5_OUTER_SECTOR_CAPTURE_CERTIFICATE") != "PASS":
        failures.append("P5 outer-sector capture did not pass")
    if d.get("all_source_handoff_branches_enter_outer_sector") is not True:
        failures.append("not all P1 handoff branches enter outer sector")
    if d.get("N_outer_words") != 0:
        failures.append("outer sector is not immediate at handoff")
    if not float(d.get("outer_sector_angle_rad", 0.0)) >= 0.80:
        failures.append("P5 outer sector regressed below 0.80 rad")
    if d.get("legacy_microscopic_inner_seed_used_as_outer_capture_target") is not False:
        failures.append("P5 returned to microscopic inner-seed capture")
    if d.get("legacy_uniform_transport_route_used") is not False:
        failures.append("P5 returned to uniform transport ceiling route")
    if d.get("P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE") is not False:
        failures.append("outer-sector certificate prematurely claims inner capture")

    ungauged = d.get("branches", {}).get("timeout_ungauged", {})
    if ungauged.get("full_heading_radius_assigned") is not False:
        failures.append("ungauged timeout assigned fictitious heading radius")
    if ungauged.get("boundary_cosine_direction_used") != "UPPER_ENCLOSURE":
        failures.append("ungauged timeout used wrong cosine enclosure direction")
    if ungauged.get("inside_outer_gravity_sector") is not True:
        failures.append("ungauged timeout tilt is outside outer sector")
    try:
        tilt_lower = float(ungauged["tilt_cosine_lower"])
        sector_upper = float(ungauged["outer_sector_cosine_upper"])
        if not tilt_lower >= sector_upper:
            failures.append("ungauged timeout cosine inequality does not prove inclusion")
    except (KeyError, TypeError, ValueError):
        failures.append("ungauged timeout cosine proof fields missing")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_OUTER_SECTOR_CAPTURE_CERTIFICATE"],
        "sector_rad": out["outer_sector_angle_rad"],
        "sector_cosine_lower": out["outer_sector_cosine_lower"],
        "sector_cosine_upper": out["outer_sector_cosine_upper"],
        "N_outer_words": out["N_outer_words"],
        "branches": out["branches"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
