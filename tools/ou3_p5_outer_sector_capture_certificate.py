#!/usr/bin/env python3
"""P5 entrance capture into the broad operation-matched outer geometry sector.

The P5 entrance is now declared independently of P4: a 45 degree SO(3)
geodesic attitude-error bound on gauged branches (tilt only on an ungauged
branch) and |delta p_i| <= 0.5 Hs.  This is intentionally not an independent
Euler-angle box.

The broad P4 operation-matched geometry sector remains 0.80 rad.  Therefore the
45 degree P5 attitude entrance is immediately inside that geometry envelope.
This object does not claim that the expensive complete 18/21-state P4 word must
be certified on the whole 0.80 rad sector.  A separate search ladder may use
narrower candidate P4 contraction sectors, and P5 must then prove finite capture
from the 45 degree entrance to the widest candidate that actually closes.

The P1 conservative physical handoff box is retained separately.  In particular,
this producer does not silently replace P1's 20 m position handoff bound with
the new 0.5 Hs entrance assumption or assume that an earlier startup interval
preserves the entrance set.  That propagation remains an explicit proof
obligation when the theorem starts before P5.

For an ungauged timeout, full yaw remains a gauge and no fictitious full-heading
radius is assigned.  P1 supplies a lower bound on the true gravity cosine.  The
legacy source-handoff inclusion test still compares it against the *upper*
validated enclosure of cos(0.80), which is the conservative direction.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_heading_handoff_contract as HEADING
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P5 outer capture must not be trajectory fitted")

    p1 = P1.build(path)
    heading = HEADING.build(path)
    sector = SECTOR.build(path)
    entrance = ENTRANCE.build(path)
    failures = [f"P1: {x}" for x in P1.validate(p1)]
    failures += [f"heading: {x}" for x in HEADING.validate(heading)]
    failures += [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]

    theta_sector = float(sector["design_full_attitude_angle_rad"])
    q_sector = float(sector["design_cayley_norm_upper"])
    sector_cos_lower = float(sector["design_full_attitude_cosine_lower"])
    sector_cos_upper = float(sector["design_full_attitude_cosine_upper"])
    p1_overlap = sector["P1_overlap"]

    p5e = entrance["P5_entrance"]
    q_entrance = float(p5e["attitude_geometry"]["cayley_norm_upper"])
    entrance_inside_outer = q_entrance <= q_sector
    if not entrance_inside_outer:
        failures.append("declared 45 degree P5 entrance misses broad outer geometry sector")

    normal_q = float(p1_overlap["normal_gauged_cayley_norm_upper"])
    timeout_q = float(p1_overlap["timeout_gauged_cayley_norm_upper"])
    normal_inside_entrance = normal_q <= q_entrance
    timeout_inside_entrance = timeout_q <= q_entrance
    if not normal_inside_entrance:
        failures.append("normal gauged P1 handoff exceeds declared 45 degree P5 entrance")
    if not timeout_inside_entrance:
        failures.append("gauged timeout P1 handoff exceeds declared 45 degree P5 entrance")

    timeout_tilt_cos = float(p1["timeout_handoff"]["combined_true_gravity_cosine_lower"])
    # P1 proves cos(theta_tilt) >= timeout_tilt_cos.  The sector producer proves
    # cos(theta_sector) <= sector_cos_upper.  Therefore timeout_tilt_cos >=
    # sector_cos_upper implies theta_tilt <= theta_sector on [0,pi].
    ungauged_tilt_inside = timeout_tilt_cos >= sector_cos_upper
    if not ungauged_tilt_inside:
        failures.append("ungauged timeout tilt misses outer gravity sector")

    physical = dict(p1["go_live"]["physical_coordinate_bounds"])
    for key, value in physical.items():
        if not (
            isinstance(value, (int, float))
            and math.isfinite(float(value))
            and float(value) >= 0.0
        ):
            failures.append(f"invalid P1 physical handoff bound {key}")

    if p1_overlap.get("normal_gauged_inside_sector") is not True:
        failures.append("normal gauged handoff misses broad outer sector")
    if p1_overlap.get("timeout_gauged_inside_sector") is not True:
        failures.append("gauged timeout handoff misses broad outer sector")

    branches = {
        "normal_gauged": {
            "attitude_representation": "FULL_SO3_CAYLEY",
            "P1_handoff_cayley_norm_upper": normal_q,
            "declared_P5_entrance_cayley_norm_upper": q_entrance,
            "inside_declared_P5_entrance": normal_inside_entrance,
            "inside_outer_sector": bool(p1_overlap["normal_gauged_inside_sector"]),
        },
        "timeout_gauged": {
            "attitude_representation": "FULL_SO3_CAYLEY",
            "P1_handoff_cayley_norm_upper": timeout_q,
            "declared_P5_entrance_cayley_norm_upper": q_entrance,
            "inside_declared_P5_entrance": timeout_inside_entrance,
            "inside_outer_sector": bool(p1_overlap["timeout_gauged_inside_sector"]),
        },
        "timeout_ungauged": {
            "attitude_representation": "GRAVITY_DIRECTION_QUOTIENT",
            "full_heading_radius_assigned": False,
            "declared_entrance_tilt_upper_deg": float(p5e["gauged_full_attitude_angle_upper_deg"]),
            "declared_entrance_tilt_inside_outer_sector": entrance_inside_outer,
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
        "qualification": "OU3_P5_DECLARED_ENTRANCE_TO_OPERATION_MATCHED_OUTER_GEOMETRY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_P5_entrance": p5e,
        "P5_entrance_is_distinct_from_P4_complete_word_sector": True,
        "P5_entrance_set_inside_outer_geometry_sector": entrance_inside_outer,
        "P1_source_gauged_handoffs_inside_declared_P5_entrance": normal_inside_entrance and timeout_inside_entrance,
        "P1_conservative_physical_handoff_product_box": physical,
        "physical_handoff_product_box": physical,
        "P1_conservative_handoff_box_replaced": False,
        "startup_propagation_of_entrance_assumed_without_proof": False,
        "sea_scaled_translation_coordinates": entrance["sea_scaled_translation_coordinates"],
        "P4_complete_word_search": entrance["P4_complete_word_search"],
        "outer_sector_angle_rad": theta_sector,
        "outer_sector_angle_deg": math.degrees(theta_sector),
        "outer_sector_cosine_lower": sector_cos_lower,
        "outer_sector_cosine_upper": sector_cos_upper,
        "outer_sector_cayley_norm_upper": q_sector,
        "validated_sector_boundary_consumed": True,
        "upper_cosine_enclosure_used_for_ungauged_inclusion": True,
        "branches": branches,
        "all_source_handoff_branches_enter_outer_sector": passed,
        "N_outer_words": 0 if passed else None,
        "legacy_microscopic_inner_seed_used_as_outer_capture_target": False,
        "legacy_uniform_transport_route_used": False,
        "P5_OUTER_SECTOR_CAPTURE_CERTIFICATE": "PASS" if passed else "FAIL",
        "P5_CAPTURE_TO_NARROWER_P4_COMPLETE_WORD_CANDIDATE_ESTABLISHED_HERE": False,
        "P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "propagate the 45 deg / 0.5 Hs entrance set through any preceding startup interval, certify complete-word full-state P4 dissipation "
            "on the widest closing candidate sector, then derive a finite P5 word count from the 45 deg entrance to that P4 set and onward to the inner stochastic localization level"
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
    if d.get("P5_entrance_is_distinct_from_P4_complete_word_sector") is not True:
        failures.append("P5 entrance was conflated with P4 complete-word sector")
    if d.get("P5_entrance_set_inside_outer_geometry_sector") is not True:
        failures.append("45 degree P5 entrance is outside broad outer geometry")
    if d.get("P1_source_gauged_handoffs_inside_declared_P5_entrance") is not True:
        failures.append("P1 gauged handoffs are outside declared P5 entrance")
    if d.get("P1_conservative_handoff_box_replaced") is not False:
        failures.append("P1 conservative handoff box was replaced")
    if d.get("startup_propagation_of_entrance_assumed_without_proof") is not False:
        failures.append("startup propagation of P5 entrance was assumed")
    p5e = d.get("declared_P5_entrance", {})
    if p5e.get("attitude_representation") != "SO3_GEODESIC":
        failures.append("P5 entrance is not SO(3) geodesic")
    if p5e.get("componentwise_euler_box_interpretation") is not False:
        failures.append("P5 entrance is an Euler component box")
    if float(p5e.get("gauged_full_attitude_angle_upper_deg", math.nan)) != 45.0:
        failures.append("P5 entrance attitude changed from 45 deg")
    if float(p5e.get("position_component_abs_error_upper_Hs_factor", math.nan)) != 0.5:
        failures.append("P5 entrance position changed from 0.5 Hs")
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
        failures.append("not all source handoff branches enter outer sector")
    if d.get("N_outer_words") != 0:
        failures.append("outer geometry sector is not immediate at entrance")
    if not float(d.get("outer_sector_angle_rad", 0.0)) >= 0.80:
        failures.append("P5 outer geometry sector regressed below 0.80 rad")
    if d.get("legacy_microscopic_inner_seed_used_as_outer_capture_target") is not False:
        failures.append("P5 returned to microscopic inner-seed capture")
    if d.get("legacy_uniform_transport_route_used") is not False:
        failures.append("P5 returned to uniform transport ceiling route")
    if d.get("P5_CAPTURE_TO_NARROWER_P4_COMPLETE_WORD_CANDIDATE_ESTABLISHED_HERE") is not False:
        failures.append("P5 prematurely claims capture to narrower P4 candidate")
    if d.get("P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE") is not False:
        failures.append("outer-sector certificate prematurely claims inner capture")

    branches = d.get("branches", {})
    for name in ("normal_gauged", "timeout_gauged"):
        if branches.get(name, {}).get("inside_declared_P5_entrance") is not True:
            failures.append(f"{name} misses declared P5 entrance")
    ungauged = branches.get("timeout_ungauged", {})
    if ungauged.get("full_heading_radius_assigned") is not False:
        failures.append("ungauged timeout assigned fictitious heading radius")
    if ungauged.get("boundary_cosine_direction_used") != "UPPER_ENCLOSURE":
        failures.append("ungauged timeout used wrong cosine enclosure direction")
    if ungauged.get("inside_outer_gravity_sector") is not True:
        failures.append("ungauged timeout tilt is outside outer sector")
    if ungauged.get("declared_entrance_tilt_inside_outer_sector") is not True:
        failures.append("declared 45 degree ungauged tilt entrance is outside outer sector")
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
        "entrance": out["declared_P5_entrance"],
        "sector_rad": out["outer_sector_angle_rad"],
        "N_outer_words": out["N_outer_words"],
        "branches": out["branches"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
