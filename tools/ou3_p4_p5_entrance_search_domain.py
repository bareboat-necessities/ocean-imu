#!/usr/bin/env python3
"""Validated P5 entrance set and reduced P4 complete-word search domain.

This producer separates three objects that must not be conflated:

1. the deployment/P5 entrance assumption;
2. the broad operation-matched finite-angle geometry envelope; and
3. narrower candidate sectors on which the expensive complete 18/21-state P4
   word can first be searched.

The P5 entrance assumption is a 45 degree SO(3) geodesic attitude error on a
gauged branch (tilt only on an ungauged branch) together with
|delta p_i| <= 0.5 Hs for i=x,y,z.  The attitude assumption is deliberately not
an independent Euler-angle box.

The position assumption is a truth-error bound, not a bound on the numerical
estimate itself.  It does not silently replace P1's conservative source handoff
box: if an earlier startup interval precedes P5, that interval still has to
propagate the entrance set source-faithfully.

For conditioning, the translation chain is represented with the dimensionless
coordinates

    p/Hs,  v Ts/Hs,  S/(Hs Ts),  a_w Ts^2/Hs.

Only p receives a new hard bound here.  No unproved sea-scaled hard bounds on
v, S, or a_w are introduced.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p5_cayley_eta_geometry as ETA
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _entrance_attitude_geometry() -> dict:
    # For theta = 45 deg, q = 2 tan(theta/2) = 2(sqrt(2)-1) exactly.
    # Reuse P1's exact-square validated sqrt enclosure rather than libm sqrt in
    # a promoted proof boundary.
    sqrt2 = P1.sqrt_interval_point(2.0)
    q_lo = down(2.0 * down(sqrt2.lo - 1.0))
    q_hi = up(2.0 * up(sqrt2.hi - 1.0))
    cos_lo = down(1.0 / sqrt2.hi)
    cos_hi = up(1.0 / sqrt2.lo)
    return {
        "full_attitude_angle_deg": 45.0,
        "full_attitude_angle_rad_nominal": math.pi / 4.0,
        "cayley_norm_lower": q_lo,
        "cayley_norm_upper": q_hi,
        "cosine_lower": cos_lo,
        "cosine_upper": cos_hi,
        "exact_identity": "q_45 = 2*(sqrt(2)-1)",
        "validated_sqrt2_interval": [sqrt2.lo, sqrt2.hi],
    }


def _candidate_row(deg: float, entrance_q_upper: float) -> dict:
    # Candidate angles are a search strategy, not a deployment theorem bound.
    # Their boundaries are nevertheless evaluated with the same validated
    # finite-angle helper used by the promoted geometry certificate.
    theta = math.radians(float(deg))
    g = SECTOR._validated_design_geometry(theta)
    q = float(g["cayley_norm_upper"])
    monotone = ETA.exact_residual_factor_lower(q)
    eta_ratio = ETA.exact_eta_to_residual_information_ratio_upper(q)
    return {
        "angle_deg": float(deg),
        "angle_rad_nominal": theta,
        "cayley_norm_upper": q,
        "strong_monotonicity_factor_lower": monotone,
        "eta_to_residual_information_ratio_upper": eta_ratio,
        "inside_45deg_entrance": q < entrance_q_upper,
        "angular_span_reduction_vs_45deg": 45.0 / float(deg),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4/P5 entrance domain must not be trajectory fitted")

    entrance = domain["initial_filter_entrance"]
    attitude = entrance["attitude"]
    position = entrance["position"]
    search = domain["certificate_search"]

    failures: list[str] = []
    if attitude.get("representation") != "SO3_GEODESIC":
        failures.append("P5 entrance attitude is not SO(3) geodesic")
    if float(attitude.get("full_attitude_error_upper_deg", math.nan)) != 45.0:
        failures.append("P5 entrance attitude is not 45 degrees")
    if attitude.get("componentwise_euler_box_interpretation") is not False:
        failures.append("P5 entrance was incorrectly turned into an Euler component box")
    if position.get("truth_error_not_estimate_magnitude") is not True:
        failures.append("P5 position entrance is not explicitly a truth-error bound")
    if float(position.get("component_abs_error_upper_Hs_factor", math.nan)) != 0.5:
        failures.append("P5 position component bound is not 0.5 Hs")
    if position.get("significant_wave_height_Hs_positive_required") is not True:
        failures.append("P5 entrance does not require Hs > 0")
    if entrance.get("additional_hard_sea_scaled_bounds_declared") is not False:
        failures.append("unproved sea-scaled v/S/a_w hard bounds were introduced")

    geom = _entrance_attitude_geometry()
    q45 = float(geom["cayley_norm_upper"])
    if not q45 < 1.0:
        failures.append("45 degree entrance left the q<1 Cayley chart")

    sqrt3 = P1.sqrt_interval_point(3.0)
    component_factor = float(position["component_abs_error_upper_Hs_factor"])
    position_norm_factor_upper = up(component_factor * sqrt3.hi)

    legacy_position_norm = float(domain["startup"]["physical_handoff_coordinate_bounds"]["position_error_norm_upper_m"])
    equivalent_Hs_lower = down(legacy_position_norm / position_norm_factor_upper)

    candidates = [float(x) for x in search["p4_complete_word_full_attitude_candidate_deg"]]
    if not candidates or any(not (0.0 < x < 45.0) for x in candidates):
        failures.append("P4 complete-word candidate angles must be strictly inside the 45 degree entrance")
    if any(candidates[i] <= candidates[i + 1] for i in range(len(candidates) - 1)):
        failures.append("P4 complete-word candidate angles must be widest-to-narrowest")
    rows = [_candidate_row(x, q45) for x in candidates]
    if not all(row["inside_45deg_entrance"] for row in rows):
        failures.append("a P4 search candidate is not inside the P5 entrance")

    sea_coords = entrance.get("sea_scaled_translation_coordinates", {})
    expected = {
        "position": "delta_p / Hs",
        "velocity": "delta_v * Ts / Hs",
        "integral_displacement": "delta_S / (Hs * Ts)",
        "latent_acceleration": "delta_a_w * Ts^2 / Hs",
    }
    for key, value in expected.items():
        if sea_coords.get(key) != value:
            failures.append(f"sea-scaled coordinate {key} changed")

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_P5_SEA_SCALED_ENTRANCE_AND_SEARCH_DOMAIN",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "P5_entrance": {
            "attitude_representation": "SO3_GEODESIC",
            "componentwise_euler_box_interpretation": False,
            "gauged_full_attitude_angle_upper_deg": 45.0,
            "ungauged_role": "TILT_BOUND_ONLY_UNTIL_MAGNETIC_REGAUGE",
            "attitude_geometry": geom,
            "position_component_abs_error_upper_Hs_factor": component_factor,
            "position_norm_upper_Hs_factor": position_norm_factor_upper,
            "position_truth_error_bound": True,
            "Hs_positive_required": True,
            "legacy_P1_position_norm_upper_m": legacy_position_norm,
            "Hs_below_this_m_guarantees_smaller_position_norm_than_legacy_P1_box": equivalent_Hs_lower,
        },
        "sea_scaled_translation_coordinates": {
            "definitions": sea_coords,
            "new_hard_bounds": {
                "position_component_abs": "<= 0.5",
                "velocity": None,
                "integral_displacement": None,
                "latent_acceleration": None,
            },
            "purpose": "CONDITION_AND_PARTITION_P4_P5_SEARCH_WITHOUT_INVENTING_NEW_PHYSICAL_BOUNDS",
        },
        "P4_complete_word_search": {
            "candidate_selection_rule": search["p4_candidate_selection_rule"],
            "candidate_rows": rows,
            "outer_geometry_sector_remains_separate": search["p4_outer_geometry_sector_remains_separate"],
            "search_strategy_not_theorem_assumption": search["search_strategy_not_theorem_assumption"],
            "promotion_rule": "PROMOTE_ONLY_THE_WIDEST_CANDIDATE_WITH_SOURCE_COMPLETE_OUTWARD_VALIDATED_FULL_STATE_DISSIPATION",
        },
        "P1_conservative_handoff_box_replaced": False,
        "startup_propagation_of_entrance_assumed_without_proof": False,
        "additional_sea_scaled_v_S_aw_hard_bounds_invented": False,
        "P4_P5_ENTRANCE_SEARCH_DOMAIN_CERTIFICATE": "PASS" if passed else "FAIL",
        "next_obligation": (
            "propagate the declared entrance set through any preceding startup/source interval, then run complete 18/21-state P4 words "
            "from the widest candidate sector downward using the sea-scaled translation coordinates; once a candidate closes, derive P5 finite capture from 45 deg to that P4 set"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("entrance/search certificate is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("entrance/search certificate uses replay")
    if d.get("filter_changed") is not False:
        failures.append("entrance/search certificate changes filter")
    if d.get("P1_conservative_handoff_box_replaced") is not False:
        failures.append("P1 handoff box was silently replaced")
    if d.get("startup_propagation_of_entrance_assumed_without_proof") is not False:
        failures.append("startup propagation of entrance was assumed")
    if d.get("additional_sea_scaled_v_S_aw_hard_bounds_invented") is not False:
        failures.append("unproved sea-scaled bounds were invented")
    p5 = d.get("P5_entrance", {})
    if p5.get("attitude_representation") != "SO3_GEODESIC":
        failures.append("P5 entrance is not SO(3) geodesic")
    if p5.get("componentwise_euler_box_interpretation") is not False:
        failures.append("P5 entrance became an Euler box")
    if float(p5.get("gauged_full_attitude_angle_upper_deg", math.nan)) != 45.0:
        failures.append("P5 entrance angle changed")
    if float(p5.get("position_component_abs_error_upper_Hs_factor", math.nan)) != 0.5:
        failures.append("P5 position factor changed")
    if not 0.86 < float(p5.get("position_norm_upper_Hs_factor", 0.0)) < 0.87:
        failures.append("P5 position norm Hs factor is invalid")
    if not float(p5.get("attitude_geometry", {}).get("cayley_norm_upper", math.inf)) < 1.0:
        failures.append("P5 entrance Cayley radius is invalid")
    rows = d.get("P4_complete_word_search", {}).get("candidate_rows", [])
    if not rows or not all(row.get("inside_45deg_entrance") is True for row in rows):
        failures.append("P4 search ladder is not strictly inside P5 entrance")
    if d.get("P4_complete_word_search", {}).get("search_strategy_not_theorem_assumption") is not True:
        failures.append("P4 search ladder was promoted to theorem assumption")
    if d.get("P4_P5_ENTRANCE_SEARCH_DOMAIN_CERTIFICATE") != "PASS":
        failures.append("P4/P5 entrance search domain did not pass")
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
        "status": out["P4_P5_ENTRANCE_SEARCH_DOMAIN_CERTIFICATE"],
        "P5_entrance": out["P5_entrance"],
        "P4_candidates": out["P4_complete_word_search"]["candidate_rows"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
