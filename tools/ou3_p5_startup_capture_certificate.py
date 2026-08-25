#!/usr/bin/env python3
"""Identify the OU-III P5 startup-to-inner-funnel capture obstruction.

P4 proves exact nonlinear source-word decrease only on a very small inner
Cayley-information sublevel.  This producer proves that the source-declared P1
handoff family is not contained in that local decrease domain, so the P4
recurrence may not be extrapolated outward.

Attitude bookkeeping is source faithful.  P1's gravity cosines are tilt-only.
Full-heading Cayley radii are taken only from the separately verified gauged
handoff contract; the timeout branch without a yaw gauge is explicitly routed
to the gravity/yaw quotient instead of being assigned a fictitious SO(3)
radius.  This file is an obstruction/requirements identifier, not the P5
completion certificate; ``ou3_p5_outer_h_bridge_certificate`` is the current
completion object.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as P4
import ou3_p5_heading_handoff_contract as HEADING
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3

HANDOFF_GROUPS = (
    ("b_g", "gyro_bias_error_norm_upper_rad_s"),
    ("v", "velocity_error_norm_upper_mps"),
    ("p", "position_error_norm_upper_m"),
    ("S", "integral_displacement_error_norm_upper_m_s"),
    ("a_w", "latent_acceleration_error_norm_upper_mps2"),
)


def _axis_witnesses(bounds: dict, m_minus: float, q_design: float,
                    W_inner: float, W_capture: float) -> list[dict]:
    rows = []
    for group, key in HANDOFF_GROUPS:
        radius = float(bounds[key])
        if not (math.isfinite(radius) and radius >= 0.0):
            raise RuntimeError(f"invalid P1 handoff bound {key}")
        r2 = P4.mul_down(radius, radius)
        W_lower = P4.mul_down(m_minus, r2)
        rows.append({
            "group": group,
            "domain_field": key,
            "axis_witness_canonical_norm": radius,
            "axis_witness_W_lower": W_lower,
            "outside_P4_nonlinear_design_radius": radius > q_design,
            "outside_P4_inner_seed": W_lower > W_inner,
            "outside_P4_strict_decrease_domain": W_lower > W_capture,
            "W_lower_over_inner_seed": W_lower / W_inner,
            "W_lower_over_strict_decrease_threshold": W_lower / W_capture,
        })
    return rows


def _cayley_norm_upper_from_cos_lower(cos_lower: float) -> float:
    """Tilt-only Cayley witness retained for the retired outer diagnostic.

    P1's input cosine is a gravity-direction/tilt cosine.  The returned value
    must not be interpreted as a full-heading handoff radius.  It remains only
    so ``ou3_p5_outer_h_word_certificate`` can demonstrate that its retired
    perturbative route already fails on an optimistic smaller attitude witness.
    The actual P5 full-heading nodes come from ``HEADING.build()`` below.
    """
    c = float(cos_lower)
    if not (-1.0 < c <= 1.0):
        raise RuntimeError("invalid strict tilt cosine lower bound")
    numerator = P4.up(1.0 - c)
    denominator = P4.down(1.0 + c)
    if not denominator > 0.0:
        raise RuntimeError("tilt witness reaches Cayley singularity")
    return P4.mul_up(2.0, P4.sqrt_up(P4.div_up(numerator, denominator)))


def _uniform_recurrence_B_cap(delta: float, W_lower: float) -> float:
    if not (delta > 0.0 and W_lower > 0.0):
        raise RuntimeError("positive delta and witness W lower bound required")
    return P4.div_down(delta, P4.mul_up(2.0, P4.sqrt_up(W_lower)))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P5 operating domain must not be trajectory fitted")

    p1 = P1.build(domain_path)
    p4 = P4.build(domain_path)
    heading = HEADING.build(domain_path)
    prereq = [f"P1: {x}" for x in P1.validate(p1)]
    prereq += [f"P4: {x}" for x in P4.validate(p4)]
    prereq += [f"heading: {x}" for x in HEADING.validate(heading)]
    if prereq:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_STARTUP_CAPTURE_IDENTIFICATION",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "P5_FINITE_CAPTURE_CERTIFICATE": "NOT_ESTABLISHED",
            "P5_OBSTRUCTION_IDENTIFIED": "NOT_EVALUATED",
            "P5_OUTER_BRIDGE_REQUIREMENTS_IDENTIFIED": "NOT_EVALUATED",
            "first_obstruction": "UPSTREAM_P1_OR_P4_OR_HEADING_FAILURE",
            "failures": prereq,
        }

    H = p4["modes"]["H"]
    delta = float(H["P3_word_endpoint_delta_lower"])
    B = float(H["transported_word_defect_B_upper"])
    m_minus = float(H["metric_lambda_min_lower"])
    W_inner = float(H["certified_level_W"])
    sqrt_W_inner = float(H["certified_level_sqrt_W"])
    q_design = float(H["correction_quadratic_bound"]["design_error_norm_radius"])
    if not (0.0 < delta < 1.0 and math.isfinite(B) and B > 0.0 and m_minus > 0.0):
        raise RuntimeError("P4 H recurrence constants are invalid")

    sqrt_W_capture = P4.div_down(delta, P4.mul_up(2.0, B))
    W_capture = P4.mul_down(sqrt_W_capture, sqrt_W_capture)
    if not (W_capture > W_inner > 0.0):
        raise RuntimeError("P4 inner seed is not inside derived strict-decrease domain")

    bounds = p1["go_live"]["physical_coordinate_bounds"]
    witnesses = _axis_witnesses(bounds, m_minus, q_design, W_inner, W_capture)
    weakest = min(witnesses, key=lambda x: x["axis_witness_W_lower"])
    strongest = max(witnesses, key=lambda x: x["axis_witness_W_lower"])

    weak_B_cap = _uniform_recurrence_B_cap(delta, float(weakest["axis_witness_W_lower"]))
    strong_B_cap = _uniform_recurrence_B_cap(delta, float(strongest["axis_witness_W_lower"]))
    obstruction = all(r["outside_P4_strict_decrease_domain"] for r in witnesses) and all(
        r["outside_P4_nonlinear_design_radius"] for r in witnesses
    )

    qn = float(heading["gauged_quality_handoff"]["full_attitude_cayley_norm_upper"])
    qt = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    promoted_limit = float(P4.PROMOTED_CAYLEY_NORM_LIMIT)
    outer = {
        "P1_gravity_cosines_are_tilt_only": True,
        "normal_gauged_full_attitude_cayley_norm_upper": qn,
        "timeout_gauged_full_attitude_cayley_norm_upper": qt,
        "timeout_ungauged_full_heading_cayley_bound_available": False,
        "timeout_ungauged_required_route": heading["ungauged_timeout_subbranch"]["required_route"],
        "normal_gauged_over_current_P4_design_radius_factor": qn / q_design,
        "timeout_gauged_over_current_P4_design_radius_factor": qt / q_design,
        "current_P4_promoted_cayley_norm_limit": promoted_limit,
        "normal_gauged_inside_current_promoted_cayley_norm_limit": qn < promoted_limit,
        "timeout_gauged_inside_current_promoted_cayley_norm_limit": qt < promoted_limit,
        "current_uniform_full_state_B_upper": B,
        "optimistic_uniform_B_cap_at_weakest_P1_axis_witness": weak_B_cap,
        "optimistic_uniform_B_cap_at_largest_P1_axis_witness": strong_B_cap,
        "uniform_B_reduction_factor_needed_at_weakest_witness": B / weak_B_cap,
        "uniform_B_reduction_factor_needed_at_largest_witness": B / strong_B_cap,
        "interpretation": (
            "The gauged P1 attitude nodes remain inside the broad Cayley chart but are many orders larger than q_design. "
            "Simply enlarging q_design while retaining the isotropic B*W perturbation recurrence is not a P5 bridge. "
            "The ungauged timeout branch is not a full-heading node at all and must use the yaw quotient until gauge acquisition."
        ),
        "required_proof_structure": [
            "branch-specific exact SO(3) finite-angle dissipation on gauged normal/timeout nodes",
            "gravity-only yaw-quotient capture on the ungauged timeout branch until a magnetic gauge hybrid event",
            "source-staged early covariance, pseudo-phase, and S-to-attitude bounds",
            "prefix-safe exact quaternion correction coverage",
            "validated outer H decrease/funnel recursion overlapping the existing P4 inner seed",
        ],
    }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_STARTUP_CAPTURE_IDENTIFICATION",
        "claim": "LOCAL_P4_CAPTURE_OBSTRUCTION_AND_SOURCE_FAITHFUL_OUTER_REQUIREMENTS",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "P1_STARTUP_CERTIFICATE": "PASS",
        "P4_EXACT_NONLINEAR_WORD_CERTIFICATE": "PASS",
        "heading_branch_contract": "PASS",
        "handoff_modes": ["normal_gauged", "timeout_gauged", "timeout_ungauged_yaw_quotient"],
        "handoff_coordinate_family": "product of source-declared goLive physical norm balls with branch-correct attitude quotient",
        "H_word_horizon_s": H["word_horizon_s"],
        "P4_H_inner_level_W": W_inner,
        "P4_H_inner_level_sqrt_W": sqrt_W_inner,
        "P4_H_delta_lower": delta,
        "P4_H_transported_word_defect_B_upper": B,
        "P4_H_metric_lambda_min_lower": m_minus,
        "P4_H_nonlinear_design_canonical_norm_radius": q_design,
        "P4_H_strict_decrease_sqrt_W_threshold_lower": sqrt_W_capture,
        "P4_H_strict_decrease_W_threshold_lower": W_capture,
        "P4_inner_seed_to_decrease_threshold_W_factor": W_capture / W_inner,
        "axis_witnesses": witnesses,
        "weakest_axis_witness": weakest,
        "largest_axis_witness": strongest,
        "outer_bridge_requirements": outer,
        "first_required_P5_inequality": "P1_handoff_subset_of_P4_certified_outer_capture_domain",
        "first_required_P5_inequality_holds": not obstruction,
        "first_obstruction": "P1_HANDOFF_OUTSIDE_P4_CERTIFIED_CAPTURE_DOMAIN" if obstruction else "NONE_AT_INITIAL_CAPTURE_DOMAIN_GATE",
        "N_H_words": None if obstruction else "PENDING_RECURRENCE_COUNT",
        "finite_capture_iteration_permitted": not obstruction,
        "reason_iteration_is_not_permitted": (
            "The local P4 nonlinear word bound is not valid/decreasing on the complete P1 handoff family; iterating it would extrapolate outside its proof domain."
            if obstruction else None
        ),
        "required_next_certificate": (
            "use the staged outer-H bridge: exact early S/covariance bounds, exact large-angle gauged-vector dissipation, "
            "and a yaw-quotient timeout route; only then compute finite H-word capture into W_*"
        ),
        "P5_FINITE_CAPTURE_CERTIFICATE": "NOT_ESTABLISHED" if obstruction else "PENDING_COUNT",
        "P5_OBSTRUCTION_IDENTIFIED": "PASS" if obstruction else "NOT_APPLICABLE",
        "P5_OUTER_BRIDGE_REQUIREMENTS_IDENTIFIED": "PASS" if obstruction else "NOT_APPLICABLE",
        "completion_object": "ou3_p5_outer_h_bridge_certificate.py",
        "next_obligation": "P5 staged outer H bridge; no finite N_H before both gauged and quotient routes close",
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("P5 identification is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("P5 identification uses replay")
    if d.get("P1_STARTUP_CERTIFICATE") != "PASS" or d.get("P4_EXACT_NONLINEAR_WORD_CERTIFICATE") != "PASS":
        failures.append("P1/P4 prerequisite did not pass")
    if d.get("heading_branch_contract") != "PASS":
        failures.append("heading branch contract did not pass")
    if d.get("P5_OBSTRUCTION_IDENTIFIED") != "PASS":
        failures.append("current P5 local-capture obstruction not identified")
    if d.get("P5_OUTER_BRIDGE_REQUIREMENTS_IDENTIFIED") != "PASS":
        failures.append("P5 outer requirements not identified")
    if d.get("P5_FINITE_CAPTURE_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("P5 identifier incorrectly promoted finite capture")
    if d.get("finite_capture_iteration_permitted") is not False:
        failures.append("P5 identifier permits local recurrence extrapolation")
    if d.get("first_required_P5_inequality_holds") is not False:
        failures.append("P1 handoff incorrectly declared inside local P4 capture")
    if d.get("first_obstruction") != "P1_HANDOFF_OUTSIDE_P4_CERTIFIED_CAPTURE_DOMAIN":
        failures.append("wrong local P5 obstruction")
    W0, Wcap = d.get("P4_H_inner_level_W"), d.get("P4_H_strict_decrease_W_threshold_lower")
    if not (isinstance(W0, (int, float)) and isinstance(Wcap, (int, float)) and 0.0 < W0 < Wcap):
        failures.append("invalid P4 inner/capture levels")
    outer = d.get("outer_bridge_requirements", {})
    if outer.get("P1_gravity_cosines_are_tilt_only") is not True:
        failures.append("gravity cosine still treated as full attitude")
    if outer.get("timeout_ungauged_full_heading_cayley_bound_available") is not False:
        failures.append("ungauged timeout assigned full-heading radius")
    if "YAW_QUOTIENT" not in str(outer.get("timeout_ungauged_required_route", "")):
        failures.append("ungauged timeout not routed to yaw quotient")
    if d.get("completion_object") != "ou3_p5_outer_h_bridge_certificate.py":
        failures.append("local obstruction identifier masquerades as completion object")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P5_FINITE_CAPTURE_CERTIFICATE": d.get("P5_FINITE_CAPTURE_CERTIFICATE"),
        "first_obstruction": d.get("first_obstruction"),
        "outer_bridge_requirements": d.get("outer_bridge_requirements"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
