#!/usr/bin/env python3
"""Identify the OU-III P5 startup-to-inner-funnel capture obligation.

P1 certifies the deployed startup/reset/handoff family and P4 certifies exact
nonlinear H/A source-word decrease only on a very small Cayley-information
sublevel.  P5 is allowed to iterate a source-word recurrence only where the
nonlinear defect estimate used by P4 is itself valid and decreasing.

For H-mode, P4 gives

    sqrt(W_next) <= (1-delta/2) sqrt(W) + B W,

as a conservative consequence of the homogeneous information-word gap and the
transported nonlinear word defect.  A sufficient strict-decrease domain is
therefore

    sqrt(W) < delta/(2 B).

This producer compares that certified domain, and the smaller P4 invariant
inner seed W_*, with the source-declared P1 physical handoff family.  It does
not extrapolate the P4 recurrence outside its proof domain and does not infer
capture from replay.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as P4
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1

# These are the H-mode error-coordinate groups present at goLive.  Each declared
# norm ball contains the axis witness with this group at its radius and the
# other groups zero, so any one witness outside the P4 capture domain disproves
# containment of the complete P1 handoff family.
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
        # W=z^T M z >= m_minus ||z||^2.  For the axis witness ||z||=radius.
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


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P5 operating domain must not be trajectory fitted")

    p1 = P1.build(domain_path)
    p1_fail = P1.validate(p1)
    p4 = P4.build(domain_path)
    p4_fail = P4.validate(p4)
    prereq_failures = [f"P1: {x}" for x in p1_fail] + [f"P4: {x}" for x in p4_fail]
    if prereq_failures:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_STARTUP_CAPTURE_IDENTIFICATION",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "prerequisite_failures": prereq_failures,
            "P5_FINITE_CAPTURE_CERTIFICATE": "NOT_ESTABLISHED",
            "P5_OBSTRUCTION_IDENTIFIED": "NOT_EVALUATED",
            "first_obstruction": "UPSTREAM_P1_OR_P4_FAILURE",
            "failures": prereq_failures,
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

    # From sqrt(W+)<=(1-delta/2)sqrt(W)+B W, strict decrease is guaranteed
    # whenever B sqrt(W)<delta/2.  Use the lower delta and upper B exactly as
    # supplied by P4 and round the threshold downward.
    sqrt_W_capture = P4.div_down(delta, P4.mul_up(2.0, B))
    W_capture = P4.mul_down(sqrt_W_capture, sqrt_W_capture)
    if not (W_capture > W_inner > 0.0):
        raise RuntimeError("P4 inner seed is not inside its derived strict-decrease domain")

    bounds = p1["go_live"]["physical_coordinate_bounds"]
    witnesses = _axis_witnesses(bounds, m_minus, q_design, W_inner, W_capture)
    weakest = min(witnesses, key=lambda x: x["axis_witness_W_lower"])
    strongest = max(witnesses, key=lambda x: x["axis_witness_W_lower"])

    # Both normal and timeout branches share these declared physical-coordinate
    # bounds.  The attitude bounds differ, but translational/bias witnesses are
    # already sufficient to refute containment, so no inverse trig is needed.
    all_outside_capture = all(r["outside_P4_strict_decrease_domain"] for r in witnesses)
    all_outside_design = all(r["outside_P4_nonlinear_design_radius"] for r in witnesses)
    obstruction = all_outside_capture and all_outside_design

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_STARTUP_CAPTURE_IDENTIFICATION",
        "claim": "FIRST_VALIDATED_P5_CAPTURE_OBSTRUCTION_FOR_DEPLOYED_STARTUP_FAMILY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "P1_STARTUP_CERTIFICATE": "PASS",
        "P4_EXACT_NONLINEAR_WORD_CERTIFICATE": "PASS",
        "handoff_modes": ["normal", "timeout"],
        "handoff_coordinate_family": "product of the source-declared goLive physical norm balls",
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
        "normal_handoff_true_gravity_cosine_lower": p1["normal_handoff"]["true_gravity_cosine_lower"],
        "timeout_handoff_true_gravity_cosine_lower": p1["timeout_handoff"]["combined_true_gravity_cosine_lower"],
        "first_required_P5_inequality": "P1_handoff_subset_of_P4_certified_outer_capture_domain",
        "first_required_P5_inequality_holds": not obstruction,
        "first_obstruction": (
            "P1_HANDOFF_OUTSIDE_P4_CERTIFIED_CAPTURE_DOMAIN" if obstruction
            else "NONE_AT_INITIAL_CAPTURE_DOMAIN_GATE"
        ),
        "N_H_words": None if obstruction else "PENDING_RECURRENCE_COUNT",
        "finite_capture_iteration_permitted": not obstruction,
        "reason_iteration_is_not_permitted": (
            "The current P4 nonlinear word bound is not valid/decreasing on the complete P1 handoff family; "
            "iterating the inner recurrence from that family would extrapolate a local certificate outside its proof domain."
            if obstruction else None
        ),
        "required_next_certificate": (
            "widen the same source-reachable Cayley-information word certificate into an outer capture bridge "
            "covering the P1 handoff family, then compute the finite H-word count to W_*"
        ),
        "P5_FINITE_CAPTURE_CERTIFICATE": "NOT_ESTABLISHED" if obstruction else "PENDING_COUNT",
        "P5_OBSTRUCTION_IDENTIFIED": "PASS" if obstruction else "NOT_APPLICABLE",
        "next_obligation": "P5 outer nonlinear capture bridge; do not proceed to a claimed finite N_H until this gate closes",
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
    if d.get("P1_STARTUP_CERTIFICATE") != "PASS":
        failures.append("P1 prerequisite did not pass")
    if d.get("P4_EXACT_NONLINEAR_WORD_CERTIFICATE") != "PASS":
        failures.append("P4 prerequisite did not pass")
    if d.get("P5_OBSTRUCTION_IDENTIFIED") != "PASS":
        failures.append("current P5 first obstruction was not identified")
    if d.get("P5_FINITE_CAPTURE_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("P5 was incorrectly promoted to finite capture")
    if d.get("finite_capture_iteration_permitted") is not False:
        failures.append("P5 permits recurrence iteration outside the certified P4 domain")
    if d.get("first_required_P5_inequality_holds") is not False:
        failures.append("P1 handoff was incorrectly declared inside P4 capture domain")
    if d.get("first_obstruction") != "P1_HANDOFF_OUTSIDE_P4_CERTIFIED_CAPTURE_DOMAIN":
        failures.append("wrong P5 first obstruction")
    W0 = d.get("P4_H_inner_level_W")
    Wcap = d.get("P4_H_strict_decrease_W_threshold_lower")
    if not (isinstance(W0, (int, float)) and isinstance(Wcap, (int, float)) and 0.0 < W0 < Wcap):
        failures.append("invalid P4 inner/capture levels")
    witnesses = d.get("axis_witnesses", [])
    if not witnesses or not all(x.get("outside_P4_strict_decrease_domain") is True for x in witnesses):
        failures.append("P1 axis witnesses do not establish the capture-domain gap")
    if not witnesses or not all(x.get("outside_P4_nonlinear_design_radius") is True for x in witnesses):
        failures.append("P1 axis witnesses do not establish the nonlinear-design-domain gap")
    return failures


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
        "P5_FINITE_CAPTURE_CERTIFICATE": out.get("P5_FINITE_CAPTURE_CERTIFICATE"),
        "P5_OBSTRUCTION_IDENTIFIED": out.get("P5_OBSTRUCTION_IDENTIFIED"),
        "first_obstruction": out.get("first_obstruction"),
        "W_inner": out.get("P4_H_inner_level_W"),
        "W_capture": out.get("P4_H_strict_decrease_W_threshold_lower"),
        "q_design": out.get("P4_H_nonlinear_design_canonical_norm_radius"),
        "weakest_axis_witness": out.get("weakest_axis_witness"),
        "largest_axis_witness": out.get("largest_axis_witness"),
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
