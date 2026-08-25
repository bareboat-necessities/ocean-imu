#!/usr/bin/env python3
"""Exact first-due S Cayley-prefix gate for the OU-III P5 outer bridge.

The source-staged first-S calculation already proves a finite full
``S -> attitude`` correction, but ``|d_S|<3`` alone is not enough for P5.  The
correction must be composed with the finite handoff attitude through the actual
deployed quaternion map and must remain in the outer Cayley proof chart.

This producer consumes the exact correction/reset algebra and evaluates

    |c+| <= (|a| + q + |a| q/2) / (1-|a|q/4),

where ``a`` is the Cayley vector of the deployed correction quaternion.  The
bound is outward and uses the source polynomial/axis-angle split through
``ou3_p5_exact_correction_transport``.

It also performs a fail-closed validated bisection for the largest correction
radius that the same formula can prove maps each handoff node into the common
``|c|<1`` outer bootstrap.  That produces a concrete target for tightening the
source-staged ``P_theta``/``K_thetaS`` enclosure; it does not alter the filter or
weaken the theorem gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_exact_correction_transport as TRANSPORT
import ou3_p5_first_s_gain_certificate as FIRST
import ou3_p5_first_s_state_prefix_certificate as SPREFIX
import ou3_p5_heading_handoff_contract as HEADING

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
TARGET_OUTER_CAYLEY_NORM = 1.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def add_up(a: float, b: float) -> float:
    return up(float(a) + float(b))


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def div_up(a: float, b: float) -> float:
    if not b > 0.0:
        raise RuntimeError("positive denominator required")
    return up(float(a) / float(b))


def post_cayley_norm_upper(q: float, delta: float) -> dict:
    tr = TRANSPORT.reset_defect_bound(q, delta)
    a = float(tr["injected_cayley_norm_upper"])
    denom = float(tr["cayley_composition_denominator_lower"])
    numerator = add_up(a, q)
    numerator = add_up(numerator, mul_up(0.5 * a, q))
    post = div_up(numerator, denom)
    return {
        **tr,
        "post_injection_cayley_norm_upper": post,
        "inside_common_outer_bootstrap": post < TARGET_OUTER_CAYLEY_NORM,
    }


def certified_delta_for_target(q: float, target: float = TARGET_OUTER_CAYLEY_NORM) -> float:
    """Largest conservative binary64 delta known to satisfy the outward gate."""
    if not (0.0 <= q < target < 2.0):
        raise ValueError("bisection requires 0<=q<target<2")
    lo = 0.0
    hi = min(2.9, math.pi - 1.0e-6)
    # The predicate is evaluated with outward source-faithful Cayley bounds.
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        try:
            row = post_cayley_norm_upper(q, mid)
            ok = row["post_injection_cayley_norm_upper"] < target
        except RuntimeError:
            ok = False
        if ok:
            lo = mid
        else:
            hi = mid
    # Move downward once so the returned design radius is strictly on the
    # certified side even if the final midpoint rounded onto the transition.
    return down(lo)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("first-S exact prefix domain must not be trajectory fitted")

    tr = TRANSPORT.build(domain_path)
    first = FIRST.build(domain_path)
    sp = SPREFIX.build(domain_path)
    heading = HEADING.build(domain_path)
    failures = [f"transport: {x}" for x in TRANSPORT.validate(tr)]
    failures += [f"first-S-gain: {x}" for x in FIRST.validate(first)]
    failures += [f"first-S-state: {x}" for x in SPREFIX.validate(sp)]
    failures += [f"heading: {x}" for x in HEADING.validate(heading)]

    current_delta = float(sp["first_due_S_induced_attitude_correction_norm_upper_rad"])
    current_ptheta = float(first["P_theta_theta_directional_lambda_max_upper"])
    nodes = {}
    for name, row in (
        ("normal_gauged", heading["gauged_quality_handoff"]),
        ("timeout_gauged", heading["gauged_timeout_subbranch"]),
    ):
        q = float(row["full_attitude_cayley_norm_upper"])
        cur = post_cayley_norm_upper(q, current_delta)
        required_delta = certified_delta_for_target(q)
        # FIRST's KthetaS is proportional to sqrt(Ptheta_upper) with every other
        # source-staged factor fixed.  This is a proof-design target, not a
        # substituted covariance claim.
        ptheta_target = down(current_ptheta * (required_delta / current_delta) ** 2)
        cur.update({
            "certified_correction_radius_for_cayley_lt_1_rad": required_delta,
            "current_correction_over_required_factor": up(current_delta / required_delta),
            "directional_Ptheta_upper_target_if_other_first_S_factors_unchanged": ptheta_target,
            "current_directional_Ptheta_upper": current_ptheta,
            "Ptheta_tightening_factor_required": up(current_ptheta / ptheta_target),
        })
        nodes[name] = cur

    all_inside = all(row["inside_common_outer_bootstrap"] for row in nodes.values())
    first_failure = (
        "NONE" if all_inside else
        "FIRST_DUE_S_EXACT_CAYLEY_PREFIX_NOT_CERTIFIED_WITH_CURRENT_STAGED_BOUND"
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_FIRST_DUE_S_CAYLEY_PREFIX_GATE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "exact_deployed_quaternion_cayley_composition_used": True,
        "full_S_to_attitude_gain_retained": True,
        "common_outer_cayley_target": TARGET_OUTER_CAYLEY_NORM,
        "current_first_due_S_correction_norm_upper_rad": current_delta,
        "nodes": nodes,
        "current_first_S_prefix_inside_outer_bootstrap": all_inside,
        "P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE": "PASS" if all_inside and not failures else "NOT_ESTABLISHED",
        "first_failure": first_failure if not failures else "UPSTREAM_PREREQUISITE_FAILURE",
        "next_obligation": (
            "tighten the source-staged pre-first-S directional attitude covariance/gain enclosure until the timeout node proves d_S below its certified correction target; then re-evaluate exact sequential vector/S prefixes"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("first-S exact prefix is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("first-S exact prefix uses replay")
    if d.get("filter_changed") is not False:
        failures.append("first-S exact prefix changes filter")
    if d.get("exact_deployed_quaternion_cayley_composition_used") is not True:
        failures.append("first-S exact prefix does not use deployed Cayley composition")
    if d.get("full_S_to_attitude_gain_retained") is not True:
        failures.append("first-S exact prefix drops S-to-attitude gain")
    for name in ("normal_gauged", "timeout_gauged"):
        row = d.get("nodes", {}).get(name, {})
        if not float(row.get("cayley_composition_denominator_lower", -1.0)) > 0.0:
            failures.append(f"{name}: exact correction reaches Cayley antipode")
        target = row.get("certified_correction_radius_for_cayley_lt_1_rad")
        if not (isinstance(target, (int, float)) and math.isfinite(float(target)) and float(target) > 0.0):
            failures.append(f"{name}: correction target is not positive")
        ptarget = row.get("directional_Ptheta_upper_target_if_other_first_S_factors_unchanged")
        if not (isinstance(ptarget, (int, float)) and math.isfinite(float(ptarget)) and float(ptarget) > 0.0):
            failures.append(f"{name}: Ptheta design target is not positive")
    # A fail-closed numerical obstruction is a valid producer output.  Validation
    # checks that it is named rather than demanding a false PASS.
    if d.get("current_first_S_prefix_inside_outer_bootstrap") is False:
        if d.get("P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE") != "NOT_ESTABLISHED":
            failures.append("failed first-S prefix was promoted")
        if d.get("first_failure") != "FIRST_DUE_S_EXACT_CAYLEY_PREFIX_NOT_CERTIFIED_WITH_CURRENT_STAGED_BOUND":
            failures.append("failed first-S prefix obstruction not named")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"],
        "first_failure": out["first_failure"],
        "current_delta": out["current_first_due_S_correction_norm_upper_rad"],
        "nodes": out["nodes"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
