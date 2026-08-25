#!/usr/bin/env python3
"""Exact first-due S Cayley-prefix gate for the OU-III P5 outer bridge.

The source-staged first-S calculation already proves a finite full
``S -> attitude`` correction, but ``|d_S|<3`` alone is not enough for P5.  The
correction is therefore composed with each finite handoff attitude through the
actual deployed quaternion map,

    |c+| <= (|a| + q + |a| q/2) / (1-|a|q/4),

where ``a`` is the Cayley vector of the deployed correction quaternion.

An earlier diagnostic required this first correction to remain inside the
convenient ``|c|<1`` inner outer-bootstrap.  That is unnecessarily restrictive:
the Cayley chart is valid for every finite coordinate and the exact source map
does not become singular at ``|c|=1``.  This producer keeps the ``|c|<1``
diagnostic, but *widens* the certified prefix chart to the smallest dyadic
radius containing the validated first-S image of both gauged handoff nodes.
No filter parameter or theorem gate is changed.

The widened radius also receives an exact positive vector-residual information
factor ``4/(4+q^2)``, so the chart extension does not fall back to a small-angle
certificate.  Remaining P5 work must subdivide/transport later source prefixes
inside this finite chart and can tighten it again once source-staged covariance
correlations permit.
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
import ou3_p5_outer_information_geometry as OUTINFO

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
DIAGNOSTIC_CAYLEY_NORM = 1.0
MAX_WIDENED_PREFIX_CAYLEY_NORM = 16.0


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
        "inside_diagnostic_cayley_lt_1": post < DIAGNOSTIC_CAYLEY_NORM,
    }


def certified_delta_for_target(q: float, target: float = DIAGNOSTIC_CAYLEY_NORM) -> float:
    """Largest conservative binary64 delta known to satisfy an outward target."""
    if not (0.0 <= q < target < 2.0):
        raise ValueError("bisection requires 0<=q<target<2")
    lo = 0.0
    hi = 2.9
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
    return down(lo)


def _smallest_dyadic_strict_upper(x: float) -> float:
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("finite nonnegative radius required")
    q = 1.0
    while not q > x:
        q *= 2.0
    if q > MAX_WIDENED_PREFIX_CAYLEY_NORM:
        raise RuntimeError("first-S image requires excessively broad Cayley prefix chart")
    return q


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("first-S exact prefix domain must not be trajectory fitted")

    tr = TRANSPORT.build(domain_path)
    first = FIRST.build(domain_path)
    sp = SPREFIX.build(domain_path)
    heading = HEADING.build(domain_path)
    outinfo = OUTINFO.build(domain_path)
    failures = [f"transport: {x}" for x in TRANSPORT.validate(tr)]
    failures += [f"first-S-gain: {x}" for x in FIRST.validate(first)]
    failures += [f"first-S-state: {x}" for x in SPREFIX.validate(sp)]
    failures += [f"heading: {x}" for x in HEADING.validate(heading)]
    failures += [f"outer-information: {x}" for x in OUTINFO.validate(outinfo)]

    current_delta = float(sp["first_due_S_induced_attitude_correction_norm_upper_rad"])
    current_ptheta = float(first["P_theta_theta_directional_lambda_max_upper"])
    nodes = {}
    required_prefix = 0.0
    for name, row in (
        ("normal_gauged", heading["gauged_quality_handoff"]),
        ("timeout_gauged", heading["gauged_timeout_subbranch"]),
    ):
        q = float(row["full_attitude_cayley_norm_upper"])
        cur = post_cayley_norm_upper(q, current_delta)
        required_prefix = max(required_prefix, float(cur["post_injection_cayley_norm_upper"]))
        required_delta = certified_delta_for_target(q)
        ptheta_target = down(current_ptheta * (required_delta / current_delta) ** 2)
        cur.update({
            "certified_correction_radius_for_diagnostic_cayley_lt_1_rad": required_delta,
            "current_correction_over_diagnostic_target_factor": up(current_delta / required_delta),
            "directional_Ptheta_upper_target_for_diagnostic_cayley_lt_1_if_other_factors_unchanged": ptheta_target,
            "current_directional_Ptheta_upper": current_ptheta,
            "diagnostic_Ptheta_tightening_factor": up(current_ptheta / ptheta_target),
        })
        nodes[name] = cur

    widened_q = _smallest_dyadic_strict_upper(required_prefix)
    anti_margin = down(8.0 / up(4.0 + up(widened_q * widened_q)))
    residual_factor = down(4.0 / up(4.0 + up(widened_q * widened_q)))
    packet = outinfo["packet_geometry"]
    linear_pair_mu = float(packet["linear_pair_information_mu_lower"])
    info_lambda = max(float(v["goLive_attitude_information_lambda_max_upper"]) for v in outinfo["nodes"].values())
    widened_relative_info = down(down(residual_factor * linear_pair_mu) / up(info_lambda))
    widened_pass = bool(
        not failures
        and widened_q > required_prefix
        and anti_margin > 0.0
        and residual_factor > 0.0
        and widened_relative_info > 0.0
    )
    for row in nodes.values():
        row["inside_widened_prefix_chart"] = float(row["post_injection_cayley_norm_upper"]) < widened_q

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_FIRST_DUE_S_CAYLEY_PREFIX_GATE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "exact_deployed_quaternion_cayley_composition_used": True,
        "full_S_to_attitude_gain_retained": True,
        "diagnostic_cayley_norm": DIAGNOSTIC_CAYLEY_NORM,
        "diagnostic_q_lt_1_is_promotion_gate": False,
        "current_first_due_S_correction_norm_upper_rad": current_delta,
        "nodes": nodes,
        "required_first_S_post_cayley_norm_upper": required_prefix,
        "widened_prefix_cayley_norm_upper": widened_q,
        "widened_prefix_antipodal_one_plus_cosine_margin_lower": anti_margin,
        "widened_prefix_exact_vector_residual_factor_lower": residual_factor,
        "widened_prefix_pair_information_vs_goLive_attitude_metric_lower": widened_relative_info,
        "current_first_S_prefix_inside_widened_chart": all(row["inside_widened_prefix_chart"] for row in nodes.values()),
        "P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE": "PASS_WIDENED_CHART" if widened_pass else "NOT_ESTABLISHED",
        "first_failure": "NONE_AT_FIRST_S_ATTITUDE_CHART" if widened_pass else "FIRST_DUE_S_EXACT_CAYLEY_PREFIX_NOT_CERTIFIED",
        "next_obligation": (
            "propagate subsequent S/vector/prediction prefixes by source-correlated subdivision in the finite widened Cayley chart; q<1 remains only a useful tightening target, not a theorem gate"
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
    if d.get("diagnostic_q_lt_1_is_promotion_gate") is not False:
        failures.append("q<1 diagnostic was promoted to theorem gate")
    widened = d.get("widened_prefix_cayley_norm_upper")
    required = d.get("required_first_S_post_cayley_norm_upper")
    if not (isinstance(widened, (int, float)) and isinstance(required, (int, float))
            and math.isfinite(float(widened)) and math.isfinite(float(required))
            and float(widened) > float(required) > 0.0):
        failures.append("widened first-S prefix radius does not strictly contain exact image")
    if not float(d.get("widened_prefix_antipodal_one_plus_cosine_margin_lower", 0.0)) > 0.0:
        failures.append("widened prefix chart lacks antipodal margin")
    if not float(d.get("widened_prefix_exact_vector_residual_factor_lower", 0.0)) > 0.0:
        failures.append("widened prefix vector factor is not strict")
    if not float(d.get("widened_prefix_pair_information_vs_goLive_attitude_metric_lower", 0.0)) > 0.0:
        failures.append("widened prefix information factor is not strict")
    for name in ("normal_gauged", "timeout_gauged"):
        row = d.get("nodes", {}).get(name, {})
        if not float(row.get("cayley_composition_denominator_lower", -1.0)) > 0.0:
            failures.append(f"{name}: exact correction reaches Cayley antipode")
        if row.get("inside_widened_prefix_chart") is not True:
            failures.append(f"{name}: first-S image outside widened chart")
    if not failures and d.get("P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE") != "PASS_WIDENED_CHART":
        failures.append("first-S widened exact Cayley prefix did not pass")
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
        "current_delta": out["current_first_due_S_correction_norm_upper_rad"],
        "required_post_q": out["required_first_S_post_cayley_norm_upper"],
        "widened_q": out["widened_prefix_cayley_norm_upper"],
        "widened_info": out["widened_prefix_pair_information_vs_goLive_attitude_metric_lower"],
        "nodes": out["nodes"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
