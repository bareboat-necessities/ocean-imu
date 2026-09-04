#!/usr/bin/env python3
"""Validated SEA3 rotational-response coupling for the P1 body-rate good event.

The #482 translational RAO family deliberately left parent rotational motion
under the separate Normal-Live body-rate assumption.  This module supplies the
missing analytical bridge without pretending that a particular hull has been
identified.

Let r(f,theta) be the complex three-axis *rotational displacement* response in
rad/m from sea-surface elevation to body attitude, and suppose a candidate
physical vessel satisfies

    ||r(f,theta)||_2 <= K min(1,(f_c/f)^q),   f>0, q>=1.

Then the body-rate response is i 2 pi f r eta.  For q>=1,

    (2 pi f)^2 ||r||^2 <= K^2 (2 pi f_c)^2

at every frequency and heading.  Since the total elevation variance is
m0=Hs^2/16,

    tr Cov[omega_rad/s] <= Hs^2 K^2 (2 pi f_c)^2 / 16.

After exact rad/s -> deg/s conversion the pi factors cancel:

    tr Cov[omega_deg/s] <= 8100 Hs^2 K^2 f_c^2.

This is a spectrum-shape-independent bound.  It therefore covers all admitted
JONSWAP gamma values, arbitrary directional spreading, arbitrary complex phase,
and cross-axis coupling once the rotational envelope itself has been validated
for the physical vessel.  The module does *not* declare a universal rotational
RAO box.  It exposes the exact monotone coupling predicate that a vessel or a
tighter response-weighted certificate must satisfy.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from ou3_interval import Interval, down
import ou3_sea3_finite_horizon_concentration as CONC
import ou3_sea3_physical_admissibility as PHYS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RESPONSE_DOMAIN = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_BODY_RATE_COVARIANCE_COUPLING_V1"
DEGREE_COEFFICIENT = 8100.0


def _outward_point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def body_rate_trace_covariance_upper_deg2_s2(
    hs_m: float,
    rotation_gain_rad_per_m: float,
    corner_hz: float,
    rolloff_power: float,
) -> float:
    """Outward upper bound ``8100 Hs^2 K^2 fc^2`` for q>=1."""
    h = float(hs_m)
    gain = float(rotation_gain_rad_per_m)
    corner = float(corner_hz)
    power = float(rolloff_power)
    if not (math.isfinite(h) and h >= 0.0):
        raise ValueError("Hs must be finite and nonnegative")
    if not (math.isfinite(gain) and gain >= 0.0):
        raise ValueError("rotational displacement gain must be finite and nonnegative")
    if not (math.isfinite(corner) and corner > 0.0):
        raise ValueError("rotational rolloff corner must be finite and positive")
    if not (math.isfinite(power) and power >= 1.0):
        raise ValueError("body-rate covariance theorem requires rotational rolloff q>=1")

    hs = _outward_point(h)
    k = _outward_point(gain)
    fc = _outward_point(corner)
    coeff = Interval.point(DEGREE_COEFFICIENT)
    return (coeff * hs.square() * k.square() * fc.square()).hi


def evaluate_tuple(
    certificate: dict,
    *,
    hs_m: float,
    rotation_gain_rad_per_m: float,
    corner_hz: float,
    rolloff_power: float,
) -> dict[str, Any]:
    """Evaluate one proposed coupled sea/rotational-response tuple."""
    failures = validate(certificate)
    if failures:
        raise ValueError(f"invalid body-rate coupling certificate: {failures}")

    reasons: list[str] = []
    h = float(hs_m)
    if not (math.isfinite(h) and 0.0 <= h <= float(certificate["repository_total_Hs_upper_m"])):
        reasons.append("Hs lies outside the declared SEA3 total-height envelope")

    try:
        bound = body_rate_trace_covariance_upper_deg2_s2(
            h,
            rotation_gain_rad_per_m,
            corner_hz,
            rolloff_power,
        )
    except (TypeError, ValueError) as exc:
        reasons.append(str(exc))
        bound = math.inf

    threshold = float(certificate["required_body_rate_trace_covariance_upper_deg2_s2"])
    if math.isfinite(bound) and bound > threshold:
        reasons.append(
            "coupled Hs/rotational-gain/corner tuple exceeds finite-horizon body-rate covariance threshold"
        )

    passed = not reasons
    return {
        "candidate_pass": passed,
        "decision": "PASS_CANDIDATE" if passed else "FAIL_CANDIDATE",
        "H_s_m": h,
        "rotation_gain_rad_per_m": float(rotation_gain_rad_per_m),
        "rotation_corner_hz": float(corner_hz),
        "rotation_rolloff_power": float(rolloff_power),
        "validated_body_rate_trace_covariance_upper_deg2_s2": bound,
        "required_body_rate_trace_covariance_upper_deg2_s2": threshold,
        "margin_deg2_s2": threshold - bound if math.isfinite(bound) else -math.inf,
        "validation_failures": reasons,
        "physical_vessel_rotational_RAO_qualified": False,
        "finite_horizon_good_event_promoted": False,
        "deterministic_left_inclusion_promoted": False,
    }


def build(
    samples: int,
    domain_path: Path = DEFAULT_DOMAIN,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
    repo: Path = REPO,
) -> dict[str, Any]:
    concentration = CONC.build(samples, domain_path, response_domain_path, repo)
    cf = CONC.validate(concentration)
    if cf:
        raise RuntimeError(f"invalid finite-horizon concentration prerequisite: {cf}")

    physical = PHYS.build(Path(domain_path))
    pf = PHYS.validate(physical)
    if pf:
        raise RuntimeError(f"invalid SEA3 physical admissibility prerequisite: {pf}")

    threshold = float(
        concentration["body_rate"]["required_trace_covariance_upper_deg2_s2"]
    )
    # This is the square of the exact monotone product bound
    # Hs*K*fc <= sqrt(threshold/8100).  Keep the proof quantity squared so no
    # validated square-root implementation is required for the decision.
    product_squared_threshold = down(threshold / DEGREE_COEFFICIENT)

    out: dict[str, Any] = {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "samples": samples,
        "trajectory_replay_used": False,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "ordinary_libm_transcendental_used_in_pass_decision": False,
        "repository_total_Hs_upper_m": float(physical["repository_total_Hs_upper_m"]),
        "required_body_rate_trace_covariance_upper_deg2_s2": threshold,
        "response_parameter_box_sha256": concentration["response_parameter_box_sha256"],
        "rotational_response_interface": {
            "definition": "||r(f,theta)||_2 <= K min(1,(f_c/f)^q), q>=1",
            "r_units": "rad per metre of sea-surface elevation",
            "body_rate_response": "omega_hat = i*2*pi*f*r*eta_hat",
            "rolloff_power_min": 1.0,
            "universal_rotational_parameter_box_declared": False,
            "actual_vessel_or_response_family_must_supply_K_fc_q": True,
        },
        "uniform_covariance_theorem": {
            "m0_identity": "m0 = Hs^2 / 16",
            "frequency_weight_bound": "(2*pi*f)^2*||r||^2 <= K^2*(2*pi*f_c)^2",
            "rad_covariance_bound": "tr Cov[omega] <= Hs^2*K^2*(2*pi*f_c)^2/16",
            "degree_covariance_bound": "tr Cov[omega_deg_s] <= 8100*Hs^2*K^2*f_c^2",
            "pi_cancels_exactly_after_degree_conversion": True,
            "covers_all_admitted_JONSWAP_gamma_conditionally_on_rotation_envelope": True,
            "covers_arbitrary_heading_phase_and_PSD_cross_axis_coupling": True,
        },
        "coupling_predicate": {
            "squared_product_form": "(Hs*K*f_c)^2 <= threshold/8100",
            "Hs_K_fc_squared_upper": product_squared_threshold,
            "validated_square_root_not_required": True,
        },
        "physical_vessel_rotational_RAO_qualified": False,
        "finite_horizon_body_rate_candidate_producer_ready": True,
        "finite_horizon_good_event_promoted": False,
        "deterministic_left_inclusion_closed": False,
        "P2_pruning_promoted": False,
        "P3_promoted": False,
        "P4_promoted": False,
        "P5_promoted": False,
        "next_obligation": (
            "attach a validated rotational RAO envelope (or a tighter direct body-rate covariance-trace enclosure) for the physical vessel population, evaluate it with this coupling predicate, then compose it with the acceleration covariance candidate in the finite-horizon concentration certificate"
        ),
    }

    # Nonzero exact-rational interface witness.  This demonstrates only that the
    # sufficient set is not the trivial zero-rotation set; it is explicitly not
    # evidence that a particular vessel has this RAO.
    witness = evaluate_tuple(
        out,
        hs_m=float(physical["repository_total_Hs_upper_m"]),
        rotation_gain_rad_per_m=0.125,
        corner_hz=0.03,
        rolloff_power=1.0,
    )
    out["constructed_nonphysical_interface_witness"] = {
        "purpose": "prove the sufficient rotational coupling set is nonempty without claiming hull identification",
        "not_a_measured_or_declared_vessel_RAO": True,
        "evaluation": witness,
    }
    return out


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    if not isinstance(d.get("samples"), int) or isinstance(d.get("samples"), bool) or d.get("samples", 0) <= 0:
        failures.append("invalid finite horizon")
    for key in (
        "validated_arithmetic",
        "outward_rounded",
        "finite_horizon_body_rate_candidate_producer_ready",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "ordinary_libm_transcendental_used_in_pass_decision",
        "physical_vessel_rotational_RAO_qualified",
        "finite_horizon_good_event_promoted",
        "deterministic_left_inclusion_closed",
        "P2_pruning_promoted",
        "P3_promoted",
        "P4_promoted",
        "P5_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    interface = d.get("rotational_response_interface", {})
    if interface.get("universal_rotational_parameter_box_declared") is not False:
        failures.append("rotational theorem silently declared a universal RAO box")
    if interface.get("actual_vessel_or_response_family_must_supply_K_fc_q") is not True:
        failures.append("physical rotational qualification obligation disappeared")
    theorem = d.get("uniform_covariance_theorem", {})
    if theorem.get("pi_cancels_exactly_after_degree_conversion") is not True:
        failures.append("exact degree-conversion simplification disappeared")
    if theorem.get("covers_all_admitted_JONSWAP_gamma_conditionally_on_rotation_envelope") is not True:
        failures.append("conditional JONSWAP coverage disappeared")
    if theorem.get("covers_arbitrary_heading_phase_and_PSD_cross_axis_coupling") is not True:
        failures.append("direction/phase/cross-axis coverage disappeared")
    threshold = d.get("required_body_rate_trace_covariance_upper_deg2_s2")
    if not isinstance(threshold, (int, float)) or isinstance(threshold, bool) or not math.isfinite(float(threshold)) or float(threshold) <= 0.0:
        failures.append("invalid body-rate covariance threshold")
    product = d.get("coupling_predicate", {}).get("Hs_K_fc_squared_upper")
    if not isinstance(product, (int, float)) or isinstance(product, bool) or not math.isfinite(float(product)) or float(product) <= 0.0:
        failures.append("invalid squared product coupling threshold")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, required=True)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--response-domain", type=Path, default=DEFAULT_RESPONSE_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    d = build(args.samples, args.domain, args.response_domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "samples": d["samples"],
        "body_rate_trace_covariance_threshold": d["required_body_rate_trace_covariance_upper_deg2_s2"],
        "Hs_K_fc_squared_upper": d["coupling_predicate"]["Hs_K_fc_squared_upper"],
        "constructed_interface_witness": d["constructed_nonphysical_interface_witness"]["evaluation"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
