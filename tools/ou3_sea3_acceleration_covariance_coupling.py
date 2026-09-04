#!/usr/bin/env python3
"""Validated SEA3 sea/RAO coupling condition for acceleration P1 good events.

For the certified translational response family

    ||h(f,theta)|| <= G min(1, (fc/f)^p), p >= 2,

the acceleration weighting obeys, for every f and heading,

    (2 pi f)^4 ||h(f,theta)||^2 <= G^2 (2 pi fc)^4.

Therefore any normalized sea spectrum with total significant height Hs and
m0 = Hs^2/16 has post-response acceleration covariance trace bounded by

    tr Cov[a] <= (Hs^2/16) G^2 (2 pi fc)^4
              = pi^4 Hs^2 G^2 fc^4.

This is spectrum-shape independent. In particular it covers every admitted
JONSWAP gamma once the tuple (Hs,G,fc) satisfies the coupling inequality. It
also makes explicit why the independently cartesianized SEA3 sea x RAO box is
too broad: high Hs, high G, and high fc cannot all be selected together.

The bound is compared to the finite-horizon Gaussian concentration threshold
from ``ou3_sea3_finite_horizon_concentration``. No square root is used in the
PASS decision; the covariance inequality is checked directly with outward
interval arithmetic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from ou3_interval import Interval
import ou3_sea3_directional_p2_ha_feasibility as SEA3
import ou3_sea3_finite_horizon_concentration as CONC
import ou3_sea3_physical_admissibility as PHYS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RESPONSE_DOMAIN = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_ACCELERATION_COVARIANCE_COUPLING_V1"


def _outward_point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def acceleration_trace_covariance_upper_m2_s4(
    hs_m: float,
    gain: float,
    corner_hz: float,
) -> float:
    """Validated upper bound pi^4 Hs^2 G^2 fc^4."""
    values = [float(hs_m), float(gain), float(corner_hz)]
    if not all(math.isfinite(x) and x >= 0.0 for x in values):
        raise ValueError("Hs, gain, and corner must be finite and nonnegative")
    hs = _outward_point(values[0])
    g = _outward_point(values[1])
    fc = _outward_point(values[2])
    pi = _outward_point(math.pi)
    bound = pi.square() * pi.square() * hs.square() * g.square() * fc.square() * fc.square()
    return bound.hi


def evaluate_tuple(certificate: dict, *, hs_m: float, gain: float, corner_hz: float) -> dict:
    """Check one coupled sea/RAO tuple against the finite-horizon threshold."""
    failures = validate(certificate)
    if failures:
        raise ValueError(f"invalid coupling certificate: {failures}")
    try:
        covariance_upper = acceleration_trace_covariance_upper_m2_s4(hs_m, gain, corner_hz)
    except (TypeError, ValueError) as exc:
        return {
            "candidate_pass": False,
            "decision": "FAIL_CANDIDATE",
            "validation_failures": [str(exc)],
            "finite_horizon_good_event_promoted": False,
            "deterministic_left_inclusion_promoted": False,
        }

    threshold = float(certificate["required_acceleration_trace_covariance_upper_m2_s4"])
    passed = covariance_upper <= threshold
    return {
        "candidate_pass": passed,
        "decision": "PASS_CANDIDATE" if passed else "FAIL_CANDIDATE",
        "H_s_m": float(hs_m),
        "RAO_gain": float(gain),
        "RAO_corner_hz": float(corner_hz),
        "validated_acceleration_trace_covariance_upper_m2_s4": covariance_upper,
        "required_acceleration_trace_covariance_upper_m2_s4": threshold,
        "margin_m2_s4": threshold - covariance_upper,
        "validation_failures": [] if passed else [
            "coupled Hs/gain/corner tuple exceeds finite-horizon acceleration covariance threshold"
        ],
        "finite_horizon_good_event_promoted": False,
        "deterministic_left_inclusion_promoted": False,
    }


def build(
    samples: int,
    domain_path: Path = DEFAULT_DOMAIN,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
    repo: Path = REPO,
) -> dict:
    concentration = CONC.build(samples, domain_path, response_domain_path, repo)
    concentration_failures = CONC.validate(concentration)
    if concentration_failures:
        raise RuntimeError(f"invalid finite-horizon concentration certificate: {concentration_failures}")

    response = SEA3.directional_response_enclosure(
        Path(repo).resolve(), Path(response_domain_path).resolve()
    )
    box = response["rao_envelope_parameter_box"]
    gain_lo, gain_hi = map(float, box["peak_translation_gain"])
    corner_lo, corner_hi = map(float, box["rolloff_corner_hz"])
    pmin = float(box["high_frequency_rolloff_power_min"])
    if pmin < 2.0:
        raise RuntimeError("acceleration covariance coupling requires RAO rolloff p >= 2")

    physical = PHYS.build(Path(domain_path))
    physical_failures = PHYS.validate(physical)
    if physical_failures:
        raise RuntimeError(f"invalid SEA3 physical admissibility contract: {physical_failures}")
    hs_upper = float(physical["repository_total_Hs_upper_m"])

    threshold = float(
        concentration["acceleration"]["required_trace_covariance_upper_m2_s4"]
    )
    out: dict[str, Any] = {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "samples": samples,
        "trajectory_replay_used": False,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "ordinary_libm_transcendental_used_in_pass_decision": False,
        "RAO_rolloff_power_min": pmin,
        "RAO_rolloff_p_at_least_2_required": True,
        "spectral_shape_quadrature_required": False,
        "JONSWAP_gamma_specific_bound_required": False,
        "temporal_independence_required": False,
        "cross_axis_independence_required": False,
        "response_parameter_box_sha256": concentration["response_parameter_box_sha256"],
        "declared_RAO_parameter_box": box,
        "repository_total_Hs_upper_m": hs_upper,
        "required_acceleration_trace_covariance_upper_m2_s4": threshold,
        "uniform_covariance_theorem": {
            "m0_identity": "m0 = Hs^2 / 16",
            "response_weight_bound": "(2*pi*f)^4*||h(f,theta)||^2 <= G^2*(2*pi*fc)^4",
            "trace_covariance_bound": "tr Cov[a] <= pi^4*Hs^2*G^2*fc^4",
            "covers_all_admitted_JONSWAP_gamma_conditionally_on_tuple": True,
            "covers_arbitrary_heading_and_complex_phase_in_response_envelope": True,
        },
        "independent_cartesian_extreme": {
            "H_s_m": hs_upper,
            "RAO_gain": gain_hi,
            "RAO_corner_hz": corner_hi,
        },
        "high_sea_low_corner_witness": {
            "H_s_m": hs_upper,
            "RAO_gain": gain_hi,
            "RAO_corner_hz": corner_lo,
        },
        "coupling_predicate_defined": True,
        "coupled_domain_nonempty_witness_required": True,
        "physical_vessel_pairing_qualified": False,
        "finite_horizon_good_event_promoted": False,
        "deterministic_left_inclusion_closed": False,
        "P2_pruning_promoted": False,
        "P3_promoted": False,
        "P4_promoted": False,
        "P5_promoted": False,
        "next_obligation": (
            "qualify the actual vessel RAO/sea pairing against this coupling predicate (or a tighter response-weighted spectral bound) and then compose the passing covariance result with the finite-horizon good-event theorem"
        ),
    }
    out["independent_cartesian_extreme"]["evaluation"] = evaluate_tuple(
        out,
        hs_m=hs_upper,
        gain=gain_hi,
        corner_hz=corner_hi,
    )
    out["high_sea_low_corner_witness"]["evaluation"] = evaluate_tuple(
        out,
        hs_m=hs_upper,
        gain=gain_hi,
        corner_hz=corner_lo,
    )
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
        "RAO_rolloff_p_at_least_2_required",
        "coupling_predicate_defined",
        "coupled_domain_nonempty_witness_required",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "ordinary_libm_transcendental_used_in_pass_decision",
        "spectral_shape_quadrature_required",
        "JONSWAP_gamma_specific_bound_required",
        "temporal_independence_required",
        "cross_axis_independence_required",
        "physical_vessel_pairing_qualified",
        "finite_horizon_good_event_promoted",
        "deterministic_left_inclusion_closed",
        "P2_pruning_promoted",
        "P3_promoted",
        "P4_promoted",
        "P5_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if float(d.get("RAO_rolloff_power_min", 0.0)) < 2.0:
        failures.append("RAO rolloff p < 2 invalidates uniform acceleration weighting")
    theorem = d.get("uniform_covariance_theorem", {})
    if theorem.get("covers_all_admitted_JONSWAP_gamma_conditionally_on_tuple") is not True:
        failures.append("conditional JONSWAP gamma coverage disappeared")
    if theorem.get("covers_arbitrary_heading_and_complex_phase_in_response_envelope") is not True:
        failures.append("direction/phase coverage disappeared")
    try:
        threshold = float(d["required_acceleration_trace_covariance_upper_m2_s4"])
    except (KeyError, TypeError, ValueError):
        threshold = math.nan
    if not (math.isfinite(threshold) and threshold > 0.0):
        failures.append("invalid acceleration covariance threshold")
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
        "threshold": d["required_acceleration_trace_covariance_upper_m2_s4"],
        "cartesian_extreme": d["independent_cartesian_extreme"]["evaluation"],
        "high_sea_low_corner": d["high_sea_low_corner_witness"]["evaluation"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
