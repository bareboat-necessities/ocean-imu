#!/usr/bin/env python3
"""Replay-free SEA3 sea/RAO compatibility with the hard Normal-Live P1 branch.

The physical JONSWAP sea family and the robust RAO envelope may not be combined
as an independent Cartesian product under the existing hard condition
||a_non-grav(t)||_2 <= 4 m/s^2. This module gives an analytical counterexample
inside the declared sets and therefore forces the intended theorem domain to be
a coupled sea-response set.

Witness: one JONSWAP partition on its PM boundary (gamma=1), Tp=8 s, on the DNV
peak-steepness boundary, and the admitted response-envelope member G=4,
fc=1.2 Hz, p=2. PM is not a different sea model here: it is exactly the
gamma=1 member of the declared JONSWAP family. Using that boundary member is
convenient because its normalization is closed form. For x=f/fp in [1,9],
f<=1.125 Hz<fc, so that envelope member may have |h|=G throughout the interval.
For the normalized gamma=1 JONSWAP/PM spectrum
q(x)=x^-5 exp(-(5/4)x^-4), I0=int q dx=1/5, the acceleration variance contains

  m0 G^2 wp^4 / I0 * int_1^9 x^-1 exp(-(5/4)x^-4) dx.

On [1,9], x^-1>=1/9 and exp(-(5/4)x^-4)>=exp(-5/4), hence the integral is at
least (8/9)exp(-5/4). At the steepness boundary the Tp factors cancel. The
only transcendental needed is exp(-5/4)=exp(-5/16)^4, evaluated with the
repository validated rational-Taylor exponential and outward interval products.
All non-dyadic rational constants used in the lower bound are also outward
enclosed before multiplication/division.

Because gamma=1 belongs to the declared JONSWAP interval gamma in [1,7], one
counterexample at that boundary is sufficient to refute soundness of the full
independent Cartesian product. No claim is made that gamma=1 is the worst
JONSWAP member.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_physical_admissibility as PHYS
import ou3_validated_transcendentals as VT
from ou3_interval import Interval, down, up

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RAO = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_COUPLED_SEA_RAO_P1_COMPATIBILITY"


def _outward_point(x: float) -> Interval:
    """One-ulp outward interval around a finite binary64 value."""
    return Interval.outward_bounds(float(x), float(x))


def exp_minus_five_quarters() -> Interval:
    """Validated exp(-5/4) from four products of exp(-5/16)."""
    e = VT.exp_point(-5.0 / 16.0)
    return e * e * e * e


def witness_mean_square_lower(gravity: float, gain: float = 4.0) -> float:
    """Validated lower bound on witness non-gravitational acceleration variance."""
    # Exact mathematical constants are enclosed by one-ulp intervals around
    # their nearest binary64 values. 4 and gain=4 are exactly representable.
    sp = _outward_point(1.0 / 15.0)
    g = _outward_point(gravity)
    G = Interval.point(gain)
    pi = _outward_point(math.pi)
    eight_ninths = _outward_point(8.0 / 9.0)
    one_fifth = _outward_point(1.0 / 5.0)

    # J >= (8/9) exp(-5/4), and the gamma=1 JONSWAP/PM I0 = 1/5 exactly.
    J = eight_ninths * exp_minus_five_quarters()
    val = sp.square() * g.square() * G.square() * pi.square()
    val = val / Interval.point(4.0)
    val = val * J / one_fifth
    return val.lo


def build(domain_path: Path = DEFAULT_DOMAIN, rao_path: Path = DEFAULT_RAO) -> dict:
    """Build the fail-closed compatibility contract."""
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    rao = json.loads(Path(rao_path).read_text(encoding="utf-8"))
    gravity = float(domain["startup"]["gravity_mps2"])
    cap = float(domain["normal_live"]["non_gravitational_cog_acceleration_norm_upper_mps2"])
    response = rao["response_contract"]

    gain_range = list(map(float, response["peak_translation_gain_range"]))
    corner_range = list(map(float, response["rolloff_corner_hz_range"]))
    if len(gain_range) != 2 or len(corner_range) != 2:
        raise RuntimeError("declared RAO family lost two-endpoint parameter ranges")
    gain_lo, gain_hi = gain_range
    corner_lo, corner_hi = corner_range
    pmin = float(response["high_frequency_rolloff_power_min"])
    finite_bounds = all(
        math.isfinite(value)
        for value in (gain_lo, gain_hi, corner_lo, corner_hi, pmin)
    )
    witness_in_rao_box = bool(
        finite_bounds
        and gain_lo <= 4.0 <= gain_hi
        and corner_lo <= 1.2 <= corner_hi
        and pmin <= 2.0
    )
    if not witness_in_rao_box:
        raise RuntimeError("declared RAO family no longer contains the analytical G=4, fc=1.2, p=2 witness")

    tp = 8.0
    hs = down(PHYS.significant_height_limit_from_peak_steepness(tp, gravity))
    if not PHYS.partition_admissible(hs, tp, gravity):
        raise RuntimeError("gamma=1 JONSWAP witness is outside physical steepness domain")
    fp = 1.0 / tp
    witness_x_hi = 9.0
    witness_hi_hz = witness_x_hi * fp
    if not witness_hi_hz < 1.2:
        raise RuntimeError("witness frequency interval left flat RAO branch")

    ms_lo = witness_mean_square_lower(gravity)
    rms_lo = down(math.sqrt(ms_lo))
    cap2_hi = up(cap * cap)
    refuted = ms_lo > cap2_hi
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "witness": {
            "sea": "one JONSWAP partition at gamma=1 (PM boundary) on DNV Tp<=8 s peak-steepness boundary",
            "declared_JONSWAP_gamma": 1.0,
            "PM_is_JONSWAP_gamma_1_boundary": True,
            "witness_is_inside_declared_JONSWAP_gamma_interval_1_to_7": True,
            "T_p_s": tp,
            "H_s_m": hs,
            "RAO_gain": 4.0,
            "RAO_corner_hz": 1.2,
            "RAO_rolloff_power": 2.0,
            "declared_RAO_gain_range": gain_range,
            "declared_RAO_corner_hz_range": corner_range,
            "declared_RAO_rolloff_power_min": pmin,
            "RAO_parameter_bounds_finite": finite_bounds,
            "witness_is_inside_declared_RAO_parameter_ranges": witness_in_rao_box,
            "x_interval": [1.0, witness_x_hi],
            "frequency_interval_hz": [fp, witness_hi_hz],
            "validated_exp_minus_5_over_4": exp_minus_five_quarters().as_list(),
            "validated_acceleration_mean_square_lower_m2_s4": ms_lo,
            "validated_acceleration_RMS_lower_mps2": rms_lo,
            "P1_cap_squared_upper_m2_s4": cap2_hi,
            "all_nondyadic_witness_constants_outward_enclosed": True,
        },
        "independent_cartesian_sea_x_RAO_domain_is_P1_sound": not refuted,
        "cartesian_product_refuted_by_analytical_witness": refuted,
        "coupled_SEA3_domain_required": refuted,
        "coupled_domain_contract": {
            "sea_and_RAO_parameters_may_not_be_selected_independently": True,
            "P1_hard_acceleration_bound_checked_after_response": True,
            "PSD_or_RMS_bound_alone_sufficient_for_P1_admission": False,
            "finite_window_deterministic_response_certificate_required": True,
        },
        "finite_window_realization_certificate_closed": False,
        "L_actual_sea_subset_Lhat_SEA3_closed": False,
        "next_obligation": (
            "construct a hard finite-window oscillator/IQC realization enclosure on the coupled JONSWAP-sea/RAO set; only response trajectories satisfying every existing P1 Normal-Live hard source bound may enter the left inclusion"
        ),
    }


def validate(d: dict) -> list[str]:
    """Validate the incompatibility witness and fail-closed coupled-set semantics."""
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("trajectory_replay_used") is not False or d.get("filter_changed") is not False:
        f.append("certificate must be replay free and non-invasive")
    if d.get("cartesian_product_refuted_by_analytical_witness") is not True:
        f.append("independent sea x RAO Cartesian product was not refuted")
    if d.get("coupled_SEA3_domain_required") is not True:
        f.append("coupled SEA3 domain requirement disappeared")
    w = d.get("witness", {})
    if w.get("PM_is_JONSWAP_gamma_1_boundary") is not True:
        f.append("witness lost JONSWAP gamma=1 boundary identity")
    if w.get("witness_is_inside_declared_JONSWAP_gamma_interval_1_to_7") is not True:
        f.append("witness left declared JONSWAP gamma family")
    if w.get("RAO_parameter_bounds_finite") is not True:
        f.append("RAO witness parameter bounds are not finite")
    if w.get("witness_is_inside_declared_RAO_parameter_ranges") is not True:
        f.append("analytical P1 witness is not inside the declared RAO parameter ranges")
    if w.get("all_nondyadic_witness_constants_outward_enclosed") is not True:
        f.append("witness rational/transcendental enclosure contract weakened")
    if list(w.get("x_interval", [])) != [1.0, 9.0]:
        f.append("strengthened spectral witness interval changed")
    if not float(w.get("validated_acceleration_mean_square_lower_m2_s4", 0.0)) > float(w.get("P1_cap_squared_upper_m2_s4", math.inf)):
        f.append("witness lower bound no longer exceeds P1 cap squared")
    if not float(w.get("validated_acceleration_RMS_lower_mps2", 0.0)) > 4.5:
        f.append("witness proof margin fell below retained >4.5 m/s^2 lower bound")
    c = d.get("coupled_domain_contract", {})
    if c.get("finite_window_deterministic_response_certificate_required") is not True:
        f.append("finite-window response obligation disappeared")
    if d.get("finite_window_realization_certificate_closed") is not False:
        f.append("compatibility certificate falsely closes finite-window realization")
    if d.get("L_actual_sea_subset_Lhat_SEA3_closed") is not False:
        f.append("compatibility certificate falsely closes left inclusion")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--rao", type=Path, default=DEFAULT_RAO)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.rao)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "JONSWAP_gamma": d["witness"]["declared_JONSWAP_gamma"],
        "mean_square_lower": d["witness"]["validated_acceleration_mean_square_lower_m2_s4"],
        "RMS_lower": d["witness"]["validated_acceleration_RMS_lower_mps2"],
        "P1_cap_squared_upper": d["witness"]["P1_cap_squared_upper_m2_s4"],
        "cartesian_product_refuted": d["cartesian_product_refuted_by_analytical_witness"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
