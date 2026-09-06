#!/usr/bin/env python3
"""Analytical SEA3 response-regularity moment bound.

This is a supporting theorem for the eventual hard finite-window SEA0 provider,
not a source generator and not a P3 promotion route.

For one directional surface component with angular frequency omega and the
continuum CoG translation RAO envelope

    ||h(omega,theta)|| <= G min(1,(omega_c/omega)^p),  p >= 2,

the acceleration response is ``a_hat=-omega^2 h eta_hat`` and therefore the
jerk response satisfies, for all frequencies,

    omega^6 ||h||^2 <= G^2 omega_c^4 omega^2.

Thus the total translational jerk second moment is bounded by

    tr M_jerk <= G^2 omega_c^4 m2_eta,omega.

The remaining surface m2 is bounded analytically over all three JONSWAP/PM
partitions without inventing independent H/T boxes.  For the normalized
JONSWAP shape, ``1 <= gamma^r <= gamma <= 7`` gives

    I2(gamma)/I0(gamma) <= 7 I2(PM)/I0(PM)
                         = 7 sqrt(5*pi/4).

The coupled physical steepness condition gives

    H_r/Tp_r^2 <= g/(30*pi)

because the maximum admitted peak steepness is 1/15.  Together with
``sum H_r^2=Hs^2``, ``M<=3`` and ``Hs<=8.5 m``, Cauchy--Schwarz gives

    sum H_r^2/Tp_r^2
      <= g/(30*pi) sum H_r
      <= g/(30*pi) sqrt(3) Hs_max.

This produces a finite theorem-wide m2/jerk-moment bound without requiring an
independent numerical Tp rectangle.  It remains a *moment* statement: it may
support a validated shaping-state/IQC construction, but it cannot by itself
materialize a deterministic 601-sample SEA3 word.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_directional_response_family as RESPONSE
import ou3_sea3_physical_admissibility as PHYSICAL

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_ANALYTICAL_RESPONSE_REGULARITY_MOMENT_V1"
PM_A = 5.0 / 4.0
GAMMA_MAX = 7.0
M_MAX = 3


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    physical = PHYSICAL.build(Path(domain_path).resolve())
    pf = PHYSICAL.validate(physical)
    response = RESPONSE.directional_response_enclosure(REPO)
    rf = RESPONSE.validate(response)
    if pf or rf:
        raise RuntimeError(f"SEA3 regularity prerequisites failed: physical={pf}, response={rf}")

    g = float(physical["gravity_mps2"])
    hs_max = float(physical["repository_total_Hs_upper_m"])
    steep_max = float(physical["peak_steepness_limit"]["T_p_le_8_s"])
    if not math.isclose(steep_max, 1.0 / 15.0, rel_tol=0.0, abs_tol=1e-15):
        raise RuntimeError("maximum SEA3 peak-steepness constant changed")

    box = response["rao_envelope_parameter_box"]
    Gmax = float(box["peak_translation_gain"][1])
    fcmax = float(box["rolloff_corner_hz"][1])
    pmin = float(box["high_frequency_rolloff_power_min"])
    if pmin < 2.0:
        raise RuntimeError("jerk regularity inequality requires p>=2")

    # I2_PM/I0_PM = sqrt(a*pi), a=5/4.  Since 1<=gamma^r<=7 pointwise,
    # I2_gamma/I0_gamma <= 7 I2_PM/I0_PM.
    shape_m2_over_m0_upper = up(GAMMA_MAX * math.sqrt(PM_A * math.pi))

    # H/Tp^2 <= g*Smax/(2*pi) with Smax=1/15.  Cauchy--Schwarz is applied to
    # the coupled partition vector, never to three independent H maxima.
    h_over_tp2_upper = up(g * steep_max / (2.0 * math.pi))
    sum_h_upper = up(math.sqrt(M_MAX) * hs_max)
    sum_h2_over_tp2_upper = up(h_over_tp2_upper * sum_h_upper)

    # m0_r=H_r^2/16 and fp=1/Tp. Convert f^2 moment to omega^2 by (2*pi)^2.
    two_pi = up(2.0 * math.pi)
    surface_m2_omega_upper = up(
        (two_pi * two_pi / 16.0)
        * shape_m2_over_m0_upper
        * sum_h2_over_tp2_upper
    )

    omega_c = up(two_pi * fcmax)
    jerk_trace_upper = up(
        Gmax * Gmax
        * omega_c**4
        * surface_m2_omega_upper
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "spectral_moment_only_source_used": False,
        "SEA3_parameter_domain_compact": True,
        "independent_H_T_rectangle_used": False,
        "independent_partition_H_maxima_used": False,
        "finite_RAO_grid_used": False,
        "same_continuum_RAO_family_consumed": True,
        "analytical_inequalities": {
            "pointwise_jerk_transfer": (
                "omega^6*||h||^2 <= G^2*omega_c^4*omega^2 for p>=2"
            ),
            "jonswap_shape_m2_over_m0_upper": (
                "I2(gamma)/I0(gamma) <= 7*sqrt(5*pi/4)"
            ),
            "partition_steepness": "H_r/Tp_r^2 <= g/(30*pi)",
            "coupled_partition_sum": (
                "sum H_r^2/Tp_r^2 <= g/(30*pi)*sqrt(3)*Hs_max"
            ),
        },
        "constants": {
            "gamma_upper": GAMMA_MAX,
            "modes_max": M_MAX,
            "Hs_upper_m": hs_max,
            "gravity_mps2": g,
            "peak_steepness_upper": steep_max,
            "RAO_gain_upper": Gmax,
            "RAO_corner_hz_upper": fcmax,
            "RAO_rolloff_power_lower": pmin,
        },
        "bounds": {
            "JONSWAP_dimensionless_I2_over_I0_upper": shape_m2_over_m0_upper,
            "H_over_Tp2_upper_mps2": h_over_tp2_upper,
            "sum_H_upper_m": sum_h_upper,
            "sum_H2_over_Tp2_upper_m2ps2": sum_h2_over_tp2_upper,
            "surface_m2_omega_upper_m2ps2": surface_m2_omega_upper,
            "translation_jerk_trace_upper_m2ps6": jerk_trace_upper,
        },
        "hard_finite_window_source_materialized": False,
        "hard_shaping_state_or_excitation_bound_closed_here": False,
        "may_support_future_hard_window_IQC": True,
        "may_substitute_for_hard_window_IQC": False,
        "P3_promoted": False,
        "next_obligation": (
            "combine this theorem-wide response regularity with the actual compact x^s/R_lambda hard "
            "state or excitation constraint; a second moment alone cannot generate a deterministic SEA3 word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "SEA3_parameter_domain_compact",
        "same_continuum_RAO_family_consumed",
        "may_support_future_hard_window_IQC",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_generator",
        "trajectory_replay_used",
        "spectral_moment_only_source_used",
        "independent_H_T_rectangle_used",
        "independent_partition_H_maxima_used",
        "finite_RAO_grid_used",
        "hard_finite_window_source_materialized",
        "hard_shaping_state_or_excitation_bound_closed_here",
        "may_substitute_for_hard_window_IQC",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    b = d.get("bounds", {})
    for key in (
        "JONSWAP_dimensionless_I2_over_I0_upper",
        "H_over_Tp2_upper_mps2",
        "sum_H_upper_m",
        "sum_H2_over_Tp2_upper_m2ps2",
        "surface_m2_omega_upper_m2ps2",
        "translation_jerk_trace_upper_m2ps6",
    ):
        x = float(b.get(key, math.nan))
        if not (math.isfinite(x) and x > 0.0):
            f.append(f"invalid analytical regularity bound {key}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "bounds": d["bounds"],
        "hard_window_materialized": d["hard_finite_window_source_materialized"],
        "P3": d["P3_promoted"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
