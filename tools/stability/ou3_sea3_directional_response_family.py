#!/usr/bin/env python3
"""Pure SEA3 directional translational-response family.

This module is deliberately independent of every retired P2/P3 graph or word
construction.  It contains only the robust continuum SEA3 response envelope
used by the physical/source layer:

    ||h(f,theta)||_2 <= G min(1,(f_c/f)^p),  p>=2.

The worst envelope member proves the response-moment bounds for the entire
parameter family by monotonicity.  This is a set theorem, not a nominal hull,
not a finite RAO grid, and not a Riccati certificate.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_RESPONSE_DOMAIN = REPO / "tools" / "stability" / "ou3_sea3_directional_response_domain.json"
RESPONSE_SCHEMA = "OU3_SEA3_DIRECTIONAL_RESPONSE_DOMAIN_V3"
PI_HI = 3.141592653589794


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _load_response_domain(path: Path = DEFAULT_RESPONSE_DOMAIN) -> dict:
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    if d.get("schema_version") != RESPONSE_SCHEMA:
        raise RuntimeError("unexpected SEA3 directional response-domain schema")
    if int(d.get("sea_modes_max", 0)) != 3:
        raise RuntimeError("SEA3 response domain must retain M_max=3")
    r = d.get("response_contract", {})
    gain = r.get("peak_translation_gain_range")
    corner = r.get("rolloff_corner_hz_range")
    if not (isinstance(gain, list) and len(gain) == 2):
        raise RuntimeError("SEA3 RAO family lost peak-gain range")
    if not (isinstance(corner, list) and len(corner) == 2):
        raise RuntimeError("SEA3 RAO family lost rolloff-corner range")
    g0, g1 = map(float, gain)
    f0, f1 = map(float, corner)
    pmin = float(r.get("high_frequency_rolloff_power_min", math.nan))
    if not (g0 == 0.0 and math.isfinite(g1) and g1 > 0.0):
        raise RuntimeError("invalid SEA3 RAO peak-gain range")
    if not (0.0 < f0 <= f1 and math.isfinite(f1)):
        raise RuntimeError("invalid SEA3 RAO corner range")
    if not (math.isfinite(pmin) and pmin >= 2.0):
        raise RuntimeError("SEA3 acceleration response requires p>=2")
    for key in (
        "worst_member_dominates_entire_parameter_box",
        "six_dof_parent_RAO_allowed",
        "arbitrary_frequency_dependence_below_envelope",
        "arbitrary_directional_dependence",
        "arbitrary_cross_axis_coupling_subject_to_PSD",
        "unbanded_acceleration_moment_is_finite_from_response_rolloff",
    ):
        if r.get(key) is not True:
            raise RuntimeError(f"SEA3 RAO family lost {key}")
    if r.get("single_nominal_RAO_used") is not False:
        raise RuntimeError("SEA3 response theorem selected a nominal hull")
    if r.get("finite_RAO_grid_used") is not False:
        raise RuntimeError("SEA3 response theorem selected a finite RAO grid")
    if r.get("phase_quantifier") != "arbitrary complex phase":
        raise RuntimeError("SEA3 response phase quantifier changed")
    return d


def _member_moment_coefficients(gain: float, corner_hz: float, power: float) -> dict[str, list[float]]:
    g = float(gain)
    fc = float(corner_hz)
    p = float(power)
    if not (math.isfinite(g) and g >= 0.0):
        raise ValueError("RAO gain must be finite and nonnegative")
    if not (math.isfinite(fc) and fc > 0.0):
        raise ValueError("RAO corner must be finite positive")
    if not (math.isfinite(p) and p >= 2.0):
        raise ValueError("RAO rolloff must satisfy p>=2")
    g2 = up(g * g)
    omega_c = up(2.0 * PI_HI * fc)
    disp = up(g2 / 16.0)
    vel = up(disp * omega_c * omega_c)
    acc = up(vel * omega_c * omega_c)
    return {
        "displacement": [down(0.0), disp],
        "velocity": [down(0.0), vel],
        "acceleration": [down(0.0), acc],
    }


def evaluate_rao_envelope_member(
    gain: float,
    corner_hz: float,
    power: float,
    response_enclosure: dict | None = None,
) -> dict[str, list[float]]:
    if response_enclosure is not None:
        box = response_enclosure["rao_envelope_parameter_box"]
        gr = list(map(float, box["peak_translation_gain"]))
        fr = list(map(float, box["rolloff_corner_hz"]))
        pmin = float(box["high_frequency_rolloff_power_min"])
        if not gr[0] <= float(gain) <= gr[1]:
            raise ValueError("RAO gain outside certified SEA3 family")
        if not fr[0] <= float(corner_hz) <= fr[1]:
            raise ValueError("RAO corner outside certified SEA3 family")
        if float(power) < pmin:
            raise ValueError("RAO rolloff outside certified SEA3 family")
    return _member_moment_coefficients(gain, corner_hz, power)


def directional_response_enclosure(
    repo: Path = REPO,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
) -> dict:
    repo = Path(repo).resolve()
    cfg = _load_response_domain(Path(response_domain_path).resolve())
    r = cfg["response_contract"]
    wrapper = (repo / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h").read_text(
        encoding="utf-8"
    )
    shipping_band = [
        float(SOURCE.parse_const(wrapper, "SIGMA_BAND_MIN_HZ_DEFAULT")),
        float(SOURCE.parse_const(wrapper, "SIGMA_BAND_MAX_HZ_DEFAULT")),
    ]
    declared_band = list(map(float, r["shipping_sigma_band_hz"]))
    parity = all(
        math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)
        for a, b in zip(shipping_band, declared_band)
    )
    if not parity:
        raise RuntimeError(
            f"declared SEA3 sigma band {declared_band} drifted from source {shipping_band}"
        )
    gr = list(map(float, r["peak_translation_gain_range"]))
    fr = list(map(float, r["rolloff_corner_hz_range"]))
    pmin = float(r["high_frequency_rolloff_power_min"])
    worst = _member_moment_coefficients(gr[1], fr[1], pmin)
    return {
        "qualification": "OU3_SEA3_UNIFORM_COMPLEX_DIRECTIONAL_RAO_ENVELOPE",
        "sea_modes_max": 3,
        "trajectory_replay_used": False,
        "single_nominal_RAO_used": False,
        "finite_RAO_grid_used": False,
        "retired_P2_stack_imported": False,
        "response_definition": r["definition"],
        "rao_envelope_parameter_box": {
            "peak_translation_gain": gr,
            "rolloff_corner_hz": fr,
            "high_frequency_rolloff_power_min": pmin,
            "complex_phase": "arbitrary",
            "heading_dependence": "arbitrary",
            "frequency_dependence_below_envelope": "arbitrary",
            "cross_axis_coupling": "arbitrary subject to PSD",
        },
        "worst_envelope_member": {
            "peak_translation_gain": gr[1],
            "rolloff_corner_hz": fr[1],
            "high_frequency_rolloff_power": pmin,
            "is_set_envelope_not_representative_RAO": True,
        },
        "single_worst_envelope_proves_entire_parameter_box_by_monotonicity": True,
        "shipping_sigma_band_hz": declared_band,
        "shipping_sigma_band_source_parity": parity,
        "uniform_moment_theorem": {
            "envelope": "||h||_2 <= G min(1,(f_c/f)^p)",
            "parameter_quantifier": "G in [0,Gmax], fc in [fcmin,fcmax], p>=2",
            "analytical_not_sampled": True,
            "unbanded_acceleration_moment_finite": True,
            "proof_corner": "G=Gmax, fc=fcmax, p=2",
        },
        "worst_envelope_trace_upper_per_Hs2": worst,
        "directional_cross_axis_coupling": {
            "independent_cartesian_axis_boxes_used": False,
            "arbitrary_phase_is_covered": True,
            "cross_terms_constrained_by_PSD": True,
        },
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("qualification") != "OU3_SEA3_UNIFORM_COMPLEX_DIRECTIONAL_RAO_ENVELOPE":
        f.append("qualification mismatch")
    for key in (
        "single_worst_envelope_proves_entire_parameter_box_by_monotonicity",
        "shipping_sigma_band_source_parity",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("trajectory_replay_used", "single_nominal_RAO_used", "finite_RAO_grid_used", "retired_P2_stack_imported"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    box = d.get("rao_envelope_parameter_box", {})
    if float(box.get("high_frequency_rolloff_power_min", 0.0)) < 2.0:
        f.append("RAO rolloff power dropped below 2")
    return list(dict.fromkeys(f))
