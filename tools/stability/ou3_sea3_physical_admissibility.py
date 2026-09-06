#!/usr/bin/env python3
"""Replay-free physical SEA3 parameter admissibility contract.

SEA3 is the compact theorem-domain sea family.  This SEA0 *parameter*
subcertificate refines that already-compact family with a coupled peak-
steepness constraint; it does not construct compactness and it must not replace
SEA3 by independent H/T boxes or another source family.  The separate open
obligation is the validated finite-window phase-continuous SEA3 realization
that drives the shipping front end, tuner, scheduler and Kalman word.

DNVGL-RP-C205 Sec. 3.5.4 defines peak steepness

    S_p = 2*pi*H_s / (g*T_p^2)

and recommends limiting S_p=1/15 for T_p<=8 s and S_p=1/25 for T_p>=15 s,
with linear interpolation between those periods when no better site-specific
source is available.  The recommendation is based on Norwegian Continental
Shelf measurements and is described there as expected to have more general
validity.

For the fixed three-partition SEA3 model we apply that criterion separately to
each active partition (H_r,T_p,r), then retain the exact energy coupling
sum_r H_r^2=H_s^2 and the repository's declared total H_s<=8.5 m envelope.
Applying the single-sea steepness criterion partitionwise is an explicit,
conservative theorem-design choice; it is not attributed to DNV as a published
multimodal partition theorem.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_PHYSICAL_STEEPNESS_ADMISSIBILITY"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def peak_steepness_limit(tp_s: float) -> float:
    """DNV peak-steepness recommendation with linear interpolation."""
    t = float(tp_s)
    if not (math.isfinite(t) and t > 0.0):
        raise ValueError("peak period must be finite and positive")
    if t <= 8.0:
        return 1.0 / 15.0
    if t >= 15.0:
        return 1.0 / 25.0
    a = (t - 8.0) / 7.0
    return (1.0 - a) / 15.0 + a / 25.0


def peak_steepness(hs_m: float, tp_s: float, gravity_mps2: float) -> float:
    h = float(hs_m)
    t = float(tp_s)
    g = float(gravity_mps2)
    if not (math.isfinite(h) and h >= 0.0):
        raise ValueError("significant height must be finite and nonnegative")
    if not (math.isfinite(t) and t > 0.0):
        raise ValueError("peak period must be finite and positive")
    if not (math.isfinite(g) and g > 0.0):
        raise ValueError("gravity must be finite and positive")
    return 2.0 * math.pi * h / (g * t * t)


def significant_height_limit_from_peak_steepness(tp_s: float, gravity_mps2: float) -> float:
    t = float(tp_s)
    g = float(gravity_mps2)
    s = peak_steepness_limit(t)
    return s * g * t * t / (2.0 * math.pi)


def partition_admissible(h_r_m: float, tp_r_s: float, gravity_mps2: float) -> bool:
    return peak_steepness(h_r_m, tp_r_s, gravity_mps2) <= peak_steepness_limit(tp_r_s)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    g = float(domain["startup"]["gravity_mps2"])
    hs_total_upper = float(
        domain["initial_filter_entrance"]["position"]["significant_wave_height_Hs_upper_m"]
    )

    # Useful exact theorem anchors.  They are evaluations of the piecewise law,
    # not replay-derived operating points.
    anchors = {}
    for t in (1.0, 8.0, 11.5, 15.0, 30.0):
        anchors[str(t)] = {
            "T_p_s": t,
            "S_p_upper": peak_steepness_limit(t),
            "H_s_from_steepness_upper_m": up(significant_height_limit_from_peak_steepness(t, g)),
            "H_s_after_repository_total_cap_m": up(
                min(hs_total_upper, significant_height_limit_from_peak_steepness(t, g))
            ),
        }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "SEA0_full_certificate_promoted": False,
        "sea_modes_max": 3,
        "gravity_mps2": g,
        "repository_total_Hs_upper_m": hs_total_upper,
        "peak_steepness_definition": "S_p = 2*pi*H/(g*T_p^2)",
        "peak_steepness_limit": {
            "T_p_le_8_s": 1.0 / 15.0,
            "T_p_ge_15_s": 1.0 / 25.0,
            "T_p_8_to_15_s": "linear interpolation of S_p between 1/15 and 1/25",
        },
        "external_basis": {
            "document": "DNVGL-RP-C205 Environmental conditions and environmental loads, Edition September 2019 amended December 2019",
            "section": "3.5.4 Steepness criterion",
            "published_scope_note": "recommended limiting peak steepness values for short-term irregular sea states in the absence of better reliable sources",
            "multimodal_partition_extension_is_our_conservative_theorem_choice": True,
        },
        "three_partition_contract": {
            "inactive_partition": "H_r = 0",
            "active_partition_constraint": "2*pi*H_r/(g*T_p,r^2) <= S_p,max(T_p,r)",
            "total_energy_coupling": "H_s^2 = sum_r H_r^2",
            "total_Hs_upper_m": hs_total_upper,
            "independent_H_r_and_T_p_rectangular_extrema_forbidden": True,
            "independent_three_partition_H_maxima_forbidden": True,
        },
        "analytical_anchor_values": anchors,
        "SEA3_parameter_domain_compact": True,
        "compact_transition_relation_is_theorem_domain": True,
        "this_subcertificate_refines_but_does_not_rectangularize_SEA3": True,
        "P3_may_not_replace_compact_SEA3_with_independent_bounds": True,
        "finite_window_realization_enclosed": False,
        "left_language_inclusion_closed": False,
        "interpretation": (
            "This refines the compact SEA3 theorem domain with the coupled sea-height/peak-period steepness condition. SEA3 compactness is retained as a theorem-domain property; the open task is a validated phase-continuous finite-window SEA3 realization, not construction of a different compact source box."
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        f.append("physical sea domain became replay fitted")
    if d.get("SEA0_full_certificate_promoted") is not False:
        f.append("parameter subcertificate incorrectly promoted SEA0")
    if int(d.get("sea_modes_max", 0)) != 3:
        f.append("M_max changed")
    for key in (
        "SEA3_parameter_domain_compact",
        "compact_transition_relation_is_theorem_domain",
        "this_subcertificate_refines_but_does_not_rectangularize_SEA3",
        "P3_may_not_replace_compact_SEA3_with_independent_bounds",
    ):
        if d.get(key) is not True:
            f.append(f"compact SEA3 contract lost {key}")
    lim = d.get("peak_steepness_limit", {})
    if not math.isclose(float(lim.get("T_p_le_8_s", math.nan)), 1.0 / 15.0):
        f.append("short-period peak-steepness limit changed")
    if not math.isclose(float(lim.get("T_p_ge_15_s", math.nan)), 1.0 / 25.0):
        f.append("long-period peak-steepness limit changed")
    c = d.get("three_partition_contract", {})
    if c.get("independent_H_r_and_T_p_rectangular_extrema_forbidden") is not True:
        f.append("independent H/T extrema were reintroduced")
    if c.get("independent_three_partition_H_maxima_forbidden") is not True:
        f.append("independent partition-height maxima were reintroduced")
    if d.get("finite_window_realization_enclosed") is not False:
        f.append("steepness contract falsely claims finite-window realization")
    if d.get("left_language_inclusion_closed") is not False:
        f.append("steepness contract falsely closes left inclusion")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "SEA3_compact": d["SEA3_parameter_domain_compact"],
        "Hs_total_upper_m": d["repository_total_Hs_upper_m"],
        "Tp8_Hs_steepness_upper_m": d["analytical_anchor_values"]["8.0"]["H_s_from_steepness_upper_m"],
        "Tp15_Hs_after_total_cap_m": d["analytical_anchor_values"]["15.0"]["H_s_after_repository_total_cap_m"],
        "left_inclusion_closed": d["left_language_inclusion_closed"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
