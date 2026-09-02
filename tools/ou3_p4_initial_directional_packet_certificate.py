#!/usr/bin/env python3
"""Initial finite-angle Joseph directional packet certificate for OU-III P4.

This producer composes three already-rigorous ingredients without scalarizing a
rank-deficient packet:

1. the source vector persistent-excitation lower bound used by P3;
2. the exact 0.80-rad Cayley vector-information retention factor; and
3. the source-correlated goLive Joseph S^-1/R^-1 attenuation from
   :mod:`ou3_p4_joseph_directional_weight`.

If alpha_6 is the P3 R^-1-weighted two-vector/gyro-bias directional information
lower, kappa_R is the exact finite-angle vector information retention, and c_J
is the minimum accelerometer/magnetometer Joseph attenuation, then

    alpha_6,J,sector >= alpha_6 * kappa_R * c_J > 0.

The result is a directional lower bound only.  The same-sample vector map has
exact rank five, so the certificate explicitly refuses a positive scalar
full-state packet margin.  Later P4 work must propagate the covariance and this
directional information through every prefix and accumulate it with S=0 and
prediction transport over the complete H/A word.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as SECTOR
import ou3_p4_effective_word_transport as WORD
import ou3_p4_joseph_directional_weight as JWEIGHT
import ou3_source_reachable_matrix_p3 as P3
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("initial directional packet must not be trajectory fitted")

    sector = SECTOR.build(path)
    word = WORD.build(path)
    weight = JWEIGHT.build(path)
    vector = VECTOR.build()
    failures = [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"word: {x}" for x in WORD.validate(word)]
    failures += [f"Joseph-weight: {x}" for x in JWEIGHT.validate(weight)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    live = domain["normal_live"]
    alpha6 = float(P3.vector_alpha6(live, vector))
    kappa = float(sector["exact_vector_information_retention_factor_lower"])
    if not (math.isfinite(alpha6) and alpha6 > 0.0):
        failures.append("P3 vector alpha6 is not finite positive")
    if not (math.isfinite(kappa) and 0.0 < kappa <= 1.0):
        failures.append("finite-angle information retention is outside (0,1]")

    modes = {}
    for mode in ("H", "A"):
        wm = word["modes"][mode]
        jm = weight["modes"][mode]
        if wm.get("vector_packet_rank_exact") != 5:
            failures.append(f"{mode}: vector packet rank is not five")
        c_acc = float(jm["accelerometer"]["joseph_vs_R_inverse_directional_attenuation_lower"])
        c_mag = float(jm["magnetometer"]["joseph_vs_R_inverse_directional_attenuation_lower"])
        c_vec = down(min(c_acc, c_mag))
        alpha_sector = down(alpha6 * kappa)
        alpha_joseph_sector = down(alpha_sector * c_vec)
        if not alpha_joseph_sector > 0.0:
            failures.append(f"{mode}: finite-angle Joseph directional information lost positivity")
        modes[mode] = {
            "dimension": 18 if mode == "H" else 21,
            "P3_R_inverse_vector_gyro_bias_alpha6_lower": alpha6,
            "finite_angle_vector_information_retention_lower": kappa,
            "initial_accelerometer_Joseph_attenuation_lower": c_acc,
            "initial_magnetometer_Joseph_attenuation_lower": c_mag,
            "initial_vector_Joseph_attenuation_lower": c_vec,
            "finite_angle_R_inverse_directional_alpha6_lower": alpha_sector,
            "initial_finite_angle_Joseph_directional_alpha6_lower": alpha_joseph_sector,
            "vector_packet_rank_exact": 5,
            "instantaneous_full_state_scalar_margin_valid": False,
            "complete_word_prefix_covariances_propagated_here": False,
            "complete_word_directional_accumulation_closed_here": False,
            "P4_PROMOTED": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_INITIAL_FINITE_ANGLE_JOSEPH_DIRECTIONAL_PACKET",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_sector_angle_rad": float(sector["outer_angle_rad"]),
        "same_information_metric_route": True,
        "standalone_eta_norm_penalty_used": False,
        "condition_number_conversion_used": False,
        "per_packet_full_state_scalarization_used": False,
        "modes": modes,
        "P4_INITIAL_DIRECTIONAL_JOSEPH_PACKET_ESTABLISHED": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "propagate source-correlated prefix covariance through prediction/Joseph/reset so the S^-1 attenuation can be re-evaluated at every accepted vector packet; "
            "pull each positive directional form to the common word endpoint, add exact S=0 directional credit, and take a scalar generalized margin only after the complete H/A word has full rank"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "same_information_metric_route",
        "P4_INITIAL_DIRECTIONAL_JOSEPH_PACKET_ESTABLISHED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "standalone_eta_norm_penalty_used",
        "condition_number_conversion_used",
        "per_packet_full_state_scalarization_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED",
        "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("declared_sector_angle_rad", 0.0)) < 0.80:
        f.append("declared P4 sector regressed below 0.80 rad")
    for mode, n in (("H", 18), ("A", 21)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != n:
            f.append(f"{mode}: dimension mismatch")
        if m.get("vector_packet_rank_exact") != 5:
            f.append(f"{mode}: vector packet rank mismatch")
        for key in (
            "P3_R_inverse_vector_gyro_bias_alpha6_lower",
            "finite_angle_vector_information_retention_lower",
            "initial_vector_Joseph_attenuation_lower",
            "finite_angle_R_inverse_directional_alpha6_lower",
            "initial_finite_angle_Joseph_directional_alpha6_lower",
        ):
            x = m.get(key)
            if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or not float(x) > 0.0:
                f.append(f"{mode}.{key} is not finite positive")
        if m.get("instantaneous_full_state_scalar_margin_valid") is not False:
            f.append(f"{mode}: impossible instantaneous scalar margin reintroduced")
        if m.get("complete_word_prefix_covariances_propagated_here") is not False:
            f.append(f"{mode}: prematurely claims prefix covariance propagation")
        if m.get("complete_word_directional_accumulation_closed_here") is not False:
            f.append(f"{mode}: prematurely claims complete directional accumulation")
        if m.get("P4_PROMOTED") is not False:
            f.append(f"{mode}: P4 prematurely promoted")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "sector_rad": d["declared_sector_angle_rad"],
        "H_alpha_J_sector": d["modes"]["H"]["initial_finite_angle_Joseph_directional_alpha6_lower"],
        "A_alpha_J_sector": d["modes"]["A"]["initial_finite_angle_Joseph_directional_alpha6_lower"],
        "P4": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
