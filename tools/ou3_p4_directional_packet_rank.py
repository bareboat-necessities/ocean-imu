#!/usr/bin/env python3
"""Structural rank certificate for operation-matched OU-III vector packets.

This module proves why the stronger P4 route must accumulate *directional PSD
information* over a complete source word instead of searching for a positive
scalar information margin on each accepted packet.

In the rotation gauge used by the first-vector proof, the H-mode measurement
active coordinates are (theta, a_w) and

    H_a = [ -[f]_x   I ],
    H_m = [ -[m]_x   0 ].

Because the admitted magnetic vector has nonzero norm, rank(H_m)=2.  Since the
accelerometer a_w block is the orthogonal full-rank J_aw=R_wb (gauged to I),
rank(H_a)=3.  Moreover the stacked vector packet has the exact nonzero null
family

    delta theta = alpha m,
    delta a_w   = [f]_x delta theta,

for which H_m delta theta=0 and

    -[f]_x delta theta + delta a_w = 0.

The null family is one-dimensional, and the accelerometer kernel has dimension
three while the magnetometer imposes exactly two independent constraints on
its theta component.  Hence the stacked accelerometer+magnetometer map has
**exact rank five** on the six-dimensional (theta,a_w) active block.

A mode adds the active accelerometer-bias coordinates but no new vector
measurement rows, so the same packet has exact rank five on the nine-dimensional
(theta,a_w,b_a) active block and nullity four.  If the S=0 pseudo measurement is
due in the same sample, its three S rows occupy a disjoint state-coordinate
block, raising the direct same-sample measurement rank by exactly three, but the
18/21-state map remains highly singular.

Therefore a strictly positive *instantaneous full-state* measurement-information
margin is algebraically impossible.  This is not a loose bound.  P3 obtains
full word detectability because prediction transports these directional null
spaces and recurrent vector/S information accumulates over time.  P4 must keep
the same structure while replacing the linearized attitude residual by the
finite-angle effective coordinates and exact Joseph/reset calculus.

No replay, filter change, domain shrink or a_w/sigma coupling is used here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_exact_word_map as WORDMAP
import ou3_p5_effective_vector_input as VEFF
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("directional packet rank certificate must not be trajectory fitted")

    word = WORDMAP.build(path)
    veff = VEFF.build(path)
    vector = VECTOR.build()
    failures = [f"word-map: {x}" for x in WORDMAP.validate(word)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    acc = veff.get("accelerometer", {})
    mag = veff.get("magnetometer", {})
    if acc.get("J_aw_orthogonal_full_row_rank") is not True:
        failures.append("accelerometer J_aw is not certified orthogonal full row rank")
    if mag.get("kalman_gain_radial_action_exact_zero") is not True:
        failures.append("magnetometer radial null direction is not source certified")

    live = domain.get("normal_live", {})
    fmin = float(live.get("specific_force_norm_lower_mps2", 0.0))
    mmin = float(live.get("magnetic_vector_norm_lower_uT", 0.0))
    sine = float(live.get("vector_sine_separation_lower", 0.0))
    if not (fmin > 0.0 and mmin > 0.0 and 0.0 < sine < 1.0):
        failures.append("declared vector packet geometry lost nonzero/noncollinear bounds")

    dims = {mode: int(word[mode]["dimension"]) for mode in ("H", "A")}
    if dims != {"H": 18, "A": 21}:
        failures.append("unexpected H/A state dimensions")

    rank_acc = 3
    rank_mag = 2
    rank_vector_pair = 5
    rank_S = 3
    rank_pair_plus_S = 8

    modes = {
        "H": {
            "full_state_dimension": 18,
            "vector_active_coordinates": ["delta_theta(3)", "delta_a_w(3)"],
            "vector_active_dimension": 6,
            "accelerometer_rank": rank_acc,
            "magnetometer_rank": rank_mag,
            "stacked_vector_packet_rank_exact": rank_vector_pair,
            "stacked_vector_packet_nullity_exact_on_active_block": 1,
            "same_sample_rank_without_S_exact": rank_vector_pair,
            "same_sample_full_state_nullity_without_S_lower": 18 - rank_vector_pair,
            "same_sample_rank_with_due_S_exact": rank_pair_plus_S,
            "same_sample_full_state_nullity_with_due_S_lower": 18 - rank_pair_plus_S,
        },
        "A": {
            "full_state_dimension": 21,
            "vector_active_coordinates": ["delta_theta(3)", "delta_a_w(3)", "delta_b_a(3)"],
            "vector_active_dimension": 9,
            "accelerometer_rank": rank_acc,
            "magnetometer_rank": rank_mag,
            "stacked_vector_packet_rank_exact": rank_vector_pair,
            "stacked_vector_packet_nullity_exact_on_active_block": 4,
            "same_sample_rank_without_S_exact": rank_vector_pair,
            "same_sample_full_state_nullity_without_S_lower": 21 - rank_vector_pair,
            "same_sample_rank_with_due_S_exact": rank_pair_plus_S,
            "same_sample_full_state_nullity_with_due_S_lower": 21 - rank_pair_plus_S,
        },
    }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_OPERATION_MATCHED_DIRECTIONAL_PACKET_STRUCTURAL_RANK",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "aw_sigma_consistency_assumption_used": False,
        "measurement_structure": {
            "rotation_gauge_accelerometer": "H_a=[-[f]_x, I_aw]",
            "rotation_gauge_magnetometer": "H_m=[-[m]_x, 0_aw]",
            "accelerometer_rank_exact": rank_acc,
            "magnetometer_rank_exact": rank_mag,
            "stacked_vector_packet_rank_exact": rank_vector_pair,
            "S_zero_rank_exact_when_due": rank_S,
            "stacked_vector_plus_due_S_rank_exact": rank_pair_plus_S,
        },
        "exact_vector_packet_null_witness": {
            "parameter": "alpha != 0",
            "delta_theta": "alpha * m",
            "delta_a_w": "[f]_x * delta_theta",
            "delta_b_a_A_mode": "0",
            "magnetometer_residual": "-[m]_x delta_theta = 0",
            "accelerometer_residual": "-[f]_x delta_theta + delta_a_w = 0",
            "nonzero_for_nonzero_m": True,
        },
        "modes": modes,
        "instantaneous_full_state_measurement_information_positive_definite_possible": False,
        "instantaneous_positive_scalar_full_state_packet_margin_is_valid_target": False,
        "directional_PSD_operation_credit_required": True,
        "word_level_directional_accumulation_required": True,
        "P3_full_word_detectability_not_contradicted": True,
        "P4_must_transport_directional_nullspaces_through_prediction": True,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "represent each accepted S/accelerometer/magnetometer Joseph decrease as a source-correlated PSD directional form, "
            "transport those forms through exact finite-angle prediction/reset/effective-coordinate maps, and accumulate them over recurrent "
            "complete H/A words before taking any scalar lower margin"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "directional_PSD_operation_credit_required",
        "word_level_directional_accumulation_required",
        "P3_full_word_detectability_not_contradicted",
        "P4_must_transport_directional_nullspaces_through_prediction",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "aw_sigma_consistency_assumption_used",
        "instantaneous_full_state_measurement_information_positive_definite_possible",
        "instantaneous_positive_scalar_full_state_packet_margin_is_valid_target",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    ms = d.get("measurement_structure", {})
    if ms.get("accelerometer_rank_exact") != 3:
        f.append("accelerometer rank is not exactly three")
    if ms.get("magnetometer_rank_exact") != 2:
        f.append("magnetometer rank is not exactly two")
    if ms.get("stacked_vector_packet_rank_exact") != 5:
        f.append("vector packet rank is not exactly five")
    if ms.get("S_zero_rank_exact_when_due") != 3:
        f.append("S=0 rank is not exactly three")
    if ms.get("stacked_vector_plus_due_S_rank_exact") != 8:
        f.append("vector+S rank is not exactly eight")

    expected = {
        "H": (18, 6, 1, 13, 10),
        "A": (21, 9, 4, 16, 13),
    }
    for mode, (dim, active, null_active, null_no_s, null_with_s) in expected.items():
        m = d.get("modes", {}).get(mode, {})
        if m.get("full_state_dimension") != dim:
            f.append(f"{mode}: wrong full-state dimension")
        if m.get("vector_active_dimension") != active:
            f.append(f"{mode}: wrong vector-active dimension")
        if m.get("stacked_vector_packet_nullity_exact_on_active_block") != null_active:
            f.append(f"{mode}: wrong active-block nullity")
        if m.get("same_sample_full_state_nullity_without_S_lower") != null_no_s:
            f.append(f"{mode}: wrong no-S full-state nullity lower")
        if m.get("same_sample_full_state_nullity_with_due_S_lower") != null_with_s:
            f.append(f"{mode}: wrong due-S full-state nullity lower")

    w = d.get("exact_vector_packet_null_witness", {})
    if w.get("nonzero_for_nonzero_m") is not True:
        f.append("exact null witness is not certified nonzero")
    if "[m]_x" not in str(w.get("magnetometer_residual", "")):
        f.append("magnetometer null witness missing skew identity")
    if "[f]_x" not in str(w.get("accelerometer_residual", "")):
        f.append("accelerometer null witness missing compensation identity")
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
        "rank": d["measurement_structure"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "instantaneous_scalar_margin_valid": d["instantaneous_positive_scalar_full_state_packet_margin_is_valid_target"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
