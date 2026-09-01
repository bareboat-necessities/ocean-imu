#!/usr/bin/env python3
"""Source-uniform two-word lift of the H-mode vector-packet null direction.

The operation-matched route established that one accepted accel+mag packet has
exact rank five on the six-dimensional (theta,a_w) active block.  Its one
nonzero null family is

    delta theta = alpha m,
    delta a_w   = [f]_x delta theta.

This producer proves a useful structural fact for the next word-level stage:
under the declared vector-PE hypotheses the same null vector necessarily has a
nonzero a_w component, and a source-complete *following* word contains the
strict four-S observation of (v,p,S,a_w).  Therefore that packet null cannot be
an unobservable direction of the two-word H observation map.

The calculation is deliberately source-uniform and conservative.  It uses the
maximum gap of one complete word between the PE packet and the beginning of the
following tile and retains only the surviving OU a_w component.  Translation
couplings can only add v/p/S components; they are not credited here.

The emitted numerical products are a structural/raw-coordinate diagnostic, not
the final P4 Lyapunov metric.  Mixed physical coordinates must still be pulled
back through the full covariance/information metric before any usable scalar
P4 margin can be promoted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_implementation_word_language as WORDS
import ou3_p4_directional_packet_rank as RANK
import ou3_p4_h18_differential_operations as DOPS
import ou3_p4_source_node_cells as SOURCE_NODES
import ou3_translational_uco_ucc as TRANS
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _exp_minus_upper(x: float) -> Interval:
    """Validated exp(-x) for finite x>=0 with range reduction if needed."""
    if not math.isfinite(x) or x < 0.0:
        raise ValueError("finite nonnegative exponent required")
    scale = 1
    while x / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    z = Interval.outward_bounds(x / scale, x / scale)
    y = VT.exp_interval(-z)
    s = scale
    while s > 1:
        y = y.square()
        s //= 2
    return y


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("two-word null lift must not be trajectory fitted")

    words = WORDS.build(path)
    wf = WORDS.validate(words)
    rank = RANK.build(path)
    rf = RANK.validate(rank)
    packet_rank = int(rank["measurement_structure"]["stacked_vector_packet_rank_exact"])
    h_active_nullity = int(rank["modes"]["H"]["stacked_vector_packet_nullity_exact_on_active_block"])
    trans = TRANS.build()
    tf = TRANS.validate(trans)
    source_nodes = SOURCE_NODES.build()
    nf = SOURCE_NODES.validate(source_nodes)

    live = domain["normal_live"]
    f_min = float(live["specific_force_norm_lower_mps2"])
    sine_min = float(live["vector_sine_separation_lower"])
    aw_per_theta_lower = math.nextafter(f_min * sine_min, -math.inf)

    # Use the configured sample upper, not the nominal recurrence time, so the
    # gap bound includes one discretization endpoint.  This safely covers a PE
    # packet at the beginning of one tile and the start of the next tile.
    wc = words["word_contract"]["conditional_word_language"]
    samples_upper = int(wc["word_samples_upper_at_configured_dt"])
    dt_hi = float(words["word_contract"]["configured_runtime"]["imu_dt_outward_interval_s"][1])
    one_word_gap_upper_s = math.nextafter(samples_upper * dt_hi, math.inf)

    # Bind the decay to the exact P2 source partition rather than a second
    # global parameter box.  The worst one-word survival occurs at the minimum
    # tau lower endpoint among the 800 source nodes.
    tau_min = min(float(n["tau_s"][0]) for n in source_nodes["nodes"])
    if tau_min <= 0.0:
        raise RuntimeError("P2 source-node tau lower endpoint is not positive")
    decay = _exp_minus_upper(one_word_gap_upper_s / tau_min)
    aw_next_per_theta_lower = math.nextafter(aw_per_theta_lower * decay.lo, -math.inf)

    tr = trans["S_observation_uco"]
    translation_info_lower = float(tr["information_gramian_lambda_min_lower"])
    # For the packet-null family, the following-word translation state has
    # norm at least its surviving a_w component.  Since the validated four-S
    # information Gramian satisfies G_S >= lambda_S I on [v,p,S,a_w],
    # x^T G_S x >= lambda_S ||x||^2 gives a concrete raw-coordinate quadratic
    # credit per ||delta theta||^2.  This is deliberately not converted into
    # the full P4 covariance/information metric here.
    raw_info_per_theta2_lower = math.nextafter(
        translation_info_lower * aw_next_per_theta_lower * aw_next_per_theta_lower,
        -math.inf,
    )

    # Exercise the exact same local H construction consumed by the #450 AD
    # route.  This is an audit witness that #449 is not maintaining a second
    # accelerometer/magnetometer derivative convention.
    f0 = [DOPS.I(0.0), DOPS.I(0.0), DOPS.I(f_min)]
    m_min = float(live["magnetic_vector_norm_lower_uT"])
    m0 = [DOPS.I(m_min * sine_min), DOPS.I(0.0), DOPS.I(0.0)]
    Ha = DOPS.H_acc_canonical(f0)
    Hm = DOPS.H_mag_canonical(m0)
    shared_H_a_aw_identity = all(Ha[i][15 + i].lo <= 1.0 <= Ha[i][15 + i].hi for i in range(3))
    shared_H_m_aw_zero = all(Hm[i][15 + j].lo <= 0.0 <= Hm[i][15 + j].hi for i in range(3) for j in range(3))

    structural = bool(
        not wf
        and not rf
        and not tf
        and not nf
        and source_nodes["partition"]["states"] == 800
        and packet_rank == 5
        and h_active_nullity == 1
        and aw_per_theta_lower > 0.0
        and decay.lo > 0.0
        and aw_next_per_theta_lower > 0.0
        and translation_info_lower > 0.0
        and raw_info_per_theta2_lower > 0.0
        and tr["aligned_firing_count"] == 4
        and shared_H_a_aw_identity
        and shared_H_m_aw_zero
    )

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_TWO_WORD_VECTOR_PACKET_NULL_STRUCTURAL_LIFT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "shared_H18_differential_operations_used": True,
        "exact_P2_source_node_partition_used": True,
        "P2_source_node_count": source_nodes["partition"]["states"],
        "packet_rank_exact": packet_rank,
        "packet_H_active_nullity": h_active_nullity,
        "packet_null_aw_per_theta_norm_lower_mps2_per_rad": aw_per_theta_lower,
        "PE_specific_force_norm_lower_mps2": f_min,
        "PE_vector_sine_separation_lower": sine_min,
        "one_word_gap_upper_s": one_word_gap_upper_s,
        "tau_aw_lower_s": tau_min,
        "aw_survival_factor_to_following_word_lower": decay.lo,
        "following_word_aw_per_theta_norm_lower": aw_next_per_theta_lower,
        "following_word_four_S_information_gramian_lambda_min_lower": translation_info_lower,
        "following_word_four_S_firing_count": tr["aligned_firing_count"],
        "packet_null_following_word_raw_information_per_theta2_lower": raw_info_per_theta2_lower,
        "packet_null_following_word_raw_information_bound_formula": "lambda_min(G_S)*(a_w_next/theta)^2",
        "shared_H_a_aw_identity_verified": shared_H_a_aw_identity,
        "shared_H_m_aw_zero_verified": shared_H_m_aw_zero,
        "two_word_packet_null_is_structurally_observed": structural,
        "raw_coordinate_product_is_final_P4_metric_margin": False,
        "full_H18_metric_directional_credit_established_here": False,
        "complete_source_branch_family_checked_here": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "failures": [
            *[f"word-language: {x}" for x in wf],
            *[f"packet-rank: {x}" for x in rf],
            *[f"translation: {x}" for x in tf],
            *[f"source-nodes: {x}" for x in nf],
        ],
        "next_obligation": (
            "attach the actual source-node covariance/information factors to the 800 P2 nodes; use the shared H18 interval-AD prediction/update/reset maps to pull back the following-word four-S PSD form to the preceding vector-packet null family; combine it with the packet Joseph directional form before scalarization"
        ),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_TWO_WORD_VECTOR_PACKET_NULL_STRUCTURAL_LIFT":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "shared_H18_differential_operations_used",
        "exact_P2_source_node_partition_used",
        "shared_H_a_aw_identity_verified",
        "shared_H_m_aw_zero_verified",
        "two_word_packet_null_is_structurally_observed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "raw_coordinate_product_is_final_P4_metric_margin",
        "full_H18_metric_directional_credit_established_here",
        "complete_source_branch_family_checked_here",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_source_node_count") != 800:
        f.append("P2 source-node count is not 800")
    if d.get("packet_rank_exact") != 5 or d.get("packet_H_active_nullity") != 1:
        f.append("packet rank/nullity contract changed")
    for key in (
        "packet_null_aw_per_theta_norm_lower_mps2_per_rad",
        "one_word_gap_upper_s",
        "tau_aw_lower_s",
        "aw_survival_factor_to_following_word_lower",
        "following_word_aw_per_theta_norm_lower",
        "following_word_four_S_information_gramian_lambda_min_lower",
        "packet_null_following_word_raw_information_per_theta2_lower",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
            f.append(f"{key} is not finite positive")
    if d.get("following_word_four_S_firing_count") != 4:
        f.append("following word is not the complete four-S observation")
    if d.get("packet_null_following_word_raw_information_bound_formula") != "lambda_min(G_S)*(a_w_next/theta)^2":
        f.append("raw packet-null information formula changed")
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
        "structural_lift": d["two_word_packet_null_is_structurally_observed"],
        "P2_source_nodes": d["P2_source_node_count"],
        "aw_per_theta_lower": d["packet_null_aw_per_theta_norm_lower_mps2_per_rad"],
        "word_gap_upper_s": d["one_word_gap_upper_s"],
        "aw_survival_lower": d["aw_survival_factor_to_following_word_lower"],
        "aw_next_per_theta_lower": d["following_word_aw_per_theta_norm_lower"],
        "four_S_information_lower": d["following_word_four_S_information_gramian_lambda_min_lower"],
        "raw_information_per_theta2_lower": d["packet_null_following_word_raw_information_per_theta2_lower"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
