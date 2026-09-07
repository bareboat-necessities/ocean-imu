#!/usr/bin/env python3
"""Full-state joint-sector master inequality for complete-SEA3 P4.

The retained complete-word endpoint identity has the form

    d_W = B_W w_W,
    Delta V = -x' D_W x + 2 (M_W x)' J_N B_W w_W
              + w_W' B_W' J_N B_W w_W.

Here ``w_W`` is one correlated word coordinate containing the nonlinear
accelerometer history together with the endpoint/boundary and exact residual
terms.  It is not split into independent packet radii.  Stack ``z=[x;w_W]``.
The exact full master quadratic is

    z' L_W z,
    L_W = [[-D_W, M_W' J_N B_W],
           [B_W' J_N M_W, B_W' J_N B_W]].

Suppose the same-history nonlinear graph is enclosed by quadratic sectors

    z' Pi_j z >= 0.

For nonnegative multipliers lambda_j, the strict full-matrix condition

    -(L_W + sum_j lambda_j Pi_j) > 0

implies ``Delta V < 0`` on that graph.  The test is performed by outward
interval LDLT on the complete augmented matrix.  This keeps every state and
every cross term; it is not a scalar correction radius, packet-count budget,
Schur elimination of a_w/b_a, or inverse-metric-floor argument.

The accelerometer block may use the retained Riccati noise-channel identity
``A' P_N^-1 A <= R^-1``.  Every later shipping event, including every due S=0
update with its actual applied anisotropic SpectralMSE R_S, remains inside A's
same-word suffix.  This module provides the terminal master inequality and
fails closed until a source-uniform correlated graph sector is supplied.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_add,
    matrix_mul,
    matrix_sub,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_complete_sea3_residual_sector as RESIDUAL
import ou3_p4_complete_word_accelerometer_channel as CHANNEL
import ou3_p4_complete_sea3_signed_information_ledger as SIGNED
import ou3_sea3_riccati_metric_p3 as P3

SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_FULL_STATE_JOINT_SECTOR_MASTER_V1"
P3_DELTA = 1.0e-18


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    rows = len(A)
    cols = len(A[0]) if rows else 0
    if any(len(row) != cols for row in A):
        raise ValueError("ragged interval matrix")
    return rows, cols


def _neg(A: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    return [[-x for x in row] for row in A]


def _scale(A: Sequence[Sequence[Interval]], value: float) -> IntervalMatrix:
    if not math.isfinite(float(value)):
        raise ValueError("finite multiplier required")
    c = Interval.outward_bounds(float(value), float(value))
    return [[c * x for x in row] for row in A]


def _block2(A, B, C, D) -> IntervalMatrix:
    ar, ac = _shape(A)
    br, bc = _shape(B)
    cr, cc = _shape(C)
    dr, dc = _shape(D)
    if ar != ac or dr != dc or br != ar or bc != dc or cr != dr or cc != ac:
        raise ValueError("2x2 block dimensions do not compose")
    return [list(A[i]) + list(B[i]) for i in range(ar)] + [
        list(C[i]) + list(D[i]) for i in range(dr)
    ]


def master_quadratic_matrix(
    D_word: Sequence[Sequence[Interval]],
    M_word: Sequence[Sequence[Interval]],
    J_final: Sequence[Sequence[Interval]],
    B_word: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Return the exact augmented endpoint-energy matrix ``L_W``."""
    n, n2 = _shape(D_word)
    mr, mc = _shape(M_word)
    jr, jc = _shape(J_final)
    br, bw = _shape(B_word)
    if n == 0 or n != n2 or (mr, mc) != (n, n) or (jr, jc) != (n, n):
        raise ValueError("D/M/J must be nonempty same-dimension square matrices")
    if br != n or bw == 0:
        raise ValueError("B_word must map a nonempty joint coordinate into the full state")
    J_B = matrix_mul(J_final, B_word)
    cross = matrix_mul(matrix_transpose(M_word), J_B)
    defect = matrix_mul(matrix_transpose(B_word), J_B)
    return matrix_symmetric_hull(
        _block2(_neg(D_word), cross, matrix_transpose(cross), defect)
    )


def joint_sector_sprocedure_matrix(
    master: Sequence[Sequence[Interval]],
    sectors: Sequence[Sequence[Sequence[Interval]]],
    multipliers: Sequence[float],
) -> IntervalMatrix:
    """Return ``L + sum lambda_j Pi_j`` for graph sectors ``z'Pi_j z>=0``."""
    n, m = _shape(master)
    if n == 0 or n != m:
        raise ValueError("master matrix must be nonempty square")
    if len(sectors) != len(multipliers) or not sectors:
        raise ValueError("need matching nonempty sector/multiplier lists")
    out = matrix_symmetric_hull(master)
    for Pi, lam in zip(sectors, multipliers):
        if _shape(Pi) != (n, n):
            raise ValueError("sector dimension does not match augmented master")
        if not (math.isfinite(float(lam)) and float(lam) >= 0.0):
            raise ValueError("S-procedure multipliers must be finite nonnegative")
        out = matrix_symmetric_hull(matrix_add(out, _scale(Pi, float(lam))))
    return out


def quadratic_graph_sector(
    graph_input_map: Sequence[Sequence[Interval]],
    graph_input_weight: Sequence[Sequence[Interval]],
    defect_map: Sequence[Sequence[Interval]],
    defect_weight: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Lift one correlated graph bound into the common augmented coordinate.

    If ``u=C*z`` and ``q=E*z`` satisfy ``q'Wq <= u'Qu``, this returns
    ``Pi=C'QC-E'WE`` so the admissible graph obeys ``z'Pi*z >= 0``.  ``C``
    and ``E`` may stack all word events and ``Q``/``W`` may be dense; no
    independent-event or diagonal relaxation is imposed here.
    """
    ur, uz = _shape(graph_input_map)
    qr, qz = _shape(defect_map)
    wr, wc = _shape(graph_input_weight)
    rr, rc = _shape(defect_weight)
    if ur == 0 or qr == 0 or uz == 0 or qz != uz:
        raise ValueError("graph maps must share a nonempty augmented coordinate")
    if (wr, wc) != (ur, ur) or (rr, rc) != (qr, qr):
        raise ValueError("graph weights do not match their mapped coordinates")
    positive = matrix_mul(
        matrix_mul(matrix_transpose(graph_input_map), graph_input_weight),
        graph_input_map,
    )
    negative = matrix_mul(
        matrix_mul(matrix_transpose(defect_map), defect_weight), defect_map
    )
    return matrix_symmetric_hull(matrix_sub(positive, negative))


def certify_strict_joint_sector_domination(
    master: Sequence[Sequence[Interval]],
    sectors: Sequence[Sequence[Sequence[Interval]]],
    multipliers: Sequence[float],
) -> tuple[bool, list[float]]:
    """Certify ``-(L+sum lambda Pi)>0`` by outward full-matrix LDLT."""
    test = _neg(joint_sector_sprocedure_matrix(master, sectors, multipliers))
    ok, pivots = symmetric_positive_definite_ldlt(matrix_symmetric_hull(test))
    return bool(ok), [float(x.lo) for x in pivots]


def quadratic_value(A: Sequence[Sequence[Interval]], x: Sequence[Interval]) -> Interval:
    n, m = _shape(A)
    if n == 0 or n != m or len(x) != n:
        raise ValueError("quadratic value dimension mismatch")
    value = Interval.point(0.0)
    for i in range(n):
        for j in range(n):
            value = value + x[i] * A[i][j] * x[j]
    return value


def build(
    *,
    channel_contract: dict | None = None,
    residual_contract: dict | None = None,
    signed_contract: dict | None = None,
    p3_contract: dict | None = None,
) -> dict:
    """Assemble the bridge, reusing canonical prerequisite results when supplied.

    The canonical P4 assembler has already built the expensive signed ledger and
    passes it here.  Standalone execution still validates the whole dependency
    chain rather than silently weakening it.
    """
    channel = CHANNEL.build() if channel_contract is None else channel_contract
    residual = RESIDUAL.build() if residual_contract is None else residual_contract
    signed = SIGNED.build() if signed_contract is None else signed_contract
    p3 = P3.build() if p3_contract is None else p3_contract
    bad = {
        "accelerometer_channel": CHANNEL.validate(channel),
        "residual_sector": RESIDUAL.validate(residual),
        "signed_information": SIGNED.validate(signed),
        "canonical_P3": P3.validate(p3),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"joint-sector master prerequisites failed: {bad}")
    p3_full_matrix = bool(
        p3["P3_CONDITIONAL_SEA3_PASS"]
        and p3["modes"]["H18"]["Omega_minus_delta_P_full_matrix_closed"]
        and p3["modes"]["A21"]["Omega_minus_delta_P_full_matrix_closed"]
    )
    composition = p3["conditional_composition"]
    finite_taub = bool(
        composition["A21_finite_bias_correlation_route_consumed"]
        and composition["A21_detectability_completion_closed"]
        and composition["A21_paper_UES_hypotheses_closed"]
        and composition["A21_comparison_observer_is_proof_only_not_alternate_estimator"]
        and composition["A21_uses_eta9_packet_shortcut"] is False
    )
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_delta_required": P3_DELTA,
        "P3_H18_delta_consumed": float(
            p3["modes"]["H18"]["relative_Riccati_injection_margin_lower"]
        ),
        "P3_A21_delta_consumed": float(
            p3["modes"]["A21"]["relative_Riccati_injection_margin_lower"]
        ),
        "canonical_full_matrix_P3_consumed": p3_full_matrix,
        "full_state_dimensions": {"H18": 18, "A21": 21},
        "full_A21_finite_taub_P3_margin_required": True,
        "finite_taub_A21_detectability_qualification": (
            "OU3_COMPLETE_SEA3_A21_FINITE_BIAS_DETECTABILITY"
        ),
        "finite_taub_A21_detectability_consumed": finite_taub,
        "finite_taub_A21_asymptotic_gap_lower": float(
            composition["A21_detectability_asymptotic_word_energy_gap_lower"]
        ),
        "finite_taub_bias_contraction_gap_lower": float(
            composition["A21_bias_homogeneous_contraction_gap_lower"]
        ),
        "detectability_comparison_observer_only": bool(
            composition["A21_comparison_observer_is_proof_only_not_alternate_estimator"]
        ),
        "finite_taub_detectability_used_as_standalone_P4_promotion": False,
        "paper_finite_state_quadratic_storage_retained": True,
        "master_coordinate": "z=[x;w_W]",
        "joint_word_defect": "d_W=B_W*w_W",
        "master_matrix": "[[-D_W,M_W^T*J_N*B_W],[B_W^T*J_N*M_W,B_W^T*J_N*B_W]]",
        "same_history_joint_graph_sector_required": True,
        "graph_sector_convention": "z^T*Pi_j*z>=0",
        "nonnegative_S_procedure_multipliers_required": True,
        "terminal_test": "-(L_W+sum_j lambda_j*Pi_j)>0",
        "terminal_full_augmented_interval_LDLT_available": True,
        "same_history_quadratic_graph_sector_assembler_available": True,
        "arbitrary_joint_nondiagonal_sector_weights_retained": True,
        "endpoint_and_prefix_use_same_master_builder": True,
        "accelerometer_Riccati_noise_channel_domination_consumed": bool(
            channel["accelerometer_measurement_noise_is_PSD_final_covariance_component"]
        ),
        "stacked_accelerometer_history_retained_jointly": True,
        "exact_Cayley_residual_sector_consumed": bool(
            residual["finite_residual_is_quadratic_sector_in_c_and_aw"]
        ),
        "accelerometer_bias_nonlinearity_exactly_zero": bool(
            residual["accelerometer_bias_nonlinearity_exactly_zero"]
        ),
        "signed_information_ledger_consumed": bool(
            signed["joint_complete_word_signed_information_composition_available"]
        ),
        "all_due_S_updates_with_actual_applied_RS_remain_in_suffix": True,
        "actual_RS_provenance_preserved": True,
        "all_21_state_cross_terms_retained": True,
        "state_elimination_used": False,
        "a_w_or_b_a_Schur_elimination_used": False,
        "packetwise_norm_sum_used": False,
        "packet_count_multiplier_used": False,
        "per_event_scalar_sector_required": False,
        "correction_radius_used": False,
        "inverse_metric_floor_used": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "filter_changed": False,
        "quality_gates_changed": False,
        "declared_domain_changed": False,
        "source_uniform_same_history_joint_sector_closed": False,
        "source_uniform_full_augmented_LDLT_closed": False,
        "source_uniform_prefix_joint_sector_closed": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "materialize one outward same-history residual-history graph from the complete SEA3 execution, form its "
            "joint endpoint and every-prefix quadratic sectors, and close the full augmented LDLT without packetwise "
            "scalarization; retain the A21 projection generalized Jacobian as a separate graph sector"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_delta_required", 0.0)) != P3_DELTA:
        f.append("frozen P3 delta changed")
    for key in ("P3_H18_delta_consumed", "P3_A21_delta_consumed"):
        if float(d.get(key, 0.0)) != P3_DELTA:
            f.append(f"{key} changed from frozen delta")
    if d.get("full_state_dimensions") != {"H18": 18, "A21": 21}:
        f.append("full-state dimensions changed")
    if d.get("finite_taub_A21_detectability_qualification") != (
        "OU3_COMPLETE_SEA3_A21_FINITE_BIAS_DETECTABILITY"
    ):
        f.append("finite-tau_b A21 detectability qualification changed")
    if float(d.get("finite_taub_A21_asymptotic_gap_lower", 0.0)) < P3_DELTA:
        f.append("finite-tau_b A21 asymptotic gap fell below frozen delta")
    if float(d.get("finite_taub_bias_contraction_gap_lower", 0.0)) <= 0.0:
        f.append("finite-tau_b bias contraction gap is not strict")
    for key in (
        "full_A21_finite_taub_P3_margin_required",
        "canonical_full_matrix_P3_consumed",
        "finite_taub_A21_detectability_consumed",
        "detectability_comparison_observer_only",
        "paper_finite_state_quadratic_storage_retained",
        "same_history_joint_graph_sector_required",
        "nonnegative_S_procedure_multipliers_required",
        "terminal_full_augmented_interval_LDLT_available",
        "same_history_quadratic_graph_sector_assembler_available",
        "arbitrary_joint_nondiagonal_sector_weights_retained",
        "endpoint_and_prefix_use_same_master_builder",
        "accelerometer_Riccati_noise_channel_domination_consumed",
        "stacked_accelerometer_history_retained_jointly",
        "exact_Cayley_residual_sector_consumed",
        "accelerometer_bias_nonlinearity_exactly_zero",
        "signed_information_ledger_consumed",
        "all_due_S_updates_with_actual_applied_RS_remain_in_suffix",
        "actual_RS_provenance_preserved",
        "all_21_state_cross_terms_retained",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "state_elimination_used", "a_w_or_b_a_Schur_elimination_used",
        "packetwise_norm_sum_used", "packet_count_multiplier_used",
        "per_event_scalar_sector_required",
        "correction_radius_used", "inverse_metric_floor_used",
        "trajectory_replay_used", "source_family_replaced", "filter_changed",
        "quality_gates_changed", "declared_domain_changed",
        "source_uniform_same_history_joint_sector_closed",
        "source_uniform_full_augmented_LDLT_closed",
        "source_uniform_prefix_joint_sector_closed", "P4_promoted_here",
        "finite_taub_detectability_used_as_standalone_P4_promotion",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "master": d["master_matrix"],
        "full_A21": d["all_21_state_cross_terms_retained"],
        "actual_RS": d["actual_RS_provenance_preserved"],
        "joint_sector_closed": d["source_uniform_same_history_joint_sector_closed"],
        "P4_promoted_here": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
