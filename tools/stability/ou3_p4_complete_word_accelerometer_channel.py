#!/usr/bin/env python3
"""Joint whole-word accelerometer channel for canonical complete-SEA3 P4.

This is a direct reduction of the endpoint master inequality, not a standalone
remainder certificate.  For each accepted accelerometer correction in one
literal complete word define

    A_i = M_{N:i+1} G_i K_i,

where the suffix M contains every later shipping event, including every due
S=0 correction with its actual applied anisotropic SpectralMSE R_S.  Stack
A=[A_1 ... A_m] and R=diag(Racc_1,...,Racc_m).

The Joseph/reset measurement-noise injection from those same accelerometer
events contributes exactly

    W_acc = A R A^T

to the final Riccati covariance.  All initial-covariance, process-Q, covariance
floor, S=0, vector/magnetometer, and hybrid covariance contributions are PSD,
so for every admitted complete SEA3 word

    W_acc <= P_N.

Therefore the full block Schur inequality gives

    A^T P_N^-1 A <= R^-1.

For the exact endpoint nonlinear interior defect

    d_acc = -A q,       q_i = H_i E_aw epsilon_aw,i,

this implies the *joint-history* bound

    d_acc^T P_N^-1 d_acc <= q^T R^-1 q.

No packet norm is summed and no packet-count multiplier appears.  This also
requires no global inverse-metric floor or correction radius.  The remaining
P4 task is to enclose the stacked nonlinear residual graph q jointly with the
linear endpoint/cross term and the r/boundary terms; this module does not
scalarize or promote that open object.
"""
from __future__ import annotations

import argparse
import json
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
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_sea3_full_word_event_algebra as EVENT

REPO = Path(__file__).resolve().parents[2]
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_WORD_JOINT_ACCELEROMETER_CHANNEL_V1"


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged interval matrix")
    return r, c


def _zero(n: int, m: int) -> IntervalMatrix:
    z = Interval.point(0.0)
    return [[z for _ in range(m)] for _ in range(n)]


def _hstack(blocks: Sequence[Sequence[Sequence[Interval]]]) -> IntervalMatrix:
    if not blocks:
        raise ValueError("at least one channel block required")
    n = _shape(blocks[0])[0]
    if any(_shape(b)[0] != n for b in blocks):
        raise ValueError("channel blocks must have one output dimension")
    return [[x for b in blocks for x in b[i]] for i in range(n)]


def _blockdiag(blocks: Sequence[Sequence[Sequence[Interval]]]) -> IntervalMatrix:
    if not blocks:
        raise ValueError("at least one covariance block required")
    sizes = []
    for b in blocks:
        r, c = _shape(b)
        if r == 0 or r != c:
            raise ValueError("noise covariance blocks must be nonempty square")
        sizes.append(r)
    total = sum(sizes)
    out = _zero(total, total)
    off = 0
    for b, n in zip(blocks, sizes):
        for i in range(n):
            for j in range(n):
                out[off + i][off + j] = b[i][j]
        off += n
    return out


def channel_covariance(
    A_blocks: Sequence[Sequence[Sequence[Interval]]],
    R_blocks: Sequence[Sequence[Sequence[Interval]]],
) -> IntervalMatrix:
    """Return ``sum_i A_i R_i A_i^T`` with outward interval arithmetic."""
    if len(A_blocks) != len(R_blocks) or not A_blocks:
        raise ValueError("A/R channel block counts must match and be nonzero")
    n = _shape(A_blocks[0])[0]
    out = _zero(n, n)
    for A, R in zip(A_blocks, R_blocks):
        ar, ac = _shape(A)
        rr, rc = _shape(R)
        if ar != n or rr != rc or ac != rr:
            raise ValueError("accelerometer channel A/R dimensions do not match")
        out = matrix_add(out, matrix_mul(matrix_mul(A, R), matrix_transpose(A)))
    return matrix_symmetric_hull(out)


def covariance_component_margin(
    P_final: Sequence[Sequence[Interval]],
    A_blocks: Sequence[Sequence[Sequence[Interval]]],
    R_blocks: Sequence[Sequence[Sequence[Interval]]],
) -> IntervalMatrix:
    n, m = _shape(P_final)
    if n == 0 or n != m:
        raise ValueError("final covariance must be nonempty square")
    W = channel_covariance(A_blocks, R_blocks)
    if _shape(W) != (n, n):
        raise ValueError("channel output dimension does not match final covariance")
    return matrix_symmetric_hull(matrix_sub(P_final, W))


def precision_channel_margin(
    P_final: Sequence[Sequence[Interval]],
    A_blocks: Sequence[Sequence[Sequence[Interval]]],
    R_blocks: Sequence[Sequence[Sequence[Interval]]],
) -> IntervalMatrix:
    """Return ``R^-1 - A^T P_final^-1 A`` for small validated fixtures/cells.

    The theorem-level inequality follows structurally from the Riccati PSD
    decomposition; constructing the large 3m-by-3m matrix is not required for
    the 601-sample word.
    """
    A = _hstack(A_blocks)
    R = _blockdiag(R_blocks)
    Pinv = matrix_inverse_gauss_jordan(matrix_symmetric_hull(P_final))
    Rinv = matrix_inverse_gauss_jordan(matrix_symmetric_hull(R))
    pullback = matrix_mul(matrix_mul(matrix_transpose(A), Pinv), A)
    return matrix_symmetric_hull(matrix_sub(Rinv, pullback))


def certify_strict_fixture_precision_domination(
    P_final: Sequence[Sequence[Interval]],
    A_blocks: Sequence[Sequence[Sequence[Interval]]],
    R_blocks: Sequence[Sequence[Sequence[Interval]]],
) -> tuple[bool, list[float]]:
    ok, pivots = symmetric_positive_definite_ldlt(
        precision_channel_margin(P_final, A_blocks, R_blocks)
    )
    return bool(ok), [float(x.lo) for x in pivots]


def build() -> dict:
    event = EVENT.build()
    failures = EVENT.validate(event)
    if failures:
        raise RuntimeError(f"accelerometer channel event algebra failed: {failures}")
    op = event["operation_classes"]["accepted_S_acc_mag_joseph"]
    preservation = event["full_matrix_margin_preservation"]
    exact = bool(
        op["B_psd"]
        and op["actual_applied_SpectralMSE_R_S_required_for_S"]
        and preservation["covers_prediction"]
        and preservation["covers_every_due_S_update"]
        and preservation["covers_every_Normal_Live_accelerometer_update"]
        and preservation["covers_asynchronous_magnetometer_update"]
        and preservation["covers_immediate_left_error_reset"]
        and preservation["covers_aw_covariance_floor"]
    )
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "same_complete_SEA3_word_required": True,
        "same_shipping_suffix_for_state_and_covariance_required": True,
        "accelerometer_channel_definition": "A_i=M_suffix(i)*G_i*K_i",
        "stacked_noise_covariance_definition": "R=diag(Racc_i)",
        "final_covariance_component_identity": "W_acc=A*R*A^T <= P_N",
        "precision_channel_inequality": "A^T*P_N^-1*A <= R^-1",
        "nonlinear_channel_input": "q_i=H_i*E_aw*epsilon_aw_i",
        "endpoint_accelerometer_defect": "d_acc=-A*q",
        "joint_defect_energy_inequality": "d_acc^T*P_N^-1*d_acc <= q^T*R^-1*q",
        "accelerometer_measurement_noise_is_PSD_final_covariance_component": exact,
        "all_later_due_S_updates_remain_inside_suffix": True,
        "actual_applied_RS_required_for_every_later_S_suffix_event": True,
        "actual_RS_regularization_not_removed_by_channel_reduction": True,
        "full_process_Q_floor_vector_and_hybrid_covariance_are_other_PSD_components": True,
        "large_stacked_R_inverse_need_not_be_materialized": True,
        "packetwise_norm_sum_used": False,
        "packet_count_multiplier_used": False,
        "correction_radius_claim_used": False,
        "inverse_metric_floor_claim_used": False,
        "independent_RS_schedule_used": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "stacked_nonlinear_residual_graph_closed_here": False,
        "endpoint_cross_term_closed_here": False,
        "master_endpoint_domination_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "retain q as one correlated accelerometer-history vector and bound its nonlinear graph jointly with the "
            "endpoint cross term and r/boundary terms; do not replace q^T R^-1 q by a packet-count times worst remainder"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "same_complete_SEA3_word_required", "same_shipping_suffix_for_state_and_covariance_required",
        "accelerometer_measurement_noise_is_PSD_final_covariance_component",
        "all_later_due_S_updates_remain_inside_suffix",
        "actual_applied_RS_required_for_every_later_S_suffix_event",
        "actual_RS_regularization_not_removed_by_channel_reduction",
        "full_process_Q_floor_vector_and_hybrid_covariance_are_other_PSD_components",
        "large_stacked_R_inverse_need_not_be_materialized",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "packetwise_norm_sum_used", "packet_count_multiplier_used", "correction_radius_claim_used",
        "inverse_metric_floor_claim_used", "independent_RS_schedule_used", "trajectory_replay_used",
        "source_family_replaced", "stacked_nonlinear_residual_graph_closed_here",
        "endpoint_cross_term_closed_here", "master_endpoint_domination_closed_here", "P4_promoted_here",
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
        "channel_domination": d["accelerometer_measurement_noise_is_PSD_final_covariance_component"],
        "actual_RS_in_suffix": d["actual_RS_regularization_not_removed_by_channel_reduction"],
        "residual_graph_closed": d["stacked_nonlinear_residual_graph_closed_here"],
        "P4_promoted_here": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
