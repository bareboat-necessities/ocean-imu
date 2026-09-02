#!/usr/bin/env python3
"""Current-head first-accelerometer 0.8-rad cover with exact replay parity.

The original fast prefix replay omitted structurally-zero columns 6..N when
updating gyro bias.  Mathematically those terms vanish, but the reference full
interval-AD predictor still evaluates the zero products/additions, so outward
rounding differs by a few ulps and the fail-closed bitwise parity guard trips.

This wrapper preserves the *exact arithmetic sequence* of the full predictor
for gyro-bias rows while retaining the fast attitude-only prefix replay.  It
also keeps inverse reporting backend-neutral: a fixed-pivot interval inverse
legitimately has no Neumann-q metadata, so None is reported rather than cast to
float.  Neither change alters the interval mathematics or proof domain.

The physical model, 0.8-rad outer set, source node, covariance, Joseph update,
and no-K-hull contract are unchanged.  This is still only a first-packet design
certificate and cannot promote complete P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_joint_first_accel_cover as C

SCHEMA = 3
DEFAULT_DOMAIN = C.DEFAULT_DOMAIN
D = C.D


def _attitude_bg_step_full_rounding(c, bg, F):
    """Bit-for-bit attitude/bg portion of the full D._prediction arithmetic."""
    n = len(c[0].der)
    Rstep = [[F[i][j] for j in range(3)] for i in range(3)]
    Bstep = [[F[i][3 + j] for j in range(3)] for i in range(3)]
    transported = D._ad_matvec(Rstep, c)
    db = D._ad_matvec(Bstep, bg)
    cnext = D.AD.deployed_correct_cayley(transported, db)

    # Mirror D._prediction rows 3:6 literally, including the structural-zero
    # columns.  This preserves outward-rounding order and therefore the exact
    # interval object seen by the slow reference replay.
    zeros = [D.AD.constant(0.0, n) for _ in range(max(0, n - 6))]
    bnext = []
    for i in range(3, 6):
        y = D.AD.constant(0.0, n)
        for j in range(3, n):
            zj = bg[j - 3] if j < 6 else zeros[j - 6]
            y = y + D.AD.constant(F[i][j], n) * zj
        bnext.append(y)
    return cnext, bnext


def _evaluate_cell_backend_neutral(mode: str, domain: dict, ctx: dict, cbox):
    """Base evaluator with optional fixed-pivot/Neumann metadata handled safely."""
    z = C._fast_state_at_first_accel(
        mode, domain, cbox, ctx["F"], ctx["first_k"], ctx["common_linear_tail"]
    )
    r = D._exact_acc_residual(mode, z, ctx["force"])
    eta = D._eta(r, D._linear_residual(ctx["Ha"], z))
    _Pout, zout, dform, meta = D._ad_joint_update(
        ctx["Pm"], z, ctx["Ha"], ctx["Racc"], r, eta
    )
    iq = meta.get("inverse_q_inf_upper")
    return {
        "correction_theta_norm_upper": float(meta["correction_theta_norm_upper"]),
        "post_cayley_norm_upper": float(D.AD._norm_upper([x.val for x in zout[:3]])),
        "signed_form_diagonal_lower": [float(dform[i][i].lo) for i in range(len(dform))],
        "inverse_backend_q_inf_upper": None if iq is None else float(iq),
        "inverse_closed_without_neumann_q": iq is None,
    }


# All imported C helpers resolve these module attributes on C at runtime.
C._attitude_bg_step = _attitude_bg_step_full_rounding
C._evaluate_cell = _evaluate_cell_backend_neutral


def _prove_structural_zero_tail(F, n):
    bad = []
    for i in range(3, 6):
        for j in range(6, n):
            x = F[i][j]
            if not (x.lo == 0.0 and x.hi == 0.0):
                bad.append((i, j, x.lo, x.hi))
    return bad


def mode_cover(mode: str, path: Path, source_node_index: int, max_depth: int, max_cells: int):
    path = Path(path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    C.CAND._configure_mode(mode)
    n = 18 if mode == "H" else 21
    src = C.NODES.h18_source_cell(source_node_index, C.NODES.build())
    F, _Q, _meta = C.G2.corrected_zero_rate_transition_process(mode, src, domain)
    structural_bad = _prove_structural_zero_tail(F, n)
    if structural_bad:
        raise RuntimeError(f"gyro-bias rows have nonzero columns >=6: {structural_bad[:10]}")

    q = 2.0 * math.tan(0.80 / 2.0)
    row = C._mode_cover(mode, path, domain, source_node_index, q, max_depth, max_cells)
    row["schema_v3"] = SCHEMA
    row["optimized_replay_full_rounding_sequence_used"] = True
    row["gyro_bias_rows_columns_6_plus_structurally_zero"] = True
    row["inverse_metadata_backend_neutral"] = True
    row["outer_angle_rad"] = 0.80
    row["outer_domain_shrunk"] = False
    row["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"] = False
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--mode", choices=("H", "A"), required=True)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--max-cells", type=int, default=5000)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    row = mode_cover(a.mode, a.domain, a.source_node_index, a.max_depth, a.max_cells)
    d = {
        "schema": SCHEMA,
        "qualification": "OU3_P4_CURRENT_HEAD_INTERVAL_IDENTICAL_FIRST_ACCEL_OUTER_BALL",
        "mode": a.mode,
        "outer_angle_rad": 0.80,
        "outer_domain_shrunk": False,
        "optimized_replay_full_rounding_sequence_used": True,
        "optimized_state_replay_verified_identical": row.get("optimized_state_replay_verified_identical"),
        "gyro_bias_rows_columns_6_plus_structurally_zero": True,
        "inverse_metadata_backend_neutral": True,
        "joint_P_H_R_r_used": True,
        "K_interval_matrix_materialized": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "result": row,
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "mode": a.mode,
        "complete": row.get("full_outer_ball_first_accel_cover_complete"),
        "roots": row.get("root_cover_cells"),
        "leaves": row.get("accepted_leaf_cells"),
        "splits": row.get("split_count"),
        "max_depth": row.get("max_depth_used"),
        "max_correction": row.get("max_correction_theta_norm_upper"),
        "max_post_q": row.get("max_post_cayley_norm_upper"),
        "replay_identical": row.get("optimized_state_replay_verified_identical"),
        "unresolved": row.get("unresolved_antipode_cells"),
        "other_failures": row.get("other_failures"),
    }, indent=2, sort_keys=True))
    ok = (
        row.get("full_outer_ball_first_accel_cover_complete") is True
        and row.get("optimized_state_replay_verified_identical") is True
        and row.get("root_cover_cells") == 96
        and not row.get("unresolved_antipode_cells")
        and not row.get("other_failures")
    )
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
