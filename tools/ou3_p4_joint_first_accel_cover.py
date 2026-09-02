#!/usr/bin/env python3
"""Adaptive full-outer-ball first-accelerometer diagnostic for OU-III P4.

The co-gauged joint-Joseph design established that the verified innovation
solve itself is viable once covariance/vector orientation is kept correlated.
Its next obstruction was an interval Cayley composition reaching the antipode
on one deliberately broad cover cell.  This diagnostic answers the immediate
question without shrinking the theorem domain:

* keep the same 0.80 rad outer attitude ball;
* build the existing outward box cover of that ball;
* propagate the source-correlated covariance only once to the first mandatory
  accelerometer packet;
* for each attitude box, propagate the interval-AD state to that packet and
  evaluate the same joint P,H,R,r Joseph solve, never materializing K;
* if and only if the direct deployed Cayley correction reaches the antipode,
  bisect that attitude box and retry its children;
* discard a child only when its minimum Euclidean norm is strictly outside the
  original Cayley ball.

The post-prediction rebase makes coordinates 3..N independent of the initial
attitude cell in this zero-rate source realization.  We therefore propagate
that linear tail once per H/A mode and, for every geometry child, replay only
the three attitude coordinates (plus the gyro-bias coordinates that drive
those attitude predictions).  This is an execution optimization only: the
resulting interval state at the first mandatory packet is identical to the
full 18/21-state replay used by the original diagnostic.

Thus successful subdivision is a dependency refinement, not a smaller basin.
This remains a design diagnostic: body rate is the admissible zero-rate source
realization, the PE vectors are the canonical floor realization, and one P2
source node is used.  It does not promote P4 or P5.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import matrix_add, matrix_mul, matrix_transpose
import ou3_p4_candidate_full_word as CAND
import ou3_p4_joint_word_dissipation_design as D
import ou3_p4_joint_word_gauge_design as G
import ou3_p4_joint_word_gauge_design_v2 as G2
import ou3_p4_source_node_cells as NODES
import ou3_p5_full_h_prefix_cells as H
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = G.DEFAULT_DOMAIN
SCHEMA = 2


def _predict_cov(P, F, Q):
    return H._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P), matrix_transpose(F)), Q))


def _slow_state_at_first_accel(mode: str, domain: dict, cbox, F, first_k: int):
    """Reference implementation used once to build/check the invariant tail."""
    zpre = D._initial_ad(mode, domain, cbox)
    z = G.POST._rebase_postprediction(D._prediction(mode, zpre, F))
    for _ in range(first_k):
        z = D._prediction(mode, z, F)
    return z


def _attitude_bg_step(c, bg, F):
    """Exact attitude/bg portion of D._prediction for the current source cell."""
    n = len(c[0].der)
    Rstep = [[F[i][j] for j in range(3)] for i in range(3)]
    Bstep = [[F[i][3 + j] for j in range(3)] for i in range(3)]
    transported = D._ad_matvec(Rstep, c)
    db = D._ad_matvec(Bstep, bg)
    cnext = D.AD.deployed_correct_cayley(transported, db)
    bnext = []
    for i in range(3, 6):
        y = D.AD.constant(0.0, n)
        for j in range(3, 6):
            y = y + D.AD.constant(F[i][j], n) * bg[j - 3]
        bnext.append(y)
    return cnext, bnext


def _fast_state_at_first_accel(mode: str, domain: dict, cbox, F, first_k: int, tail):
    """Replay only the cell-dependent attitude part after the common rebase."""
    zpre = D._initial_ad(mode, domain, cbox)
    z0 = G.POST._rebase_postprediction(D._prediction(mode, zpre, F))
    c = list(z0[:3])
    bg = list(z0[3:6])
    for _ in range(first_k):
        c, bg = _attitude_bg_step(c, bg, F)
    # The common tail already contains the same final bg variables.  Use it for
    # indices 3..N so the accelerometer residual/Jacobian sees the exact source
    # state that the slow full replay would have produced.
    return c + list(tail[3:])


def _same_ad(a, b) -> bool:
    if len(a) != len(b) or len(a[0].der) != len(b[0].der):
        return False
    for x, y in zip(a, b):
        if x.val.lo != y.val.lo or x.val.hi != y.val.hi:
            return False
        for dx, dy in zip(x.der, y.der):
            if dx.lo != dy.lo or dx.hi != dy.hi:
                return False
    return True


def _first_accel_context(mode: str, path: Path, domain: dict, source_node_index: int):
    CAND._configure_mode(mode)
    n = 18 if mode == "H" else 21
    src = NODES.h18_source_cell(source_node_index, NODES.build())
    F, Q, _meta = G2.corrected_zero_rate_transition_process(mode, src, domain)

    # G._mode starts from one already-predicted covariance/state and performs
    # another prediction at the bottom of every pre-vector sample.  Reproduce
    # that exact schedule once for P, while z remains cell-specific below.
    Ppre = G._gauged_golive_covariance(mode, src, path)
    Pm = _predict_cov(Ppre, F, Q)

    h = float(src["dt_s"])
    words = G.WORDS.build(path)
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])
    schedule = D._schedule(path, samples, h)
    first_k = int(schedule["vector_steps"][0])
    for _ in range(first_k):
        Pm = _predict_cov(Pm, F, Q)

    # Coordinates 3..N are independent of the initial attitude cell after the
    # post-prediction rebase.  Propagate them once.  Also prove the optimized
    # replay is bit-for-bit interval-identical on one nontrivial reference box.
    zero_box = [(0.0, 0.0)] * 3
    tail = _slow_state_at_first_accel(mode, domain, zero_box, F, first_k)
    qref = 2.0 * math.tan(0.80 / 2.0)
    ref_box = [(0.0, 0.25 * qref), (-0.25 * qref, 0.0), (0.0, 0.125 * qref)]
    slow_ref = _slow_state_at_first_accel(mode, domain, ref_box, F, first_k)
    fast_ref = _fast_state_at_first_accel(mode, domain, ref_box, F, first_k, tail)
    if not _same_ad(slow_ref, fast_ref):
        raise RuntimeError("optimized first-accelerometer replay is not interval-identical to full replay")

    force, _mag = D._canonical_vectors(domain)
    Ha = D._H_acc(mode, force, n)
    vc = VECTOR.build()["configured_measurement_bounds"]
    Racc = H._R_diag(float(vc["acc_measurement_std_mps2"]))
    return {
        "n": n,
        "src": src,
        "F": F,
        "Pm": Pm,
        "force": force,
        "Ha": Ha,
        "Racc": Racc,
        "first_k": first_k,
        "schedule": schedule,
        "common_linear_tail": tail,
        "optimized_state_replay_verified_identical": True,
    }


def _split_box(box):
    widths = [float(b - a) for a, b in box]
    axis = max(range(3), key=lambda i: widths[i])
    a, b = box[axis]
    m = 0.5 * (a + b)
    left = list(box)
    right = list(box)
    left[axis] = (a, m)
    right[axis] = (m, b)
    return left, right


def _evaluate_cell(mode: str, domain: dict, ctx: dict, cbox):
    z = _fast_state_at_first_accel(
        mode, domain, cbox, ctx["F"], ctx["first_k"], ctx["common_linear_tail"]
    )
    r = D._exact_acc_residual(mode, z, ctx["force"])
    eta = D._eta(r, D._linear_residual(ctx["Ha"], z))
    _Pout, zout, dform, meta = D._ad_joint_update(
        ctx["Pm"], z, ctx["Ha"], ctx["Racc"], r, eta
    )
    return {
        "correction_theta_norm_upper": float(meta["correction_theta_norm_upper"]),
        "post_cayley_norm_upper": float(D.AD._norm_upper([q.val for q in zout[:3]])),
        "signed_form_diagonal_lower": [float(dform[i][i].lo) for i in range(len(dform))],
        "inverse_backend_q_inf_upper": float(meta["inverse_q_inf_upper"]),
    }


def _mode_cover(mode: str, path: Path, domain: dict, source_node_index: int,
                q: float, max_depth: int, max_cells: int):
    ctx = _first_accel_context(mode, path, domain, source_node_index)
    roots = CAND._ball_box_cover(q, max_box_norm_factor=1.5)
    todo = [(box, 0) for box in roots]
    leaves = 0
    split_count = 0
    outside_discarded = 0
    antipode_splits = 0
    unresolved = []
    other_failures = []
    max_correction = 0.0
    max_post_q = 0.0
    max_depth_used = 0

    while todo:
        if leaves + len(todo) + split_count > max_cells:
            unresolved.append({
                "reason": "MAX_CELL_BUDGET_REACHED",
                "remaining": len(todo),
                "max_cells": max_cells,
            })
            break
        box, depth = todo.pop()
        qmin, _qmax = CAND._norm_bounds_box(box)
        if qmin > q:
            outside_discarded += 1
            continue
        max_depth_used = max(max_depth_used, depth)
        try:
            row = _evaluate_cell(mode, domain, ctx, box)
            leaves += 1
            max_correction = max(max_correction, row["correction_theta_norm_upper"])
            max_post_q = max(max_post_q, row["post_cayley_norm_upper"])
        except Exception as exc:
            msg = f"{type(exc).__name__}: {exc}"
            if "Cayley antipode" in msg and depth < max_depth:
                left, right = _split_box(box)
                todo.append((right, depth + 1))
                todo.append((left, depth + 1))
                split_count += 1
                antipode_splits += 1
                continue
            record = {"depth": depth, "box": box, "failure": msg}
            if "Cayley antipode" in msg:
                unresolved.append(record)
            else:
                other_failures.append(record)
            if len(unresolved) + len(other_failures) >= 20:
                break

    complete = not unresolved and not other_failures and not todo
    return {
        "dimension": ctx["n"],
        "source_node_index": source_node_index,
        "source_realization_body_rate_rad_s": 0.0,
        "gravity_yaw_covariance_gauge": "d=e3",
        "first_accelerometer_step": ctx["first_k"],
        "root_cover_cells": len(roots),
        "accepted_leaf_cells": leaves,
        "split_count": split_count,
        "antipode_split_count": antipode_splits,
        "outside_ball_children_discarded": outside_discarded,
        "max_depth_used": max_depth_used,
        "max_depth_allowed": max_depth,
        "max_cell_budget": max_cells,
        "max_correction_theta_norm_upper": max_correction,
        "max_post_cayley_norm_upper": max_post_q,
        "optimized_state_replay_verified_identical": ctx["optimized_state_replay_verified_identical"],
        "unresolved_antipode_cells": unresolved,
        "other_failures": other_failures,
        "full_outer_ball_first_accel_cover_complete": complete,
        "K_interval_matrix_materialized": False,
    }


def build(path: Path = DEFAULT_DOMAIN, source_node_index: int = 0,
          max_depth: int = 10, max_cells: int = 5000):
    path = Path(path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    outer_angle = 0.80
    q = 2.0 * math.tan(outer_angle / 2.0)
    modes = {}
    failures = []
    for mode in ("H", "A"):
        try:
            modes[mode] = _mode_cover(
                mode, path, domain, source_node_index, q, max_depth, max_cells
            )
        except Exception as exc:
            failures.append(f"{mode}: {type(exc).__name__}: {exc}")
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_CO_GAUGED_ADAPTIVE_FIRST_ACCEL_OUTER_BALL_DIAGNOSTIC",
        "design_diagnostic_only": True,
        "source_complete": False,
        "outer_angle_rad": outer_angle,
        "outer_cayley_norm": q,
        "outer_domain_shrunk": False,
        "adaptive_subdivision_preserves_outer_ball_cover": True,
        "trajectory_replay_used": False,
        "P1_changed": False,
        "P3_delta_used_as_physical_basin": False,
        "joint_P_H_R_r_used": True,
        "K_interval_matrix_materialized": False,
        "modes": modes,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--max-cells", type=int, default=5000)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, a.source_node_index, a.max_depth, a.max_cells)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "complete": d.get("modes", {}).get(m, {}).get("full_outer_ball_first_accel_cover_complete"),
            "roots": d.get("modes", {}).get(m, {}).get("root_cover_cells"),
            "leaves": d.get("modes", {}).get(m, {}).get("accepted_leaf_cells"),
            "splits": d.get("modes", {}).get(m, {}).get("split_count"),
            "max_depth": d.get("modes", {}).get(m, {}).get("max_depth_used"),
            "max_correction": d.get("modes", {}).get(m, {}).get("max_correction_theta_norm_upper"),
            "max_post_q": d.get("modes", {}).get(m, {}).get("max_post_cayley_norm_upper"),
            "unresolved": len(d.get("modes", {}).get(m, {}).get("unresolved_antipode_cells", [])),
            "other_failures": len(d.get("modes", {}).get(m, {}).get("other_failures", [])),
        } for m in ("H", "A")},
        "failures": d["failures"],
    }, indent=2, sort_keys=True))
    # Unresolved antipode cells remain a non-fatal diagnostic result, but any
    # other per-cell exception is an execution failure and must fail closed.
    other_failures = any(
        d.get("modes", {}).get(m, {}).get("other_failures") for m in ("H", "A")
    )
    return 0 if not d["failures"] and not other_failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
