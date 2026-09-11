#!/usr/bin/env python3
"""Non-promoting paired Live finite-increment attachment experiment for ALT.

Two estimator executions consume the SAME typed physical sample and SAME JOINT
frontend/source continuation. They may differ in covariance and physical error
state. For every matched literal Joseph event this probe extracts representative
actual P/H/R-derived innovation data and checks both inverse-free solves plus the
rank-three product-port increment identities. It does not claim a source-uniform
interval graph yet.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "tools/stability"), str(ROOT)]

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_sub, matrix_transpose
import ou3_interval_ad as AD
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_full_normal_live_word as WORD
import ou3_innovation_psd_plus_R_inverse as INNOV
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
from tools.stability.ou3_alt_contraction.rank3 import increment_product_ports

QUALIFICATION = "ALT_PAIRED_LIVE_FINITE_INCREMENT_PROBE_V1"


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _identity(n: int, scale: float = 1.0):
    return [[I(scale if i == j else 0.0) for j in range(n)] for i in range(n)]


def _perturbed_spd(n: int, scale: float):
    P = _identity(n, scale)
    P[0][3] = P[3][0] = I(0.03 * scale)
    P[2][15] = P[15][2] = I(-0.02 * scale)
    return P


def _mid(x: Interval) -> float:
    return 0.5 * (x.lo + x.hi)


def _np_matrix(A):
    return np.array([[_mid(x) for x in row] for row in A], dtype=float)


def _np_vector(v):
    return np.array([_mid(x) for x in v], dtype=float)


def _branch_key(image: JOINT.Image):
    return (image.timer_branch, image.stillness_branch, image.wpe_branch)


def _residual(cell):
    z = EVENTS._state_ad(cell.state)
    if cell.kind == "S_zero":
        y = EVENTS.residual_S(z)
    elif cell.kind == "accelerometer":
        y = EVENTS.residual_accelerometer(z, cell.f_hat, cell.R_hat)
    elif cell.kind == "magnetometer":
        y = EVENTS.residual_magnetometer(z, cell.m_body)
    else:
        raise ValueError(f"not a Joseph event: {cell.kind}")
    return _np_vector(AD.values(y))


def _residual_interval(cell):
    z = EVENTS._state_ad(cell.state)
    if cell.kind == "S_zero":
        y = EVENTS.residual_S(z)
    elif cell.kind == "accelerometer":
        y = EVENTS.residual_accelerometer(z, cell.f_hat, cell.R_hat)
    elif cell.kind == "magnetometer":
        y = EVENTS.residual_magnetometer(z, cell.m_body)
    else:
        raise ValueError(f"not a Joseph event: {cell.kind}")
    return list(AD.values(y))


def _col(v):
    return [[x] for x in v]


def _uncol(A):
    return [row[0] for row in A]


def _zero_in(v):
    return all(x.contains(0.0) for x in v)


def _interval_product_graph(a, b):
    """Outward same-pair product-port graph for one literal Joseph event."""
    cover = __import__("ou3_p4_complete_brmm_source_cover_contract")
    ea = EVENTS.source_joseph_event(**cover.joseph_event_kwargs(a))
    eb = EVENTS.source_joseph_event(**cover.joseph_event_kwargs(b))
    H0, H1 = ea["H"], eb["H"]
    N0 = matrix_mul(a.P, matrix_transpose(H0))
    N1 = matrix_mul(b.P, matrix_transpose(H1))
    S0, S1 = ea["S"], eb["S"]
    r0, r1 = _residual_interval(a), _residual_interval(b)
    Sinv0, _ = INNOV.innovation_inverse_psd_plus_R_3x3(S0, a.R)
    Sinv1, _ = INNOV.innovation_inverse_psd_plus_R_3x3(S1, b.R)
    q0 = _uncol(matrix_mul(Sinv0, _col(r0)))
    q1 = _uncol(matrix_mul(Sinv1, _col(r1)))
    dr = [r1[i] - r0[i] for i in range(3)]
    dq = [q1[i] - q0[i] for i in range(3)]
    dS = matrix_sub(S1, S0)
    dN = matrix_sub(N1, N0)
    uS = _uncol(matrix_mul(dS, _col(q0)))
    uN = _uncol(matrix_mul(dN, _col(q0)))
    solve = _uncol(matrix_add(matrix_mul(S1, _col(dq)), _col(uS)))
    solve_defect = [solve[i] - dr[i] for i in range(3)]
    corr = _uncol(matrix_add(matrix_mul(N1, _col(dq)), _col(uN)))
    direct1 = _uncol(matrix_mul(N1, _col(q1)))
    direct0 = _uncol(matrix_mul(N0, _col(q0)))
    direct = [direct1[i] - direct0[i] for i in range(len(direct1))]
    correction_defect = [corr[i] - direct[i] for i in range(len(direct))]
    return {
        "solve_graph_zero_enclosed": _zero_in(solve_defect),
        "correction_graph_zero_enclosed": _zero_in(correction_defect),
        "q0_width_max": max(x.width() for x in q0),
        "uS_width_max": max(x.width() for x in uS),
        "uN_width_max": max(x.width() for x in uN),
        "hard_product_graph_for_this_paired_cell": True,
    }


def _coefficient_record(cell):
    cover = __import__("ou3_p4_complete_brmm_source_cover_contract")
    event = EVENTS.source_joseph_event(**cover.joseph_event_kwargs(cell))
    H = _np_matrix(event["H"])
    S = _np_matrix(event["S"])
    K = _np_matrix(event["K"])
    P = _np_matrix(cell.P)
    N = P @ H.T
    r = _residual(cell)
    q = np.linalg.solve(S, r)
    correction = N @ q
    return {"H": H, "S": S, "K": K, "N": N, "r": r, "q": q,
            "correction": correction}


def _paired_event(a, b):
    if a.kind != b.kind or a.mode != b.mode or a.event_ordinal != b.event_ordinal:
        raise ValueError("paired literal events are not aligned")
    x0, x1 = _coefficient_record(a), _coefficient_record(b)
    dr = x1["r"] - x0["r"]
    dq = x1["q"] - x0["q"]
    dS = x1["S"] - x0["S"]
    dN = x1["N"] - x0["N"]
    uS = dS @ x0["q"]
    uN = dN @ x0["q"]
    lift = increment_product_ports(x1["N"], x1["S"])
    chi = np.r_[dr, dq, uS, uN]
    solve_defect = lift.equality @ chi
    correction_increment = lift.correction_difference @ chi
    actual_increment = x1["correction"] - x0["correction"]
    interval_graph = _interval_product_graph(a, b)
    return {
        "mode": a.mode,
        "kind": a.kind,
        "ordinal": a.event_ordinal,
        "residual0_norm": float(np.linalg.norm(x0["r"])),
        "residual1_norm": float(np.linalg.norm(x1["r"])),
        "delta_N_norm": float(np.linalg.norm(dN)),
        "delta_S_norm": float(np.linalg.norm(dS)),
        "uN_norm": float(np.linalg.norm(uN)),
        "uS_norm": float(np.linalg.norm(uS)),
        "solve_increment_defect_norm": float(np.linalg.norm(solve_defect)),
        "correction_increment_defect_norm": float(
            np.linalg.norm(correction_increment - actual_increment)),
        "finite_increment_identity_pass": bool(
            np.linalg.norm(solve_defect) <= 2e-10
            and np.linalg.norm(correction_increment - actual_increment) <= 2e-10
        ),
        "interval_product_graph": interval_graph,
    }


def _make_sample():
    return KERNEL.SampleCoordinates(
        gyro_measurement=KERNEL.MAHONY.Vec3(I(.01), I(-.02), I(.005)),
        omega_body_corrected=(I(.01), I(-.02), I(.005)),
        specific_force=KERNEL.MAHONY.Vec3(I(.2), I(-.1), I(-9.75)),
        f_cog_body=(I(0), I(0), I(-9.80665)),
        R_wb=_identity(3), due_S=True, aw_floor_requested=True,
        magnetometer_events_after_imu=(
            KERNEL.MagneticEvent((I(20), I(0), I(40))),
        ),
    )


def _state_pair(n: int):
    a = [I(0) for _ in range(n)]
    b = [I(0) for _ in range(n)]
    a[0], a[12], a[15] = I(.01), I(.03), I(.02)
    b[0], b[1], b[12], b[15], b[16] = (
        I(.015), I(-.007), I(.05), I(.027), I(-.01))
    if n == 21:
        a[18] = I(.02)
        b[18], b[19] = I(.03), I(-.01)
    return a, b


def build() -> dict:
    js = JOINT._smoke_state()
    sample = _make_sample()
    H0, H1 = _state_pair(18)
    A0, A1 = _state_pair(21)
    branch0 = KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(js.frontend),
        H=WORD.initialize_word("H", _identity(18)),
        A=WORD.initialize_word("A", _identity(21)),
        source_cell_id="pair-root-0")
    branch1 = KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(js.frontend),
        H=WORD.initialize_word("H", _perturbed_spd(18, 1.08)),
        A=WORD.initialize_word("A", _perturbed_spd(21, 1.05)),
        source_cell_id="pair-root-1")

    common = dict(
        joint_state=js, sample=sample, radial_scale=Interval(0.0, 1.0),
        true_bias=[I(.04), I(-.02), I(.01)], bias_projection_limit=.4,
        tau_ba=I(1800), sample_index=0)
    left = ATTACH.synchronize_sample(
        branch=branch0, state_in_H=H0, state_in_A=A0,
        next_cell_prefix="paired-left", **common)
    right = ATTACH.synchronize_sample(
        branch=branch1, state_in_H=H1, state_in_A=A1,
        next_cell_prefix="paired-right", **common)

    L = {_branch_key(h.image): (h, a) for h, a in left}
    R = {_branch_key(h.image): (h, a) for h, a in right}
    if set(L) != set(R):
        raise RuntimeError(
            "paired JOINT hybrid branch labels diverged under the same frontend/source")

    event_reports = []
    for key in sorted(L):
        for lh, rh in ((L[key][0], R[key][0]), (L[key][1], R[key][1])):
            if tuple(c.kind for c in lh.cells) != tuple(c.kind for c in rh.cells):
                raise RuntimeError("paired literal event orders diverged")
            for c0, c1 in zip(lh.cells, rh.cells):
                if c0.kind in ("S_zero", "accelerometer", "magnetometer"):
                    event_reports.append(_paired_event(c0, c1))

    all_identity = bool(event_reports) and all(
        x["finite_increment_identity_pass"] for x in event_reports)
    nontrivial_gain = any(
        x["delta_N_norm"] > 1e-9 or x["delta_S_norm"] > 1e-9
        for x in event_reports)
    nonzero_residual = any(
        x["residual0_norm"] > 1e-9 and x["residual1_norm"] > 1e-9
        for x in event_reports)
    interval_graphs = bool(event_reports) and all(
        x["interval_product_graph"]["solve_graph_zero_enclosed"]
        and x["interval_product_graph"]["correction_graph_zero_enclosed"]
        and x["interval_product_graph"]["hard_product_graph_for_this_paired_cell"]
        for x in event_reports)
    return {
        "qualification": QUALIFICATION,
        "scope": "single-sample paired Live execution; non-promoting",
        "same_COMPLETE_BRMM_typed_sample_used_by_both_executions": True,
        "same_joint_frontend_predecessor_used_by_both_executions": True,
        "paired_hybrid_branch_labels_match": set(L) == set(R),
        "literal_H18_A21_event_order_paired": True,
        "paired_r_N_S_q_materialized": True,
        "nonzero_nominal_residual_exercised": nonzero_residual,
        "endogenous_N_or_S_increment_exercised": nontrivial_gain,
        "hard_product_port_values_materialized_from_same_pair": True,
        "rank3_product_port_finite_increment_identities_hold": all_identity,
        "outward_product_graphs_close_for_all_materialized_paired_cells": interval_graphs,
        "event_reports": event_reports,
        "source_uniform_interval_product_graph_closed": False,
        "full_three_second_word_closed": False,
        "all_guard_branches_closed": False,
        "BIAS0_BIAS1_BIAS2_family_attachment_closed": False,
        "ALT_ACTUAL_SOURCE_UNIFORM_FINITE_INCREMENT_WORD_ATTACHED": False,
        "ALT_LIVE_PASS": False,
        "ALT_STARTUP_PASS": False,
        "ALT_END_TO_END_PASS": False,
        "P4_promoted": False,
        "next_falsifiable_experiment": (
            "construct the paired relation over the source-reachable Live predecessor family "
            "and compose matched hard product graphs over the literal three-second word, "
            "retaining all hybrid branch labels and BIAS ancestry without startup coupling"
        ),
    }


def validate(d: dict) -> list[str]:
    f = []
    for key in (
        "same_COMPLETE_BRMM_typed_sample_used_by_both_executions",
        "same_joint_frontend_predecessor_used_by_both_executions",
        "paired_hybrid_branch_labels_match",
        "literal_H18_A21_event_order_paired",
        "paired_r_N_S_q_materialized",
        "nonzero_nominal_residual_exercised",
        "endogenous_N_or_S_increment_exercised",
        "hard_product_port_values_materialized_from_same_pair",
        "rank3_product_port_finite_increment_identities_hold",
        "outward_product_graphs_close_for_all_materialized_paired_cells",
    ):
        if d.get(key) is not True:
            f.append(key + " not true")
    for key in (
        "source_uniform_interval_product_graph_closed",
        "full_three_second_word_closed",
        "all_guard_branches_closed",
        "BIAS0_BIAS1_BIAS2_family_attachment_closed",
        "ALT_ACTUAL_SOURCE_UNIFORM_FINITE_INCREMENT_WORD_ATTACHED",
        "ALT_LIVE_PASS", "ALT_STARTUP_PASS", "ALT_END_TO_END_PASS",
        "P4_promoted",
    ):
        if d.get(key) is not False:
            f.append(key + " not false")
    return f


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "events": len(d["event_reports"]),
        "identities": d["rank3_product_port_finite_increment_identities_hold"],
        "source_uniform": d["source_uniform_interval_product_graph_closed"],
        "failures": f,
    }, indent=2))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
