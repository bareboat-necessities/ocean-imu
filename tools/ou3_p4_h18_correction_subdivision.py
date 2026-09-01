#!/usr/bin/env python3
"""Adaptive first-cell correction-range diagnostic for the #450 H18 route.

The complete H18 interval-AD word currently reaches the deployed-quaternion
correction-range guard on the first coarse outer-ball box.  This producer asks
whether that is a genuine source obstruction or ordinary interval dependency.
It binds one exact P2 source node, propagates the same shipping H=18 covariance
and nonlinear state maps, and bisects only that failing Cayley box.

Unlike the contraction screen, this diagnostic deliberately does not construct
P3 metrics and does not compute spectral norms.  Its sole numerical question is
whether each accepted S/accelerometer/magnetometer correction stays inside the
already validated deployed-quaternion correction range.  Every correction is
attributed to its exact word sample and operation, with its interval norm upper
reported.

No theorem domain is shrunk, no filter quantity is changed, and no successful
subdivision result can promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import matrix_add, matrix_mul, matrix_transpose
import ou3_implementation_word_language as WORDS
import ou3_interval_ad as AD
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_differential_operations as DOPS
import ou3_p4_h18_interval_ad_word as SCREEN
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
import ou3_p5_deployed_quaternion_cayley_cell as QCOMP
import ou3_p5_full_h_prefix_cells as H
import ou3_p5_full_h_prefix_cells_v2 as H2
import ou3_vector_uco_certificate as VECTOR
from ou3_proof_module_state import preserve_module_bindings

DEFAULT_DOMAIN = SCREEN.DEFAULT_DOMAIN
SCHEMA = 1
N = 18


def _norm_upper(v) -> float:
    s = 0.0
    for x in v:
        a = x.val.abs_upper() if isinstance(x, AD.AD) else x.abs_upper()
        s = math.nextafter(s + math.nextafter(a * a, math.inf), math.inf)
    return math.nextafter(math.sqrt(max(0.0, s)), math.inf)


def _split_box(box) -> list[list[tuple[float, float]]]:
    children = [[]]
    for a, b in box:
        m = 0.5 * (float(a) + float(b))
        parts = [
            (float(a), math.nextafter(m, math.inf)),
            (math.nextafter(m, -math.inf), float(b)),
        ]
        children = [prefix + [part] for prefix in children for part in parts]
    return children


def _accepted(Pm, z, Hm, Rm, residual, *, step: int, operation: str):
    cell = H._measurement_cell(Pm, Hm, Rm, [x.val for x in residual])
    K = cell["K"]
    dx = DOPS.ad_matvec_interval(K, residual)
    d = [-x for x in dx[:3]]
    dnorm = _norm_upper(d)
    if dnorm > QCOMP.MAX_CORRECTION_NORM:
        return None, None, {
            "step": int(step),
            "operation": operation,
            "correction_norm_upper_rad": dnorm,
            "validated_correction_norm_limit_rad": float(QCOMP.MAX_CORRECTION_NORM),
            "inverse_backend": cell["inverse_backend"],
        }
    cp = AD.deployed_correct_cayley(z[:3], d)
    out = list(z)
    out[:3] = cp
    for i in range(3, N):
        out[i] = z[i] - dx[i]
    return cell["P_accepted"], out, {
        "step": int(step),
        "operation": operation,
        "correction_norm_upper_rad": dnorm,
        "validated_correction_norm_limit_rad": float(QCOMP.MAX_CORRECTION_NORM),
        "inverse_backend": cell["inverse_backend"],
    }


def _run_child(path: Path, domain: dict, src: dict, cbox, samples: int) -> dict:
    CAND._configure_mode("H")
    F, Q, Rstep = H2._tight_transition_and_Q(src, domain)
    Pm = H2._corrected_initial_covariance(src, path)
    z = SCREEN._initial_state(domain, cbox)
    force, mag, geometry = SCREEN._canonical_vector_cells(domain)
    Hs = H._H_S()
    Ha = DOPS.H_acc_canonical(force)
    Hm = DOPS.H_mag_canonical(mag)
    vc = VECTOR.build()["configured_measurement_bounds"]
    Racc = H._R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = H._R_diag(float(vc["mag_measurement_std_uT"]))
    RS = H._R_S(src)
    h = float(src["dt_s"])
    schedule = SCREEN._mandatory_schedule(WORDS.build(path), samples, h)
    sigma_hi = src["sigma_aw_mps2"].hi
    max_correction = 0.0
    max_record = None

    for k in range(samples):
        Pm = H._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pm), matrix_transpose(F)), Q))
        Pm = SCREEN._aw_sync_covariance_overapprox(Pm, sigma_hi)
        z = DOPS.prediction(z, F, Rstep, domain, h)

        operations = []
        if k in schedule["S_steps"]:
            operations.append(("mandatory_S", Hs, RS, DOPS.S_residual(z)))
        if k in schedule["vector_steps"]:
            operations.append(("mandatory_accelerometer", Ha, Racc, DOPS.accelerometer_residual(z, force)))
            operations.append(("mandatory_magnetometer", Hm, Rmag, DOPS.magnetometer_residual(z, mag)))

        for name, Huse, Ruse, residual in operations:
            P2, z2, rec = _accepted(Pm, z, Huse, Ruse, residual, step=k, operation=name)
            if rec["correction_norm_upper_rad"] >= max_correction:
                max_correction = rec["correction_norm_upper_rad"]
                max_record = dict(rec)
            if P2 is None:
                return {
                    "completed": False,
                    "entry_cayley_box": cbox,
                    "first_failure": rec,
                    "maximum_correction": max_record,
                    "schedule": schedule,
                    "canonical_vector_geometry": geometry,
                }
            Pm, z = P2, z2

    return {
        "completed": True,
        "entry_cayley_box": cbox,
        "first_failure": None,
        "maximum_correction": max_record,
        "schedule": schedule,
        "canonical_vector_geometry": geometry,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0,
          depth: int = 1) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("correction subdivision must not be trajectory fitted")
    if int(depth) != 1:
        raise ValueError("current focused diagnostic supports exactly one bisection depth")

    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 source-node materialization failed: {nf}")
    src = NODES.h18_source_cell(source_node_index, nodes)
    node = NODES.node(source_node_index, nodes)

    sector = SECTOR.build(path)
    sf = SECTOR.validate(sector)
    words = WORDS.build(path)
    wf = WORDS.validate(words)
    q = float(sector["design_cayley_norm_upper"])
    cover = CAND._ball_box_cover(q, max_box_norm_factor=1.5)
    parent = cover[0]
    children = [c for c in _split_box(parent) if CAND._norm_bounds_box(c)[0] <= q]
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])

    rows = []
    with preserve_module_bindings():
        H._source_cell = lambda: src
        for i, child in enumerate(children):
            try:
                row = _run_child(path, domain, src, child, samples)
            except Exception as exc:
                row = {
                    "completed": False,
                    "entry_cayley_box": child,
                    "first_failure": {
                        "step": None,
                        "operation": "unattributed_exception",
                        "error": f"{type(exc).__name__}: {exc}",
                    },
                    "maximum_correction": None,
                }
            row["child_index"] = i
            rows.append(row)

    completed = sum(1 for r in rows if r["completed"])
    finite_failures = [r["first_failure"] for r in rows if r["first_failure"] is not None]
    max_seen = max(
        (float(r["maximum_correction"]["correction_norm_upper_rad"])
         for r in rows if r.get("maximum_correction")),
        default=math.inf,
    )
    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_FAILING_CAYLEY_CELL_CORRECTION_SUBDIVISION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P3_metric_constructed_for_this_diagnostic": False,
        "spectral_norm_computed_for_this_diagnostic": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "outer_cayley_norm_upper": q,
        "source_node_index": int(source_node_index),
        "source_node": node,
        "full_word_samples": samples,
        "parent_cell": parent,
        "subdivision_depth": 1,
        "children_intersecting_outer_ball": len(children),
        "children_completed_without_correction_range_failure": completed,
        "all_children_completed_without_correction_range_failure": completed == len(children),
        "maximum_correction_norm_upper_rad": max_seen,
        "child_results": rows,
        "attributed_failures": finite_failures,
        "adaptive_state_cell_subdivision_required": True,
        "source_cell_subdivision_used_for_this_test": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "if one bisection closes every child, integrate adaptive Cayley-state subdivision into the full H18 generalized-Jacobian screen; otherwise recursively bisect only children whose attributed correction range still exceeds the deployed limit"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_FAILING_CAYLEY_CELL_CORRECTION_SUBDIVISION":
        f.append("wrong qualification")
    for key in ("source_generated_not_trajectory_fit", "adaptive_state_cell_subdivision_required"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "declared_domain_changed",
        "P3_metric_constructed_for_this_diagnostic", "spectral_norm_computed_for_this_diagnostic",
        "source_cell_subdivision_used_for_this_test", "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("outer_angle_rad", math.nan)) != 0.80:
        f.append("outer angle is not exactly 0.80 rad")
    if d.get("source_node_index") != 0:
        f.append("focused CI source node is not node zero")
    if int(d.get("children_intersecting_outer_ball", 0)) <= 1:
        f.append("failing parent cell was not meaningfully subdivided")
    rows = d.get("child_results", [])
    if len(rows) != d.get("children_intersecting_outer_ball"):
        f.append("subdivision child accounting mismatch")
    for r in rows:
        failure = r.get("first_failure")
        if failure and failure.get("operation") == "unattributed_exception":
            f.append("subdivision produced an unattributed exception")
            break
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--depth", type=int, default=1)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve(), source_node_index=args.source_node_index, depth=args.depth)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_pass": not vf,
        "source_node": d["source_node_index"],
        "samples": d["full_word_samples"],
        "children": [d["children_completed_without_correction_range_failure"], d["children_intersecting_outer_ball"]],
        "all_children_close_correction_range": d["all_children_completed_without_correction_range_failure"],
        "max_correction_norm_upper_rad": d["maximum_correction_norm_upper_rad"],
        "attributed_failures": d["attributed_failures"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
