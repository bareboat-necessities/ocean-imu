#!/usr/bin/env python3
"""Evaluate the exact centered Kalman correction on PSD innovation cells.

For every PSD-aware innovation subcell where the verified S inverse exists, use

    K r = K0 r + (P H^T - K0 S) S^-1 r

instead of first forming an entrywise interval box for K.  The residual matrix
B=P H^T-K0 S is intersected with the algebraically equivalent

    B=(I-K0 H)P H^T-K0 R.

This is an exact rearrangement for each source tuple.  It keeps the correction
center, gain-equation residual, innovation solve and measurement residual tied
until the final matrix-vector operations.  The 4^3 S partition is the same
intermediate PSD-aware partition used by the sibling diagnostic; no physical
state/source domain is reduced and P4 is never promoted here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import (
    Interval,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_point,
    matrix_sub,
)
import ou3_correlated_kalman_gain as CG
import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_correction_subdivision as SUB
import ou3_p4_h18_innovation_psd_subdivision as PSD
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
import ou3_p5_full_h_prefix_cells as H
import ou3_verified_spd_inverse as VINV
from ou3_proof_module_state import preserve_module_bindings

DEFAULT_DOMAIN = SUB.DEFAULT_DOMAIN
SCHEMA = 1


def _norm_upper(v) -> float:
    s = 0.0
    for x in v:
        a = x.abs_upper()
        s = math.nextafter(s + math.nextafter(a * a, math.inf), math.inf)
    return math.nextafter(math.sqrt(max(0.0, s)), math.inf)


def _mat_vec(A, x):
    out = []
    for row in A:
        y = Interval.point(0.0)
        for a, b in zip(row, x):
            y = y + a * b
        out.append(y)
    return out


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("centered correction must not be trajectory fitted")

    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 source nodes invalid: {nf}")
    src = NODES.h18_source_cell(source_node_index, nodes)
    sector = SECTOR.build(path)
    sf = SECTOR.validate(sector)
    words = WORDS.build(path)
    wf = WORDS.validate(words)
    q_outer = float(sector["design_cayley_norm_upper"])
    parent = CAND._ball_box_cover(q_outer, max_box_norm_factor=1.5)[0]
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])

    captured = []
    original_cell = H._measurement_cell

    def wrapped(Pm, Hm, Rm, residual):
        cell = original_cell(Pm, Hm, Rm, residual)
        if not captured and cell["inverse_backend"] == "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE":
            PHt, S = H._innovation(Pm, Hm, Rm)
            captured.append((Pm, Hm, Rm, list(residual), PHt, S))
        return cell

    with preserve_module_bindings():
        H._source_cell = lambda: src
        H._measurement_cell = wrapped
        word = SUB._run_child(path, domain, src, parent, samples)

    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]
    if not captured:
        failures.append("no limiting innovation captured")
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P4_H18_CENTERED_GAIN_CORRECTION_ON_PSD_CELLS",
            "failures": failures,
            "P4_USABLE_CERTIFICATE_PROMOTED": False,
        }

    Pm, Hm, Rm, residual, PHt, S = captured[0]
    p00 = PSD._split(S[0][0])
    p11 = PSD._split(S[1][1])
    p01 = PSD._split(S[0][1])
    rows = []
    pruned = 0
    inverse_ok = 0
    centered_ok = 0
    direct_ok = 0

    for i, a in enumerate(p00):
        for j, b in enumerate(p11):
            for k, c in enumerate(p01):
                Sc = [[S[r][q] for q in range(3)] for r in range(3)]
                Sc[0][0] = a
                Sc[1][1] = b
                Sc[0][1] = Sc[1][0] = c
                Sc = PSD._psd_clip_top(Sc, Rm)
                if Sc is None:
                    pruned += 1
                    continue
                rec = {"index": [i, j, k]}
                try:
                    Sinv, meta = VINV.inverse_enclosure(Sc)
                except Exception as exc:
                    rec["inverse_certified"] = False
                    rec["inverse_error"] = f"{type(exc).__name__}: {exc}"
                    rows.append(rec)
                    continue
                inverse_ok += 1
                rec["inverse_certified"] = True
                rec["inverse_q_inf_upper"] = meta["neumann_q_inf_upper"]

                # Direct natural extension retained for comparison.
                Kdirect = matrix_mul(PHt, Sinv)
                ddirect = _mat_vec(Kdirect, residual)
                direct_norm = _norm_upper(ddirect[:3])
                rec["direct_theta_correction_norm_upper_rad"] = direct_norm
                direct_ok += int(direct_norm <= 6.0)

                # Exact centered correction identity.
                K0f = CG._point_gain_from_PHt_S(PHt, Sc)
                K0 = matrix_point(K0f)
                B1 = matrix_sub(PHt, matrix_mul(K0, Sc))
                IminusK0H = matrix_sub(matrix_identity(len(Pm)), matrix_mul(K0, Hm))
                B2 = matrix_sub(matrix_mul(IminusK0H, PHt), matrix_mul(K0, Rm))
                B = CG._matrix_intersection(B1, B2)
                x = _mat_vec(Sinv, residual)
                d0 = _mat_vec(K0, residual)
                derr = _mat_vec(B, x)
                dcenter = [a0 + e0 for a0, e0 in zip(d0, derr)]
                centered_norm = _norm_upper(dcenter[:3])
                centered_ok += int(centered_norm <= 6.0)
                rec.update({
                    "centered_theta_correction_norm_upper_rad": centered_norm,
                    "centered_theta_component_intervals": [[z.lo, z.hi] for z in dcenter[:3]],
                    "gain_equation_residual_theta_row_norm_upper": [PSD._norm_upper(row) for row in B[:3]],
                    "innovation_solution_norm_upper": _norm_upper(x),
                })
                rows.append(rec)

    retained = len(rows)
    certified_rows = [x for x in rows if x.get("inverse_certified")]
    max_centered = max((x["centered_theta_correction_norm_upper_rad"] for x in certified_rows), default=math.inf)
    max_direct = max((x["direct_theta_correction_norm_upper_rad"] for x in certified_rows), default=math.inf)

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_CENTERED_GAIN_CORRECTION_ON_PSD_CELLS",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "source_node_index": int(source_node_index),
        "full_word_samples": samples,
        "captured_word_failure": word.get("first_failure"),
        "innovation_partition_coordinates": ["S00", "S11", "S01=S10"],
        "pieces_per_coordinate": 4,
        "cartesian_cell_count": 64,
        "PSD_impossible_cells_pruned": pruned,
        "PSD_compatible_cells_retained": retained,
        "verified_inverse_cells": inverse_ok,
        "direct_correction_cells_at_most_6rad": direct_ok,
        "centered_correction_cells_at_most_6rad": centered_ok,
        "max_direct_theta_correction_norm_upper_rad": max_direct,
        "max_centered_theta_correction_norm_upper_rad": max_centered,
        "centered_exact_identity_used": "K r = K0 r + (PHt-K0 S) S^-1 r",
        "gain_equation_residual_intersection_used": True,
        "innovation_subdivision_is_intermediate_proof_partition_only": True,
        "physical_state_or_source_domain_shrunk": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "cell_results": rows,
        "next_obligation": (
            "compare the centered correction with the four-way source-force subdivision; if the centered form closes the already invertible cells, combine it with only the force/S cells still needing inverse certification rather than refining the physical 0.80-rad state domain"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_CENTERED_GAIN_CORRECTION_ON_PSD_CELLS":
        f.append("wrong qualification")
    for key in ("source_generated_not_trajectory_fit", "gain_equation_residual_intersection_used", "innovation_subdivision_is_intermediate_proof_partition_only"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed", "declared_domain_changed", "physical_state_or_source_domain_shrunk", "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("outer_angle_rad", math.nan)) != 0.80:
        f.append("outer angle is not exactly 0.80 rad")
    if d.get("source_node_index") != 0 or d.get("full_word_samples") != 202:
        f.append("focused exact-node/full-word contract changed")
    if d.get("cartesian_cell_count") != 64:
        f.append("expected 4^3 innovation partition")
    if d.get("centered_exact_identity_used") != "K r = K0 r + (PHt-K0 S) S^-1 r":
        f.append("centered correction identity changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_node_index=a.source_node_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_failures": vf,
        "pruned": d.get("PSD_impossible_cells_pruned"),
        "retained": d.get("PSD_compatible_cells_retained"),
        "inverse_ok": d.get("verified_inverse_cells"),
        "direct_pass": d.get("direct_correction_cells_at_most_6rad"),
        "centered_pass": d.get("centered_correction_cells_at_most_6rad"),
        "max_direct": d.get("max_direct_theta_correction_norm_upper_rad"),
        "max_centered": d.get("max_centered_theta_correction_norm_upper_rad"),
        "next": d.get("next_obligation"),
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
