#!/usr/bin/env python3
"""PSD-aware subdivision of the limiting H18 accelerometer innovation block.

The exact-node anatomy shows that the Neumann obstruction lives in the top 2x2
innovation block: S00/S11 have radius about 16.17 and S01 radius about 30.90,
while the third row already has q<1.  This diagnostic therefore partitions only
(S00,S11,S01), not the physical 0.80-rad state domain.

Every exact innovation satisfies A=S-R=H P H^T >= 0.  Each Cartesian innovation
subcell is intersected with the 2x2 principal-minor consequences of A>=0; cells
that cannot contain a PSD A are rigorously discarded.  On every surviving cell
we try both:

* the verified midpoint-Neumann inverse followed by K=P H^T S^-1; and
* the shared correlation-preserving gain equation conditioned on that S cell.

The union of surviving cells covers the original innovation enclosure subject
to the exact PSD invariant.  This producer is diagnostic only and cannot promote
P4 or shrink the source/state theorem domain.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import (
    Interval,
    down,
    up,
    matrix_mul,
    symmetric_gershgorin_upper,
)
import ou3_correlated_kalman_gain as CG
import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_correction_subdivision as SUB
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
import ou3_p5_full_h_prefix_cells as H
import ou3_verified_spd_inverse as VINV
from ou3_proof_module_state import preserve_module_bindings

DEFAULT_DOMAIN = SUB.DEFAULT_DOMAIN
SCHEMA = 1
PIECES = 4


def _split(x: Interval, pieces: int = PIECES) -> list[Interval]:
    if pieces < 1:
        raise ValueError("pieces must be positive")
    width = x.hi - x.lo
    cuts = [x.lo + width * k / pieces for k in range(pieces + 1)]
    out = []
    for k in range(pieces):
        lo = x.lo if k == 0 else max(x.lo, down(cuts[k]))
        hi = x.hi if k == pieces - 1 else min(x.hi, up(cuts[k + 1]))
        out.append(Interval(lo, hi))
    return out


def _min_abs(x: Interval) -> float:
    if x.lo <= 0.0 <= x.hi:
        return 0.0
    return min(abs(x.lo), abs(x.hi))


def _clip(x: Interval, lo: float, hi: float):
    a = max(x.lo, lo)
    b = min(x.hi, hi)
    return None if a > b else Interval(a, b)


def _psd_clip_top(S, R):
    """Intersect one S cell with necessary 2x2 conditions for S-R >= 0."""
    out = [[S[i][j] for j in range(3)] for i in range(3)]
    s00 = _clip(out[0][0], R[0][0].lo, math.inf)
    s11 = _clip(out[1][1], R[1][1].lo, math.inf)
    s22 = _clip(out[2][2], R[2][2].lo, math.inf)
    if s00 is None or s11 is None or s22 is None:
        return None
    out[0][0], out[1][1], out[2][2] = s00, s11, s22

    a00_hi = max(0.0, s00.hi - R[0][0].lo)
    a11_hi = max(0.0, s11.hi - R[1][1].lo)
    off_cap = up(math.sqrt(up(a00_hi * a11_hi)))
    s01 = _clip(out[0][1], -off_cap, off_cap)
    if s01 is None:
        return None
    out[0][1] = out[1][0] = s01

    # If the off-diagonal cell stays away from zero, the principal-minor
    # inequality also supplies conditional lower bounds on both diagonals.
    amin = _min_abs(s01)
    if amin > 0.0:
        a2 = up(amin * amin)
        if a11_hi <= 0.0 or a00_hi <= 0.0:
            return None
        s00 = _clip(s00, down(R[0][0].lo + a2 / a11_hi), math.inf)
        s11 = _clip(s11, down(R[1][1].lo + a2 / a00_hi), math.inf)
        if s00 is None or s11 is None:
            return None
        out[0][0], out[1][1] = s00, s11

    return out


def _norm_upper(v) -> float:
    s = 0.0
    for x in v:
        a = x.abs_upper()
        s = up(s + up(a * a))
    return up(math.sqrt(s))


def _quad_upper(v, A) -> float:
    z = Interval.point(0.0)
    for i in range(len(v)):
        for j in range(len(v)):
            z = z + v[i] * A[i][j] * v[j]
    return up(max(0.0, z.hi))


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("innovation subdivision must not be trajectory fitted")

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
        row = SUB._run_child(path, domain, src, parent, samples)

    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]
    if not captured:
        failures.append("no limiting accelerometer innovation was captured")
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P4_H18_PSD_AWARE_INNOVATION_SUBDIVISION",
            "failures": failures,
            "P4_USABLE_CERTIFICATE_PROMOTED": False,
        }

    Pm, Hm, Rm, residual, PHt, S = captured[0]
    p00 = _split(S[0][0])
    p11 = _split(S[1][1])
    p01 = _split(S[0][1])
    theta_block = [[Pm[i][j] for j in range(3)] for i in range(3)]
    theta_lambda_upper = symmetric_gershgorin_upper(theta_block)

    cells = []
    pruned = 0
    inverse_ok = 0
    inverse_fail = 0
    direct_pass = 0
    correlated_pass = 0
    energy_pass = 0
    for i, a in enumerate(p00):
        for j, b in enumerate(p11):
            for k, c in enumerate(p01):
                Sc = [[S[r][q] for q in range(3)] for r in range(3)]
                Sc[0][0] = a
                Sc[1][1] = b
                Sc[0][1] = Sc[1][0] = c
                Sc = _psd_clip_top(Sc, Rm)
                if Sc is None:
                    pruned += 1
                    continue
                rec = {"index": [i, j, k]}
                try:
                    Sinv, meta = VINV.inverse_enclosure(Sc)
                    inverse_ok += 1
                    Kdirect = matrix_mul(PHt, Sinv)
                    dx_direct = H._mat_vec(Kdirect, residual)
                    direct_norm = _norm_upper(dx_direct[:3])
                    alpha = _quad_upper(residual, Sinv)
                    energy_norm = up(math.sqrt(up(max(0.0, theta_lambda_upper) * alpha)))
                    rec.update({
                        "inverse_certified": True,
                        "inverse_q_inf_upper": meta["neumann_q_inf_upper"],
                        "direct_theta_correction_norm_upper_rad": direct_norm,
                        "theta_energy_correction_norm_upper_rad": energy_norm,
                    })
                    direct_pass += int(direct_norm <= 6.0)
                    energy_pass += int(energy_norm <= 6.0)
                except Exception as exc:
                    inverse_fail += 1
                    rec.update({
                        "inverse_certified": False,
                        "inverse_error": f"{type(exc).__name__}: {exc}",
                    })

                try:
                    cg = CG.gain_enclosure(Pm, Hm, Rm, S_condition=Sc)
                    dx_corr = H._mat_vec(cg["K"], residual)
                    corr_norm = _norm_upper(dx_corr[:3])
                    rec.update({
                        "correlated_gain_certified": True,
                        "correlated_theta_correction_norm_upper_rad": corr_norm,
                        "correlated_theta_row_residual_radii": cg["row_gain_radius_upper"][:3],
                        "correlated_theta_row_psd_norm_upper": cg["row_psd_gain_norm_upper"][:3],
                    })
                    correlated_pass += int(corr_norm <= 6.0)
                except Exception as exc:
                    rec.update({
                        "correlated_gain_certified": False,
                        "correlated_gain_error": f"{type(exc).__name__}: {exc}",
                    })
                cells.append(rec)

    surviving = len(cells)
    max_direct = max((x.get("direct_theta_correction_norm_upper_rad", 0.0) for x in cells), default=math.inf)
    max_corr = max((x.get("correlated_theta_correction_norm_upper_rad", 0.0) for x in cells), default=math.inf)
    max_energy = max((x.get("theta_energy_correction_norm_upper_rad", 0.0) for x in cells), default=math.inf)
    max_q = max((x.get("inverse_q_inf_upper", 0.0) for x in cells), default=math.inf)

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_PSD_AWARE_INNOVATION_SUBDIVISION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "source_node_index": int(source_node_index),
        "full_word_samples": samples,
        "captured_word_failure": row.get("first_failure"),
        "innovation_subdivision_coordinates": ["S00", "S11", "S01=S10"],
        "pieces_per_coordinate": PIECES,
        "cartesian_cell_count": PIECES ** 3,
        "PSD_impossible_cells_pruned": pruned,
        "PSD_compatible_cells_retained": surviving,
        "verified_inverse_cells": inverse_ok,
        "unverified_inverse_cells": inverse_fail,
        "all_PSD_compatible_cells_have_verified_inverse": inverse_ok == surviving,
        "verified_inverse_max_q_inf_upper": max_q,
        "direct_gain_cells_with_correction_at_most_6rad": direct_pass,
        "correlated_gain_cells_with_correction_at_most_6rad": correlated_pass,
        "energy_bound_cells_with_correction_at_most_6rad": energy_pass,
        "max_direct_theta_correction_norm_upper_rad": max_direct,
        "max_correlated_theta_correction_norm_upper_rad": max_corr,
        "max_theta_energy_correction_norm_upper_rad": max_energy,
        "theta_covariance_gershgorin_lambda_upper": theta_lambda_upper,
        "innovation_subdivision_is_intermediate_proof_partition_only": True,
        "physical_state_or_source_domain_shrunk": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "cell_results": cells,
        "next_obligation": (
            "if every PSD-compatible top-block cell verifies S^-1 and the conditioned correction is within range, lift this intermediate innovation partition into the full word; otherwise recursively split only the surviving failing S00/S11/S01 cells or preserve the covariance-reset correlations that created their widths"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_PSD_AWARE_INNOVATION_SUBDIVISION":
        f.append("wrong qualification")
    for key in ("source_generated_not_trajectory_fit", "innovation_subdivision_is_intermediate_proof_partition_only"):
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
    if int(d.get("PSD_compatible_cells_retained", 0)) <= 0:
        f.append("no PSD-compatible innovation cell retained")
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
        "cells": d.get("cartesian_cell_count"),
        "pruned": d.get("PSD_impossible_cells_pruned"),
        "retained": d.get("PSD_compatible_cells_retained"),
        "inverse_ok": d.get("verified_inverse_cells"),
        "inverse_fail": d.get("unverified_inverse_cells"),
        "max_q": d.get("verified_inverse_max_q_inf_upper"),
        "direct_pass": d.get("direct_gain_cells_with_correction_at_most_6rad"),
        "correlated_pass": d.get("correlated_gain_cells_with_correction_at_most_6rad"),
        "energy_pass": d.get("energy_bound_cells_with_correction_at_most_6rad"),
        "max_direct_correction": d.get("max_direct_theta_correction_norm_upper_rad"),
        "max_correlated_correction": d.get("max_correlated_theta_correction_norm_upper_rad"),
        "max_energy_correction": d.get("max_theta_energy_correction_norm_upper_rad"),
        "next": d.get("next_obligation"),
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
