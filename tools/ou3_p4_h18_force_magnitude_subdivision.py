#!/usr/bin/env python3
"""Source-faithful H18 subdivision of the normal-Live specific-force magnitude.

The exact-node innovation anatomy shows the limiting accelerometer packet is
created by carrying one broad canonical force magnitude ``f in [5,30]`` through
both the linearized accelerometer Jacobian and the exact nonlinear residual.
This producer partitions that *source variable* into four overlapping outward
cells whose union is exactly the declared 5--30 m/s^2 theorem interval.

Each force cell is propagated through the complete 202-sample shipping H-mode
word at the same P2 source node and the same 0.80-rad Cayley entry cell.  The
shared correlation-preserving Kalman-gain enclosure is used only when the old
innovation inverse would fall back to ``S>=R``.  Fixed-pivot validated inverses
are left untouched.

This is an intermediate proof partition, not a physical-domain shrink: theorem
closure requires the union of every force cell.  No successful force cell can
promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, down, up, matrix_mul
import ou3_correlated_kalman_gain as CG
import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_correction_subdivision as SUB
import ou3_p4_h18_differential_operations as DOPS
import ou3_p4_h18_interval_ad_word as SCREEN
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
import ou3_p5_deployed_quaternion_cayley_cell as QCOMP
import ou3_p5_full_h_prefix_cells as H
import ou3_interval_ad as AD
from ou3_proof_module_state import preserve_module_bindings

DEFAULT_DOMAIN = SUB.DEFAULT_DOMAIN
SCHEMA = 1
FORCE_PIECES = 4


def _split(x: Interval, pieces: int = FORCE_PIECES) -> list[Interval]:
    width = x.hi - x.lo
    cuts = [x.lo + width * k / pieces for k in range(pieces + 1)]
    out = []
    for k in range(pieces):
        lo = x.lo if k == 0 else max(x.lo, down(cuts[k]))
        hi = x.hi if k == pieces - 1 else min(x.hi, up(cuts[k + 1]))
        out.append(Interval(lo, hi))
    return out


def _norm_upper(v) -> float:
    s = 0.0
    for x in v:
        a = x.val.abs_upper() if isinstance(x, AD.AD) else x.abs_upper()
        s = up(s + up(a * a))
    return up(math.sqrt(max(0.0, s)))


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0,
          force_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("force subdivision must not be trajectory fitted")
    if not 0 <= int(force_index) < FORCE_PIECES:
        raise ValueError("force index outside four-cell partition")

    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 source nodes invalid: {nf}")
    src = NODES.h18_source_cell(source_node_index, nodes)
    node = NODES.node(source_node_index, nodes)
    sector = SECTOR.build(path)
    sf = SECTOR.validate(sector)
    words = WORDS.build(path)
    wf = WORDS.validate(words)
    q_outer = float(sector["design_cayley_norm_upper"])
    parent = CAND._ball_box_cover(q_outer, max_box_norm_factor=1.5)[0]
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])

    original_canonical = SCREEN._canonical_vector_cells
    force_full, _mag0, _geo0 = original_canonical(domain)
    pieces = _split(force_full[2])
    fcell = pieces[int(force_index)]

    records = []
    original_accepted = SUB._accepted

    def canonical_split(d):
        force, mag, geometry = original_canonical(d)
        force = list(force)
        force[2] = fcell
        geometry = dict(geometry)
        geometry["force_magnitude_source_partition"] = [fcell.lo, fcell.hi]
        geometry["force_magnitude_partition_index"] = int(force_index)
        return force, mag, geometry

    def correlated_accepted(Pm, z, Hm, Rm, residual, *, step: int, operation: str):
        PHt, S = H._innovation(Pm, Hm, Rm)
        Sinv, backend = H._spd_inverse_enclosure(S, Rm)
        meta = None
        if backend == "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE":
            meta = CG.gain_enclosure(Pm, Hm, Rm)
            K = meta["K"]
            backend = "CORRELATED_GAIN_EQUATION_NO_S_INVERSE"
        else:
            K = matrix_mul(PHt, Sinv)
        dx_raw = H._mat_vec(K, [x.val for x in residual])
        dx = DOPS.ad_matvec_interval(K, residual)
        dtheta = [-x for x in dx[:3]]
        dnorm = _norm_upper(dtheta)
        rec = {
            "step": int(step),
            "operation": operation,
            "correction_norm_upper_rad": dnorm,
            "validated_correction_norm_limit_rad": float(QCOMP.MAX_CORRECTION_NORM),
            "inverse_backend": backend,
        }
        if meta is not None:
            rec.update({
                "gain_R_eigenvalue_lower": meta["R_eigenvalue_lower"],
                "theta_gain_residual_radius_upper": meta["row_gain_radius_upper"][:3],
                "theta_gain_psd_row_norm_upper": meta["row_psd_gain_norm_upper"][:3],
            })
            records.append(dict(rec))
        if dnorm > QCOMP.MAX_CORRECTION_NORM:
            return None, None, rec
        Pj = H._shipping_joseph(Pm, K, S, PHt)
        Pr = H._reset_covariance(Pj, dx_raw[:3])
        cp = AD.deployed_correct_cayley(z[:3], dtheta)
        out = list(z)
        out[:3] = cp
        for i in range(3, SUB.N):
            out[i] = z[i] - dx[i]
        return Pr, out, rec

    with preserve_module_bindings():
        H._source_cell = lambda: src
        SCREEN._canonical_vector_cells = canonical_split
        SUB._accepted = correlated_accepted
        try:
            row = SUB._run_child(path, domain, src, parent, samples)
        finally:
            SUB._accepted = original_accepted
            SCREEN._canonical_vector_cells = original_canonical

    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]
    # Verify the four outward pieces cover the original interval without gaps.
    if pieces[0].lo > force_full[2].lo or pieces[-1].hi < force_full[2].hi:
        failures.append("force partition does not cover declared source interval")
    for a, b in zip(pieces[:-1], pieces[1:]):
        if a.hi < b.lo:
            failures.append("force partition contains a gap")
            break

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_SOURCE_FORCE_MAGNITUDE_SUBDIVISION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "source_node_index": int(source_node_index),
        "source_node": node,
        "full_word_samples": samples,
        "force_partition_count": FORCE_PIECES,
        "force_partition_index": int(force_index),
        "declared_force_magnitude_interval_mps2": [force_full[2].lo, force_full[2].hi],
        "all_force_partition_intervals_mps2": [[x.lo, x.hi] for x in pieces],
        "force_partition_interval_mps2": [fcell.lo, fcell.hi],
        "force_partition_union_preserves_declared_domain": not failures,
        "word_completed_without_correction_range_failure": bool(row["completed"]),
        "first_failure": row.get("first_failure"),
        "maximum_correction": row.get("maximum_correction"),
        "correlated_gain_fallback_records": records,
        "source_force_subdivision_is_intermediate_proof_partition_only": True,
        "physical_state_or_source_domain_shrunk": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "take the union of all four force-magnitude cells; if any still fail at the mandatory accelerometer, combine only those failing source cells with the PSD-aware S00/S11/S01 innovation partition and the exact effective-a_w correction identity before any finer split"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_SOURCE_FORCE_MAGNITUDE_SUBDIVISION":
        f.append("wrong qualification")
    for key in ("source_generated_not_trajectory_fit", "force_partition_union_preserves_declared_domain", "source_force_subdivision_is_intermediate_proof_partition_only"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed", "declared_domain_changed", "physical_state_or_source_domain_shrunk", "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("outer_angle_rad", math.nan)) != 0.80:
        f.append("outer angle is not exactly 0.80 rad")
    if d.get("source_node_index") != 0 or d.get("full_word_samples") != 202:
        f.append("focused exact-node/full-word contract changed")
    if d.get("force_partition_count") != FORCE_PIECES:
        f.append("force partition count changed")
    idx = d.get("force_partition_index")
    if not isinstance(idx, int) or not 0 <= idx < FORCE_PIECES:
        f.append("force partition index invalid")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--force-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_node_index=a.source_node_index, force_index=a.force_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_failures": vf,
        "force_index": d["force_partition_index"],
        "force_interval": d["force_partition_interval_mps2"],
        "completed": d["word_completed_without_correction_range_failure"],
        "first_failure": d["first_failure"],
        "maximum_correction": d["maximum_correction"],
        "correlated_fallbacks": d["correlated_gain_fallback_records"],
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
