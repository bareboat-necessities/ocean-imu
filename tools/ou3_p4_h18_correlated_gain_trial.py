#!/usr/bin/env python3
"""Trial the shared correlation-preserving gain on the failing H18 word.

This is a proof-backend diagnostic only.  It preserves the exact source node,
0.80-rad outer Cayley cell, shipping Joseph covariance algebra and deployed
quaternion correction limit.  Fixed-pivot innovation inverses are left
unchanged.  Only the coarse ``S>=R`` inverse fallback is replaced by the
correlated gain-equation enclosure from ``ou3_correlated_kalman_gain``.

A successful run would justify reusing the same gain enclosure in #450's full
word Jacobian and deriving Joseph S^-1 through H K = I-R S^-1 in #449.  Failure
remains useful: the reported residual/gain radii identify whether covariance
or source subdivision is still required.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import matrix_mul
import ou3_correlated_kalman_gain as CG
import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_correction_subdivision as SUB
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
import ou3_p5_deployed_quaternion_cayley_cell as QCOMP
import ou3_p5_full_h_prefix_cells as H
import ou3_p4_h18_differential_operations as DOPS
import ou3_interval_ad as AD
from ou3_proof_module_state import preserve_module_bindings

DEFAULT_DOMAIN = SUB.DEFAULT_DOMAIN
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("correlated-gain trial must not be trajectory fitted")

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
    parent = CAND._ball_box_cover(q, max_box_norm_factor=1.5)[0]
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])

    correlated_records = []
    original_accepted = SUB._accepted

    def correlated_accepted(Pm, z, Hm, Rm, residual, *, step: int, operation: str):
        PHt, S = H._innovation(Pm, Hm, Rm)
        Sinv, inverse_backend = H._spd_inverse_enclosure(S, Rm)
        gain_meta = None
        if inverse_backend == "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE":
            gain_meta = CG.gain_enclosure(Pm, Hm, Rm)
            K = gain_meta["K"]
            inverse_backend = "CORRELATED_GAIN_EQUATION_NO_S_INVERSE"
        else:
            K = matrix_mul(PHt, Sinv)

        dx_raw = H._mat_vec(K, [x.val for x in residual])
        Pj = H._shipping_joseph(Pm, K, S, PHt)
        Pr = H._reset_covariance(Pj, dx_raw[0:3])
        dx = DOPS.ad_matvec_interval(K, residual)
        d = [-x for x in dx[:3]]
        dnorm = SUB._norm_upper(d)

        rec = {
            "step": int(step),
            "operation": operation,
            "correction_norm_upper_rad": dnorm,
            "validated_correction_norm_limit_rad": float(QCOMP.MAX_CORRECTION_NORM),
            "inverse_backend": inverse_backend,
        }
        if gain_meta is not None:
            rec.update({
                "correlated_gain_R_eigenvalue_lower": gain_meta["R_eigenvalue_lower"],
                "correlated_gain_max_row_residual_norm_upper": max(gain_meta["row_residual_norm_upper"]),
                "correlated_gain_max_row_radius_upper": max(gain_meta["row_gain_radius_upper"]),
                "correlated_gain_theta_row_radii_upper": gain_meta["row_gain_radius_upper"][0:3],
                "correlated_gain_aw_row_radii_upper": gain_meta["row_gain_radius_upper"][15:18],
            })
            correlated_records.append(dict(rec))

        if dnorm > QCOMP.MAX_CORRECTION_NORM:
            return None, None, rec
        cp = AD.deployed_correct_cayley(z[:3], d)
        out = list(z)
        out[:3] = cp
        for i in range(3, SUB.N):
            out[i] = z[i] - dx[i]
        return Pr, out, rec

    with preserve_module_bindings():
        H._source_cell = lambda: src
        SUB._accepted = correlated_accepted
        try:
            row = SUB._run_child(path, domain, src, parent, samples)
        finally:
            SUB._accepted = original_accepted

    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_CORRELATED_GAIN_FULL_WORD_TRIAL",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "source_node_index": int(source_node_index),
        "source_node": node,
        "full_word_samples": samples,
        "entry_cayley_box": parent,
        "correlated_gain_attempt_count": len(correlated_records),
        "correlated_gain_records": correlated_records,
        "word_completed_without_correction_range_failure": bool(row["completed"]),
        "first_failure": row.get("first_failure"),
        "maximum_correction": row.get("maximum_correction"),
        "shipping_filter_changed": False,
        "physical_state_or_source_domain_shrunk": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "if the correlated gain closes the sample-192 correction, integrate it into the shared H18 measurement enclosure and derive Joseph S^-1 from H K = I-R S^-1 for the #449 ledger; otherwise use its residual row radii together with the innovation anatomy to subdivide only the covariance/source factors responsible for the remaining gain radius"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_CORRELATED_GAIN_FULL_WORD_TRIAL":
        f.append("wrong qualification")
    if d.get("source_generated_not_trajectory_fit") is not True:
        f.append("trial is not source generated")
    for key in (
        "source_replay_used", "filter_changed", "declared_domain_changed",
        "shipping_filter_changed", "physical_state_or_source_domain_shrunk",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("outer_angle_rad", math.nan)) != 0.80:
        f.append("outer angle is not exactly 0.80 rad")
    if d.get("source_node_index") != 0:
        f.append("focused trial source node is not zero")
    if d.get("full_word_samples") != 202:
        f.append("focused trial is not the full 202-sample word")
    if int(d.get("correlated_gain_attempt_count", 0)) <= 0:
        f.append("correlated gain was never exercised")
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
        "attempts": d["correlated_gain_attempt_count"],
        "completed": d["word_completed_without_correction_range_failure"],
        "first_failure": d["first_failure"],
        "maximum_correction": d["maximum_correction"],
        "correlated_records": d["correlated_gain_records"],
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
