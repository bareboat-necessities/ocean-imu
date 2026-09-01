#!/usr/bin/env python3
"""Dissect the exact-node H18 innovation enclosure at the first failed vector packet.

The current complete-word screen fails at the mandatory accelerometer near the
end of the 202-sample word because fixed-pivot interval inversion falls back to
the very broad ``S >= R`` entrywise inverse enclosure.  A shared verified
midpoint-Neumann inverse also refuses the same innovation family with q>1.

This diagnostic keeps the shipping proof map, exact P2 source node zero, the
0.80-rad outer domain, Joseph covariance update and reset algebra unchanged.  It
wraps the proof measurement cell and records the first operation that uses the
spectral fallback.  For that exact innovation family it emits:

* the interval innovation covariance S and its midpoint/radii;
* the point preconditioner C and interval E=I-CS;
* row-by-row contributions to ||E||_inf and the limiting row;
* diagnostic q values if the *same S midpoint* had uniformly scaled interval
  radii (not a theorem-domain shrink and not a proof claim);
* the attitude/a_w prior covariance blocks, PH^T blocks, residual, Kalman-gain
  rows and correction rows that feed the failed accelerometer update.

The radius-scaling ladder is diagnostic only.  It estimates how much
correlation-preserving subdivision/tightening the innovation enclosure needs; it
must never be interpreted as permission to shrink the physical 0.80-rad state
set or the shipping source domain.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_identity, matrix_mul, matrix_sub
import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_correction_subdivision as SUB
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
import ou3_verified_spd_inverse as VINV
from ou3_proof_module_state import preserve_module_bindings
import ou3_p5_full_h_prefix_cells as H

DEFAULT_DOMAIN = SUB.DEFAULT_DOMAIN
SCHEMA = 1


def _iv(x: Interval) -> list[float]:
    return [float(x.lo), float(x.hi)]


def _mat(A) -> list[list[list[float]]]:
    return [[_iv(x) for x in row] for row in A]


def _mid_rad(x: Interval) -> tuple[float, float]:
    m = float(x.lo + 0.5 * (x.hi - x.lo))
    r = math.nextafter(max(m - x.lo, x.hi - m), math.inf)
    return m, r


def _scaled_about_midpoint(S, scale: float):
    out = []
    for row in S:
        r2 = []
        for x in row:
            m, r = _mid_rad(x)
            rr = math.nextafter(float(scale) * r, math.inf)
            r2.append(Interval(math.nextafter(m - rr, -math.inf), math.nextafter(m + rr, math.inf)))
        out.append(r2)
    return out


def _preconditioned_anatomy(S) -> dict:
    C = VINV._point_preconditioner(S)
    E = matrix_sub(matrix_identity(len(S)), matrix_mul(C, S))
    abs_entries = [[x.abs_upper() for x in row] for row in E]
    row_sums = [math.nextafter(sum(row), math.inf) for row in abs_entries]
    q = max(row_sums)
    return {
        "C_point_preconditioner": _mat(C),
        "E_I_minus_C_S": _mat(E),
        "E_abs_upper": abs_entries,
        "E_row_abs_sum_upper": row_sums,
        "q_inf_upper": q,
        "limiting_row": max(range(len(row_sums)), key=row_sums.__getitem__),
    }


def _block(P, rows, cols):
    return [[P[i][j] for j in cols] for i in rows]


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("innovation anatomy must not be trajectory fitted")

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

    captures = []
    original_cell = H._measurement_cell

    def wrapped(Pm, Hm, Rm, residual):
        cell = original_cell(Pm, Hm, Rm, residual)
        if cell["inverse_backend"] == "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE":
            PHt, S = H._innovation(Pm, Hm, Rm)
            anatomy = _preconditioned_anatomy(S)
            midpoint_radii = [[dict(zip(("mid", "rad"), _mid_rad(x))) for x in row] for row in S]
            scaled = {}
            for scale in (1.0, 0.5, 0.25, 0.125):
                try:
                    scaled[str(scale)] = _preconditioned_anatomy(_scaled_about_midpoint(S, scale))["q_inf_upper"]
                except Exception as exc:
                    scaled[str(scale)] = f"{type(exc).__name__}: {exc}"
            captures.append({
                "inverse_backend": cell["inverse_backend"],
                "S": _mat(S),
                "S_midpoint_and_radius": midpoint_radii,
                "R": _mat(Rm),
                "preconditioned": anatomy,
                "uniform_S_radius_scale_q_diagnostic": scaled,
                "prior_theta_covariance_block": _mat(_block(Pm, range(0, 3), range(0, 3))),
                "prior_aw_covariance_block": _mat(_block(Pm, range(15, 18), range(15, 18))),
                "prior_theta_aw_cross_block": _mat(_block(Pm, range(0, 3), range(15, 18))),
                "PHt_theta_rows": _mat(_block(PHt, range(0, 3), range(0, 3))),
                "PHt_aw_rows": _mat(_block(PHt, range(15, 18), range(0, 3))),
                "H_theta_columns": _mat(_block(Hm, range(0, 3), range(0, 3))),
                "H_aw_columns": _mat(_block(Hm, range(0, 3), range(15, 18))),
                "residual": [_iv(x) for x in residual],
                "K_theta_rows": _mat(_block(cell["K"], range(0, 3), range(0, 3))),
                "K_aw_rows": _mat(_block(cell["K"], range(15, 18), range(0, 3))),
                "dx_theta": [_iv(x) for x in cell["dx"][:3]],
                "dx_aw": [_iv(x) for x in cell["dx"][15:18]],
            })
        return cell

    with preserve_module_bindings():
        H._source_cell = lambda: src
        H._measurement_cell = wrapped
        row = SUB._run_child(path, domain, src, parent, samples)

    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]
    if not captures:
        failures.append("no spectral-fallback innovation cell was captured")
    first = captures[0] if captures else None
    if first is not None and not float(first["preconditioned"]["q_inf_upper"]) > 1.0:
        failures.append("captured fallback cell unexpectedly has q<=1")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_EXACT_NODE_INNOVATION_ENCLOSURE_ANATOMY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "source_node_index": int(source_node_index),
        "full_word_samples": samples,
        "word_completed_without_correction_range_failure": bool(row.get("completed")),
        "first_word_failure": row.get("first_failure"),
        "spectral_fallback_capture_count": len(captures),
        "first_spectral_fallback": first,
        "uniform_S_radius_scaling_is_diagnostic_only": True,
        "physical_state_or_source_domain_shrunk": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "use the limiting E row and S/P/PHt block widths to build a correlation-preserving interval solve or source/covariance subdivision; require a verified innovation solve on the union of all subdivisions and keep the original 0.80-rad physical domain"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in ("source_generated_not_trajectory_fit", "uniform_S_radius_scaling_is_diagnostic_only"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed", "declared_domain_changed", "physical_state_or_source_domain_shrunk", "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("source_node_index") != 0 or d.get("full_word_samples") != 202:
        f.append("focused exact-node/full-word contract changed")
    if float(d.get("outer_angle_rad", 0.0)) != 0.80:
        f.append("outer angle changed")
    first = d.get("first_spectral_fallback")
    if not isinstance(first, dict):
        f.append("missing first spectral fallback anatomy")
    else:
        q = first.get("preconditioned", {}).get("q_inf_upper")
        if isinstance(q, bool) or not isinstance(q, (int, float)) or not math.isfinite(float(q)) or float(q) <= 1.0:
            f.append("fallback anatomy does not reproduce q>1 obstruction")
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
    first = d.get("first_spectral_fallback") or {}
    pre = first.get("preconditioned", {})
    print(json.dumps({
        "validation_pass": not vf,
        "first_word_failure": d["first_word_failure"],
        "capture_count": d["spectral_fallback_capture_count"],
        "q_inf_upper": pre.get("q_inf_upper"),
        "limiting_row": pre.get("limiting_row"),
        "E_row_abs_sum_upper": pre.get("E_row_abs_sum_upper"),
        "uniform_S_radius_scale_q_diagnostic": first.get("uniform_S_radius_scale_q_diagnostic"),
        "S_midpoint_and_radius": first.get("S_midpoint_and_radius"),
        "prior_theta_covariance_block": first.get("prior_theta_covariance_block"),
        "prior_aw_covariance_block": first.get("prior_aw_covariance_block"),
        "prior_theta_aw_cross_block": first.get("prior_theta_aw_cross_block"),
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
