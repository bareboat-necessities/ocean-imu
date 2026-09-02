#!/usr/bin/env python3
"""Full-matrix translation transfer for one frozen P2-correlation segment.

This is the numerical primitive used by the source-varying P3 consumer.  It is
bound to ``OU3_P2_CORRELATED_STAGE_TRANSFER_V1`` and deliberately uses a fixed
physical scaling

    z = [v/h, p/h^2, S/h^3, a_w],

rather than the frozen diagnostic's sigma-dependent scaling.  Therefore a
source transition changes the OU process intensity through ``sigma_aw`` but
never causes an artificial state/metric rescaling.

For one certified V1 segment the committed ``tau``, ``sigma_aw`` and ``R_S``
come from the same physical P2 cell and remain applied for the segment.  The
one-step covariance map is outward-enclosed as

    P+ = F(tau) P F(tau)^T + sigma_aw^2 Qbar(h/tau),

then conditioned by the strongest translation measurements that can occur:
accepted accelerometer and S=0 packets every sample.  These extra accepted
measurements can only reduce the selected-process covariance, so the result is
a lower covariance for every shipping translation measurement branch.

This module proves/exports one-segment transfer semantics only.  It does not
flatten source histories and cannot promote P3/P4/P5 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_frozen_full_matrix_translation as FROZEN
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
X_SUBCELLS = 4


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _scale(A, s: Interval):
    return [[s * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _physical_measurement_variances(node: dict, h: float) -> tuple[float, float]:
    vc = BASE.VECTOR.build()["configured_measurement_bounds"]
    R_aw = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rs = Interval(*map(float, node["R_S_filter_std"]))
    axis_factor = min(BASE.source_rs_axis_std_factors())
    rS = BASE.down(rs.lo * axis_factor)
    # z_S=S/h^3, hence H_z=[0,0,h^3,0] for the physical S measurement.
    # Equivalently use coordinate z_S with R_z=R_S/h^6.
    R_S_z = BASE.down(BASE.down(rS * rS) / BASE.up(h ** 6))
    if not (R_aw > 0.0 and R_S_z > 0.0):
        raise RuntimeError("physical-scaled measurement variance lost positivity")
    return R_aw, R_S_z


def _node_x_subcells(node: dict, h: float, count: int = X_SUBCELLS):
    tau = Interval(*map(float, node["tau_s"]))
    xlo = BASE.down(h / tau.hi)
    xhi = BASE.up(h / tau.lo)
    return FROZEN._split_x(xlo, xhi, count)


def one_step(P, node: dict, x: Interval, h: float):
    """Outward enclosure of one strongest-measurement selected-process step."""
    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    Fm = FROZEN._transition(x)
    Ft = FROZEN._transpose(Fm)
    # Frozen normalized Q is in y=[v/(sigma h),...,a/sigma].  In the fixed
    # z=[v/h,...,a] scaling the same matrix is multiplied by sigma^2.
    Qz = _scale(FROZEN._scaled_Q(x), sigma.square())
    out = FROZEN._sym(FROZEN._add(FROZEN._mul(FROZEN._mul(Fm, P), Ft), Qz))
    R_aw, R_S_z = _physical_measurement_variances(node, h)
    out = FROZEN._measurement_update(out, 3, R_aw)
    out = FROZEN._measurement_update(out, 2, R_S_z)
    return FROZEN._sym(out)


def propagate_subcell(P, node: dict, samples: int, x: Interval, h: float):
    out = [[P[i][j] for j in range(4)] for i in range(4)]
    for _ in range(int(samples)):
        out = one_step(out, node, x, h)
    return out


def segment_images(P, source_node: int, gap: int,
                   runtime: dict | None = None,
                   domain_path: Path = DEFAULT_DOMAIN,
                   x_subcells: int = X_SUBCELLS):
    """Return all x-subcell enclosures for one V1 applied-source segment."""
    rt = CORR.runtime(Path(domain_path).resolve()) if runtime is None else runtime
    s = int(source_node)
    g = int(gap)
    if not 0 <= s < len(rt["nodes"]):
        raise IndexError("source node outside P2 correlation partition")
    if g not in rt["gaps"]:
        raise ValueError("gap outside certified P2 correlation alphabet")
    node = rt["nodes"][s]
    h = float(rt["clock"]["dt_binary32_s"])
    return [
        {
            "x_h_over_tau": x.as_list(),
            "posterior": propagate_subcell(P, node, g, x, h),
        }
        for x in _node_x_subcells(node, h, x_subcells)
    ]


def build(domain_path: Path = DEFAULT_DOMAIN,
          representative_nodes=(0, 137, 729, 799),
          representative_gaps=(13, 21, 26)) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("correlated segment transfer must not be trajectory fitted")
    rt = CORR.runtime(path)
    corr = CORR.build(path)
    cf = CORR.validate(corr)
    if cf:
        raise RuntimeError(f"P2 correlation interface failed: {cf}")
    if corr.get("interface_version") != CORR.INTERFACE_VERSION:
        raise RuntimeError("P2 correlation interface version mismatch")

    zero = FROZEN._mat_zero()
    rows = []
    for s in map(int, representative_nodes):
        node = rt["nodes"][s]
        for g in map(int, representative_gaps):
            images = segment_images(zero, s, g, rt, path)
            rows.append({
                "source_node": s,
                "gap_samples": g,
                "tau_s": node["tau_s"],
                "sigma_filter_committed_mps2": node["sigma_filter_committed_mps2"],
                "R_S_filter_std": node["R_S_filter_std"],
                "same_source_cell_for_tau_sigma_R_S": True,
                "x_subcells": len(images),
                "all_outputs_interval_spd": all(
                    symmetric_positive_definite_ldlt(r["posterior"])[0]
                    for r in images
                ),
            })

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_FULL_MATRIX_TRANSLATION_P2_CORRELATED_SEGMENT_TRANSFER",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "fixed_physical_scaling": "z=[v/h,p/h^2,S/h^3,a_w]",
        "sigma_dependent_state_rescaling_used": False,
        "tau_sigma_R_S_same_source_cell_per_segment": True,
        "full_4x4_translation_matrix_retained": True,
        "strongest_translation_measurements_every_sample": True,
        "gap_alphabet_samples": list(rt["gaps"]),
        "representative_rows": rows,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "next_obligation": (
            "compose this transfer along legal V1 pair-shift histories over the matrix history horizon and compare the resulting selected-process lower to a covariance upper obtained from the same source history"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_FULL_MATRIX_TRANSLATION_P2_CORRELATED_SEGMENT_TRANSFER":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed",
        "tau_sigma_R_S_same_source_cell_per_segment",
        "full_4x4_translation_matrix_retained",
        "strongest_translation_measurements_every_sample",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "sigma_dependent_state_rescaling_used", "P3_PROMOTED", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("segment transfer is not bound to frozen P2 correlation version")
    if d.get("gap_alphabet_samples") != list(range(13, 27)):
        f.append("segment transfer lost 13..26 gap alphabet")
    rows = d.get("representative_rows", [])
    if not rows:
        f.append("no representative segment transfers emitted")
    for row in rows:
        if row.get("same_source_cell_for_tau_sigma_R_S") is not True:
            f.append("representative row lost source-cell correlation")
        if int(row.get("x_subcells", 0)) <= 0:
            f.append("representative row has no x subdivisions")
        if row.get("all_outputs_interval_spd") is not True:
            f.append("representative zero-start segment lost interval SPD")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True, default=lambda x: x.as_list()), encoding="utf-8")
    print(json.dumps({
        "P2_correlation_interface_version": d["P2_correlation_interface_version"],
        "representative_segments": len(d["representative_rows"]),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
