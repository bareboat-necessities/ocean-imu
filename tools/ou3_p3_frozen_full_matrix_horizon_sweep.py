#!/usr/bin/env python3
"""Frozen-source full-matrix history-horizon sweep for OU-III P3.

This diagnostic reuses the rigorous 4x4 selected-process translation backend
from :mod:`ou3_p3_frozen_full_matrix_translation` but records the first tested
physical horizon at which the unchanged 1e-18 comparison closes.  Its purpose
is architectural: if the matrix certificate only needs a short recent suffix,
the source-complete P3 consumer can propagate the frozen P2 correlation/path
memory over that suffix instead of carrying the already-saturated full word
ancestor hull.

Every source cell is still frozen within a row, so this is NOT canonical P3 and
cannot be used to promote P3/P4/P5.  It only sizes the history that the next
source-varying correlated-matrix producer must cover.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p3_frozen_full_matrix_translation as F
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DEFAULT_HORIZONS_S = (0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00, 1.25)
X_SUBCELLS = 8


def _subcell_sweep(x: Interval, checkpoint_steps: dict[int, float],
                   R_aw_norm: float, R_S_norm: float,
                   sigma_lower: float, h: float,
                   upper: list[float], gate: float) -> dict:
    Fm = F._transition(x)
    Ft = F._transpose(Fm)
    Q = F._scaled_Q(x)
    P = F._mat_zero()
    max_step = max(checkpoint_steps)
    checkpoints: dict[str, dict] = {}
    first_closed = None
    for step in range(1, max_step + 1):
        P = F._sym(F._add(F._mul(F._mul(Fm, P), Ft), Q))
        P = F._measurement_update(P, 3, R_aw_norm)
        P = F._measurement_update(P, 2, R_S_norm)
        if step not in checkpoint_steps:
            continue
        horizon = checkpoint_steps[step]
        posterior_spd = bool(symmetric_positive_definite_ldlt(P)[0])
        gate_spd = False
        if posterior_spd:
            gate_spd = bool(symmetric_positive_definite_ldlt(
                F._gate_matrix(P, sigma_lower, h, upper, gate)
            )[0])
        checkpoints[f"{horizon:.6g}"] = {
            "step": step,
            "posterior_interval_spd": posterior_spd,
            "useful_gate_spd": gate_spd,
        }
        if gate_spd and first_closed is None:
            first_closed = horizon
    return {
        "x_h_over_tau": x.as_list(),
        "checkpoints": checkpoints,
        "first_tested_useful_horizon_s": first_closed,
    }


def build(domain_path: Path = DEFAULT_DOMAIN,
          horizons_s: tuple[float, ...] = DEFAULT_HORIZONS_S) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("frozen matrix horizon sweep must not be trajectory fitted")
    runtime = domain.get("configured_runtime", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("matrix horizon sweep requires zero lever arm")
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        raise RuntimeError("matrix horizon sweep requires transparent vibration guard")

    requested = sorted(set(float(x) for x in horizons_s))
    if not requested or any(not (math.isfinite(x) and x > 0.0) for x in requested):
        raise ValueError("positive finite horizons required")
    if requested[-1] > F.HORIZON_S + 1e-15:
        raise ValueError("horizon sweep exceeds audited frozen-matrix horizon")

    nodes_payload = NODES.build()
    nf = NODES.validate(nodes_payload)
    if nf:
        raise RuntimeError(f"source nodes invalid: {nf}")
    sched = BASE.source_schedule()
    h = float(sched["dt_s"])
    checkpoint_steps: dict[int, float] = {}
    for H in requested:
        n = int(round(H / h))
        if n <= 0 or abs(n * h - H) > 1e-12:
            raise ValueError(f"horizon {H} is not an exact deployed-sample multiple")
        checkpoint_steps[n] = H

    gate = float(BASE.MIN_USEFUL_DELTA)
    if gate != 1.0e-18:
        raise RuntimeError("history sweep requires unchanged 1e-18 useful gate")
    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    vc = BASE.VECTOR.build()["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    axis_factor = min(BASE.source_rs_axis_std_factors())

    pass_counts = {f"{H:.6g}": 0 for H in requested}
    rows = []
    worst_first = 0.0
    nonclosing = []

    for node in nodes_payload["nodes"]:
        tau = Interval(*map(float, node["tau_s"]))
        sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
        rs = Interval(*map(float, node["R_S_filter_std"]))
        upper, timing = BASE.translation_upper(tau, sigma, rs, Tpe, sched)
        upper = [float(x) for x in upper]
        if float(timing["word_horizon_s_lower"]) + 1e-15 < requested[-1]:
            raise RuntimeError("tested suffix exceeds retained P3 word horizon")

        R_aw_norm = BASE.down(ra / BASE.up(sigma.hi * sigma.hi))
        rS = BASE.down(rs.lo * axis_factor)
        denomS = BASE.up(sigma.hi * sigma.hi * h ** 6)
        R_S_norm = BASE.down(BASE.down(rS * rS) / denomS)
        if not (R_aw_norm > 0.0 and R_S_norm > 0.0):
            raise RuntimeError("normalized measurement noise lost positivity")

        xlo = BASE.down(h / tau.hi)
        xhi = BASE.up(h / tau.lo)
        subrows = [
            _subcell_sweep(
                xcell, checkpoint_steps, R_aw_norm, R_S_norm,
                sigma.lo, h, upper, gate,
            )
            for xcell in F._split_x(xlo, xhi, X_SUBCELLS)
        ]

        horizon_pass = {}
        for H in requested:
            key = f"{H:.6g}"
            ok = all(r["checkpoints"][key]["useful_gate_spd"] for r in subrows)
            horizon_pass[key] = ok
            pass_counts[key] += int(ok)
        first = next((H for H in requested if horizon_pass[f"{H:.6g}"]), None)
        if first is None:
            nonclosing.append(int(node["index"]))
        else:
            worst_first = max(worst_first, first)
        rows.append({
            "index": int(node["index"]),
            "tau_index": int(node["tau_index"]),
            "sigma_raw_index": int(node["sigma_raw_index"]),
            "R_S_index": int(node["R_S_index"]),
            "first_tested_useful_horizon_s": first,
            "horizon_pass": horizon_pass,
        })

    first_all = next(
        (H for H in requested if pass_counts[f"{H:.6g}"] == NODES.EXPECTED_STATES),
        None,
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_FROZEN_FULL_MATRIX_TRANSLATION_HISTORY_HORIZON_DIAGNOSTIC",
        "source_partition_only": True,
        "frozen_source_inside_suffix": True,
        "source_complete_time_variation_covered": False,
        "P2_correlation_interface_consumed": False,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "full_4x4_translation_matrix_retained": True,
        "determinant_e3_scalarization_used": False,
        "strongest_translation_measurements_every_sample": True,
        "useful_gate": gate,
        "dt_s": h,
        "tested_horizons_s": requested,
        "x_subdivision_target": X_SUBCELLS,
        "source_cells_scanned": len(rows),
        "pass_cells_by_horizon": pass_counts,
        "first_tested_horizon_all_800_cells_useful_s": first_all,
        "largest_first_useful_horizon_among_closing_cells_s": worst_first if rows else None,
        "nonclosing_source_cells_at_max_horizon": nonclosing,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "rows": rows,
        "next_obligation": (
            "use the first all-cell useful horizon, if one exists, to size an exact source-varying propagation over OU3_P2_CORRELATED_STAGE_TRANSFER_V1; this frozen sweep itself is not source-complete P3"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_FROZEN_FULL_MATRIX_TRANSLATION_HISTORY_HORIZON_DIAGNOSTIC":
        f.append("wrong qualification")
    for key in (
        "source_partition_only", "frozen_source_inside_suffix",
        "full_4x4_translation_matrix_retained",
        "strongest_translation_measurements_every_sample",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_complete_time_variation_covered", "P2_correlation_interface_consumed",
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "determinant_e3_scalarization_used", "P3_PROMOTED", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("useful_gate", 0.0)) != 1.0e-18:
        f.append("useful gate changed")
    if int(d.get("source_cells_scanned", 0)) != NODES.EXPECTED_STATES:
        f.append("did not scan all 800 source cells")
    horizons = d.get("tested_horizons_s", [])
    if not horizons:
        f.append("no horizons tested")
    for H in horizons:
        key = f"{float(H):.6g}"
        c = int(d.get("pass_cells_by_horizon", {}).get(key, -1))
        if not 0 <= c <= NODES.EXPECTED_STATES:
            f.append(f"invalid pass count at {H}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--horizons-s", type=float, nargs="*", default=list(DEFAULT_HORIZONS_S))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, tuple(args.horizons_s))
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "pass_cells_by_horizon": d["pass_cells_by_horizon"],
        "first_all_800_s": d["first_tested_horizon_all_800_cells_useful_s"],
        "nonclosing_at_max": len(d["nonclosing_source_cells_at_max_horizon"]),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
