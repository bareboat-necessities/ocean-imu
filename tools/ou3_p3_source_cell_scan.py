#!/usr/bin/env python3
"""Diagnostic scan of all 800 retained P2 source cells for the P3 bottleneck.

This does not promote P3.  It separates two effects that the source-uniform
LTV certificate currently combines:

* the endpoint-static covariance ceiling produced by ``translation_upper`` on
  one exact P2 cell; and
* the process floor, either with the existing source-global sigma/tau floor or
  with the same endpoint cell held static.

The first scan is still conservative in the process floor but not source-path
complete because the covariance ceiling is endpoint-static.  The second scan is
only a localization diagnostic.  Their purpose is to identify which source
coordinates must retain path dependence before the next certificate is built.
No replay values, theorem-gate changes, or filter changes are used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_ltv_translation_ucc_probe as LTV
import ou3_p3_tau_decay_budget as DECAY
import ou3_p4_source_node_cells as NODES
import ou3_p4_source_path_reachability as PATH
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
CANDIDATES = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5)


def _best(
    upper: list[float],
    tau_index: int,
    sigma_min: float,
    global_tau_bounds: tuple[float, float],
    decay_cache: dict[tuple[int, float], float],
) -> dict:
    rows = []
    for H in CANDIDATES:
        row = LTV.ltv_relative_process_floor(
            upper,
            H,
            global_tau_bounds[0],
            global_tau_bounds[1],
            sigma_min,
            decay_exponent_upper=decay_cache[(tau_index, H)],
        )
        row["clock_phase_decay_exponent_upper"] = decay_cache[(tau_index, H)]
        rows.append(row)
    return max(rows, key=lambda x: float(x["relative_process_floor_lower"]))


def _best_static(
    upper: list[float], tau_bounds: tuple[float, float], sigma_min: float
) -> dict:
    rows = []
    for H in CANDIDATES:
        row = LTV.ltv_relative_process_floor(
            upper,
            H,
            tau_bounds[0],
            tau_bounds[1],
            sigma_min,
            decay_exponent_upper=math.nextafter(H / tau_bounds[0], math.inf),
        )
        rows.append(row)
    return max(rows, key=lambda x: float(x["relative_process_floor_lower"]))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("source-cell scan must not be trajectory fitted")

    nodes_payload = NODES.build()
    nf = NODES.validate(nodes_payload)
    if nf:
        raise RuntimeError(f"P2 source nodes invalid: {nf}")

    sched = BASE.source_schedule()
    h = float(sched["dt_s"])
    Tpe = BASE.pos(
        domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence"
    )
    global_tau = tuple(map(float, sched["tau_applied_invariant_s"]))
    sigma_floor = float(nodes_payload["filter_sigma_floor_mps2"])

    # decay_budget() performs the exact staged/committed clock-phase DP.  It
    # depends on the endpoint tau cell and horizon, not sigma/R_S.  Cache the
    # 10 x 6 distinct problems instead of recomputing them for all 800 cells.
    tau_indices = sorted({int(n["tau_index"]) for n in nodes_payload["nodes"]})
    decay_cache: dict[tuple[int, float], float] = {}
    for tau_index in tau_indices:
        for H in CANDIDATES:
            samples = max(1, int(math.ceil(H / h)))
            budget = DECAY.decay_budget(tau_index, samples, path)
            decay_cache[(tau_index, H)] = float(budget["decay_exponent_upper"])

    rows = []
    for node in nodes_payload["nodes"]:
        tau = Interval(*map(float, node["tau_s"]))
        sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
        rs = Interval(*map(float, node["R_S_filter_std"]))
        upper, timing = BASE.translation_upper(tau, sigma, rs, Tpe, sched)
        upper = [float(x) for x in upper]
        global_floor = _best(
            upper,
            int(node["tau_index"]),
            sigma_floor,
            global_tau,
            decay_cache,
        )
        static_floor = _best_static(upper, (tau.lo, tau.hi), sigma.lo)
        rows.append(
            {
                "index": int(node["index"]),
                "tau_index": int(node["tau_index"]),
                "sigma_raw_index": int(node["sigma_raw_index"]),
                "R_S_index": int(node["R_S_index"]),
                "tau_s": [tau.lo, tau.hi],
                "sigma_filter_committed_mps2": [sigma.lo, sigma.hi],
                "R_S_filter_std": [rs.lo, rs.hi],
                "Sigma_translation_diagonal_upper": upper,
                "word_horizon_s_lower": float(timing["word_horizon_s_lower"]),
                "global_process_floor_best": global_floor,
                "static_cell_process_floor_best": static_floor,
            }
        )

    worst_global = min(
        rows,
        key=lambda x: float(
            x["global_process_floor_best"]["relative_process_floor_lower"]
        ),
    )
    worst_static = min(
        rows,
        key=lambda x: float(
            x["static_cell_process_floor_best"]["relative_process_floor_lower"]
        ),
    )
    useful_global = sum(
        float(x["global_process_floor_best"]["relative_process_floor_lower"])
        >= BASE.MIN_USEFUL_DELTA
        for x in rows
    )
    useful_static = sum(
        float(x["static_cell_process_floor_best"]["relative_process_floor_lower"])
        >= BASE.MIN_USEFUL_DELTA
        for x in rows
    )

    path_summary = PATH.build(path)
    pf = PATH.validate(path_summary)
    if pf:
        raise RuntimeError(f"retained P2 path graph invalid: {pf}")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_ALL_SOURCE_CELL_LOCALIZATION_SCAN",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "diagnostic_only": True,
        "P3_PROMOTED": False,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "source_cells_scanned": len(rows),
        "distinct_tau_decay_problems": len(decay_cache),
        "global_process_floor_useful_cells": useful_global,
        "static_cell_process_floor_useful_cells": useful_static,
        "worst_global_process_floor": worst_global,
        "worst_static_cell_process_floor": worst_static,
        "retained_P2_transition_edges": int(path_summary["transition_edges"]),
        "retained_P2_states": int(path_summary["partition"]["states"]),
        "retained_P2_complete_relation": (
            int(path_summary["transition_edges"])
            == int(path_summary["partition"]["states"]) ** 2
        ),
        "rows": rows,
        "next_obligation": (
            "retain source-path dependence for the coordinates identified by the scan; "
            "do not promote endpoint-static localization as source-complete P3"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_ALL_SOURCE_CELL_LOCALIZATION_SCAN":
        f.append("wrong qualification")
    for key in ("source_only", "diagnostic_only"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("trajectory_replay_used", "filter_changed", "P3_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("source_cells_scanned", 0)) != NODES.EXPECTED_STATES:
        f.append("did not scan all 800 source cells")
    if int(d.get("distinct_tau_decay_problems", 0)) != 10 * len(CANDIDATES):
        f.append("unexpected tau decay cache cardinality")
    for key in ("worst_global_process_floor", "worst_static_cell_process_floor"):
        rho = (
            d.get(key, {})
            .get(
                "global_process_floor_best"
                if key.startswith("worst_global")
                else "static_cell_process_floor_best",
                {},
            )
            .get("relative_process_floor_lower")
        )
        if not isinstance(rho, (int, float)) or not (
            math.isfinite(float(rho)) and float(rho) > 0.0
        ):
            f.append(f"{key} is not strict")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(
        json.dumps(
            {
                "source_cells_scanned": d["source_cells_scanned"],
                "distinct_tau_decay_problems": d["distinct_tau_decay_problems"],
                "global_process_floor_useful_cells": d[
                    "global_process_floor_useful_cells"
                ],
                "static_cell_process_floor_useful_cells": d[
                    "static_cell_process_floor_useful_cells"
                ],
                "worst_global": {
                    "index": d["worst_global_process_floor"]["index"],
                    "rho": d["worst_global_process_floor"][
                        "global_process_floor_best"
                    ]["relative_process_floor_lower"],
                },
                "worst_static": {
                    "index": d["worst_static_cell_process_floor"]["index"],
                    "rho": d["worst_static_cell_process_floor"][
                        "static_cell_process_floor_best"
                    ]["relative_process_floor_lower"],
                },
                "retained_P2_transition_edges": d["retained_P2_transition_edges"],
                "retained_P2_complete_relation": d[
                    "retained_P2_complete_relation"
                ],
                "validation_failures": vf,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
