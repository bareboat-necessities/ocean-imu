#!/usr/bin/env python3
"""Non-promoting multiword complete-SEA3 A21 P4 feasibility diagnostic.

This does not create a new source family.  It concatenates consecutive legal
3-second COMPLETE_SEA3_NORMAL_LIVE_WORD blocks from the same shipping history,
composes their exact 21x21 maps, and evaluates one 6 s / 9 s finite window.
Every prediction, valid accelerometer update, asynchronous vector event,
covariance floor, and every due S=0 update with its actual applied R_S remains
inside the concatenated word.

The paper's Live theorem is finite-window and explicitly permits a longer proof
window for the S regularizer.  This point experiment only asks whether a longer
complete window is a viable nonlinear P4 target.  It cannot promote P4 or stand
in for universal SEA3 coverage.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np

import ou3_p4_complete_sea3_word_feasibility as BASE
import ou3_p4_frozen_gain_shadow_feasibility as SHADOW

DONE_RE = SHADOW.DONE_RE
NX = 21


def _legal_a21(m: dict, stride: int) -> bool:
    flags = int(m["flags"])
    return bool(
        (flags & 1)
        and not (flags & (1 << 7))
        and (flags & (1 << 1))
        and (flags & (1 << 2))
        and BASE._mode(flags) == "A21"
        and int(m["acc_count"]) == stride
        and int(m["pseudo_count"]) >= 20
    )


def _compose(maps: list[dict], start: int, blocks: int) -> np.ndarray:
    M = np.eye(NX, dtype=np.float64)
    for j in range(start, start + blocks):
        M = maps[j]["M"] @ M
    return M


def _find_worst(maps: list[dict], covs: list[dict], stride: int, blocks: int) -> dict:
    rows = []
    for i in range(0, min(len(maps), len(covs)) - blocks + 1):
        if not all(_legal_a21(maps[j], stride) for j in range(i, i + blocks)):
            continue
        if not all(abs(float(maps[j]["t1"]) - float(maps[j+1]["t0"])) <= 2e-4
                   for j in range(i, i + blocks - 1)):
            continue
        M = _compose(maps, i, blocks)
        P0 = covs[i]["P0"]
        P1 = covs[i + blocks - 1]["P1"]
        rho, direction = BASE._linear_ratio(M, P0, P1)
        rows.append({
            "start_record": i,
            "end_record": i + blocks - 1,
            "blocks": blocks,
            "horizon_s": float(maps[i + blocks - 1]["t1"] - maps[i]["t0"]),
            "t0": float(maps[i]["t0"]),
            "t1": float(maps[i + blocks - 1]["t1"]),
            "rho_linear": rho,
            "distance_to_one": 1.0 - rho,
            "acc_count": sum(int(maps[j]["acc_count"]) for j in range(i, i + blocks)),
            "mag_count": sum(int(maps[j]["mag_count"]) for j in range(i, i + blocks)),
            "S_update_count": sum(int(maps[j]["pseudo_count"]) for j in range(i, i + blocks)),
            "RS_scalar_start": float(maps[i]["rs0"]),
            "RS_scalar_end": float(maps[i + blocks - 1]["rs1"]),
            "maximizing_direction": BASE._direction_json(direction, NX),
        })
    if not rows:
        raise RuntimeError(f"no legal {3*blocks}s consecutive A21 complete word")
    return max(rows, key=lambda x: x["rho_linear"]) | {"legal_windows": len(rows)}


def _scales(limit: float) -> list[float]:
    base = [0.125, 0.5, 1.0, 2.0, 4.0, 6.0, 6.5, 8.0, 12.0, 16.0, 24.0, 32.0]
    out = [x for x in base if x <= limit * (1.0 + 1e-12)]
    if not out or limit > out[-1] * (1.0 + 1e-6):
        out.append(limit)
    return out


def _run_case(sim: Path, input_path: Path, out_dir: Path, word: dict,
              direction: list[float], scale: float) -> dict:
    tag = f"{word['blocks']}b_{scale:+.8g}".replace("+", "p").replace("-", "m")
    trace = out_dir / f"multiword_{tag}.csv"
    env = os.environ.copy()
    env.update({
        "OU3_SHADOW_TRACE": str(trace),
        "OU3_SHADOW_T0": f"{float(word['t0']):.17g}",
        "OU3_SHADOW_T1": f"{float(word['t1']):.17g}",
        "OU3_SHADOW_MODE": "A21",
        "OU3_SHADOW_DIRECTION": ",".join(f"{float(x):.17g}" for x in direction),
        "OU3_SHADOW_SCALE": f"{scale:.17g}",
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": "0",
    })
    cp = subprocess.run([str(sim), "--input", str(input_path)], env=env, text=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    matches = list(DONE_RE.finditer(cp.stdout))
    result = {"scale": scale, "absolute_scale": abs(scale), "returncode": cp.returncode,
              "stdout_tail": "\n".join(cp.stdout.splitlines()[-8:]), "trace": str(trace)}
    if cp.returncode != 0 or not matches:
        result["valid"] = False
        return result
    g = matches[-1].groupdict()
    result.update({
        "valid": True,
        "rho_raw": float(g["rho"]),
        "rho_phi": float(g["rhophi"]),
        "V0_raw": float(g["V0"]), "V1_raw": float(g["V1"]),
        "V0_phi": float(g["V0phi"]), "V1_phi": float(g["V1phi"]),
        "reconstruction_max": float(g["recon"]),
        "prediction_count": int(g["pred"]), "S_count": int(g["S"]),
        "acc_count": int(g["acc"]), "vector_count": int(g["vector"]),
        "endpoint_time_s": float(g["t1"]),
    })
    result["event_counts_match"] = (
        result["prediction_count"] == int(word["acc_count"])
        and result["acc_count"] == int(word["acc_count"])
        and result["S_count"] == int(word["S_update_count"])
        and result["vector_count"] == int(word["mag_count"])
    )
    result["endpoint_matches"] = abs(result["endpoint_time_s"] - float(word["t1"])) <= 2e-6
    result["valid"] = bool(result["event_counts_match"] and result["endpoint_matches"]
                           and result["reconstruction_max"] <= 5e-5)
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", type=Path, required=True)
    ap.add_argument("--cov", type=Path, required=True)
    ap.add_argument("--domain", type=Path, required=True)
    ap.add_argument("--sim", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stride, maps = BASE.read_maps(args.map)
    cstride, covs = BASE.read_covariances(args.cov)
    if stride != 600 or cstride != stride:
        raise RuntimeError("multiword diagnostic requires the retained 600-sample base word")
    domain = json.loads(args.domain.read_text(encoding="utf-8"))

    report = {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_P4_MULTIWORD_A21_FEASIBILITY",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "concatenates_same_source_words_only": True,
        "source_family_replaced": False,
        "selected_S_subset_used": False,
        "packet_count_remainder_budget_used": False,
        "all_due_S_updates_with_actual_RS_retained": True,
        "all_valid_accelerometer_updates_retained": True,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P4_promoted": False,
        "windows": {},
    }
    for blocks in (2, 3):
        word = _find_worst(maps, covs, stride, blocks)
        direction = [float(x) for x in word["maximizing_direction"]["components"]]
        limit, limit_detail = SHADOW.max_scale(direction, "A21", domain, args.input)
        cases = []
        for mag in _scales(limit):
            cases.append(_run_case(args.sim.resolve(), args.input.resolve(), args.output_dir.resolve(), word, direction, -mag))
            cases.append(_run_case(args.sim.resolve(), args.input.resolve(), args.output_dir.resolve(), word, direction, +mag))
        valid = [c for c in cases if c.get("valid")]
        reliable = [c for c in valid if abs(float(c["scale"])) >= SHADOW.FLOAT_RESOLUTION_SCALE_FLOOR]
        report["windows"][f"{3*blocks}s"] = {
            **word,
            "declared_scale_limit": limit,
            "declared_scale_limit_detail": limit_detail,
            "cases": cases,
            "valid_cases": len(valid),
            "worst_reliable_raw_rho": max((float(c["rho_raw"]) for c in reliable), default=None),
            "worst_reliable_phi_rho": max((float(c["rho_phi"]) for c in reliable), default=None),
            "first_reliable_raw_crossing_abs_scale": min((abs(float(c["scale"])) for c in reliable if float(c["rho_raw"]) >= 1.0), default=None),
            "first_reliable_phi_crossing_abs_scale": min((abs(float(c["scale"])) for c in reliable if float(c["rho_phi"]) >= 1.0), default=None),
            "strict_raw_contraction_on_all_tested_reliable_cases": bool(reliable) and all(float(c["rho_raw"]) < 1.0 for c in reliable),
            "strict_phi_contraction_on_all_tested_reliable_cases": bool(reliable) and all(float(c["rho_phi"]) < 1.0 for c in reliable),
            "max_nominal_reconstruction_error": max((float(c["reconstruction_max"]) for c in valid), default=None),
        }

    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: {
        "rho_linear": v["rho_linear"],
        "records": [v["start_record"], v["end_record"]],
        "counts": {"acc": v["acc_count"], "S": v["S_update_count"], "vector": v["mag_count"]},
        "worst_raw": v["worst_reliable_raw_rho"],
        "worst_phi": v["worst_reliable_phi_rho"],
        "raw_cross": v["first_reliable_raw_crossing_abs_scale"],
        "phi_cross": v["first_reliable_phi_crossing_abs_scale"],
        "scale_limit": v["declared_scale_limit"],
    } for k, v in report["windows"].items()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
