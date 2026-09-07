#!/usr/bin/env python3
"""Non-promoting A21 complete-word failure-mechanism probe.

This diagnostic keeps the exact same frozen shipping word, covariance/gain
history, acceptance branches, and every actual-R_S S=0 event.  It does two
things only:

1. Schur-decomposes the *complete-word linear dissipation matrix* with respect
   to the initial a_w block.  This decides whether an a_w Schur treatment can
   materially change the limiting direction before any such proof route is
   built.
2. Replays selected subdirections of the already-observed A21 limiting
   direction through the covariance-free nonlinear shadow.  These are mechanism
   probes, not replacement theorem sources and never P4 promotion.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import struct
import subprocess
from pathlib import Path

import numpy as np

NX = 21
OFF_AW = 15
OFF_BA = 18
DONE_RE = re.compile(
    r"OU3_FROZEN_SHADOW_DONE mode=(?P<mode>H18|A21) scale=(?P<scale>[-+0-9.eE]+) "
    r"t0=(?P<t0>[-+0-9.eE]+) t1=(?P<t1>[-+0-9.eE]+) "
    r"V0=(?P<V0>[-+0-9.eE]+) V1=(?P<V1>[-+0-9.eE]+) rho=(?P<rho>[-+0-9.eE]+) "
    r"V0_phi=(?P<V0phi>[-+0-9.eE]+) V1_phi=(?P<V1phi>[-+0-9.eE]+) rho_phi=(?P<rhophi>[-+0-9.eE]+) "
    r"reconstruction_max=(?P<recon>[-+0-9.eE]+) prediction_count=(?P<pred>\d+) "
    r"S_count=(?P<S>\d+) accel_count=(?P<acc>\d+) vector_count=(?P<vector>\d+)"
)


def _read_map_record(path: Path, record: int) -> tuple[float, float, np.ndarray]:
    data = path.read_bytes()
    if data[:8] != b"OU3MAP3\0":
        raise RuntimeError("unexpected map magic")
    version, nx, _stride = struct.unpack_from("<III", data, 8)
    if version != 3 or nx != NX:
        raise RuntimeError("unexpected map schema")
    rec_bytes = 1824
    count = (len(data) - 20) // rec_bytes
    if not 0 <= record < count:
        raise RuntimeError("map record outside trace")
    off = 20 + record * rec_bytes
    t0, t1 = struct.unpack_from("<dd", data, off)
    off += 16 + 4 + 12 + 28
    M = np.frombuffer(data, dtype="<f4", count=NX*NX, offset=off).astype(float).reshape(NX, NX)
    return t0, t1, M


def _read_cov_record(path: Path, record: int) -> tuple[float, float, np.ndarray, np.ndarray]:
    data = path.read_bytes()
    if data[:8] != b"OU3COV1\0":
        raise RuntimeError("unexpected covariance magic")
    version, nx, _stride = struct.unpack_from("<III", data, 8)
    if version != 1 or nx != NX:
        raise RuntimeError("unexpected covariance schema")
    rec_bytes = 3544
    count = (len(data) - 20) // rec_bytes
    if not 0 <= record < count:
        raise RuntimeError("covariance record outside trace")
    off = 20 + record * rec_bytes
    t0, t1 = struct.unpack_from("<dd", data, off)
    off += 16
    P0 = np.frombuffer(data, dtype="<f4", count=NX*NX, offset=off).astype(float).reshape(NX, NX)
    off += NX*NX*4
    P1 = np.frombuffer(data, dtype="<f4", count=NX*NX, offset=off).astype(float).reshape(NX, NX)
    return t0, t1, 0.5*(P0+P0.T), 0.5*(P1+P1.T)


def _linear_schur(map_path: Path, cov_path: Path, worst: dict, direction: np.ndarray) -> dict:
    record = int(worst["record"])
    mt0, mt1, M = _read_map_record(map_path, record)
    ct0, ct1, P0, P1 = _read_cov_record(cov_path, record)
    for a, b, label in ((mt0, ct0, "t0"), (mt1, ct1, "t1"), (mt0, float(worst["t0"]), "selected t0"), (mt1, float(worst["t1"]), "selected t1")):
        if abs(a-b) > 2e-6:
            raise RuntimeError(f"{label} mismatch in Schur diagnostic")
    Q0 = np.linalg.inv(P0)
    Q1 = np.linalg.inv(P1)
    D = Q0 - M.T @ Q1 @ M
    D = 0.5*(D+D.T)
    aw = np.arange(OFF_AW, OFF_AW+3)
    rest = np.array([i for i in range(NX) if i not in aw])
    Daa = D[np.ix_(aw, aw)]
    Dar = D[np.ix_(aw, rest)]
    Dra = Dar.T
    Drr = D[np.ix_(rest, rest)]
    Schur = Drr - Dra @ np.linalg.solve(Daa, Dar)
    Schur = 0.5*(Schur+Schur.T)
    r = direction[rest]
    a = direction[aw]
    completion = a + np.linalg.solve(Daa, Dar @ r)
    total = float(direction @ D @ direction)
    rest_energy = float(r @ Schur @ r)
    aw_energy = float(completion @ Daa @ completion)
    return {
        "record": record,
        "word_t0": mt0,
        "word_t1": mt1,
        "full_D_lambda_min": float(np.linalg.eigvalsh(D)[0]),
        "aw_pivot_lambda_min": float(np.linalg.eigvalsh(Daa)[0]),
        "schur_rest_lambda_min": float(np.linalg.eigvalsh(Schur)[0]),
        "direction_linear_dissipation_total_at_unit_scale": total,
        "direction_linear_dissipation_schur_rest_at_unit_scale": rest_energy,
        "direction_linear_dissipation_conditional_aw_at_unit_scale": aw_energy,
        "conditional_aw_fraction_of_direction_dissipation": aw_energy / total if total > 0 else None,
        "aw_schur_materially_changes_limiting_direction": bool(total > 0 and aw_energy/total > 1e-3),
    }


def _probe_directions(full: np.ndarray) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    out["full"] = full.copy()
    x = full.copy(); x[OFF_BA:OFF_BA+3] = 0.0; out["without_ba"] = x
    x = np.zeros(NX); x[OFF_BA:OFF_BA+3] = full[OFF_BA:OFF_BA+3]; out["ba_only"] = x
    x = np.zeros(NX); x[0:3] = full[0:3]; x[OFF_BA:OFF_BA+3] = full[OFF_BA:OFF_BA+3]; out["theta_ba"] = x
    x = full.copy(); x[OFF_AW:OFF_AW+3] = 0.0; out["without_aw"] = x
    x = np.zeros(NX); x[OFF_AW:OFF_AW+3] = full[OFF_AW:OFF_AW+3]; out["aw_only"] = x
    return out


def _run(sim: Path, input_path: Path, out_dir: Path, worst: dict, name: str, direction: np.ndarray, scale: float) -> dict:
    tag = f"{scale:+.6g}".replace("+", "p").replace("-", "m")
    trace = out_dir / f"mechanism_{name}_{tag}.csv"
    env = os.environ.copy()
    env.update({
        "OU3_SHADOW_TRACE": str(trace),
        "OU3_SHADOW_T0": f"{float(worst['t0']):.17g}",
        "OU3_SHADOW_T1": f"{float(worst['t1']):.17g}",
        "OU3_SHADOW_MODE": "A21",
        "OU3_SHADOW_DIRECTION": ",".join(f"{float(x):.17g}" for x in direction),
        "OU3_SHADOW_SCALE": f"{scale:.17g}",
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": "0",
    })
    cp = subprocess.run([str(sim), "--input", str(input_path)], env=env, text=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    matches = list(DONE_RE.finditer(cp.stdout))
    result = {"probe": name, "scale": scale, "returncode": cp.returncode,
              "direction_group_norms": {
                  "theta": float(np.linalg.norm(direction[0:3])),
                  "aw": float(np.linalg.norm(direction[OFF_AW:OFF_AW+3])),
                  "ba": float(np.linalg.norm(direction[OFF_BA:OFF_BA+3])),
              },
              "stdout_tail": "\n".join(cp.stdout.splitlines()[-8:])}
    if cp.returncode != 0 or not matches:
        result["valid"] = False
        return result
    g = matches[-1].groupdict()
    result.update({
        "valid": True,
        "rho_raw": float(g["rho"]),
        "rho_phi": float(g["rhophi"]),
        "V0_raw": float(g["V0"]),
        "V1_raw": float(g["V1"]),
        "reconstruction_max": float(g["recon"]),
        "prediction_count": int(g["pred"]),
        "S_count": int(g["S"]),
        "acc_count": int(g["acc"]),
        "vector_count": int(g["vector"]),
    })
    result["event_counts_match"] = (
        result["prediction_count"] == 600
        and result["S_count"] == int(worst["S_update_count"])
        and result["acc_count"] == int(worst["acc_count"])
        and result["vector_count"] == int(worst["mag_count"])
    )
    result["valid"] = result["event_counts_match"] and result["reconstruction_max"] <= 5e-5
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--linear-json", type=Path, required=True)
    ap.add_argument("--map", type=Path, required=True)
    ap.add_argument("--cov", type=Path, required=True)
    ap.add_argument("--sim", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    linear = json.loads(args.linear_json.read_text(encoding="utf-8"))
    worst = linear["modes"]["A21"]["worst"]
    direction = np.asarray(worst["maximizing_direction"]["components"], dtype=float)
    schur = _linear_schur(args.map, args.cov, worst, direction)
    probes = _probe_directions(direction)
    scales = (4.0, 6.0, 6.5, 8.0)
    cases = []
    for name, d in probes.items():
        if not np.any(d):
            continue
        for scale in scales:
            cases.append(_run(args.sim.resolve(), args.input.resolve(), args.output_dir.resolve(), worst, name, d, -scale))
            cases.append(_run(args.sim.resolve(), args.input.resolve(), args.output_dir.resolve(), worst, name, d, +scale))
    valid = [c for c in cases if c.get("valid")]
    by_probe = {}
    for name in probes:
        xs = [c for c in valid if c["probe"] == name]
        if not xs:
            continue
        by_probe[name] = {
            "cases": xs,
            "worst_raw_rho": max(float(c["rho_raw"]) for c in xs),
            "worst_phi_rho": max(float(c["rho_phi"]) for c in xs),
            "first_raw_crossing_abs_scale": min((abs(float(c["scale"])) for c in xs if float(c["rho_raw"]) >= 1.0), default=None),
            "first_phi_crossing_abs_scale": min((abs(float(c["scale"])) for c in xs if float(c["rho_phi"]) >= 1.0), default=None),
        }
    report = {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_A21_FAILURE_MECHANISM",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "same_complete_word_retained": True,
        "same_shipping_covariance_and_gains_retained": True,
        "all_due_S_updates_with_actual_RS_retained": True,
        "source_family_replaced": False,
        "declared_domain_changed": False,
        "filter_changed": False,
        "P4_promoted": False,
        "linear_full_word_aw_schur": schur,
        "probe_scales": list(scales),
        "probes": by_probe,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "aw_schur": schur,
        "probe_summary": {k: {x: v[x] for x in ("worst_raw_rho", "worst_phi_rho", "first_raw_crossing_abs_scale", "first_phi_crossing_abs_scale")} for k, v in by_probe.items()},
    }, indent=2, sort_keys=True))
    return 0 if valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
