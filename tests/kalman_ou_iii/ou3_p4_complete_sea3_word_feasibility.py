#!/usr/bin/env python3
"""Non-promoting complete-SEA3 P4 complete-word feasibility diagnostic.

This consumes the host-only OU3MAP3 exact linear transition trace together with
OU3COV1 moving shipping covariance boundaries. Those traces come from one
genuine coupled simulation run and include the shipping prediction matrices,
every accepted accelerometer/vector correction, covariance floor/reset, and
every due S=0 update with the actual applied R_S matrix.

This file intentionally lives with the host-only simulation observers rather
than the dependency-free tools/stability certificate package. It is a point
feasibility/falsification experiment only and cannot promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path

import numpy as np

MAP_MAGIC = b"OU3MAP3\0"
COV_MAGIC = b"OU3COV1\0"
MAP_VERSION = 3
COV_VERSION = 1
EXPECTED_NX = 21
EXPECTED_HORIZON_S = 3.0


def _read_header(f, magic: bytes, version: int) -> tuple[int, int]:
    got = f.read(8)
    if got != magic:
        raise RuntimeError(f"wrong trace magic {got!r}, expected {magic!r}")
    raw = f.read(12)
    if len(raw) != 12:
        raise RuntimeError("truncated trace header")
    v, nx, stride = struct.unpack("<III", raw)
    if v != version or nx != EXPECTED_NX:
        raise RuntimeError(f"unsupported trace header version={v}, nx={nx}")
    return nx, stride


def read_maps(path: Path) -> tuple[int, list[dict]]:
    rows: list[dict] = []
    with path.open("rb") as f:
        nx, stride = _read_header(f, MAP_MAGIC, MAP_VERSION)
        prefix_fmt = "<ddIiii7f"
        prefix_size = struct.calcsize(prefix_fmt)
        matrix_size = 4 * nx * nx
        while True:
            raw = f.read(prefix_size)
            if not raw:
                break
            if len(raw) != prefix_size:
                raise RuntimeError("truncated exact-map record prefix")
            vals = struct.unpack(prefix_fmt, raw)
            mraw = f.read(matrix_size)
            if len(mraw) != matrix_size:
                raise RuntimeError("truncated exact-map matrix")
            M = np.frombuffer(mraw, dtype="<f4").astype(np.float64).reshape(nx, nx)
            rows.append({
                "t0": vals[0], "t1": vals[1], "flags": vals[2],
                "acc_count": vals[3], "mag_count": vals[4], "pseudo_count": vals[5],
                "tau0": vals[6], "sigma0": vals[7], "rs0": vals[8],
                "tau1": vals[9], "sigma1": vals[10], "rs1": vals[11],
                "linearization_residual": vals[12], "M": M,
            })
    return stride, rows


def read_covariances(path: Path) -> tuple[int, list[dict]]:
    rows: list[dict] = []
    with path.open("rb") as f:
        nx, stride = _read_header(f, COV_MAGIC, COV_VERSION)
        payload_size = 8 * nx * nx
        while True:
            raw = f.read(16)
            if not raw:
                break
            if len(raw) != 16:
                raise RuntimeError("truncated covariance record prefix")
            t0, t1 = struct.unpack("<dd", raw)
            payload = f.read(payload_size)
            if len(payload) != payload_size:
                raise RuntimeError("truncated covariance matrices")
            a = np.frombuffer(payload, dtype="<f4").astype(np.float64)
            P0 = a[: nx * nx].reshape(nx, nx)
            P1 = a[nx * nx :].reshape(nx, nx)
            rows.append({"t0": t0, "t1": t1, "P0": P0, "P1": P1})
    return stride, rows


def _sym(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.T)


def _mode(flags: int) -> str:
    start_active = bool(flags & (1 << 3))
    end_active = bool(flags & (1 << 4))
    if start_active != end_active:
        return "HYBRID"
    return "A21" if start_active else "H18"


def _linear_ratio(M: np.ndarray, P0: np.ndarray, P1: np.ndarray) -> tuple[float, np.ndarray]:
    P0 = _sym(P0)
    P1 = _sym(P1)
    L = np.linalg.cholesky(P0)
    X = np.linalg.solve(P1, M @ L)
    B = _sym(L.T @ M.T @ X)
    w, V = np.linalg.eigh(B)
    idx = int(np.argmax(w))
    rho = float(w[idx])
    x = L @ V[:, idx]
    v0 = float(x @ np.linalg.solve(P0, x))
    if not (math.isfinite(rho) and math.isfinite(v0) and v0 > 0.0):
        raise RuntimeError("non-finite generalized complete-word ratio")
    x /= math.sqrt(v0)
    return rho, x


def _direction_json(x: np.ndarray, dim: int) -> dict:
    full = np.zeros(EXPECTED_NX, dtype=np.float64)
    full[:dim] = x
    groups = {
        "theta": full[0:3], "bg": full[3:6], "v": full[6:9],
        "p": full[9:12], "S": full[12:15], "aw": full[15:18],
        "ba": full[18:21],
    }
    return {
        "information_energy_at_start": 1.0,
        "components": [float(v) for v in full],
        "csv": ",".join(f"{float(v):.17g}" for v in full),
        "group_euclidean_norms": {k: float(np.linalg.norm(v)) for k, v in groups.items()},
    }


def analyze(map_path: Path, cov_path: Path) -> dict:
    map_stride, maps = read_maps(map_path)
    cov_stride, covs = read_covariances(cov_path)
    if map_stride != cov_stride:
        raise RuntimeError(f"map/cov stride mismatch {map_stride} != {cov_stride}")
    if map_stride != 600:
        raise RuntimeError(f"complete-word feasibility requires 600 samples = 3 s, got {map_stride}")

    n = min(len(maps), len(covs))
    if n == 0:
        raise RuntimeError("empty map/covariance traces")
    rows: list[dict] = []
    rejected: dict[str, int] = {}
    for i in range(n):
        m, c = maps[i], covs[i]
        if abs(m["t0"] - c["t0"]) > 2.0e-4 or abs(m["t1"] - c["t1"]) > 2.0e-4:
            raise RuntimeError(f"map/cov block alignment lost at record {i}")
        flags = int(m["flags"])
        reasons = []
        if not (flags & 1): reasons.append("map_invalid")
        if flags & (1 << 7): reasons.append("hybrid_jump")
        if not (flags & (1 << 1)) or not (flags & (1 << 2)): reasons.append("not_live")
        mode = _mode(flags)
        if mode == "HYBRID": reasons.append("mode_change")
        horizon = float(m["t1"] - m["t0"])
        if abs(horizon - EXPECTED_HORIZON_S) > 0.02: reasons.append("not_3s")
        if int(m["acc_count"]) != map_stride: reasons.append("missing_accelerometer_update")
        if int(m["pseudo_count"]) < 20: reasons.append("missing_due_S_regularization")
        if reasons:
            for r in reasons: rejected[r] = rejected.get(r, 0) + 1
            continue
        dim = 18 if mode == "H18" else 21
        rho, direction = _linear_ratio(
            m["M"][:dim, :dim], c["P0"][:dim, :dim], c["P1"][:dim, :dim]
        )
        rows.append({
            "record": i, "mode": mode, "t0": float(m["t0"]), "t1": float(m["t1"]),
            "rho_linear": rho, "distance_to_one": 1.0 - rho,
            "acc_count": int(m["acc_count"]), "mag_count": int(m["mag_count"]),
            "S_update_count": int(m["pseudo_count"]),
            "tau_start": float(m["tau0"]), "tau_end": float(m["tau1"]),
            "sigma_start": float(m["sigma0"]), "sigma_end": float(m["sigma1"]),
            "RS_scalar_start": float(m["rs0"]), "RS_scalar_end": float(m["rs1"]),
            "exact_map_linearization_recovery_residual": float(m["linearization_residual"]),
            "maximizing_direction": _direction_json(direction, dim),
        })

    modes: dict[str, dict] = {}
    for mode in ("H18", "A21"):
        r = [x for x in rows if x["mode"] == mode]
        if not r:
            modes[mode] = {"legal_blocks": 0, "worst": None}
            continue
        worst = max(r, key=lambda x: x["rho_linear"])
        modes[mode] = {
            "legal_blocks": len(r),
            "rho_min": min(x["rho_linear"] for x in r),
            "rho_max": worst["rho_linear"],
            "worst": worst,
        }

    return {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_P4_WORD_FEASIBILITY",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "point_same_history_diagnostic_only": True,
        "P4_promoted": False,
        "source_family_materialized": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "map_stride_samples": map_stride,
        "word_horizon_s": EXPECTED_HORIZON_S,
        "shipping_prediction_matrices_consumed": True,
        "all_valid_accelerometer_updates_required": True,
        "all_due_S_updates_required": True,
        "actual_applied_RS_used_inside_each_S_gain": True,
        "RS_axis_std_factors_declared_by_complete_SEA3": [0.72, 0.72, 1.0],
        "full_shipping_covariance_metric_used": True,
        "H18_A21_separate": True,
        "H_to_A_hybrid_excluded_from_same_mode_word": True,
        "packet_count_remainder_budget_used": False,
        "selected_S_subset_used": False,
        "records_seen": n,
        "records_rejected": rejected,
        "modes": modes,
        "nonlinear_followup_required": True,
        "next_experiment": "inject each reported maximizing direction into ou3-neighborhood-sim on the identical source history and measure nonlinear rho_W over the same 3 s word",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", type=Path, required=True)
    ap.add_argument("--cov", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = analyze(args.map, args.cov)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "H18": d["modes"]["H18"],
        "A21": d["modes"]["A21"],
        "rejected": d["records_rejected"],
    }, indent=2, sort_keys=True))
    return 0 if all(d["modes"][m]["legal_blocks"] > 0 for m in ("H18", "A21")) else 2


if __name__ == "__main__":
    raise SystemExit(main())