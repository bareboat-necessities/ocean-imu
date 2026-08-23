#!/usr/bin/env python3
"""Deterministic nonzero-neighborhood diagnostic for adaptive OU-III.

This stage drives ``ou3-neighborhood-sim``.  The simulator runs a nominal and a
perturbed copy of the unchanged adaptive filter under *identical* noisy sensor
samples.  At a source point the perturbed copy receives one exact MEKF error
retraction/additive-state perturbation.  The pair is then propagated over one
certified information word.

Perturbations are normalized in the local Kalman information metric: for a raw
direction ``v`` and nominal covariance ``Sigma`` we choose ``d`` so that

    d^T Sigma^-1 d = target_W.

This removes the unit/conditioning ambiguity between attitude, biases, v, p, S
and a_w.  The resulting sampled search is a diagnostic only.  It may locate an
active nonlinear/source-guard constraint, but it can never promote the
neighborhood or deployment theorem without outward-rounded enclosure.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fnmatch
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np

import ou3_numerical_certificate as BASE
import ou3_information_certificate as INFO

DEFAULT_RESULTS = BASE.DEFAULT_OUT / "neighborhood_diagnostic"
DEFAULT_DIAGNOSTICS = BASE.REPO / "reports" / "diagnostics" / "ou3_neighborhood"
DEFAULT_SIM = BASE.TEST_DIR / "ou3-neighborhood-sim"
DEFAULT_HELD_TIME_S = 60.0
DEFAULT_ACTIVE_TIME_S = 300.0
DEFAULT_TARGET_W = (0.05,)

# Compact deterministic basis for first-pass constraint discovery.  A full
# coordinate basis is available with --full-basis.
COMPACT_H = (0, 1, 2, 3, 4, 5, 6, 9, 12, 15)
COMPACT_A = COMPACT_H + (18,)
BLOCK_NAME = {
    0: "theta_x", 1: "theta_y", 2: "theta_z",
    3: "bg_x", 4: "bg_y", 5: "bg_z",
    6: "v_x", 7: "v_y", 8: "v_z",
    9: "p_x", 10: "p_y", 11: "p_z",
    12: "S_x", 13: "S_y", 14: "S_z",
    15: "aw_x", 16: "aw_y", 17: "aw_z",
    18: "ba_x", 19: "ba_y", 20: "ba_z",
}


def record_slug(family: str, hs: float) -> str:
    return f"{family.lower().replace('-', '_')}_{hs:.2f}".replace(".", "_")


def record_index() -> dict[str, tuple[str, float, str]]:
    out = {}
    for family, hs, name in BASE.RECORDS:
        slug = record_slug(family, hs)
        out[slug] = (family, hs, name)
        out[Path(name).stem] = (family, hs, name)
    return out


def nearest_covariance(map_path: Path, t: float, dim: int) -> tuple[np.ndarray, dict]:
    maps, covs, meta = INFO.pair_map_covariance(map_path, map_path.stem.replace("_exact_maps", ""))
    candidates = []
    for i, (m, c) in enumerate(zip(maps, covs)):
        candidates.append((abs(m.t0 - t), i, "start", c.start))
        candidates.append((abs(m.t1 - t), i, "end", c.end))
    if not candidates:
        raise RuntimeError(f"no covariance blocks in {map_path}")
    dt, i, side, P = min(candidates, key=lambda x: x[0])
    if dt > 0.02:
        raise RuntimeError(f"no covariance endpoint within 20 ms of t={t}: nearest {dt}")
    P = np.asarray(P[:dim, :dim], float)
    P = 0.5 * (P + P.T)
    eig = np.linalg.eigvalsh(P)
    if not np.all(np.isfinite(eig)) or float(np.min(eig)) <= 0.0:
        raise RuntimeError(f"non-SPD covariance at {map_path} block {i} {side}")
    return P, {
        "block": i,
        "side": side,
        "time_distance_s": float(dt),
        "lambda_min": float(np.min(eig)),
        "lambda_max": float(np.max(eig)),
    }


def information_normalize(P: np.ndarray, index: int, target_W: float, sign: int) -> np.ndarray:
    if not (target_W > 0.0 and math.isfinite(target_W)):
        raise ValueError("target_W must be finite positive")
    dim = P.shape[0]
    if not (0 <= index < dim):
        raise ValueError(f"coordinate {index} outside dimension {dim}")
    v = np.zeros(dim)
    v[index] = float(sign)
    w = float(v @ np.linalg.solve(P, v))
    if not (w > 0.0 and math.isfinite(w)):
        raise RuntimeError(f"invalid information norm for coordinate {index}: {w}")
    return v * math.sqrt(target_W / w)


def to_21(delta: np.ndarray) -> np.ndarray:
    out = np.zeros(21)
    out[: len(delta)] = delta
    return out


def parse_trace(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {"status": "NO_TRACE", "pass_sampled": False}
    a = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
    if a.shape == ():
        a = np.asarray([a], dtype=a.dtype)
    if len(a) == 0:
        return {"status": "NO_ROWS", "pass_sampled": False}

    endpoint = np.flatnonzero(np.asarray(a["endpoint"], int) == 1)
    if len(endpoint) == 0:
        return {
            "status": "NO_ENDPOINT",
            "pass_sampled": False,
            "row_count": int(len(a)),
            "source_match_all": bool(np.all(np.asarray(a["source_match"], int) == 1)),
        }

    i1 = int(endpoint[-1])
    W = np.asarray(a["W_nominal"], float)
    theta = np.asarray(a["theta_rad"], float)
    covrel = np.asarray(a["covariance_rel_fro"], float)
    W0 = float(W[0])
    W1 = float(W[i1])
    ratio = W1 / W0 if W0 > 0.0 else math.inf
    source_match = bool(np.all(np.asarray(a["source_match"][: i1 + 1], int) == 1))
    acceptance_match = bool(
        np.all(np.asarray(a["acc_accept_match"][: i1 + 1], int) == 1)
        and np.all(np.asarray(a["mag_accept_match"][: i1 + 1], int) == 1)
    )
    finite = bool(
        np.all(np.isfinite(W[: i1 + 1]))
        and np.all(np.isfinite(theta[: i1 + 1]))
        and np.all(np.isfinite(covrel[: i1 + 1]))
    )
    theta_max = float(np.max(theta[: i1 + 1])) if finite else math.inf
    prefix_W_gain = float(np.max(W[: i1 + 1]) / W0) if finite and W0 > 0.0 else math.inf
    decrement = 1.0 - ratio
    passed = bool(
        finite and source_match and acceptance_match and theta_max < math.pi
        and W0 > 0.0 and W1 < W0
    )
    return {
        "status": "PASS_SAMPLED" if passed else "FAIL_SAMPLED",
        "pass_sampled": passed,
        "row_count": int(i1 + 1),
        "source_match_all": source_match,
        "measurement_acceptance_match_all": acceptance_match,
        "W0": W0,
        "W1": W1,
        "endpoint_ratio_W1_over_W0": ratio,
        "relative_decrement": decrement,
        "prefix_W_gain_max": prefix_W_gain,
        "theta_max_rad": theta_max,
        "theta_max_deg": math.degrees(theta_max) if math.isfinite(theta_max) else None,
        "covariance_relative_difference_max": float(np.max(covrel[: i1 + 1])) if finite else None,
        "actual_injection_time_s": float(a["time_s"][0]),
        "actual_endpoint_time_s": float(a["time_s"][i1]),
        "qualification": "SAMPLED_PAIRWISE_NONLINEAR_DIAGNOSTIC_ONLY",
    }


def select_records(patterns: list[str]) -> list[tuple[str, float, str]]:
    rows = list(BASE.RECORDS)
    if not patterns:
        return rows
    selected = []
    for family, hs, name in rows:
        slug = record_slug(family, hs)
        if any(fnmatch.fnmatch(slug, p) or fnmatch.fnmatch(Path(name).stem, p) for p in patterns):
            selected.append((family, hs, name))
    if not selected:
        raise ValueError(f"record patterns matched nothing: {patterns}")
    return selected


def run_case(sim: Path, data: Path, trace: Path, log: Path,
             mode: str, inject_s: float, horizon_s: float,
             delta21: np.ndarray) -> tuple[int, str]:
    trace.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({
        "OU3_NEIGHBOR_TRACE": str(trace.resolve()),
        "OU3_NEIGHBOR_INJECT_TIME_S": f"{inject_s:.9g}",
        "OU3_NEIGHBOR_HORIZON_S": f"{horizon_s:.9g}",
        "OU3_NEIGHBOR_MODE": mode,
        "OU3_NEIGHBOR_DELTA": ",".join(f"{x:.17g}" for x in delta21),
        "OU3_NEIGHBOR_TRACE_STRIDE": "50",
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": "0",
    })
    p = subprocess.run(
        [str(sim.resolve()), "--input", str(data.resolve())],
        cwd=data.parent,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    log.write_text(p.stdout)
    return p.returncode, p.stdout


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=BASE.DEFAULT_OUT)
    ap.add_argument("--data-dir", type=Path, default=BASE.DEFAULT_DATA_DIR)
    ap.add_argument("--sim", type=Path, default=DEFAULT_SIM)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--diagnostic-dir", type=Path, default=DEFAULT_DIAGNOSTICS)
    ap.add_argument("--records", action="append", default=[],
                    help="glob against certificate slug/source stem; repeatable")
    ap.add_argument("--modes", default="H,A")
    ap.add_argument("--target-W", default=",".join(str(x) for x in DEFAULT_TARGET_W))
    ap.add_argument("--held-time-s", type=float, default=DEFAULT_HELD_TIME_S)
    ap.add_argument("--active-time-s", type=float, default=DEFAULT_ACTIVE_TIME_S)
    ap.add_argument("--full-basis", action="store_true")
    ap.add_argument("--positive-only", action="store_true")
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel simulator subprocesses; result ordering remains deterministic")
    ap.add_argument("--no-build", action="store_true")
    args = ap.parse_args()

    cert = args.certificate_dir.resolve()
    data_dir = args.data_dir.resolve()
    sim = args.sim.resolve()
    out = args.output_dir.resolve()
    diag = args.diagnostic_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    diag.mkdir(parents=True, exist_ok=True)

    if args.jobs < 1:
        raise ValueError("--jobs must be >= 1")
    if not args.no_build:
        subprocess.run(["make", "-C", str(BASE.TEST_DIR), sim.name], check=True)
    if not sim.exists():
        raise FileNotFoundError(sim)

    contract = json.loads((cert / "information_enclosure_contract.json").read_text())
    modes = [x.strip() for x in args.modes.split(",") if x.strip()]
    if any(m not in ("H", "A") for m in modes):
        raise ValueError("--modes accepts H,A")
    targets = [float(x) for x in args.target_W.split(",") if x.strip()]
    if any(not (x > 0.0 and math.isfinite(x)) for x in targets):
        raise ValueError("--target-W entries must be finite positive")
    signs = (1,) if args.positive_only else (-1, 1)

    tasks = []
    for family, hs, name in select_records(args.records):
        data = (data_dir / name).resolve()
        if not data.exists():
            raise FileNotFoundError(data)
        stem = Path(name).stem
        map_path = cert / f"{stem}_exact_maps.bin"
        if not map_path.exists():
            raise FileNotFoundError(map_path)

        for mode in modes:
            dim = 18 if mode == "H" else 21
            inject_s = args.held_time_s if mode == "H" else args.active_time_s
            horizon_s = float(contract["modes"][mode]["recommended_word_horizon_s"])
            P, cov_meta = nearest_covariance(map_path, inject_s, dim)
            indices = tuple(range(dim)) if args.full_basis else (COMPACT_H if mode == "H" else COMPACT_A)
            for target in targets:
                for index in indices:
                    for sign in signs:
                        d = information_normalize(P, index, target, sign)
                        d21 = to_21(d)
                        case_id = (
                            f"{record_slug(family, hs)}_{mode}_{BLOCK_NAME[index]}_"
                            f"{'p' if sign > 0 else 'm'}_W{target:.6g}"
                        ).replace(".", "_")
                        tasks.append({
                            "case": case_id,
                            "family": family,
                            "Hs_m": hs,
                            "source_file": name,
                            "data": data,
                            "mode": mode,
                            "coordinate": index,
                            "direction": BLOCK_NAME[index],
                            "sign": sign,
                            "target_W": target,
                            "requested_injection_time_s": inject_s,
                            "word_horizon_s": horizon_s,
                            "delta_21": d21,
                            "covariance_reference": cov_meta,
                            "trace": diag / f"{case_id}.csv",
                            "log": diag / f"{case_id}.log",
                        })

    def execute(task: dict) -> dict:
        rc, stdout = run_case(
            sim, task["data"], task["trace"], task["log"], task["mode"],
            task["requested_injection_time_s"], task["word_horizon_s"],
            task["delta_21"],
        )
        result = parse_trace(task["trace"])
        result.update({
            "case": task["case"],
            "family": task["family"],
            "Hs_m": task["Hs_m"],
            "source_file": task["source_file"],
            "mode": task["mode"],
            "coordinate": task["coordinate"],
            "direction": task["direction"],
            "sign": task["sign"],
            "target_W": task["target_W"],
            "requested_injection_time_s": task["requested_injection_time_s"],
            "word_horizon_s": task["word_horizon_s"],
            "delta_21": [float(x) for x in task["delta_21"]],
            "covariance_reference": task["covariance_reference"],
            "sim_returncode": int(rc),
            "sim_completed_marker": "OU3_NEIGHBOR_DONE" in stdout,
        })
        return result

    cases = []
    if args.jobs == 1:
        cases = [execute(task) for task in tasks]
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(execute, task) for task in tasks]
            for future in as_completed(futures):
                cases.append(future.result())
    cases.sort(key=lambda c: c.get("case", ""))

    valid = [c for c in cases if c.get("status") in ("PASS_SAMPLED", "FAIL_SAMPLED")]
    all_pass = bool(valid) and len(valid) == len(cases) and all(c["pass_sampled"] for c in valid)
    decrements = [float(c["relative_decrement"]) for c in valid
                  if c.get("relative_decrement") is not None and math.isfinite(float(c["relative_decrement"]))]
    theta = [float(c["theta_max_rad"]) for c in valid
             if c.get("theta_max_rad") is not None and math.isfinite(float(c["theta_max_rad"]))]
    prefix = [float(c["prefix_W_gain_max"]) for c in valid
              if c.get("prefix_W_gain_max") is not None and math.isfinite(float(c["prefix_W_gain_max"]))]
    worst = min(valid, key=lambda c: c.get("relative_decrement", -math.inf)) if valid else None
    report = {
        "schema": 1,
        "status": "PASS_SAMPLED" if all_pass else "FAIL_OR_INCOMPLETE_SAMPLED",
        "qualification": "DENSE_OR_STRUCTURED_SAMPLING_IS_NOT_A_NEIGHBORHOOD_CERTIFICATE",
        "metric": "pairwise zeta^T Sigma_nominal^-1 zeta",
        "case_count": len(cases),
        "valid_endpoint_case_count": len(valid),
        "relative_decrement_min": min(decrements) if decrements else None,
        "theta_max_rad": max(theta) if theta else None,
        "theta_max_deg": math.degrees(max(theta)) if theta else None,
        "prefix_W_gain_max": max(prefix) if prefix else None,
        "worst_case": worst,
        "cases": cases,
        "numerical_neighborhood_certificate": "NOT_ESTABLISHED",
        "deployment_theorem_certificate": "NOT_ESTABLISHED",
    }
    (out / "neighborhood_diagnostic.json").write_text(
        json.dumps(report, indent=2, sort_keys=True)
    )
    print(json.dumps({
        "status": report["status"],
        "case_count": report["case_count"],
        "relative_decrement_min": report["relative_decrement_min"],
        "theta_max_deg": report["theta_max_deg"],
        "prefix_W_gain_max": report["prefix_W_gain_max"],
        "worst_case": None if worst is None else worst.get("case"),
    }, indent=2, sort_keys=True))
    # Sampled failure is scientific diagnostic output; only simulator/tool
    # malfunction is represented through missing/invalid cases in the JSON.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
