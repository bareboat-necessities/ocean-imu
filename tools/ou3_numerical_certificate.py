#!/usr/bin/env python3
"""Numerical source-funnel diagnostics for the deployed OU-III filter.

This tool deliberately separates two claims:
  * finite_replay_certificate: exact for the eight recorded noisy replays that
    were executed (up to floating-point arithmetic);
  * deployment_theorem_certificate: requires continuous-source validated
    enclosures and is never inferred from sampled trajectories.

The replay analysis uses the same noise generator and eight RECORDS inventory
as the OU validation tools.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ou_sweep_common import PATTERNS, RECORDS

REPO = Path(__file__).resolve().parents[1]
TEST_DIR = REPO / "tests" / "kalman_ou_iii"
DEFAULT_DATA_DIR = REPO / "plots" / "kalman_ou_ii"
DEFAULT_OUT = REPO / "reports" / "results" / "ou3_numerical_certificate"

TRACE_HZ = 20.0
WORD_SEC = 1.0
EPS = 1e-10

# These are normalization coordinates, not tuning parameters. They put the
# heterogeneous physical states on comparable numerical scales before the
# node metric is learned.
SCALE_ACTIVE = np.array(
    [1.0] * 3 +          # attitude rotvec [rad]
    [0.02] * 3 +         # gyro bias [rad/s]
    [5.0] * 3 +          # velocity [m/s]
    [10.0] * 3 +         # position [m]
    [100.0] * 3 +        # integral displacement [m s]
    [2.0] * 3 +          # acceleration [m/s2]
    [0.2] * 3,           # accel bias [m/s2]
    dtype=float,
)
SCALE_HELD = SCALE_ACTIVE[:-3]

TAU_BINS = np.array([0.0, 1.5, 2.5, 4.0, 6.0, math.inf])
SIGMA_BINS = np.array([0.0, 0.5, 1.0, 1.5, 2.5, math.inf])
RS_BINS = np.array([0.0, 0.25, 2.0, 10.0, 25.0, math.inf])


@dataclass(frozen=True)
class SourceNode:
    mode: str
    mag_lock: int
    tau_cell: int
    sigma_cell: int
    rs_cell: int

    def label(self) -> str:
        return f"{self.mode}:m{self.mag_lock}:t{self.tau_cell}:s{self.sigma_cell}:r{self.rs_cell}"


def _cell(value: float, edges: np.ndarray) -> int:
    return int(np.searchsorted(edges, value, side="right") - 1)


def source_node(row) -> SourceNode:
    return SourceNode(
        mode="A" if int(row["bias_active"]) else "H",
        mag_lock=int(row["mag_lock"]),
        tau_cell=_cell(float(row["tau_applied"]), TAU_BINS),
        sigma_cell=_cell(float(row["sigma_applied"]), SIGMA_BINS),
        rs_cell=_cell(float(row["rs_applied"]), RS_BINS),
    )


def zyx_matrix(roll_deg: float, pitch_deg: float, yaw_deg: float) -> np.ndarray:
    r, p, y = np.deg2rad([roll_deg, pitch_deg, yaw_deg])
    cr, sr = math.cos(r), math.sin(r)
    cp, sp = math.cos(p), math.sin(p)
    cy, sy = math.cos(y), math.sin(y)
    rx = np.array([[1, 0, 0], [0, cr, -sr], [0, sr, cr]], dtype=float)
    ry = np.array([[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]], dtype=float)
    rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]], dtype=float)
    return rz @ ry @ rx


def so3_log(R: np.ndarray) -> np.ndarray:
    c = float(np.clip((np.trace(R) - 1.0) * 0.5, -1.0, 1.0))
    theta = math.acos(c)
    if theta < 1e-8:
        return 0.5 * np.array(
            [R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]
        )
    if math.pi - theta < 1e-5:
        A = (R + np.eye(3)) * 0.5
        axis = np.sqrt(np.maximum(np.diag(A), 0.0))
        if R[2, 1] - R[1, 2] < 0:
            axis[0] *= -1
        if R[0, 2] - R[2, 0] < 0:
            axis[1] *= -1
        if R[1, 0] - R[0, 1] < 0:
            axis[2] *= -1
        n = np.linalg.norm(axis)
        axis = np.array([1.0, 0.0, 0.0]) if n < EPS else axis / n
        return theta * axis
    return theta / (2.0 * math.sin(theta)) * np.array(
        [R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]
    )


def group_energy(rotvec: np.ndarray) -> float:
    theta = float(np.linalg.norm(rotvec))
    return 1.0 - math.cos(min(theta, math.pi))


def zu_to_ned(v: np.ndarray) -> np.ndarray:
    return np.array([v[1], v[0], -v[2]], dtype=float)


def output_csv_for(input_path: Path) -> Path:
    name = input_path.name
    if name.startswith("wave_data_"):
        name = "w3d_" + name[len("wave_data_") :]
    stem = name[:-4] if name.endswith(".csv") else name
    return input_path.with_name(stem + "_fusion_ou3_cert.csv")


def parse_metrics(stdout: str) -> dict[str, float | None]:
    out: dict[str, float | None] = {}
    for key, pattern in PATTERNS.items():
        m = pattern.search(stdout)
        out[key] = float(m.group(1)) if m else None
    return out


def load_csv_columns(path: Path, names: list[str]) -> dict[str, np.ndarray]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = {name: [] for name in names}
        for rec in reader:
            for name in names:
                text = rec[name]
                try:
                    rows[name].append(float(text))
                except (TypeError, ValueError):
                    rows[name].append(float("nan"))
    return {k: np.asarray(v, dtype=float) for k, v in rows.items()}


def load_trace(path: Path) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)


def nearest_indices(times: np.ndarray, targets: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(times, targets, side="left")
    idx = np.clip(idx, 0, len(times) - 1)
    prev = np.clip(idx - 1, 0, len(times) - 1)
    choose_prev = np.abs(times[prev] - targets) < np.abs(times[idx] - targets)
    return np.where(choose_prev, prev, idx)


def build_error_states(trace: np.ndarray, ts_path: Path) -> tuple[np.ndarray, np.ndarray]:
    cols = [
        "time", "roll_ref", "pitch_ref", "yaw_ref", "roll_est", "pitch_est", "yaw_est",
        "disp_ref_x", "disp_ref_y", "disp_ref_z",
        "vel_ref_x", "vel_ref_y", "vel_ref_z",
        "acc_ref_x", "acc_ref_y", "acc_ref_z",
        "acc_bias_x", "acc_bias_y", "acc_bias_z",
        "gyro_bias_x", "gyro_bias_y", "gyro_bias_z",
    ]
    ts = load_csv_columns(ts_path, cols)
    t = ts["time"]
    idx = nearest_indices(t, np.asarray(trace["time_s"], dtype=float))

    # Truth S = integral of p in NED, using the full 200 Hz truth before
    # downsampling to the certificate trace.
    p_ned = np.column_stack(
        [ts["disp_ref_y"], ts["disp_ref_x"], -ts["disp_ref_z"]]
    )
    s_truth = np.zeros_like(p_ned)
    if len(t) > 1:
        increments = 0.5 * (p_ned[1:] + p_ned[:-1]) * np.diff(t)[:, None]
        s_truth[1:] = np.cumsum(increments, axis=0)

    states = np.zeros((len(trace), 21), dtype=float)
    theta = np.zeros(len(trace), dtype=float)
    for k, j in enumerate(idx):
        Re = zyx_matrix(ts["roll_est"][j], ts["pitch_est"][j], ts["yaw_est"][j])
        Rr = zyx_matrix(ts["roll_ref"][j], ts["pitch_ref"][j], ts["yaw_ref"][j])
        rv = so3_log(Re @ Rr.T)
        states[k, 0:3] = rv
        theta[k] = np.linalg.norm(rv)

        bg_true = zu_to_ned(np.array([ts["gyro_bias_x"][j], ts["gyro_bias_y"][j], ts["gyro_bias_z"][j]]))
        ba_true = zu_to_ned(np.array([ts["acc_bias_x"][j], ts["acc_bias_y"][j], ts["acc_bias_z"][j]]))
        states[k, 3:6] = np.array([trace["bg_x"][k], trace["bg_y"][k], trace["bg_z"][k]]) - bg_true

        v_true = zu_to_ned(np.array([ts["vel_ref_x"][j], ts["vel_ref_y"][j], ts["vel_ref_z"][j]]))
        p_true = p_ned[j]
        a_true = zu_to_ned(np.array([ts["acc_ref_x"][j], ts["acc_ref_y"][j], ts["acc_ref_z"][j]]))
        states[k, 6:9] = np.array([trace["v_x"][k], trace["v_y"][k], trace["v_z"][k]]) - v_true
        states[k, 9:12] = np.array([trace["p_x"][k], trace["p_y"][k], trace["p_z"][k]]) - p_true
        states[k, 12:15] = np.array([trace["S_x"][k], trace["S_y"][k], trace["S_z"][k]]) - s_truth[j]
        states[k, 15:18] = np.array([trace["aw_x"][k], trace["aw_y"][k], trace["aw_z"][k]]) - a_true
        states[k, 18:21] = np.array([trace["ba_x"][k], trace["ba_y"][k], trace["ba_z"][k]]) - ba_true
    return states, theta


def metric_from_samples(X: np.ndarray) -> np.ndarray:
    if len(X) < 2:
        return np.eye(X.shape[1])
    C = (X.T @ X) / max(len(X), 1)
    ridge = max(1e-6, 1e-4 * float(np.trace(C)) / max(C.shape[0], 1))
    return np.linalg.inv(C + ridge * np.eye(C.shape[0]))


def fit_map(X: np.ndarray, Y: np.ndarray, ridge: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    n = X.shape[1]
    gram = X.T @ X + ridge * np.eye(n)
    phi_t = np.linalg.solve(gram, X.T @ Y)
    phi = phi_t.T
    return phi, Y - X @ phi.T


def generalized_lambda(phi: np.ndarray, Pi: np.ndarray, Pj: np.ndarray) -> float:
    L = np.linalg.cholesky(Pi)
    invL = np.linalg.solve(L, np.eye(L.shape[0]))
    M = invL @ phi.T @ Pj @ phi @ invL.T
    M = 0.5 * (M + M.T)
    return float(np.max(np.linalg.eigvalsh(M)))


def analyze_mode(trace: np.ndarray, states: np.ndarray, mode: str, word_rows: int) -> dict:
    active = mode == "A"
    dim = 21 if active else 18
    scale = SCALE_ACTIVE if active else SCALE_HELD
    Z = states[:, :dim] / scale
    nodes = [source_node(row) for row in trace]
    eligible = np.array([
        bool(int(trace["live"][k])) and nodes[k].mode == mode for k in range(len(trace))
    ])

    by_node: dict[SourceNode, list[np.ndarray]] = defaultdict(list)
    for k, ok in enumerate(eligible):
        if ok:
            by_node[nodes[k]].append(Z[k])
    metrics = {node: metric_from_samples(np.asarray(vals)) for node, vals in by_node.items()}

    groups: dict[tuple, list[tuple[int, int]]] = defaultdict(list)
    covered_start = set()
    for a in range(0, len(trace) - word_rows):
        b = a + word_rows
        if not (eligible[a] and eligible[b]):
            continue
        if any(nodes[k].mode != mode for k in range(a, b + 1)):
            continue
        acc_count = int(np.sum(trace["acc_accepted"][a:b]))
        mag_count = int(np.sum(trace["mag_accepted"][a:b]))
        pseudo_count = int(np.sum(trace["pseudo_due_mirror"][a:b]))
        key = (nodes[a], nodes[b], acc_count, mag_count, pseudo_count)
        groups[key].append((a, b))
        covered_start.add(a)

    words = []
    lambdas = []
    prefix_amp = []
    mu_samples = []
    gamma_samples = []
    for key, pairs in groups.items():
        ni, nj, ac, mc, pc = key
        X = np.asarray([Z[a] for a, _ in pairs])
        Y = np.asarray([Z[b] for _, b in pairs])
        if len(X) < 2 or ni not in metrics or nj not in metrics:
            continue
        phi, residual = fit_map(X, Y)
        lam = generalized_lambda(phi, metrics[ni], metrics[nj])
        lambdas.append(lam)

        Wi = np.einsum("bi,ij,bj->b", X, metrics[ni], X)
        Wj = np.einsum("bi,ij,bj->b", Y, metrics[nj], Y)
        denom = np.maximum(np.sum(X * X, axis=1), 1e-12)
        mu = (Wi - Wj) / denom
        mu_samples.extend(mu.tolist())
        gamma = float(max(0.0, np.max(Wj - min(lam, 0.999999) * Wi)))
        gamma_samples.append(gamma)

        pmax = 1.0
        for a, b in pairs:
            base = max(float(Z[a] @ metrics[ni] @ Z[a]), 1e-12)
            for k in range(a + 1, b):
                nk = nodes[k]
                if nk in metrics:
                    wk = float(Z[k] @ metrics[nk] @ Z[k])
                    pmax = max(pmax, wk / base)
        prefix_amp.append(pmax)

        words.append({
            "start": ni.label(), "end": nj.label(),
            "acc_count": ac, "mag_count": mc, "pseudo_count": pc,
            "samples": len(pairs), "lambda_gen": lam,
            "residual_inf": float(np.max(np.abs(residual))) if residual.size else 0.0,
            "gamma_replay": gamma,
            "mu_min_replay": float(np.min(mu)),
            "mu_p05_replay": float(np.quantile(mu, 0.05)),
        })

    possible = max(0, int(np.sum(eligible)) - word_rows)
    coverage = len(covered_start) / possible if possible else 0.0
    return {
        "mode": mode, "dimension": dim,
        "node_count": len(metrics), "word_family_count": len(words),
        "word_start_coverage": coverage,
        "lambda_worst": max(lambdas) if lambdas else None,
        "B_pre_replay": max(prefix_amp) if prefix_amp else None,
        "mu_min_replay": min(mu_samples) if mu_samples else None,
        "mu_p05_replay": float(np.quantile(mu_samples, 0.05)) if mu_samples else None,
        "gamma_worst_replay": max(gamma_samples) if gamma_samples else None,
        "linear_candidate_pass": bool(lambdas) and max(lambdas) < 1.0,
        "replay_word_coverage_pass": coverage > 0.98,
        "words": words,
    }


def handoff_and_hybrid(trace: np.ndarray, states: np.ndarray, word_rows: int) -> dict:
    live = np.asarray(trace["live"], dtype=int)
    bias = np.asarray(trace["bias_active"], dtype=int)
    mag_lock = np.asarray(trace["mag_lock"], dtype=int)
    refined = np.asarray(trace["mag_refined"], dtype=int)
    live_idx = np.flatnonzero((live[1:] > live[:-1])) + 1
    bias_idx = np.flatnonzero((bias[1:] > bias[:-1])) + 1
    lock_idx = np.flatnonzero((mag_lock[1:] > mag_lock[:-1])) + 1
    refine_idx = np.flatnonzero((refined[1:] > refined[:-1])) + 1

    handoff = int(live_idx[0]) if len(live_idx) else (0 if live[0] else -1)
    theta0 = float(np.linalg.norm(states[handoff, :3])) if handoff >= 0 else None
    z = states / SCALE_ACTIVE
    energy = np.sum(z * z, axis=1)
    if len(energy):
        tail_start = max(0, int(0.8 * len(energy)))
        inner = float(np.quantile(energy[tail_start:], 0.999))
        capture = None
        if handoff >= 0:
            for k in range(handoff, len(energy)):
                if np.max(energy[k:]) <= inner * (1.0 + 1e-9):
                    capture = k
                    break
    else:
        inner, capture = None, None
    trace_dt = float(np.median(np.diff(trace["time_s"]))) if len(trace) > 1 else 1.0 / TRACE_HZ
    return {
        "live_handoff_time_s": float(trace["time_s"][handoff]) if handoff >= 0 else None,
        "handoff_theta_deg": math.degrees(theta0) if theta0 is not None else None,
        "bias_release_time_s": float(trace["time_s"][bias_idx[0]]) if len(bias_idx) else None,
        "mag_lock_time_s": float(trace["time_s"][lock_idx[0]]) if len(lock_idx) else None,
        "mag_refine_time_s": float(trace["time_s"][refine_idx[0]]) if len(refine_idx) else None,
        "inner_replay_level": inner,
        "capture_word_count_replay": (
            int(math.ceil((capture - handoff) / max(word_rows, 1)))
            if capture is not None and handoff >= 0 else None
        ),
        "capture_time_s_replay": (
            float((capture - handoff) * trace_dt)
            if capture is not None and handoff >= 0 else None
        ),
        "hybrid_events_observed": {
            "bias_release": int(len(bias_idx)),
            "mag_lock": int(len(lock_idx)),
            "mag_refine": int(len(refine_idx)),
        },
    }


def stochastic_diagnostic(word_increments: np.ndarray, horizon_words: int) -> dict:
    # Raw test-noise channels are normalized by their own standard deviations
    # before applying the Gaussian quadratic-form localization inequality.
    d = 9
    w_star_sq = float((6.0 * math.sqrt(d)) ** 2)
    lo, hi = 0.0, w_star_sq
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        q = d + 2.0 * math.sqrt(d * mid) + 2.0 * mid
        if q <= w_star_sq:
            lo = mid
        else:
            hi = mid
    t_star = lo
    localization = min(1.0, horizon_words * math.exp(-t_star))

    if len(word_increments):
        centered = word_increments - np.mean(word_increments)
        b = float(np.max(np.maximum(centered, 0.0)))
        v = float(np.var(centered))
        x = max(float(np.quantile(word_increments, 0.999) - np.mean(word_increments)), 0.0)
        if x > 0 and (v > 0 or b > 0):
            freedman = min(1.0, horizon_words * math.exp(-(x * x) / (2.0 * (v + b * x / 3.0))))
        else:
            freedman = 0.0
    else:
        b = v = freedman = None
    total = None if freedman is None else min(1.0, localization + freedman)
    return {
        "noise_coordinate": "per-sensor-sigma-normalized pre-gate Gaussian increment",
        "gaussian_dimension": d,
        "w_star_normalized": math.sqrt(w_star_sq),
        "t_star": t_star,
        "localization_union_bound": localization,
        "b_W_replay": b, "v_W_replay": v,
        "freedman_replay_bound": freedman,
        "combined_replay_bound": total,
        "qualification": "DIAGNOSTIC_EMPIRICAL_BW_VW",
    }


def analyze_record(family: str, hs: float, trace_path: Path, ts_path: Path,
                   metrics: dict, rms_pass: bool) -> dict:
    trace = load_trace(trace_path)
    states, theta = build_error_states(trace, ts_path)
    trace_dt = float(np.median(np.diff(trace["time_s"]))) if len(trace) > 1 else 1.0 / TRACE_HZ
    word_rows = max(1, int(round(WORD_SEC / trace_dt)))

    held = analyze_mode(trace, states, "H", word_rows)
    active = analyze_mode(trace, states, "A", word_rows)
    handoff = handoff_and_hybrid(trace, states, word_rows)
    V = np.asarray([group_energy(rv) for rv in states[:, :3]])
    word_inc = V[word_rows:] - V[:-word_rows] if len(V) > word_rows else np.array([])
    stochastic = stochastic_diagnostic(word_inc, max(1, len(V) // word_rows))

    theta_max = float(np.max(theta)) if len(theta) else math.nan
    source_modes = [m for m in (held, active) if m["word_family_count"]]
    linear_pass = bool(source_modes) and all(m["linear_candidate_pass"] for m in source_modes)
    coverage_pass = bool(source_modes) and all(m["replay_word_coverage_pass"] for m in source_modes)
    finite_replay_pass = bool(
        rms_pass and np.all(np.isfinite(states)) and theta_max < math.pi and coverage_pass
    )

    return {
        "family": family, "Hs_m": hs, "trace_rows": int(len(trace)),
        "theta_max_deg": math.degrees(theta_max),
        "group_energy_max": float(np.max(V)) if len(V) else None,
        "all_attitude_inside_pi": bool(theta_max < math.pi),
        "rms_regression_pass": rms_pass, "rms_metrics": metrics,
        "held_mode": held, "active_mode": active,
        "handoff_hybrid": handoff, "stochastic": stochastic,
        "linear_candidate_pass": linear_pass,
        "finite_replay_source_coverage_pass": coverage_pass,
        "finite_replay_certificate": "PASS" if finite_replay_pass else "FAIL",
        "deployment_theorem_certificate": "NOT_ESTABLISHED",
        "deployment_missing": [
            "validated continuous-source word-family enclosure",
            "robust verified path LMIs over every source cell",
            "rigorous large-angle sector theta_star per node",
            "rigorous nonlinear source-word infimum mu_W",
            "rigorous handoff/hybrid funnel inequalities",
            "non-empirical martingale b_W and v_W bounds",
        ],
    }


def run_record(exe: Path, data_path: Path, out_dir: Path) -> tuple[Path, Path, dict, bool, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    trace_path = out_dir / (data_path.stem + "_certificate_trace.csv")
    env = os.environ.copy()
    env["OU3_CERT_TRACE"] = str(trace_path)
    env.setdefault("OU3_CERT_TRACE_STRIDE", "10")
    env["W3D_WRITE_TIMESERIES"] = "1"
    env["W3D_VALIDATION_WINDOW_SEC"] = "900"
    proc = subprocess.run(
        [str(exe), "--input", str(data_path)], cwd=data_path.parent, env=env,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    ts_path = output_csv_for(data_path)
    if not trace_path.exists():
        raise RuntimeError(f"certificate trace missing for {data_path.name}\n{proc.stdout[-4000:]}")
    if not ts_path.exists():
        raise RuntimeError(f"timeseries missing for {data_path.name}\n{proc.stdout[-4000:]}")
    return trace_path, ts_path, parse_metrics(proc.stdout), proc.returncode == 0, proc.stdout


def report_markdown(report: dict) -> str:
    lines = [
        "# OU-III numerical source-funnel certificate", "",
        f"Overall eight-replay status: **{report['finite_replay_certificate']}**", "",
        "The finite replay status checks the exact eight noisy reference replays. "
        "It is deliberately not promoted to a deployment theorem certificate: "
        "continuous source cells still require validated enclosure.", "",
        "| Sea | RMS | replay | max θ | linear candidate | H words | A words | theorem |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in report["records"]:
        lines.append(
            f"| {r['family']} {r['Hs_m']:.2f} | "
            f"{'PASS' if r['rms_regression_pass'] else 'FAIL'} | "
            f"{r['finite_replay_certificate']} | {r['theta_max_deg']:.3f}° | "
            f"{'PASS' if r['linear_candidate_pass'] else 'FAIL'} | "
            f"{r['held_mode']['word_family_count']} | {r['active_mode']['word_family_count']} | "
            f"{r['deployment_theorem_certificate']} |"
        )
    lines += [
        "", "## Claim boundary", "",
        "- `finite_replay_certificate`: exact check of the executed eight noisy traces "
        "with replay-bounded word disturbances and source coverage.",
        "- `linear_candidate_pass`: path-metric diagnostic from fitted word maps. "
        "It is a candidate metric, not an interval proof.",
        "- `deployment_theorem_certificate`: remains `NOT_ESTABLISHED` until the "
        "continuous source families, nonlinear word maps, funnel jumps, and "
        "martingale constants are enclosed independently of sampled trajectories.",
        "", "The JSON report contains the worst generalized word factor, replay prefix "
        "amplification, nonlinear decrement diagnostics, startup capture, hybrid "
        "events, and Gaussian/Freedman diagnostics for each sea.", "",
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--sim", type=Path, default=TEST_DIR / "ou3-certificate-sim")
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--analysis-only", action="store_true",
                    help="analyze already-generated traces/timeseries")
    args = ap.parse_args()

    if not args.no_build and not args.analysis_only:
        subprocess.run(["make", "-C", str(TEST_DIR), "ou3-certificate-sim"], check=True)

    records = []
    logs_dir = args.output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    for family, hs, filename in RECORDS:
        data_path = args.data_dir / filename
        if not data_path.exists():
            raise FileNotFoundError(f"missing reference record: {data_path}")
        slug = f"{family.lower().replace('-', '_')}_{hs:.2f}".replace(".", "_")
        trace_path = args.output_dir / (data_path.stem + "_certificate_trace.csv")
        ts_path = output_csv_for(data_path)
        metrics: dict[str, float | None] = {}
        rms_pass = True
        if not args.analysis_only:
            trace_path, ts_path, metrics, rms_pass, log = run_record(
                args.sim, data_path, args.output_dir
            )
            (logs_dir / f"{slug}.log").write_text(log)
        records.append(analyze_record(family, hs, trace_path, ts_path, metrics, rms_pass))

    finite_pass = all(r["finite_replay_certificate"] == "PASS" for r in records)
    report = {
        "schema": 1,
        "scope": "eight_reference_noisy_replays_default_test_seeds",
        "record_count": len(records),
        "finite_replay_certificate": "PASS" if finite_pass else "FAIL",
        "deployment_theorem_certificate": "NOT_ESTABLISHED",
        "records": records,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "certificate.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    (args.output_dir / "certificate.md").write_text(report_markdown(report))
    print(report_markdown(report))
    # A failed RMS regression is a real implementation regression. Candidate
    # metric failure is reported but does not make CI red: the point of this PR
    # is to expose whether the current estimator closes the certificate.
    return 1 if any(not r["rms_regression_pass"] for r in records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
