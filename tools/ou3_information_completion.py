#!/usr/bin/env python3
"""Executed-funnel accounting in the source-varying OU-III information metric.

This stage is intentionally downstream of ``ou3_information_certificate.py``.
It does not invent a new metric: every endpoint uses the full estimator
covariance recorded at that same source point,

    W_k = e_k^T Sigma_k^{-1} e_k.

The nominal noisy replay is used only to measure disturbance allowance, startup
handoff, finite capture and observed hybrid jumps.  It is *not* promoted to a
nonzero-neighborhood or continuous-source theorem certificate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import ou3_numerical_certificate as BASE
import ou3_information_certificate as INFO


def info_energy(e: np.ndarray, Sigma: np.ndarray) -> float:
    Sigma = 0.5 * (Sigma + Sigma.T)
    return float(e @ np.linalg.solve(Sigma, e))


def nearest_row(times: np.ndarray, t: float) -> int:
    i = int(np.searchsorted(times, t))
    i = max(0, min(i, len(times) - 1))
    j = max(0, i - 1)
    return j if abs(times[j] - t) <= abs(times[i] - t) else i


def finite_capture_steps(c0: float, lam: float, gamma: float, b: float,
                         cap: int = 100000) -> int | None:
    if not (0.0 <= lam < 1.0) or gamma < 0.0 or b < 0.0:
        return None
    c = max(0.0, c0)
    for n in range(cap + 1):
        if c <= b * (1.0 + 1e-12):
            return n
        c = lam * c + gamma
    return None


def valid_word(maps: list[BASE.MapBlock], mode: str, start: int, n: int) -> bool:
    seq = maps[start:start + n]
    return (len(seq) == n and
            all(b.valid and not b.hybrid_jump and b.mode == mode for b in seq) and
            all(abs(seq[k + 1].t0 - seq[k].t1) <= 5e-3 for k in range(n - 1)))


def record_errors(trace_path: Path, timeseries: Path):
    trace = np.genfromtxt(trace_path, delimiter=",", names=True, dtype=None, encoding=None)
    E, theta = BASE.build_error_states(trace, timeseries)
    return trace, E, theta


def evaluate_mode(record_data: dict, mode: str, horizon_s: float, lam_bound: float) -> dict:
    dim = 21 if mode == "A" else 18
    residuals = []
    ratios = []
    starts = []
    endpoints = []
    worst = None

    for slug, d in record_data.items():
        maps = d["maps"]
        covs = d["covs"]
        trace = d["trace"]
        E = d["E"]
        tt = np.asarray(trace["time_s"], float)
        if not maps:
            continue
        base = float(np.median([b.t1 - b.t0 for b in maps if b.t1 > b.t0]))
        n = max(1, int(round(horizon_s / base)))
        for s in range(len(maps) - n + 1):
            if not valid_word(maps, mode, s, n):
                continue
            e0i = nearest_row(tt, maps[s].t0)
            e1i = nearest_row(tt, maps[s + n - 1].t1)
            e0 = E[e0i, :dim]
            e1 = E[e1i, :dim]
            W0 = info_energy(e0, covs[s].start[:dim, :dim])
            W1 = info_energy(e1, covs[s + n - 1].end[:dim, :dim])
            r = W1 - lam_bound * W0
            residuals.append(r)
            starts.append(W0)
            endpoints.append(W1)
            if W0 > 1e-12:
                ratios.append(W1 / W0)
            if worst is None or r > worst[0]:
                worst = (r, slug, s, s + n - 1, W0, W1)

    if not residuals:
        return {"mode": mode, "status": "NO_WORDS"}
    gamma = max(0.0, float(np.max(residuals)))
    # Smallest closed affine invariant level for W+ <= lambda W + gamma.
    b = gamma / max(1e-15, 1.0 - lam_bound)
    c0 = float(np.max(starts))
    N = finite_capture_steps(c0, lam_bound, gamma, b)
    return {
        "mode": mode,
        "status": "PASS" if N is not None else "FAIL",
        "horizon_s": horizon_s,
        "lambda_information_bound": lam_bound,
        "word_count": len(residuals),
        "gamma_replay": gamma,
        "invariant_level_b_replay": b,
        "c0_executed_word_starts": c0,
        "capture_words_N": N,
        "capture_time_s": None if N is None else N * horizon_s,
        "observed_endpoint_ratio_max": float(np.max(ratios)) if ratios else None,
        "observed_endpoint_ratio_p99": float(np.quantile(ratios, 0.99)) if ratios else None,
        "worst_replay_word": None if worst is None else {
            "record": worst[1], "start_block": worst[2], "end_block": worst[3],
            "residual": worst[0], "W0": worst[4], "W1": worst[5],
        },
        "qualification": "EXECUTED_NOISY_REPLAY_ONLY",
    }


def handoff_and_hybrid(record_data: dict) -> dict:
    handoffs = []
    jumps = []
    for slug, d in record_data.items():
        trace, E, maps, covs = d["trace"], d["E"], d["maps"], d["covs"]
        tt = np.asarray(trace["time_s"], float)
        live = np.asarray(trace["live"], int)
        rising_live = np.flatnonzero(live[1:] > live[:-1]) + 1
        if len(rising_live):
            k = int(rising_live[0]); t = float(tt[k])
            # First covariance block whose end is at/after the handoff.
            bi = min(range(len(maps)), key=lambda i: abs(maps[i].t1 - t)) if maps else None
            if bi is not None:
                dim = 21 if maps[bi].end_active else 18
                W = info_energy(E[k, :dim], covs[bi].end[:dim, :dim])
                handoffs.append({"record": slug, "time_s": t, "W_information": W,
                                 "theta_deg": math.degrees(float(np.linalg.norm(E[k, :3]))),
                                 "mode": "A" if dim == 21 else "H"})

        for i, b in enumerate(maps):
            if not b.hybrid_jump:
                continue
            k0 = nearest_row(tt, b.t0); k1 = nearest_row(tt, b.t1)
            d0 = 21 if b.start_active else 18
            d1 = 21 if b.end_active else 18
            W0 = info_energy(E[k0, :d0], covs[i].start[:d0, :d0])
            W1 = info_energy(E[k1, :d1], covs[i].end[:d1, :d1])
            jumps.append({"record": slug, "block": i, "t0": b.t0, "t1": b.t1,
                          "from_mode": "A" if b.start_active else "H",
                          "to_mode": "A" if b.end_active else "H",
                          "mag_lock_change": bool(b.start_lock != b.end_lock),
                          "mag_refine_change": bool(b.start_refined != b.end_refined),
                          "W_before": W0, "W_after": W1,
                          "amplification": W1 / max(W0, 1e-15)})
    return {
        "handoffs": handoffs,
        "handoff_W_max": max((x["W_information"] for x in handoffs), default=None),
        "hybrid_jumps": jumps,
        "hybrid_amplification_max": max((x["amplification"] for x in jumps), default=None),
        "qualification": "EXECUTED_NOISY_REPLAY_ONLY",
    }


def markdown(report: dict) -> str:
    out = ["# OU-III information-metric funnel accounting", "",
           f"Status: **{report['status']}**", "",
           "This report uses `W=e^T Sigma_KF^-1 e` at the source point. It is executed-replay accounting, not a neighborhood theorem.", ""]
    for key in ("held", "active"):
        r = report[key]
        out.append(f"{key.capitalize()}: {r.get('status')}, horizon {r.get('horizon_s')} s, "
                   f"lambda {r.get('lambda_information_bound')}, gamma {r.get('gamma_replay')}, "
                   f"b {r.get('invariant_level_b_replay')}, N {r.get('capture_words_N')}")
    hh = report["handoff_hybrid"]
    out += ["", f"Max executed handoff W: {hh.get('handoff_W_max')}",
            f"Max executed hybrid amplification: {hh.get('hybrid_amplification_max')}", "",
            "`numerical_neighborhood_certificate` remains NOT_ESTABLISHED until nonzero perturbation/validated enclosure closes the same inequalities."]
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=BASE.DEFAULT_OUT)
    ap.add_argument("--data-dir", type=Path, default=BASE.DEFAULT_DATA_DIR)
    args = ap.parse_args()
    out = args.certificate_dir.resolve(); data_dir = args.data_dir.resolve()
    info = json.loads((out / "information_certificate.json").read_text())
    if info.get("status") != "PASS":
        report = {"schema": 1, "status": "BLOCKED_AT_INFORMATION_LINEAR_GATE",
                  "information_status": info.get("status"),
                  "numerical_neighborhood_certificate": "NOT_ESTABLISHED"}
        (out / "information_completion.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        (out / "information_completion.md").write_text(markdown({
            "status": report["status"], "held": {}, "active": {},
            "handoff_hybrid": {}}))
        print(report["status"])
        return 0

    record_lookup = {f"{fam.lower().replace('-','_')}_{hs:.2f}".replace('.', '_'): name
                     for fam, hs, name in BASE.RECORDS}
    record_data = {}
    for map_path in sorted(out.glob("*_exact_maps.bin")):
        slug = map_path.name.replace("_exact_maps.bin", "")
        maps, covs, _ = INFO.pair_map_covariance(map_path, slug)
        trace_path = next(out.glob(f"*{record_lookup[slug].replace('.csv','')}*_certificate_trace.csv"), None)
        # The trace is named from the original data stem; use the exact name if glob is ambiguous.
        data_path = (data_dir / record_lookup[slug]).resolve()
        exact_trace = out / f"{data_path.stem}_certificate_trace.csv"
        if exact_trace.exists(): trace_path = exact_trace
        if trace_path is None or not trace_path.exists():
            raise FileNotFoundError(f"certificate trace for {slug}")
        timeseries = BASE.output_csv_for(data_path)
        if not timeseries.exists(): raise FileNotFoundError(timeseries)
        trace, E, theta = record_errors(trace_path, timeseries)
        record_data[slug] = {"maps": maps, "covs": covs, "trace": trace, "E": E, "theta": theta}

    hs = info["held"]["selected"]
    ac = info["active"]["selected"]
    held = evaluate_mode(record_data, "H", float(hs["horizon_s"]), float(hs["lambda_worst_information"]))
    active = evaluate_mode(record_data, "A", float(ac["horizon_s"]), float(ac["lambda_worst_information"]))
    hh = handoff_and_hybrid(record_data)
    status = "PASS_EXECUTED_REPLAY" if held.get("status") == "PASS" and active.get("status") == "PASS" else "FAIL_EXECUTED_REPLAY"
    report = {"schema": 1, "status": status, "metric": "e^T Sigma_KF^-1 e",
              "held": held, "active": active, "handoff_hybrid": hh,
              "numerical_neighborhood_certificate": "NOT_ESTABLISHED",
              "deployment_theorem_certificate": "NOT_ESTABLISHED"}
    (out / "information_completion.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    text = markdown(report); (out / "information_completion.md").write_text(text); print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
