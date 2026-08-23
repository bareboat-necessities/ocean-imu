#!/usr/bin/env python3
"""Executed-funnel accounting in the source-varying OU-III information metric.

This stage is intentionally downstream of ``ou3_information_certificate.py``.
It does not invent a new metric: every endpoint uses the full estimator
covariance recorded at that same source point,

    W_k = e_k^T Sigma_k^{-1} e_k.

The nominal noisy replay is used only to measure disturbance allowance, startup
handoff, finite capture and observed hybrid jumps. It is *not* promoted to a
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

CAPTURE_REL_SLACK = 1.0e-3
CAPTURE_ABS_SLACK = 1.0e-12


def info_energy(e: np.ndarray, Sigma: np.ndarray) -> float:
    Sigma = 0.5 * (Sigma + Sigma.T)
    return float(e @ np.linalg.solve(Sigma, e))


def nearest_row(times: np.ndarray, t: float) -> int:
    i = int(np.searchsorted(times, t))
    i = max(0, min(i, len(times) - 1))
    j = max(0, i - 1)
    return j if abs(times[j] - t) <= abs(times[i] - t) else i


def strict_capture_level(asymptotic_floor: float) -> float:
    """Choose a deterministic strict superlevel of the affine fixed point.

    For c_{n+1}=lambda*c_n+gamma with 0<=lambda<1, the minimal invariant
    fixed point b*=gamma/(1-lambda) is approached asymptotically from above.
    Finite capture is therefore claimed only into b_eta>b*.  The declared
    relative slack is certificate bookkeeping, not filter tuning.
    """
    b = max(0.0, float(asymptotic_floor))
    return b * (1.0 + CAPTURE_REL_SLACK) + CAPTURE_ABS_SLACK


def finite_capture_steps(c0: float, lam: float, gamma: float, target: float,
                         cap: int = 100000) -> int | None:
    if not (0.0 <= lam < 1.0) or gamma < 0.0 or target < 0.0:
        return None
    fixed = gamma / max(1e-15, 1.0 - lam)
    # No finite-entry claim is valid for a target at or below the fixed point
    # unless the initial level is already inside it.
    if c0 > target and target <= fixed:
        return None
    c = max(0.0, c0)
    for n in range(cap + 1):
        if c <= target:
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


def record_name_index() -> dict[str, str]:
    """Map both certificate slugs and exact replay file stems to source CSVs."""
    index: dict[str, str] = {}
    for family, hs, name in BASE.RECORDS:
        short = f"{family.lower().replace('-', '_')}_{hs:.2f}".replace(".", "_")
        index[short] = name
        index[Path(name).stem] = name
    return index


def evaluate_mode(record_data: dict, mode: str, horizon_s: float, lam_bound: float) -> dict:
    dim = 21 if mode == "A" else 18
    residuals = []
    ratios = []
    starts = []
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
            if W0 > 1e-12:
                ratios.append(W1 / W0)
            if worst is None or r > worst[0]:
                worst = (r, slug, s, s + n - 1, W0, W1)

    if not residuals:
        return {"mode": mode, "status": "NO_WORDS", "horizon_s": horizon_s}
    gamma = max(0.0, float(np.max(residuals)))
    b_star = gamma / max(1e-15, 1.0 - lam_bound)
    b_eta = strict_capture_level(b_star)
    c0 = float(np.max(starts))
    N = finite_capture_steps(c0, lam_bound, gamma, b_eta)
    return {
        "mode": mode,
        "status": "PASS" if N is not None else "FAIL",
        "horizon_s": horizon_s,
        "lambda_information_bound": lam_bound,
        "strict_margin_1_minus_lambda": 1.0 - lam_bound,
        "word_count": len(residuals),
        "gamma_replay": gamma,
        "asymptotic_floor_b_star_replay": b_star,
        "invariant_level_b_replay": b_star,
        "finite_capture_level_b_eta_replay": b_eta,
        "capture_relative_slack": CAPTURE_REL_SLACK,
        "c0_executed_word_starts": c0,
        "c0_over_b": c0 / max(b_star, 1e-15),
        "c0_over_capture_level": c0 / max(b_eta, 1e-15),
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


def evaluate_contracting_horizons(record_data: dict, mode: str,
                                  attempts: list[dict]) -> tuple[dict, list[dict]]:
    """Evaluate every strictly contracting information horizon and retain the tightest funnel.

    The linear certificate intentionally selects the *first* contracting horizon
    to show the shortest source-complete word. That is not necessarily the best
    horizon for the affine noisy funnel: when 1-lambda is tiny, gamma/(1-lambda)
    can be enormous. For replay-funnel accounting we therefore evaluate every
    already-certified contracting horizon and choose the smallest invariant
    fixed-point level b*. This changes no filter behavior and no linear
    certificate claim.
    """
    candidates = []
    for a in attempts:
        if not bool(a.get("information_pass")):
            continue
        lam = float(a["lambda_worst_information"])
        if not (0.0 <= lam < 1.0):
            continue
        candidates.append(evaluate_mode(record_data, mode, float(a["horizon_s"]), lam))
    feasible = [r for r in candidates
                if r.get("status") == "PASS" and math.isfinite(float(r.get("invariant_level_b_replay", math.inf)))]
    if not feasible:
        return {"mode": mode, "status": "NO_CONTRACTING_FUNNEL"}, candidates
    selected = min(feasible, key=lambda r: (float(r["invariant_level_b_replay"]), float(r["horizon_s"])))
    selected = dict(selected)
    selected["selection_basis"] = "MINIMUM_REPLAY_INVARIANT_LEVEL_B_OVER_CERTIFIED_HORIZONS"
    return selected, candidates


def handoff_and_hybrid(record_data: dict) -> dict:
    handoffs = []
    jumps = []
    for slug, d in record_data.items():
        trace, E, maps, covs = d["trace"], d["E"], d["maps"], d["covs"]
        tt = np.asarray(trace["time_s"], float)
        live = np.asarray(trace["live"], int)
        rising_live = np.flatnonzero(live[1:] > live[:-1]) + 1
        if len(rising_live):
            k = int(rising_live[0])
            t = float(tt[k])
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
            k0 = nearest_row(tt, b.t0)
            k1 = nearest_row(tt, b.t1)
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
           "This report uses `W=e^T Sigma_KF^-1 e` at the source point. It is executed-replay accounting, not a neighborhood theorem.",
           "For the noisy affine funnel it selects the certified horizon with the smallest replay fixed point `b*=gamma/(1-lambda)`.",
           f"Finite capture is claimed only into the declared strict superlevel `b_eta=b*(1+{CAPTURE_REL_SLACK})+{CAPTURE_ABS_SLACK}`.", ""]
    for key in ("held", "active"):
        r = report[key]
        out.append(f"{key.capitalize()}: {r.get('status')}, horizon {r.get('horizon_s')} s, "
                   f"lambda {r.get('lambda_information_bound')}, gamma {r.get('gamma_replay')}, "
                   f"b* {r.get('asymptotic_floor_b_star_replay')}, b_eta {r.get('finite_capture_level_b_eta_replay')}, "
                   f"c0/b_eta {r.get('c0_over_capture_level')}, N {r.get('capture_words_N')}")
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
    out = args.certificate_dir.resolve()
    data_dir = args.data_dir.resolve()
    info = json.loads((out / "information_certificate.json").read_text())
    if info.get("status") != "PASS":
        report = {"schema": 3, "status": "BLOCKED_AT_INFORMATION_LINEAR_GATE",
                  "information_status": info.get("status"),
                  "numerical_neighborhood_certificate": "NOT_ESTABLISHED"}
        (out / "information_completion.json").write_text(json.dumps(report, indent=2, sort_keys=True))
        (out / "information_completion.md").write_text(markdown({
            "status": report["status"], "held": {}, "active": {},
            "handoff_hybrid": {}}))
        print(report["status"])
        return 0

    record_lookup = record_name_index()
    record_data = {}
    for map_path in sorted(out.glob("*_exact_maps.bin")):
        slug = map_path.name.replace("_exact_maps.bin", "")
        record_name = record_lookup.get(slug)
        if record_name is None:
            known = ", ".join(sorted(record_lookup))
            raise KeyError(f"unknown exact-map record identity {slug!r}; known identities: {known}")
        maps, covs, _ = INFO.pair_map_covariance(map_path, slug)
        data_path = (data_dir / record_name).resolve()
        trace_path = out / f"{data_path.stem}_certificate_trace.csv"
        if not trace_path.exists():
            raise FileNotFoundError(f"certificate trace for {slug}: {trace_path}")
        timeseries = BASE.output_csv_for(data_path)
        if not timeseries.exists():
            raise FileNotFoundError(timeseries)
        trace, E, theta = record_errors(trace_path, timeseries)
        record_data[slug] = {"maps": maps, "covs": covs, "trace": trace, "E": E, "theta": theta}

    held, held_attempts = evaluate_contracting_horizons(record_data, "H", info["held"]["attempts"])
    active, active_attempts = evaluate_contracting_horizons(record_data, "A", info["active"]["attempts"])
    hh = handoff_and_hybrid(record_data)
    status = "PASS_EXECUTED_REPLAY" if held.get("status") == "PASS" and active.get("status") == "PASS" else "FAIL_EXECUTED_REPLAY"
    report = {"schema": 3, "status": status, "metric": "e^T Sigma_KF^-1 e",
              "funnel_horizon_selection": "MINIMUM_REPLAY_INVARIANT_LEVEL_B_OVER_CERTIFIED_HORIZONS",
              "finite_capture_definition": "b_eta=b_star*(1+1e-3)+1e-12, strictly above b_star=gamma/(1-lambda)",
              "held": held, "active": active,
              "held_attempts": held_attempts, "active_attempts": active_attempts,
              "handoff_hybrid": hh,
              "numerical_neighborhood_certificate": "NOT_ESTABLISHED",
              "deployment_theorem_certificate": "NOT_ESTABLISHED"}
    (out / "information_completion.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    text = markdown(report)
    (out / "information_completion.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
