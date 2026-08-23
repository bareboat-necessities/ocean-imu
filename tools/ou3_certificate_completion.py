#!/usr/bin/env python3
"""Complete the OU-III eight-replay certificate after the exact path-LMI solve.

This stage deliberately distinguishes three claims:
  1. exact linear executed-word certificate (imported from certificate.json),
  2. executed nonlinear/funnel accounting on the same eight noisy replays,
  3. neighborhood/deployment theorem certification, which requires validated
     source-cell enclosures and is never inferred from sampled trajectories.

The purpose of this tool is to make every remaining theorem obligation explicit
and machine-checkable instead of leaving a prose gap after the path-LMI stage.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

import ou3_numerical_certificate as base
from ou_sweep_common import RECORDS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DATA = REPO / "plots" / "kalman_ou_ii"
DEFAULT_CERT = REPO / "reports" / "results" / "ou3_numerical_certificate"
EPS = 1e-12


def load_metrics(path: Path) -> dict[str, dict[str, np.ndarray]]:
    if not path.exists():
        return {"H": {}, "A": {}}
    z = np.load(path, allow_pickle=False)
    ans: dict[str, dict[str, np.ndarray]] = {"H": {}, "A": {}}
    for mode in ("H", "A"):
        lk, pk = f"{mode}_labels", f"{mode}_P_local"
        if lk not in z or pk not in z:
            continue
        labels = [str(x) for x in z[lk]]
        mats = np.asarray(z[pk], float)
        ans[mode] = {label: mats[i] for i, label in enumerate(labels)}
    return ans


def nearest_index(t: np.ndarray, value: float) -> int:
    i = int(np.clip(np.searchsorted(t, value), 0, len(t) - 1))
    j = max(0, i - 1)
    return j if abs(float(t[j]) - value) < abs(float(t[i]) - value) else i


def source_node_from_trace(trace, k: int) -> base.SourceNode:
    return base.make_node(
        bool(int(trace["bias_active"][k])),
        bool(int(trace["mag_lock"][k])),
        bool(int(trace["mag_refined"][k])),
        float(trace["tau_applied"][k]),
        float(trace["sigma_applied"][k]),
        float(trace["rs_applied"][k]),
    )


def mode_scale(mode: str) -> np.ndarray:
    return base.SCALE_ACTIVE if mode == "A" else base.SCALE_HELD


def state_for_mode(e: np.ndarray, mode: str) -> np.ndarray:
    return np.asarray(e[:21 if mode == "A" else 18], float)


def group_metric_value(e: np.ndarray, P: np.ndarray, mode: str) -> float:
    x = state_for_mode(e, mode)
    s = mode_scale(mode)
    z = x / s
    a_R = 2.0 * float(P[0, 0])
    theta = float(np.linalg.norm(x[:3]))
    VR = 1.0 - math.cos(min(theta, math.pi))
    zxi = z[3:]
    Pxi = P[3:, 3:]
    return float(a_R * VR + zxi @ Pxi @ zxi)


def group_denominator(e: np.ndarray, mode: str) -> float:
    x = state_for_mode(e, mode)
    z = x / mode_scale(mode)
    VR = 1.0 - math.cos(min(float(np.linalg.norm(x[:3])), math.pi))
    return max(EPS, VR + float(z[3:] @ z[3:]))


def word_endpoint_energy(word: base.Word, trace, E: np.ndarray,
                         metrics: dict[str, np.ndarray]) -> dict | None:
    t = np.asarray(trace["time_s"], float)
    # Map blocks are [t0,t1]; word indices refer to the block inventory.
    # The exact block times are recovered from the word's start/end indices by
    # the caller and attached transiently as _t0/_t1.
    t0 = float(getattr(word, "_t0"))
    t1 = float(getattr(word, "_t1"))
    i0, i1 = nearest_index(t, t0), nearest_index(t, t1)
    li, lj = word.start_node.label(), word.end_node.label()
    if li not in metrics or lj not in metrics:
        return None
    Pi, Pj = metrics[li], metrics[lj]
    W0 = group_metric_value(E[i0], Pi, word.start_node.mode)
    W1 = group_metric_value(E[i1], Pj, word.end_node.mode)
    lam = base.generalized_lambda(word.phi, Pi, Pj)
    gamma = max(0.0, W1 - lam * W0)
    mu = (W0 - W1) / group_denominator(E[i0], word.start_node.mode)
    return {
        "W0": W0, "W1": W1, "lambda": lam, "gamma_required": gamma,
        "mu_observed": mu, "i0": i0, "i1": i1,
        "start": li, "end": lj,
    }


def attach_word_times(words: list[base.Word], blocks: list[base.MapBlock]) -> None:
    by_idx = {b.index: b for b in blocks}
    for w in words:
        object.__setattr__(w, "_t0", by_idx[w.start_index].t0)
        object.__setattr__(w, "_t1", by_idx[w.end_index].t1)


def analyze_words(blocks: list[base.MapBlock], trace, E: np.ndarray,
                  mode: str, horizon: float | None,
                  metrics: dict[str, np.ndarray]) -> tuple[dict, list[dict]]:
    if horizon is None or not metrics:
        return {"status": "NOT_AVAILABLE", "word_count": 0}, []
    words = base.compose_words(blocks, mode, float(horizon))
    attach_word_times(words, blocks)
    vals: list[dict] = []
    missing = 0
    for w in words:
        v = word_endpoint_energy(w, trace, E, metrics)
        if v is None:
            missing += 1
        else:
            vals.append(v)
    if not vals:
        return {"status": "NOT_EXERCISED", "word_count": len(words),
                "metric_missing_words": missing}, []
    lam = max(v["lambda"] for v in vals)
    gamma = max(v["gamma_required"] for v in vals)
    b = math.inf if not (lam < 1.0) else gamma / max(EPS, 1.0 - lam)
    return {
        "status": "PASS" if missing == 0 and lam < 1.0 else "FAIL",
        "word_count": len(words), "evaluated_words": len(vals),
        "metric_missing_words": missing,
        "lambda_worst": lam,
        "gamma_required_worst_replay": gamma,
        "replay_invariant_level_b": b if math.isfinite(b) else None,
        "mu_observed_min": min(v["mu_observed"] for v in vals),
        "mu_observed_p05": float(np.quantile([v["mu_observed"] for v in vals], 0.05)),
        "endpoint_W_max": max(max(v["W0"], v["W1"]) for v in vals),
        "qualification": "EXECUTED_NOISY_WORDS_ONLY",
    }, vals


def event_indices(x: np.ndarray) -> np.ndarray:
    return np.flatnonzero(x[1:] > x[:-1]) + 1


def metric_at(trace, E: np.ndarray, k: int,
              metrics: dict[str, dict[str, np.ndarray]]) -> tuple[str, str, float] | None:
    node = source_node_from_trace(trace, k)
    P = metrics[node.mode].get(node.label())
    if P is None:
        return None
    return node.mode, node.label(), group_metric_value(E[k], P, node.mode)


def analyze_handoff_and_jumps(trace, E: np.ndarray,
                              metrics: dict[str, dict[str, np.ndarray]],
                              b_mode: dict[str, float | None]) -> dict:
    live = np.asarray(trace["live"], int)
    bias = np.asarray(trace["bias_active"], int)
    lock = np.asarray(trace["mag_lock"], int)
    refined = np.asarray(trace["mag_refined"], int)
    live_i = event_indices(live)
    if len(live_i) == 0 and len(live) and live[0]:
        live_i = np.array([0])

    handoff = None
    if len(live_i):
        k = int(live_i[0])
        m = metric_at(trace, E, k, metrics)
        handoff = {
            "index": k, "time_s": float(trace["time_s"][k]),
            "theta_deg": math.degrees(float(np.linalg.norm(E[k, :3]))),
            "metric_available": m is not None,
            "mode": m[0] if m else None, "node": m[1] if m else None,
            "W": m[2] if m else None,
            "inside_replay_inner_funnel": bool(m and b_mode.get(m[0]) is not None and m[2] <= b_mode[m[0]]),
        }

    jumps = []
    for name, arr in (("bias_release", bias), ("mag_lock", lock), ("mag_refinement", refined)):
        for kk in event_indices(arr):
            k = int(kk)
            pre = metric_at(trace, E, max(0, k - 1), metrics)
            post = metric_at(trace, E, k, metrics)
            jumps.append({
                "event": name, "index": k, "time_s": float(trace["time_s"][k]),
                "pre_metric_available": pre is not None,
                "post_metric_available": post is not None,
                "pre_mode": pre[0] if pre else None,
                "post_mode": post[0] if post else None,
                "W_pre": pre[2] if pre else None,
                "W_post": post[2] if post else None,
                "jump_ratio": (post[2] / max(EPS, pre[2])) if pre and post else None,
                "inside_destination_replay_inner_funnel": bool(
                    post and b_mode.get(post[0]) is not None and post[2] <= b_mode[post[0]]
                ),
            })
    return {"handoff": handoff, "jumps": jumps}


def recurrence_capture(c0: float | None, lam: float | None,
                       gamma: float | None, b: float | None) -> dict:
    if c0 is None or lam is None or gamma is None or b is None or not (lam < 1.0):
        return {"status": "NOT_AVAILABLE"}
    target = b * (1.0 + 1e-3) + 1e-12
    c = float(c0)
    if c <= target:
        return {"status": "PASS", "N_H_words": 0, "c0": c0, "b": b}
    for n in range(1, 200001):
        c = lam * c + gamma
        if c <= target:
            return {"status": "PASS", "N_H_words": n, "c0": c0, "b": b,
                    "terminal_c": c}
    return {"status": "FAIL", "reason": "NO_FINITE_CAPTURE_WITHIN_LIMIT",
            "c0": c0, "b": b, "terminal_c": c}


def source_enclosure_contract(report: dict, metrics: dict[str, dict[str, np.ndarray]]) -> dict:
    """Describe the exact machine gate a validated backend must satisfy.

    This object is intentionally data, not a proof.  It makes the promotion
    path deterministic: a future interval/Taylor-model backend must return all
    listed bounds with strict signs for every source-reachable node/word.
    """
    lin = report["exact_linear_source_certificate"]
    required = []
    for mode, key in (("H", "held"), ("A", "active")):
        selected = lin[key]["selected"]
        required.append({
            "mode": mode,
            "horizon_s": selected.get("horizon_s"),
            "metric_nodes": sorted(metrics[mode]),
            "requirements": [
                "source_complete_word_family",
                "robust_path_LMI_lambda_gen_lt_1",
                "finite_prefix_gain",
                "theta_star_lt_pi_and_sector_positive",
                "nonlinear_word_mu_W_gt_0",
            ],
        })
    return {
        "schema": 1,
        "claim": "DEPLOYMENT_THEOREM_PROMOTION_CONTRACT",
        "required_modes": required,
        "hybrid_requirements": [
            "startup_handoff_into_certified_sublevel",
            "held_to_active_jump_into_destination_sublevel",
            "magnetic_regauge_jump_into_destination_sublevel",
            "tilt_reset_and_cooldown_jump_inequalities",
        ],
        "stochastic_requirements": [
            "source_uniform_Sigma_bar", "source_uniform_b_W", "source_uniform_v_W",
            "localized_Gaussian_and_Freedman_probability_bound",
        ],
        "promotion_rule": "all validated bounds strict; sampled replay is never sufficient",
    }


def markdown(out: dict) -> str:
    lines = [
        "# OU-III certificate completion",
        "",
        f"Exact linear executed-word gate: **{out['exact_linear_gate']}**",
        f"Executed nonlinear/funnel accounting: **{out['executed_replay_completion']}**",
        f"Neighborhood numerical certificate: **{out['neighborhood_numerical_certificate']}**",
        f"Deployment theorem certificate: **{out['deployment_theorem_certificate']}**",
        "",
        "The executed replay stage uses the exact group metric `a_R(1-cos(theta)) + xi'P_xi xi` from the solved path metrics. It does not promote observed decrements to a neighborhood theorem.",
        "",
        "| Sea | H words | A words | handoff W | handoff in inner funnel | jumps |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in out["records"]:
        h = r["modes"]["H"].get("evaluated_words", 0)
        a = r["modes"]["A"].get("evaluated_words", 0)
        ho = r["handoff_and_jumps"].get("handoff") or {}
        lines.append(
            f"| {r['family']} {r['Hs_m']:.2f} | {h} | {a} | {ho.get('W')} | "
            f"{ho.get('inside_replay_inner_funnel')} | {len(r['handoff_and_jumps']['jumps'])} |"
        )
    lines += [
        "",
        "A complete neighborhood/deployment answer is allowed only after the generated `enclosure_contract.json` is discharged by a validated source-cell backend (interval arithmetic, Taylor models, or equivalent).",
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=DEFAULT_CERT)
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    args = ap.parse_args()
    cert_dir = args.certificate_dir.resolve()
    data_dir = args.data_dir.resolve()
    src = json.loads((cert_dir / "certificate.json").read_text())
    metrics = load_metrics(cert_dir / "path_metrics.npz")
    linear_pass = src.get("exact_linear_source_certificate", {}).get("status") == "PASS"

    selected = src.get("exact_linear_source_certificate", {})
    horizons = {
        "H": selected.get("held", {}).get("selected", {}).get("horizon_s"),
        "A": selected.get("active", {}).get("selected", {}).get("horizon_s"),
    }

    records = []
    global_vals: dict[str, list[dict]] = {"H": [], "A": []}
    raw_loaded = []
    for family, hs, name in RECORDS:
        data = (data_dir / name).resolve()
        slug = f"{family.lower().replace('-','_')}_{hs:.2f}".replace(".", "_")
        tr = cert_dir / f"{data.stem}_certificate_trace.csv"
        mp = cert_dir / f"{data.stem}_exact_maps.bin"
        ts = base.output_csv_for(data)
        trace = np.genfromtxt(tr, delimiter=",", names=True, dtype=None, encoding=None)
        E, theta = base.build_error_states(trace, ts)
        blocks, _ = base.load_exact_maps(mp, slug)
        mode_results = {}
        mode_vals = {}
        for mode in ("H", "A"):
            rr, vv = analyze_words(blocks, trace, E, mode, horizons[mode], metrics[mode])
            mode_results[mode] = rr
            mode_vals[mode] = vv
            global_vals[mode].extend(vv)
        raw_loaded.append((family, hs, trace, E, mode_results, mode_vals))

    mode_envelope = {}
    for mode in ("H", "A"):
        vals = global_vals[mode]
        if vals:
            lam = max(v["lambda"] for v in vals)
            gamma = max(v["gamma_required"] for v in vals)
            b = gamma / max(EPS, 1.0 - lam) if lam < 1.0 else None
            mode_envelope[mode] = {"lambda_worst": lam, "gamma_worst_replay": gamma,
                                   "b_replay": b, "word_count": len(vals)}
        else:
            mode_envelope[mode] = {"lambda_worst": None, "gamma_worst_replay": None,
                                   "b_replay": None, "word_count": 0}
    b_mode = {m: mode_envelope[m]["b_replay"] for m in ("H", "A")}

    c0_mode: dict[str, float | None] = {"H": None, "A": None}
    temp = []
    for family, hs, trace, E, mode_results, mode_vals in raw_loaded:
        hj = analyze_handoff_and_jumps(trace, E, metrics, b_mode)
        arrivals = []
        if hj["handoff"] and hj["handoff"].get("W") is not None:
            arrivals.append((hj["handoff"].get("mode"), hj["handoff"]["W"]))
        for j in hj["jumps"]:
            if j.get("W_post") is not None:
                arrivals.append((j.get("post_mode"), j["W_post"]))
        for mode, w in arrivals:
            if mode in c0_mode:
                c0_mode[mode] = w if c0_mode[mode] is None else max(c0_mode[mode], w)
        temp.append((family, hs, mode_results, hj))

    capture = {}
    for mode in ("H", "A"):
        env = mode_envelope[mode]
        capture[mode] = recurrence_capture(c0_mode[mode], env["lambda_worst"],
                                           env["gamma_worst_replay"], env["b_replay"])
        if horizons[mode] is not None and capture[mode].get("N_H_words") is not None:
            capture[mode]["T_H_s"] = float(horizons[mode]) * int(capture[mode]["N_H_words"])

    for family, hs, mode_results, hj in temp:
        records.append({"family": family, "Hs_m": hs, "modes": mode_results,
                        "handoff_and_jumps": hj})

    replay_complete = bool(linear_pass)
    for r in records:
        for mode in ("H", "A"):
            st = r["modes"][mode]["status"]
            if st not in ("PASS", "NOT_EXERCISED"):
                replay_complete = False
        if r["handoff_and_jumps"].get("handoff") and not r["handoff_and_jumps"]["handoff"].get("metric_available"):
            replay_complete = False
        if any(not j.get("post_metric_available") for j in r["handoff_and_jumps"]["jumps"]):
            replay_complete = False

    contract = source_enclosure_contract(src, metrics)
    (cert_dir / "enclosure_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True))

    out = {
        "schema": 1,
        "exact_linear_gate": "PASS" if linear_pass else "FAIL",
        "executed_replay_completion": "PASS" if replay_complete else "FAIL",
        "metric_geometry": "W_i=a_R,i*(1-cos(theta))+xi^T P_xi,i xi",
        "mode_replay_envelopes": mode_envelope,
        "capture_recursions": capture,
        "records": records,
        "neighborhood_numerical_certificate": "NOT_ESTABLISHED",
        "neighborhood_missing": [
            "validated theta_star sector on source cells",
            "validated zero-input nonlinear word infimum mu_W>0",
            "validated hybrid jump inequalities on source sublevels",
            "source-uniform stochastic b_W and v_W",
        ],
        "deployment_theorem_certificate": "NOT_ESTABLISHED",
        "deployment_missing": "discharge enclosure_contract.json over complete continuous source families",
        "interpretation": (
            "PASS at executed_replay_completion means all eight realized trajectories are compatible with the solved exact-map metric/funnel accounting. "
            "It is not a neighborhood stability proof and cannot promote the deployment theorem."
        ),
    }
    (cert_dir / "completion.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    (cert_dir / "completion.md").write_text(markdown(out))
    print(markdown(out))
    # CI fails only on a regression/instrumentation inconsistency. An open
    # theorem obligation is a scientific result, not a workflow malfunction.
    return 1 if linear_pass and not replay_complete else 0


if __name__ == "__main__":
    raise SystemExit(main())
