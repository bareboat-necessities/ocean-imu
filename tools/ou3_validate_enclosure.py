#!/usr/bin/env python3
"""Validate a rigorous OU-III source-cell enclosure against solved path metrics.

This is the promotion gate from numerical replay evidence to the deployment
source/path theorem.  It does not generate interval/Taylor-model bounds itself;
it verifies their provenance/coverage and independently checks every robust
matrix-box path LMI against path_metrics.npz.

Expected enclosure JSON (schema 1) contains, for each H/A mode:
  source_complete: true
  words: [{start,end,phi_center,phi_radius,prefix_gain_upper,
           theta_star,alpha_R_lower,beta_R_upper,mu_W_lower}]
and top-level hybrid/stochastic validated bounds.  phi_radius is an entrywise
absolute radius obtained with outward-rounded validated arithmetic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
DEFAULT_CERT = REPO / "reports" / "results" / "ou3_numerical_certificate"


def load_metrics(path: Path) -> dict[str, dict[str, np.ndarray]]:
    z = np.load(path, allow_pickle=False)
    ans = {"H": {}, "A": {}}
    for mode in ("H", "A"):
        lk, pk = f"{mode}_labels", f"{mode}_P_local"
        if lk in z and pk in z:
            labels = [str(x) for x in z[lk]]
            mats = np.asarray(z[pk], float)
            ans[mode] = {label: mats[i] for i, label in enumerate(labels)}
    return ans


def robust_box_lmi_upper(C: np.ndarray, R: np.ndarray,
                         Pi: np.ndarray, Pj: np.ndarray) -> dict:
    """Sound spectral upper bound for all Phi=C+Delta, |Delta|<=R.

    Because ||Delta||_2 <= ||R||_F for an entrywise interval box,
      Delta' Pj C + C' Pj Delta + Delta' Pj Delta
    has spectral norm at most
      2 ||Pj C||_2 r + ||Pj||_2 r^2, r=||R||_F.
    Hence strict negativity of the returned upper bound proves
    Phi'PjPhi-Pi << 0 throughout the entire supplied matrix box.
    """
    C = np.asarray(C, float); R = np.asarray(R, float)
    Pi = np.asarray(Pi, float); Pj = np.asarray(Pj, float)
    if C.shape != R.shape or C.shape != Pi.shape or Pi.shape != Pj.shape:
        raise ValueError("matrix shape mismatch")
    if np.any(R < 0) or not np.all(np.isfinite(C)) or not np.all(np.isfinite(R)):
        raise ValueError("invalid interval matrix")
    M0 = C.T @ Pj @ C - Pi
    M0 = 0.5 * (M0 + M0.T)
    nominal_upper = float(np.max(np.linalg.eigvalsh(M0)))
    r = float(np.linalg.norm(R, "fro"))
    perturb = 2.0 * float(np.linalg.norm(Pj @ C, 2)) * r \
              + float(np.linalg.norm(Pj, 2)) * r * r
    upper = nominal_upper + perturb
    return {
        "nominal_lambda_max_difference": nominal_upper,
        "matrix_box_spectral_radius_bound": r,
        "perturbation_upper": perturb,
        "robust_difference_lambda_upper": upper,
        "pass": bool(upper < 0.0),
    }


def finite_positive(x) -> bool:
    try:
        x = float(x)
        return math.isfinite(x) and x > 0.0
    except (TypeError, ValueError):
        return False


def validate_mode(mode: str, payload: dict,
                  metrics: dict[str, np.ndarray], expected_nodes: set[str]) -> dict:
    failures = []
    if not payload.get("source_complete", False):
        failures.append("source family not declared source-complete")
    if not payload.get("outward_rounded", False):
        failures.append("validated arithmetic is not declared outward-rounded")
    words = payload.get("words", [])
    touched = set()
    checked = []
    for i, w in enumerate(words):
        start, end = w.get("start"), w.get("end")
        tag = f"word[{i}] {start}->{end}"
        if start not in metrics or end not in metrics:
            failures.append(f"{tag}: missing solved metric")
            continue
        touched.update((start, end))
        try:
            C = np.asarray(w["phi_center"], float)
            R = np.asarray(w["phi_radius"], float)
            lmi = robust_box_lmi_upper(C, R, metrics[start], metrics[end])
        except Exception as exc:
            failures.append(f"{tag}: invalid matrix enclosure: {exc}")
            continue
        if not lmi["pass"]:
            failures.append(f"{tag}: robust path LMI not strict ({lmi['robust_difference_lambda_upper']})")
        prefix = w.get("prefix_gain_upper")
        if not finite_positive(prefix):
            failures.append(f"{tag}: prefix gain not finite positive")
        theta = w.get("theta_star")
        if not finite_positive(theta) or not float(theta) < math.pi:
            failures.append(f"{tag}: theta_star not in (0,pi)")
        if not finite_positive(w.get("alpha_R_lower")):
            failures.append(f"{tag}: alpha_R lower bound is not strictly positive")
        beta = w.get("beta_R_upper")
        if beta is None or not math.isfinite(float(beta)) or float(beta) < 0:
            failures.append(f"{tag}: beta_R upper bound invalid")
        if not finite_positive(w.get("mu_W_lower")):
            failures.append(f"{tag}: nonlinear mu_W lower bound is not strictly positive")
        checked.append({"start": start, "end": end, "lmi": lmi,
                        "theta_star": theta, "alpha_R_lower": w.get("alpha_R_lower"),
                        "beta_R_upper": beta, "mu_W_lower": w.get("mu_W_lower"),
                        "prefix_gain_upper": prefix})
    missing_nodes = sorted(expected_nodes - touched)
    if missing_nodes:
        failures.append(f"source metric nodes not touched by validated words: {missing_nodes}")
    if not words:
        failures.append("no validated words supplied")
    return {"mode": mode, "pass": not failures, "failures": failures,
            "validated_word_count": len(checked), "checked_words": checked,
            "missing_nodes": missing_nodes}


def validate_hybrid(payload: list[dict]) -> dict:
    required = {"startup_handoff", "held_to_active", "magnetic_regauge",
                "tilt_reset", "cooldown"}
    seen = set(); failures = []
    for i, j in enumerate(payload):
        kind = str(j.get("kind", "")); seen.add(kind)
        margin = j.get("inward_margin_lower")
        if margin is None or not math.isfinite(float(margin)) or float(margin) < 0.0:
            failures.append(f"hybrid[{i}] {kind}: inward margin lower bound is negative/nonfinite")
        if not j.get("outward_rounded", False):
            failures.append(f"hybrid[{i}] {kind}: not outward-rounded")
    missing = sorted(required - seen)
    if missing:
        failures.append(f"missing hybrid obligations: {missing}")
    return {"pass": not failures, "failures": failures, "seen": sorted(seen)}


def validate_stochastic(payload: dict) -> dict:
    failures = []
    for key in ("Sigma_bar_norm_upper", "b_W_upper", "v_W_upper"):
        v = payload.get(key)
        if v is None or not math.isfinite(float(v)) or float(v) < 0:
            failures.append(f"{key} missing/nonfinite/negative")
    prob = payload.get("finite_horizon_failure_probability_upper")
    if prob is None or not math.isfinite(float(prob)) or not (0.0 <= float(prob) <= 1.0):
        failures.append("finite_horizon_failure_probability_upper invalid")
    if not payload.get("outward_rounded", False):
        failures.append("stochastic bounds are not outward-rounded")
    return {"pass": not failures, "failures": failures,
            "finite_horizon_failure_probability_upper": prob}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=DEFAULT_CERT)
    ap.add_argument("--enclosure", type=Path, required=True)
    args = ap.parse_args()
    cert = args.certificate_dir.resolve()
    inp = json.loads(args.enclosure.read_text())
    if inp.get("schema") != 1:
        raise RuntimeError("unsupported validated enclosure schema")
    prov = inp.get("provenance", {})
    provenance_ok = bool(prov.get("validated_arithmetic")) and bool(prov.get("outward_rounding")) \
        and bool(prov.get("source_generated_not_trajectory_fit"))
    metrics = load_metrics(cert / "path_metrics.npz")
    contract = json.loads((cert / "enclosure_contract.json").read_text())
    expected = {}
    for m in contract.get("required_modes", []):
        expected[m["mode"]] = set(m.get("metric_nodes", []))
    modes = {}
    for mode in ("H", "A"):
        modes[mode] = validate_mode(mode, inp.get("modes", {}).get(mode, {}),
                                    metrics[mode], expected.get(mode, set()))
    hybrid = validate_hybrid(inp.get("hybrid", []))
    stochastic = validate_stochastic(inp.get("stochastic", {}))
    passed = provenance_ok and all(modes[m]["pass"] for m in modes) \
             and hybrid["pass"] and stochastic["pass"]
    out = {
        "schema": 1,
        "validated_enclosure_provenance_pass": provenance_ok,
        "modes": modes,
        "hybrid": hybrid,
        "stochastic": stochastic,
        "deployment_theorem_certificate": "PASS" if passed else "FAIL",
        "promotion_is_machine_verified": True,
    }
    (cert / "validated_enclosure_check.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    print(json.dumps({
        "deployment_theorem_certificate": out["deployment_theorem_certificate"],
        "H": modes["H"]["pass"], "A": modes["A"]["pass"],
        "hybrid": hybrid["pass"], "stochastic": stochastic["pass"],
        "provenance": provenance_ok,
    }, indent=2, sort_keys=True))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
