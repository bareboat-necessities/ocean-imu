#!/usr/bin/env python3
"""Source-derived design probe for the P4 complete-word translation route.

This is deliberately NOT a theorem producer.  It answers one narrow design
question before the validated interval backend is built: does accumulating the
entire word covariance decomposition change the pathological translation margin
by orders of magnitude, or is the proposed route another dead end?

For the P3 worst H/A source cell, propagate the exact 4-state integrated-OU
translation process [v,p,S,a_w] for several finite horizons.  At every sample
apply *more* information than shipping is guaranteed to apply: an S=0 update
and an accelerometer observation of a_w.  For a covariance component this is a
conservative direction for the eventual lower-bound proof because, for any
implemented gain K,

  (I-KH) Omega (I-KH)' + K R K' >= (Omega^-1 + H'R^-1H)^-1.

The probe uses ordinary numpy/eigenvalue arithmetic and a small source grid, so
its output can select/reject an architecture but can never promote P4.  A useful
result must subsequently be reproduced with outward-rounded source-cell
matrices and adaptive subdivision.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from ou3_interval import Interval
import ou3_p4_frontier_margin_diagnostic as DIAG
import ou3_source_reachable_matrix_p3 as P3BASE
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
HORIZONS_S = (1.0, 2.0, 4.0, 8.0)


def _mid_interval_matrix(A):
    return np.asarray([[(x.lo + x.hi) * 0.5 for x in row] for row in A], dtype=float)


def _transition(tau: float, h: float) -> np.ndarray:
    x = h / tau
    a = math.exp(-x)
    em1 = math.expm1(-x)
    phi_va = -tau * em1
    if abs(x) < 1.0e-2:
        x2 = x*x; x3=x2*x; x4=x3*x; x5=x4*x
        phi_pa = tau*tau*(0.5*x2-x3/6.0+x4/24.0)
        phi_Sa = tau**3*(x3/6.0-x4/24.0+x5/120.0)
    else:
        phi_pa = tau*tau*(x+em1)
        phi_Sa = tau**3*(0.5*x*x-x-em1)
    return np.asarray([
        [1.0, 0.0, 0.0, phi_va],
        [h,   1.0, 0.0, phi_pa],
        [0.5*h*h, h, 1.0, phi_Sa],
        [0.0, 0.0, 0.0, a],
    ], dtype=float)


def _process_Q(tau: float, sigma: float, h: float) -> np.ndarray:
    x = h / tau
    qbar = _mid_interval_matrix(P3BASE.qbar_integrated_ou(Interval.outward_bounds(x, x)))
    D = np.diag([sigma*tau, sigma*tau*tau, sigma*tau**3, sigma])
    Q = D @ qbar @ D.T
    return 0.5*(Q+Q.T)


def _scalar_optimal_update(P: np.ndarray, state_index: int, variance: float) -> np.ndarray:
    c = P[:, state_index].copy()
    den = float(P[state_index, state_index] + variance)
    if not (math.isfinite(den) and den > 0.0):
        raise RuntimeError("nonpositive scalar innovation variance in design probe")
    out = P - np.outer(c, c) / den
    return 0.5*(out+out.T)


def _one_point(tau: float, sigma: float, rs: float, h: float,
               horizon: float, upper: np.ndarray, racc_var: float) -> dict:
    A = _transition(tau, h)
    Q = _process_Q(tau, sigma, h)
    n = max(1, int(math.ceil(horizon/h)))
    Omega = np.zeros((4,4), dtype=float)
    min_eig = math.inf
    for _ in range(n):
        Omega = A @ Omega @ A.T + Q
        # Deliberately over-inform the comparison: shipping cannot perform more
        # than one S correction and one accelerometer correction per IMU sample.
        Omega = _scalar_optimal_update(Omega, 2, rs*rs)
        Omega = _scalar_optimal_update(Omega, 3, racc_var)
        ev = np.linalg.eigvalsh(Omega)
        min_eig = min(min_eig, float(ev[0]))
    scale = np.diag(1.0/np.sqrt(upper))
    G = scale @ Omega @ scale
    delta = float(np.min(np.linalg.eigvalsh(0.5*(G+G.T))))
    return {
        "tau_s": tau,
        "sigma_aw_mps2": sigma,
        "R_S_std": rs,
        "horizon_s": horizon,
        "steps": n,
        "translation_complete_word_generalized_margin_design": delta,
        "Omega_lambda_min_design": float(np.min(np.linalg.eigvalsh(Omega))),
        "Omega_trace_design": float(np.trace(Omega)),
        "minimum_intermediate_Omega_lambda_min_design": min_eig,
    }


def _points(bounds):
    lo, hi = map(float, bounds)
    return (lo, math.sqrt(lo*hi), hi) if lo > 0.0 else (lo, 0.5*(lo+hi), hi)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    p = Path(domain_path).resolve()
    d = DIAG.build(p)
    vf = VECTOR.build()
    failures = [f"diagnostic: {x}" for x in DIAG.validate(d)] + [f"vector: {x}" for x in VECTOR.validate(vf)]
    racc = float(vf["configured_measurement_bounds"]["acc_measurement_std_mps2"])
    h = float(P3BASE.source_schedule()["dt_s"])
    modes = {}
    for mode in ("H", "A"):
        m = d["modes"][mode]
        w = m["p3_worst_cell"]
        xlo, xhi = map(float, w["x_h_over_tau"])
        tau_bounds = [h/xhi, h/xlo]
        # Rebuild the actual P3 worst comparison to obtain directional upper
        # variances.  H and A share the first 18 coordinate order.
        p3 = DIAG.P3.build(p)
        upper_all = list(map(float, p3["modes"][mode]["matrix_comparison"]["Sigma_diagonal_upper"]))
        upper = np.asarray([upper_all[6], upper_all[9], upper_all[12], upper_all[15]], dtype=float)
        rows=[]
        for horizon in HORIZONS_S:
            for tau, sigma, rs in itertools.product(
                _points(tau_bounds), _points(w["sigma_aw"]), _points(w["R_S"])
            ):
                rows.append(_one_point(tau, sigma, rs, h, horizon, upper, racc*racc))
        by_horizon={}
        for horizon in HORIZONS_S:
            rr=[r for r in rows if r["horizon_s"]==horizon]
            worst=min(rr, key=lambda r:r["translation_complete_word_generalized_margin_design"])
            best=max(rr, key=lambda r:r["translation_complete_word_generalized_margin_design"])
            old=float(w["delta_translation_lower"])
            by_horizon[str(horizon)]={
                "worst_grid_point": worst,
                "best_grid_point": best,
                "old_single_seed_translation_margin_lower": old,
                "design_worst_to_old_margin_ratio": worst["translation_complete_word_generalized_margin_design"]/old,
            }
        modes[mode]={
            "source_cell": w,
            "tau_s_derived_from_x_cell": tau_bounds,
            "translation_directional_covariance_upper": upper.tolist(),
            "grid_points_per_horizon": 27,
            "horizons": by_horizon,
        }
    return {
        "qualification":"OU3_P4_TRANSLATION_FULL_WORD_DESIGN_PROBE",
        "source_parameters_from_theorem_domain":True,
        "trajectory_replay_used":False,
        "ordinary_floating_point_design_only":True,
        "validated_for_theorem_promotion":False,
        "P4_USABLE_CERTIFICATE_STATUS":"NOT_ESTABLISHED",
        "purpose":"reject or select complete-word accumulation before implementing interval proof",
        "modes":modes,
        "failures":failures,
    }


def validate(d: dict) -> list[str]:
    f=list(d.get("failures",[]))
    if d.get("trajectory_replay_used") is not False: f.append("design probe uses replay")
    if d.get("ordinary_floating_point_design_only") is not True: f.append("design-only qualification missing")
    if d.get("validated_for_theorem_promotion") is not False: f.append("unvalidated design probe promoted")
    if d.get("P4_USABLE_CERTIFICATE_STATUS") != "NOT_ESTABLISHED": f.append("design probe promoted P4")
    for mode in ("H","A"):
        hs=d.get("modes",{}).get(mode,{}).get("horizons",{})
        for h in HORIZONS_S:
            r=hs.get(str(h),{})
            q=r.get("worst_grid_point",{}).get("translation_complete_word_generalized_margin_design")
            if not isinstance(q,(int,float)) or not math.isfinite(float(q)) or float(q)<=0.0:
                f.append(f"{mode}: horizon {h} lost positive design margin")
    return f


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args(); d=build(a.domain); f=validate(d); d["validation_failures"]=f
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True))
    print(json.dumps({"status":d["P4_USABLE_CERTIFICATE_STATUS"],"modes":{
        mode:{h:d["modes"][mode]["horizons"][h]["design_worst_to_old_margin_ratio"] for h in d["modes"][mode]["horizons"]}
        for mode in ("H","A")},"failures":f},indent=2,sort_keys=True))
    return 0 if not f else 2

if __name__=="__main__": raise SystemExit(main())
