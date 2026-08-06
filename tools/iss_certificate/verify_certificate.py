#!/usr/bin/env python3
"""Verify reduced-state local stability of the implemented OU-III cycle.

The implementation Jacobian is partitioned as

    [ e+ ]   [F_ee F_eb] [ e ]
    [ b+ ] = [F_be F_bb] [ b ],

where e is the 18-dimensional performance error and b is the three-dimensional
accelerometer-bias nuisance error.  The certificate verifies every F_ee is
strictly Schur stable, constructs its discrete Lyapunov matrix, independently
checks the Lyapunov residual, and reports the F_eb input-channel norm.

It also audits the two things that decide whether Theorem 1 of
doc/kalman_ou_iii/w3d-iss-stability.tex-part says anything about the running
system, rather than only about each frozen operating point:

  * the ISS constants the theorem actually delivers -- c1 = sqrt(l_bar/l_min)
    and the per-step contraction -- next to the true transient amplification
    sup_k ||F^k||, which is the physically meaningful number;

  * assumption A4.  The frozen-point argument needs the Lyapunov metric to
    drift by less than delta_L < 1/(l_bar - 1) per step.  This measures the
    metric ratio between adjacent envelope points and reports how far the
    real system is from that budget.  It is reported, not enforced: a failure
    means the certificate is a necessary-condition screen on the frozen
    dynamics and not a closed-loop stability proof.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import solve_discrete_lyapunov

PERFORMANCE_IDX = np.arange(18, dtype=int)
ACCEL_BIAS_IDX = np.arange(18, 21, dtype=int)

POINT_FIELDS = ["dt", "tau", "sigma_aw", "sigma_s", "roll", "pitch", "wx", "wy", "wz"]

# Physical scale of each 3-block, used only for reporting a conditioning number
# that is not dominated by the choice of units.  dtheta [rad], bg [rad/s],
# v [m/s], p [m], S [m s], aw [m/s^2].
STATE_SCALE = np.repeat([1e-3, 1e-3, 1.0, 1.0, 30.0, 1.0], 3)


def transient_amplification(A: np.ndarray, horizon: int = 200000) -> float:
    """sup_k ||A^k||_2, the amplification a trajectory can actually see."""
    peak = 1.0
    P = np.eye(A.shape[0])
    probe = 1
    for k in range(1, horizon + 1):
        P = P @ A
        if k == probe:
            peak = max(peak, float(np.linalg.norm(P, ord=2)))
            probe = max(probe + 1, int(probe * 1.3))
    return peak


def metric_drift_audit(metas, lyapunov):
    """Largest single-grid-step metric ratio, per swept parameter.

    delta_L is defined by L_j <= (1 + delta_L) L_i, i.e. the largest
    generalized eigenvalue of (L_j, L_i) minus one.
    """
    index = {tuple(m[f] for f in POINT_FIELDS): i for i, m in enumerate(metas)}
    axis_values = [
        sorted({m[f] for m in metas}) for f in POINT_FIELDS
    ]
    worst: dict[str, float] = {}
    pairs = 0
    for key, i in index.items():
        for axis, field in enumerate(POINT_FIELDS):
            values = axis_values[axis]
            position = values.index(key[axis])
            if position + 1 >= len(values):
                continue
            neighbour = list(key)
            neighbour[axis] = values[position + 1]
            j = index.get(tuple(neighbour))
            if j is None:
                continue
            pairs += 1
            ratio = float(np.max(np.abs(np.linalg.eigvals(
                np.linalg.solve(lyapunov[i], lyapunov[j])
            )))) - 1.0
            worst[field] = max(worst.get(field, -np.inf), ratio)
    return pairs, worst


def read_matrices(path: Path) -> tuple[list[dict[str, float]], list[np.ndarray]]:
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().split()
        if len(header) != 3 or header[0] != "OCEAN_IMU_ISS_JACOBIANS_V1":
            raise ValueError("invalid Jacobian file header")
        n, count = int(header[1]), int(header[2])
        meta: list[dict[str, float]] = []
        matrices: list[np.ndarray] = []
        for expected in range(count):
            fields = fh.readline().split()
            if len(fields) != 11 or fields[0] != "POINT" or int(fields[1]) != expected:
                raise ValueError(f"invalid point record {expected}")
            names = ["dt", "tau", "sigma_aw", "sigma_s", "roll", "pitch", "wx", "wy", "wz"]
            meta.append(dict(zip(names, map(float, fields[2:]), strict=True)))
            rows = [[float(x) for x in fh.readline().split()] for _ in range(n)]
            F = np.asarray(rows, dtype=float)
            if F.shape != (n, n) or not np.isfinite(F).all():
                raise ValueError(f"invalid matrix {expected}")
            matrices.append(F)
    return meta, matrices


def lyapunov_metrics(A: np.ndarray) -> dict[str, float]:
    n = A.shape[0]
    P = solve_discrete_lyapunov(A.T, np.eye(n))
    P = 0.5 * (P + P.T)
    eig = np.linalg.eigvalsh(P)
    if eig[0] <= 0.0 or not np.isfinite(eig).all():
        raise RuntimeError("discrete Lyapunov solution is not positive definite")

    residual = A.T @ P @ A - P + np.eye(n)
    residual_norm = float(np.linalg.norm(residual, ord=2))
    residual_relative = residual_norm / max(1.0, float(eig[-1]))

    # From A' P A = P - I, the squared contraction factor in the P metric is
    # lambda_max(P^{-1/2}(P-I)P^{-1/2}) = 1 - 1/lambda_max(P).
    lambda_metric = max(0.0, 1.0 - 1.0 / float(eig[-1]))
    return {
        "L": P,
        "p_eigen_min": float(eig[0]),
        "p_eigen_max": float(eig[-1]),
        "p_condition": float(eig[-1] / eig[0]),
        "lambda_metric": lambda_metric,
        "rho_metric": math.sqrt(lambda_metric),
        "residual_norm": residual_norm,
        "residual_relative": residual_relative,
    }


def write_tex(path: Path, result: dict[str, Any]) -> None:
    worst = result["worst_spectral_point"]
    text = rf"""% Generated by tools/iss_certificate/verify_certificate.py; do not edit.
\begin{{table}}[t]
\centering
\caption{{Implementation-generated reduced-state stability certificate for the OU--III performance error. Accelerometer bias is excluded from the certified error and retained as a bounded input channel. The lower two blocks report the constants Theorem~\ref{{thm:adaptive-local-iss}} actually delivers and the audit of assumption A4; see Sec.~\ref{{app:iss-limits}} for how to read them.}}
\label{{tab:iss-certificate}}
\begin{{tabular}}{{ll}}
\toprule
Quantity & Verified value \\
\midrule
Performance-state dimension & {result['state_dimension']} \\
Excluded nuisance-state dimension & {result['excluded_state_dimension']} \\
Generated operating points & {result['matrix_count']} \\
Maximum spectral radius of $F_{{ee}}$ & {result['spectral_radius_max']:.9f} \\
Minimum Schur margin $1-\rho(F_{{ee}})$ & {result['spectral_margin_min']:.3e} \\
Worst pointwise Lyapunov $\rho_L$ & {result['rho_metric_max']:.9f} \\
Maximum relative Lyapunov residual & {result['lyapunov_residual_relative_max']:.3e} \\
Worst $\|F_{{e b_a}}\|_2$ & {result['max_bias_coupling_2norm']:.3e} \\
\midrule
\multicolumn{{2}}{{l}}{{\emph{{Constants the ISS bound delivers}}}} \\
$c_1=\sqrt{{\overline\ell/\underline\ell}}$ & {result['iss_c1']:.3e} \\
$c_1$ in scaled coordinates & {result['iss_c1_scaled_coordinates']:.3e} \\
Certified $e$-folding time & {result['iss_efold_seconds']:.3e} s \\
True $\sup_k\|F^k\|_2$ at the worst point & {result['transient_amplification_worst_point']:.3e} \\
\midrule
\multicolumn{{2}}{{l}}{{\emph{{Assumption A4 audit}}}} \\
Budget $\delta_L<1/(\overline\ell-1)$ & {result['a4_delta_l_budget_per_step']:.3e} \\
Worst measured $\delta_L$ over one grid step & {result['a4_delta_l_worst_adjacent_step']:.3e} \\
A4 satisfied on the envelope & {'yes' if result['a4_satisfied_by_envelope'] else 'no'} \\
Worst point $(\Delta t,\tau,\sigma_a,\sigma_S)$ &
({worst['dt']:.6f}, {worst['tau']:.3f}, {worst['sigma_aw']:.3f}, {worst['sigma_s']:.3f}) \\
Worst point $(\phi,\theta,\omega_x,\omega_y,\omega_z)$ &
({worst['roll']:.3f}, {worst['pitch']:.3f}, {worst['wx']:.3f}, {worst['wy']:.3f}, {worst['wz']:.3f}) \\
\bottomrule
\end{{tabular}}
\end{{table}}
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("jacobians", type=Path)
    ap.add_argument("--json", type=Path, required=True)
    ap.add_argument("--tex", type=Path, required=True)
    args = ap.parse_args()

    meta, full_matrices = read_matrices(args.jacobians)
    performance_matrices = [F[np.ix_(PERFORMANCE_IDX, PERFORMANCE_IDX)] for F in full_matrices]
    bias_channels = [F[np.ix_(PERFORMANCE_IDX, ACCEL_BIAS_IDX)] for F in full_matrices]

    spectral_radii = np.asarray([
        float(np.max(np.abs(np.linalg.eigvals(A)))) for A in performance_matrices
    ])
    worst_spectral_index = int(np.argmax(spectral_radii))
    spectral_radius_max = float(spectral_radii[worst_spectral_index])
    spectral_margin_min = 1.0 - spectral_radius_max
    if not math.isfinite(spectral_radius_max) or spectral_margin_min <= 1.0e-8:
        raise RuntimeError(
            "reduced performance-state Jacobian is not strictly Schur stable: "
            f"rho={spectral_radius_max:.12g}, margin={spectral_margin_min:.12g}"
        )

    metrics = [lyapunov_metrics(A) for A in performance_matrices]
    residual_relative_max = max(m["residual_relative"] for m in metrics)
    if residual_relative_max > 1.0e-10:
        raise RuntimeError(
            "pointwise Lyapunov residual check failed: "
            f"max_relative_residual={residual_relative_max:.12g}"
        )

    # ---- ISS constants the theorem actually delivers -------------------
    l_bar = max(m["p_eigen_max"] for m in metrics)
    l_min = min(m["p_eigen_min"] for m in metrics)
    c1 = math.sqrt(l_bar / l_min)
    rho_iss = math.sqrt(1.0 - 1.0 / l_bar)          # per-step decay of ||e||
    dt_min = min(m["dt"] for m in meta)
    efold_s = dt_min / (1.0 - rho_iss) if rho_iss < 1.0 else math.inf

    # Same conditioning measured in physically scaled coordinates, so the
    # number is not just a statement about the choice of units.
    scale = np.diag(STATE_SCALE)
    scale_inv = np.diag(1.0 / STATE_SCALE)
    scaled = [lyapunov_metrics(scale_inv @ A @ scale) for A in performance_matrices]
    c1_scaled = math.sqrt(
        max(m["p_eigen_max"] for m in scaled) / min(m["p_eigen_min"] for m in scaled)
    )

    worst_A = performance_matrices[worst_spectral_index]
    amplification = transient_amplification(worst_A)

    # ---- A4 audit -------------------------------------------------------
    delta_l_budget = 1.0 / (l_bar - 1.0)
    pair_count, drift = metric_drift_audit(meta, [m["L"] for m in metrics])
    drift_worst = max(drift.values())

    bias_norms = np.asarray([np.linalg.norm(B, ord=2) for B in bias_channels])
    max_bias_index = int(np.argmax(bias_norms))
    rho_metric_values = np.asarray([m["rho_metric"] for m in metrics])
    worst_metric_index = int(np.argmax(rho_metric_values))

    result: dict[str, Any] = {
        "schema": "ocean-imu-iss-certificate-v3",
        "scope": "finite implementation-generated reduced performance-state envelope",
        "certificate_type": "pointwise discrete Lyapunov family",
        "state_definition": ["dtheta", "bg", "v", "p", "S", "aw"],
        "excluded_nuisance_state": ["ba"],
        "state_dimension": 18,
        "excluded_state_dimension": 3,
        "matrix_count": len(performance_matrices),
        "spectral_radius_min": float(np.min(spectral_radii)),
        "spectral_radius_max": spectral_radius_max,
        "spectral_margin_min": spectral_margin_min,
        "worst_spectral_index": worst_spectral_index,
        "worst_spectral_point": meta[worst_spectral_index],
        "rho_metric_max": float(rho_metric_values[worst_metric_index]),
        "worst_metric_index": worst_metric_index,
        "worst_metric_point": meta[worst_metric_index],
        "p_eigen_min_global": min(m["p_eigen_min"] for m in metrics),
        "p_eigen_max_global": max(m["p_eigen_max"] for m in metrics),
        "p_condition_max": max(m["p_condition"] for m in metrics),
        "lyapunov_residual_relative_max": residual_relative_max,
        "max_bias_coupling_2norm": float(bias_norms[max_bias_index]),
        "max_bias_coupling_index": max_bias_index,
        "max_bias_coupling_point": meta[max_bias_index],

        "iss_c1": c1,
        "iss_c1_scaled_coordinates": c1_scaled,
        "iss_rho_per_step": rho_iss,
        "iss_efold_seconds": efold_s,
        "transient_amplification_worst_point": amplification,

        "a4_delta_l_budget_per_step": delta_l_budget,
        "a4_adjacent_pairs": pair_count,
        "a4_delta_l_by_parameter": drift,
        "a4_delta_l_worst_adjacent_step": drift_worst,
        "a4_satisfied_by_envelope": bool(drift_worst < delta_l_budget),
    }

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.tex.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_tex(args.tex, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
