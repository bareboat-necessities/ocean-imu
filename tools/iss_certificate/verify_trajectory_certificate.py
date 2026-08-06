#!/usr/bin/env python3
"""Linear-time-varying stability certificate from trajectory transition matrices.

This is the third and least approximate of the three certificates in this
directory:

  verify_certificate.py            freezes the operating point and settles the
                                   covariance there. Measures a system that
                                   never moves; rho = 1 - 7.3e-7.

  verify_switched_certificate.py   multiplies frozen-point Jacobians along an
                                   attitude orbit. Captures the mechanism, but
                                   each factor still carries a covariance
                                   settled at its own point and held for a
                                   dwell interval.

  this one                         per-step Jacobians taken along a continuous
                                   wave trajectory against the filter's actual
                                   running covariance, accumulated into window
                                   transition matrices Phi(k+N, k).

Given those windows it reports the LTV contraction, the within-window
transient, and searches for one quadratic Lyapunov function valid for every
window. A common metric over the windows gives a genuine linear-time-varying
statement -- V decreases across every window along the trajectory -- with no
appeal to the slowly varying assumption A4.

Input comes from generate_trajectory_jacobians.cpp. Requires cvxpy for the
common-metric search; without it the spectral summary is still printed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

# Physical scale per 3-block, so the LMI is not conditioned by unit choice:
# dtheta [rad], bg [rad/s], v [m/s], p [m], S [m s], aw [m/s^2].
STATE_SCALE = np.repeat([1e-3, 1e-3, 1.0, 1.0, 30.0, 1.0], 3)


def read_windows(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().split()
        if len(header) != 6 or header[0] != "OCEAN_IMU_ISS_TRAJECTORY_V1":
            raise ValueError("invalid trajectory file header")
        n = int(header[1])
        count = int(header[2])
        period = float(header[3])
        dt = float(header[4])
        steps = int(header[5])

        times, peaks, mats = [], [], []
        for expected in range(count):
            fields = fh.readline().split()
            if len(fields) != 4 or fields[0] != "WINDOW" or int(fields[1]) != expected:
                raise ValueError(f"invalid window record {expected}")
            times.append(float(fields[2]))
            peaks.append(float(fields[3]))
            rows = [[float(x) for x in fh.readline().split()] for _ in range(n)]
            mats.append(np.asarray(rows, dtype=float))
    return {"period": period, "dt": dt, "steps": steps,
            "times": times, "peaks": peaks, "matrices": mats}


def common_metric(mats, cp, bisection_steps: int = 18):
    """Least gamma with a common L >= I satisfying Phi' L Phi <= gamma L."""
    n = mats[0].shape[0]

    def feasible(gamma: float):
        L = cp.Variable((n, n), PSD=True)
        constraints = [L >> np.eye(n)] + [P.T @ L @ P << gamma * L for P in mats]
        problem = cp.Problem(cp.Minimize(cp.trace(L)), constraints)
        try:
            problem.solve(solver=cp.CLARABEL)
        except Exception:
            return None
        if L.value is None:
            return None
        Lv = 0.5 * (L.value + L.value.T)
        tol = 1e-6 * float(np.linalg.norm(Lv, 2))
        worst = max(float(np.linalg.eigvalsh(P.T @ Lv @ P - gamma * Lv)[-1]) for P in mats)
        return Lv if worst <= tol else None

    lo, hi, best = 0.0, 1.0, None
    for _ in range(bisection_steps):
        mid = 0.5 * (lo + hi)
        L = feasible(mid)
        if L is None:
            lo = mid
        else:
            hi, best = mid, L
    return hi, best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory", type=Path)
    ap.add_argument("--skip", type=int, default=0,
                    help="discard this many leading windows as covariance settling")
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()

    data = read_windows(args.trajectory)
    period, dt = data["period"], data["dt"]
    scale = np.diag(STATE_SCALE)
    scale_inv = np.diag(1.0 / STATE_SCALE)

    mats = [scale_inv @ P @ scale for P in data["matrices"]][args.skip:]
    times = data["times"][args.skip:]
    peaks = data["peaks"][args.skip:]
    if not mats:
        print("no windows left after --skip")
        return 1

    rho = [float(abs(np.linalg.eigvals(P)).max()) for P in mats]
    norms = [float(np.linalg.norm(P, 2)) for P in mats]

    print(f"{len(mats)} windows of {period:g} s at dt={dt:g}, t = {times[0]:g}..{times[-1]:g} s")
    print(f"  rho(Phi) per window : first {rho[0]:.4f}  last {rho[-1]:.4f}  "
          f"max {max(rho):.4f}")
    print(f"  ||Phi||_2 per window: max {max(norms):.4f}")
    print(f"  within-window peak  : max {max(peaks):.1f}  <- transient amplification")

    worst = max(rho)
    if worst < 1.0:
        print(f"  spectral time constant at the worst window: "
              f"{period / -np.log(worst):.1f} s")

    result: dict[str, Any] = {
        "schema": "ocean-imu-trajectory-certificate-v1",
        "scope": "LTV: per-step Jacobians along a continuous wave trajectory",
        "period_s": period,
        "dt": dt,
        "windows": len(mats),
        "windows_skipped": args.skip,
        "rho_first": rho[0],
        "rho_last": rho[-1],
        "rho_max": worst,
        "phi_norm_max": max(norms),
        "within_window_peak_max": max(peaks),
    }

    try:
        import cvxpy as cp
    except ImportError:
        print("\ncvxpy not installed; skipping the common-metric search")
        cp = None

    if cp is not None:
        gamma, L = common_metric(mats, cp)
        if L is None:
            print("\nno common metric over the trajectory windows")
            result["common_metric_found"] = False
        else:
            ev = np.linalg.eigvalsh(L)
            c1 = float(np.sqrt(ev[-1] / ev[0]))
            per_step = float(gamma ** (1.0 / data["steps"])) ** 0.5
            print(f"\ncommon metric over every window (so delta_L = 0):")
            print(f"  V decays by gamma = {gamma:.5f} per {period:g} s window")
            print(f"  ||e|| decays by {np.sqrt(gamma):.4f} per window")
            print(f"  c1 = sqrt(cond L) = {c1:.1f}")
            print(f"  equivalent e-folding = {dt / (1.0 - per_step):.2f} s")
            result.update({
                "common_metric_found": True,
                "common_metric_gamma_per_window": gamma,
                "common_metric_c1": c1,
                "efold_seconds": dt / (1.0 - per_step),
            })

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
