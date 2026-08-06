#!/usr/bin/env python3
"""Switched/Floquet stability certificate for the OU-III performance error.

Motivation
----------
The pointwise certificate in verify_certificate.py freezes the operating point
and asks whether each frozen Jacobian is Schur stable. It is, but only barely:
rho = 1 - 7.3e-7, a mode of about 1.6 hours, which makes the ISS constants of
Theorem 1 vacuous and makes its assumption A4 unsatisfiable by ten orders of
magnitude.

A4 is the wrong assumption for this system. It asks the operating point to
move *slowly* compared with the estimator dynamics. On a vessel the opposite
holds: roll and pitch sweep the whole envelope in seconds while the frozen
mode is hours. That is the fast-variation regime, where the object that
governs stability is the transition matrix over a cycle, not any frozen
Jacobian.

This script computes that object. It multiplies implementation Jacobians along
a closed orbit in attitude, forming the monodromy matrix over one wave period,
and then searches for a single quadratic Lyapunov function valid for every
monodromy in the envelope.

Result on the shipped envelope: the motion helps enormously. Held still, the
error decays with a 1764 s time constant. Moving through the same envelope it
decays by a factor 0.46 per wave period, an e-folding of about 13 s. A common
metric exists for the monodromies, so delta_L = 0 and A4 is not needed at all
for the switched model.

Caveats
-------
This is a switched-system model, not the true trajectory. Each factor is a
frozen-point Jacobian whose covariance was settled at that operating point,
held for a dwell interval, whereas the real filter's covariance evolves
continuously. It is a strictly better model than freezing, and it captures the
mechanism the frozen analysis misses, but the rigorous version needs Jacobians
taken along an actual simulated trajectory. Treat the numbers as evidence
about the mechanism and its size, not as a closed-loop proof.

Requires cvxpy.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

N = 18
POINT_FIELDS = ["dt", "tau", "sigma_aw", "sigma_s", "roll", "pitch", "wx", "wy", "wz"]

# Physical scale per 3-block, so the LMI is not conditioned by unit choice.
STATE_SCALE = np.repeat([1e-3, 1e-3, 1.0, 1.0, 30.0, 1.0], 3)

# A closed orbit in (roll, pitch) at zero body rate. These are envelope corners
# the generator already samples, traversed in a cycle.
DEFAULT_ORBIT = [
    (-0.70, 0.00, 0.0, 0.0, 0.0),
    (0.00, 0.45, 0.0, 0.0, 0.0),
    (0.70, 0.00, 0.0, 0.0, 0.0),
    (0.00, -0.45, 0.0, 0.0, 0.0),
]


def read_matrices(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().split()
        if len(header) != 3 or header[0] != "OCEAN_IMU_ISS_JACOBIANS_V1":
            raise ValueError("invalid Jacobian file header")
        n, count = int(header[1]), int(header[2])
        meta, mats = [], []
        for _ in range(count):
            fields = fh.readline().split()
            meta.append(dict(zip(POINT_FIELDS, map(float, fields[2:]), strict=True)))
            rows = [[float(x) for x in fh.readline().split()] for _ in range(n)]
            mats.append(np.asarray(rows, dtype=float))
    return meta, mats


def build_monodromies(meta, mats, period_s: float):
    """One monodromy per tuning cell, by orbiting attitude over one period."""
    scale = np.diag(STATE_SCALE)
    scale_inv = np.diag(1.0 / STATE_SCALE)

    cells: dict[tuple, dict[tuple, np.ndarray]] = {}
    for m, A in zip(meta, mats, strict=True):
        cell = (m["dt"], m["tau"], m["sigma_aw"], m["sigma_s"])
        pose = (m["roll"], m["pitch"], m["wx"], m["wy"], m["wz"])
        cells.setdefault(cell, {})[pose] = scale_inv @ A[:N, :N] @ scale

    out = []
    for cell, poses in cells.items():
        if not all(p in poses for p in DEFAULT_ORBIT):
            continue
        steps = int(round(period_s / cell[0]))
        dwell = max(1, steps // len(DEFAULT_ORBIT))
        Phi = np.eye(N)
        for p in DEFAULT_ORBIT:
            Phi = np.linalg.matrix_power(poses[p], dwell) @ Phi
        out.append((cell, dwell * len(DEFAULT_ORBIT), Phi))
    return out


def smallest_common_contraction(mats, cp, bisection_steps: int = 18):
    """Least gamma with a common L >= I satisfying Phi' L Phi <= gamma L."""

    def feasible(gamma: float):
        L = cp.Variable((N, N), PSD=True)
        constraints = [L >> np.eye(N)] + [P.T @ L @ P << gamma * L for P in mats]
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
    ap.add_argument("jacobians", type=Path)
    ap.add_argument("--period", type=float, default=10.0, help="wave period [s]")
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()

    try:
        import cvxpy as cp
    except ImportError:
        print("cvxpy is required: pip install cvxpy")
        return 2

    meta, mats = read_matrices(args.jacobians)
    built = build_monodromies(meta, mats, args.period)
    if not built:
        print("no tuning cell contains the full attitude orbit")
        return 1

    Phis = [P for _, _, P in built]
    steps = built[0][1]
    dt = built[0][0][0]
    rho = [float(abs(np.linalg.eigvals(P)).max()) for P in Phis]

    # Frozen-point comparison: hold the first orbit pose for the same span.
    scale, scale_inv = np.diag(STATE_SCALE), np.diag(1.0 / STATE_SCALE)
    still = None
    for m, A in zip(meta, mats, strict=True):
        pose = (m["roll"], m["pitch"], m["wx"], m["wy"], m["wz"])
        if pose == (0.0, 0.0, 0.0, 0.0, 0.0) and (m["dt"], m["tau"], m["sigma_aw"], m["sigma_s"]) == built[0][0]:
            still = scale_inv @ A[:N, :N] @ scale
    rho_still = float(abs(np.linalg.eigvals(still)).max()) if still is not None else float("nan")

    print(f"{len(Phis)} monodromies over a {args.period:g} s orbit ({steps} steps at dt={dt:g})")
    print(f"  rho(Phi) in [{min(rho):.4f}, {max(rho):.4f}]")
    print(f"  frozen rho per step (held still) = {rho_still:.9f} "
          f"-> e-fold {dt / (1 - rho_still):.1f} s")

    gamma, L = smallest_common_contraction(Phis, cp)
    if L is None:
        print("no common metric found for the monodromies")
        return 1

    ev = np.linalg.eigvalsh(L)
    c1 = float(np.sqrt(ev[-1] / ev[0]))
    per_step = float(gamma ** (1.0 / steps)) ** 0.5
    efold = dt / (1.0 - per_step)

    print(f"\ncommon metric found (so delta_L = 0; A4 is not needed):")
    print(f"  V decays by gamma = {gamma:.5f} per period; ||e|| by {np.sqrt(gamma):.4f}")
    print(f"  c1 = sqrt(cond L) = {c1:.1f}")
    print(f"  equivalent e-folding = {efold:.2f} s")

    if args.json:
        result: dict[str, Any] = {
            "schema": "ocean-imu-switched-certificate-v1",
            "scope": "switched model: products of frozen-point Jacobians over an attitude orbit",
            "period_s": args.period,
            "steps_per_period": steps,
            "monodromy_count": len(Phis),
            "monodromy_rho_min": min(rho),
            "monodromy_rho_max": max(rho),
            "frozen_rho_per_step_still": rho_still,
            "frozen_efold_seconds_still": dt / (1 - rho_still),
            "common_metric_gamma_per_period": gamma,
            "common_metric_c1": c1,
            "switched_efold_seconds": efold,
            "delta_l_required": 0.0,
        }
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
