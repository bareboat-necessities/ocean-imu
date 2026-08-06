#!/usr/bin/env python3
"""Search for a common quadratic Lyapunov function over the OU-III envelope.

Assumption A4 of doc/kalman_ou_iii/w3d-iss-stability.tex-part requires the
Lyapunov metric to drift by less than 1/(l_bar - 1) per step. The certificate
measures that drift and finds it violated by more than ten orders of
magnitude, driven mostly by roll and pitch moving at wave frequency.

The standard way out is to stop using a per-point metric. If a single L > 0
satisfies

    F(xi, theta)' L F(xi, theta) - L <= -eps I

for every operating point, then delta_L = 0 identically and A4 can be dropped
from the theorem. This script tests whether such an L exists.

It probes incrementally: order the envelope by spectral radius and grow the
active set one matrix at a time, reporting the first size at which the LMI
becomes infeasible. That is faster than constraint generation and it answers
the question directly -- if no common metric exists for the three worst
points, none exists for the envelope.

Result on the shipped envelope: infeasible from three matrices. The
problem is also numerically degenerate, because the near-marginal drift mode
leaves 1 - rho^2 of order 1e-6, which is at solver tolerance; even the
single-matrix solve satisfies its own constraint only to about 1e-9. Both
facts are reported so the negative result can be read with its caveat.

This is a diagnostic, not part of the CI certificate, and it needs cvxpy.

    python3 tools/iss_certificate/search_common_lyapunov.py iss_jacobians.txt
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Physical scale of each 3-block, so the LMI is not conditioned by unit choice:
# dtheta [rad], bg [rad/s], v [m/s], p [m], S [m s], aw [m/s^2].
STATE_SCALE = np.repeat([1e-3, 1e-3, 1.0, 1.0, 30.0, 1.0], 3)
N = 18


def read_performance_matrices(path: Path) -> list[np.ndarray]:
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().split()
        if len(header) != 3 or header[0] != "OCEAN_IMU_ISS_JACOBIANS_V1":
            raise ValueError("invalid Jacobian file header")
        n, count = int(header[1]), int(header[2])
        out = []
        for _ in range(count):
            fh.readline()
            rows = [[float(x) for x in fh.readline().split()] for _ in range(n)]
            out.append(np.asarray(rows, dtype=float)[:N, :N])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("jacobians", type=Path)
    ap.add_argument("--max-subset", type=int, default=12)
    args = ap.parse_args()

    try:
        import cvxpy as cp
    except ImportError:
        print("cvxpy is required for this diagnostic: pip install cvxpy")
        return 2

    scale = np.diag(STATE_SCALE)
    scale_inv = np.diag(1.0 / STATE_SCALE)
    G = [scale_inv @ A @ scale for A in read_performance_matrices(args.jacobians)]

    rho = np.asarray([abs(np.linalg.eigvals(A)).max() for A in G])
    order = list(np.argsort(-rho))
    print(f"{len(G)} matrices, rho in [{rho.min():.9f}, {rho.max():.9f}]")
    print(
        f"worst 1 - rho^2 = {1.0 - rho.max() ** 2:.3e}  "
        "<- the LMI's numerical floor; a verdict at or below this is not conclusive"
    )

    feasible_upto = 0
    for k in range(1, args.max_subset + 1):
        subset = order[:k]
        L = cp.Variable((N, N), PSD=True)
        constraints = [L >> np.eye(N), L << 1e12 * np.eye(N)]
        constraints += [G[i].T @ L @ G[i] - L << 0 for i in subset]
        problem = cp.Problem(cp.Minimize(cp.trace(L)), constraints)
        try:
            problem.solve(solver=cp.CLARABEL)
        except cp.error.SolverError as exc:
            print(f"k={k:2d} solver error: {exc}")
            break

        if L.value is None:
            print(f"k={k:2d} {problem.status}")
            print(
                f"\nINFEASIBLE from {k} matrices: no common quadratic Lyapunov "
                "function exists over this envelope.\n"
                "A4 cannot be removed from Theorem 1 by this route while the\n"
                "near-marginal drift mode is present."
            )
            return 1

        Lv = 0.5 * (L.value + L.value.T)
        ev = np.linalg.eigvalsh(Lv)
        worst = max(float(np.linalg.eigvalsh(G[i].T @ Lv @ G[i] - Lv)[-1]) for i in subset)
        flag = "" if worst < 0 else "   (violated: at solver tolerance)"
        print(f"k={k:2d} {problem.status:10s} cond={ev[-1] / ev[0]:.3e} worst={worst:+.3e}{flag}")
        if worst < 0:
            feasible_upto = k

    print(
        f"\nA strictly satisfied common metric was found for at most "
        f"{feasible_upto} of the worst-case matrices."
    )
    return 0 if feasible_upto == args.max_subset else 1


if __name__ == "__main__":
    sys.exit(main())
