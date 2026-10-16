#!/usr/bin/env python3
"""Regenerate the deterministic example chart tables used by carma-heave.tex.

The time traces are exact-discretization samples of stationary, linear,
Gaussian SDEs; the spectral curves are theoretical continuous-time PSDs.
The figures are illustrative, not marine data or filter error measurements.

Requires NumPy and SciPy. Run from any directory:
    python3 doc/kalman_ou_iii/carma-heave-generate.py
"""
from pathlib import Path
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov

OUTPUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(20261016)
DT = 0.04
DURATION = 480.0
T = np.arange(0, DURATION, DT)
W0 = 2 * np.pi * 0.29
ZETA = 0.22

models = {}
models["ou"] = (
    np.array([[-2 * np.pi * 0.12]]),
    np.array([[1.0]]),
    np.array([[1.0]]),
)
models["sho"] = (
    np.array([[0, 1], [-W0**2, -2 * ZETA * W0]]),
    np.array([[0.0], [1.0]]),
    np.array([[1.0, 0.0]]),
)
models["colored"] = (
    np.array([
        [0, 1, 0],
        [-W0**2, -2 * ZETA * W0, 1],
        [0, 0, -2 * np.pi * 0.40],
    ]),
    np.array([[0.0], [0.0], [1.0]]),
    np.array([[1.0, 0.0, 0.0]]),
)

w1, z1 = 2 * np.pi * 0.27, 0.20
w2, z2 = 2 * np.pi * 0.56, 0.60
a = np.polymul([1, 2 * z1 * w1, w1**2],
               [1, 2 * z2 * w2, w2**2])
A4 = np.zeros((4, 4))
A4[0, 1] = A4[1, 2] = A4[2, 3] = 1.0
A4[3, :] = -a[1:][::-1]
models["carma42"] = (
    A4,
    np.array([[0.0], [0.0], [0.0], [1.0]]),
    np.array([[0.0, 0.0, 1.0, 0.0]]),
)

freqs = np.linspace(0, 1.3, 1200)
traces, spectra = {}, {}
for name, (A, G, C) in models.items():
    stationary_P = solve_continuous_lyapunov(A, -(G @ G.T))
    stationary_P = 0.5 * (stationary_P + stationary_P.T)
    sd = np.sqrt((C @ stationary_P @ C.T)[0, 0])
    F = expm(A * DT)
    Q = stationary_P - F @ stationary_P @ F.T
    Q = 0.5 * (Q + Q.T)
    chol = np.linalg.cholesky(Q + 1e-14 * np.eye(len(A)))
    innovations = RNG.standard_normal((len(T) - 1, len(A))) @ chol.T
    states = np.zeros((len(T), len(A)))
    states[0] = RNG.multivariate_normal(np.zeros(len(A)), stationary_P)
    for k in range(1, len(T)):
        states[k] = F @ states[k - 1] + innovations[k - 1]
    traces[name] = (states @ C.T).ravel() / sd

    transfer = np.array([
        (C @ np.linalg.solve(1j * (2 * np.pi * f) * np.eye(len(A)) - A,
                             G))[0, 0] for f in freqs
    ])
    power = abs(transfer)**2 / sd**2
    spectra[name] = 10 * np.log10(np.maximum(power / power.max(), 1e-12))

columns = ["ou", "sho", "colored", "carma42"]
with (OUTPUT / "carma-heave-time.csv").open("w") as file:
    file.write("t," + ",".join(columns) + "\n")
    offsets = dict(ou=9.0, sho=6.0, colored=3.0, carma42=0.0)
    for i in range(0, len(T), 5):
        if T[i] > 65 + 1e-7:
            break
        values = [traces[c][i] + offsets[c] for c in columns]
        file.write(",".join([f"{T[i]:.2f}"] +
                            [f"{v:.5f}" for v in values]) + "\n")

with (OUTPUT / "carma-heave-spectra.csv").open("w") as file:
    file.write("f," + ",".join(columns) + "\n")
    for j in range(131):
        f = j / 100.0
        i = np.abs(freqs - f).argmin()
        values = [np.clip(spectra[c][i], -70, 4) for c in columns]
        file.write(",".join([f"{f:.2f}"] +
                            [f"{v:.5f}" for v in values]) + "\n")

if __name__ == "__main__":
    print("Wrote reproducible figure tables in", OUTPUT)