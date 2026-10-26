#!/usr/bin/env python3
"""Fit CARMA(3,0), CARMA(4,1) and CARMA(5,2) displacement spectra to JONSWAP.

Objective
---------
For a model displacement PSD S_C(f; theta) = |B(i w)|^2 / |P(i w)|^2 and the
target JONSWAP elevation PSD S_J(f) (f_p = 0.29 Hz, gamma = 3.3), minimize the
RMS log-spectral error

    J(theta) = sqrt( mean_f [ d(f) - mean_f d ]^2 ),
    d(f)     = 10 log10 S_C(f; theta) - 10 log10 S_J(f),

on a uniform 0.001 Hz grid over 0.17 <= f <= 1.30 Hz.  Subtracting mean_f d
makes the fit amplitude free: only spectral *shape* is fitted, and the driving
intensity q is set at run time from the measured acceleration variance.

Constraints
-----------
* every pole damping ratio in [0.05, 1], every pole frequency positive, so all
  poles lie strictly in the left half plane;
* CARMA(5,2) numerator zeta_z >= 0.20 (the "constrained" fit).  Without it the
  optimizer places a narrow notch near the lower fitting boundary; that
  unconstrained fit is reported for reference only;
* CARMA(4,1) numerator b(s) = s + 2 pi f_z0 with f_z0 >= 0.

The relative degree p - q = 3 of every model gives a finite-variance ideal
acceleration (displacement tail ~ f^-6).  JONSWAP has an f^-5 tail, so no
model here is asymptotically equivalent to it.

Bounded trust-region least squares from a fixed seeded set of starts makes
the result deterministic.
Run:  python3 tools/carma_fit.py [--out reports/results/carma_model/carma_fit.json]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

TWO_PI = 2.0 * np.pi
FP_HZ = 0.29
GAMMA = 3.3
BAND_HZ = (0.17, 1.30)
GRID = np.round(np.arange(BAND_HZ[0], BAND_HZ[1] + 5e-4, 1e-3), 6)

# Preliminary constrained CARMA(5,2) fit supplied with the study request.
PRELIMINARY_52 = dict(f1=0.27812, z1=0.09715, f2=0.44274, z2=0.52299,
                      fl=1.86402, fz=0.15105, zz=0.20)


def jonswap(f, fp=FP_HZ, gamma=GAMMA):
    """JONSWAP shape (unnormalized); gamma = 1 gives Pierson-Moskowitz."""
    f = np.asarray(f, dtype=float)
    out = np.zeros_like(f)
    pos = f > 0
    fpos = f[pos]
    sig = np.where(fpos <= fp, 0.07, 0.09)
    r = np.exp(-((fpos - fp) ** 2) / (2.0 * sig ** 2 * fp ** 2))
    out[pos] = fpos ** -5 * np.exp(-1.25 * (fp / fpos) ** 4) * gamma ** r
    return out


def poly_quad(s, f_hz, zeta):
    w = TWO_PI * f_hz
    return s * s + 2.0 * zeta * w * s + w * w


def carma_psd(f, order, p):
    """|B(i w)|^2/|P(i w)|^2 for monic B and P (unit driving intensity)."""
    s = 1j * TWO_PI * np.asarray(f, dtype=float)
    if order == "30":
        B = np.ones_like(s)
        P = poly_quad(s, p["f1"], p["z1"]) * (s + TWO_PI * p["fl"])
    elif order == "41":
        B = s + TWO_PI * p["fz0"]
        P = poly_quad(s, p["f1"], p["z1"]) * poly_quad(s, p["f2"], p["z2"])
    elif order == "52":
        B = poly_quad(s, p["fz"], p["zz"])
        P = (poly_quad(s, p["f1"], p["z1"]) * poly_quad(s, p["f2"], p["z2"])
             * (s + TWO_PI * p["fl"]))
    else:
        raise ValueError(order)
    return np.abs(B) ** 2 / np.abs(P) ** 2


NAMES = {
    "30": ["f1", "z1", "fl"],
    "41": ["f1", "z1", "f2", "z2", "fz0"],
    "52": ["f1", "z1", "f2", "z2", "fl", "fz", "zz"],
}
BOUNDS = {
    "f1": (0.15, 0.60), "z1": (0.05, 1.0),
    "f2": (0.15, 3.00), "z2": (0.05, 1.0),
    "fl": (0.05, 5.00),
    "fz0": (0.0, 2.0),
    "fz": (0.05, 1.00), "zz": (0.20, 1.0),
}


def log_error(order, p, grid=GRID):
    d = 10.0 * np.log10(carma_psd(grid, order, p)) - 10.0 * np.log10(jonswap(grid))
    return d - d.mean()


def rms_db(order, p):
    return float(np.sqrt(np.mean(log_error(order, p) ** 2)))


def fit(order, zz_min=0.20, n_starts=48):
    """Bounded least squares on the shape residual from a fixed set of starts."""
    from scipy.optimize import least_squares

    names = NAMES[order]
    bounds = [BOUNDS[n] if n != "zz" else (zz_min, 1.0) for n in names]
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])

    def resid(x):
        return log_error(order, dict(zip(names, x)))

    rng = np.random.default_rng(20261025)
    starts = [lo + (hi - lo) * (0.1 + 0.8 * rng.random(len(names))) for _ in range(n_starts)]
    if order == "52":
        starts.append(np.clip([PRELIMINARY_52[n] for n in names], lo, hi))
    best = None
    for x0 in starts:
        res = least_squares(resid, x0, bounds=(lo, hi), method="trf",
                            x_scale="jac", ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=4000)
        if best is None or res.cost < best.cost:
            best = res
    return dict(zip(names, best.x))


def describe(order, p):
    """Diagnostics of one parameter set (unit driving intensity)."""
    f = np.linspace(1e-5, 20.0, 400001)
    df = f[1] - f[0]
    S = carma_psd(f, order, p)
    J = jonswap(f)
    w = TWO_PI * f
    m0, m2 = np.sum(S) * df, np.sum(w ** 2 * S) * df
    m0j, m2j = np.sum(J) * df, np.sum(w ** 2 * J) * df
    peak = float(f[S.argmax()])
    s0 = carma_psd(np.array([0.0]), order, p)[0]
    if order == "30":
        poles = [complex(r) for r in np.roots(np.polymul(
            [1, 2 * p["z1"] * TWO_PI * p["f1"], (TWO_PI * p["f1"]) ** 2],
            [1, TWO_PI * p["fl"]]))]
        zeros = []
    elif order == "41":
        poles = [complex(r) for r in np.roots(np.polymul(
            [1, 2 * p["z1"] * TWO_PI * p["f1"], (TWO_PI * p["f1"]) ** 2],
            [1, 2 * p["z2"] * TWO_PI * p["f2"], (TWO_PI * p["f2"]) ** 2]))]
        zeros = [complex(-TWO_PI * p["fz0"])]
    else:
        poles = [complex(r) for r in np.roots(np.polymul(np.polymul(
            [1, 2 * p["z1"] * TWO_PI * p["f1"], (TWO_PI * p["f1"]) ** 2],
            [1, 2 * p["z2"] * TWO_PI * p["f2"], (TWO_PI * p["f2"]) ** 2]),
            [1, TWO_PI * p["fl"]]))]
        zeros = [complex(r) for r in np.roots(
            [1, 2 * p["zz"] * TWO_PI * p["fz"], (TWO_PI * p["fz"]) ** 2])]
    return {
        "params": {k: float(v) for k, v in p.items()},
        "rms_db_0p17_1p30": rms_db(order, p),
        "max_abs_db_0p17_1p30": float(np.max(np.abs(log_error(order, p)))),
        "peak_hz": peak,
        "tz_times_fp_model": float(TWO_PI * np.sqrt(m0 / m2) * FP_HZ),
        "tz_times_fp_jonswap": float(TWO_PI * np.sqrt(m0j / m2j) * FP_HZ),
        "dc_level_db_re_peak": float(10 * np.log10(s0 / S.max())) if s0 > 0 else None,
        "variance_fraction_below_0p1hz": float(np.sum(S[f < 0.1]) / np.sum(S)),
        "max_pole_real_part": float(max(r.real for r in poles)),
        "poles": [[float(r.real), float(r.imag)] for r in poles],
        "zeros": [[float(r.real), float(r.imag)] for r in zeros],
        "tail_exponent_f": -6,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parents[1]
                    / "reports" / "results" / "carma_model" / "carma_fit.json")
    args = ap.parse_args(argv)

    result = {
        "target": {"spectrum": "JONSWAP", "fp_hz": FP_HZ, "gamma": GAMMA,
                   "band_hz": list(BAND_HZ), "grid_step_hz": 1e-3},
        "objective": "RMS of 10log10(S_C/S_J) minus its band mean (shape only)",
        "models": {},
    }
    result["models"]["carma30"] = describe("30", fit("30"))
    result["models"]["carma41"] = describe("41", fit("41"))
    result["models"]["carma52_preliminary"] = describe("52", PRELIMINARY_52)
    result["models"]["carma52_refit_constrained"] = describe("52", fit("52", zz_min=0.20))
    result["models"]["carma52_refit_unconstrained"] = describe("52", fit("52", zz_min=0.01))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for name, m in result["models"].items():
        print(f"{name:30s} rms={m['rms_db_0p17_1p30']:.3f} dB "
              f"peak={m['peak_hz']:.4f} Hz params=" +
              ", ".join(f"{k}={v:.5f}" for k, v in m["params"].items()))


if __name__ == "__main__":
    main()
