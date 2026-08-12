#!/usr/bin/env python3
"""Bias--variance analysis of the OU-III integral (S=0) regularizer.

Section "Riccati interpretation of the S=0 regularizer" of the OU-III paper
derives a pole-placement law for R_S and, in the strong acceleration-observation
branch that the deployed filter demonstrably occupies, predicts a schedule with
*no* leading-order dependence on sigma_aw.  The amplitude-exponent ablation
(tools/ou_rs_law_ablation.py) measures the opposite: removing the sigma_aw
factor costs +0.300 percent of Hs, and the optimum sits at the deployed p = 1.

This driver resolves that contradiction inside the paper's own reduced model.

The covariance recursion behind the pole law treats

    y_S = 0 = S + n_S,   n_S ~ N(0, R_S)

as a real observation.  It is not: the true third integral of a wave record is
a band-limited random process whose variance is set by the sea, not by R_S
(Theorem "Low-frequency-regularized triple integral of OU acceleration" gives
Var[S] = C_OU sigma_aw^2 tau^6).  The Kalman covariance is therefore the
*believed* covariance of a mis-specified filter, not its mean square error.
Restoring the missing term gives, for the drift-band triple integrator with
gain K_d = [2 w_R, 2 w_R^2, w_R^3]^T,

    MSE_p(w_R) = int |T_pe(jw)|^2 q_eff dw/pi          (drift; the Riccati term)
               + int |T_pn(jw)|^2 S_S(w) dw            (distortion; not in the
                                                        Riccati at all)

    T_pe(s) = (s + k1)/D(s),  T_pn(s) = -s(k2 s + k3)/D(s),
    D(s) = s^3 + k1 s^2 + k2 s + k3 = (s + w_R)(s^2 + w_R s + w_R^2),

with S_S = S_eta/w^2 the one-sided spectrum of the true third integral.  Only
the second term carries the sea-state amplitude, and it is the term the
covariance analysis cannot see.

Minimising the sum gives w_R* ~ (q_eff/m_{-4})^{1/7}, where

    m_{-4} = int S_eta/w^4 dw = Var[int S dt] ~ sigma_aw^2 tau^8

is the negative fourth moment of the elevation spectrum -- one integration
below Var[S], because the leakage 2 w_R^2/w is itself an integrator, but
carrying the same sea-state amplitude.  Hence

    r_S* ~ q_eff^{1/14} m_{-4}^{3/7}   and   r_S* ~ sigma_aw^p,
    p = (m + 12)/14   when   q_eff ~ sigma_aw^m.

So p = 6/7 = 0.857 even for a completely amplitude-independent q_eff -- the
paper's own strong-branch reading -- and p = 1 exactly when q_eff ~ sigma_aw^2.
Reproducing the measured p = 0 would need q_eff ~ sigma_aw^-12.  The
experimental result is not in conflict with the Riccati algebra; it is in
conflict with the assumption that a covariance recursion built on a false
measurement bounds the estimator's actual error.

Typical use:

    python3 tools/rs_law_bias_variance.py
    python3 tools/rs_law_bias_variance.py --exponents 0,1,2,3
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

GRAV = 9.80665

# Calibrated OU-III operating points, from
# doc/kalman_ou_iii/w3d-ou-validation-tuning-points-generated.tex-part, paired
# with the peak wavelength encoded in each JONSWAP scenario name.
# (label, Hs [m], peak wavelength [m], tau [s], sigma_aw [m/s^2])
SEAS = (
    ("Hs=0.27", 0.270, 14.047, 1.247, 0.354),
    ("Hs=1.50", 1.500, 50.710, 2.179, 0.724),
    ("Hs=4.00", 4.000, 112.766, 3.589, 1.124),
    ("Hs=8.50", 8.500, 202.839, 4.222, 1.420),
)

# Deployed regularizer schedule (Sec. "What the wrapper currently deploys").
C_R = 0.35
C_T = 0.015 / 1.1
TS0 = 0.015
TS_MIN, TS_MAX = 0.005, 0.250

# Validated accelerometer noise floor: R_a dt at dt = 5 ms.
DT_IMU = 0.005
SIGMA_NA = 0.0148
R_A = SIGMA_NA**2 * DT_IMU

NOMINAL = 1  # index of the Hs = 1.5 m anchor sea in SEAS


def peak_period(wavelength_m: float) -> float:
    """Deep-water peak period of a scenario's peak wavelength."""
    return math.sqrt(2.0 * math.pi * wavelength_m / GRAV)


def jonswap(w: np.ndarray, hs: float, tp: float, gamma: float = 3.3) -> np.ndarray:
    """One-sided JONSWAP elevation spectrum scaled so int S dw = (Hs/4)^2."""
    wp = 2.0 * math.pi / tp
    sigma = np.where(w <= wp, 0.07, 0.09)
    peak = np.exp(-((w - wp) ** 2) / (2.0 * sigma**2 * wp**2))
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        s = w**-5.0 * np.exp(-1.25 * (wp / w) ** 4) * gamma**peak
    s = np.nan_to_num(np.where(w > 0.0, s, 0.0))
    return s * ((hs / 4.0) ** 2 / np.trapezoid(s, w))


def transfer_gains(w: np.ndarray, w_r: float) -> tuple[np.ndarray, np.ndarray]:
    """Squared magnitudes of the two error transfers of the drift-band observer.

    Returns (|T_pe|^2, |T_pn|^2): residual acceleration error -> position error,
    and true third integral -> position error.
    """
    s = 1j * w
    k1, k2, k3 = 2.0 * w_r, 2.0 * w_r**2, w_r**3
    den = s**3 + k1 * s**2 + k2 * s + k3
    return (np.abs((s + k1) / den) ** 2,
            np.abs(-s * (k2 * s + k3) / den) ** 2)


def mse_terms(w_r: float, w: np.ndarray, s_eta: np.ndarray,
              q_eff: np.ndarray) -> tuple[float, float]:
    """Drift and distortion contributions to the mean square position error."""
    t_pe2, t_pn2 = transfer_gains(w, w_r)
    drift = float(np.trapezoid(t_pe2 * q_eff / math.pi, w))
    distortion = float(np.trapezoid(t_pn2 * s_eta / w**2, w))
    return drift, distortion


def optimal_pole(w: np.ndarray, s_eta: np.ndarray,
                 q_eff: np.ndarray) -> tuple[float, float]:
    """Minimise the total mean square position error over w_R."""
    lo, hi = -3.0, 1.0
    grid = np.logspace(lo, hi, 161)
    best_w, best_v = grid[0], math.inf
    for _ in range(4):
        vals = np.array([sum(mse_terms(x, w, s_eta, q_eff)) for x in grid])
        i = int(np.argmin(vals))
        best_w, best_v = float(grid[i]), float(vals[i])
        lo = math.log10(grid[max(i - 1, 0)])
        hi = math.log10(grid[min(i + 1, len(grid) - 1)])
        grid = np.logspace(lo, hi, 41)
    return best_w, best_v


def pseudo_period(tau: float) -> float:
    """Deployed self-similar pseudo-update cadence with its clamps."""
    return min(max(C_T * tau, TS_MIN), TS_MAX)


def deployed_rs(sigma_a: float, tau: float) -> float:
    """r_S actually reaching the filter, after the information renormalization."""
    return C_R * sigma_a * tau**3 * math.sqrt(TS0 / pseudo_period(tau))


def rs_from_pole(w_r: float, q_eff_at_pole: float, tau: float) -> float:
    """Invert w_R = (q_eff/(R_S T_S))^{1/6} for r_S = sqrt(R_S)."""
    ts = pseudo_period(tau)
    return math.sqrt(q_eff_at_pole / w_r**6 / ts)


def frequency_grid() -> np.ndarray:
    return np.logspace(-4.0, 1.2, 6001)


def self_check(w: np.ndarray) -> None:
    """Reproduce the nominal Riccati covariance from the same transfers.

    Feeding the *assumed* white pseudo-measurement noise (two-sided density
    r_S,c) back through T_pn must return the CARE solution P_pp = 3 q/w_R^3,
    which the closed form of the paper implies.  The split is exactly half and
    half: the Riccati attributes 1.5 q/w_R^3 of the steady-state position
    variance to a pseudo-measurement noise that does not exist.
    """
    print("self-check: nominal Riccati P_pp against the closed form")
    print(f"  {'w_R':>8}{'process':>12}{'pseudo':>12}{'sum':>12}"
          f"{'3q/w_R^3':>12}{'ratio':>9}")
    for w_r in (0.05, 0.2, 1.0):
        q = 1.0
        r_sc = q / w_r**6
        t_pe2, t_pn2 = transfer_gains(w, w_r)
        proc = float(np.trapezoid(t_pe2 * q / math.pi, w))
        pseudo = float(np.trapezoid(t_pn2 * r_sc / math.pi, w))
        exact = 3.0 * q / w_r**3
        print(f"  {w_r:8.2f}{proc:12.4e}{pseudo:12.4e}{proc + pseudo:12.4e}"
              f"{exact:12.4e}{(proc + pseudo) / exact:9.5f}")
    print()


def sea_table(w: np.ndarray) -> None:
    """MSE-optimal regularizer per calibrated sea, against the deployed one."""
    print("MSE-optimal regularizer, q_eff = 2 r_a (the paper's strong branch)")
    print(f"  r_a = {R_A:.3e} m^2/s^3, q_eff = {2 * R_A:.3e} m^2/s^3\n")
    print(f"  {'sea':>8}{'Tp':>7}{'tau':>7}{'sigma_aw':>10}{'w_R*':>9}"
          f"{'w_R* tau':>10}{'r_S*':>10}{'r_S dep':>10}{'r_S*/dep':>10}"
          f"{'rms %Hs':>9}")
    for label, hs, lam, tau, sigma_a in SEAS:
        tp = peak_period(lam)
        s_eta = jonswap(w, hs, tp)
        q_eff = np.full_like(w, 2.0 * R_A)
        w_r, mse = optimal_pole(w, s_eta, q_eff)
        rs = rs_from_pole(w_r, 2.0 * R_A, tau)
        dep = deployed_rs(sigma_a, tau)
        print(f"  {label:>8}{tp:7.2f}{tau:7.2f}{sigma_a:10.3f}{w_r:9.4f}"
              f"{w_r * tau:10.3f}{rs:10.3f}{dep:10.3f}{rs / dep:10.3f}"
              f"{100.0 * math.sqrt(mse) / hs:9.2f}")
    print()


def deployed_balance(w: np.ndarray) -> None:
    """Where the deployed schedule sits on the bias--variance curve."""
    print("Deployed schedule: term balance at its own operating point")
    print(f"  {'sea':>8}{'w_R dep':>10}{'drift':>12}{'distortion':>12}"
          f"{'distort/tot':>13}{'rms %Hs':>9}{'w_R*/w_R dep':>14}")
    for label, hs, lam, tau, sigma_a in SEAS:
        s_eta = jonswap(w, hs, peak_period(lam))
        q_eff = np.full_like(w, 2.0 * R_A)
        r_sc = deployed_rs(sigma_a, tau) ** 2 * pseudo_period(tau)
        w_r = (2.0 * R_A / r_sc) ** (1.0 / 6.0)
        drift, distortion = mse_terms(w_r, w, s_eta, q_eff)
        w_star, _ = optimal_pole(w, s_eta, q_eff)
        total = drift + distortion
        print(f"  {label:>8}{w_r:10.4f}{drift:12.3e}{distortion:12.3e}"
              f"{distortion / total:13.3f}{100.0 * math.sqrt(total) / hs:9.2f}"
              f"{w_star / w_r:14.3f}")
    print()


def amplitude_sweep(w: np.ndarray, exponents: list[float], span: float,
                    points: int) -> None:
    """Controlled amplitude sweep at fixed Tp and tau: measure p directly."""
    _, hs0, lam0, tau0, sigma0 = SEAS[NOMINAL]
    tp = peak_period(lam0)
    scales = np.logspace(-span, span, points)
    print(f"Controlled amplitude sweep about {SEAS[NOMINAL][0]} "
          f"(Tp = {tp:.2f} s and tau = {tau0:.2f} s held fixed)")
    print(f"  {'q_eff ~ sigma_aw^m':>20}{'p measured':>13}"
          f"{'(m+12)/14':>12}{'w_R*/w_p span':>16}")
    for m in exponents:
        logs, logr, poles = [], [], []
        for scale in scales:
            sigma_a = sigma0 * scale
            s_eta = jonswap(w, hs0 * scale, tp)
            q_at = 2.0 * R_A * (sigma_a / sigma0) ** m
            w_r, _ = optimal_pole(w, s_eta, np.full_like(w, q_at))
            logs.append(math.log(sigma_a))
            logr.append(math.log(rs_from_pole(w_r, q_at, tau0)))
            poles.append(w_r / (2.0 * math.pi / tp))
        p = float(np.polyfit(logs, logr, 1)[0])
        print(f"  {m:20.2f}{p:13.4f}{(m + 12.0) / 14.0:12.4f}"
              f"{min(poles):8.3f}-{max(poles):<7.3f}")
    print("\n  p = 0, the schedule the covariance analysis predicts, would need"
          "\n  q_eff ~ sigma_aw^-12.  The ablation minimum at p = 1 implies"
          "\n  q_eff ~ sigma_aw^2 in the drift band.\n")


def time_scale_sweep(w: np.ndarray, span: float, points: int) -> None:
    """Controlled time-scale sweep at fixed sigma_aw: measure the tau exponent.

    Tp (and with it tau) is scaled while Hs is renormalised to hold
    sigma_aw^2 = m_4 fixed, so only the time scale moves.  Exact OU similarity
    would give m_{-4} ~ tau^8 and r_S* ~ tau^(41/14); the JONSWAP shape is not
    exactly self-similar, so the measured moment exponent is used to check the
    identity d log r_S*/d log tau = (6/7 * d log m_{-4}/d log tau - 1)/2.
    """
    _, hs0, lam0, tau0, sigma0 = SEAS[NOMINAL]
    tp0 = peak_period(lam0)
    logt, logr, logm = [], [], []
    for scale in np.logspace(-span, span, points):
        tp, tau = tp0 * scale, tau0 * scale
        s_eta = jonswap(w, hs0, tp)
        s_eta = s_eta * (sigma0**2 / float(np.trapezoid(w**4 * s_eta, w)))
        q_eff = np.full_like(w, 2.0 * R_A)
        w_r, _ = optimal_pole(w, s_eta, q_eff)
        logt.append(math.log(tau))
        logr.append(math.log(rs_from_pole(w_r, 2.0 * R_A, tau)))
        logm.append(math.log(float(np.trapezoid(s_eta / w**4, w))))
    slope_r = float(np.polyfit(logt, logr, 1)[0])
    slope_m = float(np.polyfit(logt, logm, 1)[0])
    print("Controlled time-scale sweep (sigma_aw held fixed, amplitude-free q_eff)")
    print(f"  d log m_-4 / d log tau = {slope_m:8.4f}   (exact OU similarity: 8)")
    print(f"  d log r_S* / d log tau = {slope_r:8.4f}   "
          f"(identity: {(6.0 / 7.0 * slope_m - 1.0) / 2.0:.4f}, "
          f"exact similarity 41/14 = {41.0 / 14.0:.4f})")
    print(f"  deployed exponent      = {2.5:8.4f}\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exponents", default="0,1,2,3",
                    help="amplitude exponents m of q_eff to sweep (default 0,1,2,3)")
    ap.add_argument("--sweep-span", type=float, default=1.0,
                    help="half-width in decades of the amplitude sweep (default 1)")
    ap.add_argument("--sweep-points", type=int, default=21)
    ap.add_argument("--skip-self-check", action="store_true")
    args = ap.parse_args(argv)

    exponents = [float(e) for e in args.exponents.split(",") if e.strip()]
    w = frequency_grid()

    if not args.skip_self_check:
        self_check(w)
    sea_table(w)
    deployed_balance(w)
    amplitude_sweep(w, exponents, args.sweep_span, args.sweep_points)
    time_scale_sweep(w, 0.3, 13)
    return 0


if __name__ == "__main__":
    sys.exit(main())
