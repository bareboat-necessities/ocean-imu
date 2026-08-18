#!/usr/bin/env python3
"""Analytical coefficients of the OU-II dual-channel physical-MSE pseudo law.

doc/kalman_ou_ii/ou2-dual-regularization-mse.tex derives the joint
displacement-MSE optimum of the two zero pseudo-measurements (p = 0 and v = 0)
that regularize the OU-II integration chain.  Its result is a pair of *shapes*,

    r_p = C_P q_eff^(1/10) sigma_a,B^(4/5) tau^(12/5) / sqrt(T_S),
    r_v = C_V q_eff^(1/10) sigma_a,B^(4/5) tau^(7/5)  / sqrt(T_S),

against the empirical sigma_aw tau^2 and sigma_aw tau base laws.  This driver
supplies the two constants the shapes leave open, by evaluating the derivation's
own stationarity conditions on the eight reference spectra at the operating
points the OU-II tuner actually selects on them.  Nothing here is fitted to the
empirical schedule.

C_P and C_V absorb dimensionless spectral moments of the self-similar
displacement spectrum,

    M_0   = int S_eta dw          ~ sigma_a^2 tau^4,
    M_-2  = int S_eta / w^2 dw    ~ sigma_a^2 tau^6,

through Eqs. (wp-opt) and (chi-opt),

    omega_p*^5 = 3 q / (8 sqrt2 M_-2),
    chi*       = 9 q / (16 sqrt2 M_0 omega_p*^3),

with rho_p = q/omega_p^4, rho_v = q/(omega_p^2 chi) and r = sqrt(rho/T_S).

Two quantities are reported per record and in aggregate:

  C_P              the position-channel coefficient, and C_V for velocity.
  C_P/C_V          the channel ratio, which Corollary (optimal pseudo-channel
                   ratio) identifies as the better-determined half of the
                   prediction because q_eff and the cadence both cancel:
                       (r_p/r_v)^2 = (3/2) M_-2 / M_0.

The weak-regularization expansion behind Eqs. (wp-opt)-(chi-opt) replaces the
exact wave-distortion integral by M_0 chi^2 + 2 M_-2 omega_p^2.  --exact
re-solves the same objective with the exact |G(jw)-1|^2 under the finite-band
spectrum and no chi << 1 expansion, which measures how much that step costs.

Typical use:

    python3 tools/ou2_dual_mse_coefficients.py
    python3 tools/ou2_dual_mse_coefficients.py --exact
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

GRAV = 9.80665

# Density of the residual acceleration error the integration chain sees:
# r_a = R_a h and q_eff = 2 r_a.  sigma_na is the pre-band vertical acceleration
# noise floor the OU-II tuner subtracts as non-wave energy, not the
# accelerometer's bench spec: the note's q is the error left *after*
# acceleration estimation, so it carries attitude and gravity leakage and
# residual bias as well as the sensor.  Mirrors
# R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT of SeaStateFusionFilter_OU_II.h.
DT_IMU = 1.0 / 200.0
SIGMA_NA = 0.12
R_A = SIGMA_NA**2 * DT_IMU
Q_EFF = 2.0 * R_A

# Deployed pseudo-update cadence and its safety clamps.
C_T = 0.015 / 1.1
TS0 = 0.015
TS_MIN, TS_MAX = 1.0 / 200.0, 0.250

# sigma_aw = c_sigma sigma_a,B is the OU prior; the distortion penalty depends
# on the physical band RMS, so the law divides c_sigma back out.
C_SIGMA_OU_II = 0.85

# Deployed empirical coefficients, quoted for comparison only.
C_P_EMPIRICAL = 0.65
C_V_EMPIRICAL = 1.30

# The simulator's first-order spectra are band-limited to 0.02-0.8 Hz.
BAND_HZ = (0.02, 0.8)

# The eight scored records, with the (tau, sigma_aw) the OU-II tuner settles on
# for each.  Operating points are the deployed-arm means over the five seeds of
# reports/results/ou_qmekf_variances/ou_ii_raw.csv; peak wavelengths are encoded
# in the record names.
# (record, family, Hs [m], peak wavelength [m], tau [s], sigma_aw [m/s^2])
SEAS = (
    ("jonswap_H0.270", "jonswap", 0.270, 14.047, 1.2720, 0.3330),
    ("jonswap_H1.500", "jonswap", 1.500, 50.710, 2.1815, 0.6805),
    ("jonswap_H4.000", "jonswap", 4.000, 112.766, 3.5904, 1.0614),
    ("jonswap_H8.500", "jonswap", 8.500, 202.839, 4.2221, 1.3402),
    ("pmstokes_H0.270", "pm", 0.270, 14.047, 1.1643, 0.3722),
    ("pmstokes_H1.500", "pm", 1.500, 50.710, 2.0437, 0.7469),
    ("pmstokes_H4.000", "pm", 4.000, 112.766, 3.2905, 1.0617),
    ("pmstokes_H8.500", "pm", 8.500, 202.839, 4.1025, 1.4109),
)


def peak_period(wavelength_m: float) -> float:
    """Deep-water peak period of a scenario's peak wavelength."""
    return math.sqrt(2.0 * math.pi * wavelength_m / GRAV)


def jonswap(w: np.ndarray, hs: float, tp: float, gamma: float = 3.3) -> np.ndarray:
    """One-sided elevation spectrum scaled so int S dw = (Hs/4)^2.

    gamma = 1 is Pierson-Moskowitz, the first-order spectrum underlying the
    PM-Stokes records.
    """
    wp = 2.0 * math.pi / tp
    sigma = np.where(w <= wp, 0.07, 0.09)
    peak = np.exp(-((w - wp) ** 2) / (2.0 * sigma**2 * wp**2))
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        s = w**-5.0 * np.exp(-1.25 * (wp / w) ** 4) * gamma**peak
    s = np.nan_to_num(np.where(w > 0.0, s, 0.0))
    return s * ((hs / 4.0) ** 2 / np.trapezoid(s, w))


def pseudo_period(tau: float) -> float:
    """Deployed self-similar pseudo-update cadence with its clamps."""
    return min(max(C_T * tau, TS_MIN), TS_MAX)


def frequency_grid(points: int = 4001) -> np.ndarray:
    return np.logspace(math.log10(2.0 * math.pi * BAND_HZ[0]),
                       math.log10(2.0 * math.pi * BAND_HZ[1]), points)


def moments(s_eta: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """(M_0, M_-2): displacement variance and its second negative moment.

    For heave the physical displacement spectrum is the elevation spectrum, and
    the two-sided normalization of Eqs. (M0)-(Mminus2) reduces to the plain
    one-sided integrals.
    """
    return float(np.trapezoid(s_eta, w)), float(np.trapezoid(s_eta / w**2, w))


def asymptotic_optimum(m0: float, m_2: float) -> tuple[float, float]:
    """Eqs. (wp-opt) and (chi-opt): the weak-regularization stationary point."""
    wp = (3.0 * Q_EFF / (8.0 * math.sqrt(2.0) * m_2)) ** 0.2
    chi = 9.0 * Q_EFF / (16.0 * math.sqrt(2.0) * m0 * wp**3)
    return wp, chi


def j_noise(wp: float, chi: float) -> float:
    """Eq. (Jnoise-exact): the exact H_2 norm of the reduced closed loop."""
    return Q_EFF / (2.0 * wp**3 * (1.0 + chi) ** 2 * math.sqrt(2.0 + chi))


def j_wave(wp: float, chi: float, w: np.ndarray, s_eta: np.ndarray) -> float:
    """Eq. (Jwave-exact): distortion of the real wave, with no chi expansion."""
    jw = 1j * w
    g = (1.0 / (1.0 + chi)) * jw**2 / (jw**2 + wp * math.sqrt(2.0 + chi) * jw + wp**2)
    return float(np.trapezoid(np.abs(g - 1.0) ** 2 * s_eta, w))


def exact_optimum(w: np.ndarray, s_eta: np.ndarray,
                  wp0: float, chi0: float) -> tuple[float, float]:
    """Minimize J_n + J_w with the exact transfer, seeded at the asymptote."""
    from scipy.optimize import minimize

    def cost(v: np.ndarray) -> float:
        wp, chi = math.exp(v[0]), math.exp(v[1])
        return j_noise(wp, chi) + j_wave(wp, chi, w, s_eta)

    res = minimize(cost, [math.log(wp0), math.log(chi0)], method="Nelder-Mead",
                   options=dict(xatol=1e-10, fatol=1e-18,
                                maxiter=20000, maxfev=20000))
    return math.exp(res.x[0]), math.exp(res.x[1])


def coefficients(wp: float, chi: float, tau: float,
                 sigma_aw: float) -> tuple[float, float]:
    """Back out (C_P, C_V) from an optimum at a record's operating point."""
    ts = pseudo_period(tau)
    r_p = math.sqrt(Q_EFF / wp**4 / ts)
    r_v = math.sqrt(Q_EFF / (wp**2 * chi) / ts)
    scale = Q_EFF**0.1 * (sigma_aw / C_SIGMA_OU_II) ** 0.8 / math.sqrt(ts)
    return r_p / (scale * tau**2.4), r_v / (scale * tau**1.4)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exact", action="store_true",
                    help="also solve the unexpanded objective (needs SciPy)")
    ap.add_argument("--points", type=int, default=4001,
                    help="frequency grid points (default 4001)")
    args = ap.parse_args(argv)

    w = frequency_grid(args.points)
    print(f"q_eff = {Q_EFF:.4e} m^2/s^3   c_sigma = {C_SIGMA_OU_II}   "
          f"band = {BAND_HZ[0]}-{BAND_HZ[1]} Hz\n")

    head = (f"{'record':>16}{'tau':>7}{'sigma_aw':>10}{'M_0':>9}{'M_-2':>10}"
            f"{'omega_p*':>10}{'chi*':>9}{'C_P':>8}{'C_V':>8}{'C_P/C_V':>9}")
    if args.exact:
        head += f"{'C_P ex':>9}{'C_V ex':>9}{'ratio ex':>10}"
    print(head)

    cps, cvs, ratios = [], [], []
    cps_x, cvs_x, ratios_x = [], [], []
    for record, family, hs, lam, tau, sigma_aw in SEAS:
        gamma = 3.3 if family == "jonswap" else 1.0
        s_eta = jonswap(w, hs, peak_period(lam), gamma)
        m0, m_2 = moments(s_eta, w)
        wp, chi = asymptotic_optimum(m0, m_2)
        c_p, c_v = coefficients(wp, chi, tau, sigma_aw)
        cps.append(c_p); cvs.append(c_v); ratios.append(c_p / c_v)
        line = (f"{record:>16}{tau:7.3f}{sigma_aw:10.4f}{m0:9.4f}{m_2:10.4f}"
                f"{wp:10.4f}{chi:9.5f}{c_p:8.4f}{c_v:8.4f}{c_p / c_v:9.4f}")
        if args.exact:
            wpx, chix = exact_optimum(w, s_eta, wp, chi)
            cpx, cvx = coefficients(wpx, chix, tau, sigma_aw)
            cps_x.append(cpx); cvs_x.append(cvx); ratios_x.append(cpx / cvx)
            line += f"{cpx:9.4f}{cvx:9.4f}{cpx / cvx:10.4f}"
        print(line)

    def report(name: str, values: list[float]) -> float:
        mean = sum(values) / len(values)
        print(f"  {name:<26}{mean:8.4f}   spread {min(values):.4f}-{max(values):.4f}")
        return mean

    print("\nweak-regularization optimum, Eqs. (wp-opt) and (chi-opt):")
    c_p = report("C_P", cps)
    c_v = report("C_V", cvs)
    ratio = report("C_P/C_V", ratios)
    if args.exact:
        print("\nunexpanded objective, exact |G-1|^2 on the finite band:")
        c_px = report("C_P", cps_x)
        c_vx = report("C_V", cvs_x)
        ratio_x = report("C_P/C_V", ratios_x)
        print(f"\n  expansion cost: C_P {100.0 * (c_px / c_p - 1.0):+.1f} %, "
              f"C_V {100.0 * (c_vx / c_v - 1.0):+.1f} %, "
              f"ratio {100.0 * (ratio_x / ratio - 1.0):+.1f} %")

    print("\nrelative to the empirical schedule, at the same operating points:")
    print(f"  {'record':>16}{'r_p MSE/emp':>13}{'r_v MSE/emp':>13}")
    for record, _family, _hs, _lam, tau, sigma_aw in SEAS:
        ts = pseudo_period(tau)
        norm = math.sqrt(TS0 / ts)
        rp_emp = C_P_EMPIRICAL * sigma_aw * tau * tau * norm
        rv_emp = C_V_EMPIRICAL * sigma_aw * tau * norm
        scale = Q_EFF**0.1 * (sigma_aw / C_SIGMA_OU_II) ** 0.8 / math.sqrt(ts)
        rp_mse = c_p * scale * tau**2.4
        rv_mse = c_v * scale * tau**1.4
        print(f"  {record:>16}{rp_mse / rp_emp:13.3f}{rv_mse / rv_emp:13.3f}")
    print(f"\n  empirical channel ratio c_p/c_v = "
          f"{C_P_EMPIRICAL / C_V_EMPIRICAL:.3f} against the derived {ratio:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
