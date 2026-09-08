#!/usr/bin/env python3
"""Analytical coefficients of the OU-II dual-channel physical-MSE pseudo law.

doc/kalman_ou_ii/ou2-dual-regularization-mse.tex derives the joint
displacement-MSE optimum of the two zero pseudo-measurements (p = 0 and v = 0)
that regularize the OU-II integration chain.  Its result is a pair of *shapes*,

    r_p = C_P q_eff^(1/10) sigma_a,B^(4/5) tau^(12/5) / sqrt(T_S),
    r_v = C_V q_eff^(1/10) sigma_a,B^(4/5) tau^(7/5)  / sqrt(T_S),

This driver evaluates diagnostic coefficients for the shapes using the derivation's
own stationarity conditions on the eight v1.2.1 vessel CG-heave periodograms
at the operating points in the current five-draw OU-II startup study. The diagnostic
does not replace the selected deployed coefficients.

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
import csv
from pathlib import Path
from sim_dataset import input_provenance
from model_mismatch_ablation import RECORDS
from scipy.signal import periodogram

import numpy as np
trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

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
TS_MIN, TS_MAX = 1.0 / 200.0, 0.250

# sigma_aw = c_sigma sigma_a,B is the OU prior; the distortion penalty depends
# on the physical band RMS, so the law divides c_sigma back out.
C_SIGMA_OU_II = 0.85

# Selected physical-MSE coefficients, quoted for comparison only.
C_P_DEPLOYED = 0.1116
C_V_DEPLOYED = C_P_DEPLOYED / 0.4

# The simulator's first-order spectra are band-limited to 0.02-0.8 Hz.
BAND_HZ = (0.02, 0.8)

ROOT = Path(__file__).resolve().parents[1]


def reference_operating_points(path: Path):
    """Current five-draw OU-II operating points; no historical constants."""
    with path.open() as stream:
        rows = [r for r in csv.DictReader(stream)
                if r["family"] == "OU-II" and r["arm"] == "deployed"]
    for record in RECORDS:
        rr = [r for r in rows if r["input"] == record.filename]
        if len(rr) != 5:
            raise ValueError(f"{record.filename}: five current calibration draws required")
        yield (record.filename, np.mean([float(r["tau_applied_s"]) for r in rr]),
               np.mean([float(r["sigma_applied_mps2"]) for r in rr]))


def vessel_spectrum(path: Path):
    """Finite-record CG-heave periodogram, density per rad/s.

    This diagnostic uses the actual vessel output, including its response
    attenuation. It is not an incident-elevation spectrum or a continuum proof.
    """
    t, z = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(0, 3), unpack=True)
    f, density_hz = periodogram(z, fs=1.0 / np.median(np.diff(t)), detrend="constant")
    band = (f >= BAND_HZ[0]) & (f <= BAND_HZ[1])
    return 2.0 * math.pi * f[band], density_hz[band] / (2.0 * math.pi)


def pseudo_period(tau: float) -> float:
    """Deployed self-similar pseudo-update cadence with its clamps."""
    return min(max(C_T * tau, TS_MIN), TS_MAX)


def frequency_grid(points: int = 4001) -> np.ndarray:
    return np.logspace(math.log10(2.0 * math.pi * BAND_HZ[0]),
                       math.log10(2.0 * math.pi * BAND_HZ[1]), points)


def moments(s_eta: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """(M_0, M_-2): displacement variance and its second negative moment.

    S_eta denotes CG heave here. Incident elevation is filtered by the vessel
    response before these one-sided integrals are evaluated.
    """
    return float(trapezoid(s_eta, w)), float(trapezoid(s_eta / w**2, w))


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
    return float(trapezoid(np.abs(g - 1.0) ** 2 * s_eta, w))


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
    ap.add_argument("--data-dir", type=Path, default=ROOT / "plots/kalman_ou_ii")
    ap.add_argument("--operating-points", type=Path,
                    default=ROOT / "reports/results/startup_ablation/startup_runs.csv")
    args = ap.parse_args(argv)

    seas = list(reference_operating_points(args.operating_points))
    input_provenance(args.data_dir / record for record, _, _ in seas)
    print(f"q_eff = {Q_EFF:.4e} m^2/s^3   c_sigma = {C_SIGMA_OU_II}   "
          f"band = {BAND_HZ[0]}-{BAND_HZ[1]} Hz\n")

    head = (f"{'record':>16}{'tau':>7}{'sigma_aw':>10}{'M_0':>9}{'M_-2':>10}"
            f"{'omega_p*':>10}{'chi*':>9}{'C_P':>8}{'C_V':>8}{'C_P/C_V':>9}")
    if args.exact:
        head += f"{'C_P ex':>9}{'C_V ex':>9}{'ratio ex':>10}"
    print(head)

    cps, cvs, ratios = [], [], []
    cps_x, cvs_x, ratios_x = [], [], []
    for record, tau, sigma_aw in seas:
        w, s_eta = vessel_spectrum(args.data_dir / record)
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

    print("\nrelative to the selected physical-MSE schedule at the same operating points:")
    print(f"  asymptotic C_P/selected = {c_p / C_P_DEPLOYED:.3f}, "
          f"C_V/selected = {c_v / C_V_DEPLOYED:.3f}")
    print(f"  selected channel ratio = {C_P_DEPLOYED / C_V_DEPLOYED:.3f}, "
          f"diagnostic mean ratio = {ratio:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
