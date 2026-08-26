#!/usr/bin/env python3
"""Per-axis measurement of the two quantities that set the MSE-optimal R_S.

Section "Riccati interpretation of the S=0 regularizer" of the OU-III paper
carries the bias--variance optimum for the vertical channel only.  Its result,

    MSE_p(w_R) = 3 q_eff / (2 w_R^3)  +  4 m_{-4} w_R^4,
    w_R*^7 = (9/32) q_eff / m_{-4},
    r_S* = (32/9)^{3/7} q_eff^{1/14} m_{-4}^{3/7} T_S^{-1/2},

is not intrinsically vertical: nothing in the derivation of either term uses
the axis.  Reading it per axis i in {x, y, z} gives

    r_S*_i / r_S*_j = (q_i/q_j)^{1/14} (m_i/m_j)^{3/7},

so the anisotropy the filter should deploy is fixed by two measurable
quantities and no new calibration constant.  The exponents are wildly
unequal: the drift-noise intensity enters with weight 1/14 and the true
displacement moment with weight 3/7, six times larger.  Observing that x and y
carry more acceleration error is therefore *not* on its own an argument for a
tighter horizontal anchor.

This driver measures both, per axis, from simulator output rather than
inferring them from the IMU white-noise floor.

    q_eff,i = S_eps_i(0),   eps_i(t) = acc_est_i(t) - acc_ref_i(t)

is the zero-frequency intensity of the *posterior acceleration error*, which
the paper notes is dominated by attitude error through gravity and by residual
accelerometer bias, not by the sensor floor.  It is read off a Lorentzian
    S_eps(w) ~ 2 P mu / (mu^2 + w^2)
fitted over the drift band, i.e. below the sea's own energy.

    m_{-4,i} = int S_p_i^truth(w) / w^4 dw

is the negative fourth moment of the *true* displacement spectrum -- the real
motion the S_i = 0 constraint destroys.  It is taken from the reference
columns, never the estimate, and band-limited at the bottom edge of the wave
band because the 1/w^4 weight would otherwise amplify the estimator floor.

Typical use:

    python3 tools/ou_axis_rs_optimum.py
    python3 tools/ou_axis_rs_optimum.py --glob 'tests/kalman_ou_ii/w3d_*_ou2.csv'
    python3 tools/ou_axis_rs_optimum.py --csv-out axis_rs.csv
"""

from __future__ import annotations

import argparse
import csv
import glob as globmod
import math
import os
import re
import sys

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Fraction of reference displacement power left below the wave band.  The
# bottom edge found this way is where the sea's own energy starts; below it the
# synthetic records carry nothing, so it is also the honest floor for the
# 1/w^4 integrand.
BAND_FLOOR_FRACTION = 0.005
BAND_CEIL_FRACTION = 0.995

# Welch settings.  300 s segments at 75 percent overlap keep ~9 averages inside
# a 900 s window while still resolving the drift band below the sea.
SEGMENT_SECONDS = 300.0
OVERLAP = 0.75

AXES = ("x", "y", "z")


def welch_psd(x: np.ndarray, fs: float, nseg: int,
              overlap: float = OVERLAP) -> tuple[np.ndarray, np.ndarray]:
    """One-sided PSD in units per rad/s, so int S dw recovers the variance."""
    win = np.hanning(nseg)
    step = max(1, int(round(nseg * (1.0 - overlap))))
    starts = range(0, len(x) - nseg + 1, step)
    acc = np.zeros(nseg // 2 + 1)
    count = 0
    for s in starts:
        seg = x[s:s + nseg]
        spec = np.fft.rfft((seg - seg.mean()) * win)
        acc += np.abs(spec) ** 2
        count += 1
    if count == 0:
        raise ValueError("window shorter than one Welch segment")
    # Per-Hz one-sided, then per rad/s.
    psd_hz = 2.0 * acc / (count * fs * np.sum(win ** 2))
    freqs = np.fft.rfftfreq(nseg, 1.0 / fs)
    return 2.0 * math.pi * freqs, psd_hz / (2.0 * math.pi)


def wave_band(w: np.ndarray, psd: np.ndarray) -> tuple[float, float]:
    """Bottom and top edge of the band carrying the reference energy."""
    total = float(np.trapezoid(psd, w))
    if total <= 0.0:
        return float(w[1]), float(w[-1])
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (psd[1:] + psd[:-1]) * np.diff(w)))) / total
    lo = float(np.interp(BAND_FLOOR_FRACTION, cum, w))
    hi = float(np.interp(BAND_CEIL_FRACTION, cum, w))
    return lo, hi


def fit_lorentzian(w: np.ndarray, psd: np.ndarray) -> tuple[float, float]:
    """Fit S(w) = 2 P mu / (mu^2 + w^2) over a drift band; return (S(0), mu).

    1/S is linear in w^2 with slope 1/(2 P mu) and intercept mu/(2 P), so the
    fit is a straight line in w^2 and needs no iteration.  A non-positive or
    non-finite slope means the band shows no roll-off; fall back to the plateau
    mean, which is what S(0) degenerates to for mu well above the band.
    """
    good = np.isfinite(psd) & (psd > 0.0)
    if good.sum() < 4:
        return float("nan"), float("nan")
    x = w[good] ** 2
    y = 1.0 / psd[good]
    slope, intercept = np.polyfit(x, y, 1)
    if not (np.isfinite(slope) and np.isfinite(intercept)) \
            or slope <= 0.0 or intercept <= 0.0:
        return float(np.mean(psd[good])), float("inf")
    mu = math.sqrt(intercept / slope)
    return float(1.0 / intercept), float(mu)


def negative_fourth_moment(w: np.ndarray, psd: np.ndarray, w_lo: float,
                           w_hi: float) -> float:
    """int S_p(w)/w^4 dw over the band the sea actually occupies."""
    sel = (w >= w_lo) & (w <= w_hi)
    if sel.sum() < 2:
        return float("nan")
    return float(np.trapezoid(psd[sel] / w[sel] ** 4, w[sel]))


def _deployed_rho_xy() -> float:
    """The wrapper's deployed R_S_xy_factor_, so the comparison line cannot go stale."""
    src = os.path.join(REPO_ROOT, "src", "kalman_ou_iii", "SeaStateFusionFilter_OU_III.h")
    m = re.search(r"float\s+R_S_xy_factor_\s*=\s*([0-9.eE+-]+)f\s*;",
                  open(src, encoding="utf-8").read())
    return float(m.group(1)) if m else float("nan")


def scenario_label(path: str) -> str:
    base = os.path.basename(path)
    base = re.sub(r"^w3d_", "", base)
    base = re.sub(r"_fusion.*$", "", base)
    m = re.match(r"([a-z]+)_H([0-9.]+)_", base)
    return f"{m.group(1)} Hs={float(m.group(2)):.2f}" if m else base


def analyse(path: str, window_s: float) -> dict | None:
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size == 0 or "time" not in (data.dtype.names or ()):
        return None
    t = data["time"]
    mask = t >= t[-1] - window_s
    t = t[mask]
    if len(t) < 16:
        return None
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    nseg = int(round(SEGMENT_SECONDS * fs))
    if nseg > len(t):
        nseg = 1 << int(math.floor(math.log2(len(t))))

    out: dict = {"scenario": scenario_label(path), "path": path}
    for ax in AXES:
        acc_err = data["acc_est_" + ax][mask] - data["acc_ref_" + ax][mask]
        disp_ref = data["disp_ref_" + ax][mask]

        w, psd_ref = welch_psd(disp_ref, fs, nseg)
        w_lo, w_hi = wave_band(w, psd_ref)

        _, psd_err = welch_psd(acc_err, fs, nseg)
        drift = (w > 0.0) & (w < w_lo)
        q_eff, mu = fit_lorentzian(w[drift], psd_err[drift])
        plateau = float(np.mean(psd_err[drift])) if drift.sum() else float("nan")

        out[ax] = {
            "q_eff": q_eff,
            "q_plateau": plateau,
            "mu": mu,
            "m4": negative_fourth_moment(w, psd_ref, w_lo, w_hi),
            "band_lo_hz": w_lo / (2.0 * math.pi),
            "band_hi_hz": w_hi / (2.0 * math.pi),
            "acc_err_rms": float(np.sqrt(np.mean(acc_err ** 2))),
            "disp_ref_rms": float(np.sqrt(np.mean(disp_ref ** 2))),
            "n_drift_bins": int(drift.sum()),
        }
    return out


def rs_ratio(q_i: float, q_z: float, m_i: float, m_z: float) -> float:
    """r_S*_i / r_S*_z from the per-axis optimum: q^{1/14} m^{3/7}."""
    if not all(map(np.isfinite, (q_i, q_z, m_i, m_z))) or q_z <= 0 or m_z <= 0:
        return float("nan")
    return (q_i / q_z) ** (1.0 / 14.0) * (m_i / m_z) ** (3.0 / 7.0)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--glob", default="tests/kalman_ou_iii/w3d_*_fusion_ou3.csv",
                    help="fusion timeseries to analyse")
    ap.add_argument("--window-sec", type=float, default=900.0,
                    help="trailing window scored, matching the validation study")
    ap.add_argument("--csv-out", default=None, help="write the table as CSV")
    args = ap.parse_args(argv)

    paths = sorted(globmod.glob(args.glob))
    if not paths:
        print(f"no files matched {args.glob!r}", file=sys.stderr)
        return 1

    rows = []
    for path in paths:
        res = analyse(path, args.window_sec)
        if res is not None:
            rows.append(res)
    if not rows:
        print("no usable records", file=sys.stderr)
        return 1

    print("Per-axis drift-noise intensity and true-displacement moment")
    print(f"  window {args.window_sec:.0f} s, Welch {SEGMENT_SECONDS:.0f} s "
          f"segments at {OVERLAP:.0%} overlap")
    print()
    head = (f"{'scenario':22}{'band lo':>9}"
            + "".join(f"{'q_' + a:>11}" for a in AXES)
            + "".join(f"{'m4_' + a:>11}" for a in AXES))
    print(head)
    print("-" * len(head))
    for r in rows:
        print(f"{r['scenario']:22}{r['z']['band_lo_hz']:9.4f}"
              + "".join(f"{r[a]['q_eff']:11.3e}" for a in AXES)
              + "".join(f"{r[a]['m4']:11.3e}" for a in AXES))

    print()
    head2 = (f"{'scenario':22}"
             + f"{'qx/qz':>9}{'qy/qz':>9}{'m4x/m4z':>10}{'m4y/m4z':>10}"
             + f"{'rSx/rSz':>10}{'rSy/rSz':>10}")
    print(head2)
    print("-" * len(head2))
    for r in rows:
        qx = r["x"]["q_eff"] / r["z"]["q_eff"]
        qy = r["y"]["q_eff"] / r["z"]["q_eff"]
        mx = r["x"]["m4"] / r["z"]["m4"]
        my = r["y"]["m4"] / r["z"]["m4"]
        print(f"{r['scenario']:22}{qx:9.2f}{qy:9.2f}{mx:10.3f}{my:10.3f}"
              f"{rs_ratio(r['x']['q_eff'], r['z']['q_eff'], r['x']['m4'], r['z']['m4']):10.3f}"
              f"{rs_ratio(r['y']['q_eff'], r['z']['q_eff'], r['y']['m4'], r['z']['m4']):10.3f}")

    def geo(vals):
        vals = [v for v in vals if np.isfinite(v) and v > 0]
        return math.exp(sum(map(math.log, vals)) / len(vals)) if vals else float("nan")

    print()
    print("geometric mean over scenarios:")
    for ax in ("x", "y"):
        print(f"  {ax}: q/qz={geo([r[ax]['q_eff'] / r['z']['q_eff'] for r in rows]):6.2f}"
              f"  m4/m4z={geo([r[ax]['m4'] / r['z']['m4'] for r in rows]):6.3f}"
              f"  rS*/rS*z={geo([rs_ratio(r[ax]['q_eff'], r['z']['q_eff'], r[ax]['m4'], r['z']['m4']) for r in rows]):6.3f}")
    print(f"  deployed: rSx/rSz = rSy/rSz = {_deployed_rho_xy():.3f} (R_S_xy_factor)")

    if args.csv_out:
        with open(args.csv_out, "w", newline="") as fh:
            wtr = csv.writer(fh)
            wtr.writerow(["scenario", "axis", "q_eff", "q_plateau", "mu",
                          "m4", "band_lo_hz", "band_hi_hz", "acc_err_rms",
                          "disp_ref_rms", "rs_ratio_vs_z"])
            for r in rows:
                for ax in AXES:
                    wtr.writerow([r["scenario"], ax, r[ax]["q_eff"],
                                  r[ax]["q_plateau"], r[ax]["mu"], r[ax]["m4"],
                                  r[ax]["band_lo_hz"], r[ax]["band_hi_hz"],
                                  r[ax]["acc_err_rms"], r[ax]["disp_ref_rms"],
                                  rs_ratio(r[ax]["q_eff"], r["z"]["q_eff"],
                                           r[ax]["m4"], r["z"]["m4"])])
        print(f"\nwrote {args.csv_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
