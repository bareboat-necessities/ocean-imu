#!/usr/bin/env python3
"""Reproduce six *illustrative* fractional marine-motion figures.

Run: python3 doc/kalman_ou_iii/fractional-marine-motion-figures.py
Requires: numpy, scipy, matplotlib. Outputs vector PDFs and a numerical JSON
summary beside this script. No IMU/INS estimator, vessel measurements, or
hardware benchmarks are invoked.

Conventions: one-sided spectral densities S(f) in m^2/Hz, f>0, with
Var[p] = integral_0^inf S(f)df. The transfer H(i*omega) uses omega=2*pi*f;
the continuous-time two-sided PSD is |H|^2 for unit-intensity white forcing,
and S(f)=2*|H(i*2*pi*f)|^2 in the corresponding one-sided Hz convention.
Because models are rescaled by a common integrated displacement variance, the
irrelevant common white-noise intensity factor cancels for all plots.

Fitting: minimize weighted RMS log-PSD residual in dB over the explicitly
defined finite fitting band, after imposing the same integrated variance over
the numerical quadrature band. All fits are deterministic, not eye-aligned.
The existing C41J/C52J parameters are read from docs/carma-three-axis.md
(main at 2026-10-29), and are *not* re-optimized in Figure 1.

Figure 5 uses a stable/minimum-phase, *constrained* rational approximant:
    R_M(s) = k*s/L(s)^2 * product((s+z_i)/(s+r_i)), M=4,...,7,
with r_i,z_i>0 and m=M-4. This is NOT general vector fitting; it is a
reproducible, deliberately conservative state-budget experiment. It preserves
the DC zero, relative degree three, and a stable four-state oscillator core.
The metric is complex transfer relative error over the declared band, with
the corresponding PSD log error also reported.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

OUT = Path(__file__).resolve().parent
RNG_SEED = 20261029
TWOPI = 2.0 * np.pi
GRID = np.geomspace(0.008, 12.0, 13000)  # numerical variance quadrature, Hz
FIT1 = np.geomspace(0.18, 1.15, 390)  # Figure 1: energetic incident-wave band
FIT6 = np.geomspace(0.051, 0.66, 440)  # Figure 6: surrogate heave response
APPROX_BAND = (0.06, 1.50)  # Figure 5: transfer-function approximation, Hz
HS = 1.0  # Figure 1: incident significant wave height, m

COL = dict(blue="#2767a7", orange="#de792c", green="#25846a",
           purple="#8763a9", red="#c8474c", gray="#5f6773")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.6,
    "axes.labelsize": 9.0, "axes.titlesize": 9.4,
    "xtick.labelsize": 7.7, "ytick.labelsize": 7.7,
    "legend.fontsize": 7.4, "lines.linewidth": 1.65,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.alpha": 0.19, "grid.linewidth": 0.5,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.07,
})


def L(s: np.ndarray, f0: float, zeta: float) -> np.ndarray:
    w0 = TWOPI * f0
    return s**2 + 2.0 * zeta * w0 * s + w0**2


def frac(f: np.ndarray, nu: float, f0: float, zeta: float) -> np.ndarray:
    s = 1j * TWOPI * f
    return s / L(s, f0, zeta)**nu


def c41(f: np.ndarray, f1: float, z1: float, f2: float, z2: float) -> np.ndarray:
    s = 1j * TWOPI * f
    return s / (L(s, f1, z1) * L(s, f2, z2))


def c52(f: np.ndarray, f1: float, z1: float, f2: float, z2: float,
        flambda: float, fz: float, zz: float) -> np.ndarray:
    s = 1j * TWOPI * f
    wz = TWOPI * fz
    return (s**2 + 2 * zz * wz * s + wz**2) / (
        (s + TWOPI * flambda) * L(s, f1, z1) * L(s, f2, z2))


def positive_shape(h: np.ndarray) -> np.ndarray:
    return np.abs(h)**2


def to_variance(shape: np.ndarray, f: np.ndarray, variance: float) -> np.ndarray:
    return np.maximum(shape, 0.0) * variance / np.trapezoid(shape, f)


def wave_shape(f: np.ndarray, fp: float, gamma: float = 1.0) -> np.ndarray:
    # PM/JONSWAP peak-period formulation. Common coefficient cancels on m0 fit.
    f = np.maximum(f, 1e-8)
    sigma = np.where(f <= fp, 0.07, 0.09)
    r = np.exp(-0.5 * ((f / fp - 1.0) / sigma)**2)
    return f**(-5) * np.exp(-1.25 * (fp / f)**4) * gamma**r


def wave(f: np.ndarray, fp: float, gamma: float, hs: float) -> np.ndarray:
    a = (hs / 4.0)**2 / np.trapezoid(wave_shape(GRID, fp, gamma), GRID)
    return a * wave_shape(f, fp, gamma)


def scaled_model(fn, args: tuple, f: np.ndarray, variance: float) -> np.ndarray:
    scale = variance / np.trapezoid(positive_shape(fn(GRID, *args)), GRID)
    return scale * positive_shape(fn(f, *args))


def weights(target: np.ndarray) -> np.ndarray:
    # Avoid giving the near-zero tail as much influence as the physical peak.
    return 0.20 + 0.80 * np.minimum(target / np.max(target), 1.0)**0.30


def log_resid(model: np.ndarray, target: np.ndarray, w: np.ndarray) -> np.ndarray:
    return np.sqrt(w) * (10 * np.log10(np.maximum(model, 1e-100)
                                         / np.maximum(target, 1e-100)))


def rms_db(model: np.ndarray, target: np.ndarray, w: np.ndarray) -> float:
    r = log_resid(model, target, w)
    return float(np.sqrt(np.sum(r * r) / np.sum(w)))


def fit_fractional(nu: float, target_fn, band: np.ndarray,
                   variance: float, bounds=((0.07, 0.04), (0.65, 0.95))):
    target = target_fn(band)
    w = weights(target)

    def residual(x):
        return log_resid(scaled_model(frac, (nu, *x), band, variance), target, w)

    best = None
    for p0 in ((0.28, 0.15), (0.27, 0.30), (0.36, 0.12), (0.17, 0.35)):
        result = least_squares(residual, p0, bounds=bounds, max_nfev=350,
                               ftol=1e-10, xtol=1e-10, gtol=1e-10)
        if best is None or np.linalg.norm(result.fun) < np.linalg.norm(best.fun):
            best = result
    return best.x, rms_db(scaled_model(frac, (nu, *best.x), band, variance),
                          target, w)


def save(fig: plt.Figure, index: int, shortname: str):
    fig.savefig(OUT / f"fractional-fig-{index}-{shortname}.pdf")
    plt.close(fig)


def fig1(summary: dict):
    jons = lambda f: wave(f, 0.29, 3.3, HS)
    var = (HS / 4)**2
    # Actual constrained parameters from the CARMA three-axis research report.
    args41 = (0.306945, 0.104191, 0.815318, 0.478188)
    args52 = (0.27812, 0.09715, 0.44274, 0.52299, 1.86402, 0.15105, 0.20)
    allfits = {}
    for nu in (1.80, 1.85, 1.90, 2.00):
        params, err = fit_fractional(nu, jons, FIT1, var)
        allfits[f"{nu:.2f}"] = dict(f0_Hz=float(params[0]), zeta=float(params[1]),
                                    weighted_RMS_dB=err)
    summary["incident_JONSWAP_fractional_fits"] = allfits
    summary["existing_C41J"] = dict(f1_Hz=args41[0], zeta1=args41[1],
                                     f2_Hz=args41[2], zeta2=args41[3],
                                     weighted_RMS_dB=rms_db(
                                         scaled_model(c41, args41, FIT1, var),
                                         jons(FIT1), weights(jons(FIT1))))
    summary["existing_C52J"] = dict(f1_Hz=args52[0], zeta1=args52[1],
                                     f2_Hz=args52[2], zeta2=args52[3],
                                     f_lambda_Hz=args52[4], f_zero_Hz=args52[5],
                                     zeta_zero=args52[6],
                                     weighted_RMS_dB=rms_db(
                                         scaled_model(c52, args52, FIT1, var),
                                         jons(FIT1), weights(jons(FIT1))))
    f = np.linspace(0.055, 1.25, 1500)
    fig, ax = plt.subplots(1, 2, figsize=(7.12, 3.05), sharey=True)
    for a in ax:
        a.set_yscale("log"); a.set_xlim(0.06, 1.22)
        a.set_ylim(1e-5, 1.8); a.grid(True, which="major")
        a.set_xlabel("Frequency $f$ (Hz)")
    ax[0].set_ylabel(r"One-sided displacement PSD ($\mathrm{m^2/Hz}$)")
    ax[0].plot(f, jons(f), color="black", lw=2.1, label=r"JONSWAP, $\gamma=3.3$")
    ax[0].plot(f, wave(f, 0.29, 1.0, HS), color=COL["gray"], ls="--", label="PM")
    ax[0].plot(f, scaled_model(c41, args41, f, var), color=COL["blue"],
               label="Existing C41J")
    ax[0].plot(f, scaled_model(c52, args52, f, var), color=COL["red"],
               label="Existing C52J")
    ax[0].set_title("(a) Incident-wave targets and CARMA fits")
    ax[0].legend(loc="upper right", frameon=False, fontsize=7.0)
    ax[1].plot(f, jons(f), color="black", lw=2.1, label="JONSWAP target")
    for nu, c in zip((1.80, 1.85, 1.90, 2.00),
                     (COL["green"], COL["blue"], COL["purple"], COL["orange"])):
        fit = allfits[f"{nu:.2f}"]
        ax[1].plot(f, scaled_model(frac, (nu, fit["f0_Hz"], fit["zeta"]), f, var),
                   label=f"Fractional $\\nu={nu:.2f}$", color=c)
    ax[1].set_title("(b) Fractional fits; $H_s=1$ m")
    ax[1].legend(loc="upper right", frameon=False, fontsize=7.0)
    fig.subplots_adjust(wspace=0.12)
    save(fig, 1, "spectra")
    return allfits


def fig2(summary: dict):
    f = np.linspace(0.06, 1.5, 1600)
    nus = (1.80, 1.85, 1.90, 2.00)
    palette = (COL["green"], COL["blue"], COL["purple"], COL["orange"])
    f0, zeta = 0.29, 0.22
    fig, ax = plt.subplots(1, 2, figsize=(7.12, 2.9))
    bandwidth = {}
    for nu, col in zip(nus, palette):
        s = positive_shape(frac(f, nu, f0, zeta))
        peak = s.max(); db = 10 * np.log10(np.maximum(s / peak, 1e-30))
        keep = np.flatnonzero(s >= peak * 0.5)
        bw = f[keep[-1]] - f[keep[0]]
        bandwidth[f"{nu:.2f}"] = float(bw)
        ax[0].plot(f, db, color=col, label=f"$\\nu={nu:.2f}$")
    ax[0].set(xlim=(0.06, 1.45), ylim=(-65, 2),
              xlabel="Frequency (Hz)", ylabel="PSD (dB relative to own peak)",
              title=r"(a) Fixed $f_0=0.29$ Hz and $\zeta=0.22$")
    ax[0].grid(True); ax[0].legend(loc="upper right", frameon=False)
    ax[1].plot(list(map(float, bandwidth)), list(bandwidth.values()),
               color=COL["blue"], marker="o", markersize=5)
    for k, val in bandwidth.items():
        ax[1].annotate(f"{val:.3f}", (float(k), val), xytext=(0, 7),
                       textcoords="offset points", ha="center", fontsize=7)
    ax[1].set(xlabel="Fractional order $\\nu$", ylabel="Half-power bandwidth (Hz)",
              xlim=(1.77, 2.025), title="(b) Width of the principal peak")
    ax[1].grid(True)
    fig.subplots_adjust(wspace=0.33)
    save(fig, 2, "order-sensitivity")
    summary["fixed_parameter_half_power_bandwidth_Hz"] = bandwidth


def fig3(summary: dict):
    f = np.geomspace(0.001, 110, 3500)
    nu, f0, zeta = 1.85, 0.29, 0.22
    p = positive_shape(frac(f, nu, f0, zeta))
    a = (TWOPI * f)**4 * p
    fig, ax = plt.subplots(1, 2, figsize=(7.12, 2.86))
    slopes = (2.0, 2.0 - 4.0 * nu, 6.0, 6.0 - 4.0 * nu)
    for i, (v, low, high, title, col) in enumerate(((p, slopes[0], slopes[1],
        "(a) Displacement spectrum", COL["blue"]), (a, slopes[2], slopes[3],
        "(b) Acceleration spectrum", COL["red"]))):
        v = v / v.max()
        ax[i].loglog(f, v, color=col, lw=2, label="Full model")
        for fanchor, slope, interval, name in ((0.0035, low, (0.001, 0.03), "low"),
                                              (24., high, (8., 110.), "high")):
            idx = np.argmin(np.abs(f-fanchor))
            ff = np.geomspace(interval[0], interval[1], 90)
            ax[i].loglog(ff, v[idx] * (ff/f[idx])**slope,
                         color=COL["gray"], ls="--", lw=1.3,
                         label=f"${{f}}^{{{slope:g}}}$ {name}-$f$")
        ax[i].set(xlim=(0.001, 110), xlabel="Frequency (Hz)",
                  ylabel="PSD relative to its maximum", title=title)
        ax[i].grid(True, which="both"); ax[i].legend(frameon=False,
                                                      loc="lower center", fontsize=7.0)
    fig.subplots_adjust(wspace=0.27)
    save(fig, 3, "asymptotic-slopes")
    summary["asymptotic_loglog_slopes_nu_1p85"] = dict(
        displacement_low=slopes[0], displacement_high=slopes[1],
        acceleration_low=slopes[2], acceleration_high=slopes[3])


def fig4(fits: dict, summary: dict):
    fs, duration = 200.0, 240.0
    n = int(fs * duration)
    f = np.fft.rfftfreq(n, 1/fs)
    positive_f = f.copy(); positive_f[0] = 1e-12
    rng = np.random.default_rng(RNG_SEED)
    common = np.fft.rfft(rng.standard_normal(n))
    var = (HS / 4)**2
    cases = [
        ("CARMA(4,1), existing JONSWAP fit", c41,
         (0.306945, 0.104191, 0.815318, 0.478188), COL["blue"]),
        ("Fractional, $\\nu=1.85$", frac,
         (1.85, fits["1.85"]["f0_Hz"], fits["1.85"]["zeta"]), COL["green"]),
        ("CARMA(5,2), existing JONSWAP fit", c52,
         (0.27812, 0.09715, 0.44274, 0.52299, 1.86402, 0.15105, 0.20),
         COL["red"]),
    ]
    fig, ax = plt.subplots(3, 1, figsize=(7.12, 3.78), sharex=True, sharey=True)
    t = np.arange(n) / fs
    for axis, (title, fn, args, color) in zip(ax, cases):
        S = scaled_model(fn, args, positive_f, var)
        S[0] = 0.0  # all shown processes have zero DC
        # Real-white-noise FFT filter: white one-sided PSD = 2/fs, so the
        # magnitude multiplier sqrt(S*fs/2) gives the requested output PSD.
        fourier = common * np.sqrt(S * fs / 2.0)
        fourier[-1] = 0.0  # Nyquist bin, irrelevant at 200 Hz for these seas
        trace = np.fft.irfft(fourier, n=n)
        keep = (t >= 25) & (t <= 95)
        axis.plot(t[keep] - 25, trace[keep], color=color, lw=1.15)
        axis.set_ylabel("$p$ (m)"); axis.set_ylim(-0.85, 0.85)
        axis.set_title(title, loc="left", fontsize=8.7, pad=2)
        axis.grid(True)
    ax[-1].set_xlabel("Time since displayed window start (s)")
    fig.subplots_adjust(hspace=0.39, left=0.09, right=0.98)
    save(fig, 4, "realizations")
    summary["time_synthesis"] = dict(method="FFT filtering of one shared seeded Gaussian white-noise record; periodic grid",
        seed=RNG_SEED, sampling_Hz=fs, total_duration_s=duration,
        displayed_interval_s=[25, 95], theoretical_target_variance_m2=var,
        filter_approximation="spectral grid; no finite rational approximation")


def rational(f: np.ndarray, f0: float, zeta: float, gain: float,
             log_zeros: np.ndarray, log_poles: np.ndarray) -> np.ndarray:
    s = 1j * TWOPI * f
    H = gain * s / L(s, f0, zeta)**2
    for z, r in zip(log_zeros, log_poles):
        H = H * (s + np.exp(z)) / (s + np.exp(r))
    return H


def fig5(summary: dict):
    nu, f0, zeta = 1.85, 0.29, 0.22
    band = np.geomspace(APPROX_BAND[0], APPROX_BAND[1], 500)
    target = frac(band, nu, f0, zeta)
    log_g = np.log((TWOPI * f0)**(2 * (2 - nu)))
    results = {}
    fig, ax = plt.subplots(1, 2, figsize=(7.12, 2.95))
    rng = np.random.default_rng(RNG_SEED + 5)
    for m in range(4):
        order = 4 + m
        best = None
        for trial in range(5 if m else 1):
            rates = TWOPI * np.geomspace(0.10, 1.20, max(1, m))[:m]
            starts = (np.log(rates * (1.0 + 0.16*rng.standard_normal(m)))
                      if m else np.empty(0))
            x0 = np.r_[log_g, starts, starts + 0.09]

            def residual(x):
                r = rational(band, f0, zeta, np.exp(x[0]), x[1:1+m],
                             x[1+m:1+2*m]) / target - 1.0
                return np.r_[r.real, r.imag]

            lower = np.r_[-12., np.repeat(np.log(TWOPI*0.016), 2*m)]
            upper = np.r_[12., np.repeat(np.log(TWOPI*12.), 2*m)]
            fit = least_squares(residual, np.clip(x0, lower+1e-5, upper-1e-5),
                                bounds=(lower, upper), max_nfev=1300,
                                ftol=5e-11, xtol=5e-11, gtol=5e-11)
            if best is None or np.linalg.norm(fit.fun) < np.linalg.norm(best.fun):
                best = fit
        x = best.x
        ratio = rational(band, f0, zeta, np.exp(x[0]), x[1:1+m],
                         x[1+m:1+2*m]) / target
        complex_rms = 100 * float(np.sqrt(np.mean(np.abs(ratio-1)**2)))
        psd_err = 20 * np.log10(np.maximum(np.abs(ratio), 1e-15))
        results[str(order)] = dict(complex_relative_RMS_percent=complex_rms,
            log_PSD_RMS_dB=float(np.sqrt(np.mean(psd_err**2))),
            max_abs_PSD_error_dB=float(np.max(np.abs(psd_err))),
            gain=float(np.exp(x[0])), zeros_rad_s=list(np.exp(x[1:1+m])),
            poles_rad_s=list(np.exp(x[1+m:1+2*m])))
        ax[1].semilogx(band, psd_err, label=f"{order} states",
                       color=(COL["orange"], COL["green"], COL["purple"], COL["blue"])[m])
    orders = list(map(int, results))
    ax[0].plot(orders, [results[str(o)]["log_PSD_RMS_dB"] for o in orders],
               color=COL["blue"], marker="o", markersize=5, label="RMS PSD error")
    ax[0].plot(orders, [results[str(o)]["max_abs_PSD_error_dB"] for o in orders],
               color=COL["red"], marker="s", markersize=4, ls="--",
               label="Max. absolute PSD error")
    ax[0].set(xticks=orders, xlabel="Rational wave states per axis",
              ylabel="Log-spectral error (dB)",
              title="(a) Approximation error vs. state count")
    ax[0].legend(frameon=False, fontsize=7.0); ax[0].grid(True)
    ax[1].axhline(0, color="black", lw=0.8)
    ax[1].set(xlim=APPROX_BAND, xlabel="Frequency (Hz)",
              ylabel="Approx. PSD / target (dB)",
              title="(b) Finite-band shape errors")
    ax[1].legend(frameon=False, fontsize=7); ax[1].grid(True)
    fig.subplots_adjust(wspace=0.28)
    save(fig, 5, "rational-order")
    summary["rational_approximation"] = dict(target=dict(nu=nu, f0_Hz=f0,
        zeta=zeta), fitting_band_Hz=APPROX_BAND,
        objective="unweighted complex relative transfer-function least squares",
        constrained_pole_zero_family="s/L(s)^2 times m stable real pole-zero ratios",
        model_orders=results)


def fig6(summary: dict):
    # A *specified illustrative* scalar response amplitude operator, NOT an
    # identified hull response. No claim of a physically calibrated vessel.
    def heave_RAO2(f):
        return 1.0 / (1.0 + (f/0.70)**4)

    def truth(f):
        swell = wave(f, 0.095, 1.0, 1.25)
        wind = wave(f, 0.29, 3.3, 1.00)
        return (swell + wind) * heave_RAO2(f)

    variance = np.trapezoid(truth(GRID), GRID)
    target = truth(FIT6)
    w = weights(target)

    def score(y):
        return rms_db(y, target, w)

    def optimize(which):
        if which == "single_integer":
            fn = lambda x, f: scaled_model(frac, (2.0, x[0], x[1]), f, variance)
            lo, hi = [0.06, 0.035], [0.43, 0.8]
            initials = [[0.095, 0.12], [0.28, 0.2], [0.18, 0.35]]
        elif which == "single_fractional":
            fn = lambda x, f: scaled_model(frac, (x[2], x[0], x[1]), f, variance)
            lo, hi = [0.06, 0.035, 1.76], [0.43, 0.8, 2.2]
            initials = [[0.095, 0.13, 1.85], [0.28, 0.2, 1.85], [0.18, 0.3, 2.0]]
        else:
            # Independent rational CARMA(4,1) sources: *two* 4-state
            # repeated-oscillator blocks; PSDs add, not transfer functions.
            def fn(x, f):
                a = x[4]
                return variance * (a * scaled_model(frac, (2.0, x[0], x[1]), f, 1.0)
                   + (1-a) * scaled_model(frac, (2.0, x[2], x[3]), f, 1.0))
            # Assigning a fixed swell/wind identity avoids equivalent
            # permutations that can confuse optimization and reporting.
            lo = [0.065, 0.035, 0.19, 0.035, 0.04]
            hi = [0.17, 0.70, 0.42, 0.7, 0.96]
            initials = [[0.095, 0.20, 0.29, 0.16, 0.60],
                        [0.11, 0.30, 0.28, 0.25, 0.60]]
        best = None
        for x0 in initials:
            res = least_squares(lambda x: log_resid(fn(x, FIT6), target, w),
                x0, bounds=(lo, hi), max_nfev=650,
                ftol=1e-10, xtol=1e-10, gtol=1e-10)
            if best is None or np.linalg.norm(res.fun) < np.linalg.norm(best.fun):
                best = res
        return best.x, fn, score(fn(best.x, FIT6))

    results = {}
    f = np.linspace(0.045, 0.70, 1700)
    fig, ax = plt.subplots(1, 2, figsize=(7.12, 3.06))
    ax[0].semilogy(f, truth(f), lw=2.3, color="black", label="Two-sea response target")
    for name, label, col in (
        ("single_integer", "1 CARMA oscillator", COL["orange"]),
        ("single_fractional", "1 fractional oscillator", COL["purple"]),
        ("two_integer", "2 rational oscillators", COL["green"])):
        pars, fn, err = optimize(name)
        results[name] = dict(parameters=list(map(float, pars)),
                             weighted_RMS_log_PSD_dB=err)
        ax[0].semilogy(f, fn(pars, f), color=col,
                       label=f"{label} ({err:.1f} dB)")
        ratio = 10 * np.log10(np.maximum(fn(pars, FIT6), 1e-100)
                              / np.maximum(truth(FIT6), 1e-100))
        ax[1].semilogx(FIT6, ratio, color=col, label=label)
    ax[0].set(xlim=(0.052, 0.65), ylim=(1e-4, 5),
              xlabel="Frequency (Hz)", ylabel=r"Heave PSD ($\mathrm{m^2/Hz}$)",
              title="(a) Surrogate heave response, two sea peaks")
    ax[0].grid(True); ax[0].legend(frameon=False, fontsize=6.9, loc="upper right")
    ax[1].axhline(0, color="black", linewidth=0.85)
    ax[1].set(xlim=(0.051, 0.66), xlabel="Frequency (Hz)",
              ylabel="Model / target PSD (dB)",
              title="(b) Model mismatch in the fit band")
    ax[1].grid(True); ax[1].legend(frameon=False, fontsize=6.9)
    fig.subplots_adjust(wspace=0.28)
    save(fig, 6, "mixed-sea")
    summary["mixed_swell_wind"] = dict(
        swell=dict(type="PM", peak_Hz=0.095, Hs_m=1.25),
        wind=dict(type="JONSWAP", peak_Hz=0.29, gamma=3.3, Hs_m=1.0),
        surrogate_heave_RAO_squared="1/(1+(f/0.70Hz)^4), unvalidated illustration",
        target_variance_m2=float(variance), fitting_band_Hz=[float(FIT6[0]), float(FIT6[-1])],
        objective="weighted RMS log PSD residual, fixed common integrated response variance",
        fits=results)


def math_checks():
    f = np.geomspace(0.001, 10, 250)
    assert np.allclose(frac(f, 2., .29, .22), c41(f, .29, .22, .29, .22),
                       atol=1e-12, rtol=1e-12)
    nu = 1.85
    low = np.geomspace(1e-8, 3e-7, 80)
    high = np.geomspace(1e4, 1e5, 80)
    for band, expected in ((low, 2.), (high, 2.-4*nu)):
        slope = np.polyfit(np.log(band), np.log(positive_shape(frac(band, nu,
                                                                   .29, .22))), 1)[0]
        assert abs(slope - expected) < 0.003, (slope, expected)
    # Analytical implication: f^(2-4*nu) has finite k-th derivative energy
    # only if 2*k+2-4*nu < -1, equivalent to nu>(2*k+3)/4.
    assert [(2*k+3)/4 for k in range(3)] == [0.75, 1.25, 1.75]


def main():
    math_checks()
    summary = dict(
        generated_by=Path(__file__).name,
        disclaimer="Numerical spectral illustrations, not INS or hardware measurements",
        convention="one-sided S(f) in m^2/Hz; variance = integral_0^inf S(f)df",
        variance_quadrature_band_Hz=[float(GRID[0]), float(GRID[-1])],
        fit_weights="0.2+0.8*(S_target/max(S_target))^0.3",
        spectral_fitting_objective="weighted RMS of 10log10(Smodel/Starget); model fixed to target integrated variance",
    )
    fits = fig1(summary)
    fig2(summary)
    fig3(summary)
    fig4(fits, summary)
    fig5(summary)
    fig6(summary)
    with (OUT / "fractional-marine-motion-fit-summary.json").open("w") as file:
        json.dump(summary, file, indent=2, sort_keys=True)
        file.write("\n")
    print("Generated six vector PDFs and fractional-marine-motion-fit-summary.json")
    print("Fractional JONSWAP fits:", summary["incident_JONSWAP_fractional_fits"])
    print("Rational order accuracy:", {k: v["log_PSD_RMS_dB"] for k, v in
          summary["rational_approximation"]["model_orders"].items()})
    print("Mixed seas:", {k: v["weighted_RMS_log_PSD_dB"] for k, v in
          summary["mixed_swell_wind"]["fits"].items()})


if __name__ == "__main__":
    main()
