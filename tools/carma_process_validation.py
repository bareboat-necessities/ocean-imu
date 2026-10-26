#!/usr/bin/env python3
"""Independent validation of the experimental CARMA vertical priors.

Reads the long realizations written by tests/kalman_carma/carma_process-sim
(the C++ library model at kappa = 1 and unit stationary acceleration std) and
checks them against quantities computed here independently of the C++ code:

* the spectral shapes are parsed from src/kalman_carma/CarmaVerticalModel.h
  and must match reports/results/carma_model/carma_fit.json;
* stability: every root of P(s) strictly in the left half plane;
* the driving gain for unit acceleration std and the stationary variances
  come from a companion-form realization built here from the polynomial
  coefficients and scipy's continuous Lyapunov solver;
* Welch PSDs of the simulated h, v, a against g^2 |B|^2/|P|^2, w^2 of it and
  w^4 of it (sampling aliasing included), plus the empirical ratios
  S_v/S_h = w^2 and S_a/S_h = w^4, which need no normalization;
* continuity and differentiability on the 1 ms realization: h' = v and v' = a
  to discretization order, while a has Brownian increments
  (E[da^2] = g^2 dt), i.e. a is continuous but not differentiable;
* bounded sample paths and finite sample variances over 11.6 h.

With --out it also writes process_validation.json and the PSD figures
(JONSWAP, Pierson-Moskowitz, CARMA(3,0), (4,1), (5,2) on identical axes,
each normalized to its own peak).  The fractional oscillator mentioned in the
study request is not in the repository, so it is not plotted.

Run: python3 tools/carma_process_validation.py --input DIR [--check] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "src" / "kalman_carma" / "CarmaVerticalModel.h"
FIT = ROOT / "reports" / "results" / "carma_model" / "carma_fit.json"
TWO_PI = 2.0 * np.pi

FIELDS = ["f1", "z1", "f2", "z2", "fl", "fz", "zz"]


def header_shapes():
    text = HEADER.read_text()
    out = {}
    for n in (30, 41, 52):
        m = re.search(r"CARMA%d_JONSWAP_SHAPE\{([^}]*)\}" % n, text)
        if not m:
            raise SystemExit(f"CARMA{n} shape not found in {HEADER}")
        vals = [float(v) for v in m.group(1).replace("\n", " ").split(",")]
        out[n] = dict(zip(FIELDS, vals))
    return out


def polys(order, p, kappa=1.0):
    """Monic numerator and denominator coefficients (highest power first)."""
    q = lambda f, z: [1.0, 2 * z * TWO_PI * kappa * f, (TWO_PI * kappa * f) ** 2]
    if order == 30:
        return np.array([1.0]), np.polymul(q(p["f1"], p["z1"]), [1.0, TWO_PI * kappa * p["fl"]])
    if order == 41:
        return np.array([1.0, TWO_PI * kappa * p["fz"]]), np.polymul(q(p["f1"], p["z1"]), q(p["f2"], p["z2"]))
    return (np.array(q(p["fz"], p["zz"])),
            np.polymul(np.polymul(q(p["f1"], p["z1"]), q(p["f2"], p["z2"])), [1.0, TWO_PI * kappa * p["fl"]]))


def unit_psd(order, p, f):
    """|B(i w)|^2 / |P(i w)|^2, two-sided, per rad/s convention."""
    B, P = polys(order, p)
    s = 1j * TWO_PI * np.asarray(f, dtype=float)
    return np.abs(np.polyval(B, s)) ** 2 / np.abs(np.polyval(P, s)) ** 2


def companion_variances(order, p):
    """Stationary Var h, v, a for unit noise from a companion realization."""
    from scipy.linalg import solve_continuous_lyapunov

    B, P = polys(order, p)
    n = len(P) - 1
    a = P[1:] / P[0]
    A = np.zeros((n, n))
    A[:-1, 1:] = np.eye(n - 1)
    A[-1, :] = -a[::-1]
    G = np.zeros((n, 1)); G[-1, 0] = 1.0
    b = np.zeros(n); b[:len(B)] = B[::-1]          # b0 + b1 s + ...
    C = b.reshape(1, n)
    Pinf = solve_continuous_lyapunov(A, -(G @ G.T))
    rows = [C, C @ A, C @ A @ A]
    return [float((r @ Pinf @ r.T)[0, 0]) for r in rows], float(np.max(np.real(np.linalg.eigvals(A))))


def aliased_one_sided(order, p, f, fs, gain2, moment, terms=400):
    """One-sided PSD per Hz of the process sampled at fs: 2 sum_k S(f + k fs)."""
    tot = np.zeros_like(f)
    for k in range(-terms, terms + 1):
        fk = f + k * fs
        tot += (TWO_PI * fk) ** (2 * moment) * unit_psd(order, p, fk)
    return 2.0 * gain2 * tot


def jonswap(f, fp=0.29, gamma=3.3):
    f = np.asarray(f, dtype=float)
    out = np.zeros_like(f)
    m = f > 0
    sig = np.where(f[m] <= fp, 0.07, 0.09)
    r = np.exp(-((f[m] - fp) ** 2) / (2 * sig ** 2 * fp ** 2))
    out[m] = f[m] ** -5 * np.exp(-1.25 * (fp / f[m]) ** 4) * gamma ** r
    return out


def band_db_error(f, est, ref, lo, hi, nbands=24):
    """Mean dB difference est/ref in log-spaced bands over [lo, hi]."""
    edges = np.geomspace(lo, hi, nbands + 1)
    errs = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (f >= a) & (f < b)
        if m.sum() < 2:
            continue
        errs.append(10 * np.log10(np.mean(est[m]) / np.mean(ref[m])))
    errs = np.array(errs)
    return float(np.mean(np.abs(errs))), float(np.max(np.abs(errs)))


def validate(input_dir: Path, out_dir: Path | None, check: bool):
    from scipy.signal import welch

    shapes = header_shapes()
    fit = json.loads(FIT.read_text())["models"]
    fit_key = {30: "carma30", 41: "carma41", 52: "carma52_preliminary"}
    names = {"30": ["f1", "z1", "fl"], "41": ["f1", "z1", "f2", "z2", "fz"],
             "52": ["f1", "z1", "f2", "z2", "fl", "fz", "zz"]}
    json_name = {"fz": "fz0"}
    results, failures = {}, []

    for order in (30, 41, 52):
        sh = shapes[order]
        fp = fit[fit_key[order]]["params"]
        for k in names[str(order)]:
            ref = fp[json_name.get(k, k) if order == 41 else k]
            if abs(sh[k] - ref) > 5e-6:
                failures.append(f"CARMA{order}: header {k}={sh[k]} differs from carma_fit.json {ref}")
        (var_h1, var_v1, var_a1), max_re = companion_variances(order, sh)
        if not max_re < 0:
            failures.append(f"CARMA{order}: pole with non-negative real part {max_re}")
        gain2 = 1.0 / var_a1                       # unit stationary acceleration std
        var_ref = np.array([var_h1, var_v1, var_a1]) * gain2

        n = {30: 3, 41: 4, 52: 5}[order]
        coarse = np.fromfile(input_dir / f"carma{n}_coarse.f64").reshape(-1, 3)
        fine = np.fromfile(input_dir / f"carma{n}_fine.f64").reshape(-1, 3)
        fs = 50.0
        r = {"order": f"CARMA({str(order)[0]},{str(order)[1]})"}
        r["finite"] = bool(np.all(np.isfinite(coarse)) and np.all(np.isfinite(fine)))
        r["sample_var_ratio_h_v_a"] = (coarse.var(axis=0) / var_ref).tolist()
        r["max_abs_h_over_sigma"] = float(np.max(np.abs(coarse[:, 0])) / np.sqrt(var_ref[0]))
        r["max_real_pole"] = max_re

        f, Ph = welch(coarse[:, 0], fs=fs, nperseg=1 << 14)
        _, Pv = welch(coarse[:, 1], fs=fs, nperseg=1 << 14)
        _, Pa = welch(coarse[:, 2], fs=fs, nperseg=1 << 14)
        Gh = aliased_one_sided(order, sh, f, fs, gain2, 0)
        Gv = aliased_one_sided(order, sh, f, fs, gain2, 1)
        Ga = aliased_one_sided(order, sh, f, fs, gain2, 2)
        # Bands where the displacement spectrum carries energy: within 40 dB
        # of its peak, above the Welch resolution and below 2 Hz.
        pk = Gh.max()
        lo = max(f[1] * 4, float(f[np.argmax(Gh > 1e-4 * pk)]))
        hi = min(2.0, float(f[len(Gh) - 1 - np.argmax((Gh > 1e-4 * pk)[::-1])]))
        r["band_hz"] = [lo, hi]
        for key, est, ref in (("h", Ph, Gh), ("v", Pv, Gv), ("a", Pa, Ga)):
            mean_db, max_db = band_db_error(f, est, ref, lo, hi)
            r[f"psd_{key}_mean_abs_db"], r[f"psd_{key}_max_abs_db"] = mean_db, max_db
        # Normalization-free check of S_v = w^2 S_h and S_a = w^4 S_h.  The
        # sampled process aliases the slowly decaying acceleration tail into
        # every bin, so the ratio is compared where the aliased share of the
        # reference is below 1 %.
        w = TWO_PI * f
        alias_a = 1.0 - 2.0 * gain2 * w ** 4 * unit_psd(order, sh, f) / Ga
        alias_v = 1.0 - 2.0 * gain2 * w ** 2 * unit_psd(order, sh, f) / Gv
        mv = (f >= lo) & (f <= hi) & (alias_v < 0.01)
        ma = (f >= lo) & (f <= hi) & (alias_a < 0.01)
        r["ratio_band_v_hz"] = [float(f[mv].min()), float(f[mv].max())]
        r["ratio_band_a_hz"] = [float(f[ma].min()), float(f[ma].max())]
        r["ratio_v_h_over_w2_db"] = band_db_error(f[mv], Pv[mv] / Ph[mv], w[mv] ** 2,
                                                  f[mv].min(), f[mv].max(), nbands=12)
        r["ratio_a_h_over_w4_db"] = band_db_error(f[ma], Pa[ma] / Ph[ma], w[ma] ** 4,
                                                  f[ma].min(), f[ma].max(), nbands=12)

        dtf = 0.001
        h, v, a = fine[:, 0], fine[:, 1], fine[:, 2]
        r["dh_dt_minus_v_rel"] = float(np.std(np.diff(h) / dtf - 0.5 * (v[1:] + v[:-1])) / v.std())
        # The trapezoid residual of v' = a is the integral of a Brownian
        # acceleration increment: std g sqrt(dt/12).  Matching it says a is
        # continuous but not differentiable; a differentiable a would give O(dt).
        r["dv_dt_minus_a_rel"] = float(np.std(np.diff(v) / dtf - 0.5 * (a[1:] + a[:-1])) / a.std())
        r["dv_dt_minus_a_over_brownian"] = float(
            np.std(np.diff(v) / dtf - 0.5 * (a[1:] + a[:-1])) / np.sqrt(gain2 * dtf / 12.0))
        da2 = np.mean(np.diff(a) ** 2)
        r["a_increment_ratio_to_g2dt"] = float(da2 / (gain2 * dtf))
        da2_2 = np.mean((a[2:] - a[:-2]) ** 2)
        dv2 = np.mean(np.diff(v) ** 2)
        dv2_2 = np.mean((v[2:] - v[:-2]) ** 2)
        r["a_increment_scaling_exponent"] = float(np.log2(da2_2 / da2))     # 1: Brownian
        r["v_increment_scaling_exponent"] = float(np.log2(dv2_2 / dv2))     # 2: differentiable
        results[str(order)] = r
        spectra = dict(f=f, Ph=Ph, Pv=Pv, Pa=Pa, Gh=Gh, Gv=Gv, Ga=Ga)
        results[str(order)]["_spectra"] = spectra

        if check:
            name = r["order"]
            def need(cond, what, name=name):
                if not cond:
                    failures.append(f"{name}: {what}")
            need(r["finite"], "non-finite sample")
            need(r["max_abs_h_over_sigma"] < 7.0, f"|h| reached {r['max_abs_h_over_sigma']:.2f} sigma")
            for k, ratio in zip("hva", r["sample_var_ratio_h_v_a"]):
                need(abs(ratio - 1.0) < 0.08, f"Var {k} ratio {ratio:.3f}")
            for k in "hva":
                need(r[f"psd_{k}_mean_abs_db"] < 0.35, f"PSD {k} mean |dB| {r[f'psd_{k}_mean_abs_db']:.3f}")
                need(r[f"psd_{k}_max_abs_db"] < 1.2, f"PSD {k} max |dB| {r[f'psd_{k}_max_abs_db']:.3f}")
            need(r["ratio_v_h_over_w2_db"][0] < 0.3, f"S_v/S_h vs w^2 {r['ratio_v_h_over_w2_db']}")
            need(r["ratio_a_h_over_w4_db"][0] < 0.3, f"S_a/S_h vs w^4 {r['ratio_a_h_over_w4_db']}")
            need(r["dh_dt_minus_v_rel"] < 1e-3, f"h' = v residual {r['dh_dt_minus_v_rel']:.2e}")
            need(abs(r["dv_dt_minus_a_over_brownian"] - 1.0) < 0.05,
                 f"v' = a residual / g sqrt(dt/12) = {r['dv_dt_minus_a_over_brownian']:.3f}")
            need(abs(r["a_increment_ratio_to_g2dt"] - 1.0) < 0.05,
                 f"E[da^2]/(g^2 dt) = {r['a_increment_ratio_to_g2dt']:.3f}")
            need(abs(r["a_increment_scaling_exponent"] - 1.0) < 0.1,
                 f"a increment exponent {r['a_increment_scaling_exponent']:.3f}")
            need(abs(r["v_increment_scaling_exponent"] - 2.0) < 0.1,
                 f"v increment exponent {r['v_increment_scaling_exponent']:.3f}")

    summary = {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in results.items()}
    for v in summary.values():
        print(f"{v['order']}: var ratios {np.round(v['sample_var_ratio_h_v_a'], 3).tolist()} "
              f"PSD mean|dB| h={v['psd_h_mean_abs_db']:.3f} v={v['psd_v_mean_abs_db']:.3f} "
              f"a={v['psd_a_mean_abs_db']:.3f}  E[da^2]/(g^2dt)={v['a_increment_ratio_to_g2dt']:.3f} "
              f"exp(a)={v['a_increment_scaling_exponent']:.2f} exp(v)={v['v_increment_scaling_exponent']:.2f}")
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "process_validation.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        plot(results, shapes, out_dir)
    if failures:
        for msg in failures:
            print("FAIL:", msg)
        return 1
    print("carma_process_validation: all checks passed" if check else "carma_process_validation: done")
    return 0


def plot(results, shapes, out_dir: Path):
    import importlib.util
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    spec = importlib.util.spec_from_file_location("svg_determinism", ROOT / "plots" / "svg_determinism.py")
    svgd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(svgd)
    svgd.configure_svg(matplotlib)

    # Figure 1: model shapes on identical axes, each normalized to its own peak.
    f = np.linspace(0.005, 1.6, 4000)
    curves = [
        ("JONSWAP $f_p$=0.29 Hz, $\\gamma$=3.3", jonswap(f, 0.29, 3.3), "#1f77b4", "-"),
        ("Pierson-Moskowitz $f_p$=0.29 Hz", jonswap(f, 0.29, 1.0), "#7f7f7f", "--"),
        ("CARMA(3,0)", unit_psd(30, shapes[30], f), "#2ca02c", "-"),
        ("CARMA(4,1)", unit_psd(41, shapes[41], f), "#ff7f0e", "-"),
        ("CARMA(5,2)", unit_psd(52, shapes[52], f), "#d62728", "-"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for label, S, color, ls in curves:
        ax.plot(f, 10 * np.log10(np.maximum(S / S.max(), 1e-12)), color=color, ls=ls, lw=1.6, label=label)
    ax.axvspan(0.17, 1.30, color="#000000", alpha=0.05, lw=0)
    ax.text(0.18, -57, "fitted band 0.17-1.30 Hz", fontsize=8, color="#555555")
    ax.set_xlim(0, 1.6)
    ax.set_ylim(-60, 3)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Displacement PSD (dB re own peak)")
    ax.set_title("Displacement spectra, each normalized to its own peak (shape only)")
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    svgd.save_svg(fig, out_dir / "carma_psd_overlay.svg", bbox_inches="tight")
    plt.close(fig)

    # Figure 2: simulated (Welch) against analytic, h, v and a, CARMA(5,2).
    s = results["52"]["_spectra"]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    m = (s["f"] > 0.01) & (s["f"] < 5.0)
    for key, label, color in (("h", "displacement", "#1f77b4"), ("v", "velocity", "#ff7f0e"),
                              ("a", "acceleration", "#d62728")):
        ax.loglog(s["f"][m], s["P" + key][m], color=color, lw=0.9, alpha=0.8, label=f"{label}: simulated (Welch)")
        ax.loglog(s["f"][m], s["G" + key][m], color="#000000", lw=0.8, ls="--")
    ax.plot([], [], color="#000000", ls="--", lw=0.8, label="analytic $S_h$, $\\omega^2 S_h$, $\\omega^4 S_h$")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("One-sided PSD (units$^2$/Hz)")
    ax.set_title("CARMA(5,2) realization (11.6 h, unit acceleration std) vs analytic")
    ax.set_ylim(1e-8, 1e1)
    ax.grid(True, which="both", color="#eeeeee", lw=0.5)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    svgd.save_svg(fig, out_dir / "carma52_process_psd.svg", bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)
    return validate(args.input, args.out, args.check)


if __name__ == "__main__":
    sys.exit(main())
