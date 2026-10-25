#!/usr/bin/env python3
"""Steady-state error decomposition of the vertical heave filters.

For a stationary true vessel-heave spectrum S_h(f) (incident spectrum times
the benchmark truth's heave RAO, CarmaBenchTruth.h), accelerometer white noise of variance sigma_n^2 at
fs = 200 Hz and a filter that assumes measurement variance R, the error of a
linear time-invariant estimator h_hat = T(f) y_a, y_a = a + n, a = -w^2 h, is

    E(f) = |1 + w^2 T(f)|^2 S_h(f) + |T(f)|^2 N(f),     N = 2 sigma_n^2 / fs,

and Var(h - h_hat) = int_0^inf E(f) df.  T(f) is the exact causal transfer of
each filter's steady-state discrete Kalman filter (scipy DARE), so this is a
mismatched-Wiener analysis: it separates what each prior does with the true
spectrum from adaptation, bias and attitude effects.  Two parts are reported:

    signal distortion  int |1 + w^2 T|^2 S_h      (model/spectral mismatch)
    noise propagation  int |T|^2 N                (measurement observability)

OU-III's vertical chain [v, p, S, a_w] carries its S = 0 pseudo-measurement as
a second output with the information-rate-equivalent per-sample variance
r_S^2 T_S / dt (SpectralMSE r_S, T_S = c_T tau, the shipping constants).
CARMA filters are the library models in [h, v, a, latent] coordinates, built
here from the polynomial factors in their companion form; the transfer
function from y to h_hat is coordinate free.  "matched C52" is CARMA(5,2)
refitted (log-spectral, constrained as in tools/carma_fit.py) to the true
vessel spectrum of each case: the best this model class can do with a perfect
shape.

The bias state is omitted (it acts below 1e-3 Hz); the operating point is the
oracle one (tau = T_z/2, sigma_aw = 0.9 sigma_a, CARMA kappa from T_z and
sigma_a from Var a), as in the benchmark's --params oracle runs.

Run: python3 tools/carma_mismatch_analysis.py [--out reports/results/carma_benchmark]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import carma_process_validation as pv

TWO_PI = 2.0 * np.pi
FS = 200.0
DT = 1.0 / FS
SIGMA_N = 1.51e-3 * 9.80665               # injected accelerometer white noise
R_FILTER = (2.8 * SIGMA_N * 0.71) ** 2      # filter measurement variance
# OU-III constants (SeaStateFusionFilter_OU_III.h / SeaStateFusionDefaults.h)
C_J = 0.0538
R_A = 0.0148 ** 2 * 0.005
TS_RATIO = 0.015 / 1.1
TS_MIN, TS_MAX = 0.005, 0.15
RS_MIN, RS_MAX = 0.15, 100.0


def rao2(f, fn=0.55, z=0.35, avg_m=0.5):
    """|heave RAO|^2 of the benchmark truth (CarmaBenchTruth.h vessel_rao)."""
    wn, w = TWO_PI * fn, TWO_PI * f
    k = w ** 2 / 9.80665
    return np.abs(wn ** 2 * np.exp(-k * avg_m) / (wn ** 2 - w ** 2 + 2j * z * wn * w)) ** 2


CASES = {
    "jonswap_fp0.29_hs0.8_g3.3": [(0.8, 0.29, 3.3)],
    "pm_fp0.25_hs1": [(1.0, 0.25, 1.0)],
    "jonswap_fp0.20_hs1_g3.3": [(1.0, 0.20, 3.3)],
    "jonswap_fp0.45_hs0.4_g3.3": [(0.4, 0.45, 3.3)],
    "mixed_swell_wind": [(1.0, 0.09, 5.0), (0.8, 0.30, 3.3)],
    "bimodal": [(1.2, 0.13, 3.3), (0.4, 0.42, 3.3)],
}


def true_spectrum(case, f):
    """One-sided vessel heave PSD per Hz, scaled to the case's Hs.  Each part
    is band-limited to 0.55-4.5 f_p as in the benchmark's synthesis."""
    S = np.zeros_like(f)
    df = f[1] - f[0]
    for hs, fp, g in CASES[case]:
        s = np.where((f >= 0.55 * fp) & (f <= 4.5 * fp), pv.jonswap(f, fp, g), 0.0)
        S += s * (hs / 4) ** 2 / (s.sum() * df)
    return S * rao2(f)


def expm(M):
    from scipy.linalg import expm as _e
    return _e(M)


def van_loan(A, W, h):
    n = A.shape[0]
    M = np.zeros((2 * n, 2 * n))
    M[:n, :n] = -A * h
    M[:n, n:] = W * h
    M[n:, n:] = A.T * h
    E = expm(M)
    F = E[n:, n:].T
    Q = F @ E[:n, n:]
    return F, 0.5 * (Q + Q.T)


def carma_system(order, shape, kappa, accel_std):
    """Companion realization; returns (A, W, C_h, C_a)."""
    B, P = pv.polys(order, shape, kappa)
    n = len(P) - 1
    a = P[1:] / P[0]
    A = np.zeros((n, n)); A[:-1, 1:] = np.eye(n - 1); A[-1, :] = -a[::-1]
    G = np.zeros((n, 1)); G[-1, 0] = 1.0
    b = np.zeros(n); b[:len(B)] = B[::-1]
    Ch = b.reshape(1, n)
    from scipy.linalg import solve_continuous_lyapunov
    P1 = solve_continuous_lyapunov(A, -(G @ G.T))
    Ca = Ch @ A @ A
    g2 = accel_std ** 2 / float((Ca @ P1 @ Ca.T)[0, 0])
    return A, g2 * (G @ G.T), Ch, Ca


def ou3_system(tau, sigma_aw):
    """[v, p, S, a] integrated OU; outputs p (heave), a (accel), S (pseudo)."""
    A = np.zeros((4, 4))
    A[0, 3] = 1.0; A[1, 0] = 1.0; A[2, 1] = 1.0; A[3, 3] = -1.0 / tau
    W = np.zeros((4, 4)); W[3, 3] = 2.0 * sigma_aw ** 2 / tau
    Ch = np.array([[0, 1.0, 0, 0]])
    Ca = np.array([[0, 0, 0, 1.0]])
    Cs = np.array([[0, 0, 1.0, 0]])
    return A, W, Ch, Ca, Cs


def ou3_rs(tau, sigma_aw):
    TS = min(max(TS_RATIO * tau, TS_MIN), TS_MAX)
    sigma_aB = max(sigma_aw / 0.9, 1e-6)
    rs = C_J * (2 * R_A) ** (1 / 14) * (sigma_aB * tau ** 4) ** (6 / 7) / np.sqrt(TS)
    return min(max(rs, RS_MIN), RS_MAX), TS


def steady_transfer(F, Q, H, R, Ch, f, meas_index=0):
    """Transfer from measurement component meas_index to h of the steady-state
    measurement-updated Kalman filter x_k = (I - K H) F x_{k-1} + K y_k."""
    from scipy.linalg import solve_discrete_are
    P = solve_discrete_are(F.T, H.T, Q, R)          # prior covariance
    S = H @ P @ H.T + R
    K = P @ H.T @ np.linalg.inv(S)
    n = F.shape[0]
    Acl = (np.eye(n) - K @ H) @ F
    z = np.exp(1j * TWO_PI * f * DT)
    T = np.empty(len(f), dtype=complex)
    k = K[:, meas_index]
    for i, zi in enumerate(z):
        # X(z) = Acl z^-1 X + K Y  ->  X = (I - Acl z^-1)^-1 K Y
        x = np.linalg.solve(np.eye(n) - Acl / zi, k)
        T[i] = (Ch @ x)[0]
    return T, P


def error_budget(T, S_h, f):
    w = TWO_PI * f
    df = f[1] - f[0]
    N = 2.0 * SIGMA_N ** 2 / FS
    sig = np.abs(1.0 + w ** 2 * T) ** 2 * S_h
    noi = np.abs(T) ** 2 * N
    return sig, noi, float(np.sqrt((sig.sum() + noi.sum()) * df)), float(np.sqrt(sig.sum() * df)), float(np.sqrt(noi.sum() * df))


def matched_c52(S_h, f):
    """Constrained log-spectral CARMA(5,2) fit to S_h (shape only)."""
    from scipy.optimize import least_squares
    import carma_fit as cf
    pk = S_h.max()
    m = (S_h > 1e-3 * pk) & (f > 0.02) & (f < 2.0)
    fg = f[m]
    target = 10 * np.log10(S_h[m])
    names = cf.NAMES["52"]
    fp = f[S_h.argmax()]
    lo = np.array([0.3 * fp, 0.05, 0.3 * fp, 0.05, 0.05, 0.03, 0.20])
    hi = np.array([3.0 * fp, 1.0, 10 * fp, 1.0, 10.0, 3.0 * fp, 1.0])
    def resid(x):
        d = 10 * np.log10(cf.carma_psd(fg, "52", dict(zip(names, x)))) - target
        return d - d.mean()
    best = None
    rng = np.random.default_rng(7)
    for _ in range(24):
        x0 = lo + (hi - lo) * (0.1 + 0.8 * rng.random(7))
        r = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", max_nfev=3000)
        if best is None or r.cost < best.cost:
            best = r
    p = dict(zip(names, best.x))
    return {"f1": p["f1"], "z1": p["z1"], "f2": p["f2"], "z2": p["z2"], "fl": p["fl"], "fz": p["fz"], "zz": p["zz"]}, \
        float(np.sqrt(np.mean(resid(best.x) ** 2)))


def analyse(case, shapes):
    f = np.linspace(1e-4, 20.0, 200000)
    df = f[1] - f[0]
    S_h = true_spectrum(case, f)
    w = TWO_PI * f
    var_h = S_h.sum() * df
    var_v = (w ** 2 * S_h).sum() * df
    var_a = (w ** 4 * S_h).sum() * df
    fz = np.sqrt(var_v / var_h) / TWO_PI
    tau = 0.5 / fz
    sigma_a = np.sqrt(var_a)
    fe = np.linspace(1e-3, 5.0, 6000)           # evaluation grid
    S_e = np.interp(fe, f, S_h)
    out = {"case": case, "hs_heave": 4 * np.sqrt(var_h), "fz_hz": fz, "sigma_a": sigma_a, "filters": {}}
    bands = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.8), (0.8, 5.0)]

    def record(name, T):
        sig, noi, tot, s, n = error_budget(T, S_e, fe)
        e = sig + noi
        dfe = fe[1] - fe[0]
        out["filters"][name] = {
            "rms_m": tot, "signal_distortion_m": s, "noise_m": n,
            "rms_pct_hs": 100 * tot / out["hs_heave"],
            "band_rms_m": {f"{a}-{b}": float(np.sqrt(e[(fe >= a) & (fe < b)].sum() * dfe)) for a, b in bands},
        }
        out.setdefault("_curves", {})[name] = e

    # OU-III vertical chain with S pseudo-measurement at information-rate equivalence.
    sigma_aw = 0.9 * sigma_a
    A, W, Ch, Ca, Cs = ou3_system(tau, sigma_aw)
    F, Q = van_loan(A, W, DT)
    rs, TS = ou3_rs(tau, sigma_aw)
    H = np.vstack([Ca, Cs])
    R = np.diag([R_FILTER, rs ** 2 * TS / DT])
    T, _ = steady_transfer(F, Q, H, R, Ch, fe)
    record("OU3", T)

    for name, order, shape in (("C30", 30, shapes[30]), ("C41", 41, shapes[41]), ("C52", 52, shapes[52])):
        _, P1 = pv.polys(order, shape)          # base zero-crossing frequency
        (vh, vv, va), _ = pv.companion_variances(order, shape)
        fz_base = np.sqrt(vv / vh) / TWO_PI
        A, W, Ch, Ca = carma_system(order, shape, fz / fz_base, sigma_a)
        F, Q = van_loan(A, W, DT)
        T, _ = steady_transfer(F, Q, Ca, np.array([[R_FILTER]]), Ch, fe)
        record(name, T)

    mshape, mrms = matched_c52(S_h, f)
    A, W, Ch, Ca = carma_system(52, mshape, 1.0, sigma_a)
    F, Q = van_loan(A, W, DT)
    T, _ = steady_transfer(F, Q, Ca, np.array([[R_FILTER]]), Ch, fe)
    record("C52_matched", T)
    out["matched_c52_shape"] = mshape
    out["matched_c52_fit_rms_db"] = mrms
    out["_fe"] = fe
    out["_S"] = S_e
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1]
                    / "reports" / "results" / "carma_benchmark")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args(argv)
    shapes = pv.header_shapes()
    results = [analyse(c, shapes) for c in CASES]
    clean = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "mismatch_analysis.json").write_text(json.dumps(clean, indent=2, sort_keys=True) + "\n")
    for r in clean:
        print(f"{r['case']:28s} Hs_heave {r['hs_heave']:.2f} m  f_z {r['fz_hz']:.3f} Hz")
        for n, v in r["filters"].items():
            print(f"   {n:12s} RMS {v['rms_m']:.4f} m ({v['rms_pct_hs']:.2f} %Hs)  signal {v['signal_distortion_m']:.4f}"
                  f"  noise {v['noise_m']:.4f}  bands " + " ".join(f"{k}:{x:.4f}" for k, x in v["band_rms_m"].items()))
    if not args.no_plots:
        plot(results, args.out)


def plot(results, out):
    import importlib.util
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    spec = importlib.util.spec_from_file_location("svg_determinism", pv.ROOT / "plots" / "svg_determinism.py")
    svgd = importlib.util.module_from_spec(spec); spec.loader.exec_module(svgd)
    svgd.configure_svg(matplotlib)
    colors = {"OU3": "#1f77b4", "C30": "#2ca02c", "C41": "#ff7f0e", "C52": "#d62728", "C52_matched": "#9467bd"}
    picks = [r for r in results if r["case"] in ("jonswap_fp0.29_hs0.8_g3.3", "pm_fp0.25_hs1", "mixed_swell_wind")]
    fig, axes = plt.subplots(1, len(picks), figsize=(12.0, 3.6))
    for ax, r in zip(axes, picks):
        fe = r["_fe"]
        m = fe < 1.2
        ax.semilogy(fe[m], r["_S"][m], color="#000000", lw=1.2, label="true heave PSD")
        for n, e in r["_curves"].items():
            ax.semilogy(fe[m], e[m], color=colors[n], lw=0.9, label=f"{n.replace('OU3', 'OU-III')} error")
        ax.set_ylim(1e-7, 1e2)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_title(r["case"], fontsize=8)
        ax.grid(True, which="both", color="#eeeeee", lw=0.5)
    axes[0].set_ylabel("PSD (m$^2$/Hz)")
    axes[0].legend(fontsize=6, loc="upper right")
    fig.suptitle("Steady-state error PSD of each prior against the true vessel heave spectrum (oracle operating point)",
                 fontsize=9)
    fig.tight_layout()
    svgd.save_svg(fig, out / "mismatch_error_psd.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
