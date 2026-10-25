#!/usr/bin/env python3
"""Analytic study of three-axis CARMA marine-motion priors (EXPERIMENTAL).

Companion of docs/carma-three-axis.md.  Nothing here is used by a shipping
estimator.  Every quantity is computed from the formulas below and from the
benchmark truth model (tests/kalman_carma/CarmaBenchTruth.h, horizontal model
"physical"); no Monte Carlo result enters.

1. Physical spectra.  Directional random-phase sea (the benchmark's
   discretization: 5 directions with cos^2s weights, each spectral part
   band-limited to 0.55-4.5 f_p) through the vessel RAOs.  Heave: the
   benchmark's rigid-body heave RAO.  Surge/sway: a free-floating hull follows
   the orbital motion in long waves (no hydrostatic restoring, no resonance),
   with hull averaging over draft, length and beam through effective
   averaging lengths (the e^{-k l} form the truth uses for heave, roll, pitch),
       surge = e^{-k (d + l_L |cos b|)} cos b,
       sway  = e^{-k (d + l_B |sin b|)} sin b,
   in quadrature (-i) with the elevation, b = propagation direction relative
   to the bow.  The cross-spectral matrix S_xyz(f) of [x_N, y_E, z_up] is a
   sum of rank-one terms S_eta w_j t_j t_j^H, Hermitian and positive
   semidefinite by construction; both are checked numerically.
2. CARMA shapes for the horizontal axes: log-spectral shape fits (the
   objective, band and bounds of tools/carma_fit.py) to JONSWAP, vessel heave,
   surge (head sea) and sway (beam sea) targets, including CARMA(5,2) with a
   double zero at DC (B = s^2, the form proposed in stochastic-marine-motion.tex).
3. State-space specification of every candidate in the library coordinates
   [p, v, a, latent] of src/kalman_carma/CarmaVerticalModel.h (Python port of
   carma_dynamics): A, G, output rows, W, exact F and Q at 200 Hz (Van Loan),
   stationary covariance (Lyapunov), stability, controllability,
   observability from acceleration, transfer-function check, conditioning.
4. Auxiliary S = int p dt: spectrum of the augmented system, PBH tests of the
   S mode, prior growth of Var S(t), and the DC-zero criterion.
5. OU-III's S = 0 pseudo-measurement as a DC zero: the two-sided (smoothing)
   posterior of the integrated OU chain conditioned on S = 0 at information
   rate 1/(r_S^2 T_S) is a stationary CARMA(4,1) prior with B(s) = s.  Its
   poles, and the causal steady-state transfers of both.
6. Steady-state (DARE) error budget of every horizontal candidate against the
   physical surge/sway spectra (signal distortion / noise propagation), at
   the operating point the wrapper maps from the vertical channel.
7. Cross-axis information bound: the noncausal Wiener error of each axis
   using all three accelerometer axes with the TRUE cross-spectral matrix,
   against each axis from its own accelerometer only.  This upper-bounds what
   any coupled (vector) CARMA could gain over independent axes.

Run: python3 tools/carma3_analysis.py [--out reports/results/carma3_benchmark] [--no-plots]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import carma_fit as cf  # noqa: E402
import carma_mismatch_analysis as mma  # noqa: E402
import carma_process_validation as pv  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TWO_PI = 2.0 * np.pi
G = 9.80665
FS = 200.0
DT = 1.0 / FS
SIGMA_N = mma.SIGMA_N
R_FILTER = mma.R_FILTER
N_ONE_SIDED = 2.0 * SIGMA_N ** 2 / FS     # accelerometer white noise, one-sided per Hz

# Benchmark truth (CarmaBenchTruth.h).  Heave RAO: rigid-body response.
HEAVE_RAO = (0.55, 0.35, 0.5)
# Physical horizontal RAO: averaging depth and the effective averaging
# lengths of a ~7.5 m waterline and ~2.8 m beam (sinc(kL/2) matched by
# e^{-k l} at kL/2 = 1..pi/2 gives l = 0.09..0.14 L).
HORIZ_DEPTH_M = 0.5
HORIZ_LENGTH_M = 0.75
HORIZ_BEAM_M = 0.28
# OU-III horizontal regularizer factors (SeaStateFusionFilter_OU_III.h).
RS_X_FACTOR, RS_Y_FACTOR = 0.72, 0.50
KAPPA_R = 0.3627                        # OU-III pole constant (S_z pole law)

SEAS = {
    "jonswap_fp0.29_hs0.8_g3.3": [(0.8, 0.29, 3.3, 30.0, 4.0)],
    "pm_fp0.25_hs1": [(1.0, 0.25, 1.0, 30.0, 4.0)],
    "mixed_swell_wind": [(1.0, 0.09, 5.0, 10.0, 12.0), (0.8, 0.30, 3.3, 60.0, 4.0)],
    "bimodal": [(1.2, 0.13, 3.3, 0.0, 8.0), (0.4, 0.42, 3.3, 70.0, 4.0)],
}
FIELDS = ["f1", "z1", "f2", "z2", "fl", "fz", "zz"]


# ---------------------------------------------------------------------------
# 1. Physical spectra
# ---------------------------------------------------------------------------

def vessel_rao(w, fn, zeta, avg_m):
    wn = TWO_PI * fn
    k = w * w / G
    return wn * wn * np.exp(-k * avg_m) / (wn * wn - w * w + 2j * zeta * wn * w)


def horizontal_rao(w, rel):
    """Real surge and sway amplitude ratios (orbital lag applied by caller)."""
    k = w * w / G
    c, s = np.cos(rel), np.sin(rel)
    surge = np.exp(-k * (HORIZ_DEPTH_M + HORIZ_LENGTH_M * np.abs(c))) * c
    sway = np.exp(-k * (HORIZ_DEPTH_M + HORIZ_BEAM_M * np.abs(s))) * s
    return surge, sway


def part_spectrum(f, hs, fp, gamma):
    """One-sided elevation PSD per Hz, band-limited as in the benchmark."""
    s = np.where((f >= 0.55 * fp) & (f <= 4.5 * fp), pv.jonswap(f, fp, gamma), 0.0)
    fl = np.linspace(0.55 * fp, 4.5 * fp, 20001)
    m0 = np.trapezoid(pv.jonswap(fl, fp, gamma), fl)
    return s * (hs / 4.0) ** 2 / m0


def directions(spread_s, nd=5):
    th = -0.5 * np.pi + (np.arange(nd) + 0.5) * np.pi / nd
    w = np.cos(th) ** (2.0 * spread_s)
    return th, w / w.sum()


def cross_spectrum(sea, f, heading_deg=0.0, model="physical"):
    """S[f, 3, 3] (one-sided per Hz) of r = [x_N, y_E, z_up]."""
    w = TWO_PI * f
    S = np.zeros((len(f), 3, 3), dtype=complex)
    Hz = vessel_rao(w, *HEAVE_RAO)
    psi = np.deg2rad(heading_deg)
    for hs, fp, gamma, dir_deg, s in SEAS[sea] if isinstance(sea, str) else sea:
        se = part_spectrum(f, hs, fp, gamma)
        th, wt = directions(s)
        for tj, wj in zip(th, wt):
            d = np.deg2rad(dir_deg) + tj
            if model == "physical":
                su, sw = horizontal_rao(w, d - psi)
                hN = -1j * (su * np.cos(psi) - sw * np.sin(psi))
                hE = -1j * (su * np.sin(psi) + sw * np.cos(psi))
            else:   # legacy benchmark truth: 0.8 x heave RAO x orbital
                hN = -1j * 0.8 * Hz * np.cos(d)
                hE = -1j * 0.8 * Hz * np.sin(d)
            t = np.stack([hN, hE, Hz], axis=1)
            S += (se * wj)[:, None, None] * t[:, :, None] * np.conj(t[:, None, :])
    return S


def coherence(S, i, j):
    den = np.real(S[:, i, i]) * np.real(S[:, j, j])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, np.abs(S[:, i, j]) ** 2 / den, np.nan), np.angle(S[:, i, j], deg=True)


def psd_check(S):
    herm = float(np.max(np.abs(S - np.conj(np.transpose(S, (0, 2, 1))))))
    ev = np.linalg.eigvalsh(0.5 * (S + np.conj(np.transpose(S, (0, 2, 1)))))
    tr = np.real(np.trace(S, axis1=1, axis2=2))
    rel = ev[:, 0] / np.maximum(tr, 1e-300)
    return {"max_hermitian_defect": herm, "min_eig_over_trace": float(np.min(rel[tr > 0]))}


def moments(f, s):
    w = TWO_PI * f
    m0 = np.trapezoid(s, f)
    m2 = np.trapezoid(w ** 2 * s, f)
    m4 = np.trapezoid(w ** 4 * s, f)
    return m0, m2, m4


# ---------------------------------------------------------------------------
# CARMA in library coordinates (port of carma_dynamics)
# ---------------------------------------------------------------------------

def order_n(code):
    return {"30": 3, "41": 4, "41d": 4, "52": 5, "52d": 5}[code]


def carma_A(N, sh, kappa=1.0):
    w1 = TWO_PI * kappa * sh["f1"]
    z1 = sh["z1"]
    A = np.zeros((N, N))
    A[0, 1] = 1.0
    A[1, 2] = 1.0
    if N == 3:
        lam = TWO_PI * kappa * sh["fl"]
        A[2] = [-lam * w1 * w1, -(2 * lam * z1 * w1 + w1 * w1), -(lam + 2 * z1 * w1)]
    elif N == 4:
        w2 = TWO_PI * kappa * sh["f2"]
        z2 = sh["z2"]
        wz0 = TWO_PI * kappa * sh["fz"]
        d = wz0 - 2 * z2 * w2
        A[2] = [d * w1 * w1, 2 * d * z1 * w1 - w1 * w1, d - 2 * z1 * w1, -w2 * w2 - d * wz0]
        A[3] = [w1 * w1, 2 * z1 * w1, 1.0, -wz0]
    else:
        w2 = TWO_PI * kappa * sh["f2"]
        z2 = sh["z2"]
        lam = TWO_PI * kappa * sh["fl"]
        wz = TWO_PI * kappa * sh["fz"]
        zz = sh["zz"]
        b1 = 2 * zz * wz - 2 * z2 * w2
        b0 = wz * wz - w2 * w2
        gm = b1 - lam
        A[2] = [gm * w1 * w1, -w1 * w1 + 2 * gm * z1 * w1, gm - 2 * z1 * w1,
                -b1 * w2 * w2 - gm * b0, b0 - 2 * b1 * z2 * w2 - gm * b1]
        A[3, 4] = 1.0
        A[4] = [w1 * w1, 2 * z1 * w1, 1.0, -w2 * w2 - b0, -2 * z2 * w2 - b1]
    return A


def carma_polys(N, sh, kappa=1.0):
    q = lambda fr, z: [1.0, 2 * z * TWO_PI * kappa * fr, (TWO_PI * kappa * fr) ** 2]
    if N == 3:
        return np.array([1.0]), np.polymul(q(sh["f1"], sh["z1"]), [1.0, TWO_PI * kappa * sh["fl"]])
    if N == 4:
        return (np.array([1.0, TWO_PI * kappa * sh["fz"]]),
                np.polymul(q(sh["f1"], sh["z1"]), q(sh["f2"], sh["z2"])))
    return (np.array(q(sh["fz"], sh["zz"])),
            np.polymul(np.polymul(q(sh["f1"], sh["z1"]), q(sh["f2"], sh["z2"])),
                       [1.0, TWO_PI * kappa * sh["fl"]]))


def unit_psd(N, sh, f, kappa=1.0):
    """|B|^2/|P|^2 (two-sided, rad/s convention, unit noise intensity)."""
    B, P = carma_polys(N, sh, kappa)
    s = 1j * TWO_PI * np.asarray(f, dtype=float)
    return np.abs(np.polyval(B, s)) ** 2 / np.abs(np.polyval(P, s)) ** 2


def lyap(A, W):
    from scipy.linalg import solve_continuous_lyapunov
    P = solve_continuous_lyapunov(A, -W)
    return 0.5 * (P + P.T)


class Carma:
    """One configured CARMA prior: shape, frequency scale, stationary accel std."""

    def __init__(self, N, shape, kappa=1.0, accel_std=1.0):
        self.N, self.shape, self.kappa = N, dict(shape), kappa
        A1 = carma_A(N, shape, 1.0)
        W1 = np.zeros((N, N)); W1[2, 2] = 1.0
        self.base_var = np.diag(lyap(A1, W1))
        self.A = carma_A(N, shape, kappa)
        self.g2 = accel_std ** 2 * kappa / self.base_var[2]
        self.W = np.zeros((N, N)); self.W[2, 2] = self.g2

    @property
    def base_fz_hz(self):
        return np.sqrt(self.base_var[1] / self.base_var[0]) / TWO_PI

    def psd_onesided(self, f):
        return 2.0 * self.g2 * unit_psd(self.N, self.shape, f, self.kappa)

    def c0_two_sided(self):
        return self.g2 * unit_psd(self.N, self.shape, np.array([0.0]), self.kappa)[0]


def library_shapes():
    s = pv.header_shapes()
    return {"C30J": (3, s[30]), "C41J": (4, s[41]), "C52J": (5, s[52])}


# ---------------------------------------------------------------------------
# 2. Shape fits
# ---------------------------------------------------------------------------

FIT_BOUNDS = {"f1": (0.15, 0.60), "z1": (0.05, 1.0), "f2": (0.15, 3.00), "z2": (0.05, 1.0),
              "fl": (0.05, 5.00), "fz": (0.0, 2.0), "zz": (0.20, 1.0)}
FIT_FREE = {"41": ["f1", "z1", "f2", "z2", "fz"], "41d": ["f1", "z1", "f2", "z2"],
            "52": ["f1", "z1", "f2", "z2", "fl", "fz", "zz"], "52d": ["f1", "z1", "f2", "z2", "fl"]}
FIT_FIXED = {"41": {"fl": 0.0, "zz": 0.0}, "41d": {"fz": 0.0, "fl": 0.0, "zz": 0.0},
             "52": {}, "52d": {"fz": 0.0, "zz": 0.2}}


def fit_target(name, f):
    w = TWO_PI * f
    k = w * w / G
    sj = pv.jonswap(f, 0.29, 3.3)
    if name == "jonswap":
        return sj
    if name == "heave":
        return sj * np.abs(vessel_rao(w, *HEAVE_RAO)) ** 2
    if name == "surge_head":
        return sj * np.exp(-2 * k * (HORIZ_DEPTH_M + HORIZ_LENGTH_M))
    if name == "sway_beam":
        return sj * np.exp(-2 * k * (HORIZ_DEPTH_M + HORIZ_BEAM_M))
    raise ValueError(name)


def fit_grid(target_name):
    """The JONSWAP fits use carma_fit's 0.17-1.30 Hz band.  The vessel targets
    fall 30-75 dB across that band (hull averaging), where no relative-degree-3
    rational spectrum can follow them and a log-spectral fit would be decided
    by the deep tail; they are fitted on the energetic band, where the target
    is within 30 dB of its peak."""
    if target_name == "jonswap":
        return cf.GRID
    g = cf.GRID
    t = fit_target(target_name, g)
    return g[t >= 1e-3 * t.max()]


def fit_shape(code, target_name, n_starts=40, seed=20261025):
    from scipy.optimize import least_squares
    grid = fit_grid(target_name)
    target = 10 * np.log10(fit_target(target_name, grid))
    names = FIT_FREE[code]
    lo = np.array([FIT_BOUNDS[n][0] for n in names])
    hi = np.array([FIT_BOUNDS[n][1] for n in names])
    N = order_n(code)

    def shape_of(x):
        sh = dict(FIT_FIXED[code]); sh.update(zip(names, x))
        for k in FIELDS:
            sh.setdefault(k, 0.0)
        return sh

    def resid(x):
        d = 10 * np.log10(unit_psd(N, shape_of(x), grid)) - target
        return d - d.mean()

    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_starts):
        x0 = lo + (hi - lo) * (0.1 + 0.8 * rng.random(len(names)))
        r = least_squares(resid, x0, bounds=(lo, hi), method="trf", x_scale="jac",
                          ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=4000)
        if best is None or r.cost < best.cost:
            best = r
    sh = shape_of(best.x)
    e = resid(best.x)
    return sh, float(np.sqrt(np.mean(e ** 2))), float(np.max(np.abs(e)))


def describe_shape(N, sh):
    m = Carma(N, sh)
    B, P = carma_polys(N, sh)
    f = np.linspace(1e-5, 20, 400001)
    S = unit_psd(N, sh, f)
    s0 = unit_psd(N, sh, np.array([0.0]))[0]
    return {
        "shape": {k: float(sh.get(k, 0.0)) for k in FIELDS},
        "poles": [[float(r.real), float(r.imag)] for r in np.roots(P)],
        "zeros": [[float(r.real), float(r.imag)] for r in np.roots(B)] if len(B) > 1 else [],
        "dc_zero_order": int(len(B) - 1 - np.count_nonzero(np.abs(np.roots(B)) > 1e-9)) if len(B) > 1 else 0,
        "dc_level_db_re_peak": float(10 * np.log10(s0 / S.max())) if s0 > 0 else None,
        "base_zero_crossing_hz": float(m.base_fz_hz),
        "variance_fraction_below_0p1hz": float(np.sum(S[f < 0.1]) / np.sum(S)),
    }


# ---------------------------------------------------------------------------
# 3. State-space specification
# ---------------------------------------------------------------------------

def van_loan(A, W, h):
    return mma.van_loan(A, W, h)


def ctrb_rank(A, G):
    n = A.shape[0]
    M = np.hstack([np.linalg.matrix_power(A, k) @ G for k in range(n)])
    return int(np.linalg.matrix_rank(M, tol=1e-9 * np.linalg.norm(M)))


def obsv_rank(A, C):
    n = A.shape[0]
    M = np.vstack([C @ np.linalg.matrix_power(A, k) for k in range(n)])
    return int(np.linalg.matrix_rank(M, tol=1e-9 * np.linalg.norm(M)))


def spec_block(name, N, sh, fz_target_hz, accel_std):
    base = Carma(N, sh)
    kappa = fz_target_hz / base.base_fz_hz
    m = Carma(N, sh, kappa, accel_std)
    A, W = m.A, m.W
    Gv = np.zeros((N, 1)); Gv[2, 0] = np.sqrt(m.g2)
    Ch = np.eye(N)[[0]]; Cv = np.eye(N)[[1]]; Ca = np.eye(N)[[2]]
    F, Q = van_loan(A, W, DT)
    P = lyap(A, W)
    # Transfer-function check, C (iwI - A)^-1 G against g B/P.
    fchk = np.geomspace(5e-3, 20.0, 200)
    B, Pp = carma_polys(N, sh, kappa)
    err = 0.0
    for fk in fchk:
        s = 1j * TWO_PI * fk
        h = (Ch @ np.linalg.solve(s * np.eye(N) - A, Gv))[0, 0]
        ref = np.sqrt(m.g2) * np.polyval(B, s) / np.polyval(Pp, s)
        err = max(err, abs(h - ref) / abs(ref))
    ev = np.linalg.eigvals(A)
    return {
        "name": name, "states": N, "kappa": kappa, "accel_std": accel_std, "g2": m.g2,
        "A": A.tolist(), "G": Gv[:, 0].tolist(), "C_p": Ch[0].tolist(), "C_v": Cv[0].tolist(),
        "C_a": Ca[0].tolist(), "W_aa": m.g2,
        "CG": float((Ch @ Gv)[0, 0]), "CAG": float((Ch @ A @ Gv)[0, 0]), "CA2G": float((Ch @ A @ A @ Gv)[0, 0]),
        "F_200Hz": F.tolist(), "Q_200Hz": Q.tolist(), "P_stationary": P.tolist(),
        "std_p_v_a": np.sqrt(np.diag(P)[:3]).tolist(),
        "eig_A": [[float(e.real), float(e.imag)] for e in ev], "max_real_eig": float(np.max(ev.real)),
        "controllable_rank": ctrb_rank(A, Gv), "observable_rank_from_a": obsv_rank(A, Ca),
        "observable_rank_from_p": obsv_rank(A, Ch),
        "transfer_rel_err": float(err),
        "cond_P_stationary": float(np.linalg.cond(P)),
        "cond_P_scaled": float(np.linalg.cond(P / np.sqrt(np.outer(np.diag(P), np.diag(P))))),
        "Q_pp_over_Q_aa_200Hz": float(Q[0, 0] / Q[2, 2]),
        "spectral_radius_F": float(np.max(np.abs(np.linalg.eigvals(F)))),
    }


# ---------------------------------------------------------------------------
# 4. Auxiliary S state
# ---------------------------------------------------------------------------

def s_augmented(A, W):
    n = A.shape[0]
    Aa = np.zeros((n + 1, n + 1)); Aa[:n, :n] = A; Aa[n, 0] = 1.0
    Wa = np.zeros((n + 1, n + 1)); Wa[:n, :n] = W
    return Aa, Wa


def s_analysis(name, m):
    n = m.N
    Aa, Wa = s_augmented(m.A, m.W)
    Ga = np.zeros((n + 1, 1)); Ga[2, 0] = np.sqrt(m.g2)
    Ca = np.eye(n + 1)[[2]]
    Cs = np.eye(n + 1)[[n]]
    # PBH at the zero eigenvalue.
    pbh_c = int(np.linalg.matrix_rank(np.hstack([Aa, Ga]), tol=1e-9))
    pbh_o_a = int(np.linalg.matrix_rank(np.vstack([Aa, Ca]), tol=1e-9))
    pbh_o_as = int(np.linalg.matrix_rank(np.vstack([Aa, Ca, Cs]), tol=1e-9))
    # Prior growth of Var S from a stationary start, S(0) = 0.
    P = np.zeros((n + 1, n + 1)); P[:n, :n] = lyap(m.A, m.W)
    # 0.1 s steps: at 1 s the Van Loan exponential of the fast real pole
    # leaks ~1e-7 m^2 s^2 per step into the zero mode, which an
    # uncontrollable integrator then accumulates (itself the Option-B hazard:
    # every Q rounding error in the S mode is integrated without bound).
    F1, Q1 = van_loan(Aa, Wa, 0.1)
    t, vs = [0.0], [0.0]
    for k in range(1, 36001):
        P = F1 @ P @ F1.T + Q1
        if k % 100 == 0:
            t.append(k / 10.0); vs.append(float(P[n, n]))
    c0 = m.c0_two_sided()
    # DC-zero case: Var S(t) -> 2 Var S_inf, with S_inf = (output) of g B(s)/(s P(s)).
    return {
        "name": name, "eig_augmented": [[float(e.real), float(e.imag)] for e in np.linalg.eigvals(Aa)],
        "zero_mode_controllable_from_noise": pbh_c == n + 1,
        "zero_mode_observable_from_accel": pbh_o_a == n + 1,
        "zero_mode_observable_with_S_pseudo": pbh_o_as == n + 1,
        "S_h0_two_sided": float(c0), "var_S_slope_predicted_m2s2_per_s": float(c0),
        "var_S_at_600s": vs[60], "var_S_at_3600s": vs[-1],
        "var_S_slope_600_3600": float((vs[-1] - vs[60]) / 3000.0),
        "_t": t, "_vs": vs,
    }


def s_inf_variance(m):
    """Var of S_inf = int p for a DC-zero model (finite), else inf."""
    if m.c0_two_sided() > 1e-14 * m.g2:
        return float("inf")
    f = np.geomspace(1e-6, 50.0, 400000)
    w = TWO_PI * f
    S = m.psd_onesided(f) / np.maximum(w ** 2, 1e-300)
    return float(np.trapezoid(S, f))


# ---------------------------------------------------------------------------
# 5. OU-III S = 0 as a DC zero
# ---------------------------------------------------------------------------

def ou3_equivalent_carma41(tau, sigma_aw, rs, TS):
    """Spectral factor of the S-conditioned integrated-OU prior.

    S_eff(w) = 2 s^2 tau N w^2 / (2 s^2 tau + N w^6 + N tau^2 w^8), N = rs^2 T_S
             = (2 s^2 / tau) |i w|^2 / |P(i w)|^2,  P monic Hurwitz of degree 4.
    """
    N = rs ** 2 * TS
    D = np.zeros(9)                         # highest power first: s^8 .. s^0
    D[0] = N * tau ** 2
    D[2] = -N
    D[8] = 2.0 * sigma_aw ** 2 * tau
    r = np.roots(D)
    lhp = r[r.real < 0]
    P = np.real(np.poly(lhp))
    g2 = 2.0 * sigma_aw ** 2 / tau
    return P, g2, N


def s_eff_onesided(f, tau, sigma_aw, N):
    w = TWO_PI * f
    S_ou = 2 * sigma_aw ** 2 * tau / (1 + (w * tau) ** 2) / np.maximum(w, 1e-300) ** 4
    return 2.0 * S_ou * N * w ** 2 / (S_ou + N * w ** 2)


def companion(B, P, g2):
    """Companion realization (A, W, C_p, C_a) of g B/P (monic P, highest first)."""
    n = len(P) - 1
    a = P[1:] / P[0]
    A = np.zeros((n, n)); A[:-1, 1:] = np.eye(n - 1); A[-1, :] = -a[::-1]
    W = np.zeros((n, n)); W[-1, -1] = g2
    b = np.zeros(n); b[:len(B)] = np.asarray(B)[::-1]
    Cp = b.reshape(1, n)
    return A, W, Cp, Cp @ A @ A


# ---------------------------------------------------------------------------
# 6. Steady-state error budget
# ---------------------------------------------------------------------------

BANDS = [(0.0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.8), (0.8, 5.0)]


def budget(T, S, f, noise=N_ONE_SIDED):
    w = TWO_PI * f
    sig = np.abs(1.0 + w ** 2 * T) ** 2 * S
    noi = np.abs(T) ** 2 * noise
    e = sig + noi
    return {"rms_m": float(np.sqrt(np.trapezoid(e, f))),
            "signal_m": float(np.sqrt(np.trapezoid(sig, f))),
            "noise_m": float(np.sqrt(np.trapezoid(noi, f))),
            "band_rms_m": {f"{a}-{b}": float(np.sqrt(np.trapezoid(np.where((f >= a) & (f < b), e, 0), f)))
                           for a, b in BANDS}}, e


def carma_transfer(m, f, with_s=False, rs=None, TS=None, eps_s=0.0):
    if not with_s:
        F, Q = van_loan(m.A, m.W, DT)
        Ca = np.eye(m.N)[[2]]
        Cp = np.eye(m.N)[[0]]
        T, _ = mma.steady_transfer(F, Q, Ca, np.array([[R_FILTER]]), Cp, f)
        return T
    Aa, Wa = s_augmented(m.A, m.W)
    Wa[m.N, m.N] += eps_s
    F, Q = van_loan(Aa, Wa, DT)
    n = m.N + 1
    H = np.vstack([np.eye(n)[[2]], np.eye(n)[[m.N]]])
    R = np.diag([R_FILTER, rs ** 2 * TS / DT])
    T, _ = mma.steady_transfer(F, Q, H, R, np.eye(n)[[0]], f)
    return T


def ou3_transfer(tau, sigma_aw, rs, TS, f):
    A, W, Ch, Ca, Cs = mma.ou3_system(tau, sigma_aw)
    F, Q = van_loan(A, W, DT)
    H = np.vstack([Ca, Cs])
    R = np.diag([R_FILTER, rs ** 2 * TS / DT])
    T, _ = mma.steady_transfer(F, Q, H, R, Ch, f)
    return T


def operating_point(S_zz, f):
    """Oracle vertical operating point, as the benchmark's --params oracle."""
    m0, m2, m4 = moments(f, S_zz)
    fz = np.sqrt(m2 / m0) / TWO_PI
    tau = 0.5 / fz
    sigma_a = np.sqrt(m4)
    sigma_aw = 0.9 * sigma_a
    rs, TS = mma.ou3_rs(tau, sigma_aw)
    return {"fz_hz": float(fz), "tau": float(tau), "sigma_a": float(sigma_a), "sigma_aw": float(sigma_aw),
            "rs": float(rs), "TS": float(TS)}


def horizontal_budget(sea, heading_deg, shapes, f):
    S = cross_spectrum(sea, f, heading_deg)
    Sx, Sy, Sz = (np.real(S[:, i, i]) for i in range(3))
    op = operating_point(Sz, f)
    out = {"sea": sea, "heading_deg": heading_deg, "op": op, "axes": {}}
    curves = {}
    for ax, Sa, kfac in (("x", Sx, RS_X_FACTOR), ("y", Sy, RS_Y_FACTOR), ("z", Sz, 1.0)):
        m0, m2, m4 = moments(f, Sa)
        res = {"std_m": float(np.sqrt(m0)), "accel_std": float(np.sqrt(m4)),
               "fz_hz": float(np.sqrt(m2 / m0) / TWO_PI) if m0 > 0 else None, "filters": {}}
        T = ou3_transfer(op["tau"], op["sigma_aw"], kfac * op["rs"], op["TS"], f)
        res["filters"]["OU3"], curves[(ax, "OU3")] = budget(T, Sa, f)
        for name, (N, sh) in shapes.items():
            base = Carma(N, sh)
            kappa = op["fz_hz"] / base.base_fz_hz
            # The wrapper's mapping: horizontal amplitude = the vertical band
            # acceleration RMS (OU-III's isotropic S_factor = 1 prior).
            m = Carma(N, sh, kappa, op["sigma_a"])
            T = carma_transfer(m, f)
            res["filters"][name], curves[(ax, name)] = budget(T, Sa, f)
            if name == "C52J" and ax != "z":
                c0 = m.c0_two_sided()
                rsz = np.sqrt(c0) * op["tau"] / (KAPPA_R * np.sqrt(op["TS"]))
                T = carma_transfer(m, f, True, rsz, op["TS"])
                res["filters"]["C52J+S"], curves[(ax, "C52J+S")] = budget(T, Sa, f)
            if name == "C41J" and ax != "z" and m4 > 0:
                mo = Carma(N, sh, kappa, np.sqrt(m4))   # oracle per-axis amplitude
                T = carma_transfer(mo, f)
                res["filters"]["C41J_oracle_amp"], _ = budget(T, Sa, f)
                # OU-III's S = 0 update carried over unchanged onto a DC-zero
                # CARMA (its S mode is uncontrollable: a double-counted prior).
                # Not stabilizable: the S mode is uncontrollable and on the unit
                # circle, so the DARE has no stabilizing solution; evaluate the
                # limit with a vanishing S process noise.
                try:
                    carma_transfer(m, f, True, kfac * op["rs"], op["TS"])
                    res["dare_C41J+S_stabilizing"] = True
                except np.linalg.LinAlgError:
                    res["dare_C41J+S_stabilizing"] = False
                T = carma_transfer(m, f, True, kfac * op["rs"], op["TS"], eps_s=1e-10)
                res["filters"]["C41J+S(OU3 r_S, eps)"], curves[(ax, "C41J+S")] = budget(T, Sa, f)
        out["axes"][ax] = res
    return out, curves, (Sx, Sy, Sz)


# ---------------------------------------------------------------------------
# 7. Cross-axis information bound
# ---------------------------------------------------------------------------

def wiener_bound(S, f, noise_scale_h=1.0):
    """Noncausal Wiener error per axis: joint (all 3 accelerometer axes, true
    cross-spectrum) against own-axis only.  Horizontal measurement noise is
    noise_scale_h times the accelerometer white noise (tilt-coupled error)."""
    w = TWO_PI * f
    Nd = N_ONE_SIDED * np.array([noise_scale_h ** 2, noise_scale_h ** 2, 1.0])
    e_joint = np.zeros((len(f), 3))
    e_ind = np.zeros((len(f), 3))
    e_xy = np.zeros((len(f), 2))
    for k in range(len(f)):
        Sk = S[k]
        Syy = w[k] ** 4 * Sk + np.diag(Nd)
        E = Sk - w[k] ** 4 * Sk @ np.linalg.solve(Syy, Sk)
        e_joint[k] = np.real(np.diag(E))
        d = np.real(np.diag(Sk))
        e_ind[k] = d * Nd / (w[k] ** 4 * d + Nd)
        S2 = Sk[:2, :2]
        E2 = S2 - w[k] ** 4 * S2 @ np.linalg.solve(w[k] ** 4 * S2 + np.diag(Nd[:2]), S2)
        e_xy[k] = np.real(np.diag(E2))
    rj = np.sqrt(np.trapezoid(e_joint, f, axis=0))
    ri = np.sqrt(np.trapezoid(e_ind, f, axis=0))
    rxy = np.sqrt(np.trapezoid(e_xy, f, axis=0))
    return {"rms_joint_m": rj.tolist(), "rms_independent_m": ri.tolist(), "rms_xy_only_m": rxy.tolist(),
            "gain_joint_pct": (100 * (1 - rj / ri)).tolist(), "gain_xy_only_pct": (100 * (1 - rxy / ri[:2])).tolist()}


def steady_transfer_multi(F, Q, H, R, C, f):
    """Transfer T(f)[k, j] from measurement j to the output C x of the
    steady-state measurement-updated Kalman filter."""
    from scipy.linalg import solve_discrete_are
    P = solve_discrete_are(F.T, H.T, Q, R)
    K = P @ H.T @ np.linalg.inv(H @ P @ H.T + R)
    n = F.shape[0]
    Acl = (np.eye(n) - K @ H) @ F
    T = np.empty((len(f), H.shape[0]), dtype=complex)
    for i, zi in enumerate(np.exp(1j * TWO_PI * f * DT)):
        T[i] = (C @ np.linalg.solve(np.eye(n) - Acl / zi, K))[0]
    return T


def joint_error_psd(T, S, f, row, noise_diag):
    """Error PSD of r_row - T y for y = -w^2 r + n under the TRUE S(f)."""
    w = TWO_PI * f
    e = np.zeros(len(f))
    for k in range(len(f)):
        g = -w[k] ** 2 * T[k]                      # y-weights mapped to r
        g = g.astype(complex)
        g[row] += 0.0
        d = np.zeros(3, dtype=complex); d[row] = 1.0
        h = d - g                                  # e = (d - g) r - T n
        e[k] = np.real(h @ S[k] @ np.conj(h)) + np.sum(np.abs(T[k]) ** 2 * noise_diag)
    return e


def coupled_model_analysis(sea, heading_deg, f, noise_scale_h=1.0):
    """Shared-latent directional CARMA against independent axes (causal,
    steady state, true cross-spectrum).  Vertical C41J state x_z (DC zero, so
    S_inf = c' x_z = int z dt is a state functional); horizontal
        p_N = alpha_N w_z S_inf + r_N,   a_N = alpha_N w_z v_z + a(r_N),
    i.e. the narrow-band quadrature xi = w int eta of linear waves, with
    independent C41J residuals r_N, r_E.  alpha is the least-squares fit of
    the model cross-spectrum to the true one (oracle direction and coherence),
    the residual amplitude the remaining acceleration variance."""
    N, sh = library_shapes()["C41J"]
    S = cross_spectrum(sea, f, heading_deg)
    Sz = np.real(S[:, 2, 2])
    op = operating_point(Sz, f)
    w = TWO_PI * f
    wz = TWO_PI * op["fz_hz"]
    kappa = op["fz_hz"] / Carma(N, sh).base_fz_hz
    mz = Carma(N, sh, kappa, op["sigma_a"])
    c = np.linalg.solve(mz.A.T, np.eye(N)[0])     # c' A = e_h'  =>  (c'x)' = h
    assert abs(c[2]) < 1e-9 * np.max(np.abs(c)), "S_inf needs a DC zero"
    # Coherent coefficients from the true cross-spectrum (z up).
    alpha = []
    for i in (0, 1):
        m = (wz / (1j * np.maximum(w, 1e-9))) * Sz
        alpha.append(float(np.real(np.sum(np.conj(m) * S[:, i, 2])) / np.sum(np.abs(m) ** 2)))
    _, m2z, _ = moments(f, Sz)
    var_vz = m2z
    res_std = []
    for i in (0, 1):
        m4 = moments(f, np.real(S[:, i, i]))[2]
        res_std.append(float(np.sqrt(max(m4 - (alpha[i] * wz) ** 2 * var_vz, 0.05 * m4))))
    infl = []
    for i in (0, 1):
        infl.append(float(np.sqrt(max(op["sigma_a"] ** 2 - (alpha[i] * wz) ** 2 * var_vz, 0.05 * op["sigma_a"] ** 2))))
    out_variants = {}
    for vname, rstd in (("coupled", res_std), ("coupled_inflated", infl)):
        out_variants[vname] = _coupled_eval(mz, N, sh, kappa, rstd, alpha, wz, c, S, f, noise_scale_h)
    out = {"alpha": alpha, "residual_accel_std": res_std, "inflated_residual_accel_std": infl,
           "op": op, "noise_scale_h": noise_scale_h}
    for ax in ("x", "y"):
        out[ax] = {"rms_coupled_m": out_variants["coupled"][ax]["rms"],
                   "rms_coupled_inflated_m": out_variants["coupled_inflated"][ax]["rms"],
                   "rms_independent_m": out_variants["coupled"][ax]["rms_ind"],
                   "wiener_bound_joint_m": wiener_bound(S, f, noise_scale_h)["rms_joint_m"]["xy".index(ax)],
                   "wiener_bound_independent_m": wiener_bound(S, f, noise_scale_h)["rms_independent_m"]["xy".index(ax)]}
    out["model_cross_spectrum_min_eig_over_trace"] = out_variants["coupled"]["min_eig"]
    return out


def _coupled_eval(mz, N, sh, kappa, res_std, alpha, wz, c, S, f, noise_scale_h):
    op_sigma = mz.g2   # unused; kept for symmetry
    del op_sigma
    blocks = [mz] + [Carma(N, sh, kappa, sd) for sd in res_std]
    n = 3 * N
    A = np.zeros((n, n)); W = np.zeros((n, n))
    for b, m in enumerate(blocks):
        A[b * N:(b + 1) * N, b * N:(b + 1) * N] = m.A
        W[b * N:(b + 1) * N, b * N:(b + 1) * N] = m.W
    F, Q = van_loan(A, W, DT)
    H = np.zeros((3, n)); C = np.zeros((3, n))
    for i in (0, 1):
        H[i, (i + 1) * N + 2] = 1.0; H[i, 1] = alpha[i] * wz
        C[i, (i + 1) * N] = 1.0; C[i, :N] = alpha[i] * wz * c
    H[2, 2] = 1.0; C[2, 0] = 1.0
    nd = N_ONE_SIDED * np.array([noise_scale_h ** 2, noise_scale_h ** 2, 1.0])
    R = np.diag(R_FILTER * np.array([noise_scale_h ** 2, noise_scale_h ** 2, 1.0]))
    out = {}
    for i, ax in ((0, "x"), (1, "y")):
        T = steady_transfer_multi(F, Q, H, R, C[[i]], f)
        e_c = joint_error_psd(T, S, f, i, nd)
        # Independent C41J on its own axis, iso amplitude (wrapper mapping).
        mi = Carma(N, sh, kappa, np.sqrt(mz.g2 * mz.base_var[2] / kappa))
        Fi, Qi = van_loan(mi.A, mi.W, DT)
        Ti = mma.steady_transfer(Fi, Qi, np.eye(N)[[2]], np.array([[R[i, i]]]), np.eye(N)[[0]], f)[0]
        Tfull = np.zeros((len(f), 3), dtype=complex); Tfull[:, i] = Ti
        e_i = joint_error_psd(Tfull, S, f, i, nd)
        out[ax] = {"rms": float(np.sqrt(np.trapezoid(e_c, f))), "rms_ind": float(np.sqrt(np.trapezoid(e_i, f)))}
    # The coupled model's own cross-spectral matrix is Hermitian PSD by
    # construction (W >= 0); check it.
    fk = np.geomspace(0.02, 2.0, 60)
    mins = []
    for fr in fk:
        Rr = np.linalg.solve(1j * TWO_PI * fr * np.eye(n) - A, np.eye(n))
        Sm = C @ Rr @ W @ np.conj(Rr).T @ C.T
        mins.append(np.min(np.linalg.eigvalsh(0.5 * (Sm + np.conj(Sm).T))) / np.real(np.trace(Sm)))
    out["min_eig"] = float(np.min(mins))
    return out


# ---------------------------------------------------------------------------
# Cost model
# ---------------------------------------------------------------------------

def cost_model(nx, blocks, ncol_acc):
    """Multiply-adds of one predict + accelerometer update for the structured
    implementation of Kalman3D_Wave_CARMA(3) (dense counts, symmetric P)."""
    base = 6
    nl = nx - base - 3
    apply_fll = sum(b * b for b in blocks)                 # one structured F_LL application
    predict = (2 * base ** 3                                # F_AA P_AA F_AA^T
               + (2 * nl + base + 3) * apply_fll            # P_LL both sides, P_AL, P_LB
               + base * base * nl + base * base * 3 + 9 * nx)
    update = (nx * ncol_acc * 3                             # P C^T
              + ncol_acc * 9                                # S
              + nx * 9 * 2                                  # gain via LDLT solves
              + (nx * (nx + 1) // 2) * (6 + 18))            # Joseph, 3 columns
    return {"states": nx, "P_floats": nx * nx, "P_bytes_float": 4 * nx * nx,
            "predict_madds": int(predict), "acc_update_madds": int(update),
            "total_madds_per_sample": int(predict + update)}


COST_CONFIGS = {
    # name: (NX, linear-block sizes applied by F_LL, accelerometer Jacobian columns)
    "OU3 (shipping, dense 12x12 linear block)": (21, [12], 9),
    "OU3-XY + CARMA(4,1)-Z": (21, [4, 4, 4], 9),
    "OU3-XY + CARMA(5,2)-Z": (22, [4, 4, 5], 9),
    "CARMA(4,1) XYZ": (21, [4, 4, 4], 9),
    "CARMA(5,2) XYZ": (24, [5, 5, 5], 9),
    "CARMA(4,1) XY + CARMA(5,2) Z": (22, [4, 4, 5], 9),
    "CARMA(5,2) XYZ + S_x,S_y": (26, [6, 6, 5], 9),
    "CARMA(4,1) XYZ + S_x,S_y,S_z": (24, [5, 5, 5], 9),
    "coupled CARMA(4,1) XY (shared A, correlated W) + Z": (21, [4, 4, 4], 9),
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(args):
    out = {"conventions": {
        "psd": "one-sided per Hz unless stated; CARMA unit_psd two-sided per rad/s",
        "accelerometer_noise_mps2": SIGMA_N, "filter_R": R_FILTER, "fs_hz": FS,
        "horizontal_rao": {"depth_m": HORIZ_DEPTH_M, "length_m": HORIZ_LENGTH_M, "beam_m": HORIZ_BEAM_M},
        "heave_rao": dict(zip(["fn_hz", "zeta", "avg_m"], HEAVE_RAO))}}
    curves = {}
    f = np.linspace(1e-3, 5.0, 5000)

    # 1. Spectra, coherence, Hermitian PSD.
    spectra = {}
    for sea in SEAS:
        for hd in (0.0, 30.0, 120.0, 210.0) if sea.startswith("jonswap") else (0.0,):
            S = cross_spectrum(sea, f, hd)
            key = f"{sea}@heading{hd:g}"
            g_xy, ph_xy = coherence(S, 0, 1)
            g_xz, ph_xz = coherence(S, 0, 2)
            g_yz, ph_yz = coherence(S, 1, 2)
            Sd = np.real(np.stack([S[:, i, i] for i in range(3)], 1))
            mom = [moments(f, Sd[:, i]) for i in range(3)]
            wgt = np.real(S[:, 2, 2])
            spectra[key] = {
                "psd_check": psd_check(S),
                "std_m": [float(np.sqrt(m[0])) for m in mom],
                "accel_std": [float(np.sqrt(m[2])) for m in mom],
                "zero_crossing_hz": [float(np.sqrt(m[1] / m[0]) / TWO_PI) for m in mom],
                "energy_weighted_coherence_xy_xz_yz": [float(np.nansum(g * wgt) / wgt.sum())
                                                       for g in (g_xy, g_xz, g_yz)],
                "phase_xz_deg_at_heave_peak": float(ph_xz[np.argmax(wgt)]),
            }
            curves[("spec", key)] = (Sd, g_xy, g_xz, g_yz, ph_xz)
    legacy = cross_spectrum("jonswap_fp0.29_hs0.8_g3.3", f, 0.0, model="legacy")
    out["legacy_truth_horizontal_std_m"] = [float(np.sqrt(np.trapezoid(np.real(legacy[:, i, i]), f)))
                                            for i in range(2)]
    curves[("spec", "legacy")] = np.real(np.stack([legacy[:, i, i] for i in range(3)], 1))
    out["spectra"] = spectra

    # 2. Shape fits.
    shapes = library_shapes()
    fits = {}
    for target in ("jonswap", "heave", "surge_head", "sway_beam"):
        for code in ("41", "52", "52d"):
            sh, rms, mx = fit_shape(code, target, n_starts=args.starts)
            d = describe_shape(order_n(code), sh)
            d.update({"rms_db": rms, "max_abs_db": mx, "target": target, "order": code})
            fits[f"{code}:{target}"] = d
    for name, (N, sh) in shapes.items():
        d = describe_shape(N, sh)
        for target in ("jonswap", "heave", "surge_head", "sway_beam"):
            grid = fit_grid(target)
            dd = 10 * np.log10(unit_psd(N, sh, grid)) - 10 * np.log10(fit_target(target, grid))
            dd -= dd.mean()
            d[f"rms_db_vs_{target}"] = float(np.sqrt(np.mean(dd ** 2)))
        fits[name] = d
    out["fits"] = fits
    shapes["C52D"] = (5, fits["52d:jonswap"]["shape"])
    shapes["C41H"] = (4, fits["41:sway_beam"]["shape"])

    # 3. State-space specifications at f_z = 0.37 Hz (JONSWAP 0.29 Hz) and
    # accel std 0.5 m/s^2, plus conditioning across the tuner's range.
    specs = {}
    for name in ("C30J", "C41J", "C52J", "C52D", "C41H"):
        N, sh = shapes[name]
        specs[name] = spec_block(name, N, sh, 0.37, 0.5)
        specs[name]["cond_P_scaled_vs_fz"] = {
            f"{fz:g}": spec_block(name, N, sh, fz, 0.05 if fz < 0.1 else 0.5)["cond_P_scaled"]
            for fz in (0.05, 0.1, 0.2, 0.37, 0.6)}
        specs[name]["cond_P_raw_vs_fz"] = {
            f"{fz:g}": spec_block(name, N, sh, fz, 0.05 if fz < 0.1 else 0.5)["cond_P_stationary"]
            for fz in (0.05, 0.1, 0.2, 0.37, 0.6)}
    out["state_space"] = specs

    # 4. Auxiliary S.
    sres = {}
    for name in ("C41J", "C52J", "C52D"):
        N, sh = shapes[name]
        base = Carma(N, sh)
        m = Carma(N, sh, 0.37 / base.base_fz_hz, 0.5)
        r = s_analysis(name, m)
        r["var_S_inf_stationary"] = s_inf_variance(m)
        curves[("S", name)] = (r.pop("_t"), r.pop("_vs"))
        sres[name] = r
    out["auxiliary_S"] = sres

    # 5. OU-III S = 0 as a DC zero.
    eq = {}
    fz = np.linspace(1e-4, 5, 50000)
    for sea in ("jonswap_fp0.29_hs0.8_g3.3", "pm_fp0.25_hs1"):
        S = cross_spectrum(sea, f, 30.0)
        op = operating_point(np.real(S[:, 2, 2]), f)
        rec = {"op": op}
        for ax, kfac in (("z", 1.0), ("x", RS_X_FACTOR)):
            rs = kfac * op["rs"]
            P, g2, Nn = ou3_equivalent_carma41(op["tau"], op["sigma_aw"], rs, op["TS"])
            roots = np.roots(P)
            fac = [{"f_hz": float(abs(r) / TWO_PI), "zeta": float(-r.real / abs(r))} for r in roots if r.imag >= 0]
            Sfac = 2 * g2 * np.abs(1j * TWO_PI * fz) ** 2 / np.abs(np.polyval(P, 1j * TWO_PI * fz)) ** 2
            Sref = s_eff_onesided(fz, op["tau"], op["sigma_aw"], Nn)
            rel = float(np.max(np.abs(Sfac[1:] - Sref[1:]) / Sref[1:]))
            A, W, Cp, Ca = companion(np.array([1.0, 0.0]), P, g2)
            F, Q = van_loan(A, W, DT)
            T_eq, _ = mma.steady_transfer(F, Q, Ca, np.array([[R_FILTER]]), Cp, f)
            T_ou = ou3_transfer(op["tau"], op["sigma_aw"], rs, op["TS"], f)
            Sa = np.real(S[:, 2, 2]) if ax == "z" else np.real(S[:, 0, 0])
            b_eq, e_eq = budget(T_eq, Sa, f)
            b_ou, e_ou = budget(T_ou, Sa, f)
            # Causal transfer distance in the wave band.
            band = (f > 0.1) & (f < 1.0)
            dist = float(np.max(np.abs(T_eq[band] - T_ou[band]) / np.abs(T_ou[band])))
            rec[ax] = {"rs": rs, "N_S": Nn, "equivalent_carma41_factors": fac, "g2": g2,
                       "factorization_rel_err": rel, "causal_budget_ou3": b_ou,
                       "causal_budget_equivalent_carma41": b_eq,
                       "max_rel_transfer_diff_0p1_1hz": dist}
            curves[("eq", sea, ax)] = (T_ou, T_eq, e_ou, e_eq, Sa)
        eq[sea] = rec
    out["ou3_s_equivalence"] = eq

    # 6. Horizontal steady-state budgets.
    hb = {}
    use = {k: shapes[k] for k in ("C30J", "C41J", "C52J", "C52D", "C41H")}
    for sea in SEAS:
        for hd in ((0.0, 30.0, 120.0, 210.0) if sea.startswith("jonswap") else (0.0,)):
            r, c, Sxyz = horizontal_budget(sea, hd, use, f)
            hb[f"{sea}@heading{hd:g}"] = r
            curves[("hb", sea, hd)] = (c, Sxyz)
    out["horizontal_budget"] = hb

    # 7. Cross-axis bound.
    cb = {}
    for sea in SEAS:
        for hd in ((0.0, 30.0, 120.0) if sea.startswith("jonswap") else (0.0,)):
            S = cross_spectrum(sea, f, hd)
            cb[f"{sea}@heading{hd:g}"] = {f"h_noise_x{c:g}": wiener_bound(S, f, c) for c in (1.0, 3.0, 10.0)}
    out["coupling_bound"] = cb

    cm = {}
    for sea in ("jonswap_fp0.29_hs0.8_g3.3", "pm_fp0.25_hs1", "mixed_swell_wind"):
        for hd in ((30.0, 120.0) if sea.startswith("jonswap") else (0.0,)):
            for c in (1.0, 3.0):
                cm[f"{sea}@heading{hd:g}@hnoise{c:g}"] = coupled_model_analysis(sea, hd, f, c)
    out["coupled_model"] = cm

    out["cost_model"] = {k: cost_model(*v) for k, v in COST_CONFIGS.items()}
    return out, curves


def save(out, curves, outdir, plots):
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "carma3_analysis.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    if plots:
        make_plots(out, curves, outdir)


# Categorical slots of the reference palette.  A shape takes the colour of
# the benchmark filter built from it (tools/carma3_report.py); a derived
# variant is dashed in its family's colour.
COL = {"truth": "#0b0b0b", "OU3": "#2a78d6", "C41J": "#eda100", "C52J": "#e87ba4", "C52J+S": "#6250d6",
       "C52D": "#e87ba4", "C41H": "#52514e", "eq": "#2a78d6", "C30J": "#8a8a86", "C41J+S": "#e34948"}
DASH = {"C52D": "--", "eq": "--"}


def make_plots(out, curves, outdir):
    import importlib.util
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    spec = importlib.util.spec_from_file_location("svg_determinism", ROOT / "plots" / "svg_determinism.py")
    svgd = importlib.util.module_from_spec(spec); spec.loader.exec_module(svgd)
    svgd.configure_svg(matplotlib)
    f = np.linspace(1e-3, 5.0, 5000)
    grid_kw = dict(color="#e6e6e3", lw=0.5)

    def style(ax):
        ax.grid(True, which="both", **grid_kw)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    # Fig 1: axis spectra and fitted shapes.
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.7))
    for ax, hd in zip(axes, (30.0, 120.0, 210.0)):
        Sd = curves[("spec", f"jonswap_fp0.29_hs0.8_g3.3@heading{hd:g}")][0]
        m = f < 1.3
        for i, (lab, ls) in enumerate((("x (north)", "-"), ("y (east)", "--"), ("z heave", ":"))):
            ax.semilogy(f[m], Sd[m, i], color=COL["truth"], ls=ls, lw=1.4, label=lab)
        sh = curves[("spec", "legacy")]
        ax.semilogy(f[m], sh[m, 0], color="#8a8a86", lw=1.0, label="x, legacy truth (0.8 x heave RAO)")
        ax.set_ylim(1e-6, 1)
        ax.set_xlabel("Frequency (Hz)")
        rel = {30.0: "following", 120.0: "beam", 210.0: "head"}[hd]
        ax.set_title(f"JONSWAP 0.29 Hz, 0.8 m, {rel} sea (heading {hd:g} deg, waves to 30 deg)", fontsize=8)
        style(ax)
    axes[0].set_ylabel("Displacement PSD (m$^2$/Hz)")
    axes[0].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    svgd.save_svg(fig, outdir / "axis_spectra.svg", bbox_inches="tight")
    plt.close(fig)

    # Fig 2: shapes vs targets (normalized).
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.7))
    fg = np.linspace(0.02, 1.6, 3000)
    for ax, target in zip(axes, ("heave", "sway_beam")):
        t = fit_target(target, fg)
        ax.semilogy(fg, t / t.max(), color=COL["truth"], lw=1.6, label=f"target: {target.replace('_', ' ')}")
        for name, key in (("C41J", "C41J"), ("C52J", "C52J"), ("C52D", "52d:jonswap"), ("C41H", "41:sway_beam")):
            sh = out["fits"][key]["shape"]
            N = 4 if name.startswith("C41") else 5
            s = unit_psd(N, sh, fg)
            ax.semilogy(fg, s / s.max(), color=COL[name], ls=DASH.get(name, "-"), lw=1.1, label=name)
        ax.axvspan(0.17, 1.30, color="#000000", alpha=0.04, lw=0)
        ax.set_ylim(1e-5, 2)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_title(f"Shapes at kappa = 1 against the {target.replace('_', ' ')} target (peak-normalized)",
                     fontsize=8)
        style(ax)
    axes[0].set_ylabel("Normalized PSD")
    axes[0].legend(fontsize=7, loc="lower center")
    fig.tight_layout()
    svgd.save_svg(fig, outdir / "carma_shapes_horizontal.svg", bbox_inches="tight")
    plt.close(fig)

    # Fig 3: coherence and phase.
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.5))
    for ax, key, title in ((axes[0], "jonswap_fp0.29_hs0.8_g3.3@heading30", "JONSWAP, s = 4, waves on the bow"),
                           (axes[1], "mixed_swell_wind@heading0", "mixed swell (s = 12) + wind sea"),
                           (axes[2], "bimodal@heading0", "bimodal, 70 deg apart")):
        Sd, gxy, gxz, gyz, _ = curves[("spec", key)]
        # Coherence only where the sea has energy (each displacement PSD
        # within 30 dB of its peak).
        m = (f > 0.04) & (f < 1.3) & np.all(Sd > 1e-3 * Sd.max(axis=0), axis=1)
        gxy, gxz, gyz = (np.where(m, g_, np.nan) for g_ in (gxy, gxz, gyz))
        # Axis pairs are not filters: neutral inks, told apart by line style.
        ax.plot(f[m], gxy[m], color="#0b0b0b", lw=1.5, label="x-y")
        ax.plot(f[m], gxz[m], color="#52514e", lw=1.5, ls="--", label="x-z")
        ax.plot(f[m], gyz[m], color="#8a8a86", lw=1.5, ls=":", label="y-z")
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_title(title, fontsize=8)
        style(ax)
    axes[0].set_ylabel("Squared coherence")
    axes[0].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    svgd.save_svg(fig, outdir / "coherence.svg", bbox_inches="tight")
    plt.close(fig)

    # Fig 4: Var S(t) growth.
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    for name in ("C41J", "C52J", "C52D"):
        t, vs = curves[("S", name)]
        ax.plot(np.array(t) / 60, np.sqrt(vs), color=COL[name], ls=DASH.get(name, "-"), lw=1.6, label=name)
    ax.set_xlabel("Time since S = 0 (min)")
    ax.set_ylabel("Prior std of S = int p dt (m s)")
    ax.set_title("Prior spread of the displacement primitive (f_z 0.37 Hz, a 0.5 m/s$^2$)", fontsize=8)
    ax.legend(fontsize=7)
    style(ax)
    fig.tight_layout()
    svgd.save_svg(fig, outdir / "s_variance_growth.svg", bbox_inches="tight")
    plt.close(fig)

    # Fig 5: OU-III / equivalent CARMA(4,1) causal transfer and error PSD.
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    T_ou, T_eq, e_ou, e_eq, Sa = curves[("eq", "jonswap_fp0.29_hs0.8_g3.3", "z")]
    m = f < 1.3
    w = TWO_PI * f
    axes[0].loglog(f[m], np.abs(w[m] ** 2 * T_ou[m]), color=COL["OU3"], lw=1.6, label="OU-III chain + S=0")
    axes[0].loglog(f[m], np.abs(w[m] ** 2 * T_eq[m]), color=COL["eq"], lw=1.2, ls="--",
                   label="equivalent CARMA(4,1), B = s")
    axes[0].set_ylabel("|w^2 T(f)| (acceleration-to-displacement gain x w^2)")
    axes[0].set_xlabel("Frequency (Hz)")
    axes[0].set_title("Causal steady-state transfers, JONSWAP 0.29 Hz operating point", fontsize=8)
    axes[1].semilogy(f[m], Sa[m], color=COL["truth"], lw=1.4, label="true heave PSD")
    axes[1].semilogy(f[m], e_ou[m], color=COL["OU3"], lw=1.4, label="OU-III error")
    axes[1].semilogy(f[m], e_eq[m], color=COL["eq"], lw=1.2, ls="--", label="equivalent CARMA(4,1) error (causal)")
    axes[1].set_ylim(1e-7, 1)
    axes[1].set_xlabel("Frequency (Hz)")
    axes[1].set_ylabel("PSD (m$^2$/Hz)")
    for a in axes:
        style(a)
        a.legend(fontsize=7)
    fig.tight_layout()
    svgd.save_svg(fig, outdir / "ou3_s_equivalence.svg", bbox_inches="tight")
    plt.close(fig)

    # Fig 6: horizontal error PSDs.
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    for ax, (sea, hd) in zip(axes, (("jonswap_fp0.29_hs0.8_g3.3", 30.0), ("pm_fp0.25_hs1", 0.0),
                                    ("mixed_swell_wind", 0.0))):
        c, (Sx, _, _) = curves[("hb", sea, hd)]
        m = f < 1.3
        ax.semilogy(f[m], Sx[m], color=COL["truth"], lw=1.4, label="true x PSD")
        for name in ("OU3", "C41J", "C52J", "C52J+S", "C52D"):
            ax.semilogy(f[m], c[("x", name)][m], color=COL[name], ls=DASH.get(name, "-"), lw=1.1,
                        label=f"{name} error")
        ax.set_ylim(1e-7, 1)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_title(f"x axis, {sea}", fontsize=8)
        style(ax)
    axes[0].set_ylabel("PSD (m$^2$/Hz)")
    axes[0].legend(fontsize=6, loc="upper right")
    fig.tight_layout()
    svgd.save_svg(fig, outdir / "horizontal_error_psd.svg", bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=ROOT / "reports" / "results" / "carma3_benchmark")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--starts", type=int, default=40)
    args = ap.parse_args(argv)
    out, curves = run(args)
    save(out, curves, args.out, not args.no_plots)
    summarize(out)


def summarize(out):
    print("== spectra ==")
    for k, v in out["spectra"].items():
        print(f"{k:40s} std x/y/z {v['std_m'][0]:.3f} {v['std_m'][1]:.3f} {v['std_m'][2]:.3f} m  "
              f"fz {v['zero_crossing_hz'][0]:.3f} {v['zero_crossing_hz'][1]:.3f} {v['zero_crossing_hz'][2]:.3f}  "
              f"coh xy/xz/yz " + " ".join(f"{c:.2f}" for c in v["energy_weighted_coherence_xy_xz_yz"])
              + f"  minEig/tr {v['psd_check']['min_eig_over_trace']:.1e} herm {v['psd_check']['max_hermitian_defect']:.1e}")
    print("legacy truth horizontal std", out["legacy_truth_horizontal_std_m"])
    print("== fits ==")
    for k, v in out["fits"].items():
        extra = " ".join(f"{t}:{v[t]:.2f}" for t in v if t.startswith("rms_db_vs_"))
        print(f"{k:18s} rms {v.get('rms_db', float('nan')):.2f} dB fz0 {v['base_zero_crossing_hz']:.4f} "
              f"DC {v['dc_level_db_re_peak']} {extra} shape " + ",".join(f"{v['shape'][x]:.4f}" for x in FIELDS))
    print("== state space ==")
    for k, v in out["state_space"].items():
        print(f"{k}: maxRe {v['max_real_eig']:.3f} ctrb {v['controllable_rank']} obs(a) {v['observable_rank_from_a']} "
              f"tf err {v['transfer_rel_err']:.1e} CG {v['CG']:.1e} CAG {v['CAG']:.1e} CA2G {v['CA2G']:.3f} "
              f"condPscaled {v['cond_P_scaled']:.1e} raw-vs-fz " + " ".join(f"{a}:{b:.1e}" for a, b in v["cond_P_raw_vs_fz"].items()))
    print("== S ==")
    for k, v in out["auxiliary_S"].items():
        print(f"{k}: ctrl {v['zero_mode_controllable_from_noise']} obs(a) {v['zero_mode_observable_from_accel']} "
              f"obs(a,S) {v['zero_mode_observable_with_S_pseudo']} S_h(0) {v['S_h0_two_sided']:.3e} "
              f"VarS600 {v['var_S_at_600s']:.4g} VarS3600 {v['var_S_at_3600s']:.4g} slope {v['var_S_slope_600_3600']:.3e} "
              f"VarSinf {v['var_S_inf_stationary']:.4g}")
    print("== OU-III S equivalence ==")
    for sea, r in out["ou3_s_equivalence"].items():
        for ax in ("z", "x"):
            v = r[ax]
            print(f"{sea} {ax}: factors " + " ".join(f"{q['f_hz']:.3f}Hz/{q['zeta']:.3f}" for q in v["equivalent_carma41_factors"])
                  + f" fact err {v['factorization_rel_err']:.1e} causal RMS OU3 {v['causal_budget_ou3']['rms_m']:.4f} "
                  f"eqC41 {v['causal_budget_equivalent_carma41']['rms_m']:.4f} Tdiff {v['max_rel_transfer_diff_0p1_1hz']:.2f}")
    print("== horizontal budgets (rms m: signal/noise) ==")
    for k, r in out["horizontal_budget"].items():
        for ax in ("x", "y", "z"):
            a = r["axes"][ax]
            print(f"{k:38s} {ax} std {a['std_m']:.3f} " + " ".join(
                f"{n}:{v['rms_m']:.4f}({v['signal_m']:.4f}/{v['noise_m']:.4f})" for n, v in a["filters"].items()))
    print("== coupling bound (gain % joint vs independent, x/y/z) ==")
    for k, r in out["coupling_bound"].items():
        print(f"{k:38s} " + "  ".join(f"{c}: " + " ".join(f"{g:.1f}" for g in v["gain_joint_pct"])
                                       + " xy:" + " ".join(f"{g:.1f}" for g in v["gain_xy_only_pct"]) for c, v in r.items()))
    print("== coupled shared-latent model (causal) vs independent ==")
    for k, r in out["coupled_model"].items():
        print(f"{k:45s} alpha {r['alpha'][0]:.3f},{r['alpha'][1]:.3f} minEig {r['model_cross_spectrum_min_eig_over_trace']:.1e} " + "  ".join(
            f"{ax}: coupled {r[ax]['rms_coupled_m']:.4f} infl {r[ax]['rms_coupled_inflated_m']:.4f} indep {r[ax]['rms_independent_m']:.4f} "
            f"(Wiener joint {r[ax]['wiener_bound_joint_m']:.4f} / indep {r[ax]['wiener_bound_independent_m']:.4f})" for ax in ("x", "y")))
    print("== cost ==")
    for k, v in out["cost_model"].items():
        print(f"{k:55s} NX {v['states']:3d} P {v['P_bytes_float']:5d} B  madds {v['total_madds_per_sample']}")


if __name__ == "__main__":
    main()
