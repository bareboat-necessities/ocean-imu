#!/usr/bin/env python3
"""Analytic study of coupled three-axis CARMA priors (EXPERIMENTAL research).

Companion of docs/carma-coupled-three-axis.md and
doc/kalman_carma/kalman_carma_coupled.tex.  Nothing here is used by a shipping
estimator.  Every number comes from the formulas below and from the benchmark
truth model (tests/kalman_carma/CarmaBenchTruth.h, horizontal model
"physical"); no Monte Carlo result enters.  The Monte Carlo study is
tools/carma3_coupled_benchmark.py.

Models (state coordinates and noise exactly as in
src/kalman_carma/CarmaCoupledMotion.h):

  IND  baseline: an independent CARMA(4,1) C41J block [p, v, a, l] per axis.
  COR  Experiment A: the same blocks, real 3 x 3 driving covariance Sigma_w,
       Q = Sigma_w (x) Q_1 (axis-major state ordering, so the Kronecker
       factor order is Sigma (x) Q, not Q (x) Sigma).
  SHM  Experiment B: per axis the *complexified* C41J mode (8 real states),
       Hermitian 3 x 3 excitation Sigma_c.  The complexification keeps the
       upper-half-plane poles of P(s) and adds the real spectral factor D(s)
       of |P+(iw)|^2 + |P-(iw)|^2, so that the real part of the complex mode
       has EXACTLY the C41J law (theorem in the report).  SHM with a real
       Sigma_c is COR in law; with a diagonal real Sigma_c it is IND in law.
  MCA  Experiment C: a constrained (RAO-structured) MCARMA: SHM with the
       vertical output passed through a second-order heave RAO, so the
       cross-axis phase depends on frequency.
  SHR  state-reduced shared mode (analysis only): the complexified mode on z,
       x and y = Re(beta zeta_z) + independent C41J residuals (16 states).
  K2   two complexified modes per axis at two oracle frequencies (analysis
       only; multi-system seas).

Sections of the output (reports/results/carma3_coupled/analysis.json):
  1. complexification of C41J: polynomials, exactness, quadrature leakage;
  2. state-space checks of every model (stability, PSD, Kronecker identity,
     stationary covariance against Lyapunov, cross-spectra against the
     realization);
  3. truth cross-spectral matrices and the representable fraction of each
     model family (Experiment A's structural limit, quantified);
  4. oracle model parameters (the demodulated band cross-spectrum the
     deployable estimator converges to, mapped by the Sigma = c S_hat + Delta
     rule of CarmaCrossSpectrumEstimator.h);
  5. steady-state (DARE, 200 Hz) error and consistency of every model against
     the TRUE cross-spectral matrix (white accelerometer noise; horizontal
     noise x1, x3 to stand in for tilt coupling), with the noncausal Wiener
     bound;
  6. identifiability: Whittle Fisher information of each family's
     parameters from T seconds of three-axis acceleration, and of an
     unrestricted MCARMA(4,1);
  7. cost model (state dimension, multiply-adds per sample).

Run: python3 tools/carma3_coupled_analysis.py [--out reports/results/carma3_coupled] [--no-plots]
"""
from __future__ import annotations

import sys as _export_sys
from pathlib import Path as _ExportPath
_export_sys.path.insert(0, str(_ExportPath(__file__).resolve().parents[1] / "plots"))
from vector_export import save_figure


import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import carma3_analysis as ca  # noqa: E402
import carma_mismatch_analysis as mma  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TWO_PI = 2.0 * np.pi
DT = ca.DT
R_FILTER = ca.R_FILTER
N_ONE_SIDED = ca.N_ONE_SIDED
HEAVE_RAO = ca.HEAVE_RAO            # (fn Hz, zeta, averaging m) of the benchmark hull

# Seas of the analysis: the four principal seas of docs/carma-three-axis.md,
# a long-crested sea, two crossing systems at similar frequencies and the
# same JONSWAP sea in beam and head directions.  (hs, fp, gamma, dir deg, s)
SEAS = dict(ca.SEAS)
SEAS.update({
    "jonswap_longcrested": [(0.8, 0.29, 3.3, 30.0, 1.0e4)],
    "jonswap_wide_spread": [(0.8, 0.29, 3.3, 30.0, 1.0)],
    "crossing": [(0.6, 0.25, 3.3, 0.0, 6.0), (0.6, 0.31, 3.3, 90.0, 6.0)],
})
SEA_HEADING = {"jonswap_beam": 0.0}

C41J = {"f1": 0.306945, "z1": 0.104191, "f2": 0.815318, "z2": 0.478188, "fl": 0.0, "fz": 0.0, "zz": 0.0}


# ---------------------------------------------------------------------------
# 1. Complexification of a real CARMA(4,1)
# ---------------------------------------------------------------------------

def c41_poly(sh, kappa=1.0):
    """Monic real B(s) = s + w_z0 and P(s) = L1 L2 (highest power first)."""
    B, P = ca.carma_polys(4, sh, kappa)
    return np.asarray(B, float), np.asarray(P, float)


def complexify(P):
    """Complex monic P_c(s) = P+(s) D(s) for a real monic quartic P with two
    complex-conjugate pole pairs.  P+ holds the poles with Im > 0, D is the
    real Hurwitz quadratic with |D(iw)|^2 = (|P+(iw)|^2 + |P-(iw)|^2) / 2.
    Then |B/P_c|^2(iw) + |B/P_c|^2(-iw) = 2 |B/P|^2(iw) for every real
    w and real B (the real part of the complex mode is the real CARMA)."""
    r = np.roots(P)
    up = r[np.imag(r) > 0]
    if len(up) != 2:
        raise ValueError("complexify needs two complex-conjugate pole pairs")
    Pp = np.poly(up)                      # complex quadratic, monic
    Pm = np.conj(Pp)
    # Q(w) = |Pp(iw)|^2 + |Pm(iw)|^2 as a polynomial in u = w^2:
    # sample and fit (exact: degree 2 in u).
    w = np.array([0.0, 1.0, 2.0])
    q = np.abs(np.polyval(Pp, 1j * w)) ** 2 + np.abs(np.polyval(Pm, 1j * w)) ** 2
    c = np.polyfit(w ** 2, q / 2.0, 2)    # u^2 + beta u + gamma
    beta, gamma = c[1] / c[0], c[2] / c[0]
    d0 = np.sqrt(gamma)
    d1 = np.sqrt(beta + 2.0 * d0)
    D = np.array([1.0, d1, d0])
    return np.polymul(Pp, D), Pp, D


def kinematic_complex_A(B, Pc):
    """Complex 4 x 4 dynamics in [zeta, zeta', zeta'', l] coordinates for
    zeta = B(s)/Pc(s) xi (complex), B = s + w_z0; noise enters zeta'' only.
    Companion x with zeta = c x; coordinates z = T x, T = [c; cA; cA^2; e0]."""
    a = Pc[1:]                              # s^4 + a1 s^3 + a2 s^2 + a3 s + a4
    Ac = np.zeros((4, 4), complex)
    Ac[:3, 1:] = np.eye(3)
    Ac[3, :] = -a[::-1]
    b = np.zeros(4); b[3] = 1.0
    cvec = np.zeros(4, complex); cvec[:len(B)] = B[::-1]
    T = np.vstack([cvec, cvec @ Ac, cvec @ Ac @ Ac, np.eye(4)[0]])
    Az = T @ Ac @ np.linalg.inv(T)
    bz = T @ b
    assert np.allclose(bz[[0, 1, 3]], 0.0, atol=1e-12) and abs(bz[2] - 1.0) < 1e-12
    return Az


def realify(Ac):
    """Real 2n x 2n representation [[Re, -Im], [Im, Re]] (state [Re z; Im z])."""
    return np.block([[Ac.real, -Ac.imag], [Ac.imag, Ac.real]])


class ComplexMode:
    """Complexified C41J at frequency scale kappa, unit noise gain (g = 1).
    Real state [Re(p, v, a, l), Im(p, v, a, l)]; W = e_a e_a' on both halves."""

    def __init__(self, sh=C41J, kappa=1.0):
        B, P = c41_poly(sh, kappa)
        Pc, Pp, D = complexify(P)
        self.B, self.P, self.Pc, self.Pp, self.D = B, P, Pc, Pp, D
        self.Ac = kinematic_complex_A(B, Pc)
        self.A = realify(self.Ac)
        self.W = np.zeros((8, 8)); self.W[2, 2] = 1.0; self.W[6, 6] = 1.0
        self.n = 8

    def H(self, f):
        s = 1j * TWO_PI * np.asarray(f, float)
        return np.polyval(self.B, s) / np.polyval(self.Pc, s)


def complexification_report(f):
    cm = ComplexMode()
    Hp = np.abs(cm.H(f)) ** 2
    Hm = np.abs(cm.H(-f)) ** 2
    real = ca.unit_psd(4, C41J, f)
    exact = float(np.max(np.abs(0.5 * (Hp + Hm) / real - 1.0)))
    leak = Hm / Hp
    fb = [0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 1.0]
    lk = {f"{x:g}": float(np.interp(x, f, leak)) for x in fb}
    # Real-part law: stationary covariance of the realified mode, Re block,
    # against the real C41J Lyapunov solution (both unit gain).
    P8 = ca.lyap(cm.A, cm.W)
    m41 = ca.Carma(4, C41J)
    P4 = ca.lyap(m41.A, np.diag([0, 0, 1.0, 0]))
    poles = np.roots(cm.Pc)
    return {
        "Pc_coeffs_re": cm.Pc.real.tolist(), "Pc_coeffs_im": cm.Pc.imag.tolist(),
        "D_coeffs": cm.D.tolist(),
        "D_roots": [[float(r.real), float(r.imag)] for r in np.roots(cm.D)],
        "complex_poles": [[float(r.real), float(r.imag)] for r in poles],
        "max_rel_error_realpart_psd_vs_C41J": exact,
        "mirror_leakage_Hminus_over_Hplus": lk,
        # [p, v, a] only: the latent coordinates differ (C41J: the input r of
        # L2; complexified: the companion state, the complex primitive).
        "realpart_stationary_cov_pva_max_rel_dev_vs_C41J":
            float(np.max(np.abs(P8[:3, :3] - P4[:3, :3])) / np.max(np.abs(P4[:3, :3]))),
        "imag_real_cross_cov_p_a": float(P8[4, 2]),
        "hurwitz_max_real": float(np.max(np.linalg.eigvals(cm.A).real)),
    }


# ---------------------------------------------------------------------------
# 2. Three-axis systems
# ---------------------------------------------------------------------------

def rao_tf(f, fn=HEAVE_RAO[0], zeta=HEAVE_RAO[1]):
    """Second-order heave RAO wn^2/(s^2 + 2 zeta wn s + wn^2) (no averaging)."""
    wn = TWO_PI * fn
    s = 1j * TWO_PI * np.asarray(f, float)
    return wn * wn / (s * s + 2 * zeta * wn * s + wn * wn)


class System:
    """Stacked linear system: dx = A x dt + dB, Cov(dB) = W dt; p = Cp x,
    a = Ca x (world NED order x_N, y_E, z; the analysis uses z up)."""

    def __init__(self, A, W, Cp, Ca, name):
        self.A, self.W, self.Cp, self.Ca, self.name = A, W, Cp, Ca, name
        self.n = A.shape[0]

    def spectrum(self, f):
        """One-sided per-Hz cross-spectral matrix of p (complex, 3 x 3)."""
        S = np.empty((len(f), 3, 3), complex)
        I = np.eye(self.n)
        for k, fk in enumerate(f):
            Rk = np.linalg.solve(1j * TWO_PI * fk * I - self.A, I)
            Ck = self.Cp @ Rk
            S[k] = 2.0 * Ck @ self.W @ np.conj(Ck).T
        return S

    def stationary(self):
        return ca.lyap(self.A, self.W)


def herm_noise(Sigma, nb, ia_re, ia_im):
    """Real noise intensity of three complexified blocks (nb states each)
    driven by circular complex noise with E[eps eps^H] = 2 Sigma (Re-part
    variance Sigma_ii).  ia_re/ia_im: in-block indices of Re a and Im a."""
    W = np.zeros((3 * nb, 3 * nb))
    for i in range(3):
        for j in range(3):
            s = Sigma[i, j]
            W[i * nb + ia_re, j * nb + ia_re] = s.real
            W[i * nb + ia_im, j * nb + ia_im] = s.real
            W[i * nb + ia_im, j * nb + ia_re] = s.imag
            W[i * nb + ia_re, j * nb + ia_im] = -s.imag
    return W


def sys_ind_cor(kappa, Sigma_w, g2_unit):
    """IND/COR: C41J per axis; Sigma_w real 3 x 3 of acceleration variances."""
    m = ca.Carma(4, C41J, kappa, 1.0)
    A = np.kron(np.eye(3), m.A)
    W1 = np.zeros((4, 4)); W1[2, 2] = m.g2
    W = np.kron(np.real(Sigma_w), W1)
    Cp = np.kron(np.eye(3), np.eye(4)[[0]])
    Ca_ = np.kron(np.eye(3), np.eye(4)[[2]])
    return System(A, W, Cp, Ca_, "COR")


def sys_shm(kappa, Sigma_c):
    """SHM: complexified C41J per axis, Hermitian Sigma_c (acceleration
    variance units on the diagonal)."""
    cm = ComplexMode(C41J, kappa)
    m = ca.Carma(4, C41J, kappa, 1.0)       # unit-acceleration gain of the real law
    A = np.kron(np.eye(3), cm.A)
    W = herm_noise(np.asarray(Sigma_c) * m.g2, 8, 2, 6)
    Cp = np.kron(np.eye(3), np.eye(8)[[0]])
    Ca_ = np.kron(np.eye(3), np.eye(8)[[2]])
    return System(A, W, Cp, Ca_, "SHM")


def rao_block(cm, fn, zeta):
    """Complexified mode (8) followed by the heave RAO acting on Re p (2 more
    states [z, z']): z'' = wn^2 (p_mode - z) - 2 zeta wn z'.  Output z, z',
    z'' = wn^2 (x0 - z) - 2 zeta wn z'."""
    wn = TWO_PI * fn
    A = np.zeros((10, 10))
    A[:8, :8] = cm.A
    A[8, 9] = 1.0
    A[9, 0] = wn * wn
    A[9, 8] = -wn * wn
    A[9, 9] = -2 * zeta * wn
    cp = np.eye(10)[8]
    cacc = np.zeros(10); cacc[0] = wn * wn; cacc[8] = -wn * wn; cacc[9] = -2 * zeta * wn
    return A, cp, cacc


def rao_variance_gain(kappa, fn, zeta):
    """Var(z'') / Var(Re p_mode'') of the RAO block at unit mode gain."""
    cm = ComplexMode(C41J, kappa)
    A, cp, cacc = rao_block(cm, fn, zeta)
    W = np.zeros((10, 10)); W[2, 2] = 1.0; W[6, 6] = 1.0
    P = ca.lyap(A, W)
    return float(cacc @ P @ cacc / P[2, 2])


def sys_mca(kappa, Sigma_c, fn=HEAVE_RAO[0], zeta=HEAVE_RAO[1]):
    """MCA: SHM on x, y; the vertical complexified mode passes the heave RAO.
    Sigma_c[2,2] is the vertical MODE variance (the caller divides the
    requested heave acceleration variance by rao_variance_gain)."""
    cm = ComplexMode(C41J, kappa)
    m = ca.Carma(4, C41J, kappa, 1.0)
    Az, cpz, caz = rao_block(cm, fn, zeta)
    n = 8 + 8 + 10
    A = np.zeros((n, n))
    A[:8, :8] = cm.A; A[8:16, 8:16] = cm.A; A[16:, 16:] = Az
    W = np.zeros((n, n))
    offs = [0, 8, 16]
    S = np.asarray(Sigma_c) * m.g2
    for i in range(3):
        for j in range(3):
            s = S[i, j]
            oi, oj = offs[i], offs[j]
            W[oi + 2, oj + 2] = s.real; W[oi + 6, oj + 6] = s.real
            W[oi + 6, oj + 2] = s.imag; W[oi + 2, oj + 6] = -s.imag
    Cp = np.zeros((3, n)); Ca_ = np.zeros((3, n))
    Cp[0, 0] = 1; Cp[1, 8] = 1; Cp[2, 16:] = cpz
    Ca_[0, 2] = 1; Ca_[1, 10] = 1; Ca_[2, 16:] = caz
    return System(A, W, Cp, Ca_, "MCA")


def sys_shr(kappa, beta, sig_res, sig_z):
    """SHR: complexified z mode (8) + C41J residuals on x, y (4 + 4);
    p_x = Re(beta_x zeta_z) + r_x.  beta complex (2,), sig_res residual
    acceleration std (2,), sig_z vertical acceleration std."""
    cm = ComplexMode(C41J, kappa)
    m = ca.Carma(4, C41J, kappa, 1.0)
    n = 16
    A = np.zeros((n, n)); A[:8, :8] = cm.A; A[8:12, 8:12] = m.A; A[12:, 12:] = m.A
    W = np.zeros((n, n))
    W[2, 2] = W[6, 6] = m.g2 * sig_z ** 2
    W[10, 10] = m.g2 * sig_res[0] ** 2
    W[14, 14] = m.g2 * sig_res[1] ** 2
    Cp = np.zeros((3, n)); Ca_ = np.zeros((3, n))
    for i, (o, b) in enumerate(((8, beta[0]), (12, beta[1]))):
        Cp[i, 0] = b.real; Cp[i, 4] = -b.imag; Cp[i, o] = 1.0
        Ca_[i, 2] = b.real; Ca_[i, 6] = -b.imag; Ca_[i, o + 2] = 1.0
    Cp[2, 0] = 1.0; Ca_[2, 2] = 1.0
    return System(A, W, Cp, Ca_, "SHR")


def sys_k2(kappas, Sigmas):
    """K2: two complexified modes per axis (kappa_k, Hermitian Sigma_k)."""
    n = 48
    A = np.zeros((n, n)); W = np.zeros((n, n))
    Cp = np.zeros((3, n)); Ca_ = np.zeros((3, n))
    for k, (kap, Sig) in enumerate(zip(kappas, Sigmas)):
        cm = ComplexMode(C41J, kap)
        m = ca.Carma(4, C41J, kap, 1.0)
        base = 24 * k
        for i in range(3):
            o = base + 8 * i
            A[o:o + 8, o:o + 8] = cm.A
            Cp[i, o] = 1.0; Ca_[i, o + 2] = 1.0
        Wk = herm_noise(np.asarray(Sig) * m.g2, 8, 2, 6)
        W[base:base + 24, base:base + 24] = Wk
    return System(A, W, Cp, Ca_, "K2")


# ---------------------------------------------------------------------------
# 2b. Complex-OU modal decomposition of C41J (Fasen-Hartmann & Scholz)
# ---------------------------------------------------------------------------

def latent_readout(A, b):
    """Row c with c (sI - A)^-1 b = 1/P(s): Markov parameters (0, 0, 0, 1)."""
    K = np.column_stack([b, A @ b, A @ A @ b, A @ A @ A @ b])
    return np.linalg.solve(K.T, np.array([0.0, 0.0, 0.0, 1.0]))


def modal_c41(kappa=1.0, sh=C41J):
    """Real modal form of C41J (CarmaModalMotion): two complex OU modes
    z_k with dz_k = p_k z_k dt + r_k dW (one common real noise), r_k the
    residues of B/P at the upper poles, p = 2 Re(z_1 + z_2).  Real state
    [Re z1, Im z1, Re z2, Im z2]; unit noise gain."""
    B, P = c41_poly(sh, kappa)
    roots = np.roots(P)
    up = sorted(roots[np.imag(roots) > 0], key=lambda z: abs(z))
    dP = np.polyder(P)
    res = [np.polyval(B, pk) / np.polyval(dP, pk) for pk in up]
    A = np.zeros((4, 4)); G = np.zeros(4)
    for k, (pk, rk) in enumerate(zip(up, res)):
        A[2 * k:2 * k + 2, 2 * k:2 * k + 2] = [[pk.real, -pk.imag], [pk.imag, pk.real]]
        G[2 * k:2 * k + 2] = [rk.real, rk.imag]
    cp = np.array([2.0, 0.0, 2.0, 0.0])
    return {"poles": up, "residues": res, "A": A, "G": G, "cp": cp, "B": B, "P": P}


def modal_report(f):
    """Exactness, conditioning and kinematic structure of the modal form."""
    out = {}
    m = modal_c41(1.0)
    out["poles_kappa1"] = [[float(z.real), float(z.imag)] for z in m["poles"]]
    out["residues_kappa1"] = [[float(z.real), float(z.imag)] for z in m["residues"]]
    r = np.array(m["residues"]); pk = np.array(m["poles"])
    out["sum_re_r"] = float(np.sum(r.real))
    out["sum_re_rp"] = float(np.sum((r * pk).real))
    out["sum_im_r"] = float(np.sum(r.imag))
    out["sum_re_rp2"] = float(np.sum((r * pk * pk).real))
    rows = []
    for kappa in np.geomspace(0.2, 2.5, 13):
        mm = modal_c41(kappa)
        base = ca.Carma(4, C41J, kappa, 1.0)
        A, G, cp = mm["A"], mm["G"], mm["cp"]
        Ob = np.vstack([np.eye(4)[0] @ np.linalg.matrix_power(base.A, i) for i in range(4)])
        Om = np.vstack([cp @ np.linalg.matrix_power(A, i) for i in range(4)])
        T = np.linalg.solve(Ob, Om)
        sv = np.linalg.svd(T, compute_uv=False)
        # transfer function, stationary covariance, discrete F/Q, autocovariance
        s = 1j * TWO_PI * kappa * np.array([0.03, 0.1, 0.3, 0.6, 1.0, 2.0])
        Hm = np.array([cp @ np.linalg.solve(sk * np.eye(4) - A, G) for sk in s])
        Hb = np.polyval(mm["B"], s) / np.polyval(mm["P"], s)
        Wm = np.outer(G, G)
        Wb = np.zeros((4, 4)); Wb[2, 2] = 1.0
        Pm, Pb = ca.lyap(A, Wm), ca.lyap(base.A, Wb)
        Fm, Qm = mma.van_loan(A, Wm, DT)
        Fb, Qb = mma.van_loan(base.A, Wb, DT)
        Ti = np.linalg.inv(T)
        ac = []
        for lag in (0.5, 2.0, 5.0):
            Fl = mma.van_loan(A, 0 * Wm, lag)[0]
            Fbl = mma.van_loan(base.A, 0 * Wb, lag)[0]
            ac.append(abs(cp @ Fl @ Pm @ cp - (Fbl @ Pb)[0, 0]) / Pb[0, 0])
        rows.append({"kappa": float(kappa), "cond_T": float(sv[0] / sv[-1]),
                     "tf_rel_err": float(np.max(np.abs(Hm / Hb - 1.0))),
                     "stationary_rel_err": float(np.max(np.abs(T @ Pm @ T.T - Pb)) / np.max(np.abs(Pb))),
                     "F_err": float(np.max(np.abs(Ti @ Fb @ T - Fm))),
                     "Q_rel_err": float(np.max(np.abs(Ti @ Qb @ Ti.T - Qm)) / np.max(np.abs(Qm))),
                     "autocov_rel_err": float(max(ac))})
    out["equivalence_vs_kappa"] = rows
    # Kinematic structure: c G = c A G = 0; rotating the modal readout by a
    # constant phase phi (p_phi = 2 Re(e^{i phi} (z1 + z2))) breaks it.
    A, G, cp = m["A"], m["G"], m["cp"]
    norm = abs(cp @ A @ A @ G)
    rot = {}
    for deg in (5, 15, 30, 60, 90):
        c, sn = np.cos(np.radians(deg)), np.sin(np.radians(deg))
        cr = np.array([2 * c, -2 * sn, 2 * c, -2 * sn])
        rot[str(deg)] = {"cG_over_cA2G": float(abs(cr @ G) / norm), "cAG_over_cA2G": float(abs(cr @ A @ G) / norm)}
    out["kinematic_cG_cAG"] = [float(abs(cp @ G) / norm), float(abs(cp @ A @ G) / norm)]
    out["rotated_readout"] = rot
    # Independently driven oscillators (each mode its own noise, same
    # magnitudes): Q_pp slope in h, and the displacement/acceleration PSD
    # high-frequency slopes against the common-noise (exact) model.
    Wc = np.outer(G, G)
    Wi = np.zeros((4, 4)); Wi[:2, :2] = Wc[:2, :2]; Wi[2:, 2:] = Wc[2:, 2:]
    sl = {}
    for tag, W in (("common", Wc), ("independent", Wi)):
        q = [cp @ mma.van_loan(A, W, h)[1] @ cp for h in (0.005, 0.01)]
        fh = np.array([5.0, 10.0])
        Sp = [np.real(cp @ np.linalg.solve(1j * TWO_PI * x * np.eye(4) - A, np.eye(4)) @ W
                      @ np.conj(cp @ np.linalg.solve(1j * TWO_PI * x * np.eye(4) - A, np.eye(4)))) for x in fh]
        sl[tag] = {"Qpp_log2_slope": float(np.log2(q[1] / q[0])),
                   "p_psd_hf_log2_slope": float(np.log2(Sp[1] / Sp[0])),
                   "a_psd_hf_log2_slope": float(np.log2(Sp[1] / Sp[0]) + 4.0),
                   "var_p_ratio_vs_common": None}
        Pw = ca.lyap(A, W)
        sl[tag]["var_p"] = float(cp @ Pw @ cp)
    for tag in sl:
        sl[tag]["var_p_ratio_vs_common"] = sl[tag]["var_p"] / sl["common"]["var_p"]
    out["independent_driving"] = sl
    # Admissible cross-axis structure: per-mode complex gains beta_k on a
    # second axis driven by the same noise keep relative degree 3 iff
    # sum Re(beta_k r_k) = sum Re(beta_k r_k p_k) = 0; with a common phase
    # beta_k = e^{i phi} the defect is sin(phi) |sum Im r_k|.
    out["admissible_family"] = (
        "B_x(s) = b1 B(s) + b0 (real b0, b1): x = B_x/P xi; phase(S_xz) = -atan(b0/(b1 w)) "
        "for w_z0 = 0; a constant phase is not admissible; b0 = 0 leaves real coupling only")
    return out


def sys_modn(kappa, Sig_c, f_op, phase_clip_deg=70.0):
    """MODN: the 4-state C41J/modal block per axis with a real excitation and
    a per-axis numerator zero, B_x = B + b0_x (CarmaModalMotion, Form
    Numerator).  b0 from the phase of the Hermitian excitation entry and its
    coherence; Sigma~ = Re(D Sigma_c D^H), D = diag(1/u_x, 1/u_y, 1),
    u = 1 - i b0 / w_op (axis gain at the band).  Returned in baseline
    coordinates (kinematically identical to the modal form)."""
    m = ca.Carma(4, C41J, kappa, 1.0)
    b = np.zeros(4); b[2] = 1.0
    cl = latent_readout(m.A, b)
    w = TWO_PI * f_op
    Sig_c = np.asarray(Sig_c, complex)
    b0 = np.zeros(3)
    u = np.ones(3, complex)
    for i in (0, 1):
        sxz = Sig_c[i, 2]
        den = np.sqrt(max(np.real(Sig_c[i, i]) * np.real(Sig_c[2, 2]), 1e-300))
        rho = min(abs(sxz) / den, 1.0)
        if abs(sxz) == 0.0 or np.real(sxz) == 0.0 and np.imag(sxz) == 0.0:
            continue
        t = np.clip(np.imag(sxz) / np.real(sxz) if np.real(sxz) != 0 else np.sign(np.imag(sxz)) * 1e9,
                    -np.tan(np.radians(phase_clip_deg)), np.tan(np.radians(phase_clip_deg)))
        b0[i] = -w * t * rho
        u[i] = 1.0 - 1j * b0[i] / w
    D = np.diag(1.0 / u)
    Sig_t = np.real(D @ Sig_c @ np.conj(D).T)
    A = np.kron(np.eye(3), m.A)
    W1 = np.zeros((4, 4)); W1[2, 2] = m.g2
    W = np.kron(Sig_t, W1)
    Cp = np.zeros((3, 12)); Ca_ = np.zeros((3, 12))
    for i in range(3):
        c = np.eye(4)[0] + b0[i] * cl
        Cp[i, 4 * i:4 * i + 4] = c
        Ca_[i, 4 * i:4 * i + 4] = c @ m.A @ m.A
    sy = System(A, W, Cp, Ca_, "MODN")
    sy.b0, sy.Sigma_t = b0, Sig_t
    return sy


def system_checks(f):
    """Realization, Kronecker, PSD and Lyapunov checks of every model."""
    kappa = 0.33 / ca.Carma(4, C41J).base_fz_hz
    rng = np.random.default_rng(7)
    X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    Sig_c = X @ np.conj(X).T / 3.0
    Sig_w = np.real(Sig_c)
    out = {}
    cor = sys_ind_cor(kappa, Sig_w, None)
    shm = sys_shm(kappa, Sig_c)
    shm_real = sys_shm(kappa, Sig_w.astype(complex))
    mca = sys_mca(kappa, Sig_c)
    # Kronecker identity of Experiment A: joint Van Loan vs Sigma (x) Q1.
    F, Q = mma.van_loan(cor.A, cor.W, DT)
    m = ca.Carma(4, C41J, kappa, 1.0)
    F1, Q1 = mma.van_loan(m.A, np.diag([0, 0, m.g2, 0]), DT)
    out["kronecker_F_err"] = float(np.max(np.abs(F - np.kron(np.eye(3), F1))))
    out["kronecker_Q_rel_err"] = float(np.max(np.abs(Q - np.kron(Sig_w, Q1))) / np.max(np.abs(Q)))
    # Hermitian Kronecker of Experiment B: Q(i,j) = Re(S_ij) Q8 + Im(S_ij) J Q8.
    cm = ComplexMode(C41J, kappa)
    F8, Q8 = mma.van_loan(cm.A, cm.W * m.g2, DT)
    J = np.block([[np.zeros((4, 4)), -np.eye(4)], [np.eye(4), np.zeros((4, 4))]])
    Fb, Qb = mma.van_loan(shm.A, shm.W, DT)
    Qk = np.zeros_like(Qb)
    for i in range(3):
        for j in range(3):
            s = Sig_c[i, j]
            Qk[8 * i:8 * i + 8, 8 * j:8 * j + 8] = s.real * Q8 + s.imag * (J @ Q8)
    out["hermitian_kronecker_Q_rel_err"] = float(np.max(np.abs(Qb - Qk)) / np.max(np.abs(Qb)))
    out["J_commutes_with_F8"] = float(np.max(np.abs(J @ F8 - F8 @ J)))
    # SHM with a real Sigma reproduces COR's cross-spectral matrix.
    fs = np.geomspace(0.03, 3.0, 80)
    Scor, Sshr = cor.spectrum(fs), shm_real.spectrum(fs)
    out["shm_real_sigma_vs_cor_spectrum_rel_err"] = float(np.max(np.abs(Scor - Sshr)) / np.max(np.abs(Scor)))
    # PSD and Hermitian of every model's cross-spectral matrix.
    for name, sy in (("COR", cor), ("SHM", shm), ("MCA", mca)):
        S = sy.spectrum(fs)
        chk = ca.psd_check(S)
        P = sy.stationary()
        ev = np.linalg.eigvalsh(0.5 * (P + P.T))
        out[name] = {"n": sy.n, "hurwitz_max_real": float(np.max(np.linalg.eigvals(sy.A).real)),
                     "spectrum_hermitian_defect": chk["max_hermitian_defect"],
                     "spectrum_min_eig_over_trace": chk["min_eig_over_trace"],
                     "stationary_min_eig_over_max": float(ev[0] / ev[-1]),
                     # stationary covariance vs quadrature of the spectrum (p_x)
                     "var_px_lyap_vs_quadrature": float(
                         (sy.Cp @ P @ sy.Cp.T)[0, 0] /
                         np.trapezoid(np.real(sy.spectrum(np.linspace(1e-4, 20.0, 40001))[:, 0, 0]),
                                      np.linspace(1e-4, 20.0, 40001)))}
        # Kinematic consistency: v = dp/dt and a = dv/dt are noise-free
        # derivatives (Cp G = 0, Cp A G = 0, Cp A^2 G != 0).
        G = sy.W
        out[name]["CpG"] = float(np.max(np.abs(sy.Cp @ G)))
        out[name]["CpAG"] = float(np.max(np.abs(sy.Cp @ sy.A @ G)))
        out[name]["Ca_equals_CpA2"] = float(np.max(np.abs(sy.Cp @ sy.A @ sy.A - sy.Ca)))
    return out


# ---------------------------------------------------------------------------
# 3. Truth cross-spectra and representable fractions
# ---------------------------------------------------------------------------

def truth(sea, f, heading_deg=0.0):
    """[x_N, y_E, z_up] cross-spectral matrix of the benchmark truth."""
    return ca.cross_spectrum(SEAS[sea], f, heading_deg)


def turned(sea, f, heading_deg, dturn):
    """Truth of the same sea with every wave direction turned by dturn."""
    parts = [(hs, fp, g, d + dturn, sp) for hs, fp, g, d, sp in SEAS[sea]]
    return ca.cross_spectrum(parts, f, heading_deg)


def operating_point(S, f):
    return ca.operating_point(np.real(S[:, 2, 2]), f)


def demod_weight(f, f_op, bw_frac=0.25):
    """|L(f - f_op)|^2 of the estimator's complex-envelope low-pass
    (critically damped second order, corner bw_frac f_op), plus the image."""
    fc = bw_frac * f_op
    L = lambda d: 1.0 / (1.0 + (d / fc) ** 2) ** 2
    return L(f - f_op), L(f + f_op)


def band_cross(S, f, f_op, bw_frac=0.25):
    """Demodulated band cross-spectral matrix of the ACCELERATIONS (what the
    estimator converges to): sum_f w(f) w^4 S(f) df (+ image, conjugated)."""
    w4 = (TWO_PI * f) ** 4
    wp, wm = demod_weight(f, f_op, bw_frac)
    Sa = w4[:, None, None] * S
    M = np.trapezoid(wp[:, None, None] * Sa, f, axis=0) + np.trapezoid(wm[:, None, None] * np.conj(Sa), f, axis=0)
    return 0.5 * (M + np.conj(M).T)


def representable(S, f):
    """Energy-weighted fractions of the x-z, y-z and x-y cross-spectra that a
    model family can represent at all: real (COR) vs Hermitian with constant
    phase (SHM).  Weighted by acceleration band energy."""
    w4 = (TWO_PI * f) ** 4
    out = {}
    for name, (i, j) in (("xz", (0, 2)), ("yz", (1, 2)), ("xy", (0, 1))):
        s = w4 * S[:, i, j]
        e = np.trapezoid(np.abs(s) ** 2, f)
        if e <= 0:
            continue
        # COR: best real multiple of a common real shape -> only Re survives.
        re_frac = np.trapezoid(np.real(s) ** 2, f) / e
        # SHM: best constant complex coefficient times a real shape.
        mag = np.abs(s)
        c = np.trapezoid(s * mag, f) / np.trapezoid(mag ** 2, f)
        shm_frac = 1.0 - np.trapezoid(np.abs(s - c * mag) ** 2, f) / e
        coh, ph = ca.coherence(S, i, j)
        wgt = w4 * np.sqrt(np.real(S[:, i, i]) * np.real(S[:, j, j]))
        ok = np.isfinite(coh)
        out[name] = {"real_part_energy_fraction": float(re_frac),
                     "constant_phase_fit_energy_fraction": float(shm_frac),
                     "energy_weighted_coherence": float(np.trapezoid(np.where(ok, coh, 0) * wgt, f) / np.trapezoid(wgt, f)),
                     "phase_deg_at_peak": float(ph[np.argmax(np.where(ok, wgt, 0))]),
                     "phase_deg_p10_p90": [float(x) for x in np.nanpercentile(np.where(wgt > 0.1 * wgt.max(), ph, np.nan), [10, 90])]}
    return out


# ---------------------------------------------------------------------------
# 4. Oracle parameters (Sigma = c S_hat + Delta)
# ---------------------------------------------------------------------------

def map_sigma(S_hat, sig_prior, R=None):
    """Hermitian excitation from a band cross-spectrum (CarmaCoupling::excitation,
    default horizontal block, no sampling shrinkage): the baseline diagonal,
    cross terms b = c S_hat_hz with c = |R_z|^2 sig_z^2 / S_hat_zz (the measured
    regression on z), scaled down if b^H H^-1 b would exceed sig_z^2; then
    Sigma = R^-1 M R^-H for a model whose axis i responds with R_i at the band."""
    R = np.ones(3, complex) if R is None else np.asarray(R, complex)
    mzz = abs(R[2]) ** 2 * sig_prior[2] ** 2
    c = mzz / np.real(S_hat[2, 2])
    M = np.diag([abs(R[0]) ** 2 * sig_prior[0] ** 2, abs(R[1]) ** 2 * sig_prior[1] ** 2, mzz]).astype(complex)
    b = c * np.array(S_hat, complex)[:2, 2]
    q = np.real(np.conj(b) @ np.linalg.solve(M[:2, :2], b))
    if q > 0.999 * mzz:
        b *= np.sqrt(0.999 * mzz / q)
    M[:2, 2] = b
    M[2, :2] = np.conj(b)
    Rinv = np.diag(1.0 / R)
    return Rinv @ M @ np.conj(Rinv).T


def model_params(sea, f, heading_deg=0.0):
    S = truth(sea, f, heading_deg)
    op = operating_point(S, f)
    base_fz = ca.Carma(4, C41J).base_fz_hz
    kappa = op["fz_hz"] / base_fz
    sig = np.full(3, op["sigma_a"])          # wrapper mapping: isotropic
    Sh = band_cross(S, f, op["fz_hz"])
    Sig_c = map_sigma(Sh, sig)
    Sig_w = np.real(Sig_c)
    # MCA: vertical response of the RAO block at the band.
    Rz = complex(rao_tf(np.array([op["fz_hz"]]))[0])
    vg = rao_variance_gain(kappa, *HEAVE_RAO[:2])
    sig_m = sig.copy(); sig_m[2] = sig[2] / np.sqrt(vg)      # mode std so heave matches
    Sig_mca = map_sigma(Sh, sig_m, R=np.array([1.0, 1.0, Rz]))
    return S, op, kappa, sig, Sh, Sig_c, Sig_w, Sig_mca


# ---------------------------------------------------------------------------
# 5. Steady-state evaluation against the true cross-spectral matrix
# ---------------------------------------------------------------------------

def steady_eval(sy, S, f, h_noise=1.0):
    """Measurement-updated steady-state filter on y = a + n (3 axes);
    true error PSD E = (I + w^2 T) S (I + w^2 T)^H + T N T^H.  Returns RMS
    per axis, predicted std per axis and band RMS."""
    from scipy.linalg import solve_discrete_are
    F, Q = mma.van_loan(sy.A, sy.W, DT)
    Q = 0.5 * (Q + Q.T) + 1e-18 * np.eye(sy.n)
    Rm = R_FILTER * np.diag([h_noise ** 2, h_noise ** 2, 1.0])
    P = solve_discrete_are(F.T, sy.Ca.T, Q, Rm)
    Sinn = sy.Ca @ P @ sy.Ca.T + Rm
    K = P @ sy.Ca.T @ np.linalg.inv(Sinn)
    Ppost = (np.eye(sy.n) - K @ sy.Ca) @ P
    Acl = (np.eye(sy.n) - K @ sy.Ca) @ F
    w = TWO_PI * f
    Nd = N_ONE_SIDED * np.array([h_noise ** 2, h_noise ** 2, 1.0])
    E = np.empty((len(f), 3))
    I = np.eye(sy.n)
    for k, zk in enumerate(np.exp(1j * w * DT)):
        T = sy.Cp @ np.linalg.solve(I - Acl / zk, K)          # 3 x 3
        G = np.eye(3) + w[k] ** 2 * T
        Ek = G @ S[k] @ np.conj(G).T + T @ np.diag(Nd) @ np.conj(T).T
        E[k] = np.real(np.diag(Ek))
    rms = np.sqrt(np.trapezoid(E, f, axis=0))
    pred = np.sqrt(np.diag(sy.Cp @ Ppost @ sy.Cp.T))
    bands = {}
    for a, b in ca.BANDS:
        m = (f >= a) & (f < b)
        bands[f"{a}-{b}"] = np.sqrt(np.trapezoid(np.where(m[:, None], E, 0), f, axis=0)).tolist()
    return {"rms_m": rms.tolist(), "pred_std_m": pred.tolist(),
            "var_ratio": (rms ** 2 / pred ** 2).tolist(), "bands": bands}, E


def steady_study(f, seas):
    out = {}
    curves = {}
    for sea in seas:
        heading = SEA_HEADING.get(sea, 0.0)
        S, op, kappa, sig, Sh, Sig_c, Sig_w, Sig_mca = model_params(sea, f, heading)
        sys = {
            "IND": sys_ind_cor(kappa, np.diag(sig ** 2), None),
            "COR_xy": sys_ind_cor(kappa, np.block([[Sig_w[:2, :2], np.zeros((2, 1))],
                                                   [np.zeros((1, 2)), Sig_w[2:, 2:]]]), None),
            "COR": sys_ind_cor(kappa, Sig_w, None),
            "SHM": sys_shm(kappa, Sig_c),
            "SHM_diag": sys_shm(kappa, np.diag(np.real(np.diag(Sig_c))).astype(complex)),
            "MCA": sys_mca(kappa, Sig_mca),
            "MCA_diag": sys_mca(kappa, np.diag(np.real(np.diag(Sig_mca))).astype(complex)),
        }
        # Deliberately wrong coupling: the cross terms of the wave field
        # turned by 180 deg (propagation reversed: quadrature sign flipped)
        # and by 90 deg, the hull heading unchanged.
        S180 = turned(sea, f, heading, 180.0)
        S90 = turned(sea, f, heading, 90.0)
        for tag, Sx in (("SHM_wrong180", S180), ("SHM_wrong90", S90)):
            Shx = band_cross(Sx, f, op["fz_hz"])
            sys[tag] = sys_shm(kappa, map_sigma(Shx, sig))
        # Modal numerator coupling (kinematically admissible phase).
        sys["MODN"] = sys_modn(kappa, Sig_c, op["fz_hz"])
        sys["MODN_wrong180"] = sys_modn(kappa, map_sigma(band_cross(S180, f, op["fz_hz"]), sig), op["fz_hz"])
        # State-reduced shared mode: beta from the band regression on z.
        beta = np.array([Sh[0, 2] / Sh[2, 2], Sh[1, 2] / Sh[2, 2]])
        res = np.sqrt(np.maximum(sig[:2] ** 2 - np.abs(beta) ** 2 * sig[2] ** 2, 0.05 * sig[:2] ** 2))
        sys["SHR"] = sys_shr(kappa, beta, res, sig[2])
        res_sea = {"op": op, "kappa": kappa, "band_cross_accel": {"re": np.real(Sh).tolist(), "im": np.imag(Sh).tolist()},
                   "Sigma_c": {"re": np.real(Sig_c).tolist(), "im": np.imag(Sig_c).tolist()},
                   "Sigma_mca": {"re": np.real(Sig_mca).tolist(), "im": np.imag(Sig_mca).tolist()},
                   "representable": representable(S, f), "models": {},
                   "modn_b0_over_wop": (sys["MODN"].b0 / (TWO_PI * op["fz_hz"])).tolist()}
        for hn in (1.0, 3.0):
            key = f"hnoise_x{hn:g}"
            res_sea["models"][key] = {}
            for name, sy in sys.items():
                r, E = steady_eval(sy, S, f, hn)
                r["n_states"] = sy.n
                res_sea["models"][key][name] = r
                if hn == 1.0:
                    curves[(sea, name)] = E
            wb = ca.wiener_bound(S, f, hn)
            res_sea["models"][key]["wiener"] = wb
        # Two-mode oracle on two-system seas: one mode per part, each part's
        # own band cross-spectrum and frequency.
        parts = SEAS[sea]
        if len(parts) == 2:
            kap2, sig2 = [], []
            for p in parts:
                Sp = ca.cross_spectrum([p], f, heading)
                opp = operating_point(Sp, f)
                kap2.append(opp["fz_hz"] / ca.Carma(4, C41J).base_fz_hz)
                Shp = band_cross(Sp, f, opp["fz_hz"])
                sp = np.full(3, opp["sigma_a"])
                sig2.append(map_sigma(Shp, sp))
            for hn in (1.0, 3.0):
                k2 = sys_k2(kap2, sig2)
                k2d = sys_k2(kap2, [np.diag(np.real(np.diag(s))).astype(complex) for s in sig2])
                for tag, sy in (("K2", k2), ("K2_diag", k2d)):
                    r, E = steady_eval(sy, S, f, hn)
                    r["n_states"] = sy.n
                    res_sea["models"][f"hnoise_x{hn:g}"][tag] = r
                    if hn == 1.0:
                        curves[(sea, tag)] = E
        out[sea] = res_sea
        print(f"  steady state: {sea} done", flush=True)
    return out, curves


# ---------------------------------------------------------------------------
# 6. Identifiability (Whittle Fisher information)
# ---------------------------------------------------------------------------

def whittle_fim(spec_fn, theta, f, T_s, h_noise=1.0, rel_step=1e-4):
    """FIM_jk = T int_0^inf tr(S^-1 dS_j S^-1 dS_k) df for a real 3-dim
    process with one-sided cross-spectral matrix S(f) = w^4 S_p + N (per Hz,
    T seconds; the one-sided form counts each frequency once)."""
    w4 = (TWO_PI * f) ** 4
    Nd = np.diag(N_ONE_SIDED * np.array([h_noise ** 2, h_noise ** 2, 1.0]))

    def Sy(th):
        return w4[:, None, None] * spec_fn(th) + Nd[None]

    S0 = Sy(theta)
    Si = np.linalg.inv(S0)
    d = []
    for j in range(len(theta)):
        h = rel_step * max(1.0, abs(theta[j]))
        tp = np.array(theta, float); tp[j] += h
        tm = np.array(theta, float); tm[j] -= h
        d.append((Sy(tp) - Sy(tm)) / (2 * h))
    n = len(theta)
    F = np.zeros((n, n))
    for j in range(n):
        Aj = Si @ d[j]
        for k in range(j, n):
            v = np.real(np.einsum("fab,fba->f", Aj, Si @ d[k]))
            F[j, k] = F[k, j] = T_s * np.trapezoid(v, f)
    return F


def mcarma41_spectrum(theta, f, n=3):
    """Unrestricted real MCARMA(4,1): P(s) = s^4 I + A1 s^3 + A2 s^2 + A3 s + A4,
    B(s) = B0 s + B1, Sigma = I (B absorbs the noise factor).  theta packs
    A1..A4 (4 n^2), B0, B1 (2 n^2)."""
    k = n * n
    A = [np.asarray(theta[i * k:(i + 1) * k]).reshape(n, n) for i in range(4)]
    B0 = np.asarray(theta[4 * k:5 * k]).reshape(n, n)
    B1 = np.asarray(theta[5 * k:6 * k]).reshape(n, n)
    S = np.empty((len(f), n, n), complex)
    for i, fk in enumerate(f):
        s = 1j * TWO_PI * fk
        P = s ** 4 * np.eye(n) + A[0] * s ** 3 + A[1] * s ** 2 + A[2] * s + A[3]
        H = np.linalg.solve(P, B0 * s + B1)
        S[i] = 2.0 * H @ np.conj(H).T
    return S


def identifiability(f, T_s=1200.0):
    """Fisher information of each family's parameters at the JONSWAP oracle
    point, and of an unrestricted MCARMA(4,1) at the equivalent independent
    point plus a small generic perturbation."""
    sea = "jonswap_fp0.29_hs0.8_g3.3"
    S, op, kappa, sig, Sh, Sig_c, Sig_w, Sig_mca = model_params(sea, f)
    out = {"T_s": T_s, "sea": sea}
    fz = op["fz_hz"]

    def herm_pack(Sg):
        return [np.real(Sg[0, 0]), np.real(Sg[1, 1]), np.real(Sg[2, 2]),
                Sg[0, 1].real, Sg[0, 2].real, Sg[1, 2].real, Sg[0, 1].imag, Sg[0, 2].imag, Sg[1, 2].imag]

    def herm_unpack(v):
        Sg = np.diag(np.array(v[:3], complex))
        Sg[0, 1] = v[3] + 1j * v[6]; Sg[0, 2] = v[4] + 1j * v[7]; Sg[1, 2] = v[5] + 1j * v[8]
        Sg[1, 0], Sg[2, 0], Sg[2, 1] = np.conj(Sg[0, 1]), np.conj(Sg[0, 2]), np.conj(Sg[1, 2])
        return Sg

    def fam(name, spec_fn, theta, labels):
        F = whittle_fim(spec_fn, theta, f, T_s)
        ev = np.linalg.eigvalsh(F)
        try:
            crlb = np.sqrt(np.maximum(np.diag(np.linalg.pinv(F, rcond=1e-14)), 0))
        except np.linalg.LinAlgError:
            crlb = np.full(len(theta), np.nan)
        out[name] = {"n_params": len(theta), "labels": labels, "theta": [float(x) for x in theta],
                     "fim_eig_min": float(ev[0]), "fim_eig_max": float(ev[-1]),
                     "numerical_rank": int(np.sum(ev > 1e-9 * ev[-1])),
                     "crlb_std": [float(x) for x in crlb],
                     "crlb_rel_std": [float(c / abs(t)) if abs(t) > 1e-12 else None for c, t in zip(crlb, theta)]}

    # COR: log kappa, Sigma_w (6)
    def cor_spec(th):
        Sw = np.array([[th[1], th[4], th[5]], [th[4], th[2], th[6]], [th[5], th[6], th[3]]])
        return sys_ind_cor(np.exp(th[0]), Sw, None).spectrum(f)
    th = [np.log(kappa), Sig_w[0, 0], Sig_w[1, 1], Sig_w[2, 2], Sig_w[0, 1], Sig_w[0, 2], Sig_w[1, 2]]
    fam("COR", cor_spec, th, ["log_kappa", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"])

    # SHM: log kappa, Hermitian Sigma_c (9)
    def shm_spec(th):
        return sys_shm(np.exp(th[0]), herm_unpack(th[1:])).spectrum(f)
    th = [np.log(kappa)] + herm_pack(Sig_c)
    fam("SHM", shm_spec, th, ["log_kappa", "Sxx", "Syy", "Szz", "Re_xy", "Re_xz", "Re_yz", "Im_xy", "Im_xz", "Im_yz"])

    # MCA: log kappa, Sigma (9), RAO fn, zeta
    def mca_spec(th):
        return sys_mca(np.exp(th[0]), herm_unpack(th[1:10]), th[10], th[11]).spectrum(f)
    th = [np.log(kappa)] + herm_pack(Sig_mca) + [HEAVE_RAO[0], HEAVE_RAO[1]]
    fam("MCA", mca_spec, th, ["log_kappa", "Sxx", "Syy", "Szz", "Re_xy", "Re_xz", "Re_yz", "Im_xy", "Im_xz", "Im_yz",
                               "rao_fn_hz", "rao_zeta"])

    # Unrestricted MCARMA(4,1): the IND point (A_k = a_k I, B = g diag(s + 0))
    # plus a small generic perturbation so that no accidental symmetry
    # inflates the rank.
    _, P = c41_poly(C41J, kappa)
    g = np.sqrt(ca.Carma(4, C41J, kappa, op["sigma_a"]).g2)
    rng = np.random.default_rng(11)
    I3 = np.eye(3)
    th = []
    for kk in range(4):
        th += list((P[kk + 1] * (I3 + 0.05 * rng.normal(size=(3, 3)))).ravel())
    th += list((g * (I3 + 0.05 * rng.normal(size=(3, 3)))).ravel())        # B0 (s coefficient)
    th += list((g * 0.02 * rng.normal(size=(3, 3))).ravel())               # B1
    fam("MCARMA41_unrestricted", lambda t: mcarma41_spectrum(t, f), th,
        [f"A{k}_{i}{j}" for k in range(1, 5) for i in range(3) for j in range(3)] +
        [f"B0_{i}{j}" for i in range(3) for j in range(3)] + [f"B1_{i}{j}" for i in range(3) for j in range(3)])
    # The B0, B1 orthogonal rotation invariance (B -> B U, U orthogonal)
    # leaves S unchanged: 3 null directions expected.
    out["MCARMA41_unrestricted"]["expected_null_dims_rotation"] = 3
    out["operating_fz_hz"] = fz
    return out


# ---------------------------------------------------------------------------
# 7. Cost model
# ---------------------------------------------------------------------------

def cost(nl, blocks, ncol_motion):
    """Multiply-adds of one prediction + one accelerometer update with the
    generic core's block-diagonal F_LL (block sizes `blocks`) and dense Q_LL,
    nine shared states (attitude, gyro bias, accelerometer bias)."""
    nx = nl + 9
    pll = sum(2 * bi * bj * (bi + bj) / 2 for i, bi in enumerate(blocks) for j, bj in enumerate(blocks) if j >= i)
    pal = 6 * sum(b * b for b in blocks) + 6 * 6 * nl
    plb = 3 * sum(b * b for b in blocks)
    q = nl * (nl + 1) / 2
    ncol = 3 + ncol_motion + 3
    pct = nx * ncol * 3
    joseph = nx * (nx + 1) / 2 * 6
    return int(pll + pal + plb + q + pct + joseph + 400)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def make_plots(res, f, outdir, ext="svg"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "svg.hashsalt": "carma3-coupled"})
    colors = {"truth": "#222222", "COR": "#1f77b4", "SHM": "#d62728", "MCA": "#2ca02c", "IND": "#7f7f7f"}
    seas = ["jonswap_fp0.29_hs0.8_g3.3", "mixed_swell_wind", "bimodal", "crossing"]
    fig, axs = plt.subplots(3, len(seas), figsize=(12.5, 7.6), sharex=True)
    fs = np.geomspace(0.05, 1.5, 300)
    for c, sea in enumerate(seas):
        S, op, kappa, sig, Sh, Sig_c, Sig_w, Sig_mca = model_params(sea, fs)
        models = {"COR": sys_ind_cor(kappa, Sig_w, None).spectrum(fs),
                  "SHM": sys_shm(kappa, Sig_c).spectrum(fs),
                  "MCA": sys_mca(kappa, Sig_mca).spectrum(fs)}
        for r, (i, j, lab) in enumerate(((0, 2, "x-z"), (1, 2, "y-z"), (0, 1, "x-y"))):
            ax = axs[r, c]
            coh, ph = ca.coherence(S, i, j)
            ax.plot(fs, ph, color=colors["truth"], lw=2.0, label="truth")
            for name, Sm in models.items():
                _, phm = ca.coherence(Sm, i, j)
                ax.plot(fs, phm, color=colors[name], lw=1.2, ls="--" if name == "COR" else "-", label=name)
            wgt = (TWO_PI * fs) ** 4 * np.sqrt(np.real(S[:, i, i]) * np.real(S[:, j, j]))
            ax2 = ax.twinx()
            ax2.fill_between(fs, 0, wgt / wgt.max(), color="#bbbbbb", alpha=0.35, lw=0)
            ax2.set_ylim(0, 4); ax2.set_yticks([])
            ax.set_xscale("log"); ax.set_ylim(-180, 180); ax.set_yticks([-180, -90, 0, 90, 180])
            if c == 0:
                ax.set_ylabel(f"{lab} cross-spectral phase (deg)")
            if r == 0:
                ax.set_title(sea.replace("_fp0.29_hs0.8_g3.3", " 0.29 Hz").replace("_", " "))
            if r == 2:
                ax.set_xlabel("frequency (Hz)")
    axs[0, 0].legend(loc="lower left", fontsize=7)
    fig.suptitle("Cross-spectral phase: truth vs oracle-parameterized models (grey: acceleration band weight)", fontsize=10)
    fig.tight_layout()
    save_figure(fig, outdir / f"cross_phase_models.{ext}", dpi=150)
    plt.close(fig)

    fig, axs = plt.subplots(3, len(seas), figsize=(12.5, 7.6), sharex=True)
    for c, sea in enumerate(seas):
        S, op, kappa, sig, Sh, Sig_c, Sig_w, Sig_mca = model_params(sea, fs)
        models = {"COR": sys_ind_cor(kappa, Sig_w, None).spectrum(fs),
                  "SHM": sys_shm(kappa, Sig_c).spectrum(fs),
                  "MCA": sys_mca(kappa, Sig_mca).spectrum(fs)}
        for r, (i, j, lab) in enumerate(((0, 2, "x-z"), (1, 2, "y-z"), (0, 1, "x-y"))):
            ax = axs[r, c]
            coh, _ = ca.coherence(S, i, j)
            ax.plot(fs, coh, color=colors["truth"], lw=2.0, label="truth")
            for name, Sm in models.items():
                cm_, _ = ca.coherence(Sm, i, j)
                ax.plot(fs, cm_, color=colors[name], lw=1.2, ls="--" if name == "COR" else "-", label=name)
            ax.set_xscale("log"); ax.set_ylim(0, 1.05)
            if c == 0:
                ax.set_ylabel(f"{lab} coherence")
            if r == 0:
                ax.set_title(sea.replace("_fp0.29_hs0.8_g3.3", " 0.29 Hz").replace("_", " "))
            if r == 2:
                ax.set_xlabel("frequency (Hz)")
    axs[0, 0].legend(loc="lower left", fontsize=7)
    fig.suptitle("Magnitude-squared coherence: truth vs oracle-parameterized models", fontsize=10)
    fig.tight_layout()
    save_figure(fig, outdir / f"coherence_models.{ext}", dpi=150)
    plt.close(fig)

    # Steady-state RMS, x/y/z, h-noise x1 and x3.
    st = res["steady_state"]
    names = ["IND", "COR", "SHM", "SHM_diag", "MCA", "SHR", "SHM_wrong180", "SHM_wrong90"]
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 6.2))
    seas_all = list(st.keys())
    for r, hn in enumerate(("hnoise_x1", "hnoise_x3")):
        for a in range(3):
            ax = axs[r, a]
            xs = np.arange(len(seas_all))
            wdt = 0.8 / (len(names) + 1)
            for k, nm in enumerate(names + ["wiener"]):
                vals = []
                for sea in seas_all:
                    m = st[sea]["models"][hn]
                    if nm == "wiener":
                        vals.append(1000 * m["wiener"]["rms_joint_m"][a])
                    else:
                        vals.append(1000 * m[nm]["rms_m"][a])
                ax.bar(xs + (k - len(names) / 2) * wdt, vals, wdt, label=nm)
            ax.set_xticks(xs); ax.set_xticklabels([s.split("_")[0][:8] + ("" if "fp" not in s else s[s.find("fp"):s.find("fp") + 6]) for s in seas_all], rotation=30, fontsize=7)
            ax.set_yscale("log")
            ax.set_title(f"{'xyz'[a]} axis, {hn.replace('hnoise_', 'horizontal noise ')}")
            if a == 0:
                ax.set_ylabel("steady-state RMS error (mm)")
    axs[0, 0].legend(fontsize=6, ncol=3)
    fig.tight_layout()
    save_figure(fig, outdir / f"steady_state_rms.{ext}", dpi=150)
    plt.close(fig)

    # Complexification: real-part spectrum and mirror leakage.
    cm = ComplexMode()
    fq = np.geomspace(0.02, 3.0, 400)
    Hp, Hm = np.abs(cm.H(fq)) ** 2, np.abs(cm.H(-fq)) ** 2
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    axs[0].loglog(fq, ca.unit_psd(4, C41J, fq), color="k", lw=2.2, label="C41J |B/P|^2")
    axs[0].loglog(fq, 0.5 * (Hp + Hm), color="#d62728", ls="--", label="(|H+|^2+|H-|^2)/2, complexified")
    axs[0].loglog(fq, 0.5 * Hp, color="#ff9896", lw=1, label="positive-frequency part")
    axs[0].loglog(fq, 0.5 * Hm, color="#9467bd", lw=1, label="mirror part")
    axs[0].set_xlabel("frequency (Hz, kappa = 1)"); axs[0].set_ylabel("unit PSD"); axs[0].legend(fontsize=7)
    axs[1].semilogx(fq, 10 * np.log10(Hm / Hp), color="#9467bd")
    axs[1].set_xlabel("frequency (Hz, kappa = 1)"); axs[1].set_ylabel("mirror leakage |H-|^2/|H+|^2 (dB)")
    fig.tight_layout()
    save_figure(fig, outdir / f"complexification.{ext}", dpi=150)
    plt.close(fig)

    # Complex-OU modal form: kinematics of common vs independent driving,
    # conditioning of the modal transformation, admissible phase law.
    m = modal_c41(1.0)
    A, G, cp = m["A"], m["G"], m["cp"]
    Wc = np.outer(G, G)
    Wi = np.zeros((4, 4)); Wi[:2, :2] = Wc[:2, :2]; Wi[2:, 2:] = Wc[2:, 2:]
    fq = np.geomspace(0.03, 20.0, 400)
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.5))
    for tag, W, col in (("common noise (exact C41J)", Wc, "k"), ("independently driven modes", Wi, "#d62728")):
        Sa = []
        for x in fq:
            row = cp @ np.linalg.solve(1j * TWO_PI * x * np.eye(4) - A, np.eye(4))
            Sa.append(np.real(row @ W @ np.conj(row)) * (TWO_PI * x) ** 4)
        axs[0].loglog(fq, Sa, color=col, lw=1.6, label=tag)
    axs[0].set_xlabel("frequency (Hz, kappa = 1)"); axs[0].set_ylabel("acceleration PSD (unit gain)")
    axs[0].set_title("modal kinematics", fontsize=9); axs[0].legend(fontsize=7)
    eq = res["modal"]["equivalence_vs_kappa"] if res and "modal" in res else []
    if eq:
        kk = [r["kappa"] for r in eq]
        axs[1].semilogy(kk, [r["cond_T"] for r in eq], "o-", color="#8c564b", label="cond(T)")
        axs[1].semilogy(kk, [max(r["tf_rel_err"], r["stationary_rel_err"], r["Q_rel_err"]) * 1e12 for r in eq], "s--",
                        color="#1f77b4", label="worst rel. error x 1e12")
        axs[1].set_xlabel("kappa"); axs[1].legend(fontsize=7)
        axs[1].set_title("modal transformation and equivalence", fontsize=9)
    sea = "jonswap_fp0.29_hs0.8_g3.3"
    fs = np.geomspace(0.08, 1.2, 300)
    S, op, kappa, sig, Sh, Sig_c, Sig_w, Sig_mca = model_params(sea, fs)
    _, ph = ca.coherence(S, 0, 2)
    axs[2].semilogx(fs, ph, color="k", lw=2.0, label="truth")
    for name, sy, col in (("B (complexified, constant phase)", sys_shm(kappa, Sig_c), "#d62728"),
                          ("modal numerator (admissible law)", sys_modn(kappa, Sig_c, op["fz_hz"]), "#ff7f0e"),
                          ("A (real)", sys_ind_cor(kappa, Sig_w, None), "#1f77b4")):
        _, phm = ca.coherence(sy.spectrum(fs), 0, 2)
        axs[2].semilogx(fs, phm, color=col, lw=1.3, ls="--" if name.startswith("A") else "-", label=name)
    axs[2].axvline(op["fz_hz"], color="#999999", lw=0.8)
    axs[2].set_ylim(-180, 180); axs[2].set_yticks([-180, -90, 0, 90, 180])
    axs[2].set_xlabel("frequency (Hz)"); axs[2].set_ylabel("x-z phase (deg)")
    axs[2].set_title("JONSWAP x-z phase, oracle parameters", fontsize=9); axs[2].legend(fontsize=7)
    fig.tight_layout()
    save_figure(fig, outdir / f"modal_decomposition.{ext}", dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "reports/results/carma3_coupled"))
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--quick", action="store_true", help="coarser grids (development only)")
    ap.add_argument("--plots-from-json", metavar="DIR", help="only redraw the figures into DIR from <out>/analysis.json")
    ap.add_argument("--ext", default="svg", choices=["svg", "png", "pdf"])
    args = ap.parse_args(argv)
    out = Path(args.out)
    if args.plots_from_json:
        res = json.loads((out / "analysis.json").read_text())
        make_plots(res, None, Path(args.plots_from_json), args.ext)
        return 0
    out.mkdir(parents=True, exist_ok=True)
    f = np.linspace(0.002, 2.5, 1250 if args.quick else 2500)
    res = {}
    print("complexification", flush=True)
    res["complexification"] = complexification_report(np.geomspace(0.01, 5.0, 2000))
    print("system checks", flush=True)
    res["system_checks"] = system_checks(f)
    seas = ["jonswap_fp0.29_hs0.8_g3.3", "pm_fp0.25_hs1", "jonswap_longcrested", "jonswap_wide_spread",
            "mixed_swell_wind", "bimodal", "crossing"]
    print("modal decomposition", flush=True)
    res["modal"] = modal_report(f)
    print("steady state", flush=True)
    res["steady_state"], curves = steady_study(f, seas)
    print("identifiability", flush=True)
    res["identifiability"] = identifiability(np.linspace(0.01, 2.0, 800 if args.quick else 1600))
    res["cost"] = {
        "IND/COR (12 motion, blocks 4,4,4)": cost(12, [4, 4, 4], 0),
        "MOD (12 motion, modal blocks 4,4,4, dense readouts)": cost(12, [4, 4, 4], 9),
        "SHM (24 motion, blocks 8,8,8)": cost(24, [8, 8, 8], 0),
        "MCA (26 motion, blocks 8,8,10)": cost(26, [8, 8, 10], 2),
        "SHR (16 motion, blocks 8,4,4)": cost(16, [8, 4, 4], 4),
        "K2 (48 motion)": cost(48, [8] * 6, 3),
    }
    (out / "analysis.json").write_text(json.dumps(res, indent=1))
    if not args.no_plots:
        make_plots(res, f, out)
    summarize(res)
    return 0


def summarize(res):
    c = res["complexification"]
    print(f"complexification: realpart PSD max rel err {c['max_rel_error_realpart_psd_vs_C41J']:.2e}, "
          f"stationary cov dev {c['realpart_stationary_cov_pva_max_rel_dev_vs_C41J']:.2e}, leakage {c['mirror_leakage_Hminus_over_Hplus']}")
    md = res["modal"]
    print("modal:", {k: md[k] for k in ("poles_kappa1", "residues_kappa1", "sum_re_r", "sum_re_rp", "sum_im_r", "kinematic_cG_cAG")})
    print("   worst equivalence:", {k: max(r[k] for r in md["equivalence_vs_kappa"]) for k in md["equivalence_vs_kappa"][0] if k != "kappa"})
    print("   rotated:", md["rotated_readout"])
    print("   independent:", md["independent_driving"])
    sc = res["system_checks"]
    print("system checks:", {k: v for k, v in sc.items() if not isinstance(v, dict)})
    for sea, r in res["steady_state"].items():
        m = r["models"]["hnoise_x1"]
        line = " ".join(f"{k}={1000 * np.hypot(v['rms_m'][0], v['rms_m'][1]):.1f}/{1000 * v['rms_m'][2]:.1f}"
                        for k, v in m.items() if k != "wiener")
        print(f"{sea}: h/z mm {line}")
        m3 = r["models"]["hnoise_x3"]
        line3 = " ".join(f"{k}={1000 * np.hypot(v['rms_m'][0], v['rms_m'][1]):.1f}/{1000 * v['rms_m'][2]:.1f}"
                         for k, v in m3.items() if k != "wiener")
        print(f"   x3: {line3}")
        print("   wiener joint x1 (mm):", [round(1000 * x, 1) for x in m["wiener"]["rms_joint_m"]],
              "indep:", [round(1000 * x, 1) for x in m["wiener"]["rms_independent_m"]])
        print("   var ratio actual/pred x1:", {k: [round(x, 2) for x in v["var_ratio"]] for k, v in m.items() if k in ("IND", "COR", "SHM", "MCA", "SHM_wrong180")})
        print("   representable:", {k: (round(v['real_part_energy_fraction'], 3), round(v['constant_phase_fit_energy_fraction'], 3),
                                        round(v['phase_deg_at_peak'], 1)) for k, v in r["representable"].items()})
    for k, v in res["identifiability"].items():
        if isinstance(v, dict):
            print(f"identifiability {k}: n={v['n_params']} rank={v['numerical_rank']} eig {v['fim_eig_min']:.2e}..{v['fim_eig_max']:.2e}")
            if v["n_params"] < 20:
                print("    crlb rel:", {l: (round(x, 3) if x is not None else None) for l, x in zip(v["labels"], v["crlb_rel_std"])})


if __name__ == "__main__":
    raise SystemExit(main())
