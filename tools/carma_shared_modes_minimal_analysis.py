#!/usr/bin/env python3
"""Mathematical analysis of the minimal realization of the shared-mode CARMA block.

Companion of docs/kalman-carma-shared-modes-minimal.md and
doc/kalman_carma_shared_modes/kalman_carma_shared_modes-minimal.tex
(section 2-4).  Double precision, independent of the C++ code (the C++
identities are tested in tests/kalman_carma_shared_modes/carma_shared_modes_minimal-test.cpp).

  1. observability of the complexified C41J mode per axis: rank 6 of 8, the
     unobservable subspace lies in the imaginary half and equals
     {0} x span(Re w, Im w), w = [B(r), r B(r), r^2 B(r), 1], D(r) = 0
  2. its motion with kappa (principal angles between U(kappa) and U(1))
  3. cascade realization C41J * V, V = [N_r, N_i] / D co-inner: transfer
     functions, lossless identities, kappa covariance
  4. minimal stochastic order of the three-axis law: Hankel singular values
     of the positive-real part, n = 4 rank(Sigma) + 2 rank(Im Sigma)
  5. conditioning (stationary correlation matrices) of each realization
  6. motion-only Kalman experiment: fixed kappa (all realizations equal) and
     coupling off under adaptation (cascade = A41 exactly; 33/quotient differ)

Writes reports/results/kalman_carma_shared_modes_minimal/analysis.json and
figures/analysis_*.{svg,pgf}.

Run: python3 tools/carma_shared_modes_minimal_analysis.py [--ext svg,pgf]
"""
from __future__ import annotations

import sys as _export_sys
from pathlib import Path as _ExportPath
_export_sys.path.insert(0, str(_ExportPath(__file__).resolve().parents[1] / "plots"))
from vector_export import save_figure


import argparse
import json
from pathlib import Path

import numpy as np
import scipy.linalg as sl

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "results" / "kalman_carma_shared_modes_minimal"
FIG = OUT / "figures"
TWO_PI = 2.0 * np.pi
# CARMA41_JONSWAP_SHAPE (src/kalman_carma/CarmaVerticalModel.h).
SHAPE = dict(f1=0.306945, z1=0.104191, f2=0.815318, z2=0.478188, fz=0.0)
J8 = np.block([[np.zeros((4, 4)), -np.eye(4)], [np.eye(4), np.zeros((4, 4))]])
J2 = np.array([[0.0, -1.0], [1.0, 0.0]])


# ---------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------

def base_A(k):
    """C41J in [p, v, a, l] (carma_dynamics<4>)."""
    w1, z1 = TWO_PI * k * SHAPE["f1"], SHAPE["z1"]
    w2, z2 = TWO_PI * k * SHAPE["f2"], SHAPE["z2"]
    wz0 = TWO_PI * k * SHAPE["fz"]
    d = wz0 - 2 * z2 * w2
    A = np.zeros((4, 4))
    A[0, 1] = A[1, 2] = 1
    A[2] = [d * w1 * w1, 2 * d * z1 * w1 - w1 * w1, d - 2 * z1 * w1, -w2 * w2 - d * wz0]
    A[3] = [w1 * w1, 2 * z1 * w1, 1, -wz0]
    return A


def polys(k):
    """P+ = s^2 + b1 s + b0 (upper-half-plane poles of P), D = s^2 + d1 s + d0."""
    w1, z1 = TWO_PI * k * SHAPE["f1"], SHAPE["z1"]
    w2, z2 = TWO_PI * k * SHAPE["f2"], SHAPE["z2"]
    p1 = complex(-z1 * w1, w1 * np.sqrt(1 - z1 * z1))
    p2 = complex(-z2 * w2, w2 * np.sqrt(1 - z2 * z2))
    b1, b0 = -(p1 + p2), p1 * p2
    pp = lambda w: (1j * w) ** 2 + b1 * 1j * w + b0  # noqa: E731
    q = lambda w: abs(pp(w)) ** 2 + abs(pp(-w)) ** 2  # noqa: E731
    gamma = 0.5 * q(0.0)
    wr = max(w1, w2)
    beta = (0.5 * q(wr) - wr ** 4 - gamma) / wr ** 2
    d0 = np.sqrt(gamma)
    d1 = np.sqrt(beta + 2 * d0)
    return b1, b0, d1, d0


def complexified_A8(k):
    """The 33-state block's mode (coupled_detail::complexified_c41), realified."""
    b1, b0, d1, d0 = polys(k)
    wz0 = TWO_PI * k * SHAPE["fz"]
    c3, c2, c1, c0 = b1 + d1, b0 + b1 * d1 + d0, b0 * d1 + b1 * d0, b0 * d0
    Ac = np.zeros((4, 4), complex)
    Ac[0, 1] = Ac[1, 2] = Ac[2, 3] = 1
    Ac[3] = [-c0, -c1, -c2, -c3]
    cv = np.array([wz0, 1, 0, 0], complex)
    T = np.vstack([cv, cv @ Ac, cv @ Ac @ Ac, [1, 0, 0, 0]])
    Az = T @ Ac @ np.linalg.inv(T)
    return np.block([[Az.real, -Az.imag], [Az.imag, Az.real]])


def unit_g2(k):
    W = np.zeros((4, 4))
    W[2, 2] = 1
    return 1.0 / sl.solve_continuous_lyapunov(base_A(k), -W)[2, 2]


def coinner(k):
    """Input-normal realization of V = [N_r, N_i] / D (C++: coinner_input_normal)."""
    b1, b0, d1, d0 = polys(k)
    A0 = np.array([[-d1, 1.0], [-d0, 0.0]])
    B0 = np.array([[b1.real - d1, b1.imag], [b0.real - d0, b0.imag]])
    X = sl.solve_continuous_lyapunov(A0, -B0 @ B0.T)
    L = np.linalg.cholesky(X)
    Li = np.linalg.inv(L)
    A, B, C = Li @ A0 @ L, Li @ B0, np.array([[1.0, 0.0]]) @ L
    Q, R = np.linalg.qr(B)
    O = np.diag(np.sign(np.diag(R))) @ Q.T
    return O @ A @ O.T, O @ B, C @ O.T


def cascade6(k):
    Aw, Bw, Cw = coinner(k)
    A = np.zeros((6, 6))
    A[:4, :4] = base_A(k)
    A[2, 4:] = Cw[0]
    A[4:, 4:] = Aw
    G = np.zeros((6, 2))
    G[2, 0] = 1
    G[4:] = Bw
    return A, G


def quotient_N(k):
    """Rows of D(s)(alpha + beta s) on [Im zeta, Im zeta', Im zeta'', Im l]."""
    _, _, d1, d0 = polys(k)
    wz0 = TWO_PI * k * SHAPE["fz"]
    rows = []
    for al, be in ((1.0, 0.0), (0.0, 1.0)):
        n2 = be
        n1 = al + be * d1 - be * wz0
        n0 = al * d1 + be * d0 - n1 * wz0
        n3 = al * d0 - n0 * wz0
        rows.append([n0, n1, n2, n3])
    return np.array(rows)


def quotient6(k):
    A8 = complexified_A8(k)
    Pi = np.zeros((6, 8))
    Pi[:4, :4] = np.eye(4)
    Pi[4:, 4:] = quotient_N(k)
    A6 = Pi @ A8 @ np.linalg.pinv(Pi)
    G8 = np.zeros((8, 2))
    G8[2, 0] = G8[6, 1] = 1
    return A6, Pi @ G8, Pi, A8


def obsv(A, C):
    return np.vstack([C @ np.linalg.matrix_power(A, i) for i in range(A.shape[0])])


# ---------------------------------------------------------------------------
# 1-5: structure
# ---------------------------------------------------------------------------

def unobservable(k):
    A8 = complexified_A8(k)
    ca = np.zeros((1, 8))
    ca[0, 2] = 1
    # scale-aware rank: normalize the observability rows
    O = obsv(A8, ca)
    O = O / np.linalg.norm(O, axis=1, keepdims=True)
    s = np.linalg.svd(O, compute_uv=False)
    U = sl.null_space(O, rcond=1e-10)
    return U, s


def analytic_U(k):
    """{0} x span(Re w, Im w), w = [B(r), r B(r), r^2 B(r), 1] at the roots of D."""
    _, _, d1, d0 = polys(k)
    wz0 = TWO_PI * k * SHAPE["fz"]
    r = np.roots([1.0, d1, d0])
    vecs = []
    for rr in r:
        B = rr + wz0
        w = np.array([B, rr * B, rr * rr * B, 1.0])
        vecs += [np.real(w), np.imag(w)]
    S = sl.orth(np.array(vecs).T)
    return np.vstack([np.zeros((4, S.shape[1])), S])


def tf(A, G, c, w):
    return c @ np.linalg.solve(1j * w * np.eye(A.shape[0]) - A, G)


def noise_K(S):
    """Per-axis-pair 2 x 2 noise blocks [[Re, -Im], [Im, Re]] (the C++ convention)."""
    K = np.zeros((6, 6))
    for a in range(3):
        for b in range(3):
            K[2 * a:2 * a + 2, 2 * b:2 * b + 2] = [[S[a, b].real, -S[a, b].imag], [S[a, b].imag, S[a, b].real]]
    return K


def hankel_order(A1, G1, c, S, tol=1e-7):
    A = np.kron(np.eye(3), A1)
    G = np.kron(np.eye(3), G1)
    P = sl.solve_continuous_lyapunov(A, -(G @ noise_K(S) @ G.T))
    P = 0.5 * (P + P.T)
    C = np.kron(np.eye(3), c[None, :])
    B = P @ C.T                                          # positive-real part Z(s) = C (sI - A)^-1 P C^T
    Wc = sl.solve_continuous_lyapunov(A, -B @ B.T)
    Wo = sl.solve_continuous_lyapunov(A.T, -C.T @ C)
    h = np.sort(np.sqrt(np.abs(np.linalg.eigvals(Wc @ Wo))))[::-1]
    return int((h > h[0] * tol).sum()), h / h[0]


def corr_condition(P):
    d = np.sqrt(np.diag(P))
    e = np.linalg.eigvalsh(P / np.outer(d, d))
    return float(e.max() / e.min())


# ---------------------------------------------------------------------------
# 6: motion-only Kalman experiment
# ---------------------------------------------------------------------------

def disc(A, W, h):
    n = A.shape[0]
    M = np.zeros((2 * n, 2 * n))
    M[:n, :n] = -A
    M[:n, n:] = W
    M[n:, n:] = A.T
    E = sl.expm(M * h)
    F = E[n:, n:].T
    return F, F @ E[:n, n:]


class Motion:
    """Motion-only filter on y = a + n (three axes), realizations of the shared mode."""

    def __init__(self, kind, h):
        self.kind, self.h, self.k, self.scale = kind, h, None, None

    def build(self, k, sig, S):
        g2 = unit_g2(k)
        if self.kind == "33":
            A = complexified_A8(k)
            G = np.zeros((8, 2))
            G[2, 0] = G[6, 1] = 1
        elif self.kind in ("27Q", "27H"):
            A, G, self.Pi, A8 = quotient6(k)
            W8 = np.zeros((8, 8))
            W8[2, 2] = W8[6, 6] = 1
            P8 = sl.solve_continuous_lyapunov(A8, -W8)
            self.L = P8 @ self.Pi.T @ np.linalg.inv(self.Pi @ P8 @ self.Pi.T)
        elif self.kind == "27C":
            A, G = cascade6(k)
            self.scale = np.sqrt(g2) * sig
        else:   # A41
            A = base_A(k)
            G = np.zeros((4, 1))
            G[2, 0] = 1
        n = A.shape[0]
        if self.kind == "A41":
            F, QR = disc(A, g2 * G @ G.T, self.h)
            QI = 0 * QR
        else:
            F, QR = disc(A, g2 * G @ G.T, self.h)
            _, QI = disc(A, g2 * G @ J2 @ G.T, self.h)
        self.nb, self.k = n, k
        self.F = np.kron(np.eye(3), F)
        self.Q = np.zeros((3 * n, 3 * n))
        for a in range(3):
            for b in range(3):
                Sab = S[a, b] if self.kind != "A41" else (S[a, b].real if a == b else 0.0)
                self.Q[a * n:(a + 1) * n, b * n:(b + 1) * n] = np.real(Sab) * QR + np.imag(Sab) * QI
        self.Q = 0.5 * (self.Q + self.Q.T)
        self.H = np.zeros((3, 3 * n))
        for a in range(3):
            self.H[a, a * n + 2] = 1


class KF:
    def __init__(self, kind, h, k, sig, S, R):
        self.m = Motion(kind, h)
        self.m.build(k, sig, S)
        self.x = np.zeros(3 * self.m.nb)
        self.P = sl.solve_discrete_lyapunov(self.m.F, self.m.Q)
        self.R = R

    def commit(self, k, sig, S):
        old, new = self.m, Motion(self.m.kind, self.m.h)
        new.build(k, sig, S)
        M = None
        if new.kind == "27Q" and k != old.k:
            M = np.kron(np.eye(3), new.Pi @ old.L)
        elif new.kind == "27C":
            d = np.ones(18)
            for a in range(3):
                d[a * 6 + 4:a * 6 + 6] = new.scale[a] / old.scale[a]
            M = np.diag(d)
        if M is not None:
            self.x = M @ self.x
            self.P = M @ self.P @ M.T
        self.m = new

    def step(self, y):
        F, Q, H = self.m.F, self.m.Q, self.m.H
        self.x = F @ self.x
        self.P = F @ self.P @ F.T + Q
        Sv = H @ self.P @ H.T + self.R
        K = self.P @ H.T @ np.linalg.inv(Sv)
        self.x = self.x + K @ (y - H @ self.x)
        IKH = np.eye(len(self.x)) - K @ H
        self.P = IKH @ self.P @ IKH.T + K @ self.R @ K.T

    def pos(self):
        return np.array([self.x[a * self.m.nb] for a in range(3)])


def kalman_experiment(seed=11, T=6000, h=0.02):
    rng = np.random.default_rng(seed)
    nf = 40
    f0 = rng.uniform(0.22, 0.4, nf)
    ph = rng.uniform(0, 2 * np.pi, nf)
    amp = rng.rayleigh(0.06, nf)
    th = rng.uniform(-0.4, 0.4, nf) + 0.5

    def run(S0, kinds, adapt):
        fs = {kd: KF(kd, h, 1.0, np.sqrt(np.real(np.diag(S0))), S0, np.eye(3) * 0.03 ** 2) for kd in kinds}
        ktr, amp_s = 1.0, 1.0
        dmax = {kd: 0.0 for kd in kinds}
        noise = np.random.default_rng(seed + 1)
        for i in range(T):
            t = i * h
            if adapt and i % 5 == 0 and i > 0:
                ktr *= np.exp(0.01 * noise.normal())
                amp_s *= np.exp(0.05 * noise.normal())
                S = S0 * amp_s
                for f in fs.values():
                    f.commit(ktr, np.sqrt(np.real(np.diag(S))), S)
            w = TWO_PI * f0
            z, x = amp * np.sin(w * t + ph), amp * np.cos(w * t + ph)
            a = -np.array([np.sum(w * w * x * np.cos(th)), np.sum(w * w * x * np.sin(th)), np.sum(w * w * z)])
            y = a + noise.normal(size=3) * 0.03
            for f in fs.values():
                f.step(y)
            ref = fs[kinds[0]].pos()
            for kd in kinds:
                dmax[kd] = max(dmax[kd], float(np.abs(fs[kd].pos() - ref).max()))
        return dmax

    Sh = np.array([[0.30, 0.05, 0.25j], [0.05, 0.12, 0.10j], [-0.25j, -0.10j, 0.4]])
    fixed = run(Sh, ["33", "27Q", "27H", "27C"], adapt=False)
    off = run(np.diag([0.3, 0.12, 0.4]).astype(complex), ["A41", "33", "27Q", "27C"], adapt=True)
    return {"fixed_kappa_hermitian_max_abs_dp_vs_33_m": fixed, "coupling_off_adaptation_max_abs_dp_vs_A41_m": off}


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ext", default="svg,pgf")
    ap.add_argument("--no-kalman", action="store_true", help="skip the (slow) motion-only Kalman experiment")
    ap.add_argument("--figures-only", action="store_true", help="draw the figures from the committed analysis.json")
    ap.add_argument("--outdir", default=None, help="figure directory (default: the results' figures/)")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir) if args.outdir else FIG
    if args.figures_only:
        figures(json.loads((OUT / "analysis.json").read_text()), args.ext.split(","), outdir)
        return
    res = {}
    # 1. observability
    U1, s = unobservable(1.0)
    Ua = analytic_U(1.0)
    res["observability_singular_values_normalized"] = s.tolist()
    res["unobservable_dim"] = int(U1.shape[1])
    res["unobservable_real_half_norm"] = float(np.abs(U1[:4]).max())
    res["unobservable_vs_analytic_angle_deg"] = float(np.degrees(sl.subspace_angles(U1, Ua)).max())
    b1, b0, d1, d0 = polys(1.0)
    res["D_roots_kappa1"] = [complex(r).real for r in np.roots([1, d1, d0])]
    res["A8_eigenvalues_kappa1"] = [[float(e.real), float(e.imag)] for e in np.sort_complex(np.linalg.eigvals(complexified_A8(1.0)))]
    # 2. motion with kappa
    kap = np.exp(np.linspace(np.log(0.5), np.log(2.0), 61))
    ang = np.array([np.degrees(sl.subspace_angles(U1, analytic_U(k))) for k in kap])
    res["unobservable_angle_vs_kappa"] = {"kappa": kap.tolist(), "angle1_deg": ang[:, 0].tolist(),
                                          "angle2_deg": ang[:, 1].tolist()}
    res["unobservable_angle_per_percent_kappa_deg"] = float(np.degrees(sl.subspace_angles(U1, analytic_U(1.01))).max())
    # leakage gain ||Pi(k) restricted to U(1)|| (normalized quotient rows)
    leak = []
    for k in kap:
        _, _, Pi, _ = quotient6(k)
        N = Pi[4:, 4:]
        N = N / np.linalg.norm(N, axis=1, keepdims=True)
        leak.append(float(np.linalg.norm(N @ U1[4:], 2)))
    res["leak_gain_vs_kappa"] = leak
    # 3. cascade
    A8 = complexified_A8(1.0)
    G8 = np.zeros((8, 2))
    G8[2, 0] = G8[6, 1] = 1
    A6, G6 = cascade6(1.0)
    Aq, Gq, _, _ = quotient6(1.0)
    errs = {}
    for name, (A, G) in {"cascade": (A6, G6), "quotient": (Aq, Gq)}.items():
        e = 0.0
        for idx in range(3):
            c8 = np.eye(8)[idx]
            c6 = np.eye(6)[idx]
            for w in np.linspace(0.05, 30, 300):
                h8 = tf(A8, G8, c8, w)
                e = max(e, float(np.abs(h8 - tf(A, G, c6, w)).max() / np.abs(h8).max()))
        errs[name] = e
    res["tf_rel_error_vs_33"] = errs
    Aw, Bw, Cw = coinner(1.0)
    res["coinner_lossless_residual"] = float(np.abs(Aw + Aw.T + Bw @ Bw.T).max() + np.abs(Cw + Bw[:, 0]).max())
    res["coinner_kappa_covariance_residual"] = float(max(np.abs(coinner(k)[0] - k * Aw).max() for k in (0.5, 2.0)))
    # 4. stochastic order
    rng = np.random.default_rng(1)
    X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    Sh = X @ X.conj().T
    Sw = np.array([[1, 0.3, 0.6j], [0.3, 0.5, 0.2j], [-0.6j, -0.2j, 1.0]])
    v = np.array([0.8, 0.5, 1j])
    cases = {"Hermitian, full rank": Sh, "wave-like (quadrature x-z, y-z)": Sw, "real (in-phase coupling)": Sh.real.astype(complex),
             "diagonal (coupling off)": np.diag([1.0, 2.0, 0.5]).astype(complex),
             "rank 1 (one long-crested wave)": np.outer(v, v.conj())}
    cp6 = np.eye(6)[0]
    orders = {}
    for name, S in cases.items():
        n, h = hankel_order(A6, G6, cp6, S)
        orders[name] = {"rank_Sigma": int(np.linalg.matrix_rank(S, 1e-9)), "rank_Im_Sigma": int(np.linalg.matrix_rank(S.imag, 1e-9)),
                        "predicted": int(4 * np.linalg.matrix_rank(S, 1e-9) + 2 * np.linalg.matrix_rank(S.imag, 1e-9)),
                        "hankel_order": n, "hankel_sv": h[:18].tolist()}
    res["stochastic_minimal_order"] = orders
    # 5. conditioning
    g2 = unit_g2(1.0)
    W8 = np.zeros((8, 8))
    W8[2, 2] = W8[6, 6] = g2
    P8 = sl.solve_continuous_lyapunov(A8, -W8)
    _, _, Pi, _ = quotient6(1.0)
    PQ = Pi @ P8 @ Pi.T
    PC = sl.solve_continuous_lyapunov(A6, -g2 * G6 @ G6.T)
    W4 = np.zeros((4, 4))
    W4[2, 2] = g2
    P4 = sl.solve_continuous_lyapunov(base_A(1.0), -W4)
    res["stationary_correlation_condition"] = {"33 (complexified, 8)": corr_condition(P8), "27 quotient (6)": corr_condition(PQ),
                                               "27 cascade (6)": corr_condition(PC), "A41 (4)": corr_condition(P4)}
    # 6. Kalman experiment
    if not args.no_kalman:
        res["kalman"] = kalman_experiment()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "analysis.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k not in ("unobservable_angle_vs_kappa", "leak_gain_vs_kappa")}, indent=1)[:4000])
    figures(res, args.ext.split(","), outdir)


def figures(res, exts, outdir=FIG):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if "pgf" in exts:
        plt.rcParams.update({"pgf.texsystem": "lualatex", "pgf.rcfonts": False, "font.family": "serif"})
    outdir.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        for e in exts:
            save_figure(fig, outdir / f"{name}.{e}", bbox_inches="tight")
        plt.close(fig)

    d = res["unobservable_angle_vs_kappa"]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax[0].semilogx(d["kappa"], d["angle1_deg"], color="#d62728", label="largest principal angle")
    ax[0].semilogx(d["kappa"], d["angle2_deg"], color="#1f77b4", label="smallest")
    ax[0].set_xlabel(r"frequency scale $\kappa$")
    ax[0].set_ylabel(r"angle $U(\kappa)$ vs $U(1)$ (deg)")
    ax[0].legend(fontsize=7, frameon=False)
    ax[0].grid(alpha=0.3)
    ax[1].semilogx(d["kappa"], res["leak_gain_vs_kappa"], color="#2ca02c")
    ax[1].set_xlabel(r"frequency scale $\kappa$")
    ax[1].set_ylabel(r"$\|N(\kappa)\, U(1)\|_2$")
    ax[1].grid(alpha=0.3)
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    for a in ax:
        a.xaxis.set_major_locator(FixedLocator([0.5, 0.75, 1.0, 1.5, 2.0]))
        a.xaxis.set_minor_locator(NullLocator())
        a.xaxis.set_major_formatter(ScalarFormatter())
    fig.tight_layout()
    save(fig, "analysis_unobservable_subspace")

    o = res["stochastic_minimal_order"]
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    for name, v in o.items():
        h = np.array(v["hankel_sv"])
        ax.semilogy(np.arange(1, len(h) + 1), np.maximum(h, 1e-17), marker="o", ms=3,
                    label=f"{name}: {v['hankel_order']} (4 rk + 2 rk Im = {v['predicted']})")
    ax.axhline(1e-7, color="k", lw=0.6, ls=":")
    ax.set_xlabel("index")
    ax.set_ylabel("normalized Hankel s.v.")
    ax.set_xticks(range(1, 19))
    ax.legend(fontsize=6.5, frameon=False, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "analysis_stochastic_order")


if __name__ == "__main__":
    main()
