#!/usr/bin/env python3
"""Covariance/transition primitives used by the active OU-III P4 proof route.

These routines are the small reusable subset formerly embedded in experimental
P5 prefix producers.  Keeping them in a P4 utility module lets obsolete P5
search/certificate drivers be removed without changing the arithmetic used by
the current joint-Joseph/augmented complete-word construction.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from ou3_interval import Interval, matrix_identity, matrix_mul, matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_full_process_ucc as PROCESS
import ou3_p4_golive_covariance as GOLIVE
import ou3_source_reachable_matrix_p3 as P3CELL
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"


def I(x: float) -> Interval:
    return Interval.point(float(x))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def zero(rows: int, cols: int):
    return [[I(0.0) for _ in range(cols)] for _ in range(rows)]


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if lo > hi:
        raise RuntimeError(f"rigorous interval enclosures became disjoint: {a} vs {b}")
    return Interval(lo, hi)


def psd_tighten(Pm):
    """Tighten an entrywise covariance box using PSD/Cauchy-Schwarz."""
    Pm = matrix_symmetric_hull(Pm)
    n = len(Pm)
    dhi = []
    for i in range(n):
        if not math.isfinite(Pm[i][i].hi) or Pm[i][i].hi < 0.0:
            raise RuntimeError("covariance diagonal upper is invalid")
        dhi.append(up(max(0.0, Pm[i][i].hi)))
    out = [[Pm[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        out[i][i] = _intersect(out[i][i], Interval(0.0, dhi[i]))
        for j in range(i + 1, n):
            b = up(math.sqrt(up(dhi[i] * dhi[j])))
            z = _intersect(
                Interval(min(out[i][j].lo, out[j][i].lo), max(out[i][j].hi, out[j][i].hi)),
                Interval(-b, b),
            )
            out[i][j] = z
            out[j][i] = z
    return out


def reset_covariance(Pm, dx_theta):
    if len(dx_theta) != 3:
        raise ValueError("reset correction must be length three")
    n = len(Pm)
    x, y, z = dx_theta
    G = matrix_identity(n)
    h = I(0.5)
    G[0][1] = -h * z; G[0][2] = h * y
    G[1][0] = h * z;  G[1][2] = -h * x
    G[2][0] = -h * y; G[2][1] = h * x
    return psd_tighten(matrix_mul(matrix_mul(G, Pm), matrix_transpose(G)))


def R_diag(std: float):
    r = Interval.outward_bounds(float(std) ** 2, float(std) ** 2)
    out = zero(3, 3)
    for i in range(3):
        out[i][i] = r
    return out


def R_S(src: dict):
    r = src["R_S_filter_std"].square()
    factors = src.get("R_S_axis_std_factors")
    if factors is None:
        factors = P3CELL.source_rs_axis_std_factors()
    out = zero(3, 3)
    for i in range(3):
        out[i][i] = r * I(factors[i]).square()
    return out


def source_pb0() -> float:
    text = MEKF.read_text(encoding="utf-8")
    m = re.search(r"T\s+Pq0\s*=\s*T\([^)]*\),\s*T\s+Pb0\s*=\s*T\(([0-9.eE+-]+)\)", text)
    if not m:
        raise RuntimeError("cannot extract gyro-bias covariance seed")
    return float(m.group(1))


def startup_timeout_s() -> float:
    text = WRAPPER.read_text(encoding="utf-8")
    m = re.search(r"proxy_startup_timeout_sec\s*=\s*([0-9.eE+-]+)f", text)
    if not m:
        raise RuntimeError("cannot extract source startup timeout")
    t = float(m.group(1))
    if not t > 0.0:
        raise RuntimeError("invalid source startup timeout")
    return t


def rotation_step_box(rate_rad_s: float, h: float):
    theta = up(float(rate_rad_s) * float(h))
    s = VT.sin_point(theta)
    co = VT.cos_point(theta)
    sin_hi = max(abs(s.lo), abs(s.hi))
    one_minus_c = up(1.0 - co.lo)
    off = up(sin_hi + 0.5 * one_minus_c)
    R = zero(3, 3)
    B = zero(3, 3)
    for i in range(3):
        R[i][i] = Interval(co.lo, 1.0)
        B[i][i] = Interval(down(h * co.lo), up(h))
        for j in range(3):
            if i != j:
                R[i][j] = Interval(-off, off)
                B[i][j] = Interval(down(-h * off), up(h * off))
    return R, B


def _shipping_coeffs_point(tau: float, h: float):
    tau = float(tau); h = float(h)
    if not (tau > 0.0 and h > 0.0):
        raise ValueError("positive tau/h required")
    x0 = h / tau
    x = Interval.outward_bounds(x0, x0)
    tt = Interval.outward_bounds(tau, tau)
    em1 = VT.expm1_interval(-x)
    alpha = VT.exp_interval(-x)
    phi_va = -tt * em1
    if abs(x0) < 1.0e-2:
        x2 = x.square(); x3 = x2 * x; x4 = x3 * x; x5 = x4 * x
        phi_pa = tt.square() * (I(0.5) * x2 - I(1.0 / 6.0) * x3 + I(1.0 / 24.0) * x4)
        phi_Sa = (tt.square() * tt) * (
            I(1.0 / 6.0) * x3 - I(1.0 / 24.0) * x4 + I(1.0 / 120.0) * x5
        )
    else:
        phi_pa = tt.square() * (x + em1)
        phi_Sa = (tt.square() * tt) * (I(0.5) * x.square() - x - em1)
    return alpha, phi_va, phi_pa, phi_Sa


def _monotone_coeff_hull(tau: Interval, h: float):
    lo = _shipping_coeffs_point(tau.lo, h)
    hi = _shipping_coeffs_point(tau.hi, h)
    return tuple(Interval(a.lo, b.hi) for a, b in zip(lo, hi))


def _ou_process_moment_axis(tau: Interval, sigma: Interval, h: float):
    qc = up(2.0 * sigma.hi * sigma.hi / tau.lo)
    powers = (1, 2, 3, 0)
    facts = (1.0, 2.0, 6.0, 1.0)
    Q = zero(4, 4)
    for i in range(4):
        for j in range(i, 4):
            p = powers[i] + powers[j]
            den = facts[i] * facts[j] * (p + 1.0)
            ub = up(qc * (h ** (p + 1)) / den)
            z = Interval(0.0, ub)
            Q[i][j] = z
            Q[j][i] = z
    return Q


def tight_transition_and_Q(n: int, src: dict, domain: dict):
    """Dependency-preserving transition/process enclosure for H/A dimensions."""
    if n not in (18, 21):
        raise ValueError("P4 state dimension must be 18 or 21")
    h = float(src["dt_s"])
    tau = src["tau_s"]
    sigma = src["sigma_aw_mps2"]
    alpha, phi_va, phi_pa, phi_Sa = _monotone_coeff_hull(tau, h)

    F = matrix_identity(n)
    rate = float(domain["normal_live"]["body_rate_norm_upper_deg_s"]) * math.pi / 180.0
    Rstep, Bstep = rotation_step_box(rate, h)
    for i in range(3):
        for j in range(3):
            F[i][j] = Rstep[i][j]
            F[i][3 + j] = Bstep[i][j]
    for ax in range(3):
        iv, ip, iS, ia = 6 + ax, 9 + ax, 12 + ax, 15 + ax
        F[iv][iv] = I(1.0); F[iv][ia] = phi_va
        F[ip][iv] = I(h); F[ip][ip] = I(1.0); F[ip][ia] = phi_pa
        F[iS][iv] = I(0.5 * h * h); F[iS][ip] = I(h); F[iS][iS] = I(1.0); F[iS][ia] = phi_Sa
        F[ia][ia] = Interval(alpha.lo, min(1.0, alpha.hi))

    Q = zero(n, n)
    proc = PROCESS.build()["source_constants"]
    qg = float(proc["gyro_noise_density_rad_sqrt_s"]) ** 2
    qb = float(proc["gyro_bias_rw_variance_density"])
    bb = up(qb * h ** 3 / 3.0)
    cross = up(qb * h * h / 2.0)
    for i in range(3):
        Q[i][i] = Interval(down(qg * h), up(qg * h + bb))
        Q[3 + i][3 + i] = Interval.outward_bounds(qb * h, qb * h)
        for j in range(3):
            if i != j:
                Q[i][j] = Interval(-bb, bb)
            Q[i][3 + j] = Interval(-cross, cross)
            Q[3 + j][i] = Q[i][3 + j]

    qaxis = _ou_process_moment_axis(tau, sigma, h)
    groups = (6, 9, 12, 15)
    for ax in range(3):
        ids = [g + ax for g in groups]
        for i in range(4):
            for j in range(4):
                Q[ids[i]][ids[j]] = qaxis[i][j]
    return F, psd_tighten(Q), Rstep


def initial_covariance(n: int, src: dict, domain_path: Path):
    """Source goLive covariance enclosure used before P4 word propagation."""
    if n not in (18, 21):
        raise ValueError("P4 state dimension must be 18 or 21")
    seed = GOLIVE.build(domain_path)["goLive_H_covariance_seed"]
    Pm = zero(n, n)
    a = seed["attitude_covariance_seed"]
    tilt = float(a["tilt_variance"])
    yaw = float(a["gauged_yaw_variance"])
    for i in range(3):
        Pm[i][i] = Interval(down(tilt), up(yaw))
        for j in range(i + 1, 3):
            b = up(0.5 * (yaw - tilt))
            Pm[i][j] = Interval(-b, b)
            Pm[j][i] = Pm[i][j]

    proc = PROCESS.build()["source_constants"]
    qb = float(proc["gyro_bias_rw_variance_density"])
    pbg_hi = up(source_pb0() + qb * startup_timeout_s())
    for i in range(3, 6):
        Pm[i][i] = Interval(0.0, pbg_hi)
        for j in range(3, 6):
            if i != j:
                Pm[i][j] = Interval(-pbg_hi, pbg_hi)

    for i in range(6, 9):
        v = float(seed["P_vv_variance_per_axis"])
        Pm[i][i] = Interval.outward_bounds(v, v)
    for i in range(9, 12):
        v = float(seed["P_pp_variance_per_axis"])
        Pm[i][i] = Interval.outward_bounds(v, v)
    for i in range(12, 15):
        v = float(seed["P_SS_variance_per_axis"])
        Pm[i][i] = Interval.outward_bounds(v, v)
    awv = src["sigma_aw_mps2"].square()
    for i in range(15, 18):
        Pm[i][i] = awv
    return psd_tighten(Pm)
