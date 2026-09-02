#!/usr/bin/env python3
"""Shared geometry/source primitives for the active OU-III P4 proof route.

The former candidate-full-word search driver mixed an obsolete angle ladder,
retired entrance/sector diagnostics, and an old P5 prefix backend.  Only the
source/geometry primitives consumed by the current joint-Joseph/augmented route
remain here.  This module is not a certificate producer and has no CLI.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from ou3_interval import Interval
import ou3_full_process_ucc as PROCESS
import ou3_p4_covariance_primitives as COV
import ou3_validated_transcendentals as VT


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _box(a: float) -> Interval:
    a = abs(float(a))
    return Interval(down(-a), up(a))


def _indices(mode: str):
    if mode not in ("H", "A"):
        raise ValueError("mode must be H or A")
    n = 18 if mode == "H" else 21
    return {
        "n": n,
        "BG": range(3, 6),
        "V": range(6, 9),
        "P": range(9, 12),
        "SS": range(12, 15),
        "AW": range(15, 18),
        "BA": range(18, 21) if mode == "A" else range(18, 18),
    }


def _norm_bounds_box(box) -> tuple[float, float]:
    lo2 = 0.0
    hi2 = 0.0
    for a, b in box:
        if a <= 0.0 <= b:
            mn = 0.0
        else:
            mn = min(abs(a), abs(b))
        mx = max(abs(a), abs(b))
        lo2 = down(lo2 + down(mn * mn))
        hi2 = up(hi2 + up(mx * mx))
    return down(math.sqrt(max(0.0, lo2))), up(math.sqrt(max(0.0, hi2)))


def _ball_box_cover(q: float, max_box_norm_factor: float = 1.5) -> list[list[tuple[float, float]]]:
    """Outward box cover of ||c||<=q with bounded corner inflation."""
    q = float(q)
    factor = float(max_box_norm_factor)
    if not (q > 0.0 and 1.0 < factor < math.sqrt(3.0)):
        raise ValueError("invalid Cayley ball-cover parameters")
    todo = [([(-q, q), (-q, q), (-q, q)], 0)]
    out = []
    while todo:
        box, depth = todo.pop()
        qmin, qmax = _norm_bounds_box(box)
        if qmin > q:
            continue
        if qmax <= factor * q:
            out.append([(down(a), up(b)) for a, b in box])
            continue
        if depth >= 18:
            raise RuntimeError("Cayley ball cover failed to reach requested inflation")
        widths = [b - a for a, b in box]
        axis = max(range(3), key=lambda i: widths[i])
        a, b = box[axis]
        m = 0.5 * (a + b)
        left = list(box)
        right = list(box)
        left[axis] = (a, m)
        right[axis] = (m, b)
        todo.append((right, depth + 1))
        todo.append((left, depth + 1))
    out.sort(key=lambda b: _norm_bounds_box(b)[1], reverse=True)
    return out


def _configure_mode(mode: str) -> None:
    _indices(mode)


def _source_sigma_bacc0() -> float:
    text = COV.MEKF.read_text(encoding="utf-8")
    m = re.search(r"sigma_bacc0_\s*=\s*T\(([0-9.eE+-]+)\)", text)
    if not m:
        raise RuntimeError("cannot extract source sigma_bacc0")
    return float(m.group(1))


def _transition_and_Q(mode: str, src: dict, domain: dict):
    """Dependency-preserving H/A transition and process covariance."""
    idx = _indices(mode)
    F, Q, Rstep = COV.tight_transition_and_Q(idx["n"], src, domain)
    if mode == "H":
        return F, Q, Rstep, None

    h = float(src["dt_s"])
    pc = PROCESS.build()["source_constants"]
    tau_b = float(pc["accel_bias_tau_s"])
    q_b = float(pc["accel_bias_process_variance_density"])
    if not (tau_b > 0.0 and q_b > 0.0):
        raise RuntimeError("active accelerometer-bias source constants lost positivity")

    x1 = Interval.outward_bounds(h / tau_b, h / tau_b)
    phi = VT.exp_interval(-x1)
    x2 = Interval.outward_bounds(2.0 * h / tau_b, 2.0 * h / tau_b)
    em1 = VT.expm1_interval(-x2)
    qd_scale = Interval.outward_bounds(-0.5 * tau_b, -0.5 * tau_b) * em1
    qd = Interval.outward_bounds(q_b, q_b) * qd_scale
    for i in idx["BA"]:
        F[i][i] = Interval(phi.lo, min(1.0, phi.hi))
        Q[i][i] = qd
    return F, COV.psd_tighten(Q), Rstep, {
        "phi_interval": phi.as_list(),
        "Qd_variance_interval": qd.as_list(),
        "tau_bacc_s": tau_b,
        "process_variance_density": q_b,
    }


def _initial_covariance(mode: str, src: dict, domain_path: Path):
    """GoLive covariance enclosure used by the active P4 design backends."""
    idx = _indices(mode)
    Pm = COV.initial_covariance(idx["n"], src, domain_path)
    if mode == "A":
        s = _source_sigma_bacc0()
        v = Interval.outward_bounds(s * s, s * s)
        for i in idx["BA"]:
            Pm[i][i] = v
            for j in range(idx["n"]):
                if j != i:
                    Pm[i][j] = I(0.0)
                    Pm[j][i] = I(0.0)
    return COV.psd_tighten(Pm)


def _initial_error(mode: str, domain: dict) -> tuple[list[Interval], float | None, dict]:
    """Current P4 entrance error box, including the finite 0.5 Hs position bound."""
    idx = _indices(mode)
    old = domain["startup"]["physical_handoff_coordinate_bounds"]
    pos = domain["initial_filter_entrance"]["position"]
    hs = float(pos["significant_wave_height_Hs_upper_m"])
    p_factor = float(pos["component_abs_error_upper_Hs_factor"])
    if not (hs > 0.0 and p_factor == 0.5):
        raise RuntimeError("finite 0.5 Hs entrance box not declared")
    p_component = up(p_factor * hs)

    e = [I(0.0) for _ in range(idx["n"])]
    for idxs, key in (
        (idx["BG"], "gyro_bias_error_norm_upper_rad_s"),
        (idx["V"], "velocity_error_norm_upper_mps"),
        (idx["SS"], "integral_displacement_error_norm_upper_m_s"),
        (idx["AW"], "latent_acceleration_error_norm_upper_mps2"),
    ):
        a = float(old[key])
        for i in idxs:
            e[i] = _box(a)
    for i in idx["P"]:
        e[i] = _box(p_component)

    ba_cap = None
    if mode == "A":
        ba_cap = float(domain["normal_live"]["active_accelerometer_bias_state_norm_upper_mps2"])
        if not 0.0 < ba_cap < float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"]):
            raise RuntimeError("A entrance bias ball is not strictly inside projection ball")
        for i in idx["BA"]:
            e[i] = _box(ba_cap)

    return e, ba_cap, {
        "Hs_upper_m": hs,
        "position_component_abs_upper_m": p_component,
        "legacy_P1_position_norm_upper_m_not_used_as_P4_entry": float(old["position_error_norm_upper_m"]),
    }
