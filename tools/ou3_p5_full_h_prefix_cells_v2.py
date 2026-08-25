#!/usr/bin/env python3
"""Tight source-complete full-H prefix backend for OU-III P5.

Version 1 established the required 18x18/Joseph/reset/signed-Cayley plumbing,
but its single broad tau interval was evaluated by the natural interval
extension of products such as ``tau*(1-exp(-h/tau))``.  That destroys the
shared tau dependency: over the deployed invariant it enclosed the 5 ms
``a_w->v`` coefficient by values as large as roughly 0.18 s although the exact
coefficient is always below 5 ms.  The resulting first-prediction covariance
was therefore a numerical dependency artifact, not a source obstruction.

This active backend keeps the v1 full-matrix transport but replaces that one
piece by source-uniform *dependency-preserving* bounds:

* the three integrated-OU transition coefficients are positive kernel
  integrals and hence monotone in tau; evaluate the shipping formula at the two
  tau endpoints and hull those validated point evaluations;
* the integrated-OU process covariance is bounded directly from its continuous
  positive kernels.  With q_c<=2 sigma_max^2/tau_min and
  ``k_[v,p,S,a](u) <= [u,u^2/2,u^3/6,1]``, every covariance entry has an
  explicit moment upper.  No cancellation-prone broad-x formula is needed;
* the goLive gyro-bias covariance is not reset by
  ``initialize_from_attitude``.  It is therefore bounded by its constructor
  value plus the full startup-time random-walk accumulation; warmup physical
  measurements can only reduce that covariance.  This corrects the optimistic
  Pb0-only seed used by v1.

All later operations are still the v1 full 18x18 interval map: same-cell
P,H,R,S,K,r,d_eff, shipping Joseph update, immediate reset congruence, exact
magnetometer tangent reduction, effective accelerometer a_w input, and the
actual interval ``d=-E_theta K r`` fed to the signed ``a^T c`` Cayley map.
No older scalar certificate is used for promotion.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_full_process_ucc as PROCESS
import ou3_p5_full_h_prefix_cells as V1
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 2


def I(x: float) -> Interval:
    return Interval.point(float(x))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _shipping_coeffs_point(tau: float, h: float) -> tuple[Interval, Interval, Interval, Interval]:
    """Validated point evaluation of shipping alpha/phi_va/phi_pa/phi_Sa."""
    tau = float(tau)
    h = float(h)
    if not (tau > 0.0 and h > 0.0):
        raise ValueError("positive tau/h required")
    x0 = h / tau
    x = Interval.outward_bounds(x0, x0)
    tt = Interval.outward_bounds(tau, tau)
    em1 = VT.expm1_interval(-x)
    alpha = VT.exp_interval(-x)
    phi_va = -tt * em1
    if abs(x0) < 1.0e-2:
        x2 = x.square()
        x3 = x2*x
        x4 = x3*x
        x5 = x4*x
        phi_pa = tt.square() * (
            I(0.5)*x2 - I(1.0/6.0)*x3 + I(1.0/24.0)*x4
        )
        phi_Sa = (tt.square()*tt) * (
            I(1.0/6.0)*x3 - I(1.0/24.0)*x4 + I(1.0/120.0)*x5
        )
    else:
        phi_pa = tt.square() * (x + em1)
        phi_Sa = (tt.square()*tt) * (I(0.5)*x.square() - x - em1)
    return alpha, phi_va, phi_pa, phi_Sa


def _monotone_coeff_hull(tau: Interval, h: float):
    lo = _shipping_coeffs_point(tau.lo, h)
    hi = _shipping_coeffs_point(tau.hi, h)
    # alpha and every integrated initial-a coefficient are integrals of positive
    # kernels that increase with tau.  Hulling endpoint point-enclosures is
    # therefore source complete over the whole tau interval.
    return tuple(Interval(a.lo, b.hi) for a, b in zip(lo, hi))


def _ou_process_moment_axis(tau: Interval, sigma: Interval, h: float):
    """4x4 [v,p,S,a] Q enclosure from positive impulse-response moments."""
    qc = up(2.0 * sigma.hi * sigma.hi / tau.lo)
    powers = (1, 2, 3, 0)
    facts = (1.0, 2.0, 6.0, 1.0)
    Q = V1._zero(4, 4)
    for i in range(4):
        for j in range(i, 4):
            p = powers[i] + powers[j]
            den = facts[i] * facts[j] * (p + 1.0)
            ub = up(qc * (h ** (p + 1)) / den)
            z = Interval(0.0, ub)
            Q[i][j] = z
            Q[j][i] = z
    return Q


def _tight_transition_and_Q(src: dict, domain: dict):
    h = float(src["dt_s"])
    tau = src["tau_s"]
    sigma = src["sigma_aw_mps2"]
    alpha, phi_va, phi_pa, phi_Sa = _monotone_coeff_hull(tau, h)

    F = V1.matrix_identity(V1.N)
    rate = float(domain["normal_live"]["body_rate_norm_upper_deg_s"])*math.pi/180.0
    Rstep, Bstep = V1._rotation_step_box(rate, h)
    for i in range(3):
        for j in range(3):
            F[i][j] = Rstep[i][j]
            F[i][3+j] = Bstep[i][j]
    for ax in range(3):
        iv, ip, iS, ia = 6+ax, 9+ax, 12+ax, 15+ax
        F[iv][iv] = I(1.0); F[iv][ia] = phi_va
        F[ip][iv] = I(h); F[ip][ip] = I(1.0); F[ip][ia] = phi_pa
        F[iS][iv] = I(0.5*h*h); F[iS][ip] = I(h); F[iS][iS] = I(1.0); F[iS][ia] = phi_Sa
        F[ia][ia] = Interval(alpha.lo, min(1.0, alpha.hi))

    Q = V1._zero(V1.N, V1.N)
    proc = PROCESS.build()["source_constants"]
    qg = float(proc["gyro_noise_density_rad_sqrt_s"])**2
    qb = float(proc["gyro_bias_rw_variance_density"])
    bb = up(qb*h*h*h/3.0)
    cross = up(qb*h*h/2.0)
    for i in range(3):
        Q[i][i] = Interval(down(qg*h), up(qg*h+bb))
        Q[3+i][3+i] = Interval.outward_bounds(qb*h, qb*h)
        for j in range(3):
            if i != j:
                Q[i][j] = Interval(-bb, bb)
            Q[i][3+j] = Interval(-cross, cross)
            Q[3+j][i] = Q[i][3+j]

    qaxis = _ou_process_moment_axis(tau, sigma, h)
    groups = (6, 9, 12, 15)
    for ax in range(3):
        ids = [g+ax for g in groups]
        for i in range(4):
            for j in range(4):
                Q[ids[i]][ids[j]] = qaxis[i][j]
    return F, V1._psd_tighten(Q), Rstep


def _startup_timeout_s() -> float:
    text = WRAPPER.read_text(encoding="utf-8")
    m = re.search(r"proxy_startup_timeout_sec\s*=\s*([0-9.eE+-]+)f", text)
    if not m:
        raise RuntimeError("cannot extract source startup timeout")
    t = float(m.group(1))
    if not t > 0.0:
        raise RuntimeError("invalid source startup timeout")
    return t


def _corrected_initial_covariance(src: dict, domain_path: Path):
    Pm = V1._initial_covariance_original(src, domain_path)
    proc = PROCESS.build()["source_constants"]
    qb = float(proc["gyro_bias_rw_variance_density"])
    pb0 = V1._source_pb0()
    timeout = _startup_timeout_s()
    pbg_hi = up(pb0 + qb*timeout)
    for i in V1.BG:
        Pm[i][i] = Interval(0.0, pbg_hi)
        for j in V1.BG:
            if i != j:
                b = up(pbg_hi)
                Pm[i][j] = Interval(-b, b)
    return V1._psd_tighten(Pm)


def _install_backend() -> None:
    # Preserve the unpatched seed function exactly once so repeated build() calls
    # in one unittest process do not recurse through this wrapper.
    if not hasattr(V1, "_initial_covariance_original"):
        V1._initial_covariance_original = V1._initial_covariance
    V1._transition_and_Q = _tight_transition_and_Q
    V1._initial_covariance = _corrected_initial_covariance


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    _install_backend()
    out = dict(V1.build(Path(domain_path).resolve()))
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_TIGHT_FULL_18X18_H_PREFIX_INTERVAL_CELL_PROPAGATION"
    out["active_full_matrix_backend"] = "DEPENDENCY_PRESERVING_OU_KERNEL_BOUNDS"
    out["broad_tau_natural_interval_product_used"] = False
    out["integrated_ou_transition_monotone_endpoint_hull_used"] = True
    out["integrated_ou_process_positive_kernel_moment_bounds_used"] = True
    out["goLive_gyro_bias_covariance_includes_full_startup_RW_upper"] = True
    out["startup_timeout_s_used_for_gyro_bias_covariance"] = _startup_timeout_s()
    out["old_v1_natural_interval_prediction_is_promotion_route"] = False
    return out


def validate(d: dict) -> list[str]:
    # V1 validates all transport semantics but expects schema 1.  Validate a
    # schema-adjusted view, then add the dependency-preserving requirements.
    base = dict(d)
    base["schema"] = V1.SCHEMA
    failures = V1.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("active_full_matrix_backend") != "DEPENDENCY_PRESERVING_OU_KERNEL_BOUNDS":
        failures.append("tight full-matrix backend is not active")
    if d.get("broad_tau_natural_interval_product_used") is not False:
        failures.append("broad dependent tau product remains active")
    for k in (
        "integrated_ou_transition_monotone_endpoint_hull_used",
        "integrated_ou_process_positive_kernel_moment_bounds_used",
        "goLive_gyro_bias_covariance_includes_full_startup_RW_upper",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    if d.get("old_v1_natural_interval_prediction_is_promotion_route") is not False:
        failures.append("v1 dependency-artifact path remains a promotion route")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"],
        "q8_closed": out["complete_q_le_8_prefix_family_closed"],
        "max_q": out["max_reached_cayley_norm_upper"],
        "smaller_chart": out["smaller_source_reachable_chart_upper"],
        "backend": out["active_full_matrix_backend"],
        "first_failure": out["first_failure"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
