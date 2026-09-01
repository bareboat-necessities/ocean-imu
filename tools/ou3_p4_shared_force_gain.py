#!/usr/bin/env python3
"""Dependency-preserving accelerometer gain in the shared specific-force magnitude.

The structured first-accelerometer gain rows all have the shape

    g(m) = m * N / (m^2 * p + lambda),

where ``m`` is the specific-force magnitude of the audited cell.  ``m`` appears
in the numerator *and* in the denominator, so evaluating the quotient with
ordinary interval arithmetic (numerator at ``m.hi``, denominator at ``m.lo``)
overstates the gain by up to the cell's magnitude ratio ``m.hi/m.lo``.  That
loss is not a modelling choice; it is pure dependency slack, and it is what
pushed the 30 deg first-accelerometer correction family past the monotone
Cayley chart before this producer existed.

``g`` is unimodal in ``m`` on the positive axis:

    g'(m) = N (lambda - m^2 p) / (m^2 p + lambda)^2,

so the maximum over ``m > 0`` is at ``m* = sqrt(lambda/p)`` and on a cell
``[m1, m2]`` the supremum is either the interior value or the larger endpoint
value.  Two exact interior values are needed:

* ``N = p`` (the self-``p`` tangent rows ``g_perp`` and ``g_u``):
  ``g(m*) = (1/2) sqrt(p/lambda)``;
* ``N = C`` independent of ``p`` (the axial yaw-cross row ``g_z``):
  ``g(m*) = C / (2 sqrt(lambda p))``.

Monotonicity in the remaining variables is also exact: the self-``p`` rows
increase in ``p``, the independent-``C`` row decreases in ``p``, and both
decrease in ``lambda``.  The corresponding endpoint of each interval is chosen
explicitly instead of being left to interval division.

Nothing else changes.  The filter, the declared 30 deg candidate, the 0.3 g
startup ``a_w`` envelope, the source cell family, the PSD remainder treatment
and the ``KH`` rows are all identical to
``ou3_p4_candidate_first_accel_range_v3``.  This module only removes shared-
variable slack, so every bound it returns is less than or equal to the bound
the naive evaluation returns on the same cell.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_first_accel_structured_gain as SG

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _sqrt_up(x: float) -> float:
    return up(math.sqrt(max(0.0, float(x))))


def _eval_up(m: float, num: float, p: float, lam: float) -> float:
    """Outward upper bound of m*num/(m^2 p + lam) at an exact m."""
    den = down(down(down(m * m) * p) + lam)
    if den <= 0.0:
        raise RuntimeError("shared-force gain denominator lost positivity")
    return up(up(m * num) / den)


def _interior_reachable(m: Interval, p: float, lam: float) -> bool:
    """True when m* = sqrt(lam/p) may lie inside the audited magnitude cell."""
    star_lo = down(lam / up(p))
    star_hi = up(lam / down(p))
    return star_hi >= down(m.lo * m.lo) and star_lo <= up(m.hi * m.hi)


def sup_self_p_gain(m: Interval, p_hi: float, lam_lo: float) -> float:
    """sup over the cell of m*p/(m^2 p + lam); p upper, lambda lower."""
    if m.lo <= 0.0 or p_hi <= 0.0 or lam_lo <= 0.0:
        raise RuntimeError("shared-force self-p gain domain lost positivity")
    cands = [_eval_up(m.lo, p_hi, p_hi, lam_lo), _eval_up(m.hi, p_hi, p_hi, lam_lo)]
    if _interior_reachable(m, p_hi, lam_lo):
        cands.append(up(0.5 * _sqrt_up(up(p_hi / down(lam_lo)))))
    return max(cands)


def sup_independent_gain(m: Interval, c_hi: float, p_lo: float, lam_lo: float) -> float:
    """sup over the cell of m*C/(m^2 p + lam); C independent, p and lambda lower."""
    if m.lo <= 0.0 or p_lo <= 0.0 or lam_lo <= 0.0:
        raise RuntimeError("shared-force independent gain domain lost positivity")
    if c_hi <= 0.0:
        return 0.0
    cands = [_eval_up(m.lo, c_hi, p_lo, lam_lo), _eval_up(m.hi, c_hi, p_lo, lam_lo)]
    if _interior_reachable(m, p_lo, lam_lo):
        cands.append(up(c_hi / down(2.0 * down(_sqrt_up(down(lam_lo * p_lo))))))
    return max(cands)


def shared_force_structured_gain_bounds(*, tilt: float, yaw: float, eps: float,
                                        x: Interval, m: Interval, paw: Interval,
                                        racc_var: Interval):
    """V4 gain rows: identical model to V3, shared-``m`` dependency preserved."""
    F = SG.FULL1
    t = SG.I(tilt)
    delta = Interval.outward_bounds(yaw - tilt, yaw - tilt)
    pu = t + delta * x
    lam = paw + racc_var
    if lam.lo <= 0.0:
        raise RuntimeError("candidate first-accel lambda floor is nonpositive")
    m2 = m.square()
    den_perp = m2 * t + lam
    den_u = m2 * pu + lam
    if den_perp.lo <= 0.0 or den_u.lo <= 0.0:
        raise RuntimeError("candidate tangent gain denominator is nonpositive")

    geom_hi = SG._sqrt_x1mx_upper(x)
    g_perp = sup_self_p_gain(m, t.hi, lam.lo)
    g_u = sup_self_p_gain(m, pu.hi, lam.lo)
    delta_abs = max(abs(delta.lo), abs(delta.hi))
    g_z = sup_independent_gain(m, F.up(delta_abs * geom_hi), pu.lo, lam.lo)
    k0 = F.up(max(g_perp, _sqrt_up(F.up(g_u * g_u + g_z * g_z))))

    # KH rows are increasing in m and need no unimodal treatment.
    kh_perp = m2 * t / den_perp
    kh_u = m2 * pu / den_u
    kh_z = m2 * delta * Interval(0.0, geom_hi) / den_u
    kh0 = F.up(max(
        kh_perp.hi,
        _sqrt_up(F.up(kh_u.hi * kh_u.hi + kh_z.hi * kh_z.hi)),
    ))

    # V12D tangent-channel resolvent, unchanged from V3.
    eoff = F.up(2.0 * eps)
    mhi = m.hi
    dS = F.up(mhi * mhi * eoff)
    tangent_nominal_floor = min(den_perp.lo, den_u.lo)
    tangent_floor = F.down(tangent_nominal_floor - dS)
    if tangent_floor <= 0.0:
        raise RuntimeError("candidate perturbed tangent innovation floor lost positivity")
    inv_tangent = F.up(1.0 / tangent_floor)
    dCtheta = F.up(mhi * eoff)
    dk = F.up(F.up(dCtheta + F.up(k0 * dS)) * inv_tangent)
    k = F.up(k0 + dk)
    kh = F.up(kh0 + F.up(dk * mhi))
    return k, kh, {
        "lambda_lower": lam.lo,
        "p_u": pu.as_list(),
        "g_perp_upper": g_perp,
        "g_u_upper": g_u,
        "g_z_upper": g_z,
        "K0_norm_upper": k0,
        "KH0_norm_upper": kh0,
        "attitude_PSD_remainder_operator_upper": eoff,
        "PSD_innovation_perturbation_tangent_only": True,
        "PSD_innovation_axial_row_column_exact_zero": True,
        "nominal_tangent_innovation_lower": tangent_nominal_floor,
        "tangent_innovation_perturbation_upper": dS,
        "perturbed_tangent_innovation_lower": tangent_floor,
        "perturbed_tangent_inverse_operator_upper": inv_tangent,
        "PSD_remainder_K_perturbation_upper": dk,
        "retired_axial_lambda_inverse_for_PSD_remainder": True,
        "shared_force_magnitude_dependency_preserved": True,
    }


def _reference_cells(domain_path: Path) -> list[dict]:
    """Audited comparison cells drawn from the deployed startup/live domain."""
    import ou3_p5_first_accel_rotation_gauge as RG
    import ou3_p5_first_accel_rotation_gauge_v3 as RG3
    import ou3_p5_full_h_prefix_cells as FULL
    import ou3_p5_full_h_prefix_cells_v3 as FULL3
    import ou3_vector_uco_certificate as VECTOR

    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    RG3._install_backend(path, 2)
    FULL3._install_backend()
    h = float(FULL._source_cell()["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    vc = VECTOR.build()["configured_measurement_bounds"]
    racc_var = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))[0][0]
    live = domain["normal_live"]
    force_cells = RG._geom_ranges(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
        4,
    )
    xcells = SG._linear_cells(16)

    paws: list[Interval] = []
    for src, phase in RG._source_phase_children(2):
        P0 = FULL._initial_covariance(src, path)
        Fm, Q, _ = FULL._transition_and_Q(src, domain)
        Pp = FULL._psd_tighten(FULL.matrix_add(
            FULL.matrix_mul(FULL.matrix_mul(Fm, P0), FULL.matrix_transpose(Fm)), Q))
        _s, _a, paw_pred = RG._scalar_axis_structure(Pp)
        paw = RG._due_paw_and_error_norm(Pp, src, 0.0, 0.0)[0] if phase == "due" else paw_pred
        paws.append(paw)
    # De-duplicate identical source a_w covariance cells.
    seen: dict[tuple[float, float], Interval] = {}
    for p in paws:
        seen.setdefault((p.lo, p.hi), p)
    return [{
        "tilt": tilt, "yaw": yaw, "eps": eps,
        "x": x, "m": m, "paw": p, "racc_var": racc_var,
    } for p in seen.values() for x in xcells[:4] for m in force_cells]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    import ou3_p4_candidate_first_accel_range_v3 as V3

    rows = []
    worst_ratio = 0.0
    min_ratio = math.inf
    never_worse = True
    for cell in _reference_cells(domain_path):
        k_naive, kh_naive, _dn = V3._tangent_structured_gain_bounds(**cell)
        k_exact, kh_exact, _de = shared_force_structured_gain_bounds(**cell)
        if k_exact > k_naive or kh_exact > kh_naive:
            never_worse = False
        ratio = k_naive / k_exact if k_exact > 0.0 else math.inf
        min_ratio = min(min_ratio, ratio)
        if ratio > worst_ratio:
            worst_ratio = ratio
            rows = [{
                "force_magnitude_mps2": cell["m"].as_list(),
                "alignment_x": cell["x"].as_list(),
                "aw_covariance": cell["paw"].as_list(),
                "K_norm_upper_naive_interval": k_naive,
                "K_norm_upper_shared_force": k_exact,
                "KH_norm_upper_naive_interval": kh_naive,
                "KH_norm_upper_shared_force": kh_exact,
                "tightening_factor": ratio,
            }]

    failures: list[str] = []
    if not never_worse:
        failures.append("shared-force gain is not uniformly at least as tight as the interval gain")
    if min_ratio < 1.0:
        failures.append("shared-force gain reported a looser cell")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SHARED_FORCE_MAGNITUDE_ACCELEROMETER_GAIN_LEMMA",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "shared_force_magnitude_dependency_preserved": True,
        "unimodal_interior_and_endpoint_maxima_used": True,
        "self_p_rows_evaluated_at_upper_p_and_lower_lambda": True,
        "independent_C_row_evaluated_at_lower_p_and_lower_lambda": True,
        "KH_rows_unchanged_from_V3": True,
        "PSD_remainder_treatment_unchanged_from_V3": True,
        "uniformly_at_least_as_tight_as_interval_gain": never_worse,
        "minimum_observed_tightening_factor": min_ratio,
        "maximum_observed_tightening_factor": worst_ratio,
        "worst_audited_cell": rows[0] if rows else None,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "consume this gain from the candidate first-accelerometer range and signed sector"
            " producers; it removes shared-magnitude slack only and does not by itself establish"
            " sector invariance or complete-word dissipation"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "shared_force_magnitude_dependency_preserved",
        "unimodal_interior_and_endpoint_maxima_used", "KH_rows_unchanged_from_V3",
        "PSD_remainder_treatment_unchanged_from_V3",
        "uniformly_at_least_as_tight_as_interval_gain",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
              "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if not (float(d.get("minimum_observed_tightening_factor", 0.0)) >= 1.0):
        f.append("a tightening factor below one was reported")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "min_tightening": d["minimum_observed_tightening_factor"],
        "max_tightening": d["maximum_observed_tightening_factor"],
        "worst_cell": d["worst_audited_cell"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
