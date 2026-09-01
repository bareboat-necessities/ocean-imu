#!/usr/bin/env python3
"""Joint force pairing and the a_w/sigma consistency constant for P4's first accelerometer.

``ou3_p4_first_accel_sector_budget`` shows that the nuisance part of the first
deployed accelerometer correction exceeds the sector-invariance budget at every
candidate angle, and that shrinking the candidate cannot close it.  This
producer isolates the two remaining reasons and reports what each is worth.

**Unconditional part -- joint force pairing.**  The nuisance term is
``||K_theta|| * (a_w error + eta(q) |f| + b_a)``.  The gain falls roughly like
``1/|f|`` while the finite-angle force remainder grows like ``|f|``, so bounding
the two factors separately over one force cell charges the largest gain against
the largest force.  Maximising the product over subdivided magnitudes inside the
same cell is an exact dependency-preserving tightening: no domain, filter or
candidate changes, and the result is never looser.

**Conditional part -- the a_w/sigma pairing.**  The a_w covariance seeding the
gain is the tuner state, ``P_aw = sigma_applied^2`` with
``sigma_applied`` in the deployed safety range ``[0.05, 6.0] m/s^2``.  The a_w
*error* is the separately declared ``0.3 g`` startup envelope.  The two are
currently free to be chosen independently, so the worst cell pairs a tuner that
believes the sea is flat (``sigma_applied = 0.05``, hence the largest gain) with
a ``2.94 m/s^2`` latent-acceleration error -- a 60-sigma pairing.

The shipping tuner does couple them: ``sigma_target_ = min(sigma_wave *
sigma_coeff_, max_sigma_a_)`` with ``sigma_wave`` the estimated band-limited
wave-acceleration RMS, and ``sigma_applied`` an EMA toward that target.  But no
*declared domain statement* currently bounds the latent-acceleration error in
terms of the tuner state, and the EMA transient means one cannot be inferred
from the update law alone.

So this producer does not assume a coupling.  It **measures the constant a
coupling would have to supply**: the largest ``c`` such that adding

    ||delta a_w|| <= c * sigma_applied

to the declared domain brings the first-accelerometer nuisance term inside the
sector-invariance budget, reported per candidate angle.  ``c = inf`` reproduces
today's unconstrained pairing.  The result is explicitly conditional and
promotes nothing: it converts an open search into one declared-domain question
with a number attached.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_first_accel_sector_budget as BUDGET
import ou3_p4_candidate_aw_capture_budget as AWB
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p4_shared_force_gain as SHARED
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
LADDER_DEG = (15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0)
CONSISTENCY_SEARCH_UPPER = 1024.0
CONSISTENCY_BISECTION_STEPS = 50


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _force_subcells(m: Interval, pieces: int) -> list[Interval]:
    """Geometric subdivision of one audited force cell, endpoints preserved."""
    if pieces < 1:
        raise ValueError("positive force subdivision required")
    if pieces == 1:
        return [m]
    return RG._geom_ranges(m.lo, m.hi, pieces)


def _gain_table(path: Path, domain: dict, *, source_pieces: int,
                alignment_pieces: int, force_magnitude_pieces: int,
                force_subpieces: int) -> dict:
    """Angle-independent gain rows, one per subdivided force magnitude."""
    RG3._install_backend(path, source_pieces)
    FULL3._install_backend()
    h = float(FULL._source_cell()["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    vc = VECTOR.build()["configured_measurement_bounds"]
    racc_var = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))[0][0]
    startup = domain["startup"]
    handoff = startup["physical_handoff_coordinate_bounds"]
    aw0 = float(handoff["latent_acceleration_error_norm_upper_mps2"])
    ba = float(handoff["accelerometer_bias_error_norm_upper_mps2"])
    pnorm = AWB._p5_position_norm_upper(domain)
    live = domain["normal_live"]
    force_lower = float(live["specific_force_norm_lower_mps2"])
    force_upper = float(live["specific_force_norm_upper_mps2"])
    force_cells = RG._geom_ranges(force_lower, force_upper, force_magnitude_pieces)
    xcells = SG._linear_cells(alignment_pieces)

    rows = []
    sigma_lo = math.inf
    sigma_hi = 0.0
    for si, (src, phase) in enumerate(RG._source_phase_children(source_pieces)):
        sigma = src["sigma_aw_mps2"]
        sigma_lo = min(sigma_lo, sigma.lo)
        sigma_hi = max(sigma_hi, sigma.hi)
        P0 = FULL._initial_covariance(src, path)
        Fm, Q, _ = FULL._transition_and_Q(src, domain)
        Pp = FULL._psd_tighten(FULL.matrix_add(
            FULL.matrix_mul(FULL.matrix_mul(Fm, P0), FULL.matrix_transpose(Fm)), Q))
        _s, _a, paw_pred = RG._scalar_axis_structure(Pp)
        paw = RG._due_paw_and_error_norm(Pp, src, 0.0, 0.0)[0] if phase == "due" else paw_pred
        intercept, slope, _d = AWB._s_phase_affine_aw_bound(src, phase, Pp, domain, pnorm)
        aw = up(intercept + up(slope * aw0))
        for xi, x in enumerate(xcells):
            for mi, m in enumerate(force_cells):
                for mj, ms in enumerate(_force_subcells(m, force_subpieces)):
                    k, _kh, _gd = SHARED.shared_force_structured_gain_bounds(
                        tilt=tilt, yaw=yaw, eps=eps, x=x, m=ms,
                        paw=paw, racc_var=racc_var)
                    rows.append({
                        "source_phase_cell": si, "pseudo_phase": phase,
                        "alignment_cell": xi, "force_cell": mi, "force_subcell": mj,
                        "force_magnitude_mps2": ms.as_list(),
                        "force_upper_mps2": ms.hi,
                        "alignment_x": x.as_list(),
                        "tuner_sigma_aw_mps2": sigma.as_list(),
                        "tuner_sigma_aw_upper_mps2": sigma.hi,
                        "aw_after_prefix_upper_mps2": aw,
                        "Ktheta_norm_upper": k,
                    })
    return {
        "rows": rows,
        "accel_bias_error_norm_upper_mps2": ba,
        "declared_startup_aw_error_norm_upper_mps2": aw0,
        "specific_force_norm_lower_mps2": force_lower,
        "specific_force_norm_upper_mps2": force_upper,
        "tuner_sigma_aw_lower_mps2": sigma_lo,
        "tuner_sigma_aw_upper_mps2": sigma_hi,
        "prediction_step_s": h,
    }


def _worst_nuisance(table: dict, eta_per_vector: float, c: float) -> tuple[float, dict]:
    """Worst nuisance correction over all rows for a consistency constant ``c``."""
    ba = table["accel_bias_error_norm_upper_mps2"]
    worst = 0.0
    witness = None
    for r in table["rows"]:
        aw = r["aw_after_prefix_upper_mps2"]
        if math.isfinite(c):
            aw = min(aw, up(c * r["tuner_sigma_aw_upper_mps2"]))
        eta_force = up(eta_per_vector * r["force_upper_mps2"])
        rho = up(aw + up(eta_force + ba))
        n = up(r["Ktheta_norm_upper"] * rho)
        if n > worst:
            worst = n
            witness = dict(r)
            witness["applied_aw_error_upper_mps2"] = aw
            witness["force_attitude_remainder_upper_mps2"] = eta_force
            witness["effective_aw_input_upper_mps2"] = rho
            witness["nuisance_correction_norm_upper_rad"] = n
            witness["unconstrained_c_at_this_cell"] = up(
                r["aw_after_prefix_upper_mps2"] / down(r["tuner_sigma_aw_upper_mps2"]))
    return worst, witness


def _critical_consistency_constant(table: dict, eta_per_vector: float,
                                   budget: float) -> dict:
    """Largest ``c`` whose worst nuisance still fits inside ``budget``."""
    if budget <= 0.0:
        return {"critical_consistency_constant": 0.0,
                "any_finite_constant_closes_this_angle": False,
                "bracket": [0.0, 0.0]}
    zero, _w = _worst_nuisance(table, eta_per_vector, 0.0)
    if zero >= budget:
        # Even a perfect a_w estimate leaves the bias and finite-angle force
        # remainder above the budget; no consistency constant can help here.
        return {"critical_consistency_constant": 0.0,
                "any_finite_constant_closes_this_angle": False,
                "nuisance_at_zero_aw_error_rad": zero,
                "bracket": [0.0, 0.0]}
    lo = 0.0
    hi = CONSISTENCY_SEARCH_UPPER
    top, _w = _worst_nuisance(table, eta_per_vector, hi)
    if top < budget:
        return {"critical_consistency_constant": math.inf,
                "any_finite_constant_closes_this_angle": True,
                "nuisance_at_zero_aw_error_rad": zero,
                "bracket": [hi, math.inf],
                "unconstrained_pairing_already_fits": True}
    for _ in range(CONSISTENCY_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        val, _w = _worst_nuisance(table, eta_per_vector, mid)
        if val < budget:
            lo = mid
        else:
            hi = mid
    return {"critical_consistency_constant": down(lo),
            "any_finite_constant_closes_this_angle": True,
            "nuisance_at_zero_aw_error_rad": zero,
            "bracket": [down(lo), up(hi)]}


def _angle_row(deg: float, table: dict, domain: dict, outer_q: float,
               separate_budget: dict | None) -> dict:
    geom = SECTOR._validated_design_geometry(math.radians(float(deg)))
    q0 = float(geom["cayley_norm_upper"])
    q = RG._q_after_first_prediction(q0, domain, table["prediction_step_s"])
    budget = BUDGET._sector_correction_budget(q, outer_q)
    b = float(budget["sector_invariance_correction_budget_upper_rad"])
    eta_per_vector = VEFF.accel_attitude_eta_per_vector_norm_upper(q)

    joint, witness = _worst_nuisance(table, eta_per_vector, math.inf)
    crit = _critical_consistency_constant(table, eta_per_vector, b)

    row = {
        "angle_deg": float(deg),
        "candidate_q_upper": q0,
        "post_prediction_q_upper": q,
        "eta_per_vector_norm_upper": eta_per_vector,
        "sector_invariance_correction_budget_upper_rad": b,
        "nuisance_joint_force_pairing_upper_rad": joint,
        "nuisance_over_budget_ratio_joint": (math.inf if b <= 0.0 else up(joint / b)),
        "joint_pairing_fits_inside_budget": bool(b > 0.0 and joint < b),
        "worst_cell": witness,
        **crit,
    }
    unc = None if witness is None else witness.get("unconstrained_c_at_this_cell")
    row["unconstrained_c_at_worst_cell"] = unc
    cc = row["critical_consistency_constant"]
    row["consistency_tightening_needed_factor"] = (
        None if (unc is None or not (cc > 0.0)) else up(unc / cc))
    if separate_budget is not None:
        prev = float(separate_budget["nuisance_correction_norm_upper_rad"])
        row["nuisance_separate_force_bounds_upper_rad"] = prev
        row["joint_force_pairing_tightening_factor"] = (
            up(prev / joint) if joint > 0.0 else math.inf)
    return row


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4,
          force_subpieces: int = 16) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("consistency producer domain must not be trajectory fitted")

    sector = SECTOR.build(path)
    entrance = ENTRANCE.build(path)
    failures = [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    outer_q = float(sector["design_cayley_norm_upper"])

    separate = BUDGET.build(path, source_pieces=source_pieces,
                            alignment_pieces=alignment_pieces,
                            force_magnitude_pieces=force_magnitude_pieces)
    failures += [f"separate budget: {x}" for x in BUDGET.validate(separate)]
    by_angle = {float(r["angle_deg"]): r for r in separate["ladder_rows"]}

    table = _gain_table(path, domain, source_pieces=source_pieces,
                        alignment_pieces=alignment_pieces,
                        force_magnitude_pieces=force_magnitude_pieces,
                        force_subpieces=force_subpieces)

    rows = [_angle_row(d, table, domain, outer_q, by_angle.get(d))
            for d in LADDER_DEG]

    closing = [r for r in rows if r["any_finite_constant_closes_this_angle"]]
    widest = max((r["angle_deg"] for r in closing), default=None)
    widest_row = next((r for r in rows if r["angle_deg"] == widest), None)
    aw0 = table["declared_startup_aw_error_norm_upper_mps2"]
    sigma_lo = table["tuner_sigma_aw_lower_mps2"]

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_FIRST_ACCEL_JOINT_FORCE_PAIRING_AND_AW_SIGMA_CONSISTENCY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_entrance_shrunk": False,
        "aw_sigma_consistency_declared_in_domain": False,
        "consistency_constant_is_a_conditional_requirement_not_a_theorem": True,
        "joint_force_pairing_is_unconditional": True,
        "distance_only_no_verdict_emitted": True,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "operation_matched_outer_q_upper": outer_q,
        "declared_startup_aw_error_norm_upper_mps2": aw0,
        "accel_bias_error_norm_upper_mps2": table["accel_bias_error_norm_upper_mps2"],
        "specific_force_norm_lower_mps2": table["specific_force_norm_lower_mps2"],
        "specific_force_norm_upper_mps2": table["specific_force_norm_upper_mps2"],
        "tuner_sigma_aw_range_mps2": [sigma_lo, table["tuner_sigma_aw_upper_mps2"]],
        "unconstrained_pairing_sigma_ratio_at_floor": up(aw0 / down(sigma_lo)),
        "audited_gain_rows": len(table["rows"]),
        "force_subpieces": force_subpieces,
        "angle_rows": rows,
        "widest_angle_closed_by_a_finite_consistency_constant": widest,
        "consistency_constant_at_widest_closed_angle": (
            None if widest_row is None else widest_row["critical_consistency_constant"]),
        "shipping_tuner_sigma_law":
            "sigma_target = min(sigma_wave * sigma_coeff, max_sigma_a); "
            "sigma_applied is an EMA toward sigma_target; sigma_wave is the "
            "estimated band-limited wave acceleration RMS",
        "next_obligation": (
            "the joint force pairing is an unconditional tightening and may be consumed"
            " directly; the a_w/sigma consistency constant is not established here and"
            " must either be declared in the operating domain with its own justification"
            " from the tuner law and EMA transient, or replaced by carrying the tuner"
            " state jointly with the attitude in the complete-word P4 funnel"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in ("source_generated_not_trajectory_fit", "distance_only_no_verdict_emitted",
              "consistency_constant_is_a_conditional_requirement_not_a_theorem",
              "joint_force_pairing_is_unconditional"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "declared_entrance_shrunk",
              "aw_sigma_consistency_declared_in_domain",
              "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
              "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    rows = d.get("angle_rows", [])
    if [r["angle_deg"] for r in rows] != sorted(r["angle_deg"] for r in rows):
        f.append("angle rows are not ordered")
    if not rows:
        f.append("no angle rows emitted")
    for r in rows:
        if r["nuisance_joint_force_pairing_upper_rad"] <= 0.0:
            f.append(f"{r['angle_deg']}: non-positive joint nuisance")
        prev = r.get("nuisance_separate_force_bounds_upper_rad")
        if prev is not None and r["nuisance_joint_force_pairing_upper_rad"] > prev:
            f.append(f"{r['angle_deg']}: joint force pairing is looser than separate bounds")
        c = r["critical_consistency_constant"]
        if not (c >= 0.0):
            f.append(f"{r['angle_deg']}: negative consistency constant")
    if d.get("audited_gain_rows", 0) <= 0:
        f.append("no audited gain rows")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--force-subpieces", type=int, default=16)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_pieces=a.source_pieces,
              alignment_pieces=a.alignment_pieces,
              force_magnitude_pieces=a.force_magnitude_pieces,
              force_subpieces=a.force_subpieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "sigma_range": d["tuner_sigma_aw_range_mps2"],
        "aw_over_sigma_floor": d["unconstrained_pairing_sigma_ratio_at_floor"],
        "angles": [{
            "deg": r["angle_deg"],
            "budget": r["sector_invariance_correction_budget_upper_rad"],
            "nuisance_separate": r.get("nuisance_separate_force_bounds_upper_rad"),
            "nuisance_joint": r["nuisance_joint_force_pairing_upper_rad"],
            "joint_tightening": r.get("joint_force_pairing_tightening_factor"),
            "ratio": r["nuisance_over_budget_ratio_joint"],
            "fits_unconstrained": r["joint_pairing_fits_inside_budget"],
            "zero_aw_residual": r.get("nuisance_at_zero_aw_error_rad"),
            "critical_c": r["critical_consistency_constant"],
            "unconstrained_c": r["unconstrained_c_at_worst_cell"],
            "tightening_needed": r["consistency_tightening_needed_factor"],
            "closes_with_finite_c": r["any_finite_constant_closes_this_angle"],
        } for r in d["angle_rows"]],
        "widest_closed_angle": d["widest_angle_closed_by_a_finite_consistency_constant"],
        "c_at_widest": d["consistency_constant_at_widest_closed_angle"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
