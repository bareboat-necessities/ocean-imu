#!/usr/bin/env python3
"""Derive P4 latent-acceleration entrance budgets from the shipping update range.

The finite-angle candidate search must distinguish the declared startup/P5
latent-acceleration error envelope from the smaller radius, if any, actually
required by a P4 candidate.  The deployment domain now declares
``||e_aw|| <= 0.3 g = 2.941995 m/s^2`` for startup in waves; this is a physical
theorem assumption, not a replay fit or proof-specific tuning knob.

This producer independently solves the already-certified first-prefix
inequality backwards for the largest a_w entrance norm sufficient for the
unchanged 6 rad helper range, uniformly over the same source, S-pseudo phase,
force-magnitude and yaw-alignment cells used by
``ou3_p4_candidate_first_accel_range_v3``.  If the derived P4 radius contains
the declared 0.3 g startup envelope, no separate P5 a_w-only capture is needed.
If it is smaller, the difference remains an explicit P5 capture obligation.

For a not-due S phase,

    e_aw^- <= alpha_hi * A,

where A is the unknown P4 entrance radius.  For a due first S=0 update,

    e_aw^- <= alpha_hi A + |K_S| (S0 + h p0 + .5 h^2 v0 + phi_S_hi A)
            = c0 + c1 A.

The position term uses the declared P5/P4 entrance ``|delta p_i|<=0.5 Hs`` and
``Hs<=8.5 m``, hence ``||delta p||<=sqrt(3)*0.5*Hs``.  The existing v and S
bounds remain unchanged.  For each structured accelerometer child, the 6 rad
range requires

    K_theta * (rotation_residual + e_aw^- + b_a) <= 6.

Solving that inequality with directed rounding yields a certified *lower*
bound on the admissible A.  No capture, complete-word dissipation or usable P4
theorem is promoted here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_candidate_first_accel_range as BASE
import ou3_p4_candidate_first_accel_range_v3 as V3
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = 2
LIMIT = 6.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _safe_sub_lower(a: float, b: float) -> float:
    return down(float(a) - float(b))


def _p5_position_norm_upper(domain: dict) -> float:
    p = domain["initial_filter_entrance"]["position"]
    hs = float(p["significant_wave_height_Hs_upper_m"])
    fac = float(p["component_abs_error_upper_Hs_factor"])
    return up(math.sqrt(3.0) * up(fac * hs))


def _s_phase_affine_aw_bound(src: dict, phase: str, Pp, domain: dict,
                             position_norm_upper: float) -> tuple[float, float, dict]:
    """Return e_aw^- <= intercept + slope*A for unknown entrance radius A."""
    h = float(src["dt_s"])
    alpha, _pv, _pp, phiS = RG.V2._monotone_coeff_hull(src["tau_s"], h)
    slope = float(alpha.hi)
    intercept = 0.0
    detail = {
        "alpha_upper": float(alpha.hi),
        "phi_Sa_upper": float(phiS.hi),
        "S_gain_abs_upper": 0.0,
        "S_non_aw_error_part_upper_m_s": 0.0,
    }
    if phase == "due":
        pss, psa, _paw = RG._scalar_axis_structure(Pp)
        rs2 = src["R_S_filter_std"].square()
        den = pss + rs2
        if den.lo <= 0.0:
            raise RuntimeError("due-S innovation lost positive floor")
        ks = (psa / den).abs_upper()
        b = domain["startup"]["physical_handoff_coordinate_bounds"]
        v0 = float(b["velocity_error_norm_upper_mps"])
        S0 = float(b["integral_displacement_error_norm_upper_m_s"])
        const = up(S0 + up(h * position_norm_upper) + up(0.5 * h * h * v0))
        intercept = up(ks * const)
        slope = up(alpha.hi + up(ks * phiS.hi))
        detail.update({
            "S_gain_abs_upper": ks,
            "S_non_aw_error_part_upper_m_s": const,
        })
    if not (math.isfinite(intercept) and intercept >= 0.0 and math.isfinite(slope) and slope > 0.0):
        raise RuntimeError("invalid affine a_w prefix bound")
    return intercept, slope, detail


def _child_aw_radius_lower(*, k: float, rotational: float, bias: float,
                           intercept: float, slope: float) -> tuple[float, dict]:
    if not (k > 0.0 and slope > 0.0):
        raise RuntimeError("positive gain and affine slope required")
    allowed_total_residual = down(LIMIT / k)
    allowed_aw_after_prefix = _safe_sub_lower(
        _safe_sub_lower(allowed_total_residual, rotational), bias)
    numerator = _safe_sub_lower(allowed_aw_after_prefix, intercept)
    radius = 0.0 if numerator <= 0.0 else down(numerator / slope)
    radius = max(0.0, radius)
    return radius, {
        "allowed_total_residual_lower_mps2": allowed_total_residual,
        "allowed_aw_after_prefix_lower_mps2": allowed_aw_after_prefix,
        "affine_intercept_upper_mps2": intercept,
        "affine_slope_upper": slope,
        "aw_entrance_radius_lower_mps2": radius,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4 a_w budget domain must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("P4 a_w budget requires lever arm disabled")

    # The structured first-prefix proof requires the same source-preserving
    # covariance backend as the working candidate range producer.  Installing
    # only FULL3 loses exact axis symmetry through generic PSD boxing.
    RG3._install_backend(path, source_pieces)
    FULL3._install_backend()
    entrance = ENTRANCE.build(path)
    vector = VECTOR.build()
    failures = [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    src_phases = RG._source_phase_children(source_pieces)
    xcells = SG._linear_cells(alignment_pieces)
    live = domain["normal_live"]
    force_cells = RG._geom_ranges(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
        force_magnitude_pieces,
    )
    h = float(FULL._source_cell()["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    vc = vector["configured_measurement_bounds"]
    racc_var = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))[0][0]
    startup = domain["startup"]
    handoff = startup["physical_handoff_coordinate_bounds"]
    ba = float(handoff["accelerometer_bias_error_norm_upper_mps2"])
    p5_aw = float(handoff["latent_acceleration_error_norm_upper_mps2"])
    gravity = float(startup["gravity_mps2"])
    aw_fraction_g = float(startup.get("latent_acceleration_error_fraction_g", math.nan))
    pnorm = _p5_position_norm_upper(domain)

    rows = []
    for crow in entrance["P4_complete_word_search"]["candidate_rows"]:
        angle = float(crow["angle_deg"])
        q0 = float(crow["cayley_norm_upper"])
        qpred = RG._q_after_first_prediction(q0, domain, h)
        rot_gain = BASE._rotation_residual_gain_upper(qpred)
        budget = math.inf
        limiting = None
        total = 0

        for si, (src, phase) in enumerate(src_phases):
            P0 = FULL._initial_covariance(src, path)
            F, Q, _ = FULL._transition_and_Q(src, domain)
            Pp = FULL._psd_tighten(FULL.matrix_add(
                FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
            _pss, _psa, paw_pred = RG._scalar_axis_structure(Pp)
            if phase == "due":
                paw, _ignored_outer_error = RG._due_paw_and_error_norm(Pp, src, 0.0, 0.0)
            else:
                paw = paw_pred
            intercept, slope, sdetail = _s_phase_affine_aw_bound(
                src, phase, Pp, domain, pnorm)

            for xi, x in enumerate(xcells):
                for mi, m in enumerate(force_cells):
                    k, _kh, gdetail = V3._tangent_structured_gain_bounds(
                        tilt=tilt, yaw=yaw, eps=eps, x=x, m=m,
                        paw=paw, racc_var=racc_var)
                    rotational = up(rot_gain * m.hi)
                    rlower, budget_detail = _child_aw_radius_lower(
                        k=k, rotational=rotational, bias=ba,
                        intercept=intercept, slope=slope)
                    total += 1
                    if rlower < budget:
                        budget = rlower
                        limiting = {
                            "source_phase_cell": si,
                            "pseudo_phase": phase,
                            "alignment_cell": xi,
                            "alignment_x_tangent_yaw_fraction": x.as_list(),
                            "force_magnitude_cell": mi,
                            "force_magnitude_mps2": m.as_list(),
                            "tau_s": src["tau_s"].as_list(),
                            "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                            "R_S_filter_std": src["R_S_filter_std"].as_list(),
                            "Ktheta_norm_upper": k,
                            "finite_rotation_residual_norm_upper_mps2": rotational,
                            "S_phase_affine_aw_bound": sdetail,
                            "gain_detail": gdetail,
                            **budget_detail,
                        }

        if not math.isfinite(budget):
            budget = 0.0
            failures.append(f"{angle:g}deg: no finite a_w budget emitted")
        rows.append({
            "angle_deg": angle,
            "candidate_q_upper": q0,
            "post_prediction_q_upper": qpred,
            "evaluated_children": total,
            "derived_P4_aw_error_norm_upper_mps2_lower": budget,
            "P5_outer_aw_error_norm_upper_mps2": p5_aw,
            "declared_startup_aw_inside_derived_P4_radius": p5_aw <= budget,
            "finite_P5_aw_capture_required": budget < p5_aw,
            "positive_funnel_radius": budget > 0.0,
            "limiting_child": limiting,
        })

    positive = bool(rows) and all(r["positive_funnel_radius"] for r in rows)
    if not positive:
        failures.append("one or more P4 candidates has no positive derived a_w funnel")
    passed = not failures
    widest_containing_startup = next(
        (r["angle_deg"] for r in rows if r["declared_startup_aw_inside_derived_P4_radius"]), None)
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_DERIVED_LATENT_ACCELERATION_FUNNEL_BUDGET",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "derived_from_shipping_6rad_range_not_assumed": True,
        "deployed_correction_limit_rad": LIMIT,
        "deployed_correction_limit_increased": False,
        "source_preserving_rotation_gauge_v3_backend_used": True,
        "P5_position_half_Hs_norm_used_in_due_S_affine_bound": True,
        "P5_position_norm_upper_m": pnorm,
        "legacy_P1_position_20m_not_used_for_P4_budget": True,
        "P5_outer_aw_error_norm_upper_mps2": p5_aw,
        "declared_startup_gravity_mps2": gravity,
        "declared_startup_aw_error_fraction_g": aw_fraction_g,
        "declared_startup_aw_is_exactly_0p3g": p5_aw == 0.3 * gravity and aw_fraction_g == 0.3,
        "candidate_rows": rows,
        "widest_candidate_whose_derived_aw_radius_contains_declared_startup_deg": widest_containing_startup,
        "all_candidates_have_positive_derived_aw_funnel": positive,
        "P4_AW_FUNNEL_BUDGET_CERTIFICATE": "PASS" if passed else "NOT_ESTABLISHED",
        "P5_FINITE_AW_CAPTURE_ESTABLISHED_HERE": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "if the declared 0.3g startup a_w envelope is inside a candidate's derived radius, rerun that complete 18/21-state candidate word with the declared startup radius and the structured signed accelerometer/Joseph-reset map; otherwise prove source-complete finite P5 capture into the derived a_w radius before that word"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "derived_from_shipping_6rad_range_not_assumed",
        "source_preserving_rotation_gauge_v3_backend_used",
        "P5_position_half_Hs_norm_used_in_due_S_affine_bound",
        "legacy_P1_position_20m_not_used_for_P4_budget",
        "all_candidates_have_positive_derived_aw_funnel",
        "declared_startup_aw_is_exactly_0p3g",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "P5_FINITE_AW_CAPTURE_ESTABLISHED_HERE",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE", "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction range changed")
    if not (7.3 < float(d.get("P5_position_norm_upper_m", 0.0)) < 7.4):
        f.append("P5 half-Hs position norm was not used")
    rows = d.get("candidate_rows", [])
    if [float(r.get("angle_deg", -1.0)) for r in rows] != [30.0, 25.0, 20.0, 15.0]:
        f.append("candidate ladder changed")
    outer = float(d.get("P5_outer_aw_error_norm_upper_mps2", -1.0))
    gravity = float(d.get("declared_startup_gravity_mps2", -1.0))
    if not (gravity > 0.0 and outer == 0.3 * gravity):
        f.append("declared startup/P5 a_w radius is not frozen at 0.3g")
    if float(d.get("declared_startup_aw_error_fraction_g", -1.0)) != 0.3:
        f.append("declared startup a_w fraction of g changed")
    for r in rows:
        b = float(r.get("derived_P4_aw_error_norm_upper_mps2_lower", -1.0))
        if not (math.isfinite(b) and b > 0.0):
            f.append(f"{r.get('angle_deg')}deg: invalid derived a_w funnel radius")
        if r.get("declared_startup_aw_inside_derived_P4_radius") is not (outer <= b):
            f.append(f"{r.get('angle_deg')}deg: startup containment flag is inconsistent")
        if r.get("finite_P5_aw_capture_required") is not (b < outer):
            f.append(f"{r.get('angle_deg')}deg: P5 capture flag is inconsistent")
        if r.get("limiting_child") is None:
            f.append(f"{r.get('angle_deg')}deg: missing limiting child")
        if int(r.get("evaluated_children", 0)) <= 0:
            f.append(f"{r.get('angle_deg')}deg: no source children evaluated")
    if d.get("P4_AW_FUNNEL_BUDGET_CERTIFICATE") == "PASS" and f:
        f.append("PASS carries validation failures")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_pieces=a.source_pieces,
              alignment_pieces=a.alignment_pieces,
              force_magnitude_pieces=a.force_magnitude_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_AW_FUNNEL_BUDGET_CERTIFICATE"],
        "P5_outer_aw": d["P5_outer_aw_error_norm_upper_mps2"],
        "P5_outer_aw_g": d["declared_startup_aw_error_fraction_g"],
        "P5_position_norm": d["P5_position_norm_upper_m"],
        "widest_contains_startup": d["widest_candidate_whose_derived_aw_radius_contains_declared_startup_deg"],
        "rows": [{
            "angle_deg": r["angle_deg"],
            "qpred": r["post_prediction_q_upper"],
            "aw_funnel": r["derived_P4_aw_error_norm_upper_mps2_lower"],
            "startup_inside": r["declared_startup_aw_inside_derived_P4_radius"],
            "capture_required": r["finite_P5_aw_capture_required"],
            "limiting_phase": (r["limiting_child"] or {}).get("pseudo_phase"),
            "limiting_K": (r["limiting_child"] or {}).get("Ktheta_norm_upper"),
            "limiting_rotation": (r["limiting_child"] or {}).get("finite_rotation_residual_norm_upper_mps2"),
        } for r in d["candidate_rows"]],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
