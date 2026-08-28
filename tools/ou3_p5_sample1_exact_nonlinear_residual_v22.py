#!/usr/bin/env python3
"""V22: exact nonlinear accelerometer-residual audit at the first V18B q8 witness.

V21/V21B use the exact effective-input lemma but still enclose the nominal
sample-1 residual as the linear attitude term H_theta c plus four independent
norm charges: physical a_w error, accelerometer bias, the finite-angle attitude
defect, and the rotated-latent cross term.  Two pairs are algebraically the same
physical terms and need not be paid independently.

For Cayley vector c and the V7 sample-1 predicted force f, the exact rotational
residual is

    y_R(c,f) = (R(c)-I)f
             = [-2||c||^2 f + 2c(c^T f) + 4(c x f)]/(4+||c||^2).

Its small-angle term is exactly H_theta c = c x f.  Likewise the linear a_w
input plus its finite-angle latent cross term is just a rotated physical a_w
error, so orthogonality preserves its norm.  Therefore the nominal residual can
be enclosed as

    r_1 in y_R(c,f) + B(0, ||e_aw|| + ||b_a||),

without separately adding the attitude-defect and latent-cross norms.  V22 maps
that residual through the same signed V10 one-plus-two attitude gain, adds the
unchanged V12D PSD/S correction-perturbation ball, and intersects the result
with the exact V13E signed subcell at the authoritative V18B first witness.

The corrected current chart comes from V21B, which is bound to V14D's radial-
sinc quaternion semantics and reproduces q=0.6415230535178351.  The refined
correction is evaluated with the existing V16 axis cone, V15 geodesic route,
and V18 current-y/z support.  This is still a focused witness diagnostic only;
it promotes neither q8 nor sample 1 nor P5 and never sets N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_effective_input_correction_v21b as V21B
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D
import ou3_p5_sample1_signed_cayley_q8_v15 as V15
import ou3_p5_sample1_signed_cayley_q8_v16 as V16
import ou3_p5_sample1_signed_cayley_q8_v18 as V18

DEFAULT_DOMAIN = V21B.DEFAULT_DOMAIN
SCHEMA = 2200
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET


def _I(x) -> Interval:
    if isinstance(x, Interval):
        return x
    lo, hi = map(float, x)
    return Interval(lo, hi)


def _norm2_upper(a: Interval, b: Interval) -> float:
    aa = a.abs_upper(); bb = b.abs_upper()
    return FULL.up(math.sqrt(FULL.up(FULL.up(aa * aa) + FULL.up(bb * bb))))


def exact_rotation_residual(c, f):
    """Exact outward interval enclosure of (R(c)-I)f for Cayley [2,c]."""
    if len(c) != 3 or len(f) != 3:
        raise ValueError("three-component Cayley and force vectors required")
    cx, cy, cz = c
    fx, fy, fz = f
    q2 = cx.square() + cy.square() + cz.square()
    den = FULL.I(4.0) + q2
    if not den.lo > 0.0:
        raise RuntimeError("V22 Cayley denominator lost positivity")
    dot = cx * fx + cy * fy + cz * fz
    cross = [
        cy * fz - cz * fy,
        cz * fx - cx * fz,
        cx * fy - cy * fx,
    ]
    two = FULL.I(2.0); four = FULL.I(4.0)
    out = []
    for ci, fi, xi in zip(c, f, cross):
        num = -(two * q2 * fi) + two * ci * dot + four * xi
        out.append(num / den)
    return out


def _intersect_boxes(a, b):
    if len(a) != len(b):
        raise ValueError("equal correction dimensions required")
    out = []
    for x, y in zip(a, b):
        lo = max(x.lo, y.lo); hi = min(x.hi, y.hi)
        if hi < lo:
            return None
        out.append(Interval(lo, hi))
    return out


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    parent = V21B.build(
        Path(domain_path).resolve(),
        source_pieces=source_pieces,
        source_cell_index=source_cell_index,
        p_pieces=p_pieces,
        tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces,
        residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces,
    )
    failures = [f"V21B: {x}" for x in V21B.validate(parent)]
    if parent.get("P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B") != "PASS":
        failures.append("V21B authoritative witness prerequisite did not pass")

    q = float(parent.get("sample1_current_cayley_norm_upper", math.inf))
    c = [_I(x) for x in parent.get("sample1_current_component_box", [])]
    if len(c) != 3:
        failures.append("V21B current component box missing")
        c = [FULL.I(0.0), FULL.I(0.0), FULL.I(0.0)]
    fyz = parent.get("sample1_force_components_yz_mps2", [])
    if len(fyz) != 2:
        failures.append("V21B sample1 force components missing")
        fy, fz = FULL.I(0.0), FULL.I(0.0)
    else:
        fy, fz = _I(fyz[0]), _I(fyz[1])
    force = [FULL.I(0.0), fy, fz]

    gain = parent.get("gain_detail", {})
    perp = gain.get("perpendicular_gain_components", [])
    para = gain.get("parallel_gain_components", [])
    if len(perp) != 2 or len(para) != 2:
        failures.append("V21B signed one-plus-two gain components missing")
        gy = gz = kxy = kxz = FULL.I(0.0)
    else:
        gy, gz = _I(perp[0]), _I(perp[1])
        kxy, kxz = _I(para[0]), _I(para[1])

    yR = exact_rotation_residual(c, force)
    eaw = float(parent.get("post_prediction_physical_aw_error_norm_upper_mps2", math.inf))
    ba = float(parent.get("accelerometer_bias_error_norm_upper_mps2", math.inf))
    if not (math.isfinite(eaw) and eaw >= 0.0 and math.isfinite(ba) and ba >= 0.0):
        failures.append("invalid physical a_w/bias nuisance bounds")
        nuisance = math.inf
    else:
        nuisance = FULL.up(eaw + ba)

    corr_perturb = float(parent.get("V12D_correction_perturbation_norm_upper_rad", math.inf))
    if not (math.isfinite(corr_perturb) and corr_perturb >= 0.0):
        failures.append("invalid V12D correction perturbation")
        corr_perturb = math.inf

    if math.isfinite(nuisance) and math.isfinite(corr_perturb):
        n = Interval.outward_bounds(-nuisance, nuisance)
        ec = Interval.outward_bounds(-corr_perturb, corr_perturb)
        kx_norm = _norm2_upper(kxy, kxz)
        enx = Interval.outward_bounds(-FULL.up(kx_norm * nuisance),
                                      FULL.up(kx_norm * nuisance))
        source_box = [
            kxy * yR[1] + kxz * yR[2] + enx + ec,
            gy * yR[0] + gy * n + ec,
            gz * yR[0] + gz * n + ec,
        ]
    else:
        source_box = [Interval(-math.inf, math.inf) for _ in range(3)]

    baseline_box = [_I(x) for x in parent.get("baseline_correction_box_rad", [])]
    if len(baseline_box) != 3:
        failures.append("V21B baseline correction box missing")
        baseline_box = source_box
    joint = _intersect_boxes(baseline_box, source_box)
    incompatible = joint is None
    baseline_lo = float(parent.get("baseline_correction_radial_lower_rad", 0.0))
    baseline_hi = float(parent.get("baseline_correction_radial_upper_rad", math.inf))

    if incompatible:
        refined_lo = refined_hi = 0.0
        geo_q = product_q = 0.0
        product_w = math.inf
        axis_narrowed = False
        branches = ["SOURCE_INCOMPATIBLE"]
        closes = True
    else:
        box_hi = V14.CAYLEY1._norm_upper(joint)
        box_lo = V14.CAYLEY2._norm_lower(joint)
        refined_hi = min(baseline_hi, box_hi)
        refined_lo = min(refined_hi, max(baseline_lo, box_lo))
        geo = V15._geodesic_q_and_scalar_lower(q, refined_lo, refined_hi)
        geo_q = math.inf if geo is None else float(geo[0])

        wd, vd, branches, axis_narrowed = V16.axis_cone_normalized_shipping_quaternion(
            joint,
            radial_lower=refined_lo,
            radial_upper=refined_hi,
            parent=V14D.radial_sinc_normalized_shipping_quaternion,
        )
        cyz = min(q, V18._yz_norm_upper(c[1], c[2]))
        chart = {"cx": c[0], "cy": c[1], "cz": c[2], "cyz_norm_upper": cyz}
        parent_W = FULL.I(2.0) * wd - V14.CAYLEY1._dot(vd, c)
        W, _yz_box, _yz_joint = V18._support_product_scalar(
            parent_W, wd, vd, chart)
        product_w, product_q = V14._qplus_from_product_scalar(q, W)
        closes = ((math.isfinite(geo_q) and geo_q < Q_TARGET)
                  or (math.isfinite(product_q) and product_q < Q_TARGET
                      and product_w > 0.0))

    old_nuisance = float(parent.get("nominal_effective_residual_nuisance_norm_upper_mps2", math.inf))
    nuisance_strict = math.isfinite(old_nuisance) and nuisance < old_nuisance
    radial_strict = incompatible or refined_hi < baseline_hi
    next_obligation = (
        "LIFT_V22_EXACT_NONLINEAR_RESIDUAL_INTERSECTION_INTO_FULL_V18B_Q8_COVER"
        if closes else
        "SUBDIVIDE_AUTHORITATIVE_CURRENT_CHART_WITH_EXACT_NONLINEAR_RESIDUAL_AT_FIRST_Q8_WITNESS"
        if radial_strict else
        "DERIVE_SOURCE_CORRELATED_SAMPLE1_AW_ERROR_COMPONENTS_AT_FIRST_Q8_WITNESS"
    )

    status = "PASS" if not failures else "NOT_ESTABLISHED"
    return {
        **parent,
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22",
        "V21B_authoritative_current_chart_retained": True,
        "exact_Cayley_rotation_residual_used": True,
        "linear_attitude_plus_defect_double_charge_retired": True,
        "linear_aw_plus_latent_cross_double_charge_retired": True,
        "rotated_physical_aw_norm_preserved_by_orthogonality": True,
        "accelerometer_bias_norm_retained": True,
        "V12D_PSD_S_correction_perturbation_retained": True,
        "V16_axis_cone_and_V18_yz_support_used_for_product_check": True,
        "exact_rotation_residual_box_mps2": [x.as_list() for x in yR],
        "combined_rotated_aw_plus_bias_nuisance_norm_upper_mps2": nuisance,
        "V21B_previous_nuisance_norm_upper_mps2": old_nuisance,
        "nominal_nuisance_strictly_reduced": nuisance_strict,
        "exact_residual_source_correction_box_rad": [x.as_list() for x in source_box],
        "joint_exact_residual_correction_box_rad": None if joint is None else [x.as_list() for x in joint],
        "source_subcell_incompatible_under_exact_residual": incompatible,
        "refined_correction_radial_lower_rad": refined_lo,
        "refined_correction_radial_upper_rad": refined_hi,
        "correction_radial_strictly_refined_by_exact_residual": radial_strict,
        "refined_geodesic_q_upper": geo_q,
        "refined_product_abs_W_lower": product_w,
        "refined_product_q_upper": product_q,
        "V16_axis_cone_narrowed_exact_residual_quaternion": axis_narrowed,
        "exact_residual_quaternion_branches": branches,
        "first_witness_closed_inside_q8": closes,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22": status,
        "next_obligation": next_obligation,
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22":
        f.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V21B_authoritative_current_chart_retained",
        "exact_Cayley_rotation_residual_used",
        "linear_attitude_plus_defect_double_charge_retired",
        "linear_aw_plus_latent_cross_double_charge_retired",
        "rotated_physical_aw_norm_preserved_by_orthogonality",
        "accelerometer_bias_norm_retained",
        "V12D_PSD_S_correction_perturbation_retained",
        "V16_axis_cone_and_V18_yz_support_used_for_product_check",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    new_n = d.get("combined_rotated_aw_plus_bias_nuisance_norm_upper_mps2")
    old_n = d.get("V21B_previous_nuisance_norm_upper_mps2")
    if not all(isinstance(x, (int, float)) and math.isfinite(float(x)) and float(x) >= 0.0
               for x in (new_n, old_n)):
        f.append("invalid V22 nuisance accounting")
    elif float(new_n) > float(old_n):
        f.append("exact residual algebra increased nuisance norm")
    if not d.get("source_subcell_incompatible_under_exact_residual"):
        if float(d.get("refined_correction_radial_upper_rad", math.inf)) > \
                float(d.get("baseline_correction_radial_upper_rad", -math.inf)):
            f.append("V22 exact residual intersection worsened radial upper")
    st = d.get("P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22")
    if st not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V22 status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=6)
    ap.add_argument("--parallel-pieces", type=int, default=6)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(
        x.domain,
        source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index,
        p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces,
        axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces,
    )
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22"],
        "q_current": d.get("sample1_current_cayley_norm_upper"),
        "old_nuisance": d.get("V21B_previous_nuisance_norm_upper_mps2"),
        "new_nuisance": d.get("combined_rotated_aw_plus_bias_nuisance_norm_upper_mps2"),
        "baseline_radial_upper": d.get("baseline_correction_radial_upper_rad"),
        "refined_radial_upper": d.get("refined_correction_radial_upper_rad"),
        "refined_geodesic_q": d.get("refined_geodesic_q_upper"),
        "refined_product_W": d.get("refined_product_abs_W_lower"),
        "refined_product_q": d.get("refined_product_q_upper"),
        "axis_cone_narrowed": d.get("V16_axis_cone_narrowed_exact_residual_quaternion"),
        "closed_q8": d.get("first_witness_closed_inside_q8"),
        "source_incompatible": d.get("source_subcell_incompatible_under_exact_residual"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
