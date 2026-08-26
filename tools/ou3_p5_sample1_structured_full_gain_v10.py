#!/usr/bin/env python3
"""Combined perpendicular-residual transport for OU-III P5 sample 1.

V9 correctly uses the exact 1+2 gain blocks, but its x-residual cap still
triangles the finite-angle force term and a raw a_w component.  That misses an
exact source correlation and, more importantly, does not carry the latent term
in its rotated form.

Work in V7's un-Rx proof gauge.  At the first accelerometer use SO(2) symmetry to
place the useful tangent residual in the y row, so the perpendicular x row is
exactly zero.  In the ideal source-symmetric covariance block the first a_w
correction also has dx_aw,x=0 and the attitude correction is about x.  The
shipping correction is injected on the left; the simultaneous Rx(d)^T proof
gauge therefore cancels the ideal axis-angle injection and preserves the x row
of the physical attitude-error rotation.

Let E be that physical error rotation, f0=g e3 in the canonical gauge, e0 the
physical a_w error immediately before the first accelerometer, b0 the additive
accelerometer-bias error, and Delta the first a_w estimator correction.  The
perpendicular first residual is

    0 = [(E-I)f0]_x + [E e0]_x + b0_x.                 (1)

After the update and one homogeneous OU prediction the nominal force and
physical a_w error contain the same Delta with opposite signs:

    f1 = f0 + alpha Delta,
    e1 = alpha (e0-Delta).

Because Delta_x=0, their sample-1 x residual before the tiny attitude-error
transport is exactly

    [(E-I)f1]_x + [E e1]_x + b1_x
      = (1-alpha)[(E-I)f0]_x + (b1_x-alpha b0_x).      (2)

Thus the potentially large first tangent/axial a_w correction cancels exactly;
no force-slew assumption is used.  Equation (2) also carries the latent term as
E e rather than as an unrotated component.

Two small source-faithful remainders are added:

* the 5 ms attitude-error transport from gyro-bias/deterministic transport;
* below the source's 1e-2-rad quaternion threshold, the normalized polynomial
  correction differs slightly from the exact axis-angle Rx(d) used by the proof
  gauge.  The mismatch is bounded from the already validated deployed source
  branch; above the threshold it is exactly zero.

The pre-first attitude PSD cross-axis covariance remainder can spoil the ideal
correction-axis/Delta_x=0 structure and is deliberately NOT included here.
Likewise the possible sample-1 S covariance/update remains separate.  This is a
fail-closed canonical core only and does not promote sample 1, q<8, P5, or
N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_exact_correction_transport as CORR
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_structured_full_gain_v4 as V4
import ou3_p5_sample1_structured_full_gain_v8 as V8
import ou3_p5_sample1_structured_full_gain_v9 as V9

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN
SCHEMA=10
RANGE=9.0
SERIES_THRESHOLD=1.0e-2


def _series_vs_axis_rotation_mismatch_upper(d_lo:float,d_hi:float)->float:
    """Operator-norm bound between deployed series injection and ideal Rx(d)."""
    lo=max(0.0,float(d_lo)); hi=max(lo,float(d_hi))
    if lo>=SERIES_THRESHOLD or hi<=0.0:
        return 0.0
    delta=min(hi,SERIES_THRESHOLD)
    _aser,series_from_d,_=CORR._series_cayley_scalar_bounds(delta)
    _aaxis,axis_from_d=CORR._axis_cayley_scalar_bounds(delta)
    # For a common rotation axis theta(a)=2 atan(a/2) has derivative <=1,
    # hence the axis-angle mismatch is no larger than the Cayley mismatch.
    angle=FULL.up(series_from_d+axis_from_d)
    # ||R1-R2||_2 = 2 sin(|Delta theta|/2) <= |Delta theta|.
    return min(2.0,angle)


def _combined_x_residual_upper(*,alpha_lo:float,alpha_hi:float,
                               first_rot_x_upper:float,bias_upper:float,
                               error_transport_rotation_norm_upper:float,
                               series_rotation_mismatch_upper:float,
                               pre_first_aw_error_norm_upper:float,
                               gravity:float)->dict:
    vals=(alpha_lo,alpha_hi,first_rot_x_upper,bias_upper,
          error_transport_rotation_norm_upper,series_rotation_mismatch_upper,
          pre_first_aw_error_norm_upper,gravity)
    if not all(math.isfinite(float(x)) for x in vals):
        raise RuntimeError("finite combined-residual bounds required")
    if not (0.0<alpha_lo<=alpha_hi<=1.0 and min(vals[2:])>=0.0):
        raise RuntimeError("invalid combined-residual signs/ranges")

    decay=FULL.up(FULL.up(1.0-alpha_lo)*first_rot_x_upper)
    bias=FULL.up(FULL.up(1.0+alpha_hi)*bias_upper)
    mismatch=min(2.0,FULL.up(error_transport_rotation_norm_upper+series_rotation_mismatch_upper))
    vector=FULL.up(gravity+FULL.up(alpha_hi*pre_first_aw_error_norm_upper))
    geom=FULL.up(mismatch*vector)
    total=FULL.up(decay+FULL.up(bias+geom))
    return {
        "ou_decay_times_first_rotational_residual_upper_mps2":decay,
        "bias_difference_upper_mps2":bias,
        "combined_rotation_mismatch_operator_norm_upper":mismatch,
        "rotation_mismatch_vector_norm_upper_mps2":vector,
        "rotation_mismatch_residual_upper_mps2":geom,
        "combined_x_residual_upper_mps2":total,
    }


def build(domain_path:Path=DEFAULT_DOMAIN,*,source_pieces:int=4,source_cell_index:int=0,
          p_pieces:int=24,tangent_pieces:int=24,axial_pieces:int=24)->dict:
    FULL3._install_backend()
    path=Path(domain_path).resolve()
    dom=json.loads(path.read_text(encoding="utf-8"))
    core=V8.build(path,source_pieces=source_pieces,source_cell_index=source_cell_index,
                  p_pieces=p_pieces,tangent_pieces=tangent_pieces,axial_pieces=axial_pieces)
    first=FIRST.build(path,source_pieces=source_pieces)
    failures=[f"V8: {x}" for x in V8.validate(core)]
    failures += [f"first: {x}" for x in FIRST.validate(first)]

    src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
    if phase!="due": failures.append("V10 combined-residual witness requires first due source cell")
    fr=first["source_cells"][source_cell_index]
    h=float(src["dt_s"])
    g=float(dom["startup"]["gravity_mps2"])
    ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    F,_Q,_Rstep=FULL._transition_and_Q(src,dom)
    alpha=F[15][15]
    alpha_lo=float(alpha.lo); alpha_hi=float(alpha.hi)
    aw_pre=float(fr["predicted_aw_error_norm_upper_mps2"])

    cos0=float(first["post_prediction_true_gravity_cosine_lower"])
    sin0=1.0 if cos0<0.0 else FULL.up(math.sqrt(max(0.0,FULL.up(1.0-FULL.down(cos0*cos0)))))
    first_rot_x=FULL.up(g*sin0)

    # The same h*(bias+deterministic transport) bound applies to every 5 ms
    # attitude-error prediction.  Nominal body rotation is removed by V7's
    # simultaneous orthogonal gauge.
    theta_transport=float(first["first_prediction_transport_angle_upper_rad"])
    transport_chord=V4._correction_chord_upper(theta_transport)

    rows=[]; bad=None; worst=None; unclosed=0
    max_x=max_corr=max_parent=max_series=max_rho=max_kp=max_ka=0.0
    for base in core["rows"]:
        detail=base["gain_detail"]
        kperp=float(detail["scalar_gain_yz_norm_upper"])
        kpar=float(detail["two_by_two_theta_x_gain_norm_upper"])
        rho=float(base["sample1_residual_norm_upper_mps2"])
        dlo,dhi=map(float,base["first_attitude_correction_rad"])
        series=_series_vs_axis_rotation_mismatch_upper(dlo,dhi)
        comb=_combined_x_residual_upper(
            alpha_lo=alpha_lo,alpha_hi=alpha_hi,
            first_rot_x_upper=first_rot_x,bias_upper=ba,
            error_transport_rotation_norm_upper=transport_chord,
            series_rotation_mismatch_upper=series,
            pre_first_aw_error_norm_upper=aw_pre,gravity=g,
        )
        rho_x=min(rho,float(comb["combined_x_residual_upper_mps2"]))
        corr=V9._two_block_correction_upper(kperp,kpar,rho_x,rho)
        parent=float(base["correction_norm_upper_rad"])
        if corr>FULL.up(parent+1e-12):
            failures.append("V10 directional correction exceeded V8 isotropic parent")
            break
        closed=math.isfinite(corr) and corr<RANGE
        row={
            "p_cell":base["p_cell"],
            "tangent_residual_cell":base["tangent_residual_cell"],
            "axial_residual_cell":base["axial_residual_cell"],
            "first_tangent_residual_magnitude_mps2":base["first_tangent_residual_magnitude_mps2"],
            "first_axial_residual_mps2":base["first_axial_residual_mps2"],
            "first_attitude_correction_rad":base["first_attitude_correction_rad"],
            "sample1_force_norm_upper_mps2":base["sample1_force_norm_upper_mps2"],
            "sample1_full_residual_norm_upper_mps2":rho,
            "sample1_combined_source_x_residual_upper_mps2":rho_x,
            "Ktheta_perpendicular_block_upper":kperp,
            "Ktheta_parallel_block_upper":kpar,
            "series_vs_axis_rotation_mismatch_upper":series,
            **comb,
            "parent_isotropic_correction_upper_rad":parent,
            "combined_directional_correction_norm_upper_rad":corr,
            "inside_9rad_range":closed,
        }
        rows.append(row)
        max_x=max(max_x,rho_x); max_corr=max(max_corr,corr); max_parent=max(max_parent,parent)
        max_series=max(max_series,series); max_rho=max(max_rho,rho); max_kp=max(max_kp,kperp); max_ka=max(max_ka,kpar)
        if worst is None or corr>worst["combined_directional_correction_norm_upper_rad"]: worst=row
        if not closed:
            unclosed+=1
            if bad is None: bad=row

    ok=bool(rows) and unclosed==0 and bad is None and not failures
    return {
        "schema":SCHEMA,
        "qualification":"OU3_P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_TRANSPORT_V10",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "V8_positive_ratio_gain_parent_used":True,
        "analytic_one_plus_two_block_structure_retained":True,
        "first_perpendicular_residual_exact_zero_in_ideal_SO2_gauge":True,
        "first_perpendicular_aw_estimator_correction_exact_zero_in_ideal_block":True,
        "latent_term_carried_as_rotated_E_eaw_component":True,
        "first_force_and_post_update_aw_correction_cancel_in_sample1_x_residual":True,
        "combined_x_residual_identity":"r1x=(1-alpha)yR0x+(b1x-alpha*b0x)+transport/series remainder",
        "nominal_body_rotation_removed_by_V7_simultaneous_orthogonal_gauge":True,
        "same_one_step_error_transport_bound_reused_for_sample1":True,
        "deployed_series_vs_axis_gauge_mismatch_included":True,
        "above_series_threshold_axis_gauge_mismatch_exact_zero":True,
        "orthogonal_two_block_residual_ball_optimization_used":True,
        "temporal_force_slew_assumed":False,
        "raw_aw_component_used_for_x_latent_residual":False,
        "large_scalar_gain_multiplied_by_full_residual_norm":False,
        "first_attitude_PSD_cross_axis_remainder_included":False,
        "sample1_S_covariance_update_included":False,
        "sample1_S_attitude_injection_included":False,
        "complete_sample1_branch_closed_here":False,"q8_word_promoted_here":False,
        "whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,
        "alpha_interval":alpha.as_list(),
        "first_gravity_x_rotational_residual_upper_mps2":first_rot_x,
        "pre_first_aw_error_norm_upper_mps2":aw_pre,
        "one_step_error_transport_angle_upper_rad":theta_transport,
        "one_step_error_transport_rotation_norm_upper":transport_chord,
        "max_series_vs_axis_rotation_mismatch_upper":max_series,
        "evaluated_joint_cells":len(rows),"unclosed_joint_cells":unclosed,
        "max_combined_source_x_residual_upper_mps2":max_x,
        "max_full_residual_norm_upper_mps2":max_rho,
        "max_Ktheta_perpendicular_block_upper":max_kp,
        "max_Ktheta_parallel_block_upper":max_ka,
        "max_parent_isotropic_correction_upper_rad":max_parent,
        "max_combined_directional_correction_norm_upper_rad":max_corr,
        "first_unclosed_joint_cell":bad,"worst_directional_joint_cell":worst,
        "P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":(
            "ADD_FIRST_PSD_CROSS_AXIS_AND_SAMPLE1_S_BRANCH_THEN_SIGNED_CAYLEY_COMPOSE"
            if ok else "REFINE_COMBINED_X_RESIDUAL_TRANSPORT_AT_FIRST_WITNESS"
        ),
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in (
        "source_generated_not_trajectory_fit","V8_positive_ratio_gain_parent_used",
        "analytic_one_plus_two_block_structure_retained","first_perpendicular_residual_exact_zero_in_ideal_SO2_gauge",
        "first_perpendicular_aw_estimator_correction_exact_zero_in_ideal_block",
        "latent_term_carried_as_rotated_E_eaw_component",
        "first_force_and_post_update_aw_correction_cancel_in_sample1_x_residual",
        "nominal_body_rotation_removed_by_V7_simultaneous_orthogonal_gauge",
        "same_one_step_error_transport_bound_reused_for_sample1",
        "deployed_series_vs_axis_gauge_mismatch_included","above_series_threshold_axis_gauge_mismatch_exact_zero",
        "orthogonal_two_block_residual_ball_optimization_used",
    ):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in (
        "source_replay_used","filter_changed","temporal_force_slew_assumed","raw_aw_component_used_for_x_latent_residual",
        "large_scalar_gain_multiplied_by_full_residual_norm","first_attitude_PSD_cross_axis_remainder_included",
        "sample1_S_covariance_update_included","sample1_S_attitude_injection_included",
        "complete_sample1_branch_closed_here","q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
    ):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no V10 cells")
    for k in ("max_combined_source_x_residual_upper_mps2","max_combined_directional_correction_norm_upper_rad"):
        if not math.isfinite(float(d.get(k,math.nan))):f.append(f"nonfinite {k}")
    st=d.get("P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and (w is not None or int(d.get("unclosed_joint_cells",-1))!=0):f.append("PASS retains unclosed V10 cell")
    if st=="NOT_ESTABLISHED" and w is None and not f:f.append("missing V10 witness")
    return f


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces",type=int,default=4); ap.add_argument("--source-cell-index",type=int,default=0)
    ap.add_argument("--p-pieces",type=int,default=24); ap.add_argument("--tangent-pieces",type=int,default=24)
    ap.add_argument("--axial-pieces",type=int,default=24); ap.add_argument("--output",type=Path,required=True)
    x=ap.parse_args(); d=build(x.domain,source_pieces=x.source_pieces,source_cell_index=x.source_cell_index,
        p_pieces=x.p_pieces,tangent_pieces=x.tangent_pieces,axial_pieces=x.axial_pieces)
    vf=validate(d); d["validation_failures"]=vf
    x.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({
        "status":d["P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10"],"cells":d["evaluated_joint_cells"],
        "unclosed":d["unclosed_joint_cells"],"alpha":d["alpha_interval"],
        "first_rot_x":d["first_gravity_x_rotational_residual_upper_mps2"],
        "transport_chord":d["one_step_error_transport_rotation_norm_upper"],
        "series_mismatch":d["max_series_vs_axis_rotation_mismatch_upper"],
        "max_x_residual":d["max_combined_source_x_residual_upper_mps2"],
        "max_full_residual":d["max_full_residual_norm_upper_mps2"],
        "max_k_perp":d["max_Ktheta_perpendicular_block_upper"],
        "max_k_parallel":d["max_Ktheta_parallel_block_upper"],
        "max_parent":d["max_parent_isotropic_correction_upper_rad"],
        "max_directional":d["max_combined_directional_correction_norm_upper_rad"],
        "first_unclosed":d["first_unclosed_joint_cell"],"worst":d["worst_directional_joint_cell"],
        "next":d["next_obligation"],"validation_failures":vf,
    },indent=2,sort_keys=True))
    return 0 if not vf else 2

if __name__=="__main__":raise SystemExit(main())
