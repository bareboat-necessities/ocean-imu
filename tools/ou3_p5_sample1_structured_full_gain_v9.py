#!/usr/bin/env python3
"""Directional-residual V9 closure of the analytic 1+2 sample-1 gain core.

V8 preserves the gain numerator/denominator dependency, but it still multiplies
max(||K_yz,x||,||K_x,yz||) by one undirected residual norm.  In the exact
canonical proof gauge the gain is block diagonal:

    r_x   -> delta_theta_yz   with gain k_perp,
    r_yz  -> delta_theta_x    with gain k_parallel.

The first SO(2)-gauged tangent correction fixes the orthogonal first tangent
residual row to zero in the ideal source-symmetric block.  Therefore that row's
physical a_w error is source-correlated with the first gravity rotational
residual and accelerometer bias, rather than an arbitrary component of the full
latent-error ball.  The deployed reset gauge preserves this perpendicular
component; the next body-frame rotation mixes the other components only through
its certified small off-diagonal entries.

Let rho bound the complete sample-1 residual norm and rho_x bound its x
component.  Since the two attitude output blocks are orthogonal, the exact
worst allocation of the residual ball satisfies, for k_perp>=k_parallel,

    ||dtheta||^2 <= k_parallel^2 rho^2
                    +(k_perp^2-k_parallel^2) min(rho_x,rho)^2.

If k_parallel>=k_perp the worst allocation is entirely in the y/z block and is
k_parallel*rho.  This avoids both multiplying the large scalar gain by the
wrong residual direction and double-counting the same residual norm in both
blocks.

The pre-first attitude PSD cross-axis remainder can make the nominal zero
perpendicular first residual nonzero; it remains an explicit later obligation.
The possible sample-1 S covariance contraction/update also remains explicit.
This producer is a fail-closed canonical core only and does not promote sample
1, q<8, P5, or N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_structured_full_gain_v4 as V4
import ou3_p5_sample1_structured_full_gain_v8 as V8

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN
SCHEMA=9
RANGE=9.0


def _two_block_correction_upper(k_perp:float,k_parallel:float,rho_x:float,rho:float)->float:
    if not all(math.isfinite(x) and x>=0.0 for x in (k_perp,k_parallel,rho_x,rho)):
        raise RuntimeError("finite nonnegative two-block bounds required")
    R=FULL.up(rho)
    X=min(FULL.up(rho_x),R)
    kp2=FULL.up(k_perp*k_perp)
    ka2=FULL.up(k_parallel*k_parallel)
    R2=FULL.up(R*R)
    if k_perp <= k_parallel:
        return FULL.up(k_parallel*R)
    X2=FULL.up(X*X)
    diff=FULL.up(max(0.0,FULL.up(kp2-FULL.down(k_parallel*k_parallel))))
    s=FULL.up(FULL.up(ka2*R2)+FULL.up(diff*X2))
    return FULL.up(math.sqrt(max(0.0,s)))


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
    if phase!="due": failures.append("V9 directional witness requires first due source cell")
    fr=first["source_cells"][source_cell_index]
    g=float(dom["startup"]["gravity_mps2"])
    ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    F,_Q,Rstep=FULL._transition_and_Q(src,dom)
    alpha=F[15][15]
    if not alpha.lo>0.0: failures.append("positive OU decay floor lost")
    alpha_hi=float(alpha.hi)
    chord0=V4._gravity_chord_from_cos(float(first["post_prediction_true_gravity_cosine_lower"]))
    pred_chord=V4._correction_chord_upper(float(first["first_prediction_transport_angle_upper_rad"]))
    first_perp_aw=FULL.up(FULL.up(g*chord0)+ba)
    off=max(Rstep[0][1].abs_upper(),Rstep[0][2].abs_upper())

    rows=[]; bad=None; worst=None; unclosed=0
    max_x=max_corr=max_old=max_kp=max_ka=max_rho=0.0
    for base in core["rows"]:
        detail=base["gain_detail"]
        kperp=float(detail["scalar_gain_yz_norm_upper"])
        kpar=float(detail["two_by_two_theta_x_gain_norm_upper"])
        rho=float(base["sample1_residual_norm_upper_mps2"])
        eaw1=float(base["post_prediction_aw_error_norm_upper_mps2"])
        fn=float(base["sample1_force_norm_upper_mps2"])
        d_hi=max(0.0,float(base["first_attitude_correction_rad"][1]))
        chord=min(2.0,FULL.up(chord0+FULL.up(V4._correction_chord_upper(d_hi)+pred_chord)))

        # The nominal orthogonal first residual row is zero.  Its physical aw
        # component is therefore bounded by first gravity rotation plus bias.
        # The next source-enclosed body-frame rotation can mix the remaining
        # post-first aw vector into x only through the two off-diagonal entries.
        post_aw=FULL.up(eaw1/FULL.down(float(alpha.lo)))
        mixed=FULL.up(FULL.up(math.sqrt(2.0)*off)*post_aw)
        eaw_x=FULL.up(alpha_hi*FULL.up(first_perp_aw+mixed))
        rot_x=FULL.up(fn*chord)
        rho_x_source=FULL.up(rot_x+FULL.up(eaw_x+ba))
        rho_x=min(rho,rho_x_source)

        corr=_two_block_correction_upper(kperp,kpar,rho_x,rho)
        old=float(base["correction_norm_upper_rad"])
        if corr>FULL.up(old+1e-12):
            failures.append("directional correction exceeded parent isotropic bound")
            break
        closed=math.isfinite(corr) and corr<RANGE
        row={
            "p_cell":base["p_cell"],
            "tangent_residual_cell":base["tangent_residual_cell"],
            "axial_residual_cell":base["axial_residual_cell"],
            "first_tangent_residual_magnitude_mps2":base["first_tangent_residual_magnitude_mps2"],
            "first_axial_residual_mps2":base["first_axial_residual_mps2"],
            "sample1_force_norm_upper_mps2":fn,
            "sample1_full_residual_norm_upper_mps2":rho,
            "sample1_source_correlated_x_residual_upper_mps2":rho_x,
            "sample1_source_correlated_x_residual_pre_total_cap_mps2":rho_x_source,
            "post_prediction_perpendicular_aw_error_upper_mps2":eaw_x,
            "body_rotation_offdiag_abs_upper":off,
            "Ktheta_perpendicular_block_upper":kperp,
            "Ktheta_parallel_block_upper":kpar,
            "parent_isotropic_correction_upper_rad":old,
            "directional_correction_norm_upper_rad":corr,
            "inside_9rad_range":closed,
        }
        rows.append(row)
        max_x=max(max_x,rho_x); max_corr=max(max_corr,corr); max_old=max(max_old,old)
        max_kp=max(max_kp,kperp); max_ka=max(max_ka,kpar); max_rho=max(max_rho,rho)
        if worst is None or corr>worst["directional_correction_norm_upper_rad"]: worst=row
        if not closed:
            unclosed+=1
            if bad is None: bad=row

    ok=bool(rows) and unclosed==0 and bad is None and not failures
    return {
        "schema":SCHEMA,
        "qualification":"OU3_P5_SAMPLE1_DIRECTIONAL_RESIDUAL_ANALYTIC_BLOCK_GAIN_V9",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "V8_positive_ratio_gain_parent_used":True,
        "analytic_one_plus_two_block_structure_retained":True,
        "canonical_first_perpendicular_residual_zero_used":True,
        "first_exact_gravity_chord_residual_used_in_perpendicular_component":True,
        "reset_gauge_preserves_perpendicular_aw_component":True,
        "one_step_body_rotation_aw_mixing_included":True,
        "orthogonal_two_block_residual_ball_optimization_used":True,
        "large_scalar_gain_multiplied_by_full_residual_norm":False,
        "same_full_residual_norm_double_counted_across_blocks":False,
        "temporal_force_slew_assumed":False,
        "first_attitude_PSD_cross_axis_remainder_included":False,
        "sample1_S_covariance_update_included":False,
        "sample1_S_attitude_injection_included":False,
        "complete_sample1_branch_closed_here":False,"q8_word_promoted_here":False,
        "whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,
        "evaluated_joint_cells":len(rows),"unclosed_joint_cells":unclosed,
        "body_rotation_offdiag_abs_upper":off,
        "first_perpendicular_aw_component_upper_mps2":first_perp_aw,
        "max_source_correlated_x_residual_upper_mps2":max_x,
        "max_full_residual_norm_upper_mps2":max_rho,
        "max_Ktheta_perpendicular_block_upper":max_kp,
        "max_Ktheta_parallel_block_upper":max_ka,
        "max_parent_isotropic_correction_upper_rad":max_old,
        "max_directional_correction_norm_upper_rad":max_corr,
        "first_unclosed_joint_cell":bad,"worst_directional_joint_cell":worst,
        "P5_SAMPLE1_DIRECTIONAL_RESIDUAL_BLOCK_GAIN_V9":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":(
            "ADD_FIRST_PSD_CROSS_AXIS_AND_SAMPLE1_S_BRANCH_THEN_SIGNED_CAYLEY_COMPOSE"
            if ok else "REFINE_DIRECTIONAL_X_RESIDUAL_WITH_FINITE_ATTITUDE_COMPONENT_GEOMETRY"
        ),
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in (
        "source_generated_not_trajectory_fit","V8_positive_ratio_gain_parent_used",
        "analytic_one_plus_two_block_structure_retained","canonical_first_perpendicular_residual_zero_used",
        "first_exact_gravity_chord_residual_used_in_perpendicular_component",
        "reset_gauge_preserves_perpendicular_aw_component","one_step_body_rotation_aw_mixing_included",
        "orthogonal_two_block_residual_ball_optimization_used",
    ):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in (
        "source_replay_used","filter_changed","large_scalar_gain_multiplied_by_full_residual_norm",
        "same_full_residual_norm_double_counted_across_blocks","temporal_force_slew_assumed",
        "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
        "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here","q8_word_promoted_here",
        "whole_word_promoted_here","N_H_words_set_here",
    ):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no V9 cells")
    for k in ("max_source_correlated_x_residual_upper_mps2","max_directional_correction_norm_upper_rad"):
        if not math.isfinite(float(d.get(k,math.nan))):f.append(f"nonfinite {k}")
    st=d.get("P5_SAMPLE1_DIRECTIONAL_RESIDUAL_BLOCK_GAIN_V9"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and (w is not None or int(d.get("unclosed_joint_cells",-1))!=0):f.append("PASS retains unclosed cell")
    if st=="NOT_ESTABLISHED" and w is None and not f:f.append("missing V9 witness")
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
        "status":d["P5_SAMPLE1_DIRECTIONAL_RESIDUAL_BLOCK_GAIN_V9"],"cells":d["evaluated_joint_cells"],
        "unclosed":d["unclosed_joint_cells"],"max_x_residual":d["max_source_correlated_x_residual_upper_mps2"],
        "max_full_residual":d["max_full_residual_norm_upper_mps2"],"max_k_perp":d["max_Ktheta_perpendicular_block_upper"],
        "max_k_parallel":d["max_Ktheta_parallel_block_upper"],"max_parent":d["max_parent_isotropic_correction_upper_rad"],
        "max_directional":d["max_directional_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],
        "worst":d["worst_directional_joint_cell"],"next":d["next_obligation"],"validation_failures":vf,
    },indent=2,sort_keys=True))
    return 0 if not vf else 2

if __name__=="__main__":raise SystemExit(main())
