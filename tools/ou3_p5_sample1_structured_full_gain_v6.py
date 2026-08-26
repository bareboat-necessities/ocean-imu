#!/usr/bin/env python3
"""Direct first-residual-coordinate V6 sample-1 accelerometer core.

V5 still subdivided the first axial estimator correction and inverted an
interval gain to recover r_z.  That reintroduced dependency loss: a retained
child could report |r_z| larger than the already certified ||r0|| bound.  V6
uses the first residual itself as the source coordinate.

After SO(2) gravity-gauge reduction write the useful first residual as

    r0 = [r_t, 0, r_z],  r_t >= 0,  r_t^2+r_z^2 <= rho0^2.

For the ideal canonical first covariance block,

    d       = k_theta r_t,
    da_w,t  = k_aw,t r_t,
    da_w,z  = k_z r_z,

with all gains evaluated forward from the same p,t,R cell.  No gain is inverted.
The first residual ball is enforced directly, and source directional bounds are
also used.  From the startup gravity cosine lower c_g,

    ||y_R,t|| <= g sqrt(1-c_g^2),
    |y_R,z|   <= g (1-c_g).

Since the latent contribution and accelerometer bias have norms bounded by the
existing handoff contract, these give source-uniform component caps for r_t and
r_z without narrowing the theorem domain.

The physical post-first-update aw error uses the exact residual-coordinate
identity from V5,

    e_aw^+ = (I-D_aw) r0 - y_R - b_a,

with tangent/axial geometry retained before taking the norm.  The older
prior-plus-correction triangle remains a second independent upper bound.  The
complete V4 structured sample-1 nonaxial force, 3x3 innovation, and 3x3 attitude
gain are then evaluated.

The pre-first attitude PSD cross-axis remainder and sample-1 S covariance/update
remain explicit later obligations.  This producer is a fail-closed canonical
core and does not promote sample 1, q<8, the whole word, or N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_p5_sample1_structured_full_gain_v4 as V4
import ou3_p5_sample1_structured_full_gain_v5 as V5
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 6
RANGE = 9.0


def _residual_child(rt0: Interval, rz0: Interval, rho0: float):
    rt_lo=max(0.0,rt0.lo); rt_hi=max(rt_lo,rt0.hi)
    rz_min=V5._min_abs(rz0)
    if V5._norm2_lower(rt_lo,rz_min)>rho0:
        return None
    rt_cap_sq=max(0.0,FULL.up(FULL.up(rho0*rho0)-FULL.down(rz_min*rz_min)))
    rt_cap=FULL.up(math.sqrt(rt_cap_sq))
    rt_hi=min(rt_hi,rt_cap)
    if rt_hi<rt_lo:
        return None
    rt=Interval(rt_lo,rt_hi)
    rt_min=max(0.0,rt.lo)
    rz_cap_sq=max(0.0,FULL.up(FULL.up(rho0*rho0)-FULL.down(rt_min*rt_min)))
    rz_cap=FULL.up(math.sqrt(rz_cap_sq))
    rz_lo=max(rz0.lo,-rz_cap); rz_hi=min(rz0.hi,rz_cap)
    if rz_hi<rz_lo:
        return None
    rz=Interval(rz_lo,rz_hi)
    if V5._norm2_lower(max(0.0,rt.lo),V5._min_abs(rz))>rho0:
        return None
    return rt,rz


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    FULL3._install_backend()
    path=Path(domain_path).resolve()
    dom=json.loads(path.read_text(encoding="utf-8"))
    first=FIRST.build(path,source_pieces=source_pieces)
    vec=VECTOR.build()
    failures=[f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vec)]

    src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
    if phase!="due":
        failures.append("V6 direct-residual witness requires first due source cell")
    fr=first["source_cells"][source_cell_index]
    p_all=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"]))
    rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
    aw_pred=float(fr["predicted_aw_error_norm_upper_mps2"])

    h=float(src["dt_s"])
    g=float(dom["startup"]["gravity_mps2"])
    ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    tilt,yaw,eps=RG._attitude_covariance_epsilon(path,h)
    t=Interval.outward_bounds(tilt,FULL.up(tilt+eps))
    Y=Interval.outward_bounds(yaw,FULL.up(yaw+eps))
    Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))
    r=Racc[0][0]
    F,Q,_=FULL._transition_and_Q(src,dom)
    alpha=F[15][15]; alpha_hi=alpha.hi; qaw=Q[15][15]

    cos0=float(first["post_prediction_true_gravity_cosine_lower"])
    if not (-1.0 <= cos0 <= 1.0):
        failures.append("invalid first gravity cosine lower")
    sin_hi=1.0 if cos0<0.0 else FULL.up(math.sqrt(max(0.0,FULL.up(1.0-FULL.down(cos0*cos0)))))
    yRt=FULL.up(g*sin_hi)
    yRz=FULL.up(g*max(0.0,FULL.up(1.0-cos0)))
    yRnorm=FULL.up(g*math.sqrt(max(0.0,FULL.up(2.0*FULL.up(1.0-cos0)))))
    rt_source_cap=min(rho0,FULL.up(aw_pred+FULL.up(yRt+ba)))
    rz_source_cap=min(rho0,FULL.up(aw_pred+FULL.up(yRz+ba)))

    chord0=V4._gravity_chord_from_cos(cos0)
    pred_chord=V4._correction_chord_upper(float(first["first_prediction_transport_angle_upper_rad"]))

    pcells=SUB.parts(p_all.lo,p_all.hi,p_pieces)
    rtcells=SUB.parts(0.0,rt_source_cap,tangent_pieces)
    rzcells=SUB.parts(-rz_source_cap,rz_source_cap,axial_pieces)
    rows=[]; bad=None; fixed=fallback=0; pruned=0
    maxK=maxD=maxRho=maxF=maxEaw=maxEawPost=maxChord=maxD0=maxAz=0.0
    minS=math.inf

    for pi,p in enumerate(pcells):
        D=FULL.I(g*g)*t+p+r
        a=t*(p+r)/D
        c0=-(FULL.I(g)*t*p/D)
        b=p*(FULL.I(g*g)*t+r)/D
        bz=p*r/(p+r)
        Pth0=V4._diag3(a,a,Y)
        Paw0=V4._diag3(b,b,bz)
        C0=FULL._zero(3,3); C0[0][1]=-c0; C0[1][0]=c0
        ktheta=FULL.I(g)*t/D
        kaw_t=p/D
        kz=p/(p+r)
        if not (ktheta.lo>0.0 and kaw_t.lo>=0.0 and kz.lo>0.0):
            failures.append(f"first canonical gain lost positivity p={pi}")
            continue

        for ti,rt0 in enumerate(rtcells):
            for zi,rz0 in enumerate(rzcells):
                child=_residual_child(rt0,rz0,rho0)
                if child is None:
                    pruned += 1
                    continue
                rt,rz=child
                d=ktheta*rt
                awt=kaw_t*rt
                az=kz*rz
                maxD0=max(maxD0,d.abs_upper()); maxAz=max(maxAz,az.abs_upper())

                L,Rx=V4._Ltheta(d)
                Pth=matrix_mul(matrix_mul(L,Pth0),matrix_transpose(L))
                Paw_r=matrix_mul(matrix_mul(Rx,Paw0),matrix_transpose(Rx))
                C_r=matrix_mul(matrix_mul(L,C0),matrix_transpose(Rx))
                Paw=[[alpha.square()*Paw_r[i][j]+(qaw if i==j else FULL.I(0.0)) for j in range(3)] for i in range(3)]
                C=[[alpha*C_r[i][j] for j in range(3)] for i in range(3)]

                d_hi=max(0.0,d.hi)
                chord=min(2.0,FULL.up(chord0+FULL.up(V4._correction_chord_upper(d_hi)+pred_chord)))
                maxChord=max(maxChord,chord)

                raw=[FULL.I(0.0),-(alpha*awt),FULL.I(g)+alpha*az]
                f=FULL._mat_vec(Rx,raw)
                fn=V4._norm3_upper(f); maxF=max(maxF,fn)
                Ht=V4._Htheta(f)
                PHt=matrix_add(matrix_mul(Pth,matrix_transpose(Ht)),C)
                HP=matrix_mul(Ht,Pth); HC=matrix_mul(Ht,C)
                S=matrix_add(matrix_mul(HP,matrix_transpose(Ht)),HC)
                S=matrix_add(S,matrix_transpose(HC)); S=matrix_add(S,Paw)
                S=FULL.matrix_symmetric_hull(matrix_add(S,Racc))
                Sinv,backend=FULL._spd_inverse_enclosure(S,Racc)
                fixed += int(backend=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN")
                fallback += int(backend!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN")
                Kth=matrix_mul(PHt,Sinv); kn=RG._op2_upper(Kth)
                minS=min(minS,min(S[i][i].lo for i in range(3)))

                one_t=FULL.I(1.0)-kaw_t; one_z=FULL.I(1.0)-kz
                left_t=(one_t*rt).abs_upper(); left_z=(one_z*rz).abs_upper()
                post_from_residual=FULL.up(V5._norm2_upper(FULL.up(left_t+yRt),FULL.up(left_z+yRz))+ba)
                dxaw=V5._norm2_upper(awt.abs_upper(),az.abs_upper())
                post_from_triangle=FULL.up(aw_pred+dxaw)
                post_aw=min(post_from_residual,post_from_triangle)
                eaw1=FULL.up(alpha_hi*post_aw)
                maxEawPost=max(maxEawPost,post_aw); maxEaw=max(maxEaw,eaw1)

                rho=FULL.up(FULL.up(fn*chord)+FULL.up(eaw1+ba))
                corr=FULL.up(kn*rho)
                closed=backend=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN" and math.isfinite(corr) and corr<RANGE
                maxK=max(maxK,kn); maxD=max(maxD,corr); maxRho=max(maxRho,rho)
                row={
                    "p_cell":pi,"tangent_residual_cell":ti,"axial_residual_cell":zi,
                    "P_aw_variance":p.as_list(),
                    "first_tangent_residual_magnitude_mps2":rt.as_list(),
                    "first_axial_residual_mps2":rz.as_list(),
                    "first_attitude_correction_rad":d.as_list(),
                    "first_tangent_aw_correction_mps2":awt.as_list(),
                    "first_axial_aw_correction_mps2":az.as_list(),
                    "sample1_force_components_mps2":[x.as_list() for x in f],
                    "sample1_force_norm_upper_mps2":fn,
                    "attitude_rotation_chord_upper":chord,
                    "post_first_aw_error_from_residual_identity_upper_mps2":post_from_residual,
                    "post_first_aw_error_from_prior_plus_correction_upper_mps2":post_from_triangle,
                    "post_first_aw_error_upper_used_mps2":post_aw,
                    "post_prediction_aw_error_norm_upper_mps2":eaw1,
                    "sample1_residual_norm_upper_mps2":rho,
                    "inverse_backend":backend,"Ktheta_operator_norm_upper":kn,
                    "correction_norm_upper_rad":corr,"inside_9rad_range":closed,
                }
                rows.append(row)
                if not closed and bad is None: bad=row

    ok=bool(rows) and bad is None and not failures
    return {
        "schema":SCHEMA,
        "qualification":"OU3_P5_SAMPLE1_DIRECT_FIRST_RESIDUAL_COORDINATE_FULL_GAIN_V6",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "first_residual_components_are_primary_subdivision_coordinates":True,
        "first_gain_interval_never_inverted_to_reconstruct_residual":True,
        "first_residual_norm_ball_enforced_directly":True,
        "first_gravity_tangent_directional_bound_used":True,
        "first_gravity_axial_directional_bound_used":True,
        "first_latent_and_bias_source_component_caps_used":True,
        "first_attitude_and_aw_corrections_derived_forward_from_same_residual_cell":True,
        "first_posterior_aw_error_residual_identity_used":True,
        "sample1_nonaxial_force_included":True,
        "complete_3x3_accelerometer_innovation_used":True,
        "complete_3x3_attitude_gain_used":True,
        "full_propagated_aw_error_norm_used":True,
        "temporal_force_slew_assumed":False,
        "first_attitude_PSD_cross_axis_remainder_included":False,
        "sample1_S_covariance_update_included":False,
        "sample1_S_attitude_injection_included":False,
        "complete_sample1_branch_closed_here":False,"q8_word_promoted_here":False,
        "whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,
        "first_residual_norm_upper_mps2":rho0,
        "first_gravity_tangent_residual_upper_mps2":yRt,
        "first_gravity_axial_residual_upper_mps2":yRz,
        "first_gravity_total_residual_upper_mps2":yRnorm,
        "first_tangent_residual_source_cap_mps2":rt_source_cap,
        "first_axial_residual_source_cap_mps2":rz_source_cap,
        "candidate_joint_cells_before_residual_ball":len(pcells)*len(rtcells)*len(rzcells),
        "residual_incompatible_joint_cells_pruned":pruned,
        "evaluated_joint_cells":len(rows),"fixed_pivot_inverse_count":fixed,
        "spectral_fallback_inverse_count":fallback,"minimum_innovation_diagonal_lower":minS,
        "max_first_attitude_correction_norm_upper_rad":maxD0,
        "max_first_axial_aw_correction_abs_upper_mps2":maxAz,
        "max_sample1_force_norm_upper_mps2":maxF,"max_attitude_rotation_chord_upper":maxChord,
        "max_post_first_aw_error_upper_mps2":maxEawPost,
        "max_post_prediction_aw_error_norm_upper_mps2":maxEaw,
        "max_sample1_residual_norm_upper_mps2":maxRho,"max_Ktheta_operator_norm_upper":maxK,
        "max_correction_norm_upper_rad":maxD,"first_unclosed_joint_cell":bad,
        "P5_SAMPLE1_DIRECT_FIRST_RESIDUAL_FULL_GAIN_V6":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":(
            "ADD_FIRST_PSD_CROSS_AXIS_AND_SAMPLE1_S_BRANCH_THEN_SIGNED_CAYLEY_COMPOSE"
            if ok else "REMOVE_STRUCTURED_3X3_INTERVAL_INVERSE_FALLBACK_OR_REFINE_FIRST_FIXED_WITNESS"
        ),
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in (
        "source_generated_not_trajectory_fit","first_residual_components_are_primary_subdivision_coordinates",
        "first_gain_interval_never_inverted_to_reconstruct_residual","first_residual_norm_ball_enforced_directly",
        "first_gravity_tangent_directional_bound_used","first_gravity_axial_directional_bound_used",
        "first_latent_and_bias_source_component_caps_used",
        "first_attitude_and_aw_corrections_derived_forward_from_same_residual_cell",
        "first_posterior_aw_error_residual_identity_used","sample1_nonaxial_force_included",
        "complete_3x3_accelerometer_innovation_used","complete_3x3_attitude_gain_used",
        "full_propagated_aw_error_norm_used"):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in (
        "source_replay_used","filter_changed","temporal_force_slew_assumed",
        "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
        "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here",
        "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no direct-residual cells")
    if float(d.get("first_axial_residual_source_cap_mps2",math.inf))>float(d.get("first_residual_norm_upper_mps2",0.0))+1e-12:
        f.append("axial residual cap exceeds global residual norm")
    if not math.isfinite(float(d.get("max_correction_norm_upper_rad",math.nan))):f.append("nonfinite correction")
    st=d.get("P5_SAMPLE1_DIRECT_FIRST_RESIDUAL_FULL_GAIN_V6"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and w is not None:f.append("PASS retains witness")
    if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
    return f


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces",type=int,default=4); ap.add_argument("--source-cell-index",type=int,default=0)
    ap.add_argument("--p-pieces",type=int,default=24); ap.add_argument("--tangent-pieces",type=int,default=24)
    ap.add_argument("--axial-pieces",type=int,default=24); ap.add_argument("--output",type=Path,required=True)
    x=ap.parse_args(); d=build(x.domain,source_pieces=x.source_pieces,source_cell_index=x.source_cell_index,
        p_pieces=x.p_pieces,tangent_pieces=x.tangent_pieces,axial_pieces=x.axial_pieces)
    vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({
        "status":d["P5_SAMPLE1_DIRECT_FIRST_RESIDUAL_FULL_GAIN_V6"],
        "candidate_cells":d["candidate_joint_cells_before_residual_ball"],"pruned":d["residual_incompatible_joint_cells_pruned"],
        "cells":d["evaluated_joint_cells"],"fixed":d["fixed_pivot_inverse_count"],"fallback":d["spectral_fallback_inverse_count"],
        "rt_cap":d["first_tangent_residual_source_cap_mps2"],"rz_cap":d["first_axial_residual_source_cap_mps2"],
        "max_d0":d["max_first_attitude_correction_norm_upper_rad"],"max_az":d["max_first_axial_aw_correction_abs_upper_mps2"],
        "max_force":d["max_sample1_force_norm_upper_mps2"],"max_eaw_post":d["max_post_first_aw_error_upper_mps2"],
        "max_eaw1":d["max_post_prediction_aw_error_norm_upper_mps2"],"max_rho":d["max_sample1_residual_norm_upper_mps2"],
        "max_K":d["max_Ktheta_operator_norm_upper"],"max_d":d["max_correction_norm_upper_rad"],
        "first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf,
    },indent=2,sort_keys=True))
    return 0 if not vf else 2


if __name__=="__main__": raise SystemExit(main())
