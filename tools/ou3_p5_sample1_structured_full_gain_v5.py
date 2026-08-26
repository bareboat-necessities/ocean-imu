#!/usr/bin/env python3
"""Residual-coupled V5 full sample-1 accelerometer core for OU-III P5.

V4 retained the full source-correlated sample-1 force and complete 3x3 attitude
gain, but it bounded the physical latent error after the first accelerometer by

    ||e_aw^+|| <= ||e_aw^-|| + ||Delta a_w_hat||.

That inequality destroys the defining measurement dependency.  In the
canonical first gravity gauge the useful residual is

    r0 = y_R + u_aw + b_a,

where u_aw is the physical world-acceleration error expressed in the same body
proof gauge.  The ideal first Kalman aw gain is diagonal,

    D_aw = diag(k_t,k_t,k_z),
    k_t = p/(g^2 t+p+r),  k_z=p/(p+r),

and therefore exactly

    e_aw^+ = u_aw-D_aw r0 = (I-D_aw)r0-y_R-b_a.          (1)

The first attitude correction magnitude and axial aw correction are generated
by that same residual:

    d = k_theta ||r_t||,  k_theta=g t/(g^2 t+p+r),
    a_z = k_z r_z.

This producer subdivides (p,d,a_z), maps each child back to its compatible
first-residual coordinates, and removes children that cannot intersect the
already certified ||r0||<=rho0 ball.  The retained children are tightened by
that same ball before reset/force propagation.  Equation (1) then supplies a
posterior physical aw-error bound; the older triangle bound is retained only as
an independent second upper bound and the minimum of the two is used.

The complete structured sample-1 force, 3x3 innovation, and 3x3 attitude gain
from V4 are retained.  The finite-angle latent term is covered through the full
propagated aw-error norm, not one component.  As in V4, the small pre-first
attitude PSD cross-axis remainder and the sample-1 S covariance/update branch
remain explicit later obligations.  No q<8 word or P5 promotion is made here.
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
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 5
RANGE = 9.0


def _min_abs(x: Interval) -> float:
    if x.lo <= 0.0 <= x.hi:
        return 0.0
    return min(abs(x.lo), abs(x.hi))


def _norm2_upper(a: float, b: float) -> float:
    return FULL.up(math.sqrt(FULL.up(FULL.up(a*a) + FULL.up(b*b))))


def _norm2_lower(a: float, b: float) -> float:
    s = FULL.down(FULL.down(a*a) + FULL.down(b*b))
    return 0.0 if s <= 0.0 else FULL.down(math.sqrt(s))


def _compatible_residual_child(d: Interval, az: Interval, ktheta: Interval,
                               kz: Interval, rho0: float):
    """Intersect one (d,az) child with d=k_theta|r_t|, az=k_z r_z, ||r||<=rho0."""
    if not (ktheta.lo > 0.0 and kz.lo > 0.0 and rho0 > 0.0):
        raise RuntimeError("positive first residual gains/radius required")
    dlo = max(0.0, d.lo)
    dhi = max(dlo, d.hi)
    d0 = Interval(dlo, dhi)
    rt0 = d0 / ktheta
    rz0 = az / kz
    rt_min = max(0.0, rt0.lo)
    rz_min = _min_abs(rz0)
    if _norm2_lower(rt_min, rz_min) > rho0:
        return None

    # Every point in this axial cell has |rz|>=rz_min, hence any source point
    # in the residual ball must satisfy |rt|<=sqrt(rho0^2-rz_min^2).
    rt_cap_sq = max(0.0, FULL.up(FULL.up(rho0*rho0) - FULL.down(rz_min*rz_min)))
    rt_cap = FULL.up(math.sqrt(rt_cap_sq))
    d_cap = FULL.up(ktheta.hi * rt_cap)
    de_hi = min(dhi, d_cap)
    if de_hi < dlo:
        return None
    de = Interval(dlo, de_hi)
    rt = de / ktheta

    # Symmetric tightening for the axial correction from the minimum tangent
    # residual represented by the retained d child.
    rt_min = max(0.0, rt.lo)
    rz_cap_sq = max(0.0, FULL.up(FULL.up(rho0*rho0) - FULL.down(rt_min*rt_min)))
    rz_cap = FULL.up(math.sqrt(rz_cap_sq))
    az_cap = FULL.up(kz.hi * rz_cap)
    ae_lo = max(az.lo, -az_cap)
    ae_hi = min(az.hi, az_cap)
    if ae_hi < ae_lo:
        return None
    ae = Interval(ae_lo, ae_hi)
    rz = ae / kz
    rt = de / ktheta
    if _norm2_lower(max(0.0, rt.lo), _min_abs(rz)) > rho0:
        return None
    return de, ae, rt, rz


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          d_pieces: int = 24, axial_pieces: int = 24) -> dict:
    FULL3._install_backend()
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    first = FIRST.build(path, source_pieces=source_pieces)
    vec = VECTOR.build()
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vec)]

    src, phase = RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        failures.append("V5 residual-coupled witness requires first due source cell")
    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    rho0 = float(fr["combined_useful_residual_norm_upper_mps2"])
    yR0 = float(fr["rotational_residual_norm_upper_mps2"])
    dmax = float(fr["correction_norm_upper_rad"])
    aw_pred = float(fr["predicted_aw_error_norm_upper_mps2"])

    h = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    ba = float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    t = Interval.outward_bounds(tilt, FULL.up(tilt+eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw+eps))
    Racc = FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))
    r = Racc[0][0]
    F,Q,_ = FULL._transition_and_Q(src,dom)
    alpha = F[15][15]
    alpha_hi = alpha.hi
    qaw = Q[15][15]

    chord0 = V4._gravity_chord_from_cos(float(first["post_prediction_true_gravity_cosine_lower"]))
    pred_chord = V4._correction_chord_upper(float(first["first_prediction_transport_angle_upper_rad"]))

    pcells = SUB.parts(p_all.lo,p_all.hi,p_pieces)
    dcells = SUB.parts(0.0,dmax,d_pieces)
    rows=[]; bad=None; fixed=fallback=0; pruned=0
    maxK=maxD=maxRho=maxF=maxEaw=maxEawPost=maxChord=maxResidualChild=0.0
    minS=math.inf

    for pi,p in enumerate(pcells):
        D = FULL.I(g*g)*t + p + r
        a = t*(p+r)/D
        c0 = -(FULL.I(g)*t*p/D)
        b = p*(FULL.I(g*g)*t+r)/D
        bz = p*r/(p+r)
        Pth0 = V4._diag3(a,a,Y)
        Paw0 = V4._diag3(b,b,bz)
        C0 = FULL._zero(3,3)
        C0[0][1] = -c0
        C0[1][0] = c0

        ktheta = FULL.I(g)*t/D
        kaw_t = p/D
        kz = p/(p+r)
        if not (ktheta.lo > 0.0 and kaw_t.lo >= 0.0 and kz.lo > 0.0):
            failures.append(f"first canonical gain lost positivity p={pi}")
            continue
        beta = p/(FULL.I(g)*t)
        azmax = FULL.up(kz.abs_upper()*rho0)
        azcells = SUB.parts(-azmax,azmax,axial_pieces)

        for di,d0 in enumerate(dcells):
            for ai,az0 in enumerate(azcells):
                child = _compatible_residual_child(d0,az0,ktheta,kz,rho0)
                if child is None:
                    pruned += 1
                    continue
                d,az,rt,rz = child
                residual_child_norm = _norm2_upper(rt.abs_upper(),rz.abs_upper())
                maxResidualChild=max(maxResidualChild,min(residual_child_norm,rho0))

                L,Rx = V4._Ltheta(d)
                Pth = matrix_mul(matrix_mul(L,Pth0),matrix_transpose(L))
                Paw_r = matrix_mul(matrix_mul(Rx,Paw0),matrix_transpose(Rx))
                C_r = matrix_mul(matrix_mul(L,C0),matrix_transpose(Rx))
                Paw = [[alpha.square()*Paw_r[i][j] + (qaw if i==j else FULL.I(0.0)) for j in range(3)] for i in range(3)]
                C = [[alpha*C_r[i][j] for j in range(3)] for i in range(3)]

                d_hi=max(0.0,d.hi)
                chord=min(2.0,FULL.up(chord0+FULL.up(V4._correction_chord_upper(d_hi)+pred_chord)))
                maxChord=max(maxChord,chord)

                # Source-correlated filter mean after the first update.
                raw=[FULL.I(0.0), -(alpha*beta*d), FULL.I(g)+alpha*az]
                f=FULL._mat_vec(Rx,raw)
                fn=V4._norm3_upper(f)
                maxF=max(maxF,fn)

                # Complete structured sample-1 covariance/gain from V4.
                Ht=V4._Htheta(f)
                PHt=matrix_add(matrix_mul(Pth,matrix_transpose(Ht)),C)
                HP=matrix_mul(Ht,Pth)
                HC=matrix_mul(Ht,C)
                S=matrix_add(matrix_mul(HP,matrix_transpose(Ht)),HC)
                S=matrix_add(S,matrix_transpose(HC))
                S=matrix_add(S,Paw)
                S=FULL.matrix_symmetric_hull(matrix_add(S,Racc))
                Sinv,backend=FULL._spd_inverse_enclosure(S,Racc)
                fixed += int(backend=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN")
                fallback += int(backend!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN")
                Kth=matrix_mul(PHt,Sinv)
                kn=RG._op2_upper(Kth)
                minS=min(minS,min(S[i][i].lo for i in range(3)))

                # Exact first residual-coordinate identity (ideal canonical core):
                # e_aw^+=(I-D_aw)r0-y_R-b_a.  The residual cell supplies tighter
                # component ranges than the global rho0 ball.  Keep the older
                # prior+correction triangle as an independent second upper bound.
                one_t = FULL.I(1.0)-kaw_t
                one_z = FULL.I(1.0)-kz
                left_t = (one_t*rt).abs_upper()
                left_z = (one_z*rz).abs_upper()
                left_box = _norm2_upper(left_t,left_z)
                left_ball = FULL.up(max(one_t.abs_upper(),one_z.abs_upper())*rho0)
                left = min(left_box,left_ball)
                post_from_residual = FULL.up(left + FULL.up(yR0+ba))

                beta_d = (beta*d).abs_upper()
                azabs = az.abs_upper()
                dxaw = _norm2_upper(beta_d,azabs)
                post_from_triangle = FULL.up(aw_pred+dxaw)
                post_aw = min(post_from_residual,post_from_triangle)
                eaw1 = FULL.up(alpha_hi*post_aw)
                maxEawPost=max(maxEawPost,post_aw)
                maxEaw=max(maxEaw,eaw1)

                rho=FULL.up(FULL.up(fn*chord)+FULL.up(eaw1+ba))
                corr=FULL.up(kn*rho)
                closed=backend=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN" and math.isfinite(corr) and corr<RANGE
                maxK=max(maxK,kn); maxD=max(maxD,corr); maxRho=max(maxRho,rho)
                row={
                    "p_cell":pi,"d_cell":di,"axial_cell":ai,
                    "P_aw_variance":p.as_list(),
                    "first_correction_rad_after_residual_ball":d.as_list(),
                    "first_axial_aw_correction_mps2_after_residual_ball":az.as_list(),
                    "first_tangent_residual_magnitude_mps2":rt.as_list(),
                    "first_axial_residual_mps2":rz.as_list(),
                    "first_residual_norm_global_upper_mps2":rho0,
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
                if not closed and bad is None:
                    bad=row

    ok=bool(rows) and bad is None and not failures
    return {
        "schema":SCHEMA,
        "qualification":"OU3_P5_SAMPLE1_FIRST_RESIDUAL_COUPLED_FULL_GAIN_V5",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "canonical_first_full_theta_aw_Joseph_marginal_used":True,
        "first_residual_to_tangent_correction_identity_used":True,
        "first_residual_to_axial_aw_correction_identity_used":True,
        "first_residual_norm_ball_enforced_on_joint_d_axial_cells":True,
        "first_posterior_aw_error_residual_identity_used":True,
        "posterior_aw_bound_uses_min_of_residual_identity_and_triangle":True,
        "independent_prior_plus_correction_aw_bound_used_as_only_bound":False,
        "sample1_nonaxial_force_included":True,
        "complete_3x3_accelerometer_innovation_used":True,
        "complete_3x3_attitude_gain_used":True,
        "full_propagated_aw_error_norm_used":True,
        "latent_finite_rotation_covered_by_orthogonal_norm_invariance":True,
        "temporal_force_slew_assumed":False,
        "first_attitude_PSD_cross_axis_remainder_included":False,
        "sample1_S_covariance_update_included":False,
        "sample1_S_attitude_injection_included":False,
        "complete_sample1_branch_closed_here":False,"q8_word_promoted_here":False,
        "whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,
        "candidate_joint_cells_before_residual_ball":len(pcells)*len(dcells)*axial_pieces,
        "residual_incompatible_joint_cells_pruned":pruned,
        "evaluated_joint_cells":len(rows),"fixed_pivot_inverse_count":fixed,
        "spectral_fallback_inverse_count":fallback,"minimum_innovation_diagonal_lower":minS,
        "max_retained_first_residual_norm_box_upper_capped_for_report_mps2":maxResidualChild,
        "max_sample1_force_norm_upper_mps2":maxF,"max_attitude_rotation_chord_upper":maxChord,
        "max_post_first_aw_error_upper_mps2":maxEawPost,
        "max_post_prediction_aw_error_norm_upper_mps2":maxEaw,
        "max_sample1_residual_norm_upper_mps2":maxRho,"max_Ktheta_operator_norm_upper":maxK,
        "max_correction_norm_upper_rad":maxD,"first_unclosed_joint_cell":bad,
        "P5_SAMPLE1_FIRST_RESIDUAL_COUPLED_FULL_GAIN_V5":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":(
            "ADD_FIRST_PSD_CROSS_AXIS_AND_SAMPLE1_S_BRANCH_THEN_SIGNED_CAYLEY_COMPOSE"
            if ok else "REFINE_FULL_3X3_GAIN_AT_FIRST_RESIDUAL_COUPLED_WITNESS_OR_REMOVE_INTERVAL_INVERSE_FALLBACK"
        ),
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in (
        "source_generated_not_trajectory_fit","canonical_first_full_theta_aw_Joseph_marginal_used",
        "first_residual_to_tangent_correction_identity_used","first_residual_to_axial_aw_correction_identity_used",
        "first_residual_norm_ball_enforced_on_joint_d_axial_cells","first_posterior_aw_error_residual_identity_used",
        "posterior_aw_bound_uses_min_of_residual_identity_and_triangle","sample1_nonaxial_force_included",
        "complete_3x3_accelerometer_innovation_used","complete_3x3_attitude_gain_used",
        "full_propagated_aw_error_norm_used","latent_finite_rotation_covered_by_orthogonal_norm_invariance"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used","filter_changed","independent_prior_plus_correction_aw_bound_used_as_only_bound",
        "temporal_force_slew_assumed","first_attitude_PSD_cross_axis_remainder_included",
        "sample1_S_covariance_update_included","sample1_S_attitude_injection_included",
        "complete_sample1_branch_closed_here","q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:
        f.append("no residual-compatible cells")
    if int(d.get("candidate_joint_cells_before_residual_ball",0)) < int(d.get("evaluated_joint_cells",0)):
        f.append("residual-ball pruning increased cell count")
    if not math.isfinite(float(d.get("max_correction_norm_upper_rad",math.nan))):
        f.append("nonfinite correction")
    st=d.get("P5_SAMPLE1_FIRST_RESIDUAL_COUPLED_FULL_GAIN_V5")
    w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and w is not None:
        f.append("PASS retains witness")
    if st=="NOT_ESTABLISHED" and w is None:
        f.append("missing witness")
    return f


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces",type=int,default=4)
    ap.add_argument("--source-cell-index",type=int,default=0)
    ap.add_argument("--p-pieces",type=int,default=24)
    ap.add_argument("--d-pieces",type=int,default=24)
    ap.add_argument("--axial-pieces",type=int,default=24)
    ap.add_argument("--output",type=Path,required=True)
    x=ap.parse_args()
    d=build(x.domain,source_pieces=x.source_pieces,source_cell_index=x.source_cell_index,
            p_pieces=x.p_pieces,d_pieces=x.d_pieces,axial_pieces=x.axial_pieces)
    vf=validate(d); d["validation_failures"]=vf
    x.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({
        "status":d["P5_SAMPLE1_FIRST_RESIDUAL_COUPLED_FULL_GAIN_V5"],
        "candidate_cells":d["candidate_joint_cells_before_residual_ball"],
        "pruned":d["residual_incompatible_joint_cells_pruned"],
        "cells":d["evaluated_joint_cells"],"fixed":d["fixed_pivot_inverse_count"],
        "fallback":d["spectral_fallback_inverse_count"],
        "max_force":d["max_sample1_force_norm_upper_mps2"],
        "max_eaw_post":d["max_post_first_aw_error_upper_mps2"],
        "max_eaw1":d["max_post_prediction_aw_error_norm_upper_mps2"],
        "max_rho":d["max_sample1_residual_norm_upper_mps2"],
        "max_K":d["max_Ktheta_operator_norm_upper"],"max_d":d["max_correction_norm_upper_rad"],
        "first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],
        "validation_failures":vf,
    },indent=2,sort_keys=True))
    return 0 if not vf else 2


if __name__=="__main__":
    raise SystemExit(main())
