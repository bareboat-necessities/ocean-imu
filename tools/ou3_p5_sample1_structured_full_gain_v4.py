#!/usr/bin/env python3
"""Full structured sample-1 accelerometer gain after the first canonical reset.

V3 showed that keeping the first state/residual correlation closes the most
reset-sensitive scalar row, but it still treated the sample-1 predicted force
as gravity-aligned and only bounded one component of the propagated latent
acceleration error.  Neither is a final P5 assumption.

For the ideal source-symmetric first accelerometer block, rotate about gravity
so the first tangent correction is d e1.  The first Joseph posterior of
(theta,a_w) then has the exact canonical structure

    P_theta = diag(a,a,Y),
    P_aw    = diag(b,b,bz),
    P_theta,aw = [[0,-c,0],[c,0,0],[0,0,0]],

where c<0 is the theta_y/a_wx posterior cross covariance.  The shipping
left-error reset followed by the corrected-body proof gauge acts on theta with

    L_theta = diag(1, R_x(d) [[1,-d/2],[d/2,1]])

and on world-vector a_w with R_x(d).  The first tangent a_w estimate correction
is exactly -beta d e2 in this ideal block, beta=p/(g t), while the axial update
is subdivided independently.  After one OU prediction the actual predicted
specific force therefore lies in the same y-z plane:

    f1 = R_x(d) [0, -alpha beta d, g + alpha a_z].

This producer forms the complete 3x3 accelerometer innovation covariance and
all three attitude-gain rows from that structured 6x6 (theta,a_w) marginal.  It
uses the fixed-pivot validated inverse when available and fails closed on the
spectral fallback.  The residual multiplier is not a raw 30 m/s^2 packet box:
we use the exact rotation-difference chord for the source-reachable attitude
radius and the full propagated physical a_w-error norm.  Thus the finite-angle
latent rotation is covered by ||R e_aw||=||e_aw||, with no component shortcut
and no temporal force-slew hypothesis.

The small off-diagonal attitude PSD remainder present before the first
accelerometer is still omitted from the canonical first posterior, and the
possible sample-1 S=0 covariance reduction is not applied here.  Those remain
explicit obligations.  This stage does not promote the complete sample-1
branch, q<8 word, or N_H_words.
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
import ou3_p5_sample1_reset_perp_scalar_channel as V1
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 4
RANGE = 9.0


def _diag3(a: Interval, b: Interval, c: Interval):
    z = FULL.I(0.0)
    return [[a,z,z],[z,b,z],[z,z,c]]


def _norm3_upper(v) -> float:
    s = 0.0
    for x in v:
        a = x.abs_upper()
        s = FULL.up(s + FULL.up(a*a))
    return FULL.up(math.sqrt(s))


def _Ltheta(d: Interval):
    R = SUB.rx(d)
    c = R[1][1]
    s = R[2][1]
    hd = FULL.I(0.5) * d
    l00 = c - s*hd
    l01 = -(c*hd) - s
    l10 = s + c*hd
    l11 = c - s*hd
    z = FULL.I(0.0); o = FULL.I(1.0)
    return [[o,z,z],[z,l00,l01],[z,l10,l11]], R


def _Htheta(f):
    x,y,z = f
    q = FULL.I(0.0)
    return [[q,z,-y],[-z,q,x],[y,-x,q]]


def _correction_chord_upper(d_hi: float) -> float:
    if not (0.0 <= d_hi < math.pi):
        raise ValueError("correction chord requires 0<=d<pi")
    return min(2.0, FULL.up(2.0*VT.sin_point(FULL.up(0.5*d_hi)).hi))


def _gravity_chord_from_cos(cos_lower: float) -> float:
    return min(2.0, FULL.up(math.sqrt(max(0.0, FULL.up(2.0*FULL.up(1.0-cos_lower))))))


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
        failures.append("structured full-gain witness requires first due source cell")
    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    rho0 = float(fr["combined_useful_residual_norm_upper_mps2"])
    dmax = float(fr["correction_norm_upper_rad"])
    aw_pred = float(fr["predicted_aw_error_norm_upper_mps2"])

    h = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    ba = float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    # The diagonal part of the first one-step remainder is covered here.  The
    # unknown cross-axis part remains an explicit later perturbation.
    t = Interval.outward_bounds(tilt, FULL.up(tilt+eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw+eps))
    Racc = FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))
    r = Racc[0][0]
    F,Q,_Rstep = FULL._transition_and_Q(src,dom)
    alpha = F[15][15]
    alpha_hi = alpha.hi
    qaw = Q[15][15]

    chord0 = _gravity_chord_from_cos(float(first["post_prediction_true_gravity_cosine_lower"]))
    pred_chord = _correction_chord_upper(float(first["first_prediction_transport_angle_upper_rad"]))

    pcells = SUB.parts(p_all.lo,p_all.hi,p_pieces)
    dcells = SUB.parts(0.0,dmax,d_pieces)
    rows=[]; bad=None; fixed=fallback=0
    maxK=maxD=maxRho=maxF=maxEaw=maxChord=0.0
    minS=math.inf

    for pi,p in enumerate(pcells):
        D = FULL.I(g*g)*t + p + r
        a = t*(p+r)/D
        c0 = -(FULL.I(g)*t*p/D)
        b = p*(FULL.I(g*g)*t+r)/D
        bz = p*r/(p+r)
        Pth0 = _diag3(a,a,Y)
        Paw0 = _diag3(b,b,bz)
        C0 = FULL._zero(3,3)
        C0[0][1] = -c0
        C0[1][0] = c0
        beta = p/(FULL.I(g)*t)
        kgz = p/(p+r)
        azmax = FULL.up(kgz.abs_upper()*rho0)
        azcells = SUB.parts(-azmax,azmax,axial_pieces)

        for di,d in enumerate(dcells):
            L,Rx = _Ltheta(d)
            Pth = matrix_mul(matrix_mul(L,Pth0),matrix_transpose(L))
            Paw_r = matrix_mul(matrix_mul(Rx,Paw0),matrix_transpose(Rx))
            C_r = matrix_mul(matrix_mul(L,C0),matrix_transpose(Rx))
            Paw = [[alpha.square()*Paw_r[i][j] + (qaw if i==j else FULL.I(0.0)) for j in range(3)] for i in range(3)]
            C = [[alpha*C_r[i][j] for j in range(3)] for i in range(3)]

            d_hi=max(0.0,d.hi)
            chord=min(2.0,FULL.up(chord0+FULL.up(_correction_chord_upper(d_hi)+pred_chord)))
            maxChord=max(maxChord,chord)
            beta_d_abs=FULL.up(beta.abs_upper()*d_hi)

            for ai,az in enumerate(azcells):
                raw=[FULL.I(0.0), -(alpha*beta*d), FULL.I(g)+alpha*az]
                f=FULL._mat_vec(Rx,raw)
                fn=_norm3_upper(f)
                maxF=max(maxF,fn)
                Ht=_Htheta(f)
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

                azabs=az.abs_upper()
                dxaw=FULL.up(math.sqrt(FULL.up(beta_d_abs*beta_d_abs + FULL.up(azabs*azabs))))
                eaw1=FULL.up(alpha_hi*FULL.up(aw_pred+dxaw))
                rho=FULL.up(FULL.up(fn*chord)+FULL.up(eaw1+ba))
                corr=FULL.up(kn*rho)
                closed=backend=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN" and math.isfinite(corr) and corr<RANGE
                maxK=max(maxK,kn); maxD=max(maxD,corr); maxRho=max(maxRho,rho); maxEaw=max(maxEaw,eaw1)
                row={
                    "p_cell":pi,"d_cell":di,"axial_cell":ai,
                    "P_aw_variance":p.as_list(),"first_correction_rad":d.as_list(),
                    "first_axial_aw_correction_mps2":az.as_list(),
                    "sample1_force_components_mps2":[x.as_list() for x in f],
                    "sample1_force_norm_upper_mps2":fn,
                    "attitude_rotation_chord_upper":chord,
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
        "qualification":"OU3_P5_SAMPLE1_CANONICAL_RESET_FULL_FORCE_FULL_LATENT_GAIN_V4",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "canonical_first_full_theta_aw_Joseph_marginal_used":True,
        "shipping_left_error_reset_used_in_theta_marginal":True,
        "corrected_body_Rx_gauge_used_in_aw_marginal":True,
        "source_correlated_tangent_aw_mean_relation_used":True,
        "sample1_nonaxial_force_included":True,
        "complete_3x3_accelerometer_innovation_used":True,
        "complete_3x3_attitude_gain_used":True,
        "full_propagated_aw_error_norm_used":True,
        "latent_finite_rotation_covered_by_orthogonal_norm_invariance":True,
        "temporal_force_slew_assumed":False,
        "sample1_body_rotation_removed_by_simultaneous_orthogonal_gauge":True,
        "first_attitude_PSD_cross_axis_remainder_included":False,
        "sample1_S_covariance_update_included":False,
        "sample1_S_attitude_injection_included":False,
        "complete_sample1_branch_closed_here":False,"q8_word_promoted_here":False,
        "whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,
        "evaluated_joint_cells":len(rows),"fixed_pivot_inverse_count":fixed,
        "spectral_fallback_inverse_count":fallback,"minimum_innovation_diagonal_lower":minS,
        "max_sample1_force_norm_upper_mps2":maxF,"max_attitude_rotation_chord_upper":maxChord,
        "max_post_prediction_aw_error_norm_upper_mps2":maxEaw,
        "max_sample1_residual_norm_upper_mps2":maxRho,"max_Ktheta_operator_norm_upper":maxK,
        "max_correction_norm_upper_rad":maxD,"first_unclosed_joint_cell":bad,
        "P5_SAMPLE1_STRUCTURED_FULL_GAIN_V4":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":(
            "ADD_FIRST_PSD_CROSS_AXIS_AND_SAMPLE1_S_COVARIANCE_BRANCH_THEN_SIGNED_CAYLEY_COMPOSE"
            if ok else "REFINE_STRUCTURED_3X3_GAIN_OR_SOURCE_STATE_SUBDIVISION_AT_FIRST_WITNESS"
        ),
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in ("source_generated_not_trajectory_fit","canonical_first_full_theta_aw_Joseph_marginal_used",
              "shipping_left_error_reset_used_in_theta_marginal","corrected_body_Rx_gauge_used_in_aw_marginal",
              "source_correlated_tangent_aw_mean_relation_used","sample1_nonaxial_force_included",
              "complete_3x3_accelerometer_innovation_used","complete_3x3_attitude_gain_used",
              "full_propagated_aw_error_norm_used","latent_finite_rotation_covered_by_orthogonal_norm_invariance",
              "sample1_body_rotation_removed_by_simultaneous_orthogonal_gauge"):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in ("source_replay_used","filter_changed","temporal_force_slew_assumed",
              "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
              "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here","q8_word_promoted_here",
              "whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no cells")
    if not math.isfinite(float(d.get("max_correction_norm_upper_rad",math.nan))):f.append("nonfinite correction")
    st=d.get("P5_SAMPLE1_STRUCTURED_FULL_GAIN_V4"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and w is not None:f.append("PASS retains witness")
    if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
    return f


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces",type=int,default=4); ap.add_argument("--source-cell-index",type=int,default=0)
    ap.add_argument("--p-pieces",type=int,default=24); ap.add_argument("--d-pieces",type=int,default=24)
    ap.add_argument("--axial-pieces",type=int,default=24); ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args(); d=build(a.domain,source_pieces=a.source_pieces,source_cell_index=a.source_cell_index,
        p_pieces=a.p_pieces,d_pieces=a.d_pieces,axial_pieces=a.axial_pieces)
    vf=validate(d); d["validation_failures"]=vf; a.output.write_text(json.dumps(d,indent=2,sort_keys=True))
    print(json.dumps({"status":d["P5_SAMPLE1_STRUCTURED_FULL_GAIN_V4"],"cells":d["evaluated_joint_cells"],
      "fixed":d["fixed_pivot_inverse_count"],"fallback":d["spectral_fallback_inverse_count"],
      "max_force":d["max_sample1_force_norm_upper_mps2"],"max_chord":d["max_attitude_rotation_chord_upper"],
      "max_eaw":d["max_post_prediction_aw_error_norm_upper_mps2"],"max_rho":d["max_sample1_residual_norm_upper_mps2"],
      "max_K":d["max_Ktheta_operator_norm_upper"],"max_d":d["max_correction_norm_upper_rad"],
      "first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True))
    return 0 if not vf else 2

if __name__=="__main__": raise SystemExit(main())
