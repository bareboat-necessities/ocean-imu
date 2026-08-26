#!/usr/bin/env python3
"""Analytic 1+2 block sample-1 gain for the direct-residual P5 core.

V6 made the first accelerometer residual components the primary source
coordinates.  Its first unclosed cell was no longer a physical fixed-pivot
witness; it was the generic 3x3 SPD spectral inverse fallback.  That fallback is
unnecessary for the canonical source structure.

Work in the proof gauge obtained by applying R_x(d)^T simultaneously to the
post-reset attitude, world-vector a_w, measurement residual, and sample-1 force.
This orthogonal change does not alter ||K_theta||_2.  The shipping left-error
reset reduces to

    G = [[1,0,0],[0,1,-h],[0,h,1]],  h=d/2,

and the sample-1 force is exactly f=[0,f_y,f_z].  The structured covariance then
has

    P_theta = diag(a,.) with yz block G_yz diag(a,Y) G_yz^T,
    P_aw    = diag(B_t,B_t,B_z),
    C_theta,aw = [[0,-alpha c,0],[alpha c,0,0],
                  [alpha h c,0,0]],

where c<0 is the first tangent theta/aw cross covariance.  Therefore the 3x3
innovation is exactly block diagonal: one scalar x-measurement block and one
2x2 (y,z)-measurement block.

Scalar x block.
Let u=[1,h], v=[-h,1], H_x=[f_z,-f_y],

    U=f_z-h f_y,  V=-h f_z-f_y,
    N_u=a U+alpha c,
    Delta=alpha^2 (a b-c^2)+a Q_aw >0.

Then exactly

    S_x = N_u^2/a + Delta/a + Y V^2 + r,
    ||Cov(theta_yz,y_x)||^2
        = (1+h^2) [N_u^2+(YV)^2].

No cancellation-prone subtraction is used in S_x.

2x2 y/z block.
Set c_x=-alpha c, q=c_x-a f_z and A=Delta+a r.  The innovation determinant has
the positive identity

    det S_yz = f_y^2 A + (q^2+A)(B_z+r)/a >0.

The theta_x gain row is exactly

    K_xy = q(B_z+r)/det S_yz,
    K_xz = f_y A/det S_yz.

Hence the complete attitude gain has two orthogonal singular blocks and

    ||K_theta||_2 = max(sqrt(K_xy^2+K_xz^2),
                        ||K_theta_yz,x||_2).

All arithmetic below is outward interval arithmetic.  No 3x3 inverse, fixed
pivot, or spectral inverse enclosure is used.  The residual/state/force bounds
are the same direct first-residual coordinates as V6.  The pre-first attitude
PSD cross-axis remainder and sample-1 S covariance/update are still explicit
later obligations; this producer does not promote sample 1, q<8, P5, or
N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_p5_sample1_structured_full_gain_v4 as V4
import ou3_p5_sample1_structured_full_gain_v5 as V5
import ou3_p5_sample1_structured_full_gain_v6 as V6
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN
SCHEMA=7
RANGE=9.0


def _sqrt_hi(x: Interval) -> float:
    if x.hi < 0.0:
        raise RuntimeError("negative interval in norm square")
    return FULL.up(math.sqrt(max(0.0,x.hi)))


def _block_gain_bounds(a: Interval, Y: Interval, c0: Interval,
                       alpha: Interval, qaw: Interval, b: Interval, bz: Interval,
                       det_first: Interval, d: Interval, fy: Interval, fz: Interval,
                       r: Interval):
    """Return exact-structure upper bounds for the two singular gain blocks."""
    if not (a.lo>0.0 and Y.lo>0.0 and r.lo>0.0 and det_first.lo>0.0):
        raise RuntimeError("positive canonical covariance/noise floors required")
    h=FULL.I(0.5)*d
    one=FULL.I(1.0)
    h2=h.square()

    Bt=alpha.square()*b+qaw
    Bz=alpha.square()*bz+qaw
    delta=alpha.square()*det_first+a*qaw
    if not (Bt.lo>0.0 and Bz.lo>0.0 and delta.lo>0.0):
        raise RuntimeError("positive predicted aw/determinant floor required")

    # Scalar x-measurement block.  The completed-square representation is
    # positive term by term and keeps the first posterior determinant explicit.
    U=fz-h*fy
    V=-(h*fz)-fy
    cu=alpha*c0
    Nu=a*U+cu
    Sx=Nu.square()/a + delta/a + Y*V.square() + r
    if not Sx.lo>0.0:
        raise RuntimeError("scalar innovation lost positive floor")
    nperp2=(one+h2)*(Nu.square()+(Y*V).square())
    k_perp=FULL.up(_sqrt_hi(nperp2)/Sx.lo)

    # y/z measurement block.  Use its positive determinant identity rather than
    # interval-Gauss-Jordan subtraction.
    cx=-(alpha*c0)
    q=cx-a*fz
    A=delta+a*r
    if not A.lo>0.0:
        raise RuntimeError("2x2 positive A floor lost")
    det=fy.square()*A + (q.square()+A)*(Bz+r)/a
    if not det.lo>0.0:
        raise RuntimeError("2x2 innovation determinant lost positivity")
    ky_num=q*(Bz+r)
    kz_num=fy*A
    kx=FULL.up(_sqrt_hi(ky_num.square()+kz_num.square())/det.lo)
    return max(kx,k_perp), {
        "half_correction":h.as_list(),
        "Bt":Bt.as_list(),"Bz":Bz.as_list(),"Delta":delta.as_list(),
        "scalar_U":U.as_list(),"scalar_V":V.as_list(),"scalar_Nu":Nu.as_list(),
        "scalar_Sx_positive":Sx.as_list(),"scalar_gain_yz_norm_upper":k_perp,
        "two_by_two_q":q.as_list(),"two_by_two_A":A.as_list(),
        "two_by_two_det_positive":det.as_list(),"two_by_two_theta_x_gain_norm_upper":kx,
    }


def build(domain_path: Path=DEFAULT_DOMAIN, *, source_pieces:int=4,
          source_cell_index:int=0, p_pieces:int=24,
          tangent_pieces:int=24, axial_pieces:int=24)->dict:
    FULL3._install_backend()
    path=Path(domain_path).resolve(); dom=json.loads(path.read_text(encoding="utf-8"))
    first=FIRST.build(path,source_pieces=source_pieces); vec=VECTOR.build()
    failures=[f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vec)]
    src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
    if phase!="due": failures.append("V7 analytic-block witness requires first due source cell")
    fr=first["source_cells"][source_cell_index]
    p_all=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"]))
    rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
    aw_pred=float(fr["predicted_aw_error_norm_upper_mps2"])

    hstep=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"])
    ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    tilt,yaw,eps=RG._attitude_covariance_epsilon(path,hstep)
    t=Interval.outward_bounds(tilt,FULL.up(tilt+eps)); Y=Interval.outward_bounds(yaw,FULL.up(yaw+eps))
    Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); r=Racc[0][0]
    F,Q,_=FULL._transition_and_Q(src,dom); alpha=F[15][15]; alpha_hi=alpha.hi; qaw=Q[15][15]

    cos0=float(first["post_prediction_true_gravity_cosine_lower"])
    sin_hi=1.0 if cos0<0.0 else FULL.up(math.sqrt(max(0.0,FULL.up(1.0-FULL.down(cos0*cos0)))))
    yRt=FULL.up(g*sin_hi); yRz=FULL.up(g*max(0.0,FULL.up(1.0-cos0)))
    rt_cap=min(rho0,FULL.up(aw_pred+FULL.up(yRt+ba)))
    rz_cap=min(rho0,FULL.up(aw_pred+FULL.up(yRz+ba)))
    chord0=V4._gravity_chord_from_cos(cos0)
    pred_chord=V4._correction_chord_upper(float(first["first_prediction_transport_angle_upper_rad"]))

    pcells=SUB.parts(p_all.lo,p_all.hi,p_pieces)
    rtcells=SUB.parts(0.0,rt_cap,tangent_pieces)
    rzcells=SUB.parts(-rz_cap,rz_cap,axial_pieces)
    rows=[]; bad=None; pruned=0
    maxK=maxD=maxRho=maxF=maxEaw=maxEawPost=maxD0=maxAz=0.0
    minSx=minDet=math.inf

    for pi,p in enumerate(pcells):
        D=FULL.I(g*g)*t+p+r
        a=t*(p+r)/D; c0=-(FULL.I(g)*t*p/D)
        b=p*(FULL.I(g*g)*t+r)/D; bz=p*r/(p+r)
        det_first=t*p*r/D
        ktheta=FULL.I(g)*t/D; kaw_t=p/D; kz=p/(p+r)
        if not (a.lo>0.0 and det_first.lo>0.0 and ktheta.lo>0.0 and kz.lo>0.0):
            failures.append(f"canonical first block lost positivity p={pi}"); continue

        for ti,rt0 in enumerate(rtcells):
            for zi,rz0 in enumerate(rzcells):
                child=V6._residual_child(rt0,rz0,rho0)
                if child is None:
                    pruned += 1; continue
                rt,rz=child
                d=ktheta*rt; awt=kaw_t*rt; az=kz*rz
                maxD0=max(maxD0,d.abs_upper()); maxAz=max(maxAz,az.abs_upper())

                # Work directly in the un-Rx orthogonal gauge.  This is the
                # source-correlated sample-1 force before the simultaneous Rx
                # congruence used by V4/V6.
                fy=-(alpha*awt); fz=FULL.I(g)+alpha*az
                fn=V5._norm2_upper(fy.abs_upper(),fz.abs_upper()); maxF=max(maxF,fn)
                kn,detail=_block_gain_bounds(a,Y,c0,alpha,qaw,b,bz,det_first,d,fy,fz,r)
                minSx=min(minSx,float(detail["scalar_Sx_positive"][0]))
                minDet=min(minDet,float(detail["two_by_two_det_positive"][0]))

                d_hi=max(0.0,d.hi)
                chord=min(2.0,FULL.up(chord0+FULL.up(V4._correction_chord_upper(d_hi)+pred_chord)))
                one_t=FULL.I(1.0)-kaw_t; one_z=FULL.I(1.0)-kz
                left_t=(one_t*rt).abs_upper(); left_z=(one_z*rz).abs_upper()
                post_from_residual=FULL.up(V5._norm2_upper(FULL.up(left_t+yRt),FULL.up(left_z+yRz))+ba)
                dxaw=V5._norm2_upper(awt.abs_upper(),az.abs_upper())
                post_from_triangle=FULL.up(aw_pred+dxaw)
                post_aw=min(post_from_residual,post_from_triangle); eaw1=FULL.up(alpha_hi*post_aw)
                maxEawPost=max(maxEawPost,post_aw); maxEaw=max(maxEaw,eaw1)
                rho=FULL.up(FULL.up(fn*chord)+FULL.up(eaw1+ba))
                corr=FULL.up(kn*rho)
                closed=math.isfinite(corr) and corr<RANGE
                maxK=max(maxK,kn); maxD=max(maxD,corr); maxRho=max(maxRho,rho)
                row={
                    "p_cell":pi,"tangent_residual_cell":ti,"axial_residual_cell":zi,
                    "P_aw_variance":p.as_list(),"first_tangent_residual_magnitude_mps2":rt.as_list(),
                    "first_axial_residual_mps2":rz.as_list(),"first_attitude_correction_rad":d.as_list(),
                    "first_tangent_aw_correction_mps2":awt.as_list(),"first_axial_aw_correction_mps2":az.as_list(),
                    "sample1_force_unRx_components_yz_mps2":[fy.as_list(),fz.as_list()],
                    "sample1_force_norm_upper_mps2":fn,"post_prediction_aw_error_norm_upper_mps2":eaw1,
                    "sample1_residual_norm_upper_mps2":rho,"Ktheta_operator_norm_upper":kn,
                    "correction_norm_upper_rad":corr,"inside_9rad_range":closed,"gain_detail":detail,
                }
                rows.append(row)
                if not closed and bad is None: bad=row

    ok=bool(rows) and bad is None and not failures
    return {
        "schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_ANALYTIC_ONE_PLUS_TWO_BLOCK_GAIN_V7",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "direct_first_residual_coordinate_family_retained":True,
        "simultaneous_Rx_inverse_gauge_is_orthogonal":True,
        "Ktheta_operator_norm_invariant_under_block_gauge":True,
        "sample1_innovation_exactly_one_plus_two_block_diagonal":True,
        "scalar_innovation_completed_square_positive_identity_used":True,
        "two_by_two_innovation_positive_determinant_identity_used":True,
        "three_by_three_interval_inverse_used":False,"fixed_pivot_inverse_used":False,
        "spectral_inverse_fallback_used":False,"sample1_nonaxial_force_included":True,
        "full_propagated_aw_error_norm_used":True,"temporal_force_slew_assumed":False,
        "first_attitude_PSD_cross_axis_remainder_included":False,
        "sample1_S_covariance_update_included":False,"sample1_S_attitude_injection_included":False,
        "complete_sample1_branch_closed_here":False,"q8_word_promoted_here":False,
        "whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,
        "candidate_joint_cells_before_residual_ball":len(pcells)*len(rtcells)*len(rzcells),
        "residual_incompatible_joint_cells_pruned":pruned,"evaluated_joint_cells":len(rows),
        "minimum_scalar_innovation_lower":minSx,"minimum_two_by_two_determinant_lower":minDet,
        "max_first_attitude_correction_norm_upper_rad":maxD0,
        "max_first_axial_aw_correction_abs_upper_mps2":maxAz,"max_sample1_force_norm_upper_mps2":maxF,
        "max_post_first_aw_error_upper_mps2":maxEawPost,"max_post_prediction_aw_error_norm_upper_mps2":maxEaw,
        "max_sample1_residual_norm_upper_mps2":maxRho,"max_Ktheta_operator_norm_upper":maxK,
        "max_correction_norm_upper_rad":maxD,"first_unclosed_joint_cell":bad,
        "P5_SAMPLE1_ANALYTIC_BLOCK_GAIN_V7":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":(
            "ADD_FIRST_PSD_CROSS_AXIS_AND_SAMPLE1_S_BRANCH_THEN_SIGNED_CAYLEY_COMPOSE"
            if ok else "REFINE_FIRST_ANALYTIC_BLOCK_WITNESS_WITH_RESIDUAL_DIRECTION_OR_POSITIVE_RATIO_MAXIMIZATION"
        ),
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in (
        "source_generated_not_trajectory_fit","direct_first_residual_coordinate_family_retained",
        "simultaneous_Rx_inverse_gauge_is_orthogonal","Ktheta_operator_norm_invariant_under_block_gauge",
        "sample1_innovation_exactly_one_plus_two_block_diagonal",
        "scalar_innovation_completed_square_positive_identity_used",
        "two_by_two_innovation_positive_determinant_identity_used","sample1_nonaxial_force_included",
        "full_propagated_aw_error_norm_used"):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in (
        "source_replay_used","filter_changed","three_by_three_interval_inverse_used","fixed_pivot_inverse_used",
        "spectral_inverse_fallback_used","temporal_force_slew_assumed",
        "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
        "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here",
        "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no analytic block cells")
    for k in ("minimum_scalar_innovation_lower","minimum_two_by_two_determinant_lower"):
        if not float(d.get(k,0.0))>0.0:f.append(f"nonpositive {k}")
    if not math.isfinite(float(d.get("max_correction_norm_upper_rad",math.nan))):f.append("nonfinite correction")
    st=d.get("P5_SAMPLE1_ANALYTIC_BLOCK_GAIN_V7"); w=d.get("first_unclosed_joint_cell")
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
        "status":d["P5_SAMPLE1_ANALYTIC_BLOCK_GAIN_V7"],"cells":d["evaluated_joint_cells"],
        "pruned":d["residual_incompatible_joint_cells_pruned"],"min_Sx":d["minimum_scalar_innovation_lower"],
        "min_det2":d["minimum_two_by_two_determinant_lower"],"max_force":d["max_sample1_force_norm_upper_mps2"],
        "max_eaw1":d["max_post_prediction_aw_error_norm_upper_mps2"],"max_rho":d["max_sample1_residual_norm_upper_mps2"],
        "max_K":d["max_Ktheta_operator_norm_upper"],"max_d":d["max_correction_norm_upper_rad"],
        "first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf,
    },indent=2,sort_keys=True))
    return 0 if not vf else 2


if __name__=="__main__": raise SystemExit(main())
