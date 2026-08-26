#!/usr/bin/env python3
"""Positive-ratio V8 tightening of the analytic 1+2 sample-1 gain.

V7 removed every matrix-inverse fallback and exposed a genuine analytic first
witness at 9.609 rad.  Its gain bound still separated each block's numerator
maximum from its denominator minimum.  Both blocks admit the same exact
positive-ratio form, so that separation is unnecessary.

For nonnegative x,z and A>0 define

    F(x,z)=(c1*x+c2*z)/(x+z+A)^2.

On a rectangle [x0,x1]x[z0,z1], an interior stationary point exists only when
c1=c2; in that case the stationary line x+z=A intersects the rectangle boundary
whenever it intersects the rectangle.  Therefore the exact rectangle maximum is
attained at a corner or at the single one-dimensional stationary point on one
of the four edges.  We evaluate all such candidates with outward-rounded input
bounds.

V7 scalar block:

    x=N_u^2/a, z=Y V^2, A=Delta/a+r,
    ||K_yz,x||^2=(1+h^2)(a*x+Y*z)/(x+z+A)^2.

V7 2x2 block, with R=B_z+r and A2=Delta+a r:

    u=q^2/A2, v=a f_y^2/R,
    ||K_x,yz||^2=[(a^2/A2)u+(a/R)v]/(u+v+1)^2.

V8 maximizes these two forms directly instead of dividing independent interval
numerator/denominator extrema.  The direct first-residual coordinates, source
force/error bounds, and all no-promotion obligations are unchanged from V7.
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
import ou3_p5_sample1_structured_full_gain_v7 as V7
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN
SCHEMA=8
RANGE=9.0


def _rect(iv:Interval):
    return max(0.0,float(iv.lo)),max(0.0,float(iv.hi))


def _ratio(c1:float,c2:float,x:float,z:float,A:float)->float:
    den=FULL.down(FULL.down(x+z)+A)
    if not den>0.0: raise RuntimeError("positive ratio denominator lost")
    num=FULL.up(FULL.up(c1*x)+FULL.up(c2*z))
    return FULL.up(num/FULL.down(den*den))


def _linear_over_square_rect_max(c1:float,c2:float,xiv:Interval,ziv:Interval,A:float):
    """Upper max of (c1*x+c2*z)/(x+z+A)^2 on a nonnegative rectangle."""
    if not (c1>=0.0 and c2>=0.0 and A>0.0):
        raise RuntimeError("invalid positive-ratio coefficients")
    x0,x1=_rect(xiv); z0,z1=_rect(ziv)
    cand=[]
    def add(x,z,kind):
        if x0<=x<=x1 and z0<=z<=z1:
            cand.append((_ratio(c1,c2,x,z,A),x,z,kind))
    for x in (x0,x1):
        for z in (z0,z1): add(x,z,"corner")
    if c1>0.0:
        for z in (z0,z1):
            xs=FULL.up(A+z*(1.0-2.0*c2/c1))
            add(min(max(xs,x0),x1),z,"x_edge_stationary_or_clamped")
    if c2>0.0:
        for x in (x0,x1):
            zs=FULL.up(A+x*(1.0-2.0*c1/c2))
            add(x,min(max(zs,z0),z1),"z_edge_stationary_or_clamped")
    if not cand: raise RuntimeError("empty ratio candidate set")
    best=max(cand,key=lambda q:q[0])
    return FULL.up(best[0]),{"x_range":[x0,x1],"z_range":[z0,z1],"maximizer":[best[1],best[2]],"maximizer_kind":best[3],"ratio_upper":FULL.up(best[0])}


def _block_gain_bounds(a:Interval,Y:Interval,c0:Interval,alpha:Interval,qaw:Interval,
                       b:Interval,bz:Interval,det_first:Interval,d:Interval,
                       fy:Interval,fz:Interval,r:Interval):
    if not (a.lo>0.0 and Y.lo>0.0 and r.lo>0.0 and det_first.lo>0.0):
        raise RuntimeError("positive canonical covariance/noise floors required")
    h=FULL.I(0.5)*d; one=FULL.I(1.0); C=one+h.square()
    Bt=alpha.square()*b+qaw; Bz=alpha.square()*bz+qaw
    delta=alpha.square()*det_first+a*qaw
    if not (Bt.lo>0.0 and Bz.lo>0.0 and delta.lo>0.0):
        raise RuntimeError("positive predicted determinant lost")

    U=fz-h*fy; V=-(h*fz)-fy; cu=alpha*c0; Nu=a*U+cu
    x=Nu.square()/a; z=Y*V.square(); A=delta/a+r
    if not A.lo>0.0: raise RuntimeError("scalar A floor lost")
    f1,d1=_linear_over_square_rect_max(a.hi,Y.hi,x,z,A.lo)
    kperp=FULL.up(math.sqrt(max(0.0,FULL.up(C.hi*f1))))

    cx=-(alpha*c0); q=cx-a*fz; A2=delta+a*r; R=Bz+r
    if not (A2.lo>0.0 and R.lo>0.0): raise RuntimeError("2x2 A/R floor lost")
    u=q.square()/A2; v=a*fy.square()/R
    c1=FULL.up((a.square()/A2).hi); c2=FULL.up((a/R).hi)
    f2,d2=_linear_over_square_rect_max(c1,c2,u,v,1.0)
    kx=FULL.up(math.sqrt(max(0.0,f2)))

    # Retain exact positive innovation diagnostics from V7 as independent
    # algebra checks; they are not used as separated gain denominators.
    Sx=x+z+A
    det2=fy.square()*A2+(q.square()+A2)*R/a
    if not (Sx.lo>0.0 and det2.lo>0.0): raise RuntimeError("positive innovation diagnostic lost")
    return max(kperp,kx),{
        "half_correction":h.as_list(),"Bt":Bt.as_list(),"Bz":Bz.as_list(),"Delta":delta.as_list(),
        "scalar_x":x.as_list(),"scalar_z":z.as_list(),"scalar_A":A.as_list(),
        "scalar_positive_S":Sx.as_list(),"scalar_ratio_maximization":d1,
        "scalar_gain_yz_norm_upper":kperp,
        "two_by_two_u":u.as_list(),"two_by_two_v":v.as_list(),"two_by_two_A2":A2.as_list(),
        "two_by_two_R":R.as_list(),"two_by_two_positive_det":det2.as_list(),
        "two_by_two_ratio_maximization":d2,"two_by_two_theta_x_gain_norm_upper":kx,
    }


def build(domain_path:Path=DEFAULT_DOMAIN,*,source_pieces:int=4,source_cell_index:int=0,
          p_pieces:int=24,tangent_pieces:int=24,axial_pieces:int=24)->dict:
    FULL3._install_backend(); path=Path(domain_path).resolve(); dom=json.loads(path.read_text(encoding="utf-8"))
    first=FIRST.build(path,source_pieces=source_pieces); vec=VECTOR.build(); failures=[]
    failures += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
    src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
    if phase!="due": failures.append("V8 ratio witness requires first due source cell")
    fr=first["source_cells"][source_cell_index]; p_all=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"]))
    rho0=float(fr["combined_useful_residual_norm_upper_mps2"]); aw_pred=float(fr["predicted_aw_error_norm_upper_mps2"])
    hs=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"])
    ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    tilt,yaw,eps=RG._attitude_covariance_epsilon(path,hs); t=Interval.outward_bounds(tilt,FULL.up(tilt+eps)); Y=Interval.outward_bounds(yaw,FULL.up(yaw+eps))
    r=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F,Q,_=FULL._transition_and_Q(src,dom); alpha=F[15][15]; alpha_hi=alpha.hi; qaw=Q[15][15]
    cos0=float(first["post_prediction_true_gravity_cosine_lower"])
    sin_hi=1.0 if cos0<0.0 else FULL.up(math.sqrt(max(0.0,FULL.up(1.0-FULL.down(cos0*cos0)))))
    yRt=FULL.up(g*sin_hi); yRz=FULL.up(g*max(0.0,FULL.up(1.0-cos0)))
    rt_cap=min(rho0,FULL.up(aw_pred+FULL.up(yRt+ba))); rz_cap=min(rho0,FULL.up(aw_pred+FULL.up(yRz+ba)))
    chord0=V4._gravity_chord_from_cos(cos0); pred_chord=V4._correction_chord_upper(float(first["first_prediction_transport_angle_upper_rad"]))
    pcells=SUB.parts(p_all.lo,p_all.hi,p_pieces); rtcells=SUB.parts(0.0,rt_cap,tangent_pieces); rzcells=SUB.parts(-rz_cap,rz_cap,axial_pieces)
    rows=[]; bad=None; pruned=0; maxK=maxD=maxRho=maxF=maxEaw=0.0; minS=minDet=math.inf
    for pi,p in enumerate(pcells):
        D=FULL.I(g*g)*t+p+r; a=t*(p+r)/D; c0=-(FULL.I(g)*t*p/D)
        b=p*(FULL.I(g*g)*t+r)/D; bz=p*r/(p+r); det_first=t*p*r/D
        ktheta=FULL.I(g)*t/D; kaw_t=p/D; kz=p/(p+r)
        for ti,rt0 in enumerate(rtcells):
            for zi,rz0 in enumerate(rzcells):
                child=V6._residual_child(rt0,rz0,rho0)
                if child is None: pruned+=1; continue
                rt,rz=child; d=ktheta*rt; awt=kaw_t*rt; az=kz*rz
                fy=-(alpha*awt); fz=FULL.I(g)+alpha*az; fn=V5._norm2_upper(fy.abs_upper(),fz.abs_upper()); maxF=max(maxF,fn)
                kn,detail=_block_gain_bounds(a,Y,c0,alpha,qaw,b,bz,det_first,d,fy,fz,r)
                minS=min(minS,float(detail["scalar_positive_S"][0])); minDet=min(minDet,float(detail["two_by_two_positive_det"][0]))
                dhi=max(0.0,d.hi); chord=min(2.0,FULL.up(chord0+FULL.up(V4._correction_chord_upper(dhi)+pred_chord)))
                one_t=FULL.I(1.0)-kaw_t; one_z=FULL.I(1.0)-kz
                left_t=(one_t*rt).abs_upper(); left_z=(one_z*rz).abs_upper()
                post_r=FULL.up(V5._norm2_upper(FULL.up(left_t+yRt),FULL.up(left_z+yRz))+ba)
                post_tri=FULL.up(aw_pred+V5._norm2_upper(awt.abs_upper(),az.abs_upper())); post=min(post_r,post_tri)
                eaw1=FULL.up(alpha_hi*post); maxEaw=max(maxEaw,eaw1)
                rho=FULL.up(FULL.up(fn*chord)+FULL.up(eaw1+ba)); corr=FULL.up(kn*rho); closed=math.isfinite(corr) and corr<RANGE
                maxK=max(maxK,kn); maxD=max(maxD,corr); maxRho=max(maxRho,rho)
                row={"p_cell":pi,"tangent_residual_cell":ti,"axial_residual_cell":zi,
                     "first_tangent_residual_magnitude_mps2":rt.as_list(),"first_axial_residual_mps2":rz.as_list(),
                     "first_attitude_correction_rad":d.as_list(),"sample1_force_unRx_components_yz_mps2":[fy.as_list(),fz.as_list()],
                     "sample1_force_norm_upper_mps2":fn,"post_prediction_aw_error_norm_upper_mps2":eaw1,
                     "sample1_residual_norm_upper_mps2":rho,"Ktheta_operator_norm_upper":kn,
                     "correction_norm_upper_rad":corr,"inside_9rad_range":closed,"gain_detail":detail}
                rows.append(row)
                if not closed and bad is None: bad=row
    ok=bool(rows) and bad is None and not failures
    return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_POSITIVE_RATIO_ANALYTIC_BLOCK_GAIN_V8",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "direct_first_residual_coordinate_family_retained":True,"analytic_one_plus_two_block_structure_retained":True,
        "block_numerator_denominator_dependency_preserved_by_positive_ratio_maximization":True,
        "rectangle_boundary_stationary_maximization_used":True,"three_by_three_interval_inverse_used":False,
        "spectral_inverse_fallback_used":False,"sample1_nonaxial_force_included":True,"full_propagated_aw_error_norm_used":True,
        "first_attitude_PSD_cross_axis_remainder_included":False,"sample1_S_covariance_update_included":False,
        "sample1_S_attitude_injection_included":False,"complete_sample1_branch_closed_here":False,
        "q8_word_promoted_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,"evaluated_joint_cells":len(rows),
        "residual_incompatible_joint_cells_pruned":pruned,"minimum_scalar_innovation_lower":minS,
        "minimum_two_by_two_determinant_lower":minDet,"max_sample1_force_norm_upper_mps2":maxF,
        "max_post_prediction_aw_error_norm_upper_mps2":maxEaw,"max_sample1_residual_norm_upper_mps2":maxRho,
        "max_Ktheta_operator_norm_upper":maxK,"max_correction_norm_upper_rad":maxD,"first_unclosed_joint_cell":bad,
        "P5_SAMPLE1_POSITIVE_RATIO_BLOCK_GAIN_V8":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":"ADD_FIRST_PSD_CROSS_AXIS_AND_SAMPLE1_S_BRANCH_THEN_SIGNED_CAYLEY_COMPOSE" if ok else "REFINE_FIRST_RATIO_WITNESS_WITH_RESIDUAL_DIRECTION_OR_MORE_SUBDIVISION",
        "failures":failures,"rows":rows}


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in ("source_generated_not_trajectory_fit","direct_first_residual_coordinate_family_retained",
              "analytic_one_plus_two_block_structure_retained","block_numerator_denominator_dependency_preserved_by_positive_ratio_maximization",
              "rectangle_boundary_stationary_maximization_used","sample1_nonaxial_force_included","full_propagated_aw_error_norm_used"):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in ("source_replay_used","filter_changed","three_by_three_interval_inverse_used","spectral_inverse_fallback_used",
              "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
              "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here","q8_word_promoted_here",
              "whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no V8 cells")
    for k in ("minimum_scalar_innovation_lower","minimum_two_by_two_determinant_lower"):
        if not float(d.get(k,0.0))>0.0:f.append(f"nonpositive {k}")
    if not math.isfinite(float(d.get("max_correction_norm_upper_rad",math.nan))):f.append("nonfinite correction")
    st=d.get("P5_SAMPLE1_POSITIVE_RATIO_BLOCK_GAIN_V8"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and w is not None:f.append("PASS retains witness")
    if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
    return f


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces",type=int,default=4); ap.add_argument("--source-cell-index",type=int,default=0)
    ap.add_argument("--p-pieces",type=int,default=24); ap.add_argument("--tangent-pieces",type=int,default=24)
    ap.add_argument("--axial-pieces",type=int,default=24); ap.add_argument("--output",type=Path,required=True)
    x=ap.parse_args(); d=build(x.domain,source_pieces=x.source_pieces,source_cell_index=x.source_cell_index,p_pieces=x.p_pieces,tangent_pieces=x.tangent_pieces,axial_pieces=x.axial_pieces)
    vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"status":d["P5_SAMPLE1_POSITIVE_RATIO_BLOCK_GAIN_V8"],"cells":d["evaluated_joint_cells"],
        "pruned":d["residual_incompatible_joint_cells_pruned"],"min_Sx":d["minimum_scalar_innovation_lower"],
        "min_det2":d["minimum_two_by_two_determinant_lower"],"max_force":d["max_sample1_force_norm_upper_mps2"],
        "max_eaw1":d["max_post_prediction_aw_error_norm_upper_mps2"],"max_rho":d["max_sample1_residual_norm_upper_mps2"],
        "max_K":d["max_Ktheta_operator_norm_upper"],"max_d":d["max_correction_norm_upper_rad"],
        "first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True))
    return 0 if not vf else 2

if __name__=="__main__": raise SystemExit(main())
