#!/usr/bin/env python3
"""V12 refinement of the OU-III P5 sample-1 PSD/S perturbation closure.

V11 retained the correct source perturbation magnitudes but bounded the final
accelerometer gain resolvent using only the generic measurement-noise floor
S>=R_a I.  On large-force children that throws away the already certified V7/V8
1+2 innovation information: ||H||^2 DeltaP can nearly consume R_a even though
the actual nominal innovation is much better conditioned.

V12 changes only that denominator.  For every V10 child it reconstructs the
same nominal reduced (theta,a_w) covariance and 1+2 accelerometer innovation
used by V11, then certifies a strict lower eigenvalue with outward interval LDLT
and bisection of S-lambda I.  The exact same V11 PSD cross-axis, second-prediction
process, sample-1 S covariance/reset, S attitude, S a_w and residual perturbation
bounds are retained.  The resolvent now uses

    lambda_pert >= lambda_min(S_nominal) - ||H||^2 DeltaP,

rather than R_a-||H||^2 DeltaP.

This stage is still fail-closed.  It does not compose the signed Cayley state,
promote sample 1 or a whole word, set N_H_words, change the filter, or widen the
deployed six-radian correction primitive.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p5_sample1_structured_full_gain_v11 as V11

DEFAULT_DOMAIN=V11.DEFAULT_DOMAIN
SCHEMA=12
RANGE=V11.RANGE
FULL=V11.FULL


def _shift(A, lam:float):
    z=FULL.I(lam)
    return [[A[i][j]-(z if i==j else FULL.I(0.0)) for j in range(len(A))] for i in range(len(A))]


def _nominal_lambda_min_lower(S)->float:
    """Validated lower bound on lambda_min(S) by interval LDLT bisection."""
    ok,_=symmetric_positive_definite_ldlt(S)
    if not ok:
        return 0.0
    hi=min(float(S[i][i].lo) for i in range(len(S)))
    if not hi>0.0:
        return 0.0
    lo=0.0
    for _ in range(52):
        mid=0.5*(lo+hi)
        ok,_=symmetric_positive_definite_ldlt(_shift(S,mid))
        if ok:
            lo=mid
        else:
            hi=mid
    return FULL.down(lo)


def build(domain_path:Path=DEFAULT_DOMAIN,*,source_pieces:int=4,source_cell_index:int=0,
          p_pieces:int=24,tangent_pieces:int=24,axial_pieces:int=24)->dict:
    # Certification-only refinement: proof domain and shipping correction limit are unchanged.
    V11.FULL3._install_backend(); path=Path(domain_path).resolve(); dom=json.loads(path.read_text(encoding="utf-8"))
    core=V11.V10.build(path,source_pieces=source_pieces,source_cell_index=source_cell_index,
                       p_pieces=p_pieces,tangent_pieces=tangent_pieces,axial_pieces=axial_pieces)
    first=V11.FIRST.build(path,source_pieces=source_pieces); vec=V11.VECTOR.build(); failures=[]
    failures += [f"V10: {x}" for x in V11.V10.validate(core)]
    failures += [f"first: {x}" for x in V11.FIRST.validate(first)]
    failures += [f"vector: {x}" for x in V11.VECTOR.validate(vec)]
    if core.get("P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10")!="PASS":
        failures.append("V10 prerequisite did not pass")

    src,phase=V11.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase!="due": failures.append("V12 focused perturbation requires first due source cell")
    if source_cell_index!=0: failures.append("V12 S perturbation helper currently certified for source cell 0")

    fr=first["source_cells"][source_cell_index]
    p_all=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"]))
    rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
    aw_pre=float(fr["predicted_aw_error_norm_upper_mps2"])
    h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"])
    tilt,yaw,eps=V11.RG._attitude_covariance_epsilon(path,h)
    t=Interval.outward_bounds(tilt,FULL.up(tilt+eps)); Y=Interval.outward_bounds(yaw,FULL.up(yaw+eps))
    r=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F,Q,_=FULL._transition_and_Q(src,dom); alpha=F[15][15]; alpha_hi=alpha.hi; qaw=Q[15][15]
    pcells=V11.SUB.parts(p_all.lo,p_all.hi,p_pieces)
    sb=V11._sample1_s_bounds(path,src,dom,first,vec)

    rows=[]; bad=None; worst=None; unclosed=0
    maxd=maxdp=maxdr=maxdk=maxlamloss=0.0
    min_nominal=math.inf; min_perturbed=math.inf
    for base in core["rows"]:
        pi=int(base["p_cell"]); p=pcells[pi]
        rt=Interval.outward_bounds(*map(float,base["first_tangent_residual_magnitude_mps2"]))
        rz=Interval.outward_bounds(*map(float,base["first_axial_residual_mps2"]))
        d=Interval.outward_bounds(*map(float,base["first_attitude_correction_rad"]))
        D=FULL.I(g*g)*t+p+r
        fy=-(alpha*(p/D)*rt); fz=FULL.I(g)+alpha*(p/(p+r))*rz
        psd=V11._first_psd_perturbation(
            t=t,Y=Y,p=p,r=r,g=g,eps=eps,rho0=rho0,dhi=max(0.0,d.hi),
            rt=rt,rz=rz,alpha_hi=alpha_hi,aw_pre=aw_pre)

        vecnorm=FULL.up(g+FULL.up(alpha_hi*(aw_pre+psd["first_nominal_aw_correction_norm_upper_mps2"]+psd["first_aw_x_correction_upper_mps2"])))
        drS=FULL.up(sb["sample1_S_aw_correction_upper_mps2"]+FULL.up(sb["sample1_S_attitude_correction_upper_rad"]*vecnorm))
        drho=FULL.up(psd["PSD_induced_sample1_residual_perturbation_upper_mps2"]+drS)
        dP=FULL.up(psd["sample1_reduced_covariance_PSD_perturbation_upper"]+sb["sample1_S_total_reduced_covariance_perturbation_upper"])

        Pn,Hn,Sn=V11._nominal_sample1_matrices(t=t,Y=Y,p=p,r=r,g=g,alpha=alpha,qaw=qaw,d=d,fy=fy,fz=fz)
        h1=V11._op(Hn); pnorm=V11._op(Pn)
        lam=_nominal_lambda_min_lower(Sn)
        dS=FULL.up(FULL.up(h1*h1)*dP)
        lamp=FULL.down(lam-dS)
        min_nominal=min(min_nominal,lam)
        if not lamp>0.0:
            row={
                "reason":"perturbed innovation lower bound nonpositive",
                "nominal_innovation_lambda_lower":lam,"innovation_perturbation_upper":dS,
                "perturbed_innovation_lambda_lower":lamp,"p_cell":pi,
                "tangent_residual_cell":base["tangent_residual_cell"],"axial_residual_cell":base["axial_residual_cell"],
            }
            rows.append(row); unclosed+=1
            if bad is None: bad=row
            continue

        k0=max(float(base["Ktheta_perpendicular_block_upper"]),float(base["Ktheta_parallel_block_upper"]))
        dC=FULL.up(dP*h1)
        dk=FULL.up(FULL.up(dC+FULL.up(k0*dS))/lamp)
        rho=float(base["sample1_full_residual_norm_upper_mps2"])
        d10=float(base["combined_directional_correction_norm_upper_rad"])
        d12=FULL.up(d10+FULL.up(FULL.up(k0*drho)+FULL.up(dk*FULL.up(rho+drho))))
        closed=math.isfinite(d12) and d12<RANGE
        row={
            "p_cell":pi,"tangent_residual_cell":base["tangent_residual_cell"],"axial_residual_cell":base["axial_residual_cell"],
            "V10_directional_correction_upper_rad":d10,
            "PSD_residual_perturbation_upper_mps2":psd["PSD_induced_sample1_residual_perturbation_upper_mps2"],
            "S_residual_perturbation_upper_mps2":drS,"total_residual_perturbation_upper_mps2":drho,
            "PSD_reduced_covariance_perturbation_upper":psd["sample1_reduced_covariance_PSD_perturbation_upper"],
            "S_reduced_covariance_perturbation_upper":sb["sample1_S_total_reduced_covariance_perturbation_upper"],
            "total_reduced_covariance_perturbation_upper":dP,"nominal_reduced_covariance_operator_upper":pnorm,
            "sample1_H_operator_upper":h1,"nominal_innovation_lambda_lower":lam,
            "measurement_noise_floor_lower":r.lo,"innovation_perturbation_upper":dS,
            "perturbed_innovation_lambda_lower":lamp,"sample1_gain_operator_perturbation_upper":dk,
            "V12_correction_norm_upper_rad":d12,"inside_9rad_range":closed,**psd,
        }
        rows.append(row)
        maxd=max(maxd,d12); maxdp=max(maxdp,dP); maxdr=max(maxdr,drho); maxdk=max(maxdk,dk); maxlamloss=max(maxlamloss,dS)
        min_perturbed=min(min_perturbed,lamp)
        if worst is None or d12>worst.get("V12_correction_norm_upper_rad",-math.inf): worst=row
        if not closed:
            unclosed+=1
            if bad is None: bad=row

    ok=bool(rows) and unclosed==0 and bad is None and not failures
    return {
        "schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_RESOLVENT_V12",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "V10_canonical_core_retained":True,"V11_PSD_and_S_perturbation_magnitudes_retained":True,
        "actual_nominal_innovation_lambda_certified_by_interval_LDLT":True,
        "measurement_noise_only_resolvent_floor_used":False,
        "first_attitude_PSD_cross_axis_remainder_included":True,
        "second_prediction_attitude_process_remainder_included":True,
        "sample1_S_covariance_update_included":True,"sample1_S_attitude_injection_included":True,
        "sample1_S_aw_mean_correction_included":True,"sample1_S_solver_identity_branch_contained_as_zero_perturbation":True,
        "broad_sample1_3x3_interval_inverse_reintroduced":False,"temporal_force_slew_assumed":False,
        "complete_sample1_branch_closed_here":False,"signed_cayley_q8_composed_here":False,
        "q8_word_promoted_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,"deployed_correction_limit_rad":6.0,"deployed_correction_limit_increased":False,
        "attitude_covariance_remainder_spectral_upper":eps,"sample1_S_perturbation_bounds":sb,
        "evaluated_joint_cells":len(rows),"unclosed_joint_cells":unclosed,
        "minimum_nominal_innovation_lambda_lower":min_nominal,"max_innovation_perturbation_upper":maxlamloss,
        "minimum_perturbed_innovation_lambda_lower":min_perturbed,
        "max_total_residual_perturbation_upper_mps2":maxdr,"max_total_reduced_covariance_perturbation_upper":maxdp,
        "max_sample1_gain_operator_perturbation_upper":maxdk,"max_V12_correction_norm_upper_rad":maxd,
        "first_unclosed_joint_cell":bad,"worst_joint_cell":worst,
        "P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":"SIGNED_CAYLEY_COMPOSE_SAMPLE1_S_AND_ACCELERATOR_INSIDE_Q8" if ok else "REFINE_BLOCKWISE_PERTURBATION_GAIN_NUMERATOR_AT_FIRST_WITNESS",
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in (
        "source_generated_not_trajectory_fit","V10_canonical_core_retained","V11_PSD_and_S_perturbation_magnitudes_retained",
        "actual_nominal_innovation_lambda_certified_by_interval_LDLT","first_attitude_PSD_cross_axis_remainder_included",
        "second_prediction_attitude_process_remainder_included","sample1_S_covariance_update_included",
        "sample1_S_attitude_injection_included","sample1_S_aw_mean_correction_included",
        "sample1_S_solver_identity_branch_contained_as_zero_perturbation",
    ):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in (
        "source_replay_used","filter_changed","measurement_noise_only_resolvent_floor_used",
        "broad_sample1_3x3_interval_inverse_reintroduced","temporal_force_slew_assumed",
        "complete_sample1_branch_closed_here","signed_cayley_q8_composed_here","q8_word_promoted_here",
        "whole_word_promoted_here","N_H_words_set_here","deployed_correction_limit_increased",
    ):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad",0.0))!=6.0:f.append("deployed correction limit changed")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no V12 cells")
    for k in ("minimum_nominal_innovation_lambda_lower","max_total_residual_perturbation_upper_mps2",
              "max_total_reduced_covariance_perturbation_upper","max_sample1_gain_operator_perturbation_upper","max_V12_correction_norm_upper_rad"):
        if not math.isfinite(float(d.get(k,math.nan))):f.append(f"nonfinite {k}")
    st=d.get("P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and (w is not None or int(d.get("unclosed_joint_cells",-1))!=0):f.append("PASS retains V12 unclosed cell")
    if st=="NOT_ESTABLISHED" and w is None and not f:f.append("missing V12 witness")
    return f


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces",type=int,default=4); ap.add_argument("--source-cell-index",type=int,default=0)
    ap.add_argument("--p-pieces",type=int,default=24); ap.add_argument("--tangent-pieces",type=int,default=24); ap.add_argument("--axial-pieces",type=int,default=24)
    ap.add_argument("--output",type=Path,required=True); x=ap.parse_args()
    d=build(x.domain,source_pieces=x.source_pieces,source_cell_index=x.source_cell_index,
            p_pieces=x.p_pieces,tangent_pieces=x.tangent_pieces,axial_pieces=x.axial_pieces)
    vf=validate(d); d["validation_failures"]=vf
    x.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({
        "status":d["P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12"],"cells":d["evaluated_joint_cells"],"unclosed":d["unclosed_joint_cells"],
        "min_nominal_lambda":d["minimum_nominal_innovation_lambda_lower"],"max_dS":d["max_innovation_perturbation_upper"],
        "min_perturbed_lambda":d["minimum_perturbed_innovation_lambda_lower"],"max_drho":d["max_total_residual_perturbation_upper_mps2"],
        "max_dP":d["max_total_reduced_covariance_perturbation_upper"],"max_dK":d["max_sample1_gain_operator_perturbation_upper"],
        "max_d":d["max_V12_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],
        "worst":d["worst_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True))
    return 0 if not vf else 2

if __name__=="__main__": raise SystemExit(main())
