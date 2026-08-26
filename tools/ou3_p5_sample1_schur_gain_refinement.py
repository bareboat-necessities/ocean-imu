#!/usr/bin/env python3
"""Bound the first sample-1 P5 witness without interval inversion.

For an accepted shipping Kalman update K=P H' S^-1,
    K S K' <= P,  S >= R > 0.
Hence for any state group g,
    ||K_g||_2 <= sqrt(lambda_max(P_gg)/lambda_min(R)).
If safe_ldlt3_ fails, shipping performs the identity/no-correction branch, so
zero correction is already contained.  This producer uses that dichotomy for
the sample-1 S and accelerometer updates and never needs an interval inverse.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval,matrix_add,matrix_mul,matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_entry as ENTRY
import ou3_p5_sample1_prefix_v2 as PREFIX2
import ou3_p5_sample1_rotation_gauge_refinement as BASE
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=BASE.DEFAULT_DOMAIN; SCHEMA=1; N=FULL.N; LIMIT=6.0

def lam_upper(P,idxs):
    ans=0.0
    for i in idxs:
        row=max(0.0,P[i][i].hi)
        for j in idxs:
            if j!=i: row+=max(abs(P[i][j].lo),abs(P[i][j].hi))
        ans=max(ans,FULL.up(row))
    return ans

def gain_bound(P,idxs,rmin):
    if not rmin>0: raise RuntimeError("measurement floor lost")
    return FULL.up(math.sqrt(FULL.up(lam_upper(P,idxs)/rmin)))

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2,source_cell_index=0,delta_pieces=8):
    FULL3._install_backend(); p=Path(domain_path).resolve(); dom=json.loads(p.read_text())
    first=FIRST.build(p,source_pieces=source_pieces); vec=VECTOR.build(); fails=[]
    fails += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
    src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
    if phase0!="due": fails.append("expected first due witness")
    fr=first["source_cells"][source_cell_index]; dmax=float(fr["correction_norm_upper_rad"]); rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
    h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tf=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); racc=min(Racc[i][i].lo for i in range(3)); RS=FULL._R_S(src); rsmin=min(RS[i][i].lo for i in range(3)); H0=ENTRY._canonical_first_H(g)
    F,Q,Rstep=FULL._transition_and_Q(src,dom); P0=FULL._initial_covariance(src,p); Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,p,h)
    ep=FULL._predict_error(FULL._initial_error(dom),F); x0=[FULL.I(0) for _ in range(N)]; Ppre,_=ENTRY._zero_residual_S_covariance(Pp,src); qpre=float(first["post_prediction_full_cayley_norm_upper"])
    rows=[]; bad=None; maxds=maxdaw=maxrho=maxd1=0.0
    for di,d in enumerate(SUB.parts(0,dmax,delta_pieces)):
        Pa,ea,xa,_=SUB.first_child(Ppre,ep,x0,H0,Racc,rho0,d); q0=PREFIX2._post_correction_q_upper(qpre,d.hi)
        Pm=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pa),matrix_transpose(F)),Q)); em=BASE._predict_state(ea,F); xm=BASE._predict_state(xa,F); q1=RG._q_after_first_prediction(q0,dom,h)
        P1=BASE._transform_covariance(Pm,Rstep); e1=BASE._transform_state(em,Rstep); x1=BASE._transform_state(xm,Rstep)
        rS=FULL._norm_upper([-x1[12+i] for i in range(3)]); kSth=gain_bound(P1,FULL.TH,rsmin); kSaw=gain_bound(P1,FULL.AW,rsmin); ds=FULL.up(kSth*rS); daw=FULL.up(kSaw*rS)
        maxds=max(maxds,ds); maxdaw=max(maxdaw,daw); q2=PREFIX2._post_correction_q_upper(q1,ds)
        xaw=FULL.up(BASE._norm(x1,FULL.AW)+daw); eaw=FULL.up(BASE._norm(e1,FULL.AW)+daw); rho=FULL.up(FULL.up(BASE._rot_diff(q2)*tf)+FULL.up(eaw+ba)); maxrho=max(maxrho,rho)
        # Accepted S posterior is <= P1.  Only the attitude reset can enlarge
        # lambda_max(P_theta); ||G(d)||^2 = 1+||d||^2/4.
        lam1=lam_upper(P1,FULL.TH); lam2=FULL.up(lam1*FULL.up(1.0+0.25*ds*ds)); kA=FULL.up(math.sqrt(FULL.up(lam2/racc))); d1=FULL.up(kA*rho); qa=PREFIX2._post_correction_q_upper(q2,d1); maxd1=max(maxd1,d1)
        closed=ds<=LIMIT and d1<=LIMIT and math.isfinite(qa) and qa<8
        r={"delta_cell":di,"delta_rad":d.as_list(),"sample1_S_residual_norm_upper":rS,"sample1_S_Ktheta_norm_upper":kSth,"sample1_S_attitude_correction_norm_upper_rad":ds,"sample1_S_aw_correction_norm_upper_mps2":daw,"sample1_acc_residual_norm_upper_mps2":rho,"sample1_acc_Ktheta_norm_upper":kA,"sample1_acc_correction_norm_upper_rad":d1,"sample1_post_accel_q_upper":qa,"closed":closed}; rows.append(r)
        if not closed and bad is None: bad=r
    ok=bool(rows) and bad is None and not fails
    return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_SCHUR_GAIN_WITNESS_REFINEMENT","source_generated_not_trajectory_fit":True,"filter_changed":False,"accepted_gain_uses_KSK_le_P":True,"solver_failure_identity_branch_included":True,"no_sample1_interval_inverse_used":True,"S_posterior_loewner_upper_by_prior_used":True,"reset_operator_exact_norm_formula_used":True,"delta_pieces":delta_pieces,"evaluated_delta_cells":len(rows),"max_sample1_S_attitude_correction_norm_upper_rad":maxds,"max_sample1_S_aw_correction_norm_upper_mps2":maxdaw,"max_sample1_acc_residual_norm_upper_mps2":maxrho,"max_sample1_acc_correction_norm_upper_rad":maxd1,"first_unclosed_delta_cell":bad,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"complete_source_cell_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"P5_SAMPLE1_SCHUR_GAIN_WITNESS_REFINEMENT":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"LIFT_SCHUR_GAIN_REFINEMENT_TO_ALL_SAMPLE1_BRANCHES_AND_SOURCE_CELLS" if ok else "TIGHTEN_POST_RESET_ATTITUDE_COVARIANCE_OR_RESIDUAL_CHART","failures":fails,"rows":rows}

def validate(d):
    f=list(d.get("failures",[]))
    for k in ("source_generated_not_trajectory_fit","accepted_gain_uses_KSK_le_P","solver_failure_identity_branch_included","no_sample1_interval_inverse_used","S_posterior_loewner_upper_by_prior_used","reset_operator_exact_norm_formula_used"):
        if d.get(k) is not True:f.append(k)
    for k in ("filter_changed","deployed_correction_limit_increased","complete_source_cell_refined_here","whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False:f.append(k)
    st=d.get("P5_SAMPLE1_SCHUR_GAIN_WITNESS_REFINEMENT"); w=d.get("first_unclosed_delta_cell")
    if st=="PASS" and w is not None:f.append("PASS retains witness")
    if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
    return f

def main():
    a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=2); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--delta-pieces",type=int,default=8); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.delta_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_SCHUR_GAIN_WITNESS_REFINEMENT"],"delta_cells":d["evaluated_delta_cells"],"max_S_d":d["max_sample1_S_attitude_correction_norm_upper_rad"],"max_S_aw":d["max_sample1_S_aw_correction_norm_upper_mps2"],"max_rho":d["max_sample1_acc_residual_norm_upper_mps2"],"max_acc_d":d["max_sample1_acc_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_delta_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
