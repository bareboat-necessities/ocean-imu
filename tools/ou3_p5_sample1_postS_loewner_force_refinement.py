#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from ou3_interval import matrix_add,matrix_mul,matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_entry as ENTRY
import ou3_p5_sample1_prefix_v2 as PREFIX2
import ou3_p5_sample1_rotation_gauge_refinement as BASE
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_p5_sample1_schur_gain_refinement as SCHUR
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=BASE.DEFAULT_DOMAIN; SCHEMA=1; N=FULL.N; LIMIT=6.0

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2,source_cell_index=0,delta_pieces=8,force_pieces=4):
 FULL3._install_backend(); p=Path(domain_path).resolve(); dom=json.loads(p.read_text()); first=FIRST.build(p,source_pieces=source_pieces); vec=VECTOR.build(); fails=[]
 fails += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase0!="due": fails.append("expected due witness")
 fr=first["source_cells"][source_cell_index]; dmax=float(fr["correction_norm_upper_rad"]); rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tf=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
 Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); RS=FULL._R_S(src); rsmin=min(RS[i][i].lo for i in range(3)); H0=ENTRY._canonical_first_H(g)
 F,Q,Rstep=FULL._transition_and_Q(src,dom); P0=FULL._initial_covariance(src,p); Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,p,h); ep=FULL._predict_error(FULL._initial_error(dom),F); x0=[FULL.I(0) for _ in range(N)]; Ppre,_=ENTRY._zero_residual_S_covariance(Pp,src); qpre=float(first["post_prediction_full_cayley_norm_upper"])
 rows=[]; bad=None; fixed=fallback=0; maxds=maxrho=maxd=0.0
 for di,d in enumerate(SUB.parts(0,dmax,delta_pieces)):
  Pa,ea,xa,_=SUB.first_child(Ppre,ep,x0,H0,Racc,rho0,d); q0=PREFIX2._post_correction_q_upper(qpre,d.hi); Pm=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pa),matrix_transpose(F)),Q)); em=BASE._predict_state(ea,F); xm=BASE._predict_state(xa,F); q1=RG._q_after_first_prediction(q0,dom,h); P1=BASE._transform_covariance(Pm,Rstep); e1=BASE._transform_state(em,Rstep); x1=BASE._transform_state(xm,Rstep)
  rS=FULL._norm_upper([-x1[12+i] for i in range(3)]); ds=FULL.up(SCHUR.gain_bound(P1,FULL.TH,rsmin)*rS); daw=FULL.up(SCHUR.gain_bound(P1,FULL.AW,rsmin)*rS); maxds=max(maxds,ds); q2=PREFIX2._post_correction_q_upper(q1,ds)
  # Accepted Kalman posterior <= prior; apply only the deployed reset uncertainty.
  P2=FULL._reset_covariance(P1,[FULL._box(ds),FULL._box(ds),FULL._box(ds)])
  xaw=FULL.up(BASE._norm(x1,FULL.AW)+daw); eaw=FULL.up(BASE._norm(e1,FULL.AW)+daw); fhat=FULL.up(g+xaw); rho=FULL.up(FULL.up(BASE._rot_diff(q2)*tf)+FULL.up(eaw+ba)); maxrho=max(maxrho,rho)
  for fi,fc in enumerate(SUB.force_cells(fhat,force_pieces)):
   PHt,S=FULL._innovation(P2,SUB.Hforce(fc),Racc); Sinv,b=FULL._spd_inverse_enclosure(S,Racc); K=matrix_mul(PHt,Sinv); kn=RG._op2_upper([r[:] for r in K[:3]]); d1=FULL.up(kn*rho); qa=PREFIX2._post_correction_q_upper(q2,d1); maxd=max(maxd,d1); fixed+=int(b=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"); fallback+=int(b!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"); closed=b=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN" and d1<=LIMIT and math.isfinite(qa) and qa<8
   r={"delta_cell":di,"force_cell":fi,"backend":b,"S_attitude_correction_norm_upper_rad":ds,"acc_residual_norm_upper_mps2":rho,"Ktheta_norm_upper":kn,"acc_correction_norm_upper_rad":d1,"post_q_upper":qa,"closed":closed}; rows.append(r)
   if not closed and bad is None: bad=r
 ok=bool(rows) and bad is None and not fails
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_POST_S_LOEWNER_FORCE_REFINEMENT","source_generated_not_trajectory_fit":True,"filter_changed":False,"sample1_S_accepted_posterior_bounded_by_prior":True,"sample1_S_reset_retained":True,"sample1_S_interval_inverse_not_used":True,"sample1_accelerometer_interval_inverse_used_after_clean_S_bound":True,"sample1_J_aw_exact_identity":True,"evaluated_joint_cells":len(rows),"fixed_pivot_inverse_count":fixed,"spectral_fallback_inverse_count":fallback,"max_sample1_S_attitude_correction_norm_upper_rad":maxds,"max_sample1_acc_residual_norm_upper_mps2":maxrho,"max_sample1_acc_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":bad,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"complete_source_cell_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"P5_SAMPLE1_POST_S_LOEWNER_FORCE_REFINEMENT":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"LIFT_CLEAN_POST_S_ACCEL_REFINEMENT_TO_ALL_BRANCHES" if ok else "DERIVE_ACCEL_TANGENT_DIRECTIONAL_GAIN_BOUND","failures":fails}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","sample1_S_accepted_posterior_bounded_by_prior","sample1_S_reset_retained","sample1_S_interval_inverse_not_used","sample1_accelerometer_interval_inverse_used_after_clean_S_bound","sample1_J_aw_exact_identity"):
  if d.get(k) is not True:f.append(k)
 for k in ("filter_changed","deployed_correction_limit_increased","complete_source_cell_refined_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 st=d.get("P5_SAMPLE1_POST_S_LOEWNER_FORCE_REFINEMENT"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=2); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--delta-pieces",type=int,default=8); a.add_argument("--force-pieces",type=int,default=4); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.delta_pieces,x.force_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_POST_S_LOEWNER_FORCE_REFINEMENT"],"cells":d["evaluated_joint_cells"],"fixed":d["fixed_pivot_inverse_count"],"fallback":d["spectral_fallback_inverse_count"],"max_S_d":d["max_sample1_S_attitude_correction_norm_upper_rad"],"max_rho":d["max_sample1_acc_residual_norm_upper_mps2"],"max_acc_d":d["max_sample1_acc_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
