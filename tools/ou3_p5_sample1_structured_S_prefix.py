#!/usr/bin/env python3
"""Structured consecutive sample-1 S=0 prefix diagnostic for OU-III P5.

This stage consumes the exact first-accelerometer Joseph posterior from the
transported-crosscovariance backend. For each canonical first correction cell
it applies the shipping reset, exact body-gauge coordinate transform, one 5-ms
prediction, and source-enclosed body rotation, then evaluates the sample-1 due
S=0 measurement with the shipping 18x18 covariance and actual estimator S
residual.

H_S is constant and sparse. Fixed-pivot inversion is used where certified; the
existing rigorous SPD S>=R inverse enclosure is also admissible and reported.
Unlike the second-accelerometer diagnostic, use of that fallback is not itself
a failure: the gate is whether its resulting S attitude correction remains
finite and inside q<8. The due solver-identity branch has zero correction. The
second accelerometer is not evaluated here, so no complete sample-1/P5
promotion is made.
"""
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
import ou3_p5_sample1_transported_crosscov_refinement as TC
import ou3_vector_uco_certificate as VECTOR
DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN; SCHEMA=2

def build(domain_path=DEFAULT_DOMAIN,source_pieces=4,source_cell_index=0,d_pieces=32):
 FULL3._install_backend(); path=Path(domain_path).resolve(); domain=json.loads(path.read_text()); first=FIRST.build(path,source_pieces=source_pieces); vector=VECTOR.build(); failures=[]
 failures += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vector)]
 src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase0!="due": failures.append("structured consecutive-S witness requires sample-0 due source cell")
 phase1=PREFIX2.BASE._phase1_options(src,phase0)
 if phase1 != ["due"]: failures.append(f"sample-0 due source cell did not force sample-1 due: {phase1!r}")
 fr=first["source_cells"][source_cell_index]; dmax=float(fr["correction_norm_upper_rad"]); rho0=float(fr["combined_useful_residual_norm_upper_mps2"]); h=float(src["dt_s"]); g=float(domain["startup"]["gravity_mps2"])
 Racc=FULL._R_diag(float(vector["configured_measurement_bounds"]["acc_measurement_std_mps2"])); H0=ENTRY._canonical_first_H(g); F,Q,Rstep=FULL._transition_and_Q(src,domain)
 P0=FULL._initial_covariance(src,path); Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,path,h); Ppre,bS0=ENTRY._zero_residual_S_covariance(Pp,src)
 if bS0!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": failures.append("sample-0 S prerequisite did not use fixed pivots")
 qpre=float(first["post_prediction_full_cayley_norm_upper"]); rows=[]; max_ds=max_rs=max_q=0.0; fixed=fallback=0; first_bad=None
 for di,d0 in enumerate(SUB.parts(0.0,dmax,d_pieces)):
  try:
   _Pj,Pa,_K0,_Rx,b0,dx0=TC.first_posterior(Ppre,H0,Racc,rho0,d0); fixed += b0=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"; fallback += b0!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"
   x0=[FULL.I(0.0) for _ in range(FULL.N)]
   for i in range(3,FULL.N): x0[i]=dx0[i]
   q0=PREFIX2._post_correction_q_upper(qpre,d0.hi); Pm=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pa),matrix_transpose(F)),Q)); xm=BASE._predict_state(x0,F); q1=RG._q_after_first_prediction(q0,domain,h); P1=BASE._transform_covariance(Pm,Rstep); x1=BASE._transform_state(xm,Rstep)
   rS=[-x1[12+i] for i in range(3)]; rSn=FULL._norm_upper(rS); Scell=FULL._measurement_cell(P1,FULL._H_S(),FULL._R_S(src),rS); fixed += Scell["inverse_backend"]=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"; fallback += Scell["inverse_backend"]!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"
   ds=FULL._norm_upper(Scell["dx"][0:3]); qS=PREFIX2._post_correction_q_upper(q1,ds); finite=math.isfinite(ds) and math.isfinite(qS) and qS<8.0
   row={"d_cell":di,"first_correction_rad":d0.as_list(),"sample1_S_residual_norm_upper_m_s":rSn,"sample1_S_inverse_backend":Scell["inverse_backend"],"sample1_S_attitude_correction_norm_upper_rad":ds,"sample1_q_before_S_upper":q1,"sample1_q_after_S_upper":qS,"finite_q8":finite}; rows.append(row); max_ds=max(max_ds,ds); max_rs=max(max_rs,rSn); max_q=max(max_q,qS if math.isfinite(qS) else math.inf)
   if not finite and first_bad is None:first_bad=row
  except Exception as exc:
   row={"d_cell":di,"exception":f"{type(exc).__name__}: {exc}"}; rows.append(row)
   if first_bad is None:first_bad=row
 ok=bool(rows) and first_bad is None and not failures
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_STRUCTURED_CONSECUTIVE_S_PREFIX_DIAGNOSTIC_V2","source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,"sample0_due_implies_sample1_due_for_refined_source_cell":phase1==["due"],"exact_first_Joseph_reset_gauge_child_used":True,"one_step_process_prediction_included":True,"sample1_body_rotation_gauge_included":True,"actual_estimator_S_residual_used":True,"sample1_S_shipping_gain_and_Joseph_reset_used":True,"rigorous_S_ge_R_spectral_inverse_is_admissible":True,"due_S_solver_identity_branch_has_zero_correction":True,"second_accelerometer_evaluated_here":False,"complete_sample1_branch_closed_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"evaluated_d_cells":len(rows),"fixed_pivot_inverse_count":fixed,"spectral_fallback_inverse_count":fallback,"max_sample1_S_residual_norm_upper_m_s":max_rs,"max_sample1_S_attitude_correction_norm_upper_rad":max_ds,"max_sample1_q_after_S_upper":max_q,"first_unclosed_cell":first_bad,"P5_SAMPLE1_STRUCTURED_S_PREFIX_DIAGNOSTIC":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"RECOMPUTE_SECOND_ACCELEROMETER_DIRECTIONAL_CHANNEL_AFTER_STRUCTURED_RESET_PROCESS_S_PREFIX","failures":failures,"rows":rows}

def validate(d):
 f=list(d.get("failures",[]))
 if d.get("schema")!=SCHEMA:f.append("schema mismatch")
 for k in ("source_generated_not_trajectory_fit","sample0_due_implies_sample1_due_for_refined_source_cell","exact_first_Joseph_reset_gauge_child_used","one_step_process_prediction_included","sample1_body_rotation_gauge_included","actual_estimator_S_residual_used","sample1_S_shipping_gain_and_Joseph_reset_used","rigorous_S_ge_R_spectral_inverse_is_admissible","due_S_solver_identity_branch_has_zero_correction"):
  if d.get(k) is not True:f.append(k)
 for k in ("source_replay_used","filter_changed","second_accelerometer_evaluated_here","complete_sample1_branch_closed_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 n=int(d.get("evaluated_d_cells",0));
 if n<=0:f.append("no d cells")
 if int(d.get("fixed_pivot_inverse_count",0))+int(d.get("spectral_fallback_inverse_count",0)) != 2*n:f.append("inverse count mismatch")
 for k in ("max_sample1_S_attitude_correction_norm_upper_rad","max_sample1_q_after_S_upper"):
  if not math.isfinite(float(d.get(k,math.nan))):f.append(f"nonfinite {k}")
 if d.get("P5_SAMPLE1_STRUCTURED_S_PREFIX_DIAGNOSTIC")=="PASS" and d.get("first_unclosed_cell") is not None:f.append("PASS retains unclosed cell")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=4); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--d-pieces",type=int,default=32); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.d_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_STRUCTURED_S_PREFIX_DIAGNOSTIC"],"cells":d["evaluated_d_cells"],"fixed":d["fixed_pivot_inverse_count"],"fallback":d["spectral_fallback_inverse_count"],"max_rS":d["max_sample1_S_residual_norm_upper_m_s"],"max_dS":d["max_sample1_S_attitude_correction_norm_upper_rad"],"max_qS":d["max_sample1_q_after_S_upper"],"first_unclosed":d["first_unclosed_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
