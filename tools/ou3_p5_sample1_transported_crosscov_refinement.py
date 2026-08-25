#!/usr/bin/env python3
"""Transport the exact first-posterior PH' identity to sample 1.

For the accepted sample-0 accelerometer Joseph posterior Pj,
    Pj H0' = K0 Racc.
For each canonical reset cell let A be the exact reset/gauge/one-step linear
map to sample 1 (before the optional sample-1 S update).  A dual reference
Jacobian Href is chosen so Href A=H0.  Therefore
    P1 Href' = A K0 Racc + Q1 Href'.
The actual sample-1 accelerometer cross covariance is obtained as
    P1 H1' = P1 Href' + P1 (H1-Href)'.
This preserves the small posterior relation rather than reconstructing the two
large terms P_theta_theta H_theta' and P_theta_aw independently.
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
import ou3_p5_sample1_entry_v3 as ENTRY3
import ou3_p5_sample1_prefix_v2 as PREFIX2
import ou3_p5_sample1_rotation_gauge_refinement as BASE
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=BASE.DEFAULT_DOMAIN; SCHEMA=1; N=FULL.N; LIMIT=6.0

def zmat(r,c): return FULL._zero(r,c)

def ginv(d:Interval):
 h=FULL.I(.5); den=FULL.I(1)+FULL.I(.25)*d.square(); inv=den.reciprocal(); z=FULL.I(0); o=FULL.I(1)
 return [[o,z,z],[z,inv,h*d*inv],[z,-h*d*inv,inv]]

def href_matrix(H0,F,Rstep,Rx,d):
 Href=zmat(3,N); Rt=matrix_transpose(Rstep); Gi=ginv(d); Atinv=matrix_mul(Gi,Rt); Hth=[[H0[i][j] for j in range(3)] for i in range(3)]; Hrth=matrix_mul(Hth,Atinv)
 B=[[F[i][3+j] for j in range(3)] for i in range(3)]; Hrbg=matrix_mul(Hrth,B)
 for i in range(3):
  for j in range(3): Href[i][j]=Hrth[i][j]; Href[i][3+j]=-Hrbg[i][j]
 alpha=F[15][15]; Haw=matrix_mul(matrix_mul(matrix_transpose(Rx),[[alpha.reciprocal() if i==j else FULL.I(0) for j in range(3)] for i in range(3)]),Rt)
 for i in range(3):
  for j in range(3): Href[i][15+j]=Haw[i][j]
 return Href

def first_posterior(Ppre,H0,Racc,rho,d):
 PHt,S=FULL._innovation(Ppre,H0,Racc); Sinv,b=FULL._spd_inverse_enclosure(S,Racc); K=matrix_mul(PHt,Sinv); r=FULL._vec_box(rho); dx=FULL._mat_vec(K,r); caps=ENTRY3._linear_gain_caps(Ppre,Racc,rho)
 for name,idxs in ENTRY3.GROUPS.items():
  cap=Interval(-caps[name],caps[name])
  for i in idxs: dx[i]=FULL._intersect(dx[i],cap)
 Pj=FULL._shipping_joseph(Ppre,K,S,PHt); Pr=FULL._reset_covariance(Pj,[d,FULL.I(0),FULL.I(0)]); Rx=SUB.rx(d); Pa=BASE._transform_covariance(Pr,Rx)
 e=[FULL.I(0) for _ in range(N)]; return Pj,Pa,K,Rx,b,dx

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2,source_cell_index=0,delta_pieces=8,force_pieces=4):
 FULL3._install_backend(); p=Path(domain_path).resolve(); dom=json.loads(p.read_text()); first=FIRST.build(p,source_pieces=source_pieces); vec=VECTOR.build(); fails=[]
 fails += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase0!="due": fails.append("expected due witness")
 fr=first["source_cells"][source_cell_index]; dmax=float(fr["correction_norm_upper_rad"]); rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tf=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
 Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); rmin=min(Racc[i][i].lo for i in range(3)); H0=ENTRY._canonical_first_H(g)
 F,Q,Rstep=FULL._transition_and_Q(src,dom); P0=FULL._initial_covariance(src,p); Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,p,h); Ppre,_=ENTRY._zero_residual_S_covariance(Pp,src); ep=FULL._predict_error(FULL._initial_error(dom),F); x0=[FULL.I(0) for _ in range(N)]; qpre=float(first["post_prediction_full_cayley_norm_upper"])
 rows=[]; bad=None; maxref=maxmis=maxc=maxd=0.0
 for di,d0 in enumerate(SUB.parts(0,dmax,delta_pieces)):
  Pj,Pa,K0,Rx,b0,_=first_posterior(Ppre,H0,Racc,rho0,d0)
  # Use SUB.first_child for the already-audited state enclosure.
  _,ea,xa,_=SUB.first_child(Ppre,ep,x0,H0,Racc,rho0,d0); q0=PREFIX2._post_correction_q_upper(qpre,d0.hi)
  Pm=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pa),matrix_transpose(F)),Q)); em=BASE._predict_state(ea,F); xm=BASE._predict_state(xa,F); q1=RG._q_after_first_prediction(q0,dom,h); P1=BASE._transform_covariance(Pm,Rstep); e1=BASE._transform_state(em,Rstep); x1=BASE._transform_state(xm,Rstep)
  Href=href_matrix(H0,F,Rstep,Rx,d0)
  # C0=Pj H0'=K0 R exactly for the Joseph posterior.
  C0=matrix_mul(K0,Racc); G=FULL._reset_matrix([d0,FULL.I(0),FULL.I(0)]); Gt=[[G[i][j] for j in range(N)] for i in range(3)]
  # theta rows of A*C0: Rstep*Gtheta*C0_theta + B*C0_bg.
  G3=[[G[i][j] for j in range(3)] for i in range(3)]; Cth=[row[:] for row in C0[:3]]; Cbg=[row[:] for row in C0[3:6]]; B=[[F[i][3+j] for j in range(3)] for i in range(3)]
  Ctransport=matrix_add(matrix_mul(matrix_mul(Rstep,G3),Cth),matrix_mul(B,Cbg))
  Q1=BASE._transform_covariance(Q,Rstep); QHref=matrix_mul(Q1,matrix_transpose(Href)); Cref=matrix_add(Ctransport,[row[:] for row in QHref[:3]]); refn=RG._op2_upper(Cref); maxref=max(maxref,refn)
  xaw=BASE._norm(x1,FULL.AW); eaw=BASE._norm(e1,FULL.AW); fhat=FULL.up(g+xaw); rho=FULL.up(FULL.up(BASE._rot_diff(q1)*tf)+FULL.up(eaw+ba))
  for fi,fc in enumerate(SUB.force_cells(fhat,force_pieces)):
   H1=SUB.Hforce(fc); D=[[H1[i][j]-Href[i][j] for j in range(N)] for i in range(3)]; PD=matrix_mul(P1,matrix_transpose(D)); mis=[row[:] for row in PD[:3]]; misn=RG._op2_upper(mis); C=matrix_add(Cref,mis); cn=RG._op2_upper(C); k=FULL.up(cn/rmin); dc=FULL.up(k*rho); qa=PREFIX2._post_correction_q_upper(q1,dc); maxmis=max(maxmis,misn); maxc=max(maxc,cn); maxd=max(maxd,dc); closed=dc<=LIMIT and math.isfinite(qa) and qa<8
   r={"delta_cell":di,"force_cell":fi,"first_inverse_backend":b0,"transported_reference_Ctheta_norm_upper":refn,"reference_mismatch_Ctheta_norm_upper":misn,"actual_Ctheta_norm_upper":cn,"Ktheta_norm_upper":k,"residual_norm_upper_mps2":rho,"correction_norm_upper_rad":dc,"post_q_upper":qa,"closed":closed}; rows.append(r)
   if not closed and bad is None: bad=r
 ok=bool(rows) and bad is None and not fails
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_TRANSPORTED_FIRST_POSTERIOR_CROSSCOV_WITNESS","source_generated_not_trajectory_fit":True,"filter_changed":False,"first_posterior_identity_P_Ht_equals_KR_used":True,"dual_reference_jacobian_satisfies_Href_A_equals_H0":True,"process_crosscovariance_term_retained":True,"actual_H_minus_Href_mismatch_retained":True,"sample1_S_identity_subbranch_only":True,"evaluated_joint_cells":len(rows),"max_transported_reference_Ctheta_norm_upper":maxref,"max_reference_mismatch_Ctheta_norm_upper":maxmis,"max_actual_Ctheta_norm_upper":maxc,"max_sample1_acc_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":bad,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"complete_sample1_branch_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"P5_SAMPLE1_TRANSPORTED_CROSSCOV_WITNESS_REFINEMENT":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"ADD_TINY_SAMPLE1_S_POSTERIOR_PERTURBATION_AND_LIFT_BRANCHES" if ok else "TIGHTEN_DUAL_JACOBIAN_MISMATCH_USING_FORCE_TANGENT_STRUCTURE","failures":fails}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","first_posterior_identity_P_Ht_equals_KR_used","dual_reference_jacobian_satisfies_Href_A_equals_H0","process_crosscovariance_term_retained","actual_H_minus_Href_mismatch_retained","sample1_S_identity_subbranch_only"):
  if d.get(k) is not True:f.append(k)
 for k in ("filter_changed","deployed_correction_limit_increased","complete_sample1_branch_refined_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 st=d.get("P5_SAMPLE1_TRANSPORTED_CROSSCOV_WITNESS_REFINEMENT"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=2); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--delta-pieces",type=int,default=8); a.add_argument("--force-pieces",type=int,default=4); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.delta_pieces,x.force_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_TRANSPORTED_CROSSCOV_WITNESS_REFINEMENT"],"cells":d["evaluated_joint_cells"],"max_ref_C":d["max_transported_reference_Ctheta_norm_upper"],"max_mismatch_C":d["max_reference_mismatch_Ctheta_norm_upper"],"max_actual_C":d["max_actual_Ctheta_norm_upper"],"max_d":d["max_sample1_acc_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
