#!/usr/bin/env python3
"""Tighten the transported sample-1 PH' mismatch with the source force cone.

At sample 0, before the first accelerometer, P_aw=pI and P_theta,aw=0 while
H_theta=-[g e3]_x.  Hence S is exactly block diagonal between the gravity-axis
measurement channel and its tangent plane.  Since P_theta,t >= t I,
  ||K_aw,t|| <= p_hi/(p_lo+R_lo+g^2 t).
The axial gain may be much larger but remains parallel to gravity.  Reset/gauge
rotations preserve this tangent/axial decomposition and the next OU prediction
only multiplies a_w by scalar alpha.  Thus the sample-1 predicted force is
covered by an axial interval along the transported gravity axis plus a small
tangent ball, rather than an arbitrary 3-D component cube.
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
import ou3_p5_sample1_transported_crosscov_refinement as TC
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=TC.DEFAULT_DOMAIN; SCHEMA=1; N=FULL.N; LIMIT=6.0

def axial_parts(lo,hi,n): return SUB.parts(lo,hi,n)

def force_cone_cells(u,zint,tan,n):
 for z in axial_parts(zint.lo,zint.hi,n):
  yield tuple(z*u[i]+Interval(-FULL.up(tan),FULL.up(tan)) for i in range(3))

def block_contribution(P,D,cols):
 Db=[[D[i][j] for j in cols] for i in range(3)]; Pb=[[P[i][j] for j in cols] for i in range(3)]; return matrix_mul(Pb,matrix_transpose(Db))

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2,source_cell_index=0,delta_pieces=8,axial_pieces=8):
 FULL3._install_backend(); pth=Path(domain_path).resolve(); dom=json.loads(pth.read_text()); first=FIRST.build(pth,source_pieces=source_pieces); vec=VECTOR.build(); fails=[]
 fails += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase0!="due": fails.append("expected due witness")
 fr=first["source_cells"][source_cell_index]; dmax=float(fr["correction_norm_upper_rad"]); rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tf=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
 Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); rmin=min(Racc[i][i].lo for i in range(3)); H0=ENTRY._canonical_first_H(g)
 F,Q,Rstep=FULL._transition_and_Q(src,dom); alpha=F[15][15]; P0=FULL._initial_covariance(src,pth); Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,pth,h); Ppre,_=ENTRY._zero_residual_S_covariance(Pp,src); ep=FULL._predict_error(FULL._initial_error(dom),F); x0=[FULL.I(0) for _ in range(N)]; qpre=float(first["post_prediction_full_cayley_norm_upper"])
 tilt,_,_=RG._attitude_covariance_epsilon(pth,h); paw=Ppre[15][15]; den_t=FULL.down(paw.lo+rmin+g*g*tilt); den_z=FULL.down(paw.lo+rmin)
 if not (den_t>0 and den_z>0): fails.append("first aw gain denominator lost positivity")
 kaw_t=FULL.up(paw.hi/den_t); kaw_z=FULL.up(paw.hi/den_z); aw_t0=FULL.up(kaw_t*rho0); aw_z0=FULL.up(kaw_z*rho0); aw_t1=FULL.up(alpha.hi*aw_t0); aw_z1=FULL.up(alpha.hi*aw_z0); zint=Interval.outward_bounds(-g-aw_z1,-g+aw_z1)
 rows=[]; bad=None; maxref=maxth=maxbg=maxaw=maxmis=maxc=maxd=0.0
 e3=[FULL.I(0),FULL.I(0),FULL.I(1)]
 for di,d0 in enumerate(SUB.parts(0,dmax,delta_pieces)):
  Pj,Pa,K0,Rx,b0,_=TC.first_posterior(Ppre,H0,Racc,rho0,d0); _,ea,xa,_=SUB.first_child(Ppre,ep,x0,H0,Racc,rho0,d0); q0=PREFIX2._post_correction_q_upper(qpre,d0.hi)
  Pm=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pa),matrix_transpose(F)),Q)); em=BASE._predict_state(ea,F); xm=BASE._predict_state(xa,F); q1=RG._q_after_first_prediction(q0,dom,h); P1=BASE._transform_covariance(Pm,Rstep); e1=BASE._transform_state(em,Rstep); x1=BASE._transform_state(xm,Rstep)
  Href=TC.href_matrix(H0,F,Rstep,Rx,d0); C0=matrix_mul(K0,Racc); G=FULL._reset_matrix([d0,FULL.I(0),FULL.I(0)]); G3=[[G[i][j] for j in range(3)] for i in range(3)]; Cth=[r[:] for r in C0[:3]]; Cbg=[r[:] for r in C0[3:6]]; B=[[F[i][3+j] for j in range(3)] for i in range(3)]; Ctransport=matrix_add(matrix_mul(matrix_mul(Rstep,G3),Cth),matrix_mul(B,Cbg)); Q1=BASE._transform_covariance(Q,Rstep); QHref=matrix_mul(Q1,matrix_transpose(Href)); Cref=matrix_add(Ctransport,[r[:] for r in QHref[:3]]); refn=RG._op2_upper(Cref); maxref=max(maxref,refn)
  u=FULL._mat_vec(matrix_mul(Rstep,Rx),e3); eaw=BASE._norm(e1,FULL.AW); rho=FULL.up(FULL.up(BASE._rot_diff(q1)*tf)+FULL.up(eaw+ba))
  for zi,fc in enumerate(force_cone_cells(u,zint,aw_t1,axial_pieces)):
   H1=SUB.Hforce(fc); D=[[H1[i][j]-Href[i][j] for j in range(N)] for i in range(3)]; cth=block_contribution(P1,D,range(0,3)); cbg=block_contribution(P1,D,range(3,6)); caw=block_contribution(P1,D,range(15,18)); nth=RG._op2_upper(cth); nbg=RG._op2_upper(cbg); naw=RG._op2_upper(caw); mis=matrix_add(matrix_add(cth,cbg),caw); misn=RG._op2_upper(mis); C=matrix_add(Cref,mis); cn=RG._op2_upper(C); k=FULL.up(cn/rmin); dc=FULL.up(k*rho); qa=PREFIX2._post_correction_q_upper(q1,dc); maxth=max(maxth,nth); maxbg=max(maxbg,nbg); maxaw=max(maxaw,naw); maxmis=max(maxmis,misn); maxc=max(maxc,cn); maxd=max(maxd,dc); closed=dc<=LIMIT and math.isfinite(qa) and qa<8
   r={"delta_cell":di,"axial_cell":zi,"force_components":[x.as_list() for x in fc],"transported_reference_Ctheta_norm_upper":refn,"theta_mismatch_contribution_norm_upper":nth,"bg_mismatch_contribution_norm_upper":nbg,"aw_mismatch_contribution_norm_upper":naw,"total_mismatch_Ctheta_norm_upper":misn,"actual_Ctheta_norm_upper":cn,"Ktheta_norm_upper":k,"residual_norm_upper_mps2":rho,"correction_norm_upper_rad":dc,"post_q_upper":qa,"closed":closed}; rows.append(r)
   if not closed and bad is None: bad=r
 ok=bool(rows) and bad is None and not fails
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_FORCE_CONE_TRANSPORTED_CROSSCOV_WITNESS","source_generated_not_trajectory_fit":True,"filter_changed":False,"first_aw_tangent_gain_uses_gravity_information_denominator":True,"first_aw_axial_gain_kept_separate":True,"reset_and_prediction_preserve_force_cone_decomposition":True,"transported_first_posterior_crosscovariance_used":True,"sample1_S_identity_subbranch_only":True,"first_aw_tangent_gain_norm_upper":kaw_t,"first_aw_axial_gain_abs_upper":kaw_z,"sample0_aw_tangent_mean_norm_upper_mps2":aw_t0,"sample0_aw_axial_mean_abs_upper_mps2":aw_z0,"sample1_aw_tangent_mean_norm_upper_mps2":aw_t1,"sample1_aw_axial_mean_abs_upper_mps2":aw_z1,"evaluated_joint_cells":len(rows),"max_transported_reference_Ctheta_norm_upper":maxref,"max_theta_mismatch_contribution_norm_upper":maxth,"max_bg_mismatch_contribution_norm_upper":maxbg,"max_aw_mismatch_contribution_norm_upper":maxaw,"max_total_mismatch_Ctheta_norm_upper":maxmis,"max_actual_Ctheta_norm_upper":maxc,"max_sample1_acc_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":bad,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"complete_sample1_branch_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"P5_SAMPLE1_FORCE_CONE_CROSSCOV_WITNESS_REFINEMENT":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"ADD_SAMPLE1_S_PERTURBATION_AND_LIFT_SOURCE_FAMILY" if ok else "TIGHTEN_DOMINANT_REPORTED_MISMATCH_COMPONENT","failures":fails}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","first_aw_tangent_gain_uses_gravity_information_denominator","first_aw_axial_gain_kept_separate","reset_and_prediction_preserve_force_cone_decomposition","transported_first_posterior_crosscovariance_used","sample1_S_identity_subbranch_only"):
  if d.get(k) is not True:f.append(k)
 for k in ("filter_changed","deployed_correction_limit_increased","complete_sample1_branch_refined_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 st=d.get("P5_SAMPLE1_FORCE_CONE_CROSSCOV_WITNESS_REFINEMENT"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=2); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--delta-pieces",type=int,default=8); a.add_argument("--axial-pieces",type=int,default=8); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.delta_pieces,x.axial_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_FORCE_CONE_CROSSCOV_WITNESS_REFINEMENT"],"cells":d["evaluated_joint_cells"],"Kaw_t":d["first_aw_tangent_gain_norm_upper"],"aw_t1":d["sample1_aw_tangent_mean_norm_upper_mps2"],"max_ref":d["max_transported_reference_Ctheta_norm_upper"],"max_theta_mismatch":d["max_theta_mismatch_contribution_norm_upper"],"max_bg_mismatch":d["max_bg_mismatch_contribution_norm_upper"],"max_aw_mismatch":d["max_aw_mismatch_contribution_norm_upper"],"max_mismatch":d["max_total_mismatch_Ctheta_norm_upper"],"max_C":d["max_actual_Ctheta_norm_upper"],"max_d":d["max_sample1_acc_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
