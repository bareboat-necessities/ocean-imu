#!/usr/bin/env python3
"""Preserve first-posterior cancellation by pushing H1 through the forward map.

For the accepted first accelerometer posterior Pj, Pj H0'=K0 R.  Let A map
that pre-reset error coordinate through the deployed reset, canonical linear
gauge, one prediction, and sample-1 body gauge.  Before the sample-1 S update,
  P1 = A Pj A' + Q1.
Hence for the actual sample-1 accelerometer Jacobian H1,
  P1 H1' = A [K0 R + Pj E'] + Q1 H1',
  E = H1 A - H0.
The reset and force-direction mismatch is therefore formed before multiplication
by the large posterior covariance.  This retains cancellations lost by the
dual-reference P1(H1-Href)' interval product.
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
import ou3_p5_sample1_force_cone_crosscov_refinement as FC
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=TC.DEFAULT_DOMAIN; SCHEMA=1; N=FULL.N; LIMIT=6.0

def op(M): return RG._op2_upper(M)

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2,source_cell_index=0,delta_pieces=8,axial_pieces=8):
 FULL3._install_backend(); pth=Path(domain_path).resolve(); dom=json.loads(pth.read_text()); first=FIRST.build(pth,source_pieces=source_pieces); vec=VECTOR.build(); fails=[]
 fails += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase0!="due": fails.append("expected due witness")
 fr=first["source_cells"][source_cell_index]; dmax=float(fr["correction_norm_upper_rad"]); rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tf=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
 Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); rmin=min(Racc[i][i].lo for i in range(3)); H0=ENTRY._canonical_first_H(g)
 F,Q,Rstep=FULL._transition_and_Q(src,dom); alpha=F[15][15]; B=[[F[i][3+j] for j in range(3)] for i in range(3)]; P0=FULL._initial_covariance(src,pth); Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,pth,h); Ppre,_=ENTRY._zero_residual_S_covariance(Pp,src); ep=FULL._predict_error(FULL._initial_error(dom),F); x0=[FULL.I(0) for _ in range(N)]; qpre=float(first["post_prediction_full_cayley_norm_upper"])
 tilt,_,_=RG._attitude_covariance_epsilon(pth,h); paw=Ppre[15][15]; kaw_t=FULL.up(paw.hi/FULL.down(paw.hi+rmin+g*g*tilt)); kaw_z=FULL.up(paw.hi/FULL.down(paw.hi+rmin)); aw_t1=FULL.up(alpha.hi*FULL.up(kaw_t*rho0)); aw_z1=FULL.up(alpha.hi*FULL.up(kaw_z*rho0)); zint=FULL.Interval.outward_bounds(-g-aw_z1,-g+aw_z1) if hasattr(FULL,'Interval') else None
 # FULL does not export Interval as an attribute on every import path.
 if zint is None:
  from ou3_interval import Interval
  zint=Interval.outward_bounds(-g-aw_z1,-g+aw_z1)
 e3=[FULL.I(0),FULL.I(0),FULL.I(1)]; Q1=BASE._transform_covariance(Q,Rstep)
 rows=[]; bad=None; maxEth=maxEbg=maxEaw=maxPE=maxC=maxD=0.0
 for di,d0 in enumerate(SUB.parts(0,dmax,delta_pieces)):
  Pj,Pa,K0,Rx,b0,_=TC.first_posterior(Ppre,H0,Racc,rho0,d0); _,ea,xa,_=SUB.first_child(Ppre,ep,x0,H0,Racc,rho0,d0); q0=PREFIX2._post_correction_q_upper(qpre,d0.hi); Pm=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pa),matrix_transpose(F)),Q)); em=BASE._predict_state(ea,F); q1=RG._q_after_first_prediction(q0,dom,h); e1=BASE._transform_state(em,Rstep)
  Gfull=FULL._reset_matrix([d0,FULL.I(0),FULL.I(0)]); G=[[Gfull[i][j] for j in range(3)] for i in range(3)]; At=matrix_mul(Rstep,G); Aaw=matrix_mul(matrix_mul(Rstep,[[alpha if i==j else FULL.I(0) for j in range(3)] for i in range(3)]),Rx); u=FULL._mat_vec(matrix_mul(Rstep,Rx),e3); eaw=BASE._norm(e1,FULL.AW); rho=FULL.up(FULL.up(BASE._rot_diff(q1)*tf)+FULL.up(eaw+ba)); C0=matrix_mul(K0,Racc)
  for zi,fc in enumerate(FC.force_cone_cells(u,zint,aw_t1,axial_pieces)):
   H1=SUB.Hforce(fc); Hth=[[H1[i][j] for j in range(3)] for i in range(3)]; H0th=[[H0[i][j] for j in range(3)] for i in range(3)]; Eth0=matrix_mul(Hth,At); Eth=[[Eth0[i][j]-H0th[i][j] for j in range(3)] for i in range(3)]; Ebg=matrix_mul(Hth,B); Eaw0=Aaw; Eaw=[[Eaw0[i][j]-FULL.I(1 if i==j else 0) for j in range(3)] for i in range(3)]; Ethn=op(Eth); Ebgn=op(Ebg); Eawn=op(Eaw); maxEth=max(maxEth,Ethn); maxEbg=max(maxEbg,Ebgn); maxEaw=max(maxEaw,Eawn)
   E=FULL._zero(3,N)
   for i in range(3):
    for j in range(3): E[i][j]=Eth[i][j]; E[i][3+j]=Ebg[i][j]; E[i][15+j]=Eaw[i][j]
   PE=matrix_mul(Pj,matrix_transpose(E)); PEn=op([r[:] for r in PE[:6]]); maxPE=max(maxPE,PEn); M=matrix_add(C0,PE); Mth=[r[:] for r in M[:3]]; Mbg=[r[:] for r in M[3:6]]; CA=matrix_add(matrix_mul(At,Mth),matrix_mul(B,Mbg)); QH=matrix_mul(Q1,matrix_transpose(H1)); C=matrix_add(CA,[r[:] for r in QH[:3]]); cn=op(C); k=FULL.up(cn/rmin); dc=FULL.up(k*rho); qa=PREFIX2._post_correction_q_upper(q1,dc); maxC=max(maxC,cn); maxD=max(maxD,dc); closed=dc<=LIMIT and math.isfinite(qa) and qa<8
   r={"delta_cell":di,"axial_cell":zi,"first_inverse_backend":b0,"E_theta_norm_upper":Ethn,"E_bg_norm_upper":Ebgn,"E_aw_norm_upper":Eawn,"Pj_Et_first6_norm_upper":PEn,"actual_Ctheta_norm_upper":cn,"Ktheta_norm_upper":k,"residual_norm_upper_mps2":rho,"correction_norm_upper_rad":dc,"post_q_upper":qa,"closed":closed}; rows.append(r)
   if not closed and bad is None: bad=r
 ok=bool(rows) and bad is None and not fails
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_FORWARD_MAP_CROSSCOV_WITNESS","source_generated_not_trajectory_fit":True,"filter_changed":False,"first_posterior_identity_P_Ht_equals_KR_used":True,"actual_H1_pushed_through_forward_reset_prediction_map":True,"mismatch_E_formed_before_covariance_multiplication":True,"source_reachable_force_cone_used":True,"sample1_S_identity_subbranch_only":True,"first_aw_tangent_gain_norm_upper":kaw_t,"first_aw_axial_gain_abs_upper":kaw_z,"sample1_aw_tangent_mean_norm_upper_mps2":aw_t1,"sample1_aw_axial_mean_abs_upper_mps2":aw_z1,"evaluated_joint_cells":len(rows),"max_E_theta_norm_upper":maxEth,"max_E_bg_norm_upper":maxEbg,"max_E_aw_norm_upper":maxEaw,"max_Pj_Et_first6_norm_upper":maxPE,"max_actual_Ctheta_norm_upper":maxC,"max_sample1_acc_correction_norm_upper_rad":maxD,"first_unclosed_joint_cell":bad,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"complete_sample1_branch_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"P5_SAMPLE1_FORWARD_CROSSCOV_WITNESS_REFINEMENT":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"ADD_SAMPLE1_S_PERTURBATION_AND_LIFT_SOURCE_FAMILY" if ok else "TIGHTEN_REPORTED_FORWARD_MISMATCH_COMPONENT","failures":fails}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","first_posterior_identity_P_Ht_equals_KR_used","actual_H1_pushed_through_forward_reset_prediction_map","mismatch_E_formed_before_covariance_multiplication","source_reachable_force_cone_used","sample1_S_identity_subbranch_only"):
  if d.get(k) is not True:f.append(k)
 for k in ("filter_changed","deployed_correction_limit_increased","complete_sample1_branch_refined_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 st=d.get("P5_SAMPLE1_FORWARD_CROSSCOV_WITNESS_REFINEMENT"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=2); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--delta-pieces",type=int,default=8); a.add_argument("--axial-pieces",type=int,default=8); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.delta_pieces,x.axial_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_FORWARD_CROSSCOV_WITNESS_REFINEMENT"],"cells":d["evaluated_joint_cells"],"max_Etheta":d["max_E_theta_norm_upper"],"max_Ebg":d["max_E_bg_norm_upper"],"max_Eaw":d["max_E_aw_norm_upper"],"max_PjE":d["max_Pj_Et_first6_norm_upper"],"max_C":d["max_actual_Ctheta_norm_upper"],"max_d":d["max_sample1_acc_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
