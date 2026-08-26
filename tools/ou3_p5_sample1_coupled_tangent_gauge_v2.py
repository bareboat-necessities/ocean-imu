#!/usr/bin/env python3
"""Gauge-consistent sample-1 coupled-tangent/axial P5 witness.

The established first accelerometer gauge maps the predicted specific force to
+g e3; ENTRY._canonical_first_H(g) is exactly -[+g e3]x.  A follow-on diagnostic
must therefore reconstruct the no-correction sample-1 force from +g e3, not
-g e3.  This producer repeats the schema-3 coupled-tangent/axial witness with
that sign fixed and records the gauge contract explicitly.

No filter or theorem bound is changed.  This remains source-cell-0 and the
sample-1 S-identity subbranch only.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_entry as ENTRY
import ou3_p5_sample1_prefix_v2 as PREFIX2
import ou3_p5_sample1_rotation_gauge_refinement as RBASE
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_p5_sample1_transported_crosscov_refinement as TC
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=TC.DEFAULT_DOMAIN; SCHEMA=4; N=FULL.N; LIMIT=6.0

def op(M): return RG._op2_upper(M)

def build(domain_path=DEFAULT_DOMAIN,source_pieces=4,source_cell_index=0,delta_pieces=16,axial_pieces=16):
 FULL3._install_backend(); p=Path(domain_path).resolve(); dom=json.loads(p.read_text()); first=FIRST.build(p,source_pieces=source_pieces); vec=VECTOR.build(); failures=[]
 failures += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase0!="due": failures.append("expected due source witness")
 fr=first["source_cells"][source_cell_index]; dmax=float(fr["correction_norm_upper_rad"]); rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tf=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
 Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); rvar=Racc[0][0]; H0=ENTRY._canonical_first_H(g)
 F,Q,Rstep=FULL._transition_and_Q(src,dom); alpha=F[15][15]; B=[[F[i][3+j] for j in range(3)] for i in range(3)]; Q1=RBASE._transform_covariance(Q,Rstep)
 P0=FULL._initial_covariance(src,p); Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,p,h); Ppre,_=ENTRY._zero_residual_S_covariance(Pp,src)
 ep=FULL._predict_error(FULL._initial_error(dom),F); x0=[FULL.I(0) for _ in range(N)]; qpre=float(first["post_prediction_full_cayley_norm_upper"])
 paw=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"])); tilt,_,eps=RG._attitude_covariance_epsilon(p,h); beta=Interval(FULL.down(paw.lo/(g*tilt)),FULL.up(paw.hi/(g*tilt)))
 lam=FULL.down(paw.lo+rvar.lo+g*g*tilt); ds=FULL.up(g*g*eps)
 if not lam>0: failures.append("tangent innovation floor not positive")
 dinv=FULL.up(ds/(lam*lam)); dkaw=FULL.up(paw.hi*dinv); dkth=FULL.up(FULL.up(eps*g/lam)+FULL.up(g*tilt*dinv)); defect=FULL.up(dkaw+FULL.up(beta.hi*dkth)); rem0=FULL.up(defect*rho0)
 axial_gain=FULL.up(paw.hi/(paw.hi+rvar.lo)); axial0=FULL.up(axial_gain*rho0); axial=SUB.parts(-axial0,axial0,axial_pieces); deltas=SUB.parts(0,dmax,delta_pieces)
 rows=[]; first_bad=None; fixed=fallback=0; maxe=maxpe=maxc=maxd=0.0
 for di,d0 in enumerate(deltas):
  Pj,_,K0,Rx,b0,_=TC.first_posterior(Ppre,H0,Racc,rho0,d0); _,ea,_,_=SUB.first_child(Ppre,ep,x0,H0,Racc,rho0,d0); q0=PREFIX2._post_correction_q_upper(qpre,d0.hi); em=RBASE._predict_state(ea,F); q1=RG._q_after_first_prediction(q0,dom,h); e1=RBASE._transform_state(em,Rstep); rho=FULL.up(FULL.up(RBASE._rot_diff(q1)*tf)+FULL.up(RBASE._norm(e1,FULL.AW)+ba))
  Gf=FULL._reset_matrix([d0,FULL.I(0),FULL.I(0)]); G=[[Gf[i][j] for j in range(3)] for i in range(3)]; At=matrix_mul(Rstep,G); AI=[[alpha if i==j else FULL.I(0) for j in range(3)] for i in range(3)]; Aaw=matrix_mul(matrix_mul(Rstep,AI),Rx); C0=matrix_mul(K0,Racc); rb=Interval(-rem0,rem0)
  for ai,az in enumerate(axial):
   aw0=[rb,-(beta*d0)+rb,az]; w=[alpha*aw0[i] for i in range(3)]
   # Gauge contract: sample-0 H0 is for +g e3, so the unchanged sample-1
   # reference force is +g e3 in the same canonical proof coordinates.
   w[2]=w[2]+FULL.I(g)
   fc=FULL._mat_vec(matrix_mul(Rstep,Rx),w); H1=SUB.Hforce(fc)
   Hth=[[H1[i][j] for j in range(3)] for i in range(3)]; H0th=[[H0[i][j] for j in range(3)] for i in range(3)]; Eth0=matrix_mul(Hth,At); Eth=[[Eth0[i][j]-H0th[i][j] for j in range(3)] for i in range(3)]; Ebg=matrix_mul(Hth,B); Eaw=[[Aaw[i][j]-FULL.I(1 if i==j else 0) for j in range(3)] for i in range(3)]; E=FULL._zero(3,N)
   for i in range(3):
    for j in range(3): E[i][j]=Eth[i][j]; E[i][3+j]=Ebg[i][j]; E[i][15+j]=Eaw[i][j]
   HA=[[H0[i][j]+E[i][j] for j in range(N)] for i in range(3)]; PE=matrix_mul(Pj,matrix_transpose(E)); M=matrix_add(C0,PE); Mth=[r[:] for r in M[:3]]; Mbg=[r[:] for r in M[3:6]]; CA=matrix_add(matrix_mul(At,Mth),matrix_mul(B,Mbg)); QH=matrix_mul(Q1,matrix_transpose(H1)); C=matrix_add(CA,[r[:] for r in QH[:3]])
   Ssig=matrix_mul(matrix_mul(HA,Pj),matrix_transpose(HA)); Sq=matrix_mul(matrix_mul(H1,Q1),matrix_transpose(H1)); S1=FULL.matrix_symmetric_hull(matrix_add(matrix_add(Ssig,Sq),Racc)); Sinv,backend=FULL._spd_inverse_enclosure(S1,Racc); Kth=matrix_mul(C,Sinv); kn=op(Kth); corr=FULL.up(kn*rho); qa=PREFIX2._post_correction_q_upper(q1,corr)
   fixed += backend=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"; fallback += backend!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"; en=op(Eth); pen=op([r[:] for r in PE[:6]]); cn=op(C); maxe=max(maxe,en); maxpe=max(maxpe,pen); maxc=max(maxc,cn); maxd=max(maxd,corr); closed=backend=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN" and corr<=LIMIT and math.isfinite(qa) and qa<8
   row={"delta_cell":di,"axial_cell":ai,"delta_rad":d0.as_list(),"axial_aw_correction_mps2":az.as_list(),"force_components_mps2":[x.as_list() for x in fc],"E_theta_norm_upper":en,"Pj_Et_first6_norm_upper":pen,"actual_Ctheta_norm_upper":cn,"inverse_backend":backend,"Ktheta_norm_upper":kn,"residual_norm_upper_mps2":rho,"correction_norm_upper_rad":corr,"post_q_upper":qa,"closed":closed}; rows.append(row)
   if not closed and first_bad is None:first_bad=row
 passed=bool(rows) and first_bad is None and not failures
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_CANONICAL_PLUS_GRAVITY_GAUGE_COUPLED_AXIAL_WITNESS","source_generated_not_trajectory_fit":True,"filter_changed":False,"sample0_canonical_specific_force_is_plus_g_e3":True,"sample1_zero_correction_reference_specific_force_is_plus_g_e3":True,"gravity_sign_consistent_with_canonical_H0":True,"same_first_residual_attitude_aw_coupling_used":True,"structured_tangent_gain_cancellation_used_before_interval_inversion":True,"axial_first_aw_correction_subdivided":True,"sample1_S_identity_subbranch_only":True,"beta_interval_mps2_per_rad":beta.as_list(),"first_tangent_relation_remainder_norm_upper_mps2":rem0,"first_axial_aw_correction_abs_upper_mps2":axial0,"evaluated_joint_cells":len(rows),"fixed_pivot_inverse_count":fixed,"spectral_fallback_inverse_count":fallback,"max_E_theta_norm_upper":maxe,"max_Pj_Et_first6_norm_upper":maxpe,"max_actual_Ctheta_norm_upper":maxc,"max_sample1_acc_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":first_bad,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"complete_sample1_branch_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"P5_SAMPLE1_COUPLED_TANGENT_GAUGE_V2":"PASS" if passed else "NOT_ESTABLISHED","next_obligation":"ADD_SAMPLE1_S_PERTURBATION_AND_LIFT_SOURCE_FAMILY" if passed else "IF_PIVOTS_STILL_LOST_DERIVE_REPEATED_SCALAR_TANGENT_CHANNEL","failures":failures,"rows":rows}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","sample0_canonical_specific_force_is_plus_g_e3","sample1_zero_correction_reference_specific_force_is_plus_g_e3","gravity_sign_consistent_with_canonical_H0","same_first_residual_attitude_aw_coupling_used","structured_tangent_gain_cancellation_used_before_interval_inversion","axial_first_aw_correction_subdivided","sample1_S_identity_subbranch_only"):
  if d.get(k) is not True:f.append(k)
 for k in ("filter_changed","deployed_correction_limit_increased","complete_sample1_branch_refined_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 if d.get("deployed_correction_limit_rad")!=6.0:f.append("correction limit changed")
 if int(d.get("fixed_pivot_inverse_count",0))+int(d.get("spectral_fallback_inverse_count",0))!=int(d.get("evaluated_joint_cells",0)):f.append("inverse count mismatch")
 st=d.get("P5_SAMPLE1_COUPLED_TANGENT_GAUGE_V2"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=4); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--delta-pieces",type=int,default=16); a.add_argument("--axial-pieces",type=int,default=16); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.delta_pieces,x.axial_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_COUPLED_TANGENT_GAUGE_V2"],"cells":d["evaluated_joint_cells"],"fixed":d["fixed_pivot_inverse_count"],"fallback":d["spectral_fallback_inverse_count"],"rem0":d["first_tangent_relation_remainder_norm_upper_mps2"],"max_Etheta":d["max_E_theta_norm_upper"],"max_PjE":d["max_Pj_Et_first6_norm_upper"],"max_C":d["max_actual_Ctheta_norm_upper"],"max_d":d["max_sample1_acc_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
