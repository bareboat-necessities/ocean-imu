#!/usr/bin/env python3
"""Dependency-preserving V2 of the reset-sensitive sample-1 scalar channel.

V1 formed the scalar innovation as m^2 T + 2 m C + B + r.  Although exact,
its natural interval extension loses the covariance dependency when C<0 and
can drive the lower bound back to r.  The same 2x2 covariance has determinant
Delta=T B-C^2>=0 and the exact identity

    T (S-r) = (m T+C)^2 + Delta.

For the ideal structured first tangent/aw posterior

    det(P+) = t p r / D,  D=g^2 t+p+r.

After reset+proof-gauge yaw mixing, with theta'=l theta+h yaw,

    Delta_r = l^2 t p r/D + h^2 Y b.

After one OU prediction, T1=T_r+q_theta, C1=alpha C_r,
B1=alpha^2 b+q_aw, so

    Delta_1 = alpha^2 Delta_r + T_r q_aw
              + q_theta (alpha^2 b+q_aw).

All terms are nonnegative.  Setting N=m T1+C1 and A=r T1+Delta_1 gives

    K_theta = T1 N / (A+N^2).

This producer bounds that ratio directly, including the exact one-dimensional
maximum in |N|, rather than independently dividing an interval numerator by a
cancellation-damaged innovation interval.

The modeled core includes the dominant reset/yaw mixing and one-step OU/process
variances.  As in V1, source attitude-remainder cross terms, body-rate mixing,
non-axial sample-1 force, and the already certified 1.8e-12-rad sample-1 S
attitude correction remain explicit perturbation obligations.  No complete
sample-1/P5 promotion is made.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_p5_sample1_reset_perp_scalar_channel as V1
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN; SCHEMA=2; RANGE=9.0

def _abs_range(x:Interval):
 if x.lo<=0.0<=x.hi:return 0.0,max(abs(x.lo),abs(x.hi))
 return min(abs(x.lo),abs(x.hi)),max(abs(x.lo),abs(x.hi))

def _ratio_upper(T:Interval,A:Interval,N:Interval):
 """Bound T*|N|/(A+N^2) from T<=Thi, A>=Alo>0."""
 if not (T.hi>0.0 and A.lo>0.0): raise RuntimeError("positive T/A floor required")
 nlo,nhi=_abs_range(N); candidates=[]
 def fup(n):
  if n<=0.0:return 0.0
  den=FULL.down(A.lo+FULL.down(n*n))
  if not den>0.0:raise RuntimeError("ratio denominator lost positivity")
  return FULL.up(n/den)
 candidates.extend((fup(nlo),fup(nhi)))
 root=math.sqrt(A.lo)
 # If the exact stationary point might lie in the admissible absolute range,
 # include its closed-form maximum 1/(2 sqrt(Alo)).
 if nlo<=root<=nhi:
  rlo=FULL.down(root)
  candidates.append(FULL.up(1.0/FULL.down(2.0*rlo)))
 return FULL.up(T.hi*max(candidates)),[nlo,nhi]

def build(domain_path=DEFAULT_DOMAIN,source_pieces=4,source_cell_index=0,p_pieces=32,d_pieces=32,axial_pieces=32):
 FULL3._install_backend(); path=Path(domain_path).resolve(); dom=json.loads(path.read_text()); first=FIRST.build(path,source_pieces=source_pieces); vec=VECTOR.build(); failures=[]
 failures += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase!="due":failures.append("expected due source witness")
 fr=first["source_cells"][source_cell_index]; p_all=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"])); rho0=float(fr["combined_useful_residual_norm_upper_mps2"]); dmax=float(fr["correction_norm_upper_rad"])
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tilt,yaw,eps=RG._attitude_covariance_epsilon(path,h); t=Interval.outward_bounds(tilt,FULL.up(tilt+eps)); Y=Interval.outward_bounds(yaw,FULL.up(yaw+eps)); Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); r=Racc[0][0]
 F,Q,_=FULL._transition_and_Q(src,dom); alpha=F[15][15]; qaw=Q[15][15]; qth=Interval(0.0,FULL.up(eps)); rho1=FULL.up(float(dom["normal_live"]["specific_force_norm_upper_mps2"])+float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"]))
 rows=[]; bad=None; maxk=maxd=maxT=0.0; minDelta=math.inf; minA=math.inf; pcells=SUB.parts(p_all.lo,p_all.hi,p_pieces); dcells=SUB.parts(0,dmax,d_pieces)
 for pi,p in enumerate(pcells):
  D=FULL.I(g*g)*t+p+r; a=t*(p+r)/D; c=-(FULL.I(g)*t*p/D); b=p*(FULL.I(g*g)*t+r)/D; det_first=t*p*r/D
  kg=Interval(FULL.down(p.lo/(p.lo+r.hi)),FULL.up(p.hi/(p.hi+r.lo))); azmax=FULL.up(kg.hi*rho0); azcells=SUB.parts(-azmax,azmax,axial_pieces)
  for di,d in enumerate(dcells):
   l,hc=V1._L_coeffs(d); Tr=l.square()*a+hc.square()*Y; Cr=l*c; det_r=l.square()*det_first+hc.square()*Y*b
   T1=Tr+qth; C1=alpha*Cr; B1=alpha.square()*b+qaw; det1=alpha.square()*det_r+Tr*qaw+qth*(alpha.square()*b+qaw)
   if not (T1.lo>0.0 and det1.lo>=0.0): failures.append(f"lost positive covariance structure p={pi} d={di}"); continue
   A=r*T1+det1
   if not A.lo>0.0: failures.append(f"lost A floor p={pi} d={di}"); continue
   maxT=max(maxT,T1.hi); minDelta=min(minDelta,det1.lo); minA=min(minA,A.lo)
   for ai,az in enumerate(azcells):
    m=FULL.I(g)+alpha*az; N=m*T1+C1; k,nabs=_ratio_upper(T1,A,N); corr=FULL.up(k*rho1); closed=math.isfinite(corr) and corr<RANGE; maxk=max(maxk,k); maxd=max(maxd,corr)
    # Positive innovation enclosure is reported for diagnosis, but K is bounded
    # from the equivalent ratio directly.
    Spos=r+(N.square()+det1)/T1
    row={"p_cell":pi,"d_cell":di,"axial_cell":ai,"P_aw_variance":p.as_list(),"first_correction_rad":d.as_list(),"first_axial_aw_correction_mps2":az.as_list(),"sample1_aligned_signed_force_mps2":m.as_list(),"T1":T1.as_list(),"C1":C1.as_list(),"B1":B1.as_list(),"determinant_lower":det1.lo,"A_lower":A.lo,"N_interval":N.as_list(),"N_abs_range":nabs,"positive_innovation_interval":Spos.as_list(),"Ktheta_abs_upper":k,"tangent_residual_norm_upper_mps2":rho1,"correction_norm_upper_rad":corr,"inside_9rad_range":closed}; rows.append(row)
    if not closed and bad is None:bad=row
 ok=bool(rows) and bad is None and not failures
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_POSITIVE_DETERMINANT_V2","source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,"first_scalar_Joseph_posterior_used":True,"exact_reset_plus_body_gauge_yaw_mixing_used":True,"first_posterior_determinant_identity_used":True,"reset_determinant_identity_used":True,"one_step_process_determinant_identity_used":True,"positive_innovation_identity_used":True,"direct_scalar_ratio_maximization_used":True,"dominant_reset_yaw_and_process_core_included":True,"source_attitude_remainder_cross_terms_included":False,"sample1_body_rate_rotation_perturbation_included":False,"sample1_tangent_force_perturbation_included":False,"sample1_S_attitude_correction_included":False,"complete_sample1_branch_closed_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"validated_deployed_quaternion_range_rad":RANGE,"evaluated_joint_cells":len(rows),"tangent_residual_norm_upper_mps2":rho1,"minimum_predicted_determinant_lower":minDelta,"minimum_A_lower":minA,"max_predicted_tangent_variance":maxT,"max_Ktheta_abs_upper":maxk,"max_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":bad,"P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V2":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"ADD_REMAINDER_CROSS_BODY_RATE_TANGENT_FORCE_AND_TINY_S_PERTURBATIONS" if ok else "SUBDIVIDE_SOURCE_PARAMETERS_OR_TIGHTEN_REMAINDER_SEPARATION","failures":failures,"rows":rows}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","first_scalar_Joseph_posterior_used","exact_reset_plus_body_gauge_yaw_mixing_used","first_posterior_determinant_identity_used","reset_determinant_identity_used","one_step_process_determinant_identity_used","positive_innovation_identity_used","direct_scalar_ratio_maximization_used","dominant_reset_yaw_and_process_core_included"):
  if d.get(k) is not True:f.append(k)
 for k in ("source_replay_used","filter_changed","source_attitude_remainder_cross_terms_included","sample1_body_rate_rotation_perturbation_included","sample1_tangent_force_perturbation_included","sample1_S_attitude_correction_included","complete_sample1_branch_closed_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 if int(d.get("evaluated_joint_cells",0))<=0:f.append("no cells")
 if not float(d.get("minimum_A_lower",0.0))>0.0:f.append("nonpositive A lower")
 st=d.get("P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V2"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=4); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--p-pieces",type=int,default=32); a.add_argument("--d-pieces",type=int,default=32); a.add_argument("--axial-pieces",type=int,default=32); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.p_pieces,x.d_pieces,x.axial_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V2"],"cells":d["evaluated_joint_cells"],"rho":d["tangent_residual_norm_upper_mps2"],"min_det":d["minimum_predicted_determinant_lower"],"min_A":d["minimum_A_lower"],"max_T":d["max_predicted_tangent_variance"],"max_K":d["max_Ktheta_abs_upper"],"max_d":d["max_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
