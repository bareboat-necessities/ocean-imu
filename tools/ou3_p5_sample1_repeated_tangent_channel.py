#!/usr/bin/env python3
"""Dependency-preserving repeated tangent channel witness for P5 sample 1.

For one gravity-tangent axis before the first accelerometer,
  P0=diag(t,p), H0=[g,1], R=r.
After the first accepted measurement the exact structured posterior is
  a=t(p+r)/D, c=-gtp/D, b=p(g^2 t+r)/D, D=g^2t+p+r.
For a later aligned signed force magnitude m, H1=[m,1], hence
  C_theta=t[p(m-g)+r m]/D
and
  S1=r+[t p (m-g)^2+t r m^2+p r]/D.
This positive form preserves the repeated-measurement dependencies that the
3x3 interval inverse loses.

The correction considered here is one tangent attitude channel.  In the aligned
core its predicted tangent specific-force component is zero, so the raw tangent
accelerometer innovation is bounded directly by the declared measured-force
norm upper plus the H-mode accelerometer-bias error bound.  It is incorrect to
charge this scalar row with twice the full 3-D force norm plus an independent
latent-a_w error: the a_w coordinate is already part of H1 and the scalar
posterior above.  No temporal force-slew assumption is introduced.

Reset/process/tangent-force perturbations remain deliberately omitted; PASS is
only a repeated scalar-core range statement and cannot promote P5/N_H_words.
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
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN; SCHEMA=2; LIMIT=6.0

def parts(x:Interval,n:int): return SUB.parts(x.lo,x.hi,n)

def build(domain_path=DEFAULT_DOMAIN,source_pieces=4,source_cell_index=0,p_pieces=16,axial_pieces=16):
 FULL3._install_backend(); path=Path(domain_path).resolve(); dom=json.loads(path.read_text()); first=FIRST.build(path,source_pieces=source_pieces); vec=VECTOR.build(); failures=[]
 failures += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase!="due": failures.append("expected due source witness")
 fr=first["source_cells"][source_cell_index]; p0=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"])); rho0=float(fr["combined_useful_residual_norm_upper_mps2"])
 live_force=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"]); rho1=FULL.up(live_force+ba)
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); t,_,_=RG._attitude_covariance_epsilon(path,h); r=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]; F,_,_=FULL._transition_and_Q(src,dom); alpha=F[15][15]
 rows=[]; bad=None; maxk=maxd=0.0; minS=math.inf
 for pi,p in enumerate(parts(p0,p_pieces)):
  llo=FULL.down(p.lo/(p.lo+r.hi)); lhi=FULL.up(p.hi/(p.hi+r.lo)); l=Interval(llo,lhi); azmax=FULL.up(l.hi*rho0)
  for ai,az in enumerate(SUB.parts(-azmax,azmax,axial_pieces)):
   m=FULL.I(g)+alpha*az; D=FULL.I(g*g*t)+p+r
   C=FULL.I(t)*(p*(m-FULL.I(g))+r*m)/D
   S=r+(FULL.I(t)*p*(m-FULL.I(g)).square()+FULL.I(t)*r*m.square()+p*r)/D
   if S.lo<=0: failures.append(f"nonpositive repeated scalar innovation floor p={pi} axial={ai}"); continue
   k=FULL.up(C.abs_upper()/S.lo); corr=FULL.up(k*rho1); safe=math.isfinite(corr) and corr<=LIMIT; maxk=max(maxk,k); maxd=max(maxd,corr); minS=min(minS,S.lo)
   row={"p_cell":pi,"axial_cell":ai,"P_aw_variance":p.as_list(),"first_axial_gain":l.as_list(),"first_axial_correction_mps2":az.as_list(),"sample1_aligned_signed_force_mps2":m.as_list(),"Ctheta_interval":C.as_list(),"innovation_variance_interval":S.as_list(),"Ktheta_abs_upper":k,"tangent_measured_force_component_abs_upper_mps2":live_force,"accelerometer_bias_error_component_abs_upper_mps2":ba,"tangent_residual_norm_upper_mps2":rho1,"correction_norm_upper_rad":corr,"inside_deployed_range":safe}; rows.append(row)
   if not safe and bad is None: bad=row
 passed=bool(rows) and bad is None and not failures
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_REPEATED_SCALAR_TANGENT_CORE_ROWWISE_SOURCE_WITNESS","source_generated_not_trajectory_fit":True,"filter_changed":False,"exact_first_structured_posterior_used":True,"dependency_preserving_positive_repeated_innovation_formula_used":True,"signed_aligned_sample1_force_used":True,"rowwise_tangent_residual_bound_used":True,"temporal_force_slew_assumed":False,"latent_aw_error_double_counted_in_raw_residual":False,"reset_process_and_tangent_force_perturbations_included":False,"complete_sample1_branch_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"source_pieces":source_pieces,"p_cell_count":p_pieces,"axial_cell_count":axial_pieces,"evaluated_joint_cells":len(rows),"tangent_measured_force_component_abs_upper_mps2":live_force,"accelerometer_bias_error_component_abs_upper_mps2":ba,"sample1_tangent_residual_norm_upper_mps2":rho1,"minimum_scalar_innovation_variance_lower":minS,"max_scalar_Ktheta_abs_upper":maxk,"max_scalar_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":bad,"P5_SAMPLE1_REPEATED_TANGENT_CORE_WITNESS":"PASS" if passed else "NOT_ESTABLISHED","next_obligation":"ADD_RESET_PROCESS_AND_TANGENT_FORCE_PERTURBATIONS_TO_SCALAR_CORE" if passed else "EXTEND_DEPLOYED_QUATERNION_RANGE_OR_TIGHTEN_SOURCE_CORRELATED_TANGENT_RESIDUAL","failures":failures,"rows":rows}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","exact_first_structured_posterior_used","dependency_preserving_positive_repeated_innovation_formula_used","signed_aligned_sample1_force_used","rowwise_tangent_residual_bound_used"):
  if d.get(k) is not True:f.append(k)
 for k in ("filter_changed","temporal_force_slew_assumed","latent_aw_error_double_counted_in_raw_residual","reset_process_and_tangent_force_perturbations_included","complete_sample1_branch_refined_here","whole_word_promoted_here","N_H_words_set_here","deployed_correction_limit_increased"):
  if d.get(k) is not False:f.append(k)
 if d.get("deployed_correction_limit_rad")!=6.0:f.append("correction limit changed")
 if int(d.get("evaluated_joint_cells",0))<=0:f.append("no cells")
 if not (math.isfinite(float(d.get("minimum_scalar_innovation_variance_lower",math.nan))) and float(d.get("minimum_scalar_innovation_variance_lower",0))>0):f.append("invalid innovation floor")
 st=d.get("P5_SAMPLE1_REPEATED_TANGENT_CORE_WITNESS"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=4); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--p-pieces",type=int,default=16); a.add_argument("--axial-pieces",type=int,default=16); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.p_pieces,x.axial_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_REPEATED_TANGENT_CORE_WITNESS"],"cells":d["evaluated_joint_cells"],"rho_tangent":d["sample1_tangent_residual_norm_upper_mps2"],"min_S":d["minimum_scalar_innovation_variance_lower"],"max_K":d["max_scalar_Ktheta_abs_upper"],"max_d":d["max_scalar_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
