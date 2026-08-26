#!/usr/bin/env python3
"""Reset-sensitive perpendicular scalar accelerometer channel for P5 sample 1.

The first accelerometer is gravity-gauged with force +g e3.  Rotational
symmetry places its attitude correction on +x.  The second tangent channel
perpendicular to that correction is the one most exposed to yaw covariance.
For the shipping reset G(d)=I+0.5[d e1]x followed by the corrected-body proof
gauge R_x(d),

    [theta_y'; theta_z'] = L(d)[theta_y;theta_z],
    L = R_x(d) [[1,-d/2],[d/2,1]].

Before reset, one first tangent/aw scalar posterior from prior diag(t,p) and
H=[g,1] is

    a=t(p+r)/D, c=-g t p/D, b=p(g^2 t+r)/D,
    D=g^2 t+p+r.

Yaw is unobserved with variance Y.  Since aw_x is unchanged by R_x(d), the
perpendicular post-reset pair is exactly

    t_r=L00^2 a + L01^2 Y, c_r=L00 c, b_r=b.

One next OU prediction contracts aw by alpha.  A separately proved one-step
attitude PSD remainder eps is added to t_r and the exact Q_aw diagonal is added
to b.  For an aligned sample-1 signed force m e3 the repeated scalar channel is

    K_theta=(m t_1+c_1)/(m^2 t_1+2m c_1+b_1+r).

The denominator is additionally intersected with the exact SPD floor r.  This
producer subdivides p, first correction d, and the first axial aw correction.
It includes the dominant reset/yaw mixing and OU process terms.  Tiny body-rate
rotation of this channel, sample-1 tangent force, and the already-certified
1.8e-12-rad S correction are left as explicit later perturbations.  It neither
promotes the whole sample-1 branch nor P5.
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
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN; SCHEMA=1; RANGE=9.0

def _L_coeffs(d:Interval):
 s=Interval(VT.sin_point(d.lo).lo,VT.sin_point(d.hi).hi); c=Interval(VT.cos_point(d.hi).lo,VT.cos_point(d.lo).hi); h=FULL.I(.5)
 return c-s*h*d, -(c*h*d)-s

def build(domain_path=DEFAULT_DOMAIN,source_pieces=4,source_cell_index=0,p_pieces=32,d_pieces=32,axial_pieces=32):
 FULL3._install_backend(); path=Path(domain_path).resolve(); dom=json.loads(path.read_text()); first=FIRST.build(path,source_pieces=source_pieces); vec=VECTOR.build(); failures=[]
 failures += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
 src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
 if phase!="due": failures.append("expected due source witness")
 fr=first["source_cells"][source_cell_index]; p_all=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"])); rho0=float(fr["combined_useful_residual_norm_upper_mps2"]); dmax=float(fr["correction_norm_upper_rad"])
 h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tilt,yaw,eps=RG._attitude_covariance_epsilon(path,h); t=Interval.outward_bounds(tilt,FULL.up(tilt+eps)); Y=Interval.outward_bounds(yaw,FULL.up(yaw+eps)); Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); r=Racc[0][0]
 F,Q,_=FULL._transition_and_Q(src,dom); alpha=F[15][15]; qaw=Q[15][15]; rho1=FULL.up(float(dom["normal_live"]["specific_force_norm_upper_mps2"])+float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"]))
 pcells=SUB.parts(p_all.lo,p_all.hi,p_pieces); dcells=SUB.parts(0,dmax,d_pieces); rows=[]; bad=None; maxk=maxd=0.0; minS=math.inf; maxtr=0.0
 for pi,p in enumerate(pcells):
  D=FULL.I(g*g)*t+p+r; a=t*(p+r)/D; c0=-(FULL.I(g)*t*p/D); b=p*(FULL.I(g*g)*t+r)/D
  # First axial aw correction uses the exact scalar axial gain p/(p+r).
  kg=Interval(FULL.down(p.lo/(p.lo+r.hi)),FULL.up(p.hi/(p.hi+r.lo))); azmax=FULL.up(kg.hi*rho0); azcells=SUB.parts(-azmax,azmax,axial_pieces)
  for di,d in enumerate(dcells):
   l00,l01=_L_coeffs(d); tr=l00.square()*a+l01.square()*Y+Interval(0.0,FULL.up(eps)); cr=alpha*l00*c0; br=alpha.square()*b+qaw; maxtr=max(maxtr,tr.hi)
   for ai,az in enumerate(azcells):
    m=FULL.I(g)+alpha*az; num=m*tr+cr; S=m.square()*tr+FULL.I(2.0)*m*cr+br+r; slo=max(r.lo,S.lo)
    if not slo>0:
     failures.append(f"lost scalar SPD floor p={pi} d={di} a={ai}"); continue
    k=FULL.up(num.abs_upper()/slo); corr=FULL.up(k*rho1); closed=math.isfinite(corr) and corr<RANGE; maxk=max(maxk,k); maxd=max(maxd,corr); minS=min(minS,slo)
    row={"p_cell":pi,"d_cell":di,"axial_cell":ai,"P_aw_variance":p.as_list(),"first_correction_rad":d.as_list(),"first_axial_aw_correction_mps2":az.as_list(),"sample1_aligned_signed_force_mps2":m.as_list(),"post_reset_predicted_tangent_variance":tr.as_list(),"post_reset_predicted_theta_aw_covariance":cr.as_list(),"post_reset_predicted_aw_variance":br.as_list(),"innovation_variance_raw":S.as_list(),"innovation_variance_lower_used":slo,"Ktheta_abs_upper":k,"tangent_residual_norm_upper_mps2":rho1,"correction_norm_upper_rad":corr,"inside_9rad_range":closed}; rows.append(row)
    if not closed and bad is None:bad=row
 ok=bool(rows) and bad is None and not failures
 return {"schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_RESET_PERPENDICULAR_SCALAR_CHANNEL_WITNESS","source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,"first_scalar_Joseph_posterior_used":True,"exact_reset_plus_body_gauge_yaw_mixing_used":True,"one_step_attitude_PSD_remainder_included":True,"one_step_OU_process_variance_included":True,"aligned_force_core_only":True,"sample1_body_rate_rotation_perturbation_included":False,"sample1_tangent_force_perturbation_included":False,"sample1_S_attitude_correction_included":False,"complete_sample1_branch_closed_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"validated_deployed_quaternion_range_rad":RANGE,"evaluated_joint_cells":len(rows),"tangent_residual_norm_upper_mps2":rho1,"minimum_scalar_innovation_variance_lower":minS,"max_post_reset_predicted_tangent_variance":maxtr,"max_Ktheta_abs_upper":maxk,"max_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":bad,"P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"ADD_BODY_RATE_TANGENT_FORCE_AND_TINY_S_PERTURBATIONS_TO_RESET_PERP_SCALAR_CHANNEL" if ok else "TIGHTEN_POSITIVE_SCALAR_INNOVATION_FORMULA_AFTER_RESET","failures":failures,"rows":rows}

def validate(d):
 f=list(d.get("failures",[]))
 for k in ("source_generated_not_trajectory_fit","first_scalar_Joseph_posterior_used","exact_reset_plus_body_gauge_yaw_mixing_used","one_step_attitude_PSD_remainder_included","one_step_OU_process_variance_included","aligned_force_core_only"):
  if d.get(k) is not True:f.append(k)
 for k in ("source_replay_used","filter_changed","sample1_body_rate_rotation_perturbation_included","sample1_tangent_force_perturbation_included","sample1_S_attitude_correction_included","complete_sample1_branch_closed_here","whole_word_promoted_here","N_H_words_set_here"):
  if d.get(k) is not False:f.append(k)
 if int(d.get("evaluated_joint_cells",0))<=0:f.append("no cells")
 if not float(d.get("minimum_scalar_innovation_variance_lower",0))>0:f.append("nonpositive innovation floor")
 st=d.get("P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL"); w=d.get("first_unclosed_joint_cell")
 if st=="PASS" and w is not None:f.append("PASS retains witness")
 if st=="NOT_ESTABLISHED" and w is None:f.append("missing witness")
 return f

def main():
 a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=4); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--p-pieces",type=int,default=32); a.add_argument("--d-pieces",type=int,default=32); a.add_argument("--axial-pieces",type=int,default=32); a.add_argument("--output",type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces,x.source_cell_index,x.p_pieces,x.d_pieces,x.axial_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL"],"cells":d["evaluated_joint_cells"],"rho":d["tangent_residual_norm_upper_mps2"],"min_S":d["minimum_scalar_innovation_variance_lower"],"max_tvar":d["max_post_reset_predicted_tangent_variance"],"max_K":d["max_Ktheta_abs_upper"],"max_d":d["max_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__":raise SystemExit(main())
