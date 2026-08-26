#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_entry as ENTRY
import ou3_p5_sample1_entry_v3 as ENTRY3
import ou3_p5_sample1_prefix_v2 as PREFIX2
import ou3_p5_sample1_rotation_gauge_refinement as BASE
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=BASE.DEFAULT_DOMAIN; SCHEMA=2; N=FULL.N; LIMIT=6.0

def parts(lo,hi,n):
    return [Interval.outward_bounds(lo+(hi-lo)*i/n,lo+(hi-lo)*(i+1)/n) for i in range(n)]

def rx(d):
    s=Interval(VT.sin_point(d.lo).lo,VT.sin_point(d.hi).hi)
    c=Interval(VT.cos_point(d.hi).lo,VT.cos_point(d.lo).hi)
    z=FULL.I(0); o=FULL.I(1)
    return [[o,z,z],[z,c,-s],[z,s,c]]

def first_child(P,e,x,H,R,rho,d):
    PHt,S=FULL._innovation(P,H,R); Sinv,b=FULL._spd_inverse_enclosure(S,R); K=matrix_mul(PHt,Sinv)
    dx=FULL._mat_vec(K,FULL._vec_box(rho)); caps=ENTRY3._linear_gain_caps(P,R,rho)
    for name,idxs in ENTRY3.GROUPS.items():
        cap=Interval(-caps[name],caps[name])
        for i in idxs: dx[i]=FULL._intersect(dx[i],cap)
    Pj=FULL._shipping_joseph(P,K,S,PHt); Pr=FULL._reset_covariance(Pj,[d,FULL.I(0),FULL.I(0)])
    ee=list(e); xx=list(x)
    for i in range(3,N): ee[i]=e[i]-dx[i]; xx[i]=x[i]+dx[i]
    for i in FULL.TH: xx[i]=FULL.I(0)
    Rg=rx(d)
    return BASE._transform_covariance(Pr,Rg),BASE._transform_state(ee,Rg),BASE._transform_state(xx,Rg),b

def minabs(a):
    return 0.0 if a.lo<=0<=a.hi else min(abs(a.lo),abs(a.hi))

def force_cells(f,n):
    aa=parts(-f,f,n); r2=FULL.up(f*f)
    for x in aa:
      for y in aa:
       for z in aa:
        if FULL.down(minabs(x)**2+minabs(y)**2+minabs(z)**2)<=r2: yield x,y,z

def Hforce(v):
    x,y,z=v; H=FULL._zero(3,N)
    H[0][1]=z; H[0][2]=-y; H[1][0]=-z; H[1][2]=x; H[2][0]=y; H[2][1]=-x
    for i in range(3): H[i][15+i]=FULL.I(1)
    return H

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2,source_cell_index=0,delta_pieces=4,force_pieces=4):
    FULL3._install_backend(); p=Path(domain_path).resolve(); dom=json.loads(p.read_text())
    first=FIRST.build(p,source_pieces=source_pieces); vec=VECTOR.build(); fails=[]
    fails += [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
    src,phase0=RG._source_phase_children(source_pieces)[source_cell_index]
    if phase0!="due": fails.append("expected first due witness")
    row0=first["source_cells"][source_cell_index]; dmax=float(row0["correction_norm_upper_rad"]); rho0=float(row0["combined_useful_residual_norm_upper_mps2"])
    h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"]); tf=float(dom["normal_live"]["specific_force_norm_upper_mps2"]); ba=float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    Racc=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])); H0=ENTRY._canonical_first_H(g)
    F,Q,Rstep=FULL._transition_and_Q(src,dom); P0=FULL._initial_covariance(src,p)
    Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q)); Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,p,h)
    ep=FULL._predict_error(FULL._initial_error(dom),F); x0=[FULL.I(0) for _ in range(N)]
    Ppre,bS0=ENTRY._zero_residual_S_covariance(Pp,src); fixed=int(bS0=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"); fallback=1-fixed
    qpre=float(first["post_prediction_full_cayley_norm_upper"]); rows=[]; bad=None; maxd=maxrho=0.0
    for di,d in enumerate(parts(0,dmax,delta_pieces)):
        Pa,ea,xa,b0=first_child(Ppre,ep,x0,H0,Racc,rho0,d); fixed+=int(b0=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"); fallback+=int(b0!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN")
        q0=PREFIX2._post_correction_q_upper(qpre,d.hi); Pm=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pa),matrix_transpose(F)),Q)); em=BASE._predict_state(ea,F); xm=BASE._predict_state(xa,F); q1=RG._q_after_first_prediction(q0,dom,h)
        P1=BASE._transform_covariance(Pm,Rstep); e1=BASE._transform_state(em,Rstep); x1=BASE._transform_state(xm,Rstep)
        Sc=FULL._measurement_cell(P1,FULL._H_S(),FULL._R_S(src),[-x1[12+i] for i in range(3)]); fixed+=int(Sc["inverse_backend"]=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"); fallback+=int(Sc["inverse_backend"]!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN")
        P2=Sc["P_accepted"]; e2=list(e1); x2=list(x1)
        for i in range(3,N): e2[i]=e1[i]-Sc["dx"][i]; x2[i]=x1[i]+Sc["dx"][i]
        q2=PREFIX2._post_correction_q_upper(q1,FULL._norm_upper(Sc["dx"][0:3])); xaw=BASE._norm(x2,FULL.AW); eaw=BASE._norm(e2,FULL.AW); fhat=FULL.up(g+xaw); rho=FULL.up(FULL.up(BASE._rot_diff(q2)*tf)+FULL.up(eaw+ba)); maxrho=max(maxrho,rho)
        for fi,fc in enumerate(force_cells(fhat,force_pieces)):
            PHt,S=FULL._innovation(P2,Hforce(fc),Racc); Sinv,b=FULL._spd_inverse_enclosure(S,Racc); K=matrix_mul(PHt,Sinv); kn=RG._op2_upper([r[:] for r in K[:3]]); d1=FULL.up(kn*rho); qa=PREFIX2._post_correction_q_upper(q2,d1); maxd=max(maxd,d1)
            fixed+=int(b=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"); fallback+=int(b!="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN")
            closed=b=="FIXED_PIVOT_INTERVAL_GAUSS_JORDAN" and d1<=LIMIT and math.isfinite(qa) and qa<8
            r={"delta_cell":di,"delta_rad":d.as_list(),"force_cell":fi,"force_components":[a.as_list() for a in fc],"backend":b,"Ktheta_norm_upper":kn,"residual_norm_upper_mps2":rho,"correction_norm_upper_rad":d1,"post_q_upper":qa,"closed":closed}; rows.append(r)
            if not closed and bad is None: bad=r
    ok=bool(rows) and bad is None and not fails
    return {"schema":SCHEMA,"source_generated_not_trajectory_fit":True,"filter_changed":False,"sample0_canonical_correction_magnitude_subdivided":True,"sample1_force_component_cube_subdivided":True,"sample1_H_theta_exact_skew_structure_retained":True,"sample1_J_aw_exact_identity_in_transported_body_gauge":True,"delta_pieces":delta_pieces,"force_pieces_per_axis":force_pieces,"evaluated_joint_cells":len(rows),"fixed_pivot_inverse_count":fixed,"spectral_fallback_inverse_count":fallback,"max_sample1_residual_norm_upper_mps2":maxrho,"max_sample1_correction_norm_upper_rad":maxd,"first_unclosed_joint_cell":bad,"deployed_correction_limit_rad":LIMIT,"deployed_correction_limit_increased":False,"complete_source_cell_refined_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,"P5_SAMPLE1_DELTA_FORCE_SUBDIVIDED_WITNESS_REFINEMENT":"PASS" if ok else "NOT_ESTABLISHED","next_obligation":"LIFT_SUBDIVIDED_GAUGE_TO_REMAINING_SAMPLE1_BRANCHES" if ok else "INCREASE_SUBDIVISION_OR_DERIVE_DIRECTIONAL_SCHUR_GAIN_BOUND","failures":fails}

def validate(d):
    f=list(d.get("failures",[]))
    for k in ("source_generated_not_trajectory_fit","sample0_canonical_correction_magnitude_subdivided","sample1_force_component_cube_subdivided","sample1_H_theta_exact_skew_structure_retained","sample1_J_aw_exact_identity_in_transported_body_gauge"):
        if d.get(k) is not True: f.append(k)
    for k in ("filter_changed","deployed_correction_limit_increased","complete_source_cell_refined_here","whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False: f.append(k)
    if d.get("evaluated_joint_cells",0)<=0: f.append("no cells")
    st=d.get("P5_SAMPLE1_DELTA_FORCE_SUBDIVIDED_WITNESS_REFINEMENT"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and w is not None: f.append("PASS retains witness")
    if st=="NOT_ESTABLISHED" and w is None: f.append("missing witness")
    return f

def main():
    a=argparse.ArgumentParser(); a.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); a.add_argument("--source-pieces",type=int,default=2); a.add_argument("--source-cell-index",type=int,default=0); a.add_argument("--delta-pieces",type=int,default=4); a.add_argument("--force-pieces",type=int,default=4); a.add_argument("--output",type=Path,required=True); x=a.parse_args()
    d=build(x.domain,x.source_pieces,x.source_cell_index,x.delta_pieces,x.force_pieces); vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"status":d["P5_SAMPLE1_DELTA_FORCE_SUBDIVIDED_WITNESS_REFINEMENT"],"cells":d["evaluated_joint_cells"],"fixed":d["fixed_pivot_inverse_count"],"fallback":d["spectral_fallback_inverse_count"],"max_rho":d["max_sample1_residual_norm_upper_mps2"],"max_d":d["max_sample1_correction_norm_upper_rad"],"first_unclosed":d["first_unclosed_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=="__main__": raise SystemExit(main())
