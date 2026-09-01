#!/usr/bin/env python3
"""Propagate the 45 deg P5 first-packet bridge to sample 1 in H/A modes.

Starts from the sign-complete 45 deg -> q<8 bridge, not the old global q8
cube.  Every first-prefix source/phase cell is propagated with the shipping
Joseph update, immediate left-error reset, accepted/rejected hull, 0.5 Hs
position entrance, and fixed dimensions H=18 / A=21.  The first accelerometer
attitude correction uses the source-certified norm cap; linear/bias state
corrections use group operator-norm caps from P C^T and the first innovation
Loewner floor.  The producer stops at sample-1 pre-measurement entry and does
not claim recapture into the 30 deg P4 sector.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval,matrix_abs_col_sum_upper,matrix_abs_row_sum_upper
import ou3_p4_candidate_first_accel_exact_source as FIRST
import ou3_p4_candidate_full_word as CAND
import ou3_p5_45deg_first_accel_q8_bridge as BRIDGE
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN; SCHEMA=1

def up(x): return math.nextafter(float(x),math.inf)

def _op2(A): return up(math.sqrt(up(matrix_abs_col_sum_upper(A)*matrix_abs_row_sum_upper(A))))

def _canon_att(P,path,h):
 tilt,yaw,eps=RG._attitude_covariance_epsilon(path,h); Pt=RG._ptheta_cell([FULL.I(0),FULL.I(0),FULL.I(1)],tilt,yaw,eps)
 out=[[P[i][j] for j in range(FULL.N)] for i in range(FULL.N)]
 for i in range(3):
  for j in range(3): out[i][j]=Pt[i][j]
 return FULL._psd_tighten(out)

def _zero_S_cov(P,src):
 H=FULL._H_S(); R=FULL._R_S(src); PHt,S=FULL._innovation(P,H,R); Sinv,b=FULL._spd_inverse_enclosure(S,R); K=FULL.matrix_mul(PHt,Sinv)
 return FULL._shipping_joseph(P,K,S,PHt),b

def _Hacc(mode,g):
 H=FULL._zero(3,FULL.N); gg=FULL.I(g); H[0][1]=gg; H[1][0]=-gg
 for i in range(3): H[i][15+i]=FULL.I(1)
 if mode=='A':
  for i in range(3): H[i][18+i]=FULL.I(1)
 return H

def _floor(P,R,mode):
 r=min(R[i][i].lo for i in range(3)); paw=min(max(0.0,P[15+i][15+i].lo) for i in range(3)); pba=0.0
 if mode=='A': pba=min(max(0.0,P[18+i][18+i].lo) for i in range(3))
 x=FULL.down(paw+pba+r)
 if not x>0: raise RuntimeError('first accelerometer innovation Loewner floor lost')
 return x

def _caps(PHt,P,R,mode,rho):
 fl=_floor(P,R,mode); groups={'gyro_bias':tuple(FULL.BG),'velocity':tuple(FULL.V),'position':tuple(FULL.P),'S':tuple(FULL.SS),'aw':tuple(FULL.AW)}
 if mode=='A': groups['accel_bias']=tuple(FULL.BA)
 out={}
 for name,idxs in groups.items():
  B=[[PHt[i][j] for j in range(3)] for i in idxs]; k=up(_op2(B)/fl)
  if name=='aw': k=min(1.0,k)
  out[name]=up(k*rho)
 return out

def _acc_child(P,e,mode,H,R,rho,dcap,ba_cap,proj):
 PHt,S=FULL._innovation(P,H,R); Sinv,b=FULL._spd_inverse_enclosure(S,R); K=FULL.matrix_mul(PHt,Sinv); dx=FULL._mat_vec(K,FULL._vec_box(rho)); caps=_caps(PHt,P,R,mode,rho); dx2=list(dx)
 dc=Interval(-up(dcap),up(dcap))
 for i in range(3): dx2[i]=FULL._intersect(dx2[i],dc)
 gm={'gyro_bias':tuple(FULL.BG),'velocity':tuple(FULL.V),'position':tuple(FULL.P),'S':tuple(FULL.SS),'aw':tuple(FULL.AW)}
 if mode=='A': gm['accel_bias']=tuple(FULL.BA)
 for name,idxs in gm.items():
  c=Interval(-caps[name],caps[name])
  for i in idxs: dx2[i]=FULL._intersect(dx2[i],c)
 Pj=FULL._shipping_joseph(P,K,S,PHt); Pr=FULL._reset_covariance(Pj,dx2[:3]); ea=list(e)
 for i in range(3,FULL.N): ea[i]=e[i]-dx2[i]
 ba2=ba_cap
 if mode=='A':
  ba2=up(float(ba_cap)+caps['accel_bias'])
  if not ba2<proj: raise RuntimeError(f'A bias enclosure reaches projection surface: {ba2} >= {proj}')
  for i in FULL.BA: ea[i]=FULL._intersect(ea[i],Interval(-ba2,ba2))
 return FULL._psd_tighten(FULL._mat_hull(P,Pr)),FULL._vec_hull(e,ea),ba2,{'backend':b,'floor':_floor(P,R,mode),'caps':caps}

def _gn(e,mode):
 d={'gyro_bias':FULL._norm_upper([e[i] for i in FULL.BG]),'velocity':FULL._norm_upper([e[i] for i in FULL.V]),'position':FULL._norm_upper([e[i] for i in FULL.P]),'S':FULL._norm_upper([e[i] for i in FULL.SS]),'aw':FULL._norm_upper([e[i] for i in FULL.AW])}
 if mode=='A': d['accel_bias']=FULL._norm_upper([e[i] for i in FULL.BA])
 return d

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2):
 FULL3._install_backend(); path=Path(domain_path).resolve(); dom=json.loads(path.read_text()); first=FIRST.build(path,source_pieces=source_pieces); bridge=BRIDGE.build(path,source_pieces=source_pieces); vec=VECTOR.build(); failures=[]
 failures += [f'first: {x}' for x in FIRST.validate(first)] + [f'bridge: {x}' for x in BRIDGE.validate(bridge)] + [f'vector: {x}' for x in VECTOR.validate(vec)]
 srcs=RG._source_phase_children(source_pieces); base=first['candidate_rows'][0]['source_rows']; h=float(FULL._source_cell()['dt_s']); g=float(dom['startup']['gravity_mps2']); baH=float(dom['startup']['physical_handoff_coordinate_bounds']['accelerometer_bias_error_norm_upper_mps2']); qpred=float(bridge['P5_45deg_entrance_first_accel']['post_prediction_q_upper']); qpost=float(bridge['P5_45deg_entrance_first_accel']['product']['post_update_q_upper_from_scalar']); q1=RG._q_after_first_prediction(qpost,dom,h); qsafe=math.isfinite(q1) and q1<8.0
 if len(srcs)!=len(base): failures.append('source-phase ordering/count mismatch')
 if not qsafe: failures.append('sample1 prediction leaves q8')
 Racc=FULL._R_diag(float(vec['configured_measurement_bounds']['acc_measurement_std_mps2'])); proj=float(dom['normal_live']['active_accelerometer_bias_projection_limit_mps2']); modes={}
 for mode in ('H','A'):
  CAND._configure_mode(mode); rows=[]; bad=None; maxstate={}; maxdiag=maxcross=0.0
  for si,((src,phase),sr) in enumerate(zip(srcs,base)):
   try:
    F,Q,_R,_bp=CAND._transition_and_Q(mode,src,dom); P0=CAND._initial_covariance(mode,src,path); e0,ba,pos=CAND._initial_error(mode,dom); Pp=FULL._psd_tighten(FULL.matrix_add(FULL.matrix_mul(FULL.matrix_mul(F,P0),FULL.matrix_transpose(F)),Q)); ep,ba=CAND._predict_error(mode,e0,F,ba); Pp=_canon_att(Pp,path,h); sb='NOT_DUE_IDENTITY'; Ppre=Pp
    if phase=='due': Ppre,sb=_zero_S_cov(Pp,src)
    rot=up(g*qpred); aw=float(sr['predicted_aw_error_norm_upper_mps2']); rho=up(rot+up(aw+baH)); dcap=up(float(sr['Ktheta_norm_upper'])*rho); Pa,ea,ba2,a=_acc_child(Ppre,ep,mode,_Hacc(mode,g),Racc,rho,dcap,ba,proj); P1=FULL._psd_tighten(FULL.matrix_add(FULL.matrix_mul(FULL.matrix_mul(F,Pa),FULL.matrix_transpose(F)),Q)); e1,ba1=CAND._predict_error(mode,ea,F,ba2); sm=FULL._matrix_summary(P1); gn=_gn(e1,mode); maxdiag=max(maxdiag,max(x[1] for x in sm['diagonal_intervals'])); maxcross=max(maxcross,sm['max_offdiagonal_abs_upper'])
    for k,v in gn.items(): maxstate[k]=max(maxstate.get(k,0.0),v)
    rows.append({'source_phase_cell':si,'pseudo_phase':phase,'position_entrance':pos,'residual_norm_upper_mps2':rho,'attitude_correction_norm_upper_rad':dcap,'S_backend':sb,'acc_backend':a['backend'],'innovation_floor':a['floor'],'group_correction_norm_caps':a['caps'],'sample1_state_group_norm_uppers':gn,'sample1_covariance':sm,'sample1_active_bias_norm_cap':ba1})
   except Exception as exc:
    bad={'mode':mode,'source_phase_cell':si,'pseudo_phase':phase,'reason':f'{type(exc).__name__}: {exc}'}; break
  complete=len(rows)==len(srcs) and bad is None; modes[mode]={'dimension':18 if mode=='H' else 21,'evaluated':len(rows),'expected':len(srcs),'complete':complete,'max_state':maxstate,'max_cov_diag':maxdiag,'max_cov_cross':maxcross,'first_failure':bad,'rows':rows}
  if not complete: failures.append(f'{mode} did not reach sample1 for every source-phase cell')
 ok=qsafe and all(modes[m]['complete'] for m in ('H','A')) and not failures
 return {'schema':SCHEMA,'qualification':'OU3_P5_45DEG_FIRST_ACCEL_JOSEPH_RESET_TO_SAMPLE1_HA','source_generated_not_trajectory_fit':True,'source_replay_used':False,'filter_changed':False,'starts_from_45deg_sign_complete_q8_bridge':True,'position_entrance_uses_half_Hs':True,'shipping_Joseph_update_used':True,'shipping_left_error_reset_congruence_used':True,'accepted_and_identity_branches_hulled':True,'first_due_S_zero_mean_covariance_update_retained':True,'H_dimension':18,'A_dimension':21,'A_Jba_identity_retained':True,'state_corrections_use_group_operator_norm_caps':True,'raw_residual_component_cube_not_used_as_group_norm':True,'sample0_post_update_q_upper':qpost,'sample1_pre_measurement_q_upper':q1,'sample1_entry_inside_q8':qsafe,'modes':modes,'sample1_measurement_prefix_evaluated_here':False,'returned_to_30deg_P4_sector_here':False,'P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE':False,'P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE':False,'P5_45DEG_SAMPLE1_ENTRY_HA_CERTIFICATE':'PASS' if ok else 'NOT_ESTABLISHED','next_obligation':'evaluate sample1 source-correlated S/accelerometer/magnetometer branches from these H/A children and test recapture toward 30deg' if ok else 'refine the first reported H/A child without returning to the global q8 cube','failures':list(dict.fromkeys(failures))}

def validate(d):
 f=list(d.get('failures',[]))
 for k in ('source_generated_not_trajectory_fit','starts_from_45deg_sign_complete_q8_bridge','position_entrance_uses_half_Hs','shipping_Joseph_update_used','shipping_left_error_reset_congruence_used','accepted_and_identity_branches_hulled','first_due_S_zero_mean_covariance_update_retained','A_Jba_identity_retained','state_corrections_use_group_operator_norm_caps','raw_residual_component_cube_not_used_as_group_norm'):
  if d.get(k) is not True:f.append(k)
 for k in ('source_replay_used','filter_changed','sample1_measurement_prefix_evaluated_here','returned_to_30deg_P4_sector_here','P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE','P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE'):
  if d.get(k) is not False:f.append(k)
 if d.get('H_dimension')!=18 or d.get('A_dimension')!=21:f.append('H/A dimensions')
 if d.get('sample1_entry_inside_q8') is not True:f.append('sample1 q8')
 for mode,dim in (('H',18),('A',21)):
  m=d.get('modes',{}).get(mode,{})
  if m.get('dimension')!=dim or m.get('complete') is not True or m.get('evaluated')!=m.get('expected') or m.get('first_failure') is not None:f.append(f'{mode} incomplete')
 return list(dict.fromkeys(f))

def main():
 a=argparse.ArgumentParser(); a.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); a.add_argument('--source-pieces',type=int,default=2); a.add_argument('--output',type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces); vf=validate(d); d['validation_pass']=not vf; d['validation_failures']=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({'status':d['P5_45DEG_SAMPLE1_ENTRY_HA_CERTIFICATE'],'q1':d['sample1_pre_measurement_q_upper'],'H':{k:d['modes']['H'][k] for k in ('complete','max_state','first_failure')},'A':{k:d['modes']['A'][k] for k in ('complete','max_state','first_failure')},'next':d['next_obligation'],'validation_failures':vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=='__main__': raise SystemExit(main())
