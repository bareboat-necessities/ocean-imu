#!/usr/bin/env python3
"""Physical group-norm transport for the deployed H-first P5 sample-1 bridge.

The interval state boxes used by older P5 matrix code conservatively project
vector norm balls onto coordinates.  Reading those coordinate boxes back as a
Euclidean group norm pays an artificial sqrt(3): e.g. ||v||<=5 becomes the
cube [-5,5]^3 and is later reported as 8.66.  That is not a physical theorem
bound.

This producer carries the declared group norms in parallel with the interval
matrix enclosure.  The integrated-OU chain is axis-isotropic, hence the vector
norms obey the exact scalar 4x4 comparison

  V+ <= V + |phi_va| A
  P+ <= h V + P + |phi_pa| A
  S+ <= .5 h^2 V + h P + S + |phi_Sa| A
  A+ <= |alpha| A.

At the first accelerometer, the accepted branch adds at most the already
source-certified group correction operator cap; the identity branch adds zero,
so their hull is bounded by the accepted sum.  The optional first S=0 update
has exact zero estimator residual and therefore adds no state correction on
this first prefix.  A second prediction gives the sample-1 pre-measurement
physical group-norm caps.

Position is intentionally different: the new P5 assumption is componentwise
|delta p_i|<=0.5 Hs, so its physical norm cap is sqrt(3)/2 Hs.  No new hard
bounds are invented.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_p4_candidate_full_word as CAND
import ou3_p5_45deg_sample1_entry_ha as ENTRY
import ou3_p5_45deg_sample1_entry_ha_v2 as ENTRY2
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_startup_stability_certificate as P1

DEFAULT_DOMAIN=ENTRY.DEFAULT_DOMAIN
SCHEMA=1


def up(x): return math.nextafter(float(x),math.inf)

def _coeffs(F):
 # One representative axis; the deployed integrated-OU coefficients are shared
 # identically by x/y/z.
 return {
  'h':F[9][6].abs_upper(),
  'half_h2':F[12][6].abs_upper(),
  'phi_va':F[6][15].abs_upper(),
  'phi_pa':F[9][15].abs_upper(),
  'phi_Sa':F[12][15].abs_upper(),
  'alpha':F[15][15].abs_upper(),
 }

def _predict(c,k):
 V,P,S,A=(float(c[x]) for x in ('velocity','position','S','aw'))
 return {
  'velocity':up(V+up(k['phi_va']*A)),
  'position':up(up(k['h']*V)+up(P+up(k['phi_pa']*A))),
  'S':up(up(k['half_h2']*V)+up(k['h']*P)+up(S+up(k['phi_Sa']*A))),
  'aw':up(k['alpha']*A),
  'gyro_bias':float(c['gyro_bias']),
 }

def _correct(c,corr):
 out=dict(c)
 for name in ('gyro_bias','velocity','position','S','aw'):
  out[name]=up(float(c[name])+float(corr[name]))
 return out

def build(domain_path=DEFAULT_DOMAIN,source_pieces=2):
 path=Path(domain_path).resolve(); dom=json.loads(path.read_text()); raw=ENTRY.build(path,source_pieces); reach=ENTRY2._source_reachability(); failures=[]
 if not raw['modes']['H']['complete']: failures.append('H interval bridge did not reach sample1')
 if not all(reach.values()): failures.append('deployed H-first reachability markers failed')
 srcs=RG._source_phase_children(source_pieces); hrows=raw['modes']['H']['rows']
 if len(srcs)!=len(hrows): failures.append('source/H row count mismatch')
 old=dom['startup']['physical_handoff_coordinate_bounds']; pos=dom['initial_filter_entrance']['position']; hs=float(pos['significant_wave_height_Hs_upper_m']); pf=float(pos['component_abs_error_upper_Hs_factor']); sqrt3=P1.sqrt_interval_point(3.0).hi
 initial={
  'gyro_bias':float(old['gyro_bias_error_norm_upper_rad_s']),
  'velocity':float(old['velocity_error_norm_upper_mps']),
  'position':up(sqrt3*up(pf*hs)),
  'S':float(old['integral_displacement_error_norm_upper_m_s']),
  'aw':float(old['latent_acceleration_error_norm_upper_mps2']),
 }
 rows=[]; maxima={k:0.0 for k in initial}; max_box={k:0.0 for k in initial}; max_ratio={k:0.0 for k in initial}
 for ((src,phase),r) in zip(srcs,hrows):
  CAND._configure_mode('H'); F,_Q,_R,_bp=CAND._transition_and_Q('H',src,dom); k=_coeffs(F); pred0=_predict(initial,k); corr=r['group_correction_norm_caps']; post=_correct(pred0,corr); sample1=_predict(post,k); box=r['sample1_state_group_norm_uppers']
  ratios={name:(float(box[name])/sample1[name] if sample1[name]>0 else 1.0) for name in initial}
  for name in initial:
   maxima[name]=max(maxima[name],sample1[name]); max_box[name]=max(max_box[name],float(box[name])); max_ratio[name]=max(max_ratio[name],ratios[name])
  rows.append({'source_phase_cell':r['source_phase_cell'],'pseudo_phase':phase,'coefficients':k,'initial_group_norm_caps':initial,'after_first_prediction_group_norm_caps':pred0,'first_accel_group_correction_norm_caps':corr,'after_accepted_or_identity_hull_group_norm_caps':post,'sample1_group_norm_caps':sample1,'old_interval_box_reported_group_norms':box,'old_box_to_physical_norm_ratio':ratios})
 for name in ('velocity','S','aw'):
  if not maxima[name] < max_box[name]: failures.append(f'{name} physical norm transport did not improve old box readback')
 ok=bool(rows) and not failures
 return {'schema':SCHEMA,'qualification':'OU3_P5_45DEG_H_FIRST_SAMPLE1_PHYSICAL_GROUP_NORM_TRANSPORT','source_generated_not_trajectory_fit':True,'source_replay_used':False,'filter_changed':False,'deployed_first_live_mode':'H','norm_balls_not_reinterpreted_as_cartesian_cubes':True,'position_componentwise_half_Hs_converted_to_sqrt3_norm':True,'integrated_OU_axis_isotropy_used':True,'first_due_S_mean_correction_exact_zero':True,'accepted_or_identity_state_norm_hull_used':True,'new_hard_bounds_invented':False,'initial_group_norm_caps':initial,'sample1_max_group_norm_caps':maxima,'old_interval_box_max_reported_group_norms':max_box,'max_old_box_to_physical_norm_ratio':max_ratio,'rows':rows,'P5_45DEG_H_SAMPLE1_GROUP_NORM_CERTIFICATE':'PASS' if ok else 'NOT_ESTABLISHED','next_obligation':'use these physical group-norm caps, together with the source-correlated covariance child, in the signed sample1 S/accelerometer/magnetometer recapture calculation toward the 30deg P4 sector','failures':failures}

def validate(d):
 f=list(d.get('failures',[]))
 for k in ('source_generated_not_trajectory_fit','norm_balls_not_reinterpreted_as_cartesian_cubes','position_componentwise_half_Hs_converted_to_sqrt3_norm','integrated_OU_axis_isotropy_used','first_due_S_mean_correction_exact_zero','accepted_or_identity_state_norm_hull_used'):
  if d.get(k) is not True:f.append(k)
 for k in ('source_replay_used','filter_changed','new_hard_bounds_invented'):
  if d.get(k) is not False:f.append(k)
 if d.get('deployed_first_live_mode')!='H':f.append('first mode')
 if d.get('P5_45DEG_H_SAMPLE1_GROUP_NORM_CERTIFICATE')=='PASS' and not d.get('rows'):f.append('no rows')
 return list(dict.fromkeys(f))

def main():
 a=argparse.ArgumentParser(); a.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); a.add_argument('--source-pieces',type=int,default=2); a.add_argument('--output',type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces); vf=validate(d); d['validation_pass']=not vf; d['validation_failures']=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({'status':d['P5_45DEG_H_SAMPLE1_GROUP_NORM_CERTIFICATE'],'initial':d['initial_group_norm_caps'],'sample1':d['sample1_max_group_norm_caps'],'old_box':d['old_interval_box_max_reported_group_norms'],'ratios':d['max_old_box_to_physical_norm_ratio'],'next':d['next_obligation'],'validation_failures':vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=='__main__': raise SystemExit(main())
