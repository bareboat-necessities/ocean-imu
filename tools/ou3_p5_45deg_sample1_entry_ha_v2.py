#!/usr/bin/env python3
"""Projection-aware H/A sample-1 entry from the 45 deg P5 bridge.

V1 deliberately failed if the active accelerometer-bias enclosure crossed the
0.5 m/s^2 shipping projection surface.  That is appropriate for the smooth P4
inner word, whose A source node is restricted to 0.45, but it is too strong for
outer P5: the shipping operation at an accepted injection is the exact Euclidean
projection onto the closed 0.5 ball.  Projection is nonexpansive and gives the
post-injection state norm cap min(pre_projection_cap, 0.5).

This wrapper changes only that proof enclosure.  The rejected branch retains its
pre-update cap, and the accepted/rejected hull is therefore bounded by the
maximum of the rejected cap and the projected accepted cap.

Source reachability is recorded separately.  The deployed wrapper defaults to
with_mag=true and goLive(..., allow_acc_bias=false).  enterCold() therefore
locks the accelerometer bias; enterLive() keeps it held; updateMag() can release
it only after the configured magnetometer-update count and >1 s guard.  Thus the
deployed first-Live packet is H=18.  The A=21 sample-zero calculation emitted
here is a projection-aware configuration diagnostic and a bound useful for the
later H->A hybrid transition; it is not claimed to be a deployed first-Live
branch.
"""
from __future__ import annotations
import argparse,json,math,re
from pathlib import Path
from ou3_interval import Interval
import ou3_p5_45deg_sample1_entry_ha as V1

REPO=Path(__file__).resolve().parents[1]
WRAPPER=REPO/'src'/'kalman_ou_iii'/'SeaStateFusionFilter_OU_III.h'
MEKF=REPO/'src'/'kalman_ou_iii'/'Kalman3D_Wave_OU_III.h'
DEFAULT_DOMAIN=V1.DEFAULT_DOMAIN
SCHEMA=2


def up(x): return math.nextafter(float(x),math.inf)


def _source_reachability():
 w=WRAPPER.read_text(encoding='utf-8'); k=MEKF.read_text(encoding='utf-8')
 checks={
  'deployed_config_with_mag_default_true':'bool with_mag = true;' in w,
  'goLive_allow_acc_bias_default_false':re.search(r'void\s+goLive\([^)]*bool\s+allow_acc_bias\s*=\s*false',w,re.S) is not None,
  'cold_lock_follows_with_mag':'accel_bias_locked_   = with_mag_;' in w,
  'live_enable_requires_lock_cleared':'const bool allow_bias = !accel_bias_locked_ && !acc_bias_hold_;' in w,
  'mag_unlock_requires_update_count':'mag_updates_applied_ >= mag_updates_to_unlock_' in w,
  'mag_unlock_has_one_second_guard':re.search(r'time_\)\s*-\s*first_mag_update_time_\)\s*>\s*1\.0f',w) is not None,
  'deployed_unlock_count_250':'int acc_bias_unlock_mag_updates = 250;' in w and 'MAG_UPDATES_TO_UNLOCK = 250' in w,
  'shipping_projection_is_radial_closed_ball':'if (n > acc_bias_limit_)' in k and 'b *= (acc_bias_limit_ / n);' in k,
 }
 return checks


def _acc_child_projection(P,e,mode,H,R,rho,dcap,ba_cap,proj):
 PHt,S=V1.FULL._innovation(P,H,R); Sinv,b=V1.FULL._spd_inverse_enclosure(S,R); K=V1.FULL.matrix_mul(PHt,Sinv); dx=V1.FULL._mat_vec(K,V1.FULL._vec_box(rho)); caps=V1._caps(PHt,P,R,mode,rho); dx2=list(dx)
 dc=Interval(-up(dcap),up(dcap))
 for i in range(3): dx2[i]=V1.FULL._intersect(dx2[i],dc)
 gm={'gyro_bias':tuple(V1.FULL.BG),'velocity':tuple(V1.FULL.V),'position':tuple(V1.FULL.P),'S':tuple(V1.FULL.SS),'aw':tuple(V1.FULL.AW)}
 if mode=='A': gm['accel_bias']=tuple(V1.FULL.BA)
 for name,idxs in gm.items():
  c=Interval(-caps[name],caps[name])
  for i in idxs: dx2[i]=V1.FULL._intersect(dx2[i],c)
 Pj=V1.FULL._shipping_joseph(P,K,S,PHt); Pr=V1.FULL._reset_covariance(Pj,dx2[:3]); ea=list(e)
 for i in range(3,V1.FULL.N): ea[i]=e[i]-dx2[i]
 ba2=ba_cap; preproj=None; active=False
 if mode=='A':
  preproj=up(float(ba_cap)+caps['accel_bias'])
  accepted_cap=min(preproj,float(proj))
  rejected_cap=float(ba_cap)
  ba2=up(max(accepted_cap,rejected_cap))
  active=preproj>float(proj)
  for i in V1.FULL.BA: ea[i]=V1.FULL._intersect(ea[i],Interval(-ba2,ba2))
 return V1.FULL._psd_tighten(V1.FULL._mat_hull(P,Pr)),V1.FULL._vec_hull(e,ea),ba2,{'backend':b,'floor':V1._floor(P,R,mode),'caps':caps,'A_pre_projection_norm_cap':preproj,'A_projection_active_possible':active,'A_post_projection_or_rejected_hull_norm_cap':ba2}


def build(domain_path=DEFAULT_DOMAIN,source_pieces=2):
 checks=_source_reachability(); old=V1._acc_child; V1._acc_child=_acc_child_projection
 try: out=dict(V1.build(Path(domain_path).resolve(),source_pieces))
 finally: V1._acc_child=old
 failures=[]
 if not all(checks.values()): failures += [f'source reachability missing: {k}' for k,v in checks.items() if not v]
 # V1 failures caused solely by its old A projection refusal disappear through
 # the patched child.  Any remaining failure is still real and is retained.
 failures += list(out.get('failures',[]))
 out.update({
  'schema':SCHEMA,
  'qualification':'OU3_P5_45DEG_SAMPLE1_ENTRY_HA_PROJECTION_AWARE_SOURCE_REACHABILITY',
  'source_reachability':checks,
  'deployed_default_first_live_mode':'H',
  'deployed_default_A_first_live_reachable':False,
  'A_sample0_evaluation_role':'CONFIGURATION_DIAGNOSTIC_AND_LATER_H_TO_A_BOUND',
  'A_projection_surface_forbidden_in_outer_P5':False,
  'A_shipping_closed_ball_projection_used':True,
  'A_projection_nonexpansive_used':True,
  'A_P4_inner_045_ball_promoted_here':False,
  'failures':list(dict.fromkeys(failures)),
 })
 ok=out['sample1_entry_inside_q8'] and out['modes']['H']['complete'] and out['modes']['A']['complete'] and not out['failures']
 out['P5_45DEG_SAMPLE1_ENTRY_HA_V2_CERTIFICATE']='PASS' if ok else 'NOT_ESTABLISHED'
 out['next_obligation']='continue the deployed H child through source-correlated sample1 operations and magnetic unlock; at the actual H->A jump use the projection-aware A outer cap and prove finite recapture from <=0.5 into the 0.45 smooth P4 A node'
 return out


def validate(d):
 f=list(d.get('failures',[]))
 if d.get('schema')!=SCHEMA:f.append('schema')
 if d.get('deployed_default_first_live_mode')!='H':f.append('first mode')
 if d.get('deployed_default_A_first_live_reachable') is not False:f.append('A first-live reachability')
 for k in ('A_shipping_closed_ball_projection_used','A_projection_nonexpansive_used'):
  if d.get(k) is not True:f.append(k)
 for k in ('A_projection_surface_forbidden_in_outer_P5','A_P4_inner_045_ball_promoted_here'):
  if d.get(k) is not False:f.append(k)
 if not all(d.get('source_reachability',{}).values()):f.append('source reachability')
 for m in ('H','A'):
  if d.get('modes',{}).get(m,{}).get('complete') is not True:f.append(f'{m} incomplete')
 if d.get('sample1_entry_inside_q8') is not True:f.append('q8')
 return list(dict.fromkeys(f))


def main():
 a=argparse.ArgumentParser(); a.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); a.add_argument('--source-pieces',type=int,default=2); a.add_argument('--output',type=Path,required=True); x=a.parse_args(); d=build(x.domain,x.source_pieces); vf=validate(d); d['validation_pass']=not vf; d['validation_failures']=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({'status':d['P5_45DEG_SAMPLE1_ENTRY_HA_V2_CERTIFICATE'],'first_mode':d['deployed_default_first_live_mode'],'q1':d['sample1_pre_measurement_q_upper'],'H_complete':d['modes']['H']['complete'],'A_diagnostic_complete':d['modes']['A']['complete'],'A_max_state':d['modes']['A']['max_state'],'A_failure':d['modes']['A']['first_failure'],'next':d['next_obligation'],'validation_failures':vf},indent=2,sort_keys=True)); return 0 if not vf else 2
if __name__=='__main__': raise SystemExit(main())
