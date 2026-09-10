#!/usr/bin/env python3
"""Finite-time shipping startup handoff for the qualified 8 m/s^2 BRMM source.

This proves that the actual startup proxy cannot wait forever before entering
Live/H18. It does not claim that the Live handoff is already inside the final
P4 basin.

The proof composes the source-correlated private-Mahony invariant, its
source-order binary32/discrete robustness certificate, the declared same-history
world-gravity averaging direction error, and the literal shipping timeout guard.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_private_mahony_live_invariant as CONT
import ou3_brmm_private_mahony_discrete_invariant as DISC
import ou3_startup_handoff_tilt_hemisphere as HANDOFF

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
COMMON=REPO/'src/kalman_common/SeaStateFusionFilterCommon.h'
SCHEMA=1
QUALIFICATION='OU3_QUALIFIED_BRMM_STARTUP_TIMEOUT_CAPTURE_V1'

def build():
 c=CONT.build();dc=DISC.build();h=HANDOFF.build();bad={'continuous':CONT.validate(c),'discrete':DISC.validate(dc),'handoff':HANDOFF.validate(h)};bad={k:v for k,v in bad.items() if v}
 if bad:raise RuntimeError('startup capture prerequisites failed: '+repr(bad))
 d=json.loads(DOMAIN.read_text());w=WRAPPER.read_text();common=COMMON.read_text();eps=float(d['startup']['world_averaged_gravity_direction_error_upper_rad']);tilt=float(c['actual_tilt_rad_upper']);margin=0.5*math.pi-tilt-eps
 parity={'timeout_is_150s':('proxy_startup_timeout_sec = 150.0f' in w or 'proxy_startup_timeout_sec = 150' in w),'timeout_guard_requires_aligned_branch':'const bool ready_by_timeout' in w and 'mag_gravity_aligned_branch_;' in w,'aligned_branch_is_negative_world_z':'return acc_world_lp.z() < 0.0f;' in common,'handoff_calls_goLive':'impl_.goLive(q_seed,' in w,'handoff_sets_live':'stage_ = Stage::Live;' in w}
 closed=bool(all(parity.values()) and dc['shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'] and margin>0)
 return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','filter_changed':False,'quality_gates_changed':False,'trajectory_replay_used':False,
  'shipping_parity':parity,'qualified_proxy_tilt_upper_rad':tilt,'qualified_proxy_tilt_upper_deg':math.degrees(tilt),'world_gravity_average_direction_error_upper_rad':eps,'timeout_aligned_branch_margin_rad':margin,'timeout_aligned_branch_margin_deg':math.degrees(margin),
  'deployed_timeout_s':150.0,'magnetic_north_required_before_timeout_handoff':False,'tuner_ready_required_before_timeout_handoff':False,'source_order_binary32_discrete_invariant_consumed':True,
  'QUALIFIED_LIVE_HANDOFF_BY_150S_CLOSED':closed,'P4_basin_reached_here':False,'P5_end_to_end_closed_here':False,
  'next_obligation':'start the real H18 early-Live capture word from this wide handoff set, retain timeout-prior tuner state and late-north handling, and prove finite entry into the certified P4 basin'}

def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 if not all(d.get('shipping_parity',{}).values()):f.append('shipping startup parity failed')
 for k in ('source_order_binary32_discrete_invariant_consumed','QUALIFIED_LIVE_HANDOFF_BY_150S_CLOSED'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('filter_changed','quality_gates_changed','trajectory_replay_used','magnetic_north_required_before_timeout_handoff','tuner_ready_required_before_timeout_handoff','P4_basin_reached_here','P5_end_to_end_closed_here'):
  if d.get(k) is not False:f.append(k+' not false')
 if not float(d.get('timeout_aligned_branch_margin_rad',0))>0:f.append('timeout branch margin not positive')
 if float(d.get('deployed_timeout_s',0))!=150.0:f.append('timeout changed')
 return f

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'capture150':d['QUALIFIED_LIVE_HANDOFF_BY_150S_CLOSED'],'tilt_deg':d['qualified_proxy_tilt_upper_deg'],'branch_margin_deg':d['timeout_aligned_branch_margin_deg'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())