#!/usr/bin/env python3
"""Shipping-bound scope certificate for the >70 deg hard attitude reinitialization.

This module deliberately distinguishes the wrapper's hybrid hard-reinit branch
from the ordinary MEKF Joseph update/injection/recentering machinery.

Out of theorem scope:
  * tilt > 70 deg held for 0.35 s with cooldown expired;
  * Live initialize_from_acc_preserve_yaw();
  * pre-Live initialize_from_acc()+enterCold_()+resetTrackingState_();
  * the associated 3 s cooldown/re-entry hybrid history.

Still in theorem scope:
  * every ordinary prediction and accepted Joseph update;
  * quaternion correction injection / MEKF error-coordinate recentering;
  * covariance/Joseph transport;
  * S=0 events with actual applied anisotropic R_S;
  * A21 accelerometer-bias radial projection;
  * startup handoff, H18->A21, and late magnetic yaw regauging.

The shipping implementation is unchanged.  P4/P5 therefore quantify over the
explicit no-hard-attitude-reinit execution family.  Once P4 prefix retention is
closed below the shipping 70 deg guard, the exclusion is automatically
invariant after P4 entry; startup/capture remains responsible for reaching P4
without taking this excluded branch.
"""
from __future__ import annotations
import argparse,json,re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
WRAPPER=ROOT/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
DOMAIN=ROOT/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_NO_HARD_ATTITUDE_REINIT_SCOPE_V1'

def _constant(text,name):
    m=re.search(rf'constexpr\s+float\s+{re.escape(name)}\s*=\s*([0-9.]+)f\s*;',text)
    if not m: raise RuntimeError('shipping hard-reinit constant not found: '+name)
    return float(m.group(1))

def build():
    w=WRAPPER.read_text();d=json.loads(DOMAIN.read_text());live=d['normal_live']
    deg=_constant(w,'TILT_RESET_DEG');hold=_constant(w,'TILT_RESET_HOLD_SEC');cool=_constant(w,'TILT_RESET_COOLDOWN_SEC')
    parity={
      'guard_uses_hold_and_cooldown':'tilt_over_limit_sec_ >= TILT_RESET_HOLD_SEC && tilt_reset_cooldown_sec_ <= 0.0f' in w,
      'live_hard_reinit_is_preserve_yaw':'mekf_->initialize_from_acc_preserve_yaw(acc_in);' in w,
      'startup_hard_reinit_is_acc_only':'mekf_->initialize_from_acc(acc_in);' in w,
      'startup_hard_reinit_returns_cold':'enterCold_();' in w and 'resetTrackingState_();' in w,
      'normal_live_domain_declares_no_hard_rewrite':live['hard_attitude_rewrite_inside_word'] is False,
      'tilt_reset_declared_separate_hybrid':'tilt_reset' in live['hybrid_events_separate_from_P3_word'],
    }
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'shipping_hard_reinit_threshold_deg':deg,'shipping_hard_reinit_hold_sec':hold,'shipping_hard_reinit_cooldown_sec':cool,
      'shipping_parity':parity,'filter_changed':False,
      'hard_attitude_reinitialization_branch_in_P4_scope':False,
      'hard_attitude_reinitialization_branch_in_P5_capture_scope':False,
      'ordinary_MEKF_Joseph_updates_in_scope':True,
      'ordinary_MEKF_quaternion_injection_and_error_recentering_in_scope':True,
      'ordinary_MEKF_covariance_transport_in_scope':True,
      'S_zero_actual_RS_events_in_scope':True,
      'A21_accelerometer_bias_projection_in_scope':True,
      'late_magnetic_yaw_regauging_in_startup_capture_scope':True,
      'startup_handoff_and_H18_to_A21_in_scope':True,
      'P4_must_prove_prefix_attitude_below_hard_reinit_threshold':True,
      'P5_must_reach_P4_without_excluded_hard_reinit':True,
      'scope_closed':all(parity.values()) and deg==70.0 and hold==0.35 and cool==3.0,
      'P4_promoted_here':False,'P5_promoted_here':False,
    }

def validate(x):
    f=[]
    if x.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if not all(x.get('shipping_parity',{}).values()):f.append('shipping hard-reinit parity failed')
    for k in ('ordinary_MEKF_Joseph_updates_in_scope','ordinary_MEKF_quaternion_injection_and_error_recentering_in_scope','ordinary_MEKF_covariance_transport_in_scope','S_zero_actual_RS_events_in_scope','A21_accelerometer_bias_projection_in_scope','late_magnetic_yaw_regauging_in_startup_capture_scope','startup_handoff_and_H18_to_A21_in_scope','P4_must_prove_prefix_attitude_below_hard_reinit_threshold','P5_must_reach_P4_without_excluded_hard_reinit','scope_closed'):
        if x.get(k) is not True:f.append(k+' not true')
    for k in ('filter_changed','hard_attitude_reinitialization_branch_in_P4_scope','hard_attitude_reinitialization_branch_in_P5_capture_scope','P4_promoted_here','P5_promoted_here'):
        if x.get(k) is not False:f.append(k+' not false')
    if x.get('shipping_hard_reinit_threshold_deg')!=70.0:f.append('hard reinit threshold changed')
    if x.get('shipping_hard_reinit_hold_sec')!=0.35:f.append('hard reinit hold changed')
    if x.get('shipping_hard_reinit_cooldown_sec')!=3.0:f.append('hard reinit cooldown changed')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();x=build();f=validate(x);x['validation_pass']=not f;x['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n');print(json.dumps({'scope_closed':x['scope_closed'],'threshold_deg':x['shipping_hard_reinit_threshold_deg'],'ordinary_mekf_in_scope':x['ordinary_MEKF_Joseph_updates_in_scope'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
