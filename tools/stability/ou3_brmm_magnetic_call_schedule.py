#!/usr/bin/env python3
"""Deterministic magnetometer call schedule and H18->A21 release timing.

The non-ALT BRMM route declares `hybrid_events_separate_from_P3_word` but has
no timing assumption that makes the shipping accelerometer-bias release
reachable in finite time. Without one, `H18` may persist forever and the
`A21` half of the required `H=18` / `A=21` pair is unreachable by assumption
rather than by proof.

`MAG-CALL-SCHEDULE-v1` is the deployment/source timing class

    first post-Live `updateMag()` call within 0.04 s of Live,
    every subsequent `updateMag()` call gap at most 0.04 s,

i.e. only 25 Hz. It is a call-cadence assumption on the host integration, kept
strictly separate from the magnetic *value* class `MAG-BMM150-DET-v1`.

Shipping `updateMag()` increments `mag_updates_applied_` immediately after the
MEKF magnetometer update and does not gate that counter on innovation
acceptance, so the release predicate

    accel_bias_locked_ && Live && mag_updates_applied_ >= 250 &&
    (time_ - first_mag_update_time_) > 1.0

is driven by call count and wall time alone.

## Two things the schedule alone does not give

**The counter is behind the outer north-lock gate.** The only call site of the
counter-owning inner method is

    if (mag_ref_set_ && stage_ == Stage::Live) {
        impl_.updateMag(mag_body_ned - mag_hard_iron_body_uT_);
    }

so a host call schedule does *not* advance `mag_updates_applied_` until the
startup magnetic reference is locked. On the explicitly admitted ungauged
timeout path `mag_ref_set_` stays false, and the release is then unreachable
however fast the host calls. The release time below is therefore conditional on
north lock, and north lock is exactly what
`ou3_brmm_magnetic_startup_yaw_capture` shows this route cannot supply. The two
obligations share one prerequisite.

**The count bound is an upper bound, so it cannot prove the strict guard.**
`(n-1)*gap_max` bounds the elapsed time from above; calls faster than 25 Hz
reach 250 sooner, and 250 calls 1 ms apart clear the count at `0.249 s` with the
shipping `> 1.0f` predicate still false. The release time is therefore derived
by a case split on whether the count or the guard binds last, never by comparing
an upper bound against the guard.

The external `acc_bias_hold_` branch can also keep
`set_acc_bias_updates_enabled(true)` from firing, so eventual A21 under an
arbitrary external hold is not claimed either.
"""
from __future__ import annotations
import argparse,json,re
from fractions import Fraction as F
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
SCHEMA=1
QUALIFICATION='OU3_BRMM_MAGNETIC_CALL_SCHEDULE_V1'
ASSUMPTION_ID='MAG-CALL-SCHEDULE-v1'

FIRST_CALL_AFTER_LIVE_MAX_S=F(1,25)
CALL_GAP_MAX_S=F(1,25)


def _one(pattern:str,text:str,label:str)->str:
 m=re.search(pattern,text,flags=re.MULTILINE)
 if not m:raise RuntimeError(f'cannot extract shipping constant {label}')
 return m.group(1)


def shipping_constants()->dict:
 w=WRAPPER.read_text(encoding='utf-8')
 return {
  'MAG_UPDATES_TO_UNLOCK':int(_one(r'static\s+constexpr\s+int\s+MAG_UPDATES_TO_UNLOCK\s*=\s*([0-9]+)\s*;',w,'MAG_UPDATES_TO_UNLOCK')),
  'acc_bias_unlock_mag_updates':int(_one(r'int\s+acc_bias_unlock_mag_updates\s*=\s*([0-9]+)\s*;',w,'acc_bias_unlock_mag_updates')),
  'MAG_DELAY_SEC':float(_one(r'constexpr\s+float\s+MAG_DELAY_SEC\s*=\s*([0-9.eE+-]+)f',w,'MAG_DELAY_SEC')),
  'strict_guard_sec':float(_one(r'first_mag_update_time_\)\s*>\s*([0-9.]+)f\)',w,'one-second guard')),
  'counter_increment_is_unconditional':'mekf_->measurement_update_mag_only(mag_body_ned);\n        mag_updates_applied_++;' in w,
  'release_requires_live_stage':'startup_stage_ == StartupStage::Live &&' in w,
  'release_requires_count':'mag_updates_applied_ >= mag_updates_to_unlock_ &&' in w,
  'external_hold_gates_enable':'if (!acc_bias_hold_) {' in w,
  # The counter-owning inner call is behind the outer north-lock gate; the exact
  # guarded call is checked, not merely the presence of its name.
  'counter_call_is_behind_north_lock_gate':(
   'if (mag_ref_set_ && stage_ == Stage::Live) {\n            impl_.updateMag(mag_body_ned - mag_hard_iron_body_uT_);' in w),
  'proxy_startup_timeout_sec':float(_one(r'float\s+proxy_startup_timeout_sec\s*=\s*([0-9.eE+-]+)f',w,'proxy_startup_timeout_sec')),
 }


def release_reachability(n:int,gap:F=CALL_GAP_MAX_S,first:F=FIRST_CALL_AFTER_LIVE_MAX_S,guard:F=F(1))->dict:
 """Finite release time from the schedule, by a case split on the binding term.

 `gap` bounds call spacing only from above, so `(n-1)*gap` is an upper bound on
 the elapsed time to the n-th call and cannot establish the strict
 `elapsed > guard` predicate. Instead:

 * the count condition `applied >= n` holds from the n-th call onward, and that
   call arrives no later than `first + (n-1)*gap`;
 * the guard condition `time - first_call > guard` holds from the first call
   strictly after `first_call + guard`, and with spacing at most `gap` such a
   call arrives no later than `first_call + guard + gap`.

 Both conditions are monotone once true, so the release fires no later than the
 maximum of the two, which is finite without any lower bound on spacing.
 """
 if n<1:raise ValueError('positive unlock count required')
 if gap<=0 or guard<0:raise ValueError('positive gap and nonnegative guard required')
 count_upper=first+(n-1)*gap              # Live -> n-th call, upper bound
 guard_upper=first+guard+gap              # Live -> first call past the guard
 total=max(count_upper,guard_upper)
 return {'unlock_count':n,'first_call_after_live_upper_s':first,
  'count_condition_satisfied_by_upper_s':count_upper,
  'guard_condition_satisfied_by_upper_s':guard_upper,
  'binding_condition':'count' if count_upper>=guard_upper else 'strict_guard',
  'live_to_unlock_upper_s':total,'strict_guard_s':guard,
  'release_derived_by_case_split_not_upper_bound_comparison':True,
  'lower_bound_on_call_spacing_assumed':False}


def build()->dict:
 c=shipping_constants();domain=json.loads(DOMAIN.read_text(encoding='utf-8'))
 n=c['acc_bias_unlock_mag_updates']
 if n!=c['MAG_UPDATES_TO_UNLOCK']:
  raise RuntimeError('shipping unlock count is not the wrapper default; schedule certificate needs the configured value')
 r=release_reachability(n)
 events=list(domain['normal_live']['hybrid_events_separate_from_P3_word'])
 parity_ok=all(v is True for k,v in c.items() if isinstance(v,bool))
 # Finite release time is established, but only from north lock onward: the
 # counter-owning call sits behind the outer mag_ref_set_ gate.
 release_after_north_lock=bool(parity_ok and r['release_derived_by_case_split_not_upper_bound_comparison'])
 return {
  'schema':SCHEMA,'qualification':QUALIFICATION,'assumption_id':ASSUMPTION_ID,
  'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','filter_changed':False,'quality_gates_changed':False,
  'trajectory_replay_used':False,'magnetic_value_class_used_here':False,
  'schedule':{'first_post_live_call_upper_s':float(FIRST_CALL_AFTER_LIVE_MAX_S),
   'call_gap_upper_s':float(CALL_GAP_MAX_S),'effective_rate_hz':float(1/CALL_GAP_MAX_S),
   'is_deployment_timing_assumption_not_sensor_guarantee':True},
  'shipping_constants':c,
  'shipping_parity':{k:v for k,v in c.items() if isinstance(v,bool)},
  'release':{'unlock_count':r['unlock_count'],
   'live_to_unlock_upper_s':float(r['live_to_unlock_upper_s']),
   'count_condition_satisfied_by_upper_s':float(r['count_condition_satisfied_by_upper_s']),
   'guard_condition_satisfied_by_upper_s':float(r['guard_condition_satisfied_by_upper_s']),
   'binding_condition':r['binding_condition'],
   'strict_one_second_guard_s':float(r['strict_guard_s']),
   'release_derived_by_case_split_not_upper_bound_comparison':r['release_derived_by_case_split_not_upper_bound_comparison'],
   'lower_bound_on_call_spacing_assumed':r['lower_bound_on_call_spacing_assumed'],
   'measured_from_north_lock_not_from_Live':True},
  # The internal MAG_DELAY_SEC gate is measured on the impl clock, which starts
  # at construction; both handoff paths are strictly later, so the delay never
  # postpones the first post-Live call.
  'mag_delay_precedes_live_on_timeout_path':c['MAG_DELAY_SEC']<c['proxy_startup_timeout_sec'],
  'H18_A21_release_time_after_north_lock_upper_s':float(r['live_to_unlock_upper_s']),
  'H18_A21_RELEASE_TIME_CLOSED_AFTER_NORTH_LOCK':release_after_north_lock,
  'H18_A21_RELEASE_TIME_CLOSED_FROM_LIVE':False,
  'counter_advances_without_north_lock':False,
  'north_lock_is_an_unmet_prerequisite_on_this_route':True,
  'ungauged_timeout_path_reaches_the_release':False,
  'hybrid_events_outside_P3_word':events,
  'innovation_acceptance_used_for_counter':False,
  'eventual_A21_under_arbitrary_external_acc_bias_hold_claimed':False,
  'H18_A21_covariance_transport_closed_here':False,
  'P4_promoted_here':False,'P5_promoted_here':False,
  'next_obligation':'supply north lock -- the same prerequisite the startup magnetic capture cannot meet on this route -- then attach this release time to the actual H18->A21 joint24 covariance transport edge and the reverse hold edge with BA cross-covariances zeroed',
 }


def validate(d:dict)->list:
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 if d.get('assumption_id')!=ASSUMPTION_ID:f.append('wrong call-schedule assumption id')
 if not all(d.get('shipping_parity',{}).values()):f.append('shipping magnetometer release parity failed')
 for k in ('H18_A21_RELEASE_TIME_CLOSED_AFTER_NORTH_LOCK','mag_delay_precedes_live_on_timeout_path',
           'north_lock_is_an_unmet_prerequisite_on_this_route'):
  if d.get(k) is not True:f.append(k+' not true')
 r=d.get('release',{})
 if r.get('release_derived_by_case_split_not_upper_bound_comparison') is not True:
  f.append('release time still derived by comparing an upper bound with the guard')
 if r.get('lower_bound_on_call_spacing_assumed') is not False:f.append('a call-spacing lower bound was assumed')
 if r.get('measured_from_north_lock_not_from_Live') is not True:f.append('release clock not anchored at north lock')
 for k in ('filter_changed','quality_gates_changed','trajectory_replay_used','magnetic_value_class_used_here',
           'innovation_acceptance_used_for_counter','eventual_A21_under_arbitrary_external_acc_bias_hold_claimed',
           'H18_A21_covariance_transport_closed_here','H18_A21_RELEASE_TIME_CLOSED_FROM_LIVE',
           'counter_advances_without_north_lock','ungauged_timeout_path_reaches_the_release',
           'P4_promoted_here','P5_promoted_here'):
  if d.get(k) is not False:f.append(k+' not false')
 if not 0.0<float(d.get('H18_A21_release_time_after_north_lock_upper_s',0.0))<=10.5:f.append('release time not finite and inside 10.5 s')
 if 'magnetic_lock' not in d.get('hybrid_events_outside_P3_word',[]):f.append('magnetic hybrid events detached from declared domain')
 return f


def main()->int:
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
 a.output.parent.mkdir(parents=True,exist_ok=True)
 a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'assumption':d['assumption_id'],
  'release_after_north_lock_s':d['H18_A21_release_time_after_north_lock_upper_s'],
  'binding':d['release']['binding_condition'],
  'closed_after_north_lock':d['H18_A21_RELEASE_TIME_CLOSED_AFTER_NORTH_LOCK'],
  'closed_from_live':d['H18_A21_RELEASE_TIME_CLOSED_FROM_LIVE'],'failures':f},sort_keys=True))
 return int(bool(f))


if __name__=='__main__':raise SystemExit(main())
