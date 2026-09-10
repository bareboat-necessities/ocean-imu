#!/usr/bin/env python3
"""Minimal finite magnetic-progress admission for OU-III P5.

The timeout path can enter Live without north. Full-yaw capture therefore needs
finite progress of the literal MagAutoTuner accepted state, not merely packet
recurrence. This adds only the observability premise needed by shipping:
accepted_count and accepted_window_sec plus the tuner's own finite/horizontal
mean checks. It does not alter BRMM motion caps, BIAS families, filter code or
quality gates.
"""
from __future__ import annotations
import json,math,re
from pathlib import Path
REPO=Path(__file__).resolve().parents[2]
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
MAG=REPO/'src/tuner/MagAutoTuner.h'
QUALIFICATION='OU3_P5_MAGNETIC_PROGRESS_ADMISSION_V1'

def _num(text,name):
 m=re.search(rf'(?:float|int)\s+{re.escape(name)}\s*=\s*([0-9.eE+\-/ ]+)f?\s*;',text)
 if not m:raise RuntimeError('cannot extract '+name)
 raw=m.group(1).strip()
 if '/' in raw:
  a,b=raw.split('/',1);return float(a)/float(b)
 return float(raw)

def build(progress_horizon_s=None):
 w=WRAPPER.read_text();m=MAG.read_text();timeout=_num(w,'proxy_startup_timeout_sec');n=int(_num(w,'mag_min_samples'));win=_num(w,'mag_min_window_sec');ew=_num(w,'mag_min_effective_weight')
 horizon=timeout if progress_horizon_s is None else float(progress_horizon_s)
 if not(math.isfinite(horizon) and horizon>0):raise ValueError('positive finite horizon required')
 parity={
  'timeout_can_enter_live_without_north':'const bool ready_by_timeout' in w and 'mag_gravity_aligned_branch_;' in w,
  'postlive_first_north_keeps_accumulating':'if (!mag_ref_set_)' in w and 'mag_auto_tuner_.addSampleWithTiltQuatDt(' in w and 'impl_.mekf().set_quaternion_boat(q_new);' in w,
  'accepted_count_gate':'accepted_count_ < std::max(1, cfg_.min_samples)' in m,
  'accepted_window_gate':'accepted_window_sec_ < cfg_.min_window_sec' in m,
  'accepted_window_advances_only_after_acceptance':'accepted_window_sec_ += dt_use;' in m and m.index('accepted_window_sec_ += dt_use;')>m.index('if (w < cfg_.min_sample_weight)'),
  'mean_finite_and_nonzero':'if (!mean.allFinite())' in m and 'mean_norm > cfg_.mag_norm_min' in m,
  'horizontal_fraction_gate':'horiz_frac < cfg_.min_horizontal_fraction' in m,
 }
 default_ok=all(parity.values()) and ew==0.0 and 'bool  mag_estimate_hard_iron = false;' in w and 'bool  mag_enable_quality_weighting = false;' in w
 return {'qualification':QUALIFICATION,'filter_changed':False,'quality_gates_changed':False,'BRMM_motion_caps_changed':False,'BIAS_family_changed':False,'canonical_motion_source_remains_COMPLETE_BRMM':True,'shipping':{'startup_timeout_s':timeout,'min_accepted_samples':n,'min_accepted_window_s':win,'parity':parity},'admission':{'with_mag_false_requires_progress':False,'with_mag_true_progress_horizon_s':horizon,'same_history_required':True,'accepted_count_lower':n,'accepted_window_sec_lower':win,'accepted_mean_finite':True,'accepted_mean_norm_strictly_above_mag_norm_min':True,'accepted_mean_horizontal_fraction_at_least_shipping_min':True,'packet_recurrence_alone_is_not_substitute':True},'FINITE_NORTH_EVENT_UNDER_ADMISSION_CLOSED':default_ok,'finite_north_event_time_upper_from_progress_origin_s':horizon if default_ok else None,'global_physical_deployment_left_inclusion_claimed':False,'P4_PASS':False,'P5_PASS':False}
def validate(d):
 f=[]
 if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
 for k in ('canonical_motion_source_remains_COMPLETE_BRMM','FINITE_NORTH_EVENT_UNDER_ADMISSION_CLOSED'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('filter_changed','quality_gates_changed','BRMM_motion_caps_changed','BIAS_family_changed','global_physical_deployment_left_inclusion_claimed','P4_PASS','P5_PASS'):
  if d.get(k) is not False:f.append(k+' not false')
 s=d.get('shipping',{})
 if s.get('min_accepted_samples')!=128:f.append('min samples changed')
 if s.get('min_accepted_window_s')!=15.0:f.append('min window changed')
 if s.get('startup_timeout_s')!=150.0:f.append('timeout changed')
 if not all(s.get('parity',{}).values()):f.append('shipping magnetic parity failed')
 return f
if __name__=='__main__':
 d=build();f=validate(d);print(json.dumps({**d,'validation_pass':not f,'validation_failures':f},indent=2,sort_keys=True));raise SystemExit(bool(f))
