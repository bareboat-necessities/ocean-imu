#!/usr/bin/env python3
"""Same-history low-frequency qualification for BRMM gravity-direction forcing.

The padded COMPLETE-BRMM family permits ||a_m||<=8.8 m/s^2 instantaneously.
A pointwise bound alone permits adversarial rectification of the normalized
specific-force direction; no PI gravity observer can uniformly reject an
arbitrary persistent false-gravity direction. Marine wave forcing is
oscillatory, so the theorem keeps an explicit low-frequency qualification on
the same physical direction history consumed by the private Mahony observer.

Let r(t) be the normalized world-gravity direction correction disturbance.
Require one same-history decomposition

    r(t) = m(t) + d xi(t)/dt,
    ||m(t)|| <= 0.05,
    ||xi(t)|| <= 1.0 s.

This does not reduce the instantaneous padded acceleration cap. The witness
xi/m are source coordinates carried with the same BRMM history; they are not
independently selected per sample.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
SCHEMA=2
QUALIFICATION='OU3_BRMM_GRAVITY_DIRECTION_FORCING_V2'
MEAN_CHORD_NORM_UPPER=0.05
PRIMITIVE_NORM_UPPER_S=1.0
PROOF_CHART_DEG=87.0
EXPECTED_PADDED_ACCEL_MPS2=8.8


def build(domain_path:Path=DOMAIN):
 d=json.loads(Path(domain_path).read_text());g=float(d['startup']['gravity_mps2']);a=float(d['normal_live']['non_gravitational_cog_acceleration_norm_upper_mps2'])
 if not 0<a<g:raise RuntimeError('gravity-direction qualification requires 0<a<g')
 env=d.get('complete_brmm_physical_envelope',{})
 if float(env.get('wave_acceleration_norm_upper_mps2',math.nan))!=a:raise RuntimeError('gravity-direction cap detached from COMPLETE-BRMM physical envelope')
 rho=a/g;alpha=math.asin(rho);chord=2*math.sin(alpha/2)
 return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
  'same_history_required':True,'independent_sample_direction_boxes_used':False,'trajectory_replay_used':False,
  'instantaneous_acceleration_norm_upper_mps2':a,'instantaneous_direction_angle_upper_rad':alpha,
  'instantaneous_direction_angle_upper_deg':math.degrees(alpha),'instantaneous_direction_chord_upper':chord,
  'direction_forcing_decomposition':'r=m+d(xi)/dt','mean_direction_chord_norm_upper':MEAN_CHORD_NORM_UPPER,
  'direction_primitive_norm_upper_s':PRIMITIVE_NORM_UPPER_S,'private_mahony_proof_chart_deg':PROOF_CHART_DEG,
  'pointwise_padded_acceleration_cap_retained':a==EXPECTED_PADDED_ACCEL_MPS2,
  'padded_acceleration_cap_mps2':EXPECTED_PADDED_ACCEL_MPS2,
  'low_frequency_false_gravity_is_explicitly_bounded':True,
  'physical_interpretation':'instantaneous orbital acceleration may be large; persistent/rectified normalized false-gravity direction must remain small and oscillatory direction has a bounded primitive',
  'source_applicability_refinement':True,'P4_promoted_here':False,'P5_promoted_here':False}


def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 for k in ('same_history_required','pointwise_padded_acceleration_cap_retained','low_frequency_false_gravity_is_explicitly_bounded','source_applicability_refinement'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('independent_sample_direction_boxes_used','trajectory_replay_used','P4_promoted_here','P5_promoted_here'):
  if d.get(k) is not False:f.append(k+' not false')
 if float(d.get('instantaneous_acceleration_norm_upper_mps2',math.nan))!=EXPECTED_PADDED_ACCEL_MPS2:f.append('padded acceleration cap mismatch')
 if not 0<float(d.get('mean_direction_chord_norm_upper',1))<=0.05:f.append('mean direction chord bound invalid')
 if not 0<float(d.get('direction_primitive_norm_upper_s',0))<=1.0:f.append('direction primitive bound invalid')
 if not 63<float(d.get('instantaneous_direction_angle_upper_deg',0))<65:f.append('8.8 mps2 direction angle mismatch')
 if float(d.get('private_mahony_proof_chart_deg',0))!=87.0:f.append('proof chart changed')
 return f


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'accel_cap':d['instantaneous_acceleration_norm_upper_mps2'],'angle_deg':d['instantaneous_direction_angle_upper_deg'],'chord':d['instantaneous_direction_chord_upper'],'mean':d['mean_direction_chord_norm_upper'],'primitive_s':d['direction_primitive_norm_upper_s'],'chart_deg':d['private_mahony_proof_chart_deg'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
