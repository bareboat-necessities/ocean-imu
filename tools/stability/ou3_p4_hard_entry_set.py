#!/usr/bin/env python3
"""Deterministic hard P4 entry set, independent of covariance consistency.

The entry set is a physical-coordinate box equal to a fixed fraction of the
already-declared startup/handoff coordinate envelope.  Membership is a hard
state-error hypothesis, not a probabilistic covariance ellipsoid and not a
replay-fit radius.  Reachability of this entry set from startup is a separate
hybrid/P5 obligation; P4 itself is the regional implication from this set.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'


def build(path: Path=CLOSURE):
    c=json.loads(Path(path).read_text())['hard_entry_search']
    scale=float(c['minimum_certified_scale'])
    b={k:float(v) for k,v in c['base_coordinate_radii'].items()}
    if not (0<scale<=1): raise RuntimeError('invalid hard entry scale')
    r={k:math.nextafter(scale*v,math.inf) for k,v in b.items()}
    return {
      'qualification':'OU3_P4_DETERMINISTIC_HARD_ENTRY_SET_V1',
      'source':'predeclared fraction of physical startup/handoff envelope',
      'trajectory_fit':False,'covariance_ellipsoid_used':False,
      'shipping_covariance_membership_used':False,
      'hard_state_error_membership_required':True,
      'scale_of_declared_handoff_envelope':scale,
      'coordinate_radii':r,
      'H18_groups':['attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m','integral_displacement_norm_m_s','latent_acceleration_norm_mps2'],
      'A21_additional_group':'accelerometer_bias_error_norm_mps2',
      'nontrivial_minimum_scale_enforced':scale>=0.01,
      'entry_set_declared_and_membership_checkable':True,
      'startup_reachability_proved_here':False,
      'P4_may_use_as_regional_entry_hypothesis':True,
    }


def validate(d):
    f=[]
    for k in ('hard_state_error_membership_required','nontrivial_minimum_scale_enforced','entry_set_declared_and_membership_checkable','P4_may_use_as_regional_entry_hypothesis'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('trajectory_fit','covariance_ellipsoid_used','shipping_covariance_membership_used','startup_reachability_proved_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('scale_of_declared_handoff_envelope',0))<0.01:f.append('hard entry scale below floor')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
