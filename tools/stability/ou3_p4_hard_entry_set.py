#!/usr/bin/env python3
"""Legacy conditional product and the exact fresh-Live source entry graph.

The legacy radii remain available to reproduce earlier non-promoting diagnostics
and as candidate working-domain scales. They are not startup reachability facts.
Fresh entry uses the shipping initialization and shared physical-S origin graph;
its uniform reachable attitude/covariance/tuner cover is still an obligation.
No covariance-confidence ellipsoid is used as a true-error set.
"""
from __future__ import annotations
import argparse, json, math
import ou3_p4_live_entry_graph as LIVE_ENTRY
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'


def build(path: Path=CLOSURE):
    c=json.loads(Path(path).read_text())['hard_entry_search']
    scale=float(c['minimum_certified_scale'])
    b={k:float(v) for k,v in c['base_coordinate_radii'].items()}
    if scale != 1.0 or list(map(float,c['candidate_scale_factors'])) != [1.0]:
        raise RuntimeError('P4 hard entry set may not shrink the declared operating domain')
    r={k:math.nextafter(v,math.inf) for k,v in b.items()}
    return {
      'qualification':'OU3_P4_DETERMINISTIC_HARD_ENTRY_SET_V2',
      'source':'legacy conditional diagnostic/working product; not fresh-Live reachability',
      'legacy_box_is_certified_fresh_live_entry':False,
      'fresh_live_entry_graph':LIVE_ENTRY.build(),
      'S_radius_role':'legacy conditional diagnostic only; not independent fresh-entry uncertainty',
      'trajectory_fit':False,'covariance_ellipsoid_used':False,
      'shipping_covariance_membership_used':False,
      'hard_state_error_membership_required':True,
      'scale_of_declared_handoff_envelope':scale,
      'coordinate_radii':r,
      'H18_groups':['attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m','integral_displacement_norm_m_s','latent_acceleration_norm_mps2'],
      'A21_additional_group':'accelerometer_bias_error_norm_mps2',
      'full_declared_scale_enforced':True,
      'entry_set_declared_and_membership_checkable':True,
      'startup_reachability_proved_here':False,
      'P4_may_use_as_regional_entry_hypothesis':True,
      'proof_domain_shrunk':False,
    }


def validate(d):
    f=[]
    for k in ('hard_state_error_membership_required','full_declared_scale_enforced','entry_set_declared_and_membership_checkable','P4_may_use_as_regional_entry_hypothesis'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('trajectory_fit','covariance_ellipsoid_used','shipping_covariance_membership_used','startup_reachability_proved_here','proof_domain_shrunk'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('scale_of_declared_handoff_envelope',0)) != 1.0:f.append('hard entry scale is not full declared envelope')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
