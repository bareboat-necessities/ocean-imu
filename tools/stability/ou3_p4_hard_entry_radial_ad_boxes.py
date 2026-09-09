#!/usr/bin/env python3
"""Continuous hard-entry radial cells -> outward AD state boxes for P4.

The theorem entry geometry is the deterministic product of 3-vector Euclidean
balls from ``ou3_p4_hard_entry_set``.  Nonlinear interval AD additionally needs
a rectangular state enclosure.  For a closed radial interval a=[a0,a1] subset
[0,1], every state t*x with t in a and ||x_g||<=r_g satisfies componentwise

    |(t x_g)_j| <= a1 r_g.

This module emits that rectangular *outer enclosure* while retaining the exact
ball premises separately.  The box is never claimed equal to the hard set and
never replaces the affine hard-ball IQCs in an S-procedure.  Recursive radial
partitioning uses the exact closed-interval lineage API and therefore covers
all scales without point sampling.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_brmm_event_lineage_cover as LINEAGE

SCHEMA=1
QUALIFICATION='OU3_P4_HARD_ENTRY_CONTINUOUS_RADIAL_AD_BOX_V1'
GROUP_ORDER=(
 'attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m',
 'integral_displacement_norm_m_s','latent_acceleration_norm_mps2','accelerometer_bias_error_norm_mps2')

def radial_state_box(mode:str,radial:Interval,entry:dict|None=None):
    if mode not in ('H','A'):raise ValueError('mode must be H/A')
    if not isinstance(radial,Interval) or radial.lo<0 or radial.hi>1:raise ValueError('radial interval must lie in [0,1]')
    e=ENTRY.build() if entry is None else entry
    ef=ENTRY.validate(e)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    names=GROUP_ORDER[:6] if mode=='H' else GROUP_ORDER
    out=[]
    for name in names:
        r=float(e['coordinate_radii'][name]);b=math.nextafter(radial.hi*r,math.inf)
        out.extend(Interval(-b,b) for _ in range(3))
    return out

def box_contains_radial_ball_components(mode,radial,box,entry):
    names=GROUP_ORDER[:6] if mode=='H' else GROUP_ORDER
    if len(box)!=3*len(names):return False
    for g,name in enumerate(names):
        b=radial.hi*float(entry['coordinate_radii'][name])
        for j in range(3):
            x=box[3*g+j]
            if not(x.lo<=-b and x.hi>=b):return False
    return True

def build():
    e=ENTRY.build();ef=ENTRY.validate(e)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    root=LINEAGE.root_radial_cell();parts=LINEAGE.partition(root,4)
    hboxes=[radial_state_box('H',p.radial_scale,e) for p in parts]
    aboxes=[radial_state_box('A',p.radial_scale,e) for p in parts]
    exact_cover=LINEAGE.exact_partition_cover(parts)
    contains=all(box_contains_radial_ball_components(m,p.radial_scale,b,e) for p,b in zip(parts,hboxes) for m in ('H',) ) and all(box_contains_radial_ball_components('A',p.radial_scale,b,e) for p,b in zip(parts,aboxes))
    nested=all(hboxes[i][j].abs_upper()<=hboxes[i+1][j].abs_upper() for i in range(len(hboxes)-1) for j in range(18))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'hard_entry_full_declared_scale_consumed':bool(e['full_declared_scale_enforced']),
      'continuous_radial_interval_not_point_sample':True,'recursive_partition_exactly_covers_unit_interval':exact_cover,
      'rectangular_box_is_outer_AD_enclosure_only':True,'rectangular_box_replaces_hard_ball_IQCs':False,
      'component_bound_formula':'abs(x_gj)<=radial_hi*r_g','all_partition_boxes_contain_corresponding_radial_ball_components':contains,
      'nested_component_boxes_with_increasing_radial_upper':nested,'H18_box_dimension':18,'A21_box_dimension':21,
      'radial_partition_count':len(parts),'radial_partition_bounds':[p.radial_scale.as_list() for p in parts],
      'covariance_membership_used':False,'trajectory_replay_used':False,'domain_shrunk':False,
      'production_event_state_reachability_proved_here':False,'production_every_prefix_hard_domain_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'use each radial_state_box only to outward-enclose nonlinear event Jacobians; retain the original hard-ball IQCs as premises and prove every-prefix retention/first-exit closure before calling the event state reachable'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('hard_entry_full_declared_scale_consumed','continuous_radial_interval_not_point_sample','recursive_partition_exactly_covers_unit_interval','rectangular_box_is_outer_AD_enclosure_only','all_partition_boxes_contain_corresponding_radial_ball_components','nested_component_boxes_with_increasing_radial_upper'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rectangular_box_replaces_hard_ball_IQCs','covariance_membership_used','trajectory_replay_used','domain_shrunk','production_event_state_reachability_proved_here','production_every_prefix_hard_domain_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('H18_box_dimension')!=18 or d.get('A21_box_dimension')!=21:f.append('state box dimensions changed')
    if int(d.get('radial_partition_count',0))<=0:f.append('empty radial partition')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'radial_cells':d['radial_partition_count'],'exact_cover':d['recursive_partition_exactly_covers_unit_interval'],'AD_outer_only':d['rectangular_box_is_outer_AD_enclosure_only'],'retention':d['production_every_prefix_hard_domain_retention_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
