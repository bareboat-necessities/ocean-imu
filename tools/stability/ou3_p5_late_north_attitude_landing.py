#!/usr/bin/env python3
"""Late-north yaw rewrite landing into the widest current P4 attitude sector.

The shipping rewrite is yaw-only and leaves tilt invariant. P5 targets a 19 deg
tilt core; the existing 10 deg heading-gauge bound then gives <29 deg full SO(3)
error by triangle inequality, strictly inside the 30 deg P4 sector. This does not
shrink P4 and closes only the attitude coordinate of the late-north landing.
"""
from __future__ import annotations
import json,math
from pathlib import Path
import ou3_p5_late_north_yaw_reset as RESET
REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_P5_LATE_NORTH_ATTITUDE_LANDING_V1'
CAPTURE_TILT_CORE_DEG=19.0

def build():
 d=json.loads(DOMAIN.read_text());r=RESET.build()
 if RESET.validate(r):raise RuntimeError('late-north reset prerequisite failed')
 widest=max(map(float,d['certificate_search']['p4_complete_word_full_attitude_candidate_deg']))
 yaw=math.degrees(float(d['startup']['internal_heading_gauge_error_upper_rad']))
 landing=CAPTURE_TILT_CORE_DEG+yaw
 structure=bool(r['LATE_NORTH_HYBRID_RESET_MAP_STRUCTURE_CLOSED'] and r['gravity_quotient_tilt_error_exactly_invariant'] and r['linear_navigation_state_exactly_unchanged'] and r['covariance_exactly_unchanged_by_setter'])
 closed=bool(structure and landing<widest)
 return {'qualification':QUALIFICATION,'filter_changed':False,'P4_attitude_basin_shrunk':False,'P4_widest_full_attitude_candidate_deg':widest,'P5_tilt_capture_core_deg':CAPTURE_TILT_CORE_DEG,'P5_tilt_core_is_not_P4_basin_radius':True,'declared_internal_heading_gauge_error_upper_deg':yaw,'SO3_triangle_inequality_consumed':True,'late_north_reset_structure_closed':structure,'post_reset_full_attitude_strict_upper_deg':landing,'strict_margin_to_widest_P4_attitude_sector_deg':widest-landing,'LATE_NORTH_ATTITUDE_LANDING_INTO_WIDEST_P4_SECTOR_CLOSED':closed,'other_P4_coordinates_landing_closed_here':False,'finite_time_reach_of_19deg_tilt_core_closed_here':False,'full_P4_membership_closed_here':False,'P4_PASS':False,'P5_PASS':False}
def validate(d):
 f=[]
 if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
 for k in ('P5_tilt_core_is_not_P4_basin_radius','SO3_triangle_inequality_consumed','late_north_reset_structure_closed','LATE_NORTH_ATTITUDE_LANDING_INTO_WIDEST_P4_SECTOR_CLOSED'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('filter_changed','P4_attitude_basin_shrunk','other_P4_coordinates_landing_closed_here','finite_time_reach_of_19deg_tilt_core_closed_here','full_P4_membership_closed_here','P4_PASS','P5_PASS'):
  if d.get(k) is not False:f.append(k+' not false')
 if not float(d.get('post_reset_full_attitude_strict_upper_deg',99))<30:f.append('landing outside P4')
 return f
if __name__=='__main__':
 d=build();f=validate(d);print(json.dumps({**d,'validation_pass':not f,'validation_failures':f},indent=2,sort_keys=True));raise SystemExit(bool(f))
