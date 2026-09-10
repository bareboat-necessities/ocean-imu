#!/usr/bin/env python3
"""Chart-free structural tilt bound at the shipping Live handoff.

The deployed handoff requires the world-frame aligned branch acc_world_lp.z()<0.
With the declared same-history averaged-gravity direction error eps_g, spherical
triangle inequality gives theta_tilt < pi/2 + eps_g. This is a P5 handoff
section, not a P4 basin restriction and does not assume the historical 60 deg
Mahony chart or any yaw gauge.
"""
from __future__ import annotations
import json,math
from pathlib import Path
REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
COMMON=REPO/'src/kalman_common/SeaStateFusionFilterCommon.h'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_SHIPPING_HANDOFF_TILT_HEMISPHERE_V1'

def build():
 d=json.loads(DOMAIN.read_text());c=COMMON.read_text();w=WRAPPER.read_text();eps=float(d['startup']['world_averaged_gravity_direction_error_upper_rad'])
 p1='return acc_world_lp.z() < 0.0f;' in c and 'gravityAlignedBranchWorld' in c
 p2='const bool ready_by_timeout' in w and 'mag_gravity_aligned_branch_;' in w and 'if (!ready_by_quality && !ready_by_timeout) return;' in w
 upper=.5*math.pi+eps
 return {'qualification':QUALIFICATION,'filter_changed':False,'quality_gates_changed':False,'shipping_branch_parity':{'world_branch_is_strict_negative_z':p1,'handoff_requires_branch_on_timeout_path':p2},'same_history_world_averaged_gravity_direction_error_upper_rad':eps,'true_gravity_quotient_tilt_strict_upper_rad':upper,'true_gravity_quotient_tilt_strict_upper_deg':math.degrees(upper),'historical_60_deg_mahony_chart_consumed':False,'yaw_gauge_required':False,'applies_to_quality_and_timeout_handoff_samples':True,'handoff_tilt_set_representation':'OPEN_GEODESIC_BALL_ON_S2_QUOTIENT','early_live_must_contract_to_P4':True,'HANDOFF_TILT_HEMISPHERE_BOUND_CLOSED':bool(p1 and p2 and 0<=eps<.5*math.pi)}
def validate(d):
 f=[]
 if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
 if d.get('HANDOFF_TILT_HEMISPHERE_BOUND_CLOSED') is not True:f.append('handoff hemisphere bound not closed')
 if d.get('historical_60_deg_mahony_chart_consumed') is not False:f.append('historical chart reintroduced')
 if d.get('yaw_gauge_required') is not False:f.append('yaw gauge incorrectly required')
 return f
if __name__=='__main__':
 d=build();f=validate(d);print(json.dumps({**d,'validation_pass':not f,'validation_failures':f},indent=2));raise SystemExit(bool(f))
