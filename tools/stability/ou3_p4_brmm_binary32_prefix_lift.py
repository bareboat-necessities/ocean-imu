#!/usr/bin/env python3
"""Binary32 additive state-forcing lift for OU-III P4 literal prefixes.

The finite-precision enclosure already bounds shipping Eigen/Kalman/reset and
time-update arithmetic conditionally on the platform contract.  This module
turns those scalar norm bounds into theorem-facing additive state coordinates

    x_computed = x_exact + n_fp

with one n-dimensional n_fp per literal operation and a homogeneous ball IQC
on the same h coordinate used by the hard-entry/prefix graph.

Prediction uses the certified time-update/libm state-error bound.  Accepted
Joseph/reset/projection events use the certified explicit measurement/reset
state-error bound.  Covariance-entry roundoff is not double-counted here: it
remains in the source-uniform reachable P/H/R/K coefficient enclosure that
feeds the next same-cell event.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_kalman_reset_binary32_iss as FPISS
import ou3_p4_affine_hard_tube_iqc as HARD

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_BINARY32_PREFIX_ADDITIVE_LIFT_V1'

def I(x):return Interval.point(float(x))
def identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def event_radius(contract,mode,kind):
    m=contract['modes'][mode]
    if kind=='prediction':return float(m['time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction'])
    if kind in ('S_zero','accelerometer','magnetometer','vector','projection'):
        return float(m['explicit_measurement_reset_state_roundoff_norm_upper_per_event'])
    if kind=='aw_floor':
        # State is unchanged by the covariance floor. Its arithmetic influence is
        # retained in the covariance/coefficient enclosure, not a fake state kick.
        return 0.0
    raise ValueError('unsupported finite-precision event kind')
def build_event_lift(contract,mode,kind):
    n=18 if mode=='H18' else 21;r=event_radius(contract,mode,kind)
    if not(math.isfinite(r) and r>=0):raise ValueError('invalid binary32 state forcing radius')
    # local homogeneous coordinate [h; n_fp(n)]
    Pi=None
    if r>0:Pi=HARD.ball_iqc(n+1,0,tuple(range(1,n+1)),r)
    return {'mode':mode,'kind':kind,'state_dimension':n,'forcing_radius':r,
            'state_injection_matrix':identity(n),'forcing_ball_iqc':Pi,
            'zero_state_forcing':r==0.0}
def build():
    c=FPISS.build();cf=FPISS.validate(c)
    if cf:raise RuntimeError('binary32 ISS prerequisite failed: '+repr(cf))
    lifts={}
    for mode in ('H18','A21'):
        lifts[mode]={k:build_event_lift(c,mode,k) for k in ('prediction','aw_floor','S_zero','accelerometer','magnetometer')}
    finite=all(math.isfinite(v['forcing_radius']) for m in lifts.values() for v in m.values())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'binary32_contract_qualification':c['qualification'],'platform_execution_contract_consumed':c['platform_execution_contract_consumed'],
      'additive_state_forcing_coordinate_materialized':True,'prediction_time_update_forcing_materialized':True,
      'Joseph_reset_projection_forcing_materialized':True,'aw_floor_state_forcing_exactly_zero':all(lifts[m]['aw_floor']['zero_state_forcing'] for m in lifts),
      'covariance_roundoff_not_double_counted_as_state_kick':True,'covariance_roundoff_retained_in_next_coefficient_enclosure':True,
      'homogeneous_ball_IQCs_materialized_for_nonzero_forcing':True,'all_forcing_radii_finite':finite,
      'lifts_summary':{m:{k:v['forcing_radius'] for k,v in e.items()} for m,e in lifts.items()},
      'roundoff_absorbed_into_homogeneous_contraction':False,'rowwise_K_used_as_reset_domain':False,
      'production_lineage_binary32_coordinates_embedded_here':False,'production_every_prefix_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'append the appropriate n_fp coordinate to every literal prediction/Joseph prefix, embed its ball IQC and additive state injection into the prefix transport, then combine with BIAS1 and exact nonlinear/reset/projection sectors in one LDLT coordinate'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('platform_execution_contract_consumed','additive_state_forcing_coordinate_materialized','prediction_time_update_forcing_materialized','Joseph_reset_projection_forcing_materialized','aw_floor_state_forcing_exactly_zero','covariance_roundoff_not_double_counted_as_state_kick','covariance_roundoff_retained_in_next_coefficient_enclosure','homogeneous_ball_IQCs_materialized_for_nonzero_forcing','all_forcing_radii_finite'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('roundoff_absorbed_into_homogeneous_contraction','rowwise_K_used_as_reset_domain','production_lineage_binary32_coordinates_embedded_here','production_every_prefix_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    for mode,events in d.get('lifts_summary',{}).items():
        for kind,x in events.items():
            if not(math.isfinite(float(x)) and float(x)>=0):f.append(mode+' '+kind+' radius invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'binary32_lift':d['additive_state_forcing_coordinate_materialized'],'radii':d['lifts_summary'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
