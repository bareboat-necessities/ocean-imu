#!/usr/bin/env python3
"""Discrete 5 ms closure of the qualified BRMM private-Mahony invariant.

The former comparison treated the gravity-direction correction as one arbitrary
pointwise chord and therefore encoded the old 4 m/s^2 source.  The widened
8 m/s^2 theorem uses the transformed same-history decomposition r=m+d xi/dt.
The production source-order binary32 certificate already includes both the
forward-Euler h^2 term and quaternion/integral roundoff, so this module is now
the compatibility/theorem interface consumed by the connected frontend proof.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_brmm_private_mahony_live_invariant as LIVE
import ou3_brmm_private_mahony_discrete_invariant as B32

SCHEMA=2
QUALIFICATION='OU3_BRMM_PRIVATE_MAHONY_DISCRETE_PI_COMPARISON_V2'

def build():
    live=LIVE.build();lf=LIVE.validate(live);b=B32.build();bf=B32.validate(b)
    if lf or bf:raise RuntimeError(f'discrete Mahony prerequisites failed live={lf} binary32={bf}')
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','source_generator':False,
      'same_BRMM_forcing_as_continuous_invariant':True,'same_history_direction_primitive_retained':True,'sample_period_s':b['dt_s'],
      'continuous_invariant_consumed':True,'old_independent_pointwise_chord_model_used':False,
      'continuous_normalized_boundary_margin_lower':live['boundary_validation']['strict_inward_margin_lower'],
      'continuous_boundary_margin_after_binary32':b['continuous_boundary_margin_after_binary32'],
      'h2_metric_increase_upper':b['forward_euler_quadratic_V_upper'],'first_order_step_decrease_lower':b['guaranteed_first_order_V_decrease_lower'],
      'discrete_metric_decrease_lower':b['discrete_V_margin_lower'],'ideal_5ms_discrete_PI_invariant_closed':b['discrete_V_margin_lower']>0,
      'shipping_binary32_quaternion_map_error_composed':b['shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'],
      'shipping_source_order_binary32_discrete_invariant_closed':b['shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'],
      'toolchain_independent_binary32_invariant_closed':b['toolchain_independent_binary32_invariant_closed'],
      'complete_BRMM_family_materialized_here':False,'P3_promoted':False,
      'next_obligation':'propagate the connected source-order binary32 Mahony state through WPE/tuner and separately discharge compiler/FMA realization qualification'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_BRMM_forcing_as_continuous_invariant','same_history_direction_primitive_retained','continuous_invariant_consumed','ideal_5ms_discrete_PI_invariant_closed','shipping_binary32_quaternion_map_error_composed','shipping_source_order_binary32_discrete_invariant_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_generator','old_independent_pointwise_chord_model_used','toolchain_independent_binary32_invariant_closed','complete_BRMM_family_materialized_here','P3_promoted'):
        if d.get(k) is not False:f.append(k+' not false')
    for k in ('continuous_boundary_margin_after_binary32','discrete_metric_decrease_lower'):
        if not float(d.get(k,0))>0:f.append(k+' not positive')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'ideal_closed':d['ideal_5ms_discrete_PI_invariant_closed'],'binary32_closed':d['shipping_source_order_binary32_discrete_invariant_closed'],'metric_margin':d['discrete_metric_decrease_lower'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
