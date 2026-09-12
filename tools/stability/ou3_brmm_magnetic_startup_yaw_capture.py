#!/usr/bin/env python3
"""Startup magnetic yaw capture and the non-ALT fresh-attitude entrance.

The non-ALT BRMM startup certificate closes a *hemisphere* handoff: the shipping
`ready_by_timeout` branch needs only `mag_gravity_aligned_branch_`, so
`ou3_startup_timeout_capture` reports
`magnetic_north_required_before_timeout_handoff = False` and
`ou3_startup_handoff_tilt_hemisphere` reports `yaw_gauge_required = False`.
Nothing in the route bounded yaw at all, while
`initial_filter_entrance.attitude` declares a 45 deg *full* attitude error set
whose gauged branch bound the route never supplied.

This module composes the named magnetic classes into that gap and states the
result fail-closed. At the shipping handoff

    q_seed = boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_),

so the seed's tilt is the private observer's tilt and its yaw is the one-time
magnetic gauge. The geodesic triangle inequality gives

    ||log(q_true q_seed^-1)|| <= theta_yaw + theta_proxy_tilt,

and the *same* private-observer tilt also drives the accumulation frame that
produced the gauge, so both terms are functions of one supplied tilt bound.

Two supplies are compared:

* the operating domain's declared `world_averaged_gravity_direction_error`
  (0.02 rad). With MAG-BMM150-DET-v1 this yields `|sin(theta_yaw)| <= 17/30`,
  `theta_yaw < 0.61` rad by the alternating Taylor bound
  `sin(0.61) >= 0.61 - 0.61^3/6`, and a full attitude error `< 0.63` rad
  `= 36.0962 deg`, strictly inside the declared 45 deg entrance set.

* the route's own *certified* private-Mahony tilt bound, which is the level-set
  radius of the 87 deg-chart PI invariant. It is far too weak both for the
  magnetic gates and, on its own, for the 45 deg entrance set.

The declared 0.02 rad quantity is the low-passed world-gravity direction error,
not the private observer's tilt error, so it is recorded as an unsupplied
hypothesis rather than as a discharged one.
"""
from __future__ import annotations
import argparse,json,math
from fractions import Fraction as F
from pathlib import Path

import ou3_brmm_magnetic_source_envelope as ENV
import ou3_brmm_private_mahony_live_invariant as CONT
import ou3_brmm_gravity_direction_forcing_qualification as DIR

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
SCHEMA=1
QUALIFICATION='OU3_BRMM_MAGNETIC_STARTUP_YAW_CAPTURE_V1'
YAW_CERT_RAD=F(61,100)


def _taylor_sin_lower(x:F)->F:
 """sin(x) >= x - x^3/6 on [0,1], exact over the rationals."""
 if x<0 or x>1:raise ValueError('Taylor lower bound declared on [0,1]')
 return x-x*x*x/6


def yaw_certificate(sin_yaw_upper:F)->dict:
 x=YAW_CERT_RAD;lower=_taylor_sin_lower(x)
 return {'sin_yaw_error_upper':sin_yaw_upper,'yaw_test_angle_rad':x,
  'sin_test_angle_lower':lower,'yaw_strictly_below_test_angle':lower>sin_yaw_upper}


def _proxy_tilt_supply()->dict:
 c=CONT.build();failures=CONT.validate(c)
 return {'certified_tilt_rad_upper':float(c['actual_tilt_rad_upper']),
  'certified_tilt_deg_upper':float(c['actual_tilt_deg_upper']),
  'invariant_level_sqrt_C':float(c['sqrt_C']),'proof_chart_deg':float(c['chart_deg']),
  'invariant_closed_at_padded_envelope':bool(c['continuous_all_live_PI_invariant_closed']),
  'invariant_validation_failures':failures,
  'bound_is_level_set_radius_not_accuracy_bound':True}


def _forcing_scaling(required_rad:float)->dict:
 """Scaling argument for the best tilt this forcing qualification can support.

 The PI observer's own steady state under a constant admitted mean chord m is
 theta_ss = -m/s, and the transformed state carries the primitive term 0.1*xi
 directly into theta = z_theta - 0.1*xi. Neither term depends on the metric, so
 they bound from below what *any* storage/metric choice can certify.
 """
 q=DIR.build();m=float(q['mean_direction_chord_norm_upper']);xi=float(q['direction_primitive_norm_upper_s'])
 chart=math.radians(float(q['private_mahony_proof_chart_deg']));s_min=math.sin(chart)/chart
 dc=m/s_min;prim=0.1*xi;total=dc+prim
 return {'declared_mean_direction_chord_norm_upper':m,'declared_direction_primitive_norm_upper_s':xi,
  'endpoint_sector_s_lower':s_min,'metric_free_mean_chord_tilt_floor_rad':dc,
  'metric_free_primitive_tilt_floor_rad':prim,'metric_free_tilt_floor_rad':total,
  'metric_free_tilt_floor_deg':math.degrees(total),'required_tilt_rad':required_rad,
  'shortfall_factor':total/required_rad if required_rad>0 else math.inf,
  'metric_refinement_alone_cannot_reach_requirement':total>required_rad,
  'is_scaling_argument_not_certificate':True}


def build()->dict:
 domain=json.loads(DOMAIN.read_text(encoding='utf-8'));w=WRAPPER.read_text(encoding='utf-8')
 env=ENV.build();ef=ENV.validate(env)
 if ef:raise RuntimeError('magnetic source envelope invalid: '+repr(ef))
 entrance=domain['initial_filter_entrance']['attitude']
 declared_entrance_deg=float(entrance['full_attitude_error_upper_deg'])
 declared_tilt=F(str(domain['startup']['world_averaged_gravity_direction_error_upper_rad'])).limit_denominator(10**6)
 c=ENV.shipping_constants()

 # Declared-supply branch, exact over the rationals.
 declared=ENV.derive(declared_tilt/2,c)
 yaw=yaw_certificate(declared['sin_yaw_error_upper'])
 full=YAW_CERT_RAD+declared_tilt
 declared_branch={
  'supplied_tilt_rad_upper':float(declared_tilt),'supplied_tilt_deg_upper':math.degrees(float(declared_tilt)),
  'mean_perturbation_upper_uT':float(declared['mean_perturbation_upper_uT']),
  'sin_yaw_error_upper_exact':f"{declared['sin_yaw_error_upper'].numerator}/{declared['sin_yaw_error_upper'].denominator}",
  'sin_test_angle_lower':float(yaw['sin_test_angle_lower']),
  'yaw_rad_upper':float(YAW_CERT_RAD),'yaw_deg_upper':math.degrees(float(YAW_CERT_RAD)),
  'full_attitude_rad_upper':float(full),'full_attitude_deg_upper':math.degrees(float(full)),
  'north_nonvanishing':declared['north_nonvanishing'],
  'horizontal_fraction_gate_satisfied':declared['horizontal_fraction_gate_satisfied'],
  'inside_declared_entrance_set':yaw['yaw_strictly_below_test_angle'] and math.degrees(float(full))<declared_entrance_deg,
  'sqrtN_statistical_reduction_used':False}

 # Certified-supply branch.
 proxy=_proxy_tilt_supply()
 req=env['derived_tilt_requirements']
 binding=float(req['max_tilt_rad_for_horizontal_fraction_gate'])
 certified_branch={
  **proxy,
  'required_tilt_rad_for_binding_gate':binding,
  'required_tilt_deg_for_binding_gate':float(req['max_tilt_deg_for_horizontal_fraction_gate']),
  'required_tilt_deg_for_north_capture':float(req['max_tilt_deg_for_north_capture']),
  'binding_gate':req['binding_requirement'],
  'gate_shortfall_factor':proxy['certified_tilt_rad_upper']/binding,
  'north_capture_shortfall_factor':proxy['certified_tilt_rad_upper']/float(req['max_tilt_rad_for_north_capture']),
  'entrance_shortfall_factor':proxy['certified_tilt_deg_upper']/declared_entrance_deg,
  'magnetic_capture_forced_from_certified_supply':proxy['certified_tilt_rad_upper']<binding,
  'certified_tilt_alone_inside_declared_entrance_set':proxy['certified_tilt_deg_upper']<declared_entrance_deg,
  'perfect_yaw_gauge_would_repair_entrance':False}

 parity={'timeout_branch_needs_only_aligned_branch':'mag_gravity_aligned_branch_;' in w and 'const bool ready_by_timeout' in w,
  'quality_branch_needs_north_ready':'const bool north_ready = !cfg_.with_mag || mag_ref_set_;' in w,
  'seed_composes_proxy_tilt_with_gauge_yaw':'boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_)' in w,
  'ungauged_seed_uses_free_yaw_sigma':'cfg_.proxy_handoff_yaw_sigma_free_rad' in w}

 scaling=_forcing_scaling(binding)
 return {
  'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
  'magnetic_value_assumption':env['assumption_id'],'filter_changed':False,'quality_gates_changed':False,
  'trajectory_replay_used':False,'shipping_parity':parity,
  'declared_entrance_full_attitude_deg':declared_entrance_deg,
  'declared_entrance_gauged_branch_bound':bool(entrance['gauged_branch_full_attitude_bound']),
  'declared_entrance_ungauged_branch_is_tilt_only':bool(entrance['ungauged_branch_tilt_bound_only_until_regauge']),
  'declared_supply_branch':declared_branch,
  'certified_supply_branch':certified_branch,
  'forcing_scaling_argument':scaling,
  'DECLARED_SUPPLY_ENTRANCE_CLOSED':bool(declared_branch['inside_declared_entrance_set']),
  'CERTIFIED_SUPPLY_ENTRANCE_CLOSED':False,
  'accumulation_frame_tilt_supply_discharged':False,
  'timeout_branch_yaw_gauge_forced':False,
  'P4_basin_reached_here':False,'P5_end_to_end_closed_here':False,
  'limiting_quantity':'private-observer accumulation tilt-frame error bound',
  'alternatives':[
   'tighten the gravity-direction forcing qualification (mean chord 0.05, primitive 1.0 s) so a metric-free ultimate tilt below the binding gate exists at all',
   'narrow MAG-BMM150-DET-v1 to a commissioned geomagnetic band (higher horizontal floor, lower total ceiling) so the shipping gates are forced at the tilt this route can certify',
   'prove that the shipping quality handoff cannot be preceded by the timeout handoff under the declared schedule, making north_ready a theorem consequence rather than an assumption'],
  'next_obligation':'supply a private-observer accumulation tilt-frame accuracy bound, not a level-set radius, or requalify one of the three recorded alternatives',
 }


def validate(d:dict)->list:
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 if not all(d.get('shipping_parity',{}).values()):f.append('shipping handoff parity failed')
 if d.get('DECLARED_SUPPLY_ENTRANCE_CLOSED') is not True:f.append('declared-supply entrance certificate not closed')
 for k in ('CERTIFIED_SUPPLY_ENTRANCE_CLOSED','accumulation_frame_tilt_supply_discharged',
           'timeout_branch_yaw_gauge_forced','filter_changed','quality_gates_changed','trajectory_replay_used',
           'P4_basin_reached_here','P5_end_to_end_closed_here'):
  if d.get(k) is not False:f.append(k+' not false')
 b=d.get('declared_supply_branch',{})
 if b.get('sqrtN_statistical_reduction_used') is not False:f.append('statistical averaging gain used')
 if not 0.0<float(b.get('full_attitude_deg_upper',180.0))<float(d.get('declared_entrance_full_attitude_deg',0.0)):
  f.append('declared-supply full attitude bound not inside the declared entrance set')
 c=d.get('certified_supply_branch',{})
 if c.get('magnetic_capture_forced_from_certified_supply') is not False:
  f.append('certified-supply capture claim must stay fail-closed until the tilt supply is proved')
 if not float(c.get('gate_shortfall_factor',0.0))>1.0:f.append('certified-supply shortfall factor not reported')
 s=d.get('forcing_scaling_argument',{})
 if s.get('is_scaling_argument_not_certificate') is not True:f.append('scaling argument mislabelled as a certificate')
 if len(d.get('alternatives',[]))<3:f.append('two-strike rule requires at least three qualitatively different alternatives')
 return f


def main()->int:
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
 a.output.parent.mkdir(parents=True,exist_ok=True)
 a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'declared_supply_full_attitude_deg':d['declared_supply_branch']['full_attitude_deg_upper'],
  'declared_entrance_deg':d['declared_entrance_full_attitude_deg'],
  'declared_supply_closed':d['DECLARED_SUPPLY_ENTRANCE_CLOSED'],
  'certified_tilt_deg':d['certified_supply_branch']['certified_tilt_deg_upper'],
  'required_tilt_deg':d['certified_supply_branch']['required_tilt_deg_for_binding_gate'],
  'gate_shortfall':d['certified_supply_branch']['gate_shortfall_factor'],
  'metric_free_floor_deg':d['forcing_scaling_argument']['metric_free_tilt_floor_deg'],
  'failures':f},sort_keys=True))
 return int(bool(f))


if __name__=='__main__':raise SystemExit(main())
