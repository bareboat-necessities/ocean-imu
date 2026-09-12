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

* the route's own *certified* private-Mahony tilt bounds from the two-phase
  certificate: the all-time outer level, the bound at the deployed 150 s
  startup horizon, and the asymptotic bound. All three are far too weak for the
  magnetic gates and, on their own, for the 45 deg entrance set.

The declared 0.02 rad quantity is the low-passed world-gravity direction error,
not the private observer's tilt error, so it is recorded as an unsupplied
hypothesis rather than as a discharged one.

## Why the private observer can never supply it

The obstruction is not the choice of metric. The declared forcing qualification
admits `r = m + d xi/dt` with `||m||<=0.05` and `||xi||<=1.0 s`, and the
certificate quantifies over the sector value `s` as an independent parameter in
`[sinc(87 deg),1]`. For the member `s=1` the tilt transfer function is

    theta/r = -(0.01 + 0.1 j w) / (0.01 s - w^2 + 0.1 s j w),

so an admitted primitive-bounded sinusoid `r = j w xi` produces
`|theta| = w |theta/r| |xi|`. That response peaks at `w = 0.117 rad/s`
(`0.0186 Hz`), which lies *inside* the declared physical band
`[0.018,0.88] Hz`, giving `|theta| >= 0.146781` rad from the primitive channel
alone; superposing an admitted DC `m` adds `0.05` rad. Hence no certificate in
this formulation -- any metric, any level, any number of subdivisions -- can
bound the tilt below

    0.196781 rad = 11.2748 deg,

against the `6.1145 deg` north-capture and `2.8378 deg` horizontal-fraction
requirements. This is a property of the declared forcing pair and the deployed
gains, so the limiter is `(mean chord 0.05, primitive 1.0 s)` and the magnetic
envelope, not the storage geometry.
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
  'certified_tilt_at_deployed_horizon_deg_upper':float(c['tilt_at_deployed_horizon_deg_upper']),
  'certified_ultimate_tilt_deg_upper':float(c['ultimate_tilt_deg_upper']),
  'invariant_level_sqrt_C':float(c['sqrt_C']),'proof_chart_deg':float(c['chart_deg']),
  'invariant_closed_at_padded_envelope':bool(c['continuous_all_live_PI_invariant_closed']),
  'invariant_validation_failures':failures,
  'two_phase_certificate_consumed':True,
  'accumulation_frame_may_start_before_the_horizon':True,
  'sound_supply_for_accumulation_frame_is_all_time_bound':True}


PEAK_OMEGA_RAD_S=F(117,1000)


def _sqrt_lower(x:F,digits:int=30)->F:
 """Rational lower bound for sqrt(x), x>=0, by integer square root."""
 if x<0:raise ValueError('nonnegative argument required')
 scale=10**digits
 n=(x*scale*scale).numerator//(x*scale*scale).denominator if x else 0
 return F(math.isqrt(n),scale)


def _declared_formulation_tilt_floor(required_rad:float)->dict:
 """Lower bound on any tilt the DECLARED formulation can certify.

 The certificate quantifies over the sector value s as an independent parameter
 in [sinc(87 deg),1], so its conclusion must also hold for the member s=1. For
 that member the deployed PI gains give

     theta/r = -(0.01 + 0.1 j w)/(0.01 - w^2 + 0.1 j w),

 and an admitted primitive-bounded sinusoid r = j w xi produces
 |theta| = w |theta/r| |xi|. Evaluated at the declared rational peak frequency
 with exact rational arithmetic, plus a superposed admitted DC mean chord. This
 is a statement about what the declared formulation can yield, not a claim
 about the physical observer.
 """
 q=DIR.build();m=F(str(q['mean_direction_chord_norm_upper'])).limit_denominator(10**6)
 xi=F(str(q['direction_primitive_norm_upper_s'])).limit_denominator(10**6)
 w=PEAK_OMEGA_RAD_S;w2=w*w
 num_sq=F(1,10000)+F(1,100)*w2                      # |0.01 + 0.1 j w|^2
 den_sq=(F(1,100)-w2)**2+F(1,100)*w2                # |0.01 - w^2 + 0.1 j w|^2
 prim=w*_sqrt_lower(num_sq/den_sq)*xi
 total=prim+m
 f=float(total)
 # The peak must lie inside the declared COMPLETE-BRMM frequency support, or the
 # floor would rest on a frequency the physical envelope does not admit.
 band=[float(x) for x in json.loads(DOMAIN.read_text(encoding='utf-8'))['complete_brmm_physical_envelope']['frequency_support_hz']]
 peak_hz=float(w)/(2*math.pi)
 in_band=band[0]<=peak_hz<=band[1]
 return {'declared_mean_direction_chord_norm_upper':float(m),
  'declared_direction_primitive_norm_upper_s':float(xi),
  'peak_omega_rad_s':float(w),'peak_frequency_hz':peak_hz,
  'declared_wave_band_hz':band,
  'peak_frequency_inside_declared_wave_band':in_band,
  'primitive_channel_tilt_floor_rad':float(prim),
  'mean_chord_tilt_floor_rad':float(m),
  'declared_formulation_tilt_floor_rad':f,
  'declared_formulation_tilt_floor_deg':math.degrees(f),
  'required_tilt_rad':required_rad,
  'shortfall_factor':f/required_rad if required_rad>0 else math.inf,
  'no_metric_or_level_can_reach_requirement':f>required_rad,
  'frozen_sector_member_s_equals_one':True,
  'is_formulation_lower_bound_not_physical_claim':True}


def _required_envelope_narrowing(tilt_rad:float,c:dict)->dict:
 """Commissioned band that would force the shipping gates at a supplied tilt.

 The horizontal-fraction gate needs E <= (H - f Bmax)/(1+f) with
 E = P + 2 Bmax sin(delta/2), so for a supplied delta the requirement is a
 linear lower bound on the horizontal floor H as a function of Bmax.
 """
 f=F(str(c['min_horizontal_fraction'])).limit_denominator(10**6)
 P=ENV.HARD_IRON_NORM_MAX_UT+ENV.RESIDUAL_NORM_MAX_UT
 half=math.sin(tilt_rad/2.0)
 rows=[]
 for bmax in (45,55,65,75):
  B=F(bmax)
  need=(1+f)*(P+2*B*F(half).limit_denominator(10**9))+f*B
  rows.append({'world_field_norm_upper_uT':bmax,'required_horizontal_floor_uT':float(need),
   'feasible_within_total_field':float(need)<=bmax})
 return {'supplied_tilt_rad':tilt_rad,'supplied_tilt_deg':math.degrees(tilt_rad),
  'sin_half_supplied_tilt':half,'declared_combined_perturbation_uT':float(P),
  'shipping_min_horizontal_fraction':float(f),'rows':rows,
  'declared_envelope_row':{'world_field_norm_upper_uT':float(ENV.WORLD_FIELD_NORM_MAX_UT),
   'declared_horizontal_floor_uT':float(ENV.WORLD_FIELD_HORIZONTAL_MIN_UT)}}


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
 north=float(req['max_tilt_rad_for_north_capture'])
 certified_branch={
  **proxy,
  'required_tilt_rad_for_binding_gate':binding,
  'required_tilt_deg_for_binding_gate':float(req['max_tilt_deg_for_horizontal_fraction_gate']),
  'required_tilt_deg_for_north_capture':float(req['max_tilt_deg_for_north_capture']),
  'binding_gate':req['binding_requirement'],
  'gate_shortfall_factor':proxy['certified_tilt_rad_upper']/binding,
  'north_capture_shortfall_factor':proxy['certified_tilt_rad_upper']/north,
  'ultimate_gate_shortfall_factor':math.radians(proxy['certified_ultimate_tilt_deg_upper'])/binding,
  'ultimate_north_capture_shortfall_factor':math.radians(proxy['certified_ultimate_tilt_deg_upper'])/north,
  'entrance_shortfall_factor':proxy['certified_tilt_deg_upper']/declared_entrance_deg,
  'magnetic_capture_forced_from_certified_supply':proxy['certified_tilt_rad_upper']<binding,
  'magnetic_capture_forced_from_ultimate_supply':math.radians(proxy['certified_ultimate_tilt_deg_upper'])<binding,
  'certified_tilt_alone_inside_declared_entrance_set':proxy['certified_tilt_deg_upper']<declared_entrance_deg,
  'perfect_yaw_gauge_would_repair_entrance':False}

 parity={'timeout_branch_needs_only_aligned_branch':'mag_gravity_aligned_branch_;' in w and 'const bool ready_by_timeout' in w,
  'quality_branch_needs_north_ready':'const bool north_ready = !cfg_.with_mag || mag_ref_set_;' in w,
  'seed_composes_proxy_tilt_with_gauge_yaw':'boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_)' in w,
  'ungauged_seed_uses_free_yaw_sigma':'cfg_.proxy_handoff_yaw_sigma_free_rad' in w}

 floor=_declared_formulation_tilt_floor(binding)
 narrowing=_required_envelope_narrowing(floor['declared_formulation_tilt_floor_rad'],c)
 return {
  'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
  'magnetic_value_assumption':env['assumption_id'],'filter_changed':False,'quality_gates_changed':False,
  'trajectory_replay_used':False,'shipping_parity':parity,
  'declared_entrance_full_attitude_deg':declared_entrance_deg,
  'declared_entrance_gauged_branch_bound':bool(entrance['gauged_branch_full_attitude_bound']),
  'declared_entrance_ungauged_branch_is_tilt_only':bool(entrance['ungauged_branch_tilt_bound_only_until_regauge']),
  'declared_supply_branch':declared_branch,
  'certified_supply_branch':certified_branch,
  'declared_formulation_tilt_floor':floor,
  'required_envelope_narrowing_at_the_floor':narrowing,
  'DECLARED_SUPPLY_ENTRANCE_CLOSED':bool(declared_branch['inside_declared_entrance_set']),
  'CERTIFIED_SUPPLY_ENTRANCE_CLOSED':False,
  'accumulation_frame_tilt_supply_discharged':False,
  'timeout_branch_yaw_gauge_forced':False,
  'P4_basin_reached_here':False,'P5_end_to_end_closed_here':False,
  'limiting_quantity':'declared gravity-direction forcing pair (mean chord, primitive) and the commissioned magnetic band',
  'private_observer_can_supply_the_declared_gates':False,
  'quadratic_metric_class_ruled_out':True,
  'alternatives':[
   'tighten the gravity-direction forcing qualification: the primitive channel alone contributes w|theta/r||xi| at the in-band peak, so the declared (0.05, 1.0 s) pair must come down together and must be argued from the COMPLETE-BRMM spectral support',
   'narrow MAG-BMM150-DET-v1 to a commissioned geomagnetic band: the gates depend on Bmax and H_min, not tilt alone, and required_envelope_narrowing_at_the_floor gives the horizontal floor needed at each total-field ceiling',
   'replace the static quadratic storage with a frequency-dependent multiplier/IQC on the primitive channel: the obstruction is a frequency-response peak, which a static metric cannot see',
   'prove that the shipping quality handoff cannot be preceded by the timeout handoff under the declared schedule, making north_ready a theorem consequence rather than an assumption'],
  'next_obligation':'the private observer cannot supply the declared gates in this formulation; close either the forcing qualification or the commissioned band, and use an IQC rather than a static metric for the tilt certificate',
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
 cb=d.get('certified_supply_branch',{})
 if cb.get('magnetic_capture_forced_from_certified_supply') is not False:
  f.append('certified-supply capture claim must stay fail-closed until the tilt supply is proved')
 if not float(cb.get('gate_shortfall_factor',0.0))>1.0:f.append('certified-supply shortfall factor not reported')
 s=d.get('declared_formulation_tilt_floor',{})
 if s.get('is_formulation_lower_bound_not_physical_claim') is not True:f.append('formulation floor mislabelled as a physical claim')
 if s.get('no_metric_or_level_can_reach_requirement') is not True:f.append('formulation floor no longer dominates the requirement')
 if s.get('peak_frequency_inside_declared_wave_band') is not True:f.append('peak primitive frequency outside the declared wave band')
 for k in ('private_observer_can_supply_the_declared_gates','quadratic_metric_class_ruled_out'):
  if k=='quadratic_metric_class_ruled_out':
   if d.get(k) is not True:f.append(k+' not true')
  elif d.get(k) is not False:f.append(k+' not false')
 c=d.get('certified_supply_branch',{})
 if c.get('magnetic_capture_forced_from_ultimate_supply') is not False:
  f.append('ultimate-supply capture claim must stay fail-closed')
 if c.get('two_phase_certificate_consumed') is not True:f.append('two-phase certificate not consumed')
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
  'formulation_floor_deg':d['declared_formulation_tilt_floor']['declared_formulation_tilt_floor_deg'],
  'failures':f},sort_keys=True))
 return int(bool(f))


if __name__=='__main__':raise SystemExit(main())
