#!/usr/bin/env python3
"""Named deterministic magnetic source envelope for the non-ALT BRMM route.

`MAG-BMM150-DET-v1` is the commissioned-installation admission class

    20 uT <= ||B_W||_2 <= 75 uT,
    ||(B_Wx,B_Wy)||_2 >= 15 uT,
    ||b_HI||_2 <= 5 uT,
    ||n_m||_2 <= 2 uT per theorem sample,

for the exact startup magnetic identity

    m_body = R_true B_W + b_HI + n_m.

It is an engineering source requirement on a commissioned installation, not a
sensor datasheet guarantee, and `Rmag` is never read as a deterministic bound.

The non-ALT route had no magnetic source class at all: its magnetic PE numbers
(`normal_live.magnetic_vector_norm_lower_uT = 10`, upper 200) were bare declared
hypotheses. This module attaches the named class to the actual shipping startup
acquisition path and derives, exactly over the rationals:

1. the measured body-field norm band, hence whether the declared PE magnetic
   floor/ceiling are consequences of the class rather than extra hypotheses;
2. whether the class clears the shipping `mag_init_min_mag_norm` guard;
3. the combined hard-iron + residual budget the shipping
   `MagAutoTuner::max_sample_norm_ratio_from_mean` gate needs before *every*
   qualified sample is admitted;
4. the accumulation tilt-frame error the shipping
   `MagAutoTuner::min_horizontal_fraction` gate and the non-vanishing-north
   condition need before the startup reference/yaw gauge can be forced.

Items 3 and 4 are derived *requirements*, not assumptions: they are reported
whether or not the declared class and the certified startup tilt satisfy them.
"""
from __future__ import annotations
import argparse,json,math,re
from fractions import Fraction as F
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
TUNER=REPO/'src/tuner/MagAutoTuner.h'
SCHEMA=1
QUALIFICATION='OU3_BRMM_MAGNETIC_SOURCE_ENVELOPE_V1'
ASSUMPTION_ID='MAG-BMM150-DET-v1'

# MAG-BMM150-DET-v1, exactly as the named class declares it.
WORLD_FIELD_NORM_MIN_UT=F(20)
WORLD_FIELD_NORM_MAX_UT=F(75)
WORLD_FIELD_HORIZONTAL_MIN_UT=F(15)
HARD_IRON_NORM_MAX_UT=F(5)
RESIDUAL_NORM_MAX_UT=F(2)


def _one(pattern:str,text:str,label:str)->str:
 m=re.search(pattern,text,flags=re.MULTILINE)
 if not m:raise RuntimeError(f'cannot extract shipping constant {label}')
 return m.group(1)


def shipping_constants()->dict:
 w=WRAPPER.read_text(encoding='utf-8');t=TUNER.read_text(encoding='utf-8')
 return {
  'mag_init_min_mag_norm_uT':float(_one(r'float\s+mag_init_min_mag_norm\s*=\s*([0-9.eE+-]+)f',w,'mag_init_min_mag_norm')),
  'mag_min_samples':int(_one(r'int\s+mag_min_samples\s*=\s*([0-9]+)\s*;',w,'mag_min_samples')),
  'mag_min_window_sec':float(_one(r'float\s+mag_min_window_sec\s*=\s*([0-9.eE+-]+)f',w,'mag_min_window_sec')),
  'mag_max_window_sec':float(_one(r'float\s+mag_max_window_sec\s*=\s*([0-9.eE+-]+)f',w,'mag_max_window_sec')),
  'mag_tilt_fallback_sec':float(_one(r'float\s+mag_tilt_fallback_sec\s*=\s*([0-9.eE+-]+)f',w,'mag_tilt_fallback_sec')),
  'proxy_mag_settle_sec':float(_one(r'float\s+proxy_mag_settle_sec\s*=\s*([0-9.eE+-]+)f',w,'proxy_mag_settle_sec')),
  'proxy_startup_timeout_sec':float(_one(r'float\s+proxy_startup_timeout_sec\s*=\s*([0-9.eE+-]+)f',w,'proxy_startup_timeout_sec')),
  'mag_gravity_align_hold_sec':float(_one(r'float\s+mag_gravity_align_hold_sec\s*=\s*([0-9.eE+-]+)f',w,'mag_gravity_align_hold_sec')),
  'mag_enable_quality_weighting':_one(r'bool\s+mag_enable_quality_weighting\s*=\s*(true|false)\s*;',w,'mag_enable_quality_weighting')=='true',
  'mag_estimate_hard_iron':_one(r'bool\s+mag_estimate_hard_iron\s*=\s*(true|false)\s*;',w,'mag_estimate_hard_iron')=='true',
  'max_sample_norm_ratio_from_mean':float(_one(r'float\s+max_sample_norm_ratio_from_mean\s*=\s*([0-9.eE+-]+)f',t,'max_sample_norm_ratio_from_mean')),
  'min_horizontal_fraction':float(_one(r'float\s+min_horizontal_fraction\s*=\s*([0-9.eE+-]+)f',t,'min_horizontal_fraction')),
  'tilt_frame_is_private_observer_tilt':'tiltOnlyQuatFromBoatQuat_(attitudeReferenceQuat_())' in w,
  'world_sample_is_rotated_body_sample':'const Eigen::Vector3f mag_world_i = q * mag_body_ned;' in t,
  'reference_gauge_fixed_to_horizontal_x':'mag_world_ref_ = Eigen::Vector3f(horiz, 0.0f, mean.z());' in t,
 }


def _sin_half_upper_from_angle(theta:F)->F:
 """Rational upper bound for sin(theta/2) using sin x <= x."""
 if theta<0:raise ValueError('nonnegative tilt angle required')
 return theta/2


def derive(sin_half_tilt_upper:F,c:dict)->dict:
 """Deterministic magnetic derivations at a supplied sin(delta_tilt/2) bound."""
 s=F(sin_half_tilt_upper)
 if s<0 or s>1:raise ValueError('sin(delta/2) bound must lie in [0,1]')
 P=HARD_IRON_NORM_MAX_UT+RESIDUAL_NORM_MAX_UT            # combined deterministic body perturbation
 rot=2*WORLD_FIELD_NORM_MAX_UT*s                          # ||(Rhat Rtrue' - I) B|| chord bound
 E=P+rot                                                  # accumulation-mean perturbation bound
 H=WORLD_FIELD_HORIZONTAL_MIN_UT
 f=F(c['min_horizontal_fraction']).limit_denominator(10**6)
 ratio=F(c['max_sample_norm_ratio_from_mean']).limit_denominator(10**6)
 # Universal per-sample admission at the norm-ratio gate:
 #   every ||m_body|| and the running mean lie in [Bmin-P, Bmax+P]; a sufficient
 #   condition is 2P/(Bmin-P) <= ratio, i.e. P <= ratio*Bmin/(2+ratio).
 P_budget=ratio*WORLD_FIELD_NORM_MIN_UT/(2+ratio)
 # min_horizontal_fraction gate on the finalized mean:
 #   (H-E)/(Bmax+E) >= f  <=>  E <= (H - f*Bmax)/(1+f).
 E_horiz=(H-f*WORLD_FIELD_NORM_MAX_UT)/(1+f)
 s_horiz=(E_horiz-P)/(2*WORLD_FIELD_NORM_MAX_UT) if E_horiz>P else F(0)
 # Non-vanishing north condition: E < H.
 s_capture=(H-P)/(2*WORLD_FIELD_NORM_MAX_UT) if H>P else F(0)
 return {
  'sin_half_tilt_upper':s,'combined_body_perturbation_upper_uT':P,
  'earth_field_rotation_chord_upper_uT':rot,'mean_perturbation_upper_uT':E,
  'horizontal_true_lower_uT':H,
  'measured_body_norm_lower_uT':WORLD_FIELD_NORM_MIN_UT-P,
  'measured_body_norm_upper_uT':WORLD_FIELD_NORM_MAX_UT+P,
  'mean_horizontal_lower_uT':H-E,'mean_norm_upper_uT':WORLD_FIELD_NORM_MAX_UT+E,
  'horizontal_fraction_lower':(H-E)/(WORLD_FIELD_NORM_MAX_UT+E) if E<H else F(0),
  'sin_yaw_error_upper':E/H,
  'north_nonvanishing':E<H,
  'horizontal_fraction_gate_satisfied':E<H and (H-E)/(WORLD_FIELD_NORM_MAX_UT+E)>=f,
  'universal_sample_admission_perturbation_budget_uT':P_budget,
  'universal_sample_admission_forced':P<=P_budget,
  'max_sin_half_tilt_for_horizontal_gate':s_horiz,
  'max_sin_half_tilt_for_north_capture':s_capture,
 }


def _asin_upper(x:F)->float:
 return math.asin(float(x)) if 0<=float(x)<=1 else float('nan')


def build()->dict:
 d=json.loads(DOMAIN.read_text(encoding='utf-8'));live=d['normal_live'];c=shipping_constants()
 declared_tilt=F(str(d['startup']['world_averaged_gravity_direction_error_upper_rad'])).limit_denominator(10**6)
 # Unconditional part: no tilt bound is used.
 P=HARD_IRON_NORM_MAX_UT+RESIDUAL_NORM_MAX_UT
 measured_lo=WORLD_FIELD_NORM_MIN_UT-P;measured_hi=WORLD_FIELD_NORM_MAX_UT+P
 pe_floor=F(str(live['magnetic_vector_norm_lower_uT'])).limit_denominator(10**6)
 pe_ceiling=F(str(live['magnetic_vector_norm_upper_uT'])).limit_denominator(10**6)
 guard=F(str(c['mag_init_min_mag_norm_uT'])).limit_denominator(10**9)
 # Tilt-parameterised part, evaluated at the declared startup direction error.
 declared=derive(_sin_half_upper_from_angle(declared_tilt),c)
 thresholds=derive(F(0),c)
 s_horiz=thresholds['max_sin_half_tilt_for_horizontal_gate']
 s_capture=thresholds['max_sin_half_tilt_for_north_capture']
 tilt_horiz_rad=2*_asin_upper(s_horiz);tilt_capture_rad=2*_asin_upper(s_capture)
 return {
  'schema':SCHEMA,'qualification':QUALIFICATION,'assumption_id':ASSUMPTION_ID,
  'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','filter_changed':False,'quality_gates_changed':False,
  'trajectory_replay_used':False,'Rmag_used_as_deterministic_bound':False,
  'source_identity':'m_body = R_true B_W + b_HI + n_m',
  'envelope':{'world_field_norm_lower_uT':float(WORLD_FIELD_NORM_MIN_UT),'world_field_norm_upper_uT':float(WORLD_FIELD_NORM_MAX_UT),
   'world_field_horizontal_lower_uT':float(WORLD_FIELD_HORIZONTAL_MIN_UT),'hard_iron_norm_upper_uT':float(HARD_IRON_NORM_MAX_UT),
   'residual_norm_upper_uT':float(RESIDUAL_NORM_MAX_UT),'combined_body_perturbation_upper_uT':float(P)},
  'shipping_constants':c,
  'shipping_parity':{'accumulation_tilt_frame_is_private_observer':c['tilt_frame_is_private_observer_tilt'],
   'world_sample_is_norm_preserving_rotation':c['world_sample_is_rotated_body_sample'],
   'reference_is_gauge_fixed_to_horizontal_x':c['reference_gauge_fixed_to_horizontal_x'],
   'startup_quality_weighting_disabled':c['mag_enable_quality_weighting'] is False,
   'startup_hard_iron_solve_disabled':c['mag_estimate_hard_iron'] is False},
  'measured_body_field_norm_lower_uT':float(measured_lo),'measured_body_field_norm_upper_uT':float(measured_hi),
  'declared_PE_magnetic_floor_uT':float(pe_floor),'declared_PE_magnetic_ceiling_uT':float(pe_ceiling),
  'declared_PE_magnetic_floor_is_now_derived':measured_lo>=pe_floor,
  'declared_PE_magnetic_ceiling_is_now_derived':measured_hi<=pe_ceiling,
  'shipping_mag_norm_guard_cleared_unconditionally':measured_lo>guard,
  'norm_ratio_gate':{'shipping_ratio':c['max_sample_norm_ratio_from_mean'],
   'universal_admission_perturbation_budget_uT':float(thresholds['universal_sample_admission_perturbation_budget_uT']),
   'declared_combined_perturbation_uT':float(P),
   'universal_sample_admission_forced':thresholds['universal_sample_admission_forced'],
   'note':'sufficient condition 2P/(Bmin-P) <= ratio; a weaker necessary condition needs an accepted-sample supply assumption that is not declared'},
  'derived_tilt_requirements':{
   'max_sin_half_tilt_for_north_capture':float(s_capture),'max_tilt_rad_for_north_capture':tilt_capture_rad,
   'max_tilt_deg_for_north_capture':math.degrees(tilt_capture_rad),
   'max_sin_half_tilt_for_horizontal_fraction_gate':float(s_horiz),
   'max_tilt_rad_for_horizontal_fraction_gate':tilt_horiz_rad,
   'max_tilt_deg_for_horizontal_fraction_gate':math.degrees(tilt_horiz_rad),
   'binding_requirement':'min_horizontal_fraction' if s_horiz<s_capture else 'north_nonvanishing'},
  'at_declared_startup_direction_error':{
   'declared_tilt_rad_upper':float(declared_tilt),'declared_tilt_deg_upper':math.degrees(float(declared_tilt)),
   'mean_perturbation_upper_uT':float(declared['mean_perturbation_upper_uT']),
   'sin_yaw_error_upper':float(declared['sin_yaw_error_upper']),
   'sin_yaw_error_upper_exact':f"{declared['sin_yaw_error_upper'].numerator}/{declared['sin_yaw_error_upper'].denominator}",
   'horizontal_fraction_lower':float(declared['horizontal_fraction_lower']),
   'north_nonvanishing':declared['north_nonvanishing'],
   'horizontal_fraction_gate_satisfied':declared['horizontal_fraction_gate_satisfied']},
  'declared_startup_direction_error_is_the_accumulation_frame_error':False,
  'accumulation_frame_error_supply_is_open_obligation':True,
  'magnetic_vector_sine_separation_derivable_from_this_class':False,
  'sine_separation_remains_declared_PE_hypothesis':True,
  'P4_promoted_here':False,'P5_promoted_here':False,
  'next_obligation':'supply the actual startup accumulation tilt-frame error for the private observer below the derived binding requirement, or widen MAG-BMM150-DET-v1 so the shipping gates are forced at the tilt error the non-ALT route can certify',
 }


def validate(d:dict)->list:
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 if d.get('assumption_id')!=ASSUMPTION_ID:f.append('wrong magnetic assumption id')
 if not all(d.get('shipping_parity',{}).values()):f.append('shipping magnetic acquisition parity failed')
 for k in ('declared_PE_magnetic_floor_is_now_derived','declared_PE_magnetic_ceiling_is_now_derived',
           'shipping_mag_norm_guard_cleared_unconditionally','accumulation_frame_error_supply_is_open_obligation',
           'sine_separation_remains_declared_PE_hypothesis'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('filter_changed','quality_gates_changed','trajectory_replay_used','Rmag_used_as_deterministic_bound',
           'declared_startup_direction_error_is_the_accumulation_frame_error',
           'magnetic_vector_sine_separation_derivable_from_this_class','P4_promoted_here','P5_promoted_here'):
  if d.get(k) is not False:f.append(k+' not false')
 t=d.get('derived_tilt_requirements',{})
 if not 0.0<float(t.get('max_tilt_rad_for_horizontal_fraction_gate',0.0))<float(t.get('max_tilt_rad_for_north_capture',0.0)):
  f.append('derived tilt requirement ordering invalid')
 a=d.get('at_declared_startup_direction_error',{})
 if a.get('north_nonvanishing') is not True or a.get('horizontal_fraction_gate_satisfied') is not True:
  f.append('declared startup direction error does not clear the shipping magnetic gates')
 return f


def main()->int:
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
 a.output.parent.mkdir(parents=True,exist_ok=True)
 a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'assumption':d['assumption_id'],
  'measured_norm_uT':[d['measured_body_field_norm_lower_uT'],d['measured_body_field_norm_upper_uT']],
  'pe_floor_derived':d['declared_PE_magnetic_floor_is_now_derived'],
  'universal_admission_forced':d['norm_ratio_gate']['universal_sample_admission_forced'],
  'tilt_required_deg':d['derived_tilt_requirements']['max_tilt_deg_for_horizontal_fraction_gate'],
  'failures':f},sort_keys=True))
 return int(bool(f))


if __name__=='__main__':raise SystemExit(main())
