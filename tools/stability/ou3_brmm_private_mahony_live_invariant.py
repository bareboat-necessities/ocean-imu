#!/usr/bin/env python3
"""Validated all-Live PI invariant for the private Mahony observer.

The widened BRMM source permits large instantaneous orbital acceleration, so the
old proof that treated the normalized gravity-direction error as an arbitrary
pointwise disturbance inside 60 degrees is no longer valid. The qualified
marine source instead carries the same-history decomposition

    r = m + d xi/dt,   ||m||<=0.05,   ||xi||<=1.0 s.

For x=(theta,beta), deployed gains give

    xdot = A_s x + B r + d,
    A_s=[[-0.1 s,1],[-0.01 s,0]], B=[-0.1,-0.01]^T.

Set z=x-B xi. Then

    zdot=A_s z + A_s B xi + B m + d,

so the oscillatory primitive is retained as one bounded same-history source
coordinate rather than replaced by an arbitrary r box. The existing exact
quadratic metric P=R^T R is used at sqrt(C)=1.21. A rational stereographic
cover of the boundary, outward interval arithmetic, and the endpoint sector
s in [sinc(87 deg),1] prove strict inward flow. The actual tilt is recovered
as theta=z_theta-0.1 xi and remains strictly below the 87 degree proof chart.

The first-sample accelerometer seed angle is the source's own instantaneous
gravity-direction angle `asin(||a_ng||/g)`, taken from the qualification module
rather than frozen as a constant. At the padded 8.8 m/s^2 envelope that seed is
1.1137 rad, and this certificate reports fail-closed: the level needed to
contain the seed and the largest level that still fits inside the 87 degree
chart form an empty window, so no level of this fixed metric closes the
invariant for the padded family. See `admissible_level_window` in the output.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval
import ou3_validated_transcendentals as VT
import ou3_brmm_gravity_direction_forcing_qualification as DIR

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
SCHEMA=2
QUALIFICATION='OU3_BRMM_PRIVATE_MAHONY_ALL_LIVE_PI_INVARIANT_V2'
P11=1.0;P12=-6.5;P22=163.25;DET_P=121.0
SQRT_C=1.21;C_LEVEL=SQRT_C*SQRT_C
CELLS_PER_CHART=8192

def I(x):return Interval.outward_bounds(float(x),float(x))
def abs_upper(x):return max(abs(x.lo),abs(x.hi))

def _source_constants(domain):
 s=domain['startup'];live=domain['normal_live'];q=DIR.build();f=DIR.validate(q)
 if f:raise RuntimeError('gravity-direction source qualification invalid: '+repr(f))
 return {'g':float(s['gravity_mps2']),'a_ng':float(live['non_gravitational_cog_acceleration_norm_upper_mps2']),
  'gyro_transport':float(s['effective_deterministic_gyro_transport_disturbance_upper_rad_s']),
  'bias_transport':float(s['effective_deterministic_bias_transport_disturbance_upper_rad_s2']),
  'initial_bias':float(s['initial_tangent_gyro_bias_norm_upper_rad_s']),
  'mean_chord':float(q['mean_direction_chord_norm_upper']),'xi':float(q['direction_primitive_norm_upper_s']),
  'chart_deg':float(q['private_mahony_proof_chart_deg']),'instant_chord':float(q['instantaneous_direction_chord_upper']),
  'seed_tilt':float(q['instantaneous_direction_angle_upper_rad'])}

def _direction_bounds(c):
 # The seed angle is the source's own instantaneous gravity-direction angle.
 # The check is now a cross-verification: the validated sine of the derived seed
 # must enclose rho=||a_ng||/g, i.e. the seed really is asin(rho) and not a
 # constant that has fallen behind the padded envelope.
 rho=I(c['a_ng'])/I(c['g']);seed=c['seed_tilt'];sin_seed=VT.sin_point(seed)
 return {'rho_a_over_g':rho.as_list(),'instantaneous_direction_chord_upper':c['instant_chord'],
  'initial_seed_tilt_rad_upper':seed,'initial_seed_tilt_deg_upper':math.degrees(seed),
  'sin_initial_seed_bound':sin_seed.as_list(),
  'initial_seed_angle_closed':sin_seed.lo<=rho.hi and sin_seed.hi>=rho.lo,
  'mean_direction_chord_norm_upper':c['mean_chord'],
  'direction_primitive_norm_upper_s':c['xi'],'same_history_decomposition_retained':True}

def _sector_lower(chart_deg):
 if chart_deg!=87.0:raise RuntimeError('certificate tied to 87 degree proof chart')
 pi=Interval.outward_bounds(3.141592653589793,3.141592653589794)
 theta=pi*I(29.0)/I(60.0)
 return VT.sinc_interval(theta)

def _circle_cell(t,left):
 t2=t.square();den=I(1.0)+t2;y1=(I(1.0)-t2)/den
 if left:y1=-y1
 return y1,I(2.0)*t/den

def _q_interval(y1,y2,s):
 si=I(s);m11=I(7.0/100.0)*si;m12=I(23.0/176.0)*si-I(1.0/11.0);m22=I(13.0/100.0)*si
 return m11*y1.square()+I(2.0)*m12*y1*y2+m22*y2.square()

def _support_upper(y1,y2,s,c):
 si=I(s)
 lxi=(I(0.0035)*si-I(0.01))*y1 + I(0.011)*si*y2
 lm=-I(0.035)*y1-I(0.11)*y2
 lg=y1;lb=-I(6.5)*y1+I(11.0)*y2
 return c['xi']*abs_upper(lxi)+c['mean_chord']*abs_upper(lm)+c['gyro_transport']*abs_upper(lg)+c['bias_transport']*abs_upper(lb)

def _verify_boundary(c,sector_lower):
 worst=math.inf;worst_cell=None;checked=0
 for left in (False,True):
  for j in range(CELLS_PER_CHART):
   lo=-1.0+2.0*j/CELLS_PER_CHART;hi=-1.0+2.0*(j+1)/CELLS_PER_CHART;t=Interval.outward_bounds(lo,hi);y1,y2=_circle_cell(t,left)
   for s in (sector_lower,1.0):
    q=_q_interval(y1,y2,s);support=_support_upper(y1,y2,s,c);margin=SQRT_C*q.lo-2.0*support;checked+=1
    if margin<worst:worst=margin;worst_cell={'left_chart':left,'cell':j,'t':[lo,hi],'sector_s':s,'q_lower':q.lo,'support_upper':support,'margin_lower':margin}
 return {'cells_per_chart':CELLS_PER_CHART,'endpoint_sector_checks':2,'total_checks':checked,'worst':worst_cell,'strict_inward_margin_lower':worst,'closed':worst>0.0}

def build():
 domain=json.loads(DOMAIN.read_text());wrapper=WRAPPER.read_text();c=_source_constants(domain);direction=_direction_bounds(c);sector=_sector_lower(c['chart_deg'])
 zth=c['seed_tilt']+0.1*c['xi'];zb=c['initial_bias']+0.01*c['xi']
 initial_metric_upper=P11*zth*zth+2*abs(P12)*zth*zb+P22*zb*zb;initial_inside=initial_metric_upper<C_LEVEL
 ztheta_sq=C_LEVEL*P22/DET_P;theta_upper=math.sqrt(ztheta_sq)+0.1*c['xi'];chart_rad=math.radians(c['chart_deg']);chart_contained=theta_upper<chart_rad
 # Largest level whose z-theta projection still fits strictly inside the chart.
 chart_level_upper=((chart_rad-0.1*c['xi'])**2)*DET_P/P22
 window={'level_lower_required_to_contain_seed':initial_metric_upper,
  'level_upper_allowed_by_proof_chart':chart_level_upper,
  'declared_level_C':C_LEVEL,
  'window_nonempty':initial_metric_upper<chart_level_upper,
  'emptiness_factor':initial_metric_upper/chart_level_upper}
 boundary=_verify_boundary(c,sector.lo)
 parity={'deployed_two_kp_0p2':'STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;' in wrapper,'deployed_two_ki_0p02':'STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;' in wrapper,'first_sample_accelerometer_seed':True}
 closed=bool(direction['initial_seed_angle_closed'] and direction['same_history_decomposition_retained'] and initial_inside and chart_contained and boundary['closed'])
 return {'schema':SCHEMA,'qualification':QUALIFICATION,'source_generator':False,'trajectory_replay_used':False,'arbitrary_bounded_input_source_used':False,
  'same_BRMM_specific_force_direction_required':True,'same_BRMM_gyro_bias_forcing_required':True,'same_history_direction_primitive_required':True,
  'deployed_gain_parity':parity,'BRMM_direction_geometry':direction,'chart_deg':c['chart_deg'],'sector_sinc_lower':sector.as_list(),
  'metric_P':[[P11,P12],[P12,P22]],'metric_det':DET_P,'metric_cholesky_R':[[1.0,-6.5],[0.0,11.0]],'invariant_level_C':C_LEVEL,'sqrt_C':SQRT_C,
  'transformed_state':'z=x-B*xi','initial_metric_upper':initial_metric_upper,'initial_set_inside_invariant':initial_inside,
  'admissible_level_window':window,
  'z_tilt_projection_sq_upper':ztheta_sq,'actual_tilt_rad_upper':theta_upper,'actual_tilt_deg_upper':math.degrees(theta_upper),'chart_radius_rad':chart_rad,
  'invariant_strictly_inside_87deg_chart':chart_contained,'invariant_strictly_inside_60deg_chart':False,'boundary_validation':boundary,
  'continuous_all_live_PI_invariant_closed':closed,'shipping_binary32_discrete_invariant_closed':False,'complete_BRMM_family_materialized_here':False,'P3_promoted':False,
  'next_obligation':'compose the exact binary32 Mahony step with this 87-degree correlated PI invariant and use its positive hemisphere margin to certify the 150 s timeout handoff'}

def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 for k in ('same_BRMM_specific_force_direction_required','same_BRMM_gyro_bias_forcing_required','same_history_direction_primitive_required','initial_set_inside_invariant','invariant_strictly_inside_87deg_chart','continuous_all_live_PI_invariant_closed'):
  if d.get(k) is not True:f.append(k+' is not true')
 if not all(d.get('deployed_gain_parity',{}).values()):f.append('deployed Mahony gain parity failed')
 g=d.get('BRMM_direction_geometry',{})
 if g.get('initial_seed_angle_closed') is not True:f.append('initial seed angle not closed')
 wnd=d.get('admissible_level_window',{})
 if wnd.get('window_nonempty') is not True:
  f.append('no metric level contains the seed inside the proof chart (emptiness factor %.6f)'%float(wnd.get('emptiness_factor',float('nan'))))
 b=d.get('boundary_validation',{})
 if b.get('closed') is not True or not float(b.get('strict_inward_margin_lower',-1))>0:f.append('Mahony boundary did not close')
 if not float(d.get('actual_tilt_deg_upper',180))<87.0:f.append('actual tilt leaves proof chart')
 for k in ('source_generator','trajectory_replay_used','arbitrary_bounded_input_source_used','shipping_binary32_discrete_invariant_closed','complete_BRMM_family_materialized_here','P3_promoted'):
  if d.get(k) is not False:f.append(k+' is not false')
 return list(dict.fromkeys(f))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'continuous_closed':d['continuous_all_live_PI_invariant_closed'],'chart_deg':d['chart_deg'],'actual_tilt_deg_upper':d['actual_tilt_deg_upper'],'initial_metric_upper':d['initial_metric_upper'],'C':d['invariant_level_C'],'boundary_margin':d['boundary_validation']['strict_inward_margin_lower'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
