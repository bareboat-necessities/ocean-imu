#!/usr/bin/env python3
"""Two-phase all-Live PI certificate for the private Mahony observer.

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
coordinate rather than replaced by an arbitrary r box. The actual tilt is
recovered as theta=z_theta-0.1 xi.

## Why a single level failed and what replaced it

The metric is `P=R^T R` with `R=[[1,-p],[0,c]]`, so `r_22=p^2+c^2` and
`det P=c^2`. Writing `u=Rz`, `y=u/||u||` and `M=R A_s R^-1`,

    Vdot = C y^T(M+M^T)y + 2 sqrt(C) y^T R w = -sqrt(C) [ sqrt(C) q - 2 sup ],

with `q(y,s)=-y^T(M+M^T)y` and `sup` the support of the admitted forcing. Both
`q` and `sup` are level-independent, so the boundary condition is a *lower*
bound on the level,

    sqrt(C) >= W_in := 2 sup_max / q_min,

while the 87 degree proof chart imposes an *upper* bound `C r_22/c^2 <
(chart - 0.1 xi)^2`, and containment of the first-sample accelerometer seed
imposes a second lower bound `C >= C_seed`.

At the padded 8.8 m/s^2 envelope the seed angle is `asin(8.8/g)=1.1137` rad and
the previous metric `p=6.5, c=11` left those three bounds mutually
unsatisfiable: the level needed for the seed was `1.8540056` against a chart
ceiling of `1.4912551`. That is a failure of the *single-level* formulation, not
of the enclosure -- the boundary flow still closed with margin `0.0357940`.

This certificate replaces it with two nested levels of one better-conditioned
metric. Because `Vdot<0` holds on *every* level with `sqrt(C)>=W_in`, not only
on the certified boundary:

* the outer level `C_out` contains the seed and lies inside the chart, so
  `{V<=C_out}` is forward invariant and the observer never leaves the chart;
* `W=sqrt(V)` obeys `Wdot <= -(q_min/2)(W-W_in)`, so the state enters the inner
  level `C_in` exponentially with time constant `2/q_min`, giving a genuine
  ultimate tilt bound instead of an outer level radius.

Both levels are verified with the same rational stereographic boundary cover,
outward interval arithmetic and the endpoint sector `s in [sinc(87 deg),1]`;
`q` is affine and `sup` convex in `s`, so the margin is concave in `s` and
endpoint checking is rigorous.

The metric is chosen for seed containment and chart retention. It is *not*
tuned toward the magnetic capture requirement: a non-promoting scan of the whole
`(p,c)` class puts the best certifiable ultimate tilt near 41 deg, far above
both the `6.1145 deg` north-capture and `2.8378 deg` horizontal-fraction
requirements, so no member of this class supplies the magnetic gates. See
`docs/ou3-proof-research-state.md`.
"""
from __future__ import annotations
import argparse,json,math
from fractions import Fraction as F
from pathlib import Path
from ou3_interval import Interval,down,up
import ou3_validated_transcendentals as VT
import ou3_brmm_gravity_direction_forcing_qualification as DIR

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
SCHEMA=3
QUALIFICATION='OU3_BRMM_PRIVATE_MAHONY_TWO_PHASE_PI_INVARIANT_V3'

# Metric R=[[1,-P_SKEW],[0,C_SCALE]] chosen for seed containment inside the
# 87 degree chart. det P = C_SCALE^2 is unchanged from the retired metric.
P_SKEW=F(3)
C_SCALE=F(11)
SQRT_C_OUT=F(133,100)
SQRT_C_IN=F(1)
# Rational upper bound for the seed angle asin(8.8/g), cross-verified below
# against the validated sine rather than assumed.
SEED_TILT_RAD_UPPER=F(111373,100000)
CAPTURE_HORIZON_S=F(150)
CELLS_PER_CHART=8192

def I(x):return Interval.outward_bounds(float(x),float(x))
def abs_upper(x):return max(abs(x.lo),abs(x.hi))
def _r22():return P_SKEW*P_SKEW+C_SCALE*C_SCALE
def _det():return C_SCALE*C_SCALE

def _source_constants(domain):
 s=domain['startup'];live=domain['normal_live'];q=DIR.build();f=DIR.validate(q)
 if f:raise RuntimeError('gravity-direction source qualification invalid: '+repr(f))
 return {'g':float(s['gravity_mps2']),'a_ng':float(live['non_gravitational_cog_acceleration_norm_upper_mps2']),
  'gyro_transport':float(s['effective_deterministic_gyro_transport_disturbance_upper_rad_s']),
  'bias_transport':float(s['effective_deterministic_bias_transport_disturbance_upper_rad_s2']),
  'initial_bias':float(s['initial_tangent_gyro_bias_norm_upper_rad_s']),
  'mean_chord':float(q['mean_direction_chord_norm_upper']),'xi':float(q['direction_primitive_norm_upper_s']),
  'chart_deg':float(q['private_mahony_proof_chart_deg']),'instant_chord':float(q['instantaneous_direction_chord_upper'])}

def _direction_bounds(c):
 # The declared rational seed bound must actually dominate the source's own
 # instantaneous gravity-direction angle: sin is increasing on [0,pi/2], so
 # sin(SEED_TILT_RAD_UPPER) >= rho = ||a_ng||/g certifies asin(rho) <= bound.
 rho=I(c['a_ng'])/I(c['g']);seed=float(SEED_TILT_RAD_UPPER);sin_seed=VT.sin_point(seed)
 return {'rho_a_over_g':rho.as_list(),'instantaneous_direction_chord_upper':c['instant_chord'],
  'seed_tilt_rad_upper':seed,'seed_tilt_deg_upper':math.degrees(seed),
  'sin_seed_tilt_bound':sin_seed.as_list(),
  'seed_tilt_dominates_source_angle':sin_seed.lo>=rho.hi,
  'mean_direction_chord_norm_upper':c['mean_chord'],
  'direction_primitive_norm_upper_s':c['xi'],'same_history_decomposition_retained':True}

def _sector_lower(chart_deg):
 if chart_deg!=87.0:raise RuntimeError('certificate tied to 87 degree proof chart')
 theta=_chart_interval(chart_deg)
 return VT.sinc_interval(theta)

def _chart_interval(chart_deg):
 if chart_deg!=87.0:raise RuntimeError('certificate tied to 87 degree proof chart')
 pi=Interval.outward_bounds(3.141592653589793,3.141592653589794)
 return pi*I(29.0)/I(60.0)

def _circle_cell(t,left):
 t2=t.square();den=I(1.0)+t2;y1=(I(1.0)-t2)/den
 if left:y1=-y1
 return y1,I(2.0)*t/den

def _q_interval(y1,y2,s):
 """q = -y'(M+M')y with M = R A_s R^-1; affine in the sector value s."""
 si=I(s);p=I(float(P_SKEW));cc=I(float(C_SCALE))
 m11=I(0.02)*si*(I(10.0)-p)
 m22=I(0.02)*p*si
 m12=(p*si*(I(0.1)-I(0.01)*p)-I(1.0))/cc+I(0.01)*cc*si
 return m11*y1.square()+I(2.0)*m12*y1*y2+m22*y2.square()

def _support_upper(y1,y2,s,c):
 """Support of R*(A_s B xi + B m + d) over the admitted same-history forcing."""
 si=I(s);p=I(float(P_SKEW));cc=I(float(C_SCALE))
 lxi=(I(0.01)*si-I(0.001)*p*si-I(0.01))*y1 + I(0.001)*cc*si*y2
 lm=(-I(0.1)+I(0.01)*p)*y1 + (-I(0.01)*cc)*y2
 lg=y1;lb=-p*y1+cc*y2
 return c['xi']*abs_upper(lxi)+c['mean_chord']*abs_upper(lm)+c['gyro_transport']*abs_upper(lg)+c['bias_transport']*abs_upper(lb)

def _verify_boundary(c,sector_lower,sqrt_c):
 """Strict inward flow on {V = sqrt_c^2}: sqrt_c*q - 2*support > 0."""
 worst=math.inf;worst_cell=None;checked=0;q_floor=math.inf;support_ceiling=0.0
 for left in (False,True):
  for j in range(CELLS_PER_CHART):
   lo=-1.0+2.0*j/CELLS_PER_CHART;hi=-1.0+2.0*(j+1)/CELLS_PER_CHART;t=Interval.outward_bounds(lo,hi);y1,y2=_circle_cell(t,left)
   for s in (sector_lower,1.0):
    q=_q_interval(y1,y2,s);support=_support_upper(y1,y2,s,c);margin=float(sqrt_c)*q.lo-2.0*support;checked+=1
    if q.lo<q_floor:q_floor=q.lo
    if support>support_ceiling:support_ceiling=support
    if margin<worst:worst=margin;worst_cell={'left_chart':left,'cell':j,'t':[lo,hi],'sector_s':s,'q_lower':q.lo,'support_upper':support,'margin_lower':margin}
 return {'cells_per_chart':CELLS_PER_CHART,'endpoint_sector_checks':2,'total_checks':checked,
  'level_sqrt_C':float(sqrt_c),'q_lower_floor':q_floor,'support_upper_ceiling':support_ceiling,
  'worst':worst_cell,'strict_inward_margin_lower':worst,'closed':worst>0.0}

def _exp_neg_upper(x:float)->float:
 """Outward upper bound for exp(-x), x>=0, by composing validated half-steps."""
 if x<0:raise ValueError('nonnegative argument required')
 n=max(1,int(math.ceil(x/0.5)))
 step=VT.exp_point(-x/n)
 acc=Interval.outward_bounds(1.0,1.0)
 for _ in range(n):acc=acc*step
 return acc.hi

def _tilt_upper(level_sqrt:float,xi:float)->float:
 """theta bound from a level: sqrt(C*r22/det) + 0.1*xi, outward rounded."""
 scale=up(math.sqrt(up(float(_r22())/float(_det()))))
 return up(up(level_sqrt*scale)+up(0.1*xi))

def build():
 domain=json.loads(DOMAIN.read_text());wrapper=WRAPPER.read_text()
 c=_source_constants(domain);direction=_direction_bounds(c);sector=_sector_lower(c['chart_deg'])
 r22=_r22();det=_det();xi=c['xi']
 # Seed level, worst relative sign between the seed tilt and the initial bias.
 zth=SEED_TILT_RAD_UPPER+F(1,10)*F(xi).limit_denominator(10**6)
 zb=F(c['initial_bias']).limit_denominator(10**9)+F(1,100)*F(xi).limit_denominator(10**6)
 c_seed=zth*zth+2*P_SKEW*zth*zb+r22*zb*zb
 c_out=SQRT_C_OUT*SQRT_C_OUT;c_in=SQRT_C_IN*SQRT_C_IN
 # Chart ceiling: theta-extent + 0.1*xi must stay strictly inside the chart.
 chart=_chart_interval(c['chart_deg'])
 chart_room=chart-I(0.1*xi)
 c_chart_upper=down(chart_room.lo*chart_room.lo*float(det)/float(r22))
 seed_contained=float(c_seed)<=float(c_out)
 chart_contained=float(c_out)<c_chart_upper
 outer=_verify_boundary(c,sector.lo,SQRT_C_OUT)
 inner=_verify_boundary(c,sector.lo,SQRT_C_IN)
 # Capture: Wdot <= -(q/2)(W - W_in) with W_in = 2*sup/q.
 q_min=min(outer['q_lower_floor'],inner['q_lower_floor'])
 support_max=max(outer['support_upper_ceiling'],inner['support_upper_ceiling'])
 w_in=up(2.0*support_max/down(q_min))
 rate=down(q_min/2.0)
 w_out=float(SQRT_C_OUT)
 decay=_exp_neg_upper(down(rate*float(CAPTURE_HORIZON_S)))
 w_horizon=up(w_in+up(up(w_out-w_in)*decay))
 inner_reached=w_in<=float(SQRT_C_IN)
 outer_tilt=_tilt_upper(w_out,xi)
 horizon_tilt=_tilt_upper(w_horizon,xi)
 # W(t) -> W_in, so the asymptotic bound comes from the minimal certifiable
 # level, not from the nested level whose boundary is separately verified.
 ultimate_tilt=_tilt_upper(w_in,xi)
 inner_level_tilt=_tilt_upper(float(SQRT_C_IN),xi)
 chart_rad=chart.hi
 parity={'deployed_two_kp_0p2':'STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;' in wrapper,'deployed_two_ki_0p02':'STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;' in wrapper,'first_sample_accelerometer_seed':True}
 closed=bool(direction['seed_tilt_dominates_source_angle'] and direction['same_history_decomposition_retained']
  and seed_contained and chart_contained and outer['closed'] and inner['closed'] and inner_reached
  and outer_tilt<chart_rad)
 return {'schema':SCHEMA,'qualification':QUALIFICATION,'source_generator':False,'trajectory_replay_used':False,'arbitrary_bounded_input_source_used':False,
  'same_BRMM_specific_force_direction_required':True,'same_BRMM_gyro_bias_forcing_required':True,'same_history_direction_primitive_required':True,
  'deployed_gain_parity':parity,'BRMM_direction_geometry':direction,'chart_deg':c['chart_deg'],'sector_sinc_lower':sector.as_list(),
  'metric_P':[[1.0,-float(P_SKEW)],[-float(P_SKEW),float(r22)]],'metric_det':float(det),
  'metric_cholesky_R':[[1.0,-float(P_SKEW)],[0.0,float(C_SCALE)]],
  'metric_skew_p':float(P_SKEW),'metric_scale_c':float(C_SCALE),
  'invariant_level_C':float(c_out),'sqrt_C':float(SQRT_C_OUT),
  'transformed_state':'z=x-B*xi',
  'two_phase_levels':{'seed_level_required':float(c_seed),'outer_level_C':float(c_out),
   'chart_level_ceiling':c_chart_upper,'inner_level_C':float(c_in),
   'seed_contained_in_outer_level':seed_contained,'outer_level_inside_chart':chart_contained,
   'chart_headroom_relative':down((c_chart_upper-float(c_out))/c_chart_upper),
   'outer_boundary_margin_lower':outer['strict_inward_margin_lower'],
   'inner_boundary_margin_lower':inner['strict_inward_margin_lower']},
  'initial_metric_upper':float(c_seed),'initial_set_inside_invariant':seed_contained,
  'capture':{'q_lower_floor':q_min,'support_upper_ceiling':support_max,
   'minimal_certifiable_sqrt_level':w_in,'inner_level_above_minimal':inner_reached,
   'contraction_rate_lower_per_s':rate,'time_constant_upper_s':up(2.0/down(q_min)),
   'deployed_horizon_s':float(CAPTURE_HORIZON_S),'horizon_decay_upper':decay,
   'sqrt_level_at_horizon_upper':w_horizon,
   'monotone_decrease_above_inner_level':True,'capture_uses_deployed_timeout_as_hypothesis':False},
  'actual_tilt_rad_upper':outer_tilt,'actual_tilt_deg_upper':math.degrees(outer_tilt),
  'tilt_at_deployed_horizon_rad_upper':horizon_tilt,'tilt_at_deployed_horizon_deg_upper':math.degrees(horizon_tilt),
  'ultimate_tilt_rad_upper':ultimate_tilt,'ultimate_tilt_deg_upper':math.degrees(ultimate_tilt),
  'certified_inner_level_tilt_rad_upper':inner_level_tilt,'certified_inner_level_tilt_deg_upper':math.degrees(inner_level_tilt),
  'chart_radius_rad':chart_rad,
  'invariant_strictly_inside_87deg_chart':chart_contained and outer_tilt<chart_rad,
  'invariant_strictly_inside_60deg_chart':False,
  'boundary_validation':outer,'inner_boundary_validation':inner,
  'continuous_all_live_PI_invariant_closed':closed,'shipping_binary32_discrete_invariant_closed':False,
  'complete_BRMM_family_materialized_here':False,'P3_promoted':False,
  'metric_tuned_toward_magnetic_requirement':False,
  'next_obligation':'compose the exact binary32 Mahony step with these two levels and supply an accumulation tilt-frame bound the magnetic capture can consume; the quadratic-metric class cannot reach it'}

def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 for k in ('same_BRMM_specific_force_direction_required','same_BRMM_gyro_bias_forcing_required','same_history_direction_primitive_required','initial_set_inside_invariant','invariant_strictly_inside_87deg_chart','continuous_all_live_PI_invariant_closed'):
  if d.get(k) is not True:f.append(k+' is not true')
 if not all(d.get('deployed_gain_parity',{}).values()):f.append('deployed Mahony gain parity failed')
 g=d.get('BRMM_direction_geometry',{})
 if g.get('seed_tilt_dominates_source_angle') is not True:f.append('declared seed angle does not dominate the source angle')
 lv=d.get('two_phase_levels',{})
 for k in ('seed_contained_in_outer_level','outer_level_inside_chart'):
  if lv.get(k) is not True:f.append(k+' is not true')
 if not float(lv.get('chart_headroom_relative',-1.0))>0.0:f.append('outer level has no chart headroom')
 for b in ('boundary_validation','inner_boundary_validation'):
  x=d.get(b,{})
  if x.get('closed') is not True or not float(x.get('strict_inward_margin_lower',-1))>0:f.append(b+' did not close')
 cap=d.get('capture',{})
 if cap.get('inner_level_above_minimal') is not True:f.append('inner level below the minimal certifiable level')
 if not float(cap.get('contraction_rate_lower_per_s',0.0))>0.0:f.append('capture rate not positive')
 if not float(cap.get('sqrt_level_at_horizon_upper',math.inf))<float(d.get('sqrt_C',0.0)):f.append('no strict decrease over the deployed horizon')
 if not float(d.get('actual_tilt_deg_upper',180))<87.0:f.append('actual tilt leaves proof chart')
 if not float(d.get('ultimate_tilt_deg_upper',180))<float(d.get('tilt_at_deployed_horizon_deg_upper',180))<float(d.get('actual_tilt_deg_upper',180)):f.append('tilt bounds are not ordered ultimate < horizon < all-time')
 for k in ('source_generator','trajectory_replay_used','arbitrary_bounded_input_source_used','shipping_binary32_discrete_invariant_closed','complete_BRMM_family_materialized_here','P3_promoted','metric_tuned_toward_magnetic_requirement'):
  if d.get(k) is not False:f.append(k+' is not false')
 return list(dict.fromkeys(f))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'closed':d['continuous_all_live_PI_invariant_closed'],
  'seed_deg':d['BRMM_direction_geometry']['seed_tilt_deg_upper'],
  'levels':[d['two_phase_levels']['seed_level_required'],d['two_phase_levels']['outer_level_C'],d['two_phase_levels']['chart_level_ceiling']],
  'chart_headroom':d['two_phase_levels']['chart_headroom_relative'],
  'outer_margin':d['boundary_validation']['strict_inward_margin_lower'],
  'inner_margin':d['inner_boundary_validation']['strict_inward_margin_lower'],
  'tilt_all_time_deg':d['actual_tilt_deg_upper'],
  'tilt_150s_deg':d['tilt_at_deployed_horizon_deg_upper'],
  'tilt_ultimate_deg':d['ultimate_tilt_deg_upper'],
  'failures':f},sort_keys=True))
 return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
