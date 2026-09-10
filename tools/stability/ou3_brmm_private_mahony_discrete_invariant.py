#!/usr/bin/env python3
"""Source-order binary32/discrete charge for the widened Mahony invariant.

This certificate does not replay a trajectory. It takes the already validated
87-degree continuous transformed PI invariant and charges the literal shipping
Euler/quaternion evaluation at dt=5 ms. The fast-inverse-square-root scale is a
common positive quaternion scale and therefore does not rotate the attitude;
only componentwise binary32 rounding after the common scaling is charged as an
orientation perturbation.

The proof uses conservative gamma_n absolute-error budgets for the raw
quaternion component expression and Ki accumulation. It also adds the exact
forward-Euler h^2 f^T P f term to the quadratic boundary inequality. Compiler
reassociation/FMA is intentionally separate: this closes the declared source
order binary32 model, not yet every toolchain realization.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_private_mahony_live_invariant as CONT
import ou3_fast_inv_sqrt_interval as FINV

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
SCHEMA=1
QUALIFICATION='OU3_BRMM_PRIVATE_MAHONY_SOURCE_ORDER_BINARY32_INVARIANT_V1'
U=2.0**-24

def gamma(n):
 x=n*U
 if not x<1:raise RuntimeError('gamma_n overflow')
 return x/(1-x)

def build():
 c=CONT.build();cf=CONT.validate(c)
 if cf:raise RuntimeError('continuous Mahony invariant invalid: '+repr(cf))
 d=json.loads(DOMAIN.read_text());dt=float(d['configured_runtime']['imu_dt_s']);body=math.radians(float(d['normal_live']['body_rate_norm_upper_deg_s']))
 shell=FINV.all_positive_normal_normalized_norm2_enclosure();qn_lo=math.sqrt(shell.lo);qn_hi=math.sqrt(shell.hi)
 C=float(c['invariant_level_C']);P22=float(c['metric_P'][1][1]);det=float(c['metric_det']);xi=float(c['BRMM_direction_geometry']['direction_primitive_norm_upper_s']);mean=float(c['BRMM_direction_geometry']['mean_direction_chord_norm_upper'])
 ztheta=math.sqrt(C*P22/det);zbeta=math.sqrt(C/det);theta=ztheta+0.1*xi;beta=zbeta+0.01*xi
 half_g=0.5*qn_hi*qn_hi;half_err=qn_hi*half_g
 omega=body+beta+0.2*half_err
 half_rate=0.5*dt*omega
 # q_i <- q_i + three signed q_j*(0.5*dt*omega_j) terms. 24 ops is
 # deliberately above the literal unfused count including the preceding rate
 # scaling and three-term sum; cancellation is not used.
 raw_sum_abs=qn_hi*(1.0+3.0*half_rate)
 raw_component_round=gamma(24)*raw_sum_abs
 raw_vector_round=2.0*raw_component_round
 raw_angle_round=2.0*raw_vector_round/qn_lo
 # Common invSqrt scale does not rotate; only four final multiplies round.
 norm_mult_vector_round=2.0*U*qn_hi
 norm_angle_round=2.0*norm_mult_vector_round/qn_lo
 # Error-vector/gyro arithmetic is charged separately with a conservative 24-op
 # gamma against the complete corrected-rate magnitude.
 feedback_rate_round=gamma(24)*omega
 attitude_rate_fp=(raw_angle_round+norm_angle_round)/dt+feedback_rate_round
 # Ki recurrence: old integral magnitude is dominated by beta plus initial
 # physical gyro-bias. Three literal operations plus one accumulation; use 8.
 initial_bias=float(d['startup']['initial_tangent_gyro_bias_norm_upper_rad_s'])
 integral_mag=beta+initial_bias;ki_increment=0.02*half_err*dt
 integral_step_round=gamma(8)*(integral_mag+ki_increment)
 bias_rate_fp=integral_step_round/dt
 # Robustify continuous boundary margin against the extra d_g,d_b support.
 bias_row_norm=math.hypot(6.5,11.0)
 fp_support=attitude_rate_fp+bias_row_norm*bias_rate_fp
 continuous_margin=float(c['boundary_validation']['strict_inward_margin_lower'])
 robust_margin=continuous_margin-2.0*fp_support
 # Bound transformed tangent derivative and the forward-Euler quadratic term.
 smin=float(c['sector_sinc_lower'][0]);dg=float(d['startup']['effective_deterministic_gyro_transport_disturbance_upper_rad_s'])+attitude_rate_fp;db=float(d['startup']['effective_deterministic_bias_transport_disturbance_upper_rad_s2'])+bias_rate_fp
 ftheta=0.1*ztheta+zbeta+max(abs(0.01*smin-0.01),0.0)*xi+0.1*mean+dg
 fbeta=0.01*ztheta+0.001*xi+0.01*mean+db
 Rf1=ftheta+6.5*fbeta;Rf2=11.0*fbeta
 euler_quadratic=dt*dt*(Rf1*Rf1+Rf2*Rf2)
 # From the continuous certificate convention dot(V)<=-2*sqrt(C)*margin.
 guaranteed_first_order_decrease=2.0*float(c['sqrt_C'])*robust_margin*dt
 discrete_margin=guaranteed_first_order_decrease-euler_quadratic
 chart_margin=math.radians(float(c['chart_deg']))-float(c['actual_tilt_rad_upper'])
 # One-step angle rounding must fit comfortably inside geometric margin; the
 # invariant inequality, rather than accumulation, handles indefinite steps.
 one_step_chart_round_margin=chart_margin-(raw_angle_round+norm_angle_round)
 closed=bool(robust_margin>0 and discrete_margin>0 and one_step_chart_round_margin>0)
 return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','source_order_binary32_model':True,'trajectory_replay_used':False,'compiler_reassociation_or_FMA_closed':False,
  'continuous_invariant_consumed':True,'dt_s':dt,'quaternion_norm_shell':[qn_lo,qn_hi],'half_error_norm_upper':half_err,'corrected_gyro_norm_upper_rad_s':omega,
  'gamma24':gamma(24),'gamma8':gamma(8),'raw_quaternion_component_roundoff_upper':raw_component_round,'raw_orientation_roundoff_upper_rad':raw_angle_round,'normalization_multiply_orientation_roundoff_upper_rad':norm_angle_round,
  'equivalent_attitude_rate_roundoff_upper_rad_s':attitude_rate_fp,'equivalent_integral_bias_rate_roundoff_upper_rad_s2':bias_rate_fp,
  'continuous_boundary_margin_before_binary32':continuous_margin,'binary32_support_charge':fp_support,'continuous_boundary_margin_after_binary32':robust_margin,
  'forward_euler_quadratic_V_upper':euler_quadratic,'guaranteed_first_order_V_decrease_lower':guaranteed_first_order_decrease,'discrete_V_margin_lower':discrete_margin,
  'continuous_chart_margin_rad':chart_margin,'one_step_chart_round_margin_rad':one_step_chart_round_margin,
  'shipping_binary32_discrete_invariant_closed_conditionally_on_source_order':closed,'toolchain_independent_binary32_invariant_closed':False,'P4_promoted_here':False,'P5_promoted_here':False,
  'next_obligation':'bind the shipping Mahony state-step/source lineage to this source-order discrete invariant and separately qualify compiler contraction/FMA behavior'}

def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 for k in ('source_order_binary32_model','continuous_invariant_consumed','shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('trajectory_replay_used','compiler_reassociation_or_FMA_closed','toolchain_independent_binary32_invariant_closed','P4_promoted_here','P5_promoted_here'):
  if d.get(k) is not False:f.append(k+' not false')
 for k in ('continuous_boundary_margin_after_binary32','discrete_V_margin_lower','one_step_chart_round_margin_rad'):
  if not float(d.get(k,0))>0:f.append(k+' not positive')
 return f

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'binary32_closed':d['shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'],'attitude_fp_rate':d['equivalent_attitude_rate_roundoff_upper_rad_s'],'bias_fp_rate':d['equivalent_integral_bias_rate_roundoff_upper_rad_s2'],'robust_margin':d['continuous_boundary_margin_after_binary32'],'discrete_margin':d['discrete_V_margin_lower'],'chart_margin':d['one_step_chart_round_margin_rad'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
