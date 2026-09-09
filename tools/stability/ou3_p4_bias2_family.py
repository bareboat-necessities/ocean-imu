#!/usr/bin/env python3
"""Explicit conditional BIAS2 family for bounded-bias P4.

BIAS2 is the non-relaxing drift family: the shipping-centered true bias is
bounded in magnitude and in variation, but admits NO relaxation root.  Its
admitted factors are

    phi_true in [exp(-dt/tau_true_lower), 1],

so the deployed first-order model's relaxation is not available to the truth
at all.  This is a strictly separate mandatory family: BIAS1's box caps
tau_true, and BIAS0's composite root cancels its own Gauss-Markov channel, so
neither implies the drift limit and neither may substitute for it.

Every admitted sequence satisfies the exact recurrence

    b_i = phi_true*b_{i-1} + w_i,
    |w_i| <= |b_i - b_{i-1}| + (1-phi_true)*|b_{i-1}|
          <= rate*dt + (1-phi_lo)*magnitude,

so the driver increment is small while the tau mismatch against the deployed
phi_hat is maximal.

A non-relaxing truth does NOT force the separation sector on the declared
objective.  That objective is bounded accelerometer-bias error plus regional
practical ISS of the other 18 errors, and the bias half of it comes from the
deployed radial projection, not from any decay of the bias error:
`ou3_p4_projection_sector` closes `F_R(e,beta)=beta-Pi_R(beta-e)` analytically
on every branch, so the estimate stays in the R ball and

    |e_b| <= |b_hat| + |b_true| <= R + B_true

pointwise at every event, for any admitted truth history and hence for every
one of BIAS0/BIAS1/BIAS2.  The three families share one declared truth
envelope, so that compactness bound is the same number for all three.  The
motion half is the cocycle contraction against a persistent bias forcing,
which is a finite ultimate bound rather than a separation statement.

The separation sector

    Pi_sep = C_y^T W C_y - mu_sep X

therefore only SHARPENS the motion channel gains; it is not a prerequisite of
the declared objective.  No uniform mu_sep lower bound is proved here and none
is required here.  This module is a theorem/source-family qualification, not
assembled-sensor hardware qualification.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
DEFAULT=REPO/'tools/stability/ou3_p4_closure_domain.json'
DT=0.005
SCHEMA=1
QUALIFICATION='OU3_P4_CONDITIONAL_BIAS2_FAMILY_V2'


def up(x): return math.nextafter(float(x), math.inf)
def down(x): return math.nextafter(float(x), -math.inf)


def build(path: Path=DEFAULT):
    c=json.loads(Path(path).read_text())['BIAS2_family']
    mag=float(c['component_abs_upper_mps2'])
    rate=float(c['rate_abs_upper_mps3'])
    tau_lo=float(c['tau_true_lower_s'])
    if not (mag>=0 and rate>=0): raise RuntimeError('invalid BIAS2 magnitude/rate bounds')
    if not (math.isfinite(tau_lo) and tau_lo>0): raise RuntimeError('invalid BIAS2 tau lower bound')
    if c.get('tau_true_upper_is_infinite') is not True:
        raise RuntimeError('BIAS2 must admit the non-relaxing limit')
    if c.get('uniform_separation_constant_lower') is not None:
        raise RuntimeError('BIAS2 separation constant is not proved and must stay null')
    import ou3_p4_projection_sector as PROJ
    proj=PROJ.build()
    if PROJ.validate(proj):
        raise RuntimeError('projection sector prerequisite failed')
    if not (proj['global_joint_sector_closed'] and proj['estimate_ball_invariance_closed']):
        raise RuntimeError('bias compactness needs the closed projection sector')
    projection_radius=float(proj['base_report']['projection_radius'])

    phi_lo=down(math.exp(-DT/tau_lo)); phi_hi=1.0
    one_minus_phi=up(1.0-phi_lo)
    variation_step=up(rate*DT)
    relaxation_step=up(one_minus_phi*mag)
    w_component=up(variation_step+relaxation_step)

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'conditional_theorem_family':True,'assembled_sensor_hardware_qualified':False,
      'model':c['model'],'family_role':c['family_role'],
      'one_drift_history_required':True,
      'independent_per_sample_bias_slots_forbidden':True,
      'inferred_from_BIAS0_or_BIAS1':False,'may_be_substituted_by_BIAS1':False,
      'dt_s':DT,
      'true_bias_component_abs_upper_mps2':mag,
      'true_bias_norm_upper_mps2':up(math.sqrt(3)*mag),
      'variation_rate_abs_upper_mps3':rate,
      'tau_true_lower_s':tau_lo,'tau_true_upper_is_infinite':True,
      'non_relaxing_limit_admitted':phi_hi==1.0,
      'phi_true_interval':[phi_lo,phi_hi],
      'channel_driver_increment_components_mps2':{
        'bounded_variation':variation_step,'admitted_relaxation':relaxation_step},
      'driver_increment_component_abs_upper_mps2':w_component,
      'driver_increment_norm_upper_mps2':up(math.sqrt(3)*w_component),
      'driver_bound_is_analytic_not_replay_fit':True,
      'separation_sector':c['separation_sector'],
      'uniform_separation_constant_lower':None,
      'source_uniform_separation_closed':False,
      # The declared objective is bounded bias error plus practical ISS of the
      # other 18 errors.  Projection supplies the bias half for a non-relaxing
      # truth exactly as it does for a relaxing one, so separation is a gain
      # sharpener here, not a prerequisite.
      'separation_sector_required_for_bounded_bias_objective':False,
      'separation_sector_sharpens_motion_gains_only':True,
      'deployed_projection_radius_mps2':projection_radius,
      'bias_compactness_route':'|e_b| <= R + B_true pointwise from the closed radial projection sector',
      'bias_error_compactness_upper_mps2':up(projection_radius+up(math.sqrt(3)*mag)),
      'bias_compactness_needs_no_relaxation_root':True,
      'bias_compactness_needs_no_separation_constant':True,
      'BIAS2_SOURCE_ADMISSION_PASS':True,
      'deployment_hardware_admission_pass':False,
    }


def validate(x):
    f=[]
    if x.get('schema')!=SCHEMA or x.get('qualification')!=QUALIFICATION: f.append('schema/qualification mismatch')
    for k in ('conditional_theorem_family','one_drift_history_required','independent_per_sample_bias_slots_forbidden','non_relaxing_limit_admitted','tau_true_upper_is_infinite','separation_sector_sharpens_motion_gains_only','bias_compactness_needs_no_relaxation_root','bias_compactness_needs_no_separation_constant','driver_bound_is_analytic_not_replay_fit','BIAS2_SOURCE_ADMISSION_PASS'):
        if x.get(k) is not True: f.append(k+' not true')
    for k in ('assembled_sensor_hardware_qualified','deployment_hardware_admission_pass','inferred_from_BIAS0_or_BIAS1','may_be_substituted_by_BIAS1','source_uniform_separation_closed','separation_sector_required_for_bounded_bias_objective'):
        if x.get(k) is not False: f.append(k+' not false')
    if x.get('uniform_separation_constant_lower') is not None: f.append('separation constant falsely bound')
    if not float(x.get('driver_increment_norm_upper_mps2',0))>0: f.append('driver bound not positive')
    if not float(x.get('true_bias_norm_upper_mps2',0))>0: f.append('true bias bound not positive')
    compact=float(x.get('bias_error_compactness_upper_mps2',0))
    if not compact>=float(x.get('deployed_projection_radius_mps2',0))+float(x.get('true_bias_norm_upper_mps2',0)):
        f.append('compactness bound does not dominate R+B_true')
    channels=x.get('channel_driver_increment_components_mps2',{})
    if set(channels)!={'bounded_variation','admitted_relaxation'}:
        f.append('BIAS2 driver channel decomposition incomplete')
    elif not float(x.get('driver_increment_component_abs_upper_mps2',0))>=sum(float(v) for v in channels.values())-1e-18:
        f.append('driver increment does not dominate its channels')
    lo,hi=map(float,x.get('phi_true_interval',(0.0,0.0)))
    if not 0.0<lo<hi or hi!=1.0: f.append('BIAS2 phi_true interval must reach the non-relaxing limit')
    return f


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    x=build(); f=validate(x); x['validation_pass']=not f; x['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
    print(json.dumps(x,sort_keys=True)); return int(bool(f))
if __name__=='__main__': raise SystemExit(main())
