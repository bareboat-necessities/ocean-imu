#!/usr/bin/env python3
"""Explicit conditional BIAS0 family for bounded-bias P4.

BIAS0 is the composite physical-composition accelerometer-bias family

    b_true(t)=beta_gm(t)+b0+K_T*dT(t)+K_strain*strain(t)+d_nonGM(t),

i.e. the shipping-centered truth with every uncompensated deterministic term
retained.  It is a mandatory family in its own right: it is never inferred
from, substituted by, or reduced to BIAS1's one-root one-parameter graph.

Every admitted sequence satisfies the exact recurrence

    b_i = phi_true*b_{i-1} + w_i,    w_i := b_i - phi_true*b_{i-1},

for any factor, so the content of this module is the outward bound on w.  The
declared Gauss-Markov channel root is used as phi_true, which cancels that
channel's homogeneous part exactly and charges every non-GM channel to the
driver increment through its declared magnitude and rate:

    |c_i - phi*c_{i-1}| <= |c_i - c_{i-1}| + (1-phi)*|c_{i-1}|
                        <= rate*dt + (1-phi_lo)*magnitude.

The Gauss-Markov channel contributes a declared PATHWISE increment cap.  A
stationary OU PSD supplies no hard pathwise cap, so that cap is an explicit
conditional hypothesis of the family, not a consequence of the sigma/tau pair;
the sigma-multiple below is descriptive only.  This is a theorem/source-family
qualification, not assembled-sensor hardware qualification.  One composite
physical history is retained over a word; there are no independent per-sample
bias slots.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
DEFAULT=REPO/'tools/stability/ou3_p4_closure_domain.json'
DT=0.005
SCHEMA=1
QUALIFICATION='OU3_P4_CONDITIONAL_BIAS0_FAMILY_V1'


def up(x): return math.nextafter(float(x), math.inf)
def down(x): return math.nextafter(float(x), -math.inf)


def _deterministic_channel(magnitude,rate,one_minus_phi):
    """Outward |c_i-phi c_{i-1}| for a bounded-magnitude bounded-rate channel."""
    if magnitude<0 or rate<0: raise RuntimeError('deterministic channel bounds must be nonnegative')
    return up(up(rate*DT)+up(one_minus_phi*magnitude))


def build(path: Path=DEFAULT):
    c=json.loads(Path(path).read_text())['BIAS0_family']
    gm=float(c['gauss_markov_component_abs_upper_mps2'])
    gm_w=float(c['gauss_markov_pathwise_increment_component_abs_upper_mps2'])
    gm_sigma=float(c['gauss_markov_stationary_sigma_descriptive_mps2'])
    tlo,thi=map(float,c['gauss_markov_tau_true_s'])
    b0=float(c['turn_on_offset_component_abs_upper_mps2'])
    kt=float(c['thermal_coefficient_abs_upper_mps2_per_K'])
    dT=float(c['thermal_deviation_abs_upper_K'])
    dTdot=float(c['thermal_deviation_rate_abs_upper_K_s'])
    st=float(c['strain_component_abs_upper_mps2'])
    stdot=float(c['strain_component_rate_abs_upper_mps3'])
    dn=float(c['non_gauss_markov_component_abs_upper_mps2'])
    dndot=float(c['non_gauss_markov_rate_abs_upper_mps3'])
    if not (0<tlo<=thi): raise RuntimeError('invalid BIAS0 Gauss-Markov tau interval')
    if min(gm,gm_w,gm_sigma,b0,kt,dT,dTdot,st,stdot,dn,dndot)<0: raise RuntimeError('invalid BIAS0 family bounds')
    if c.get('gauss_markov_pathwise_increment_cap_is_declared_hypothesis') is not True:
        raise RuntimeError('BIAS0 pathwise increment cap must stay a declared hypothesis')

    phi_lo=down(math.exp(-DT/tlo)); phi_hi=up(math.exp(-DT/thi))
    one_minus_phi=up(1.0-phi_lo)

    thermal_mag=up(kt*dT); thermal_rate=up(kt*dTdot)
    turn_on_step=_deterministic_channel(b0,0.0,one_minus_phi)
    thermal_step=_deterministic_channel(thermal_mag,thermal_rate,one_minus_phi)
    strain_step=_deterministic_channel(st,stdot,one_minus_phi)
    non_gm_step=_deterministic_channel(dn,dndot,one_minus_phi)
    w_component=up(gm_w+up(turn_on_step+up(thermal_step+up(strain_step+non_gm_step))))
    beta_component=up(gm+up(b0+up(thermal_mag+up(st+dn))))

    # Descriptive only: the one-step OU innovation standard deviation implied by
    # the declared sigma/tau pair, taken at the shortest admitted tau so the
    # reported multiple is the least favourable one.  It does not derive the
    # pathwise cap, which stays a declared hypothesis of the family.
    innovation_sigma=up(gm_sigma*math.sqrt(up(1.0-math.exp(-2.0*DT/tlo))))
    sigma_multiple=up(gm_w/innovation_sigma) if innovation_sigma>0 else math.inf

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'conditional_theorem_family':True,'assembled_sensor_hardware_qualified':False,
      'model':c['model'],'family_role':c['family_role'],
      'one_composite_history_required':True,
      'independent_per_sample_bias_slots_forbidden':True,
      'inferred_from_BIAS1':False,'may_be_substituted_by_BIAS1':False,
      'dt_s':DT,
      'gauss_markov_component_abs_upper_mps2':gm,
      'gauss_markov_tau_true_s':[tlo,thi],
      'gauss_markov_pathwise_increment_component_abs_upper_mps2':gm_w,
      'gauss_markov_pathwise_increment_cap_is_declared_hypothesis':True,
      'stationary_OU_PSD_used_to_derive_pathwise_cap':False,
      'descriptive_largest_one_step_innovation_sigma_mps2':innovation_sigma,
      'descriptive_pathwise_cap_sigma_multiple':sigma_multiple,
      'turn_on_offset_component_abs_upper_mps2':b0,
      'thermal_component_abs_upper_mps2':thermal_mag,
      'thermal_component_rate_abs_upper_mps3':thermal_rate,
      'strain_component_abs_upper_mps2':st,
      'strain_component_rate_abs_upper_mps3':stdot,
      'non_gauss_markov_component_abs_upper_mps2':dn,
      'non_gauss_markov_rate_abs_upper_mps3':dndot,
      'phi_true_interval':[phi_lo,phi_hi],
      'channel_driver_increment_components_mps2':{
        'gauss_markov':gm_w,'turn_on_offset':turn_on_step,'thermal':thermal_step,
        'strain':strain_step,'non_gauss_markov':non_gm_step},
      'true_bias_component_abs_upper_mps2':beta_component,
      'true_bias_norm_upper_mps2':up(math.sqrt(3)*beta_component),
      'driver_increment_component_abs_upper_mps2':w_component,
      'driver_increment_norm_upper_mps2':up(math.sqrt(3)*w_component),
      'driver_bound_is_analytic_not_replay_fit':True,
      'BIAS0_SOURCE_ADMISSION_PASS':True,
      'deployment_hardware_admission_pass':False,
    }


def validate(x):
    f=[]
    if x.get('schema')!=SCHEMA or x.get('qualification')!=QUALIFICATION: f.append('schema/qualification mismatch')
    for k in ('conditional_theorem_family','one_composite_history_required','independent_per_sample_bias_slots_forbidden','gauss_markov_pathwise_increment_cap_is_declared_hypothesis','driver_bound_is_analytic_not_replay_fit','BIAS0_SOURCE_ADMISSION_PASS'):
        if x.get(k) is not True: f.append(k+' not true')
    for k in ('assembled_sensor_hardware_qualified','deployment_hardware_admission_pass','inferred_from_BIAS1','may_be_substituted_by_BIAS1','stationary_OU_PSD_used_to_derive_pathwise_cap'):
        if x.get(k) is not False: f.append(k+' not false')
    if not float(x.get('driver_increment_norm_upper_mps2',0))>0: f.append('driver bound not positive')
    if not float(x.get('true_bias_norm_upper_mps2',0))>0: f.append('true bias bound not positive')
    channels=x.get('channel_driver_increment_components_mps2',{})
    if set(channels)!={'gauss_markov','turn_on_offset','thermal','strain','non_gauss_markov'}:
        f.append('BIAS0 driver channel decomposition incomplete')
    elif not float(x.get('driver_increment_component_abs_upper_mps2',0))>=sum(float(v) for v in channels.values())-1e-18:
        f.append('driver increment does not dominate its channels')
    lo,hi=map(float,x.get('phi_true_interval',(0.0,0.0)))
    if not 0.0<lo<=hi<=1.0: f.append('invalid phi_true interval')
    return f


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    x=build(); f=validate(x); x['validation_pass']=not f; x['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
    print(json.dumps(x,sort_keys=True)); return int(bool(f))
if __name__=='__main__': raise SystemExit(main())
