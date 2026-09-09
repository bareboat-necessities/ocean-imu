#!/usr/bin/env python3
"""Explicit conditional BIAS1 family for bounded-bias P4.

This is a theorem/source-family qualification, not assembled-sensor hardware
qualification.  One root and one parameter history are retained over a word;
there are no independent per-sample bias slots.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
DEFAULT=REPO/'tools/stability/ou3_p4_closure_domain.json'
DT=0.005


def up(x): return math.nextafter(float(x), math.inf)
def down(x): return math.nextafter(float(x), -math.inf)


def build(path: Path=DEFAULT):
    c=json.loads(Path(path).read_text())['BIAS1_family']
    b=float(c['root_component_abs_upper_mps2'])
    a=float(c['sinusoid_component_amplitude_abs_upper_mps2'])
    d=float(c['deterministic_mismatch_component_abs_upper_mps2'])
    dl=float(c['deterministic_mismatch_derivative_component_abs_upper_mps3'])
    tlo,thi=map(float,c['tau_true_s']); plo,phi=map(float,c['sinusoid_period_s'])
    if not (0<tlo<=thi and 0<plo<=phi and min(b,a,d,dl)>=0): raise RuntimeError('invalid BIAS1 family')
    omega_hi=up(2*math.pi/plo)
    phi_lo=down(math.exp(-DT/tlo)); phi_hi=up(math.exp(-DT/thi))
    # beta=b0 exp(-t/tau)+a sin(omega t+phase)+d_det(t).
    # In beta_i=phi beta_{i-1}+w_i, the homogeneous root cancels exactly.
    # |sin q_i-phi sin q_{i-1}| <= 2 sin(omega*h/2)+(1-phi).
    sine_step=up(2*math.sin(up(omega_hi*DT/2)) + (1-phi_lo))
    # |d_i-phi d_{i-1}| <= |d_i-d_{i-1}|+(1-phi)|d_{i-1}|.
    det_step=up(dl*DT + (1-phi_lo)*d)
    w_component=up(a*sine_step + det_step)
    beta_component=up(b+a+d)
    return {
      'qualification':'OU3_P4_CONDITIONAL_BIAS1_FAMILY_V1',
      'conditional_theorem_family':True,'assembled_sensor_hardware_qualified':False,
      'model':c['model'],'one_root_one_parameter_history_required':True,
      'independent_per_sample_bias_slots_forbidden':True,'dt_s':DT,
      'root_component_abs_upper_mps2':b,'tau_true_s':[tlo,thi],
      'sinusoid_component_amplitude_abs_upper_mps2':a,'sinusoid_period_s':[plo,phi],
      'deterministic_mismatch_component_abs_upper_mps2':d,
      'deterministic_mismatch_derivative_component_abs_upper_mps3':dl,
      'phi_true_interval':[phi_lo,phi_hi],
      'true_bias_component_abs_upper_mps2':beta_component,
      'true_bias_norm_upper_mps2':up(math.sqrt(3)*beta_component),
      'driver_increment_component_abs_upper_mps2':w_component,
      'driver_increment_norm_upper_mps2':up(math.sqrt(3)*w_component),
      'driver_bound_is_analytic_not_replay_fit':True,
      'existing_probe_member': b>=0.08 and a>=0.015 and tlo<=1200<=thi and plo<=600<=phi,
      'BIAS1_SOURCE_ADMISSION_PASS':True,
      'deployment_hardware_admission_pass':False,
    }


def validate(x):
    f=[]
    for k in ('conditional_theorem_family','one_root_one_parameter_history_required','independent_per_sample_bias_slots_forbidden','driver_bound_is_analytic_not_replay_fit','existing_probe_member','BIAS1_SOURCE_ADMISSION_PASS'):
        if x.get(k) is not True: f.append(k+' not true')
    if x.get('assembled_sensor_hardware_qualified') is not False or x.get('deployment_hardware_admission_pass') is not False: f.append('hardware falsely qualified')
    if not float(x.get('driver_increment_norm_upper_mps2',0))>0: f.append('driver bound not positive')
    return f


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    x=build(); f=validate(x); x['validation_pass']=not f; x['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
    print(json.dumps(x,sort_keys=True)); return int(bool(f))
if __name__=='__main__': raise SystemExit(main())
