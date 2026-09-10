#!/usr/bin/env python3
"""Correct accelerometer attitude-covariance cap for the deployed residual.

A prior-independent r/|f|^2 transverse attitude cap is false for the shipping
filter because the accelerometer residual also contains latent acceleration and
accelerometer bias. What is valid for every prior and every unit e orthogonal to
f is the variational bound

 e^T P+_theta e <= (sigma_a^2 + 2 lambda_max(P_(a_w,b_a)))/|f|^2.

This module records that correction and derives the conditional cap from the
certified A21 covariance ceiling. It promotes nothing; the remaining obligation
is window observability/detectability strong enough to reduce the coupled
latent/bias covariance rather than a one-shot measurement argument.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_riccati_tube_factored as TUBE
REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_P4_ATTITUDE_MEASUREMENT_CAP_V1'

def up(x): return math.nextafter(float(x),math.inf)
def down(x): return math.nextafter(float(x),-math.inf)
def conditional_cap(sigma_a,f_lower,lam):
 if not(sigma_a>0 and f_lower>0 and lam>=0): raise ValueError('invalid cap inputs')
 return up(up(sigma_a*sigma_a+2*lam)/down(f_lower*f_lower))
def build(domain_path=DOMAIN):
 d=json.loads(Path(domain_path).read_text());tube=TUBE.build();tf=TUBE.validate_covariance_ceiling(tube)
 if tf: raise RuntimeError('covariance ceiling invalid: '+repr(tf))
 std=min(map(float,d['configured_runtime']['measurement_noise_std']['accelerometer_mps2']))
 fmin=float(d['normal_live']['specific_force_norm_lower_mps2'])
 diag=list(map(float,tube['modes']['A']['Pbar_diagonal_variance_upper']))
 latent=max(diag[15:18]);bias=max(diag[18:21]);joint=up(latent+bias);cap=conditional_cap(std,fmin,joint)
 retracted=down(std*std/up(fmin*fmin))
 return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','deployed_attitude_jacobian':'J_att=-skew(f_cog_b)','residual_also_carries_latent_and_bias_blocks':True,'prior_independent_transverse_cap_holds_for_deployed_filter':False,'retracted_prior_independent_cap_rad2':retracted,'conditional_cap_statement':'e^T P+_theta e <= (sigma_a^2 + 2 lambda_max(P_(a_w,b_a)))/|f|^2 for e perpendicular to f','variational_witness':'k=-(f x e)/|f|^2','conditional_on_joint_latent_bias_block':True,'latent_lambda_max':latent,'bias_lambda_max':bias,'joint_latent_bias_lambda_max':joint,'sigma_accelerometer_mps2':std,'specific_force_norm_lower_mps2':fmin,'conditional_transverse_cap_rad2':cap,'yaw_about_specific_force_capped_here':False,'one_shot_information_route_closes_correction_domain':False,'attitude_covariance_bound_needs_window_observability':True,'P4_promoted_here':False}
def validate(d):
 f=[]
 if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
 for k in ('residual_also_carries_latent_and_bias_blocks','conditional_on_joint_latent_bias_block','attitude_covariance_bound_needs_window_observability'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('prior_independent_transverse_cap_holds_for_deployed_filter','yaw_about_specific_force_capped_here','one_shot_information_route_closes_correction_domain','P4_promoted_here'):
  if d.get(k) is not False:f.append(k+' not false')
 if not float(d.get('conditional_transverse_cap_rad2',0))>0:f.append('conditional cap invalid')
 if not float(d.get('conditional_transverse_cap_rad2',0))>float(d.get('retracted_prior_independent_cap_rad2',0)):f.append('conditional cap failed to expose coupling')
 return f
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args();d=build();f=validate(d)
 if a.output:a.output.write_text(json.dumps({**d,'validation_pass':not f,'validation_failures':f},indent=2,sort_keys=True)+'\n')
 print(json.dumps({**d,'validation_pass':not f,'validation_failures':f},sort_keys=True));raise SystemExit(bool(f))
