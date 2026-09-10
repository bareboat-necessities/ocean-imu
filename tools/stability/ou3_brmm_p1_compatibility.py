#!/usr/bin/env python3
"""Replay-free BRMM sea/RAO compatibility with the hard Normal-Live P1 branch.

The canonical theorem requires a *coupled same-history* sea/RAO realization as a
source-semantic requirement.  Earlier versions additionally used one analytical
G=4, gamma=1 witness to refute an independent sea x RAO Cartesian product under
the former 4 m/s^2 acceleration cap.  After widening the physical BRMM cap to
8 m/s^2 that witness (RMS >4.5 m/s^2) no longer exceeds the hard cap, so it may
no longer be used as a refutation.  This does not legalize detached sea/RAO
coordinates: the coupled source relation is required directly by the theorem.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_physical_admissibility as PHYS
import ou3_validated_transcendentals as VT
from ou3_interval import Interval,down,up

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/'tools'/'stability'/'ou3_proof_operating_domain.json'
DEFAULT_RAO=REPO/'tools'/'stability'/'ou3_brmm_directional_response_domain.json'
SCHEMA=2
QUALIFICATION='OU3_BRMM_COUPLED_SEA_RAO_P1_COMPATIBILITY_V2'

def _outward_point(x):return Interval.outward_bounds(float(x),float(x))
def exp_minus_five_quarters():
 e=VT.exp_point(-5.0/16.0);return e*e*e*e

def witness_mean_square_lower(gravity:float,gain:float=4.0)->float:
 sp=_outward_point(1.0/15.0);g=_outward_point(gravity);G=Interval.point(gain);pi=_outward_point(math.pi)
 J=_outward_point(8.0/9.0)*exp_minus_five_quarters();one_fifth=_outward_point(1.0/5.0)
 val=sp.square()*g.square()*G.square()*pi.square()/Interval.point(4.0);val=val*J/one_fifth
 return val.lo

def build(domain_path:Path=DEFAULT_DOMAIN,rao_path:Path=DEFAULT_RAO)->dict:
 domain=json.loads(Path(domain_path).read_text());rao=json.loads(Path(rao_path).read_text());gravity=float(domain['startup']['gravity_mps2']);cap=float(domain['normal_live']['non_gravitational_cog_acceleration_norm_upper_mps2']);response=rao['response_contract']
 gain_range=list(map(float,response['peak_translation_gain_range']));corner_range=list(map(float,response['rolloff_corner_hz_range']));pmin=float(response['high_frequency_rolloff_power_min'])
 finite=all(math.isfinite(v) for v in (*gain_range,*corner_range,pmin));inside=bool(finite and gain_range[0]<=4.0<=gain_range[1] and corner_range[0]<=1.2<=corner_range[1] and pmin<=2.0)
 if not inside:raise RuntimeError('declared RAO family no longer contains analytical witness')
 tp=8.0;hs=down(PHYS.significant_height_limit_from_peak_steepness(tp,gravity));
 if not PHYS.partition_admissible(hs,tp,gravity):raise RuntimeError('gamma=1 witness outside physical domain')
 ms=witness_mean_square_lower(gravity);rms=down(math.sqrt(ms));cap2=up(cap*cap);refuted=ms>cap2
 return {'schema':SCHEMA,'qualification':QUALIFICATION,'trajectory_replay_used':False,'filter_changed':False,
  'witness':{'declared_JONSWAP_gamma':1.0,'PM_is_JONSWAP_gamma_1_boundary':True,'witness_is_inside_declared_JONSWAP_gamma_interval_1_to_7':True,'T_p_s':tp,'H_s_m':hs,'RAO_gain':4.0,'RAO_corner_hz':1.2,'RAO_rolloff_power':2.0,'declared_RAO_gain_range':gain_range,'declared_RAO_corner_hz_range':corner_range,'declared_RAO_rolloff_power_min':pmin,'RAO_parameter_bounds_finite':finite,'witness_is_inside_declared_RAO_parameter_ranges':inside,'validated_acceleration_mean_square_lower_m2_s4':ms,'validated_acceleration_RMS_lower_mps2':rms,'P1_cap_squared_upper_m2_s4':cap2,'analytical_witness_exceeds_current_P1_cap':refuted,'all_nondyadic_witness_constants_outward_enclosed':True},
  'independent_cartesian_sea_x_RAO_domain_is_P1_sound':None,
  'cartesian_product_refuted_by_analytical_witness':refuted,
  'cartesian_product_refutation_required_for_canonical_BRMM':False,
  'coupled_BRMM_domain_required':True,
  'coupled_domain_requirement_source':'canonical same-history BRMM source semantics, independent of whether one counterexample exceeds the current cap',
  'coupled_domain_contract':{'sea_and_RAO_parameters_may_not_be_selected_independently':True,'P1_hard_acceleration_bound_checked_after_response':True,'PSD_or_RMS_bound_alone_sufficient_for_P1_admission':False,'finite_window_deterministic_response_certificate_required':True},
  'finite_window_realization_certificate_closed':False,'L_actual_sea_subset_Lhat_BRMM_closed':False,
  'next_obligation':'propagate the coupled same-history finite-window source relation; the old 4.635 m/s^2 witness is diagnostic only at the widened 8 m/s^2 cap'}

def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 if d.get('trajectory_replay_used') is not False or d.get('filter_changed') is not False:f.append('certificate must be replay free/non-invasive')
 if d.get('coupled_BRMM_domain_required') is not True:f.append('coupled BRMM domain requirement disappeared')
 if d.get('cartesian_product_refutation_required_for_canonical_BRMM') is not False:f.append('obsolete witness refutation made mandatory')
 w=d.get('witness',{})
 for k in ('PM_is_JONSWAP_gamma_1_boundary','witness_is_inside_declared_JONSWAP_gamma_interval_1_to_7','RAO_parameter_bounds_finite','witness_is_inside_declared_RAO_parameter_ranges','all_nondyadic_witness_constants_outward_enclosed'):
  if w.get(k) is not True:f.append(k+' not true')
 if not float(w.get('validated_acceleration_RMS_lower_mps2',0))>4.5:f.append('analytical witness regression lost')
 c=d.get('coupled_domain_contract',{})
 if c.get('sea_and_RAO_parameters_may_not_be_selected_independently') is not True:f.append('detached sea/RAO coordinates reintroduced')
 if c.get('finite_window_deterministic_response_certificate_required') is not True:f.append('finite-window response obligation disappeared')
 if d.get('finite_window_realization_certificate_closed') is not False or d.get('L_actual_sea_subset_Lhat_BRMM_closed') is not False:f.append('compatibility falsely closes deployment inclusion')
 return list(dict.fromkeys(f))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);ap.add_argument('--rao',type=Path,default=DEFAULT_RAO);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build(a.domain,a.rao);f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True));print(json.dumps({'cap':math.sqrt(d['witness']['P1_cap_squared_upper_m2_s4']),'witness_rms_lower':d['witness']['validated_acceleration_RMS_lower_mps2'],'witness_refutes_current_cap':d['cartesian_product_refuted_by_analytical_witness'],'coupled_required':d['coupled_BRMM_domain_required'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
