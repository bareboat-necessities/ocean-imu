#!/usr/bin/env python3
"""Source-uniform same-cell Joseph correction-domain certificate for P4.

For accepted vector Joseph updates, the exact same-cell identity
K S K^T=P^- - P^+ <= P^- bounds the attitude correction without a rowwise K
box.  The source-uniform covariance ceiling and exact residual geometry then
bound accelerometer and magnetometer correction magnitudes.

The S=0 pseudo measurement is different. Its physical residual is

    y_S = -S_hat = e_S - S_true,

not e_S alone.  The large session integration constant is a shared source/gauge
column and cancels at every literal prefix.  Therefore this module explicitly
REFUSES to substitute the legacy hard e_S radius (300 m*s) for y_S.  The S-event
correction/chart target must be certified in the production augmented prefix
graph using the centered joint selector e_S-S_true and the same-cell d=K y.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_exact_reset_transport as RESET
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_S_origin_literal_prefix_invariant as ORIGIN

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_P4_SAME_CELL_JOSEPH_CORRECTION_DOMAIN_V3'

def up(x):return math.nextafter(float(x),math.inf)
def down(x):return math.nextafter(float(x),-math.inf)

def _residual_bounds(domain,entry,dynamic):
    live=domain['normal_live'];q=float(entry['coordinate_radii']['attitude_cayley_norm']);theta=2*math.atan(.5*q)
    rmi=up(2*math.sin(.5*theta));fmax=float(live['specific_force_norm_upper_mps2']);mmax=float(live['magnetic_vector_norm_upper_uT'])
    aw=float(entry['coordinate_radii']['latent_acceleration_norm_mps2']);ba=float(entry['coordinate_radii'].get('accelerometer_bias_error_norm_mps2',0.0))
    noise=domain['configured_runtime']['measurement_noise_std'];acc_std=min(map(float,noise['accelerometer_mps2']));mag_std=min(map(float,noise['magnetometer_uT']))
    if not(acc_std>0 and mag_std>0):raise RuntimeError('configured vector measurement std lost positivity')
    inv=dynamic['dynamic_invariant'];rs_lo=float(inv['R_S_applied'][0]);rs_std_min=down(rs_lo*.72)
    if not rs_std_min>0:raise RuntimeError('actual applied R_S lower lost positivity')
    acc=up(rmi*fmax+aw+ba);mag=up(rmi*mmax)
    return {'attitude_cayley_norm_upper':q,'attitude_angle_rad':theta,'rotation_minus_identity_norm_upper':rmi,
      'accelerometer_residual_norm_upper':acc,'magnetometer_residual_norm_upper':mag,
      'accelerometer_measurement_std_lower':acc_std,'magnetometer_measurement_std_lower':mag_std,
      'accelerometer_Rinv_energy_upper':up(acc*acc/down(acc_std*acc_std)),
      'magnetometer_Rinv_energy_upper':up(mag*mag/down(mag_std*mag_std)),
      'S_zero_Rinv_energy_upper':None,'S_zero_residual_static_error_ball_used':False,
      'S_zero_residual_definition':'e_S-S_true=-S_hat in the centered joint prefix graph',
      'S_zero_requires_centered_prefix_graph':True,'actual_RS_std_lower_with_horizontal_factor':rs_std_min}

def build(domain_path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve();domain=json.loads(path.read_text());tube=TUBE.build(path);entry=ENTRY.build();dynamic=DYNAMIC.build(path);cayley=CAYLEY.build(path);origin=ORIGIN.build()
    bad={'endpoint_covariance_envelope':TUBE.validate_covariance_ceiling(tube),'entry':ENTRY.validate(entry),'dynamic':DYNAMIC.validate(dynamic),'cayley':CAYLEY.validate(cayley),'origin':ORIGIN.validate(origin)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('same-cell correction-domain prerequisites failed: '+repr(bad))
    residual=_residual_bounds(domain,entry,dynamic);modes={}
    for mode,key in (('H18','H'),('A21','A')):
        p=list(map(float,tube['modes'][key]['Pbar_diagonal_variance_upper']));patt=up(sum(p[:3]));events={}
        for name,energy in (('accelerometer',residual['accelerometer_Rinv_energy_upper']),('magnetometer',residual['magnetometer_Rinv_energy_upper'])):
            delta=up(math.sqrt(max(0.0,up(patt*float(energy)))))
            events[name]={'attitude_covariance_trace_upper':patt,'residual_Rinv_energy_upper':float(energy),
              'same_cell_attitude_correction_norm_upper':delta,'inside_reset_utility_domain':delta<=RESET.CAYLEY_MONOTONE_NORM_MAX,
              'inside_one_radian':delta<=1.0,'requires_centered_prefix_graph':False}
        events['S_zero']={'attitude_covariance_trace_upper':patt,'residual_Rinv_energy_upper':None,
          'same_cell_attitude_correction_norm_upper':None,'inside_reset_utility_domain':None,'inside_one_radian':None,
          'requires_centered_prefix_graph':True,'static_hard_eS_radius_used':False}
        bounded=[(n,e['same_cell_attitude_correction_norm_upper']) for n,e in events.items() if e['same_cell_attitude_correction_norm_upper'] is not None]
        worst_name,worst=max(bounded,key=lambda x:x[1]);vector_closed=all(events[n]['inside_reset_utility_domain'] for n in ('accelerometer','magnetometer'))
        modes[mode]={'attitude_covariance_trace_upper':patt,'events':events,'limiting_bounded_vector_event':worst_name,
          'bounded_vector_event_correction_norm_upper':worst,'vector_events_reset_utility_domain_closed':vector_closed,
          'S_zero_centered_prefix_graph_required':True,'correction_norm_upper':None,'limiting_event':'S_zero:centered-prefix-required',
          'reset_utility_domain_closed':False}
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_Joseph_covariance_identity_consumed':True,'identity':'K*S*K^T=Pminus-Pplus<=Pminus','attitude_block_only_used':True,
      'rowwise_K_bound_used':False,'independent_K_box_used':False,'hard_entry_residual_geometry_consumed_for_vector_events':True,
      'source_uniform_endpoint_covariance_envelope_consumed':True,'failed_moving_Riccati_relative_margin_consumed':False,
      'actual_applied_RS_lower_available_for_integrated_S_event':True,'shared_S_origin_prefix_invariant_consumed':origin['finite_literal_prefix_composition_closed_by_induction'],
      'legacy_300m_s_eS_radius_used_as_S_zero_residual':False,'S_zero_reset_domain_is_integrated_prefix_obligation':True,
      'residual_bounds':residual,'modes':modes,
      'all_bounded_vector_events_inside_exact_reset_utility_domain':all(m['vector_events_reset_utility_domain_closed'] for m in modes.values()),
      'all_events_inside_exact_reset_utility_domain':False,'production_same_graph_correction_domain_promoted_here':False,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_Joseph_covariance_identity_consumed','attitude_block_only_used','hard_entry_residual_geometry_consumed_for_vector_events','source_uniform_endpoint_covariance_envelope_consumed','actual_applied_RS_lower_available_for_integrated_S_event','shared_S_origin_prefix_invariant_consumed','S_zero_reset_domain_is_integrated_prefix_obligation'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rowwise_K_bound_used','independent_K_box_used','failed_moving_Riccati_relative_margin_consumed','legacy_300m_s_eS_radius_used_as_S_zero_residual','all_events_inside_exact_reset_utility_domain','production_same_graph_correction_domain_promoted_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('residual_bounds',{}).get('S_zero_Rinv_energy_upper') is not None:f.append('S zero scalar energy reintroduced')
    for mode,m in d.get('modes',{}).items():
        s=m.get('events',{}).get('S_zero',{})
        if s.get('requires_centered_prefix_graph') is not True or s.get('same_cell_attitude_correction_norm_upper') is not None:f.append(mode+' S event not fail-closed on centered graph')
        for ev in ('accelerometer','magnetometer'):
            x=float(m['events'][ev]['same_cell_attitude_correction_norm_upper'])
            if not(math.isfinite(x) and x>=0):f.append(mode+' '+ev+' correction invalid')
    return f

def main():
    p=argparse.ArgumentParser();p.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build(a.domain);f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'vector_events_closed':d['all_bounded_vector_events_inside_exact_reset_utility_domain'],'S_event':'centered-prefix-required','legacy_S_radius_used':d['legacy_300m_s_eS_radius_used_as_S_zero_residual'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
