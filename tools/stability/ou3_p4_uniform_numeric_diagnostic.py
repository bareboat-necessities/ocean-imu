#!/usr/bin/env python3
"""Numerical diagnostic for the source-uniform P4 whole-word route.

No promotion. Reports the actual theorem-level Live covariance seed, canonical
endpoint-referenced covariance envelope, complete-word H18 information lower,
scheduler recurrence and finite-bias A21 headroom.  The historical scalar
moving-Riccati injection margin is retained as diagnostic data only; P4 no
longer treats that abandoned scalar architecture as a prerequisite.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_h18_information_composition as HINFO
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_brmm_riccati_metric_p3 as P3

QUALIFICATION='OU3_P4_SOURCE_UNIFORM_NUMERIC_DIAGNOSTIC_V2'

def build():
    live=LIVE.build(); tube=TUBE.build(); info=HINFO.build(); sched=SCHED.build(); p3=P3.build()
    bad={
      'live':LIVE.validate(live),'endpoint_covariance_envelope':TUBE.validate_covariance_ceiling(tube),'info':HINFO.validate(info),
      'scheduler':SCHED.validate(sched),'P3':P3.validate(p3)}
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError(str(bad))
    tri=info['triangular_information_composition']
    modes={}
    for mode,key in [('H18','H'),('A21','A')]:
        row=tube['modes'][key]
        diag=[float(x) for x in row['Pbar_diagonal_variance_upper']]
        modes[mode]={
          'Pbar_diagonal_variance_upper':diag,
          'Pbar_trace_diagonal_sum_upper':math.nextafter(sum(diag),math.inf),
          'Pbar_max_diagonal_upper':max(diag),
          'historical_scalar_relative_Riccati_injection_margin_lower_diagnostic_only':float(row['relative_Riccati_injection_margin_lower']),
          'historical_scalar_useful_margin_pass':bool(row['useful_margin_pass'])}
    return {
      'qualification':QUALIFICATION,'promotes_P4':False,'trajectory_replay_used':False,
      'canonical_source':p3['canonical_source'],
      'canonical_endpoint_referenced_covariance_envelope_consumed':True,
      'failed_moving_Riccati_relative_margin_consumed_as_P4_premise':False,
      'historical_scalar_moving_Riccati_architecture_abandoned':True,
      'scheduler_uniform_max_gap_samples':sched['certified_uniform_max_gap_samples'],
      'scheduler_uniform_max_gap_s':sched['certified_uniform_max_gap_s'],
      'minimum_S_updates_in_600_sample_word':600//int(sched['certified_uniform_max_gap_samples']),
      'H18_eta6_information_lower':info['eta6_information_lower'],
      'H18_directional_translation_information_lower':info['directional_translation_information_lower'],
      'H18_complete_word_information_lambda_min_lower':tri['D_H18_lambda_min_lower'],
      'H18_complete_word_information_legacy_scalar_lower':tri['legacy_scalarized_D_H18_lambda_min_lower_diagnostic'],
      'H18_accel_aw_cross_norm_sq_upper':info['accelerometer_translation_cross_norm_squared_upper'],
      'modes':modes,
      'A21_bias_homogeneous_contraction_gap_lower':p3['conditional_composition']['A21_bias_homogeneous_contraction_gap_lower'],
      'A21_detectability_asymptotic_word_energy_gap_lower':p3['conditional_composition']['A21_detectability_asymptotic_word_energy_gap_lower'],
    }

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('canonical_endpoint_referenced_covariance_envelope_consumed','historical_scalar_moving_Riccati_architecture_abandoned'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('promotes_P4','trajectory_replay_used','failed_moving_Riccati_relative_margin_consumed_as_P4_premise'):
        if d.get(k) is not False:f.append(k+' not false')
    for mode,m in d.get('modes',{}).items():
        x=float(m.get('historical_scalar_relative_Riccati_injection_margin_lower_diagnostic_only',math.nan))
        if not(math.isfinite(x) and x>=0):f.append(mode+' historical scalar margin invalid')
        if not m.get('Pbar_diagonal_variance_upper'):f.append(mode+' Pbar missing')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,indent=2,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
