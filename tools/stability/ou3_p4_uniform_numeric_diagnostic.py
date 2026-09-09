#!/usr/bin/env python3
"""Numerical diagnostic for the source-uniform P4 whole-word route.

No promotion.  Reports the actual theorem-level Live covariance seed/tube,
complete-word H18 information lower, scheduler recurrence and finite-bias A21
headroom so the nonlinear master can be scaled against the right quantities
instead of the deliberately prior-free P3 trace upper.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_riccati_tube as TUBE
import ou3_brmm_h18_information_composition as HINFO
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_brmm_riccati_metric_p3 as P3

QUALIFICATION='OU3_P4_SOURCE_UNIFORM_NUMERIC_DIAGNOSTIC_V1'

def build():
    live=LIVE.build(); tube=TUBE.build(); info=HINFO.build(); sched=SCHED.build(); p3=P3.build()
    bad={
      'live':LIVE.validate(live),'tube':TUBE.validate(tube),'info':HINFO.validate(info),
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
          'relative_Riccati_injection_margin_lower':float(row['relative_Riccati_injection_margin_lower'])}
    return {
      'qualification':QUALIFICATION,'promotes_P4':False,'trajectory_replay_used':False,
      'canonical_source':p3['canonical_source'],
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

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,indent=2,sort_keys=True));return 0
if __name__=='__main__':raise SystemExit(main())
