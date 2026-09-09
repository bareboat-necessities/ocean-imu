#!/usr/bin/env python3
"""Non-promoting exact-reset absorption diagnostic for the signed P4 ledger.

Consumes the full declared 45-degree hard entry chart and the source-uniform
Joseph correction magnitude ceiling.  The exact reset transport then bounds
rho in

    e_reset = G(d) t + rho,   ||G(d)^-1||_2 = 1.

We report dimensionless ratios ||rho||/||d|| and ||rho||/(q ||d||).  Small
ratios support a local signed-ledger treatment; large ratios reject that route
without changing the theorem domain.  No covariance inverse or packet-count
bound is introduced here, so this diagnostic cannot promote P4.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_exact_reset_transport as RESET

QUALIFICATION='OU3_P4_RESET_SIGNED_ABSORPTION_DIAGNOSTIC_V1'

def build():
    c=COEFF.build();e=ENTRY.build();bad={'coeff':COEFF.validate(c),'entry':ENTRY.validate(e)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('reset diagnostic prerequisites failed: '+repr(bad))
    q=float(e['coordinate_radii']['attitude_cayley_norm']);modes={}
    for mode in ('H18','A21'):
        delta=float(c['modes'][mode]['attitude_correction_norm_upper'])
        row=RESET.reset_defect_bound(q,delta)
        rho=float(row['reset_attitude_defect_norm_upper'])
        modes[mode]={
          'correction_norm_upper':delta,'reset_attitude_defect_norm_upper':rho,
          'rho_over_delta':rho/delta if delta>0 else 0.0,
          'rho_over_q_delta':rho/(q*delta) if q>0 and delta>0 else 0.0,
          'cayley_composition_denominator_lower':row['cayley_composition_denominator_lower'],
          'injected_cayley_norm_upper':row['injected_cayley_norm_upper'],
          'chart_safe':row['chart_safe'],'limiting_correction_event':c['modes'][mode]['limiting_correction_event']}
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','full_declared_entry_cayley_norm':q,
      'same_source_uniform_correction_magnitude_envelope_consumed':True,'exact_reset_transport_consumed':True,
      'reset_inverse_operator_norm_upper':1.0,'modes':modes,'covariance_inverse_scalarized_here':False,
      'packet_count_multiplier_used':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_source_uniform_correction_magnitude_envelope_consumed','exact_reset_transport_consumed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('covariance_inverse_scalarized_here','packet_count_multiplier_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('reset_inverse_operator_norm_upper')!=1.0:f.append('reset inverse norm changed')
    for mode,m in d.get('modes',{}).items():
        if m.get('chart_safe') is not True:f.append(mode+' chart unsafe')
        for k in ('correction_norm_upper','reset_attitude_defect_norm_upper','rho_over_delta','rho_over_q_delta','cayley_composition_denominator_lower'):
            if not (math.isfinite(float(m.get(k,math.nan))) and float(m[k])>=0):f.append(mode+' '+k+' invalid')
        if not float(m['cayley_composition_denominator_lower'])>0:f.append(mode+' denominator nonpositive')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'modes':d['modes'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
